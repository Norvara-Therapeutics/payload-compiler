"""Independent consistency checks for persisted optimization audit trails."""

import hashlib
import re

from Bio.Data import CodonTable
from Bio.Seq import Seq

from payload_compiler.analysis import scan_splice_signals
from payload_compiler.exceptions import InputError
from payload_compiler.models import CodingRegion, OptimizeReport, Payload
from payload_compiler.profiles import HOST_PROFILE_VERSIONS, HOST_PROFILES
from payload_compiler.validation import coding_region_positions, extract_cds, validate_payload


def _cpg_positions(sequence: str) -> tuple[int, ...]:
    return tuple(match.start() + 1 for match in re.finditer(r"(?=CG)", sequence))


def _translation(payload: Payload, region: CodingRegion) -> str:
    coding = extract_cds(payload, region)
    protein = str(Seq(coding).translate(table=region.translation_table))
    table = CodonTable.unambiguous_dna_by_id[region.translation_table]
    if coding[:3] in table.start_codons:
        protein = "M" + protein[1:]
    return protein


def _codon(sequence: str, positions: tuple[int, ...], strand: int) -> str:
    bases = "".join(sequence[position] for position in positions)
    return bases if strand == 1 else str(Seq(bases).complement())


def _mutation_positions(
    region: CodingRegion, start: int, end: int
) -> tuple[tuple[int, ...], int]:
    transcript_positions = coding_region_positions(region)
    for offset in range(0, len(transcript_positions), 3):
        positions = transcript_positions[offset : offset + 3]
        if len(positions) == 3 and min(positions) + 1 == start and max(positions) + 1 == end:
            return positions, offset // 3
    raise InputError(f"Mutation coordinates {start}:{end} do not identify a codon in {region.name}")


def verify_optimization_audit(
    source: Payload, optimized: Payload, audit: OptimizeReport
) -> None:
    """Recompute hashes, CpGs, translations, and every base edit in an audit."""
    if source.identifier != audit.identifier or optimized.identifier != audit.identifier:
        raise InputError("Sequence identifiers do not match the optimization audit")
    if audit.host not in HOST_PROFILES or audit.host_profile != HOST_PROFILE_VERSIONS[audit.host]:
        raise InputError("Host profile name or version does not match a bundled profile")
    if len(source.sequence) != len(optimized.sequence):
        raise InputError("Source and optimized sequence lengths differ")
    source_digest = hashlib.sha256(source.sequence.encode("utf-8")).hexdigest()
    optimized_digest = hashlib.sha256(optimized.sequence.encode("utf-8")).hexdigest()
    if source_digest != audit.source_sha256:
        raise InputError("Source sequence does not match the audit SHA-256")
    if optimized_digest != audit.optimized_sha256:
        raise InputError("Optimized sequence does not match the audit SHA-256")

    if source.coding_regions and source.coding_regions != audit.coding_regions:
        raise InputError("Annotated source CDS definitions do not match the optimization audit")
    audited_source = source.model_copy(update={"coding_regions": audit.coding_regions})
    audited_optimized = optimized.model_copy(update={"coding_regions": audit.coding_regions})
    for label, payload in (("source", audited_source), ("optimized", audited_optimized)):
        validation = validate_payload(payload, require_cds=True)
        if not validation.valid:
            raise InputError(
                f"{label.capitalize()} sequence does not satisfy the audited CDS model"
            )

    source_cpgs = _cpg_positions(source.sequence)
    optimized_cpgs = _cpg_positions(optimized.sequence)
    if len(source_cpgs) != audit.cpg_count_before:
        raise InputError("Source CpG count does not match the optimization audit")
    if optimized_cpgs != audit.residual_cpg_positions:
        raise InputError("Optimized CpG positions do not match the optimization audit")

    regions = {region.name: region for region in audit.coding_regions}
    replay = list(source.sequence)
    for mutation in audit.mutations:
        region = regions.get(mutation.region)
        if region is None:
            raise InputError(f"Mutation references unknown CDS '{mutation.region}'")
        positions, codon_index = _mutation_positions(region, mutation.start, mutation.end)
        current = "".join(replay)
        if _codon(current, positions, region.strand) != mutation.codon_before:
            raise InputError(f"Mutation {mutation.index} codon-before value is inconsistent")
        if len(_cpg_positions(current)) != mutation.cpg_count_before:
            raise InputError(f"Mutation {mutation.index} CpG-before value is inconsistent")
        if not set(change.position - 1 for change in mutation.base_changes).issubset(positions):
            raise InputError(f"Mutation {mutation.index} changes a base outside its codon")
        for change in mutation.base_changes:
            position = change.position - 1
            if not 0 <= position < len(replay) or replay[position] != change.before:
                raise InputError(f"Mutation {mutation.index} base-change trace is inconsistent")
            replay[position] = change.after
        current = "".join(replay)
        if _codon(current, positions, region.strand) != mutation.codon_after:
            raise InputError(f"Mutation {mutation.index} codon-after value is inconsistent")
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        before_aa = (
            "M"
            if codon_index == 0 and mutation.codon_before in table.start_codons
            else table.forward_table.get(mutation.codon_before)
        )
        after_aa = (
            "M"
            if codon_index == 0 and mutation.codon_after in table.start_codons
            else table.forward_table.get(mutation.codon_after)
        )
        if before_aa != mutation.amino_acid or after_aa != mutation.amino_acid:
            raise InputError(f"Mutation {mutation.index} is not synonymous as audited")
        profile = HOST_PROFILES[audit.host]
        if mutation.host_codon_frequency_before != profile.get(mutation.codon_before, 0.0):
            raise InputError(f"Mutation {mutation.index} source codon frequency is inconsistent")
        if mutation.host_codon_frequency_after != profile.get(mutation.codon_after, 0.0):
            raise InputError(f"Mutation {mutation.index} optimized codon frequency is inconsistent")
        if len(_cpg_positions(current)) != mutation.cpg_count_after:
            raise InputError(f"Mutation {mutation.index} CpG-after value is inconsistent")

    if "".join(replay) != optimized.sequence:
        raise InputError("Mutation trace does not reproduce the optimized sequence")
    source_translations = {
        region.name: _translation(audited_source, region) for region in audit.coding_regions
    }
    optimized_translations = {
        region.name: _translation(audited_optimized, region) for region in audit.coding_regions
    }
    translations_preserved = source_translations == optimized_translations
    if translations_preserved != audit.protein_translation_preserved:
        raise InputError("Translation-preservation claim does not match the audited sequences")
    if not translations_preserved:
        raise InputError("Optimized sequence changes an audited CDS translation")

    expected_regressions: list[str] = []
    if audit.after_lint.configuration.splice_scan:
        before_splice = {
            (signal.kind, signal.start, signal.end, signal.strand)
            for signal in scan_splice_signals(source.sequence)
        }
        after_splice = {
            (signal.kind, signal.start, signal.end, signal.strand)
            for signal in scan_splice_signals(optimized.sequence)
        }
        if new_splice := after_splice - before_splice:
            expected_regressions.append(
                f"optimization introduced {len(new_splice)} splice-consensus signal(s)"
            )
    before_folds = {
        (window.start, window.end): window for window in audit.before_lint.metrics.folding_windows
    }
    newly_flagged = tuple(
        window
        for window in audit.after_lint.metrics.folding_windows
        if window.flagged
        and (
            (previous := before_folds.get((window.start, window.end))) is None
            or not previous.flagged
        )
    )
    if newly_flagged:
        expected_regressions.append(
            f"optimization introduced {len(newly_flagged)} newly flagged RNAfold window(s)"
        )
    if tuple(expected_regressions) != audit.safety_regressions:
        raise InputError("Biological safety regression claims do not match the audited sequences")
