"""Deterministic synonymous CpG depletion with an auditable mutation trace."""

import hashlib
import platform
import re
import sys
from importlib.metadata import version

from Bio.Data import CodonTable
from Bio.Seq import Seq

from payload_compiler import __version__
from payload_compiler.analysis import scan_splice_signals
from payload_compiler.design import (
    codon_constraints_pass,
    codon_metrics,
    get_profile,
    profile_digest,
    protected_positions,
    resolve_protection,
    vector_budget,
)
from payload_compiler.exceptions import InputError
from payload_compiler.lint import _motif_matches, lint_payload
from payload_compiler.models import (
    BaseChange,
    CodingRegion,
    CodonProfile,
    LintConfig,
    Mutation,
    OptimizeReport,
    OptimizeSettings,
    Payload,
    RuntimeMetadata,
    SearchResult,
)
from payload_compiler.screening import hit_keys, screen_sequence
from payload_compiler.validation import (
    DNA_BASES,
    coding_region_positions,
    extract_cds,
    validate_payload,
)


def _cpg_positions(sequence: str) -> tuple[int, ...]:
    return tuple(match.start() + 1 for match in re.finditer(r"(?=CG)", sequence))


def _cpg_count(sequence: str) -> int:
    return len(_cpg_positions(sequence))


def _coding_translations(
    payload: Payload, sequence: str, regions: tuple[CodingRegion, ...] | None = None
) -> dict[str, str]:
    result = {}
    for region in regions if regions is not None else payload.coding_regions:
        if region.end > len(sequence):
            continue
        fragment = Seq(extract_cds(payload.model_copy(update={"sequence": sequence}), region))
        if len(fragment) % 3 or set(str(fragment)) - DNA_BASES:
            raise InputError(f"CDS '{region.name}' cannot be translated unambiguously")
        translation = str(fragment.translate(table=region.translation_table))
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        if str(fragment[:3]) in table.start_codons:
            translation = "M" + translation[1:]
        result[region.name] = translation
    return result


def _local_cpg_count(sequence: str, changed_positions: tuple[int, ...]) -> int:
    starts = {
        start
        for position in changed_positions
        for start in (position - 1, position)
        if 0 <= start < len(sequence) - 1
    }
    return sum(sequence[start : start + 2] == "CG" for start in starts)


def _codon_positions(region: CodingRegion, codon_offset: int) -> tuple[int, int, int]:
    positions = coding_region_positions(region)[codon_offset : codon_offset + 3]
    if len(positions) != 3:
        raise InputError(f"CDS '{region.name}' contains an incomplete codon")
    return positions


def _codon_in_sequence(sequence: str, positions: tuple[int, int, int], strand: int) -> str:
    bases = "".join(sequence[position] for position in positions)
    return bases if strand == 1 else str(Seq(bases).complement())


def _write_codon(sequence: str, positions: tuple[int, int, int], strand: int, codon: str) -> str:
    updated = list(sequence)
    if strand == 1:
        genomic_bases = codon
    else:
        genomic_bases = str(Seq(codon).complement())
    for position, base in zip(positions, genomic_bases, strict=True):
        updated[position] = base
    return "".join(updated)


def _translate_codon(codon: str, table: CodonTable.CodonTable, *, first: bool) -> str | None:
    if first and codon in table.start_codons:
        return "M"
    if codon in table.stop_codons:
        return "*"
    return table.forward_table.get(codon)


def _codon_amino_acid(sequence: str, region: CodingRegion, codon_index: int) -> str | None:
    positions = _codon_positions(region, codon_index * 3)
    codon = _codon_in_sequence(sequence, positions, region.strand)
    table = CodonTable.unambiguous_dna_by_id[region.translation_table]
    return _translate_codon(codon, table, first=codon_index == 0)


def _translations_preserved(
    payload: Payload,
    before: str,
    after: str,
    changed_positions: tuple[int, ...],
) -> bool:
    """Check only codons intersecting an edit, including all overlapping CDSs."""
    for region in payload.coding_regions:
        region_positions = coding_region_positions(region)
        position_offsets = {position: offset for offset, position in enumerate(region_positions)}
        offsets = [
            position_offsets[position]
            for position in changed_positions
            if position in position_offsets
        ]
        for codon_index in {offset // 3 for offset in offsets}:
            if _codon_amino_acid(before, region, codon_index) != _codon_amino_acid(
                after, region, codon_index
            ):
                return False
    return True


def _motif_keys(sequence: str, config: LintConfig) -> set[tuple[str, int, int, str]]:
    return {
        (match.name, match.start, match.end, match.strand)
        for match in _motif_matches(sequence, config.forbidden_motifs)
    }


def _splice_keys(sequence: str) -> set[tuple[str, int, int, str]]:
    return {
        (signal.kind, signal.start, signal.end, signal.strand)
        for signal in scan_splice_signals(sequence)
    }


def _candidate_codons(
    payload: Payload,
    sequence: str,
    cpg_positions: tuple[int, ...],
    current_cpg_count: int,
    before_motifs: set[tuple[str, int, int, str]],
    before_splice_signals: set[tuple[str, int, int, str]],
    lint_config: LintConfig,
    profile: CodonProfile,
    settings: OptimizeSettings,
    before_screen_keys: set,
) -> list[tuple[tuple, str, Mutation]]:
    frequencies = profile.frequencies
    protected = protected_positions(settings)
    current_motif_keys = before_motifs
    candidate_codons: set[tuple[str, int]] = set()

    # A CpG can only be changed by a codon containing its C or G. Map those
    # genomic bases into each unambiguous, in-frame CDS; overlapping CDS are
    # permitted only when the resulting sequence preserves every translation.
    for position in cpg_positions:
        for genomic_position in (position - 1, position):
            for region in payload.coding_regions:
                region_positions = coding_region_positions(region)
                if genomic_position in region_positions:
                    transcript_offset = region_positions.index(genomic_position)
                    codon_index = transcript_offset // 3
                    candidate_codons.add((region.name, codon_index))

    regions = {region.name: region for region in payload.coding_regions}
    candidates: list[tuple[tuple, str, Mutation]] = []
    for region_name, codon_index in sorted(candidate_codons):
        region = regions[region_name]
        if codon_index == 0:
            continue
        codon_offset = codon_index * 3
        positions = _codon_positions(region, codon_offset)
        original_codon = _codon_in_sequence(sequence, positions, region.strand)
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        amino_acid = _codon_amino_acid(sequence, region, codon_index)
        if amino_acid is None or amino_acid == "*":
            continue
        synonymous = sorted(
            codon
            for codon in (set(table.forward_table) | set(table.start_codons))
            if codon != original_codon
            and _translate_codon(codon, table, first=codon_index == 0) == amino_acid
        )
        for replacement in synonymous:
            updated = _write_codon(sequence, positions, region.strand, replacement)
            changed_positions = tuple(
                position for position in positions if sequence[position] != updated[position]
            )
            if protected.intersection(changed_positions):
                continue
            if settings.cpg_policy == "strict_subset" and not set(
                _cpg_positions(updated)
            ).issubset(cpg_positions):
                continue
            next_cpg_count = (
                current_cpg_count
                - _local_cpg_count(sequence, changed_positions)
                + _local_cpg_count(updated, changed_positions)
            )
            reduction = current_cpg_count - next_cpg_count
            if reduction <= 0:
                continue
            if not _translations_preserved(payload, sequence, updated, changed_positions):
                continue
            if not current_motif_keys.issuperset(_motif_keys(updated, lint_config)):
                continue
            if lint_config.splice_scan and not before_splice_signals.issuperset(
                _splice_keys(updated)
            ):
                continue
            candidate = payload.model_copy(update={"sequence": updated})
            if not codon_constraints_pass(
                codon_metrics(candidate, profile, settings.rare_codon_threshold), settings
            ):
                continue
            if settings.prevent_screen_regressions and not hit_keys(
                screen_sequence(candidate, lint_config.screens)
            ).issubset(before_screen_keys):
                continue
            base_changes = []
            for position in changed_positions:
                if sequence[position] != updated[position]:
                    base_changes.append(
                        BaseChange(
                            position=position + 1,
                            before=sequence[position],
                            after=updated[position],
                        )
                    )
            host_before = frequencies[original_codon]
            host_after = frequencies[replacement]
            mutation = Mutation(
                index=1,
                region=region.name,
                amino_acid=amino_acid,
                start=min(positions) + 1,
                end=max(positions) + 1,
                strand=region.strand,
                codon_before=original_codon,
                codon_after=replacement,
                host_codon_frequency_before=host_before,
                host_codon_frequency_after=host_after,
                cpg_count_before=current_cpg_count,
                cpg_count_after=next_cpg_count,
                reason="synonymous codon replacement reduces whole-payload CpG count",
                base_changes=tuple(base_changes),
            )
            # Lower tuples sort first: maximize CpG reduction and preferred host
            # codon frequency, then minimize changes and use stable lexical ties.
            score = (
                -reduction,
                -host_after,
                len(base_changes),
                min(positions),
                region.name,
                replacement,
            )
            candidates.append((score, updated, mutation))
    return candidates


def _oversize_message(payload: Payload, lint_config: LintConfig) -> str:
    """Explain an over-budget input in component terms instead of refusing blankly.

    Synonymous substitution is length-preserving, so an oversized payload can only
    be fixed by changing what the vector contains. When a vector architecture is
    declared the message names the accounted unit, states whether both ITRs are
    inside the submitted sequence, and lists the shortest component swap plan the
    declared catalog permits.
    """
    length = len(payload.sequence)
    lines = [
        f"Input is {length} nt, over the {lint_config.max_length} nt limit; "
        "synonymous substitution preserves length and cannot shorten a payload."
    ]
    if lint_config.vector is None:
        lines.append(
            "No --vector architecture was supplied, so no component-level size accounting "
            "is available. Supply one to see ITR accounting and shorter-component options."
        )
        return " ".join(lines)
    budget = vector_budget(payload, lint_config.vector)
    if lint_config.vector.representation == "itr_to_itr_unit":
        lines.append("Both ITRs are inside the submitted sequence and are counted in the unit.")
    else:
        lines.append(
            f"Submitted sequence is a cassette excluding ITRs; "
            f"{budget['omitted_itr_overhead']} nt of declared ITR overhead was added, "
            f"giving an accounted unit of {budget['accounted_unit_length']} nt."
        )
    lines.append(
        f"{budget['genome_type']} budget {budget['limit']} nt, headroom {budget['headroom']} nt."
    )
    plan = budget["plan"]
    if plan["changes"] and plan["fits"]:
        swaps = "; ".join(
            f"{change['component']} -> {change['alternative']} ({change['length']} nt)"
            for change in plan["changes"]
        )
        lines.append(
            f"Shortest declared plan that fits ({plan['estimated_length']} nt): {swaps}. "
            "These are size estimates only; this tool does not perform the substitution."
        )
    else:
        replaceable = [
            component["name"]
            for component in budget["components"]
            if component["shorter_alternatives"]
        ]
        lines.append(
            "No declared alternative catalog brings the unit within budget"
            + (f"; replaceable components: {', '.join(replaceable)}." if replaceable else ".")
        )
    return " ".join(lines)


def optimize_payload(
    payload: Payload,
    *,
    host: str = "human",
    lint_config: LintConfig | None = None,
    settings: OptimizeSettings | None = None,
    profile: CodonProfile | None = None,
) -> tuple[Payload, OptimizeReport]:
    """Remove CpGs using single-codon synonymous changes with deterministic ranking."""
    profile = get_profile(host, profile)
    host = profile.name
    lint_config = lint_config or LintConfig()
    settings = resolve_protection(payload, settings or OptimizeSettings(), lint_config.vector)
    baseline_validation = validate_payload(payload, require_cds=True)
    if not baseline_validation.valid:
        problems = "; ".join(
            f"{finding.code}: {finding.message}"
            for finding in baseline_validation.findings
            if finding.severity == "error"
        )
        raise InputError(f"Input must pass sequence validation before optimization ({problems})")
    if len(payload.sequence) > lint_config.max_length:
        raise InputError(_oversize_message(payload, lint_config))

    initial_translations = _coding_translations(payload, payload.sequence)
    before_lint = lint_payload(payload, lint_config)
    before_metrics = codon_metrics(payload, profile, settings.rare_codon_threshold)
    sequence = payload.sequence
    original_cpg_count = _cpg_count(sequence)
    current_cpg_count = original_cpg_count
    mutations: list[Mutation] = []
    before_motifs = _motif_keys(sequence, lint_config)
    before_splice_signals = _splice_keys(sequence) if lint_config.splice_scan else set()
    before_screen_keys = hit_keys(before_lint.metrics.screen_hits)
    steps = 0

    while True:
        positions = _cpg_positions(sequence)
        if not positions:
            break
        candidates = _candidate_codons(
            payload,
            sequence,
            positions,
            current_cpg_count,
            before_motifs,
            before_splice_signals,
            lint_config,
            profile,
            settings,
            before_screen_keys,
        )
        steps += 1
        if not candidates:
            break
        candidates.sort(key=lambda item: item[0])
        _, sequence, mutation = candidates[0]
        mutations.append(mutation.model_copy(update={"index": len(mutations) + 1}))
        current_cpg_count = mutation.cpg_count_after

    search = SearchResult(
        solver="greedy", status="locally_optimal", explored_states=steps,
        objective="CpG reduction per step, descending codon preference, changed bases, coordinates",
        proof_scope=(
            "No improving permitted single-codon move under the "
            f"{settings.cpg_policy} CpG policy; no global optimum guarantee"
        ),
    )
    edit_semantics = "sequential"
    if settings.solver == "global":
        from payload_compiler.search import global_search

        original_positions = set(_cpg_positions(payload.sequence))
        protected = protected_positions(settings)
        def feasible(candidate):
            if settings.cpg_policy == "strict_subset" and not set(
                _cpg_positions(candidate)
            ).issubset(original_positions):
                return False
            if settings.cpg_policy == "count_only" and len(
                _cpg_positions(candidate)
            ) > len(original_positions):
                return False
            if any(candidate[p] != payload.sequence[p] for p in protected):
                return False
            if _coding_translations(payload, candidate) != initial_translations:
                return False
            if not _motif_keys(candidate, lint_config).issubset(before_motifs):
                return False
            if lint_config.splice_scan and not _splice_keys(candidate).issubset(
                before_splice_signals
            ):
                return False
            proposed = payload.model_copy(update={"sequence": candidate})
            if not codon_constraints_pass(
                codon_metrics(proposed, profile, settings.rare_codon_threshold), settings
            ):
                return False
            if settings.prevent_screen_regressions and not hit_keys(
                screen_sequence(proposed, lint_config.screens)
            ).issubset(before_screen_keys):
                return False
            candidate_lint = lint_payload(proposed, lint_config)
            old_folds = {(w.start, w.end) for w in before_lint.metrics.folding_windows if w.flagged}
            return candidate_lint.valid and not any(
                w.flagged and (w.start, w.end) not in old_folds
                for w in candidate_lint.metrics.folding_windows
            )

        seed = sequence if feasible(sequence) else (
            payload.sequence if feasible(payload.sequence) else None
        )
        result, search = global_search(payload, profile, settings, feasible, seed)
        sequence = result if result is not None else payload.sequence
        mutations = _joint_mutations(payload, sequence, profile)
        edit_semantics = "joint"

    optimized_payload = payload.model_copy(update={"sequence": sequence})
    final_translations = _coding_translations(optimized_payload, sequence)
    if final_translations != initial_translations:
        raise InputError("Internal safety check failed: CDS translation changed")
    after_lint = lint_payload(optimized_payload, lint_config)
    remaining_positions = _cpg_positions(sequence)
    warnings = []
    safety_regressions: list[str] = []
    after_metrics = codon_metrics(optimized_payload, profile, settings.rare_codon_threshold)
    if not codon_constraints_pass(after_metrics, settings):
        safety_regressions.append("optimized sequence does not meet configured codon constraints")
    if settings.prevent_screen_regressions and not hit_keys(
        after_lint.metrics.screen_hits
    ).issubset(before_screen_keys):
        safety_regressions.append("optimization introduced new sequence-screen findings")
    introduced_cpg = tuple(
        sorted(set(_cpg_positions(sequence)) - set(_cpg_positions(payload.sequence)))
    )
    if settings.cpg_policy == "strict_subset" and introduced_cpg:
        raise InputError("Internal safety check failed: introduced a new CpG position")
    if introduced_cpg:
        warnings.append(
            f"cpg-policy count-only: {len(introduced_cpg)} CpG(s) introduced at new position(s) "
            f"{', '.join(str(position) for position in introduced_cpg)}; "
            "total count still decreased"
        )
    if lint_config.splice_scan:
        new_splice_signals = _splice_keys(sequence) - before_splice_signals
        if new_splice_signals:
            safety_regressions.append(
                f"optimization introduced {len(new_splice_signals)} splice-consensus signal(s)"
            )
    before_folds = {
        (window.start, window.end): window for window in before_lint.metrics.folding_windows
    }
    newly_flagged_folds = tuple(
        window
        for window in after_lint.metrics.folding_windows
        if window.flagged
        and (
            (previous := before_folds.get((window.start, window.end))) is None
            or not previous.flagged
        )
    )
    if newly_flagged_folds:
        safety_regressions.append(
            f"optimization introduced {len(newly_flagged_folds)} newly flagged RNAfold window(s)"
        )
    if remaining_positions:
        warnings.append(
            f"{len(remaining_positions)} CpG site(s) remain outside editable CDS contexts "
            "or were retained under the selected constraints and search budget"
        )
    if search.status in ("best_found", "no_feasible_found"):
        warnings.append("Search budget exhausted; no global optimality or infeasibility claim")
    if lint_config.vector is None:
        warnings.append(
            "No vector architecture supplied: length was checked against a flat "
            f"{lint_config.max_length}-nt limit with no ITR or component accounting"
        )
    else:
        if lint_config.vector.representation == "cassette":
            warnings.append(
                "Submitted sequence excludes both ITRs; "
                f"{sum(lint_config.vector.omitted_itr_lengths or ())} nt of declared ITR "
                "overhead was added to the accounted unit"
            )
        else:
            warnings.append(
                "Submitted sequence is an ITR-to-ITR unit; both ITRs are counted "
                "and held immutable"
            )

    report = OptimizeReport(
        tool_version=__version__,
        identifier=payload.identifier,
        source_sha256=baseline_validation.sequence_sha256,
        optimized_sha256=hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
        host=host,
        host_profile=profile.version,
        max_length=lint_config.max_length,
        cpg_depletion_requested=True,
        protein_translation_preserved=True,
        cpg_count_before=original_cpg_count,
        cpg_count_after=len(remaining_positions),
        cpg_reduction=original_cpg_count - len(remaining_positions),
        coding_regions=payload.coding_regions,
        mutations=tuple(mutations),
        residual_cpg_positions=remaining_positions,
        introduced_cpg_positions=introduced_cpg,
        before_lint=before_lint,
        after_lint=after_lint,
        safety_checks_passed=not safety_regressions,
        safety_regressions=tuple(safety_regressions),
        runtime=RuntimeMetadata(
            python=sys.version.split()[0],
            platform=platform.platform(),
            biopython=version("biopython"),
            pydantic=version("pydantic"),
            seqfold=version("seqfold"),
            typer=version("typer"),
        ),
        warnings=tuple(warnings),
        settings=settings,
        profile=profile,
        profile_sha256=profile_digest(profile),
        search=search,
        edit_semantics=edit_semantics,
        codon_metrics_before=before_metrics,
        codon_metrics_after=after_metrics,
    )
    return optimized_payload, report


def _joint_mutations(payload: Payload, sequence: str, profile: CodonProfile) -> list[Mutation]:
    """Describe every changed CDS codon against the source, applied simultaneously.

    Overlapping CDS rows may share base edits. They must agree; these rows are
    explanations of one final transaction, not independently applicable steps.
    """
    mutations = []
    for region in payload.coding_regions:
        positions = coding_region_positions(region)
        for offset in range(0, len(positions), 3):
            triplet = positions[offset:offset + 3]
            before = _codon_in_sequence(payload.sequence, triplet, region.strand)
            after = _codon_in_sequence(sequence, triplet, region.strand)
            if before == after:
                continue
            mutations.append(Mutation(
                index=len(mutations) + 1, region=region.name,
                amino_acid=_codon_amino_acid(payload.sequence, region, offset // 3),
                start=min(triplet) + 1, end=max(triplet) + 1, strand=region.strand,
                codon_before=before, codon_after=after,
                host_codon_frequency_before=profile.frequencies[before],
                host_codon_frequency_after=profile.frequencies[after],
                cpg_count_before=_cpg_count(payload.sequence), cpg_count_after=_cpg_count(sequence),
                reason="joint synonymous solution; apply all rows simultaneously",
                base_changes=tuple(BaseChange(position=p + 1, before=payload.sequence[p],
                                              after=sequence[p])
                                   for p in triplet if payload.sequence[p] != sequence[p]),
            ))
    return mutations
