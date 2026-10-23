"""Structural checks only; this module never changes the input sequence."""

import hashlib

from Bio.Data import CodonTable
from Bio.Seq import Seq

from payload_compiler import __version__
from payload_compiler.models import CodingRegion, Finding, Payload, ValidationReport

DNA_BASES = frozenset("ACGT")
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


def coding_region_positions(region: CodingRegion) -> tuple[int, ...]:
    """Return genomic positions in CDS transcript order."""
    parts = region.parts or ((region.start, region.end),)
    ordered = sorted(parts, reverse=region.strand == -1)
    if region.strand == 1:
        return tuple(position for start, end in ordered for position in range(start, end))
    return tuple(position for start, end in ordered for position in range(end - 1, start - 1, -1))


def extract_cds(payload: Payload, region: CodingRegion) -> str:
    if not 0 <= region.start < region.end <= len(payload.sequence):
        raise ValueError("CDS coordinates are outside the sequence")
    positions = coding_region_positions(region)
    sequence = Seq("".join(payload.sequence[position] for position in positions))
    return str(sequence if region.strand == 1 else sequence.complement())


def validate_payload(payload: Payload, *, require_cds: bool = False) -> ValidationReport:
    findings: list[Finding] = []

    def add(code, message, *, severity="error", region=None, positions=()):
        findings.append(
            Finding(
                code=code, message=message, severity=severity, region=region, positions=positions
            )
        )

    seq = payload.sequence
    invalid = tuple(i + 1 for i, base in enumerate(seq) if base not in IUPAC_DNA)
    ambiguous = tuple(i + 1 for i, base in enumerate(seq) if base in IUPAC_DNA - DNA_BASES)
    if invalid:
        add(
            "INVALID_BASE",
            "Sequence contains characters outside the IUPAC DNA alphabet",
            positions=invalid,
        )
    if ambiguous:
        add(
            "AMBIGUOUS_BASE",
            "Sequence contains ambiguous IUPAC DNA bases",
            severity="warning",
            positions=ambiguous,
        )
    if not payload.coding_regions:
        add(
            "NO_CDS",
            "No CDS annotated; no reading frame was inferred",
            severity="error" if require_cds else "warning",
        )
    seen_names: set[str] = set()
    for region in payload.coding_regions:
        if region.name in seen_names:
            add("DUPLICATE_CDS_NAME", "CDS names must be unique", region=region.name)
        seen_names.add(region.name)
        if region.end > len(seq):
            add("CDS_OUT_OF_BOUNDS", "CDS extends beyond the sequence", region=region.name)
            continue
        positions = coding_region_positions(region)
        if any(position >= len(seq) for position in positions):
            add("CDS_OUT_OF_BOUNDS", "CDS part extends beyond the sequence", region=region.name)
            continue
        coding = extract_cds(payload, region)
        if len(coding) % 3:
            add("CDS_FRAME", "CDS length must be divisible by three", region=region.name)
            continue
        if set(coding) - DNA_BASES:
            add("CDS_UNRESOLVED_BASE", "CDS requires unambiguous A/C/G/T bases", region=region.name)
            continue
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        if set(table.stop_codons) & set(table.forward_table):
            add(
                "CONTEXT_DEPENDENT_TABLE",
                "Context-dependent stop codons are unsupported",
                region=region.name,
            )
            continue
        protein = str(Seq(coding).translate(table=region.translation_table))
        if "*" in protein[:-1]:
            add("CDS_INTERNAL_STOP", "CDS contains an internal stop codon", region=region.name)
        if region.annotated_translation is not None:
            # GenBank complete CDS annotations translate alternative starts to methionine.
            annotated_protein = protein
            if coding[:3] in table.start_codons:
                annotated_protein = "M" + protein[1:]
            if annotated_protein.removesuffix("*") != region.annotated_translation.removesuffix(
                "*"
            ):
                add(
                    "TRANSLATION_MISMATCH",
                    "CDS does not match its /translation annotation",
                    region=region.name,
                )
    region_position_sets = {
        region.name: set(coding_region_positions(region)) for region in payload.coding_regions
    }
    for i, left in enumerate(payload.coding_regions):
        for right in payload.coding_regions[i + 1 :]:
            if region_position_sets[left.name] & region_position_sets[right.name]:
                add("OVERLAPPING_CDS", f"{left.name} overlaps {right.name}", severity="warning")
    return ValidationReport(
        tool_version=__version__,
        identifier=payload.identifier,
        source_format=payload.source_format,
        length=len(seq),
        sequence_sha256=hashlib.sha256(seq.encode("utf-8")).hexdigest(),
        cds_count=len(payload.coding_regions),
        findings=tuple(findings),
    )
