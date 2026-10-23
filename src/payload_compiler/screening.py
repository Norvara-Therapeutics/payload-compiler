"""Transparent sequence heuristics; scores do not establish biological activity."""

import re

from Bio.Seq import Seq

from payload_compiler.models import MotifDefinition, Payload, ScreenConfig, ScreenHit
from payload_compiler.validation import coding_region_positions

ASSEMBLY_MOTIFS = {
    "golden-gate-bsai": (("BsaI", "GGTCTC"),),
    "golden-gate-bsmbi": (("BsmBI", "CGTCTC"),),
    "biobrick": (("EcoRI", "GAATTC"), ("XbaI", "TCTAGA"),
                 ("SpeI", "ACTAGT"), ("PstI", "CTGCAG")),
}


def screen_sequence(payload: Payload, config: ScreenConfig) -> tuple[ScreenHit, ...]:
    if not config.enabled:
        return ()
    # Local import avoids a module initialization cycle with the linter.
    from payload_compiler.lint import _motif_matches

    hits = []
    sequence = payload.sequence
    length = len(sequence)
    intended_starts = {
        (coding_region_positions(r)[0], "+" if r.strand == 1 else "-")
        for r in payload.coding_regions
    }
    for strand, oriented in (("+", sequence), ("-", str(Seq(sequence).reverse_complement()))):
        # `strand` is bound as a default so the closure cannot capture the loop
        # variable by reference and mislabel reverse-strand hits.
        def add(kind, name, start, end, evidence, score=None, strand=strand):
            left, right = (start + 1, end) if strand == "+" else (length - end + 1, length - start)
            hits.append(ScreenHit(kind=kind, name=name, start=left, end=right,
                                  strand=strand, evidence=evidence, score=score))

        for match in re.finditer(r"(?=(AATAAA|ATTAAA))", oriented):
            add("polya", "poly(A) signal-like hexamer", match.start(), match.start() + 6,
                "Consensus hexamer only; transcript context and cleavage are not predicted")
        for match in re.finditer(r"(?=(TATA[AT]A[AT]))", oriented):
            add("promoter", "TATA-like motif", match.start(), match.start() + 7,
                "Motif only; does not establish internal promoter activity")
        for match in re.finditer(r"(?=(ATG|CTG|GTG|TTG))", oriented):
            start = match.start()
            genomic_start = start if strand == "+" else length - start - 1
            if (genomic_start, strand) in intended_starts:
                continue
            strong = start >= 3 and start + 3 < length and (
                oriented[start - 3] in "AG" and oriented[start + 3] == "G"
            )
            if strong:
                add("alternative_start", "Kozak-context start-like site", start, start + 3,
                    "ATG or near-cognate codon with purine at -3 and G at +4; "
                    "not initiation evidence")
            if match.group(1) != "ATG":
                continue
            for end in range(start + 3, length - 2, 3):
                codon = oriented[end:end + 3]
                if set(codon) - set("ACGT"):
                    break
                if codon in ("TAA", "TAG", "TGA"):
                    if (end - start) // 3 >= config.orf_min_codons:
                        add("orf", "unannotated-start ORF", start, end + 3,
                            "Standard-code ATG-to-stop ORF; annotation-relative, "
                            "not expression evidence")
                    break
        for matrix in config.matrices:
            width = len(matrix.weights)
            for start in range(length - width + 1):
                fragment = oriented[start:start + width]
                if set(fragment) - set("ACGT"):
                    continue
                score = sum(row[base] for row, base in zip(matrix.weights, fragment, strict=True))
                if score >= matrix.threshold:
                    add(matrix.kind, matrix.name, start, start + width,
                        f"User matrix {matrix.version}; {matrix.source}; {matrix.citation}", score)
    for match in re.finditer(r"([ACGT])\1{" + str(config.homopolymer_min - 1) + r",}", sequence):
        hits.append(ScreenHit(
            kind="homopolymer", name="homopolymer", start=match.start() + 1,
            end=match.end(), strand="both", evidence="Synthesis-complexity heuristic"))
    for width in range(2, config.tandem_repeat_max_unit + 1):
        pattern = (r"([ACGT]{" + str(width) + r"})\1{"
                   + str(config.tandem_repeat_min_copies - 1) + r",}")
        for match in re.finditer(pattern, sequence):
            if len(set(match.group(1))) == 1:
                continue
            hits.append(ScreenHit(kind="tandem_repeat", name=f"{width}-nt tandem repeat",
                                  start=match.start() + 1, end=match.end(), strand="both",
                                  evidence="Synthesis-complexity heuristic, not vendor acceptance"))
    for standard in config.assembly_standards:
        for enzyme, site in ASSEMBLY_MOTIFS[standard]:
            for match in _motif_matches(sequence, (MotifDefinition(name=enzyme, sequence=site),)):
                hits.append(ScreenHit(
                    kind="assembly_site", name=f"{standard}:{enzyme}",
                    start=match.start, end=match.end, strand=match.strand,
                    evidence=(f"{enzyme} site ({site}) is disallowed inside a {standard} insert; "
                              "cloning constraint, not a biological finding")))
    for match in _motif_matches(sequence, config.motifs):
        hits.append(ScreenHit(
            kind="custom_motif", name=match.name, start=match.start,
            end=match.end, strand=match.strand, evidence="User-configured IUPAC motif"))
    return tuple(sorted(hits, key=lambda h: (h.start, h.end, h.kind, h.name, h.strand)))


def hit_keys(hits):
    return {(h.kind, h.name, h.start, h.end, h.strand) for h in hits}
