"""Typed contracts. Internal coordinates are zero-based, half-open."""

import re
from typing import Literal

from Bio.Data import CodonTable
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class Model(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class CodingRegion(Model):
    name: str = Field(min_length=1)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    strand: Literal[1, -1] = 1
    translation_table: int = 1
    annotated_translation: str | None = None
    # Empty means one contiguous interval spanning start:end. For joined CDS
    # features these are zero-based, half-open genomic intervals.
    parts: tuple[tuple[int, int], ...] = ()

    @model_validator(mode="after")
    def check_region(self) -> "CodingRegion":
        if self.end <= self.start:
            raise ValueError("CDS end must be greater than start")
        if self.translation_table not in CodonTable.unambiguous_dna_by_id:
            raise ValueError(f"Unknown translation table: {self.translation_table}")
        if self.parts:
            ordered = sorted(self.parts)
            if any(part_end <= part_start for part_start, part_end in ordered):
                raise ValueError("CDS parts must have positive length")
            if any(
                left[1] > right[0]
                for left, right in zip(ordered, ordered[1:], strict=False)
            ):
                raise ValueError("CDS parts must not overlap")
            if self.start != ordered[0][0] or self.end != ordered[-1][1]:
                raise ValueError("CDS start/end must bound all CDS parts")
        return self


class FeatureAnnotation(Model):
    """Retained GenBank feature metadata; location uses Biopython's notation."""

    kind: str
    location: str | None
    qualifiers: dict[str, tuple[str, ...]] = Field(default_factory=dict)


class Payload(Model):
    identifier: str = Field(min_length=1)
    description: str = ""
    sequence: str = Field(min_length=1)
    source_format: Literal["fasta", "genbank", "embl"]
    topology: Literal["linear", "circular"] = "linear"
    coding_regions: tuple[CodingRegion, ...] = ()
    features: tuple[FeatureAnnotation, ...] = ()


class Finding(Model):
    code: str
    severity: Literal["error", "warning"]
    message: str
    region: str | None = None
    # Findings are user-facing: positions are one-based.
    positions: tuple[int, ...] = ()


class ValidationReport(Model):
    schema_version: str = "1.0"
    tool_version: str
    identifier: str
    source_format: Literal["fasta", "genbank", "embl"]
    length: int
    sequence_sha256: str
    cds_count: int
    findings: tuple[Finding, ...]

    @property
    def valid(self) -> bool:
        return not any(item.severity == "error" for item in self.findings)


class MotifDefinition(Model):
    name: str = Field(min_length=1)
    sequence: str = Field(min_length=1)

    @field_validator("name")
    @classmethod
    def normalize_name(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("motif name cannot be blank")
        return value

    @field_validator("sequence")
    @classmethod
    def normalize_sequence(cls, value: str) -> str:
        value = value.strip().upper()
        invalid = sorted(set(value) - set("ACGTRYSWKMBDHVN"))
        if not value:
            raise ValueError("motif sequence cannot be blank")
        if invalid:
            raise ValueError(f"unsupported motif DNA symbols: {', '.join(invalid)}")
        return value


class ProtectedRegion(Model):
    """User-facing one-based inclusive immutable interval."""

    start: int = Field(gt=0)
    end: int = Field(gt=0)
    name: str = "protected"

    @model_validator(mode="after")
    def ordered(self):
        if self.end < self.start:
            raise ValueError("protected interval end must not precede start")
        return self


class ScoringMatrix(Model):
    name: str = Field(min_length=1)
    kind: Literal["splice", "promoter", "transcription_factor"]
    source: str = Field(min_length=1)
    version: str = Field(min_length=1)
    citation: str = Field(min_length=1)
    # User-supplied additive weights, not inferred probabilities.
    weights: tuple[dict[str, float], ...] = Field(min_length=1)
    threshold: float = Field(allow_inf_nan=False)

    @model_validator(mode="after")
    def complete_weights(self):
        import math

        if any(set(row) != set("ACGT") for row in self.weights):
            raise ValueError("each matrix position requires A, C, G, T weights")
        if any(not math.isfinite(v) for row in self.weights for v in row.values()):
            raise ValueError("matrix weights must be finite")
        return self


class ScreenConfig(Model):
    enabled: bool = False
    homopolymer_min: int = Field(default=8, ge=3)
    tandem_repeat_min_copies: int = Field(default=4, ge=3)
    tandem_repeat_max_unit: int = Field(default=6, ge=2, le=12)
    orf_min_codons: int = Field(default=30, ge=2)
    matrices: tuple[ScoringMatrix, ...] = ()
    motifs: tuple[MotifDefinition, ...] = ()
    #: Named restriction-site libraries screened as assembly constraints. Each entry
    #: reports the enzyme sites that standard forbids inside an insert.
    assembly_standards: tuple[
        Literal["golden-gate-bsai", "golden-gate-bsmbi", "biobrick"], ...
    ] = ()

    @model_validator(mode="after")
    def unique_standards(self) -> "ScreenConfig":
        if len(set(self.assembly_standards)) != len(self.assembly_standards):
            raise ValueError("assembly standards must be unique")
        return self


class ScreenHit(Model):
    kind: str
    name: str
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    strand: Literal["+", "-", "both"]
    score: float | None = None
    evidence: str


class CodonProfile(Model):
    name: str = Field(min_length=1)
    version: str = Field(min_length=1)
    source: str = Field(min_length=1)
    citation: str = Field(min_length=1)
    description: str = Field(min_length=1)
    translation_table: int = 1
    frequencies: dict[str, float]
    pair_scores: dict[str, float] = Field(default_factory=dict)

    @model_validator(mode="after")
    def check_frequencies(self):
        import itertools
        import math

        codons = {"".join(c) for c in itertools.product("ACGT", repeat=3)}
        if set(self.frequencies) != codons:
            raise ValueError("profile requires all 64 DNA codons")
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in self.frequencies.values()):
            raise ValueError("codon frequencies must be finite fractions between zero and one")
        table = CodonTable.unambiguous_dna_by_id.get(self.translation_table)
        if table is None or set(table.stop_codons) & set(table.forward_table):
            raise ValueError("unsupported profile translation table")
        for aa in set(table.forward_table.values()):
            total = sum(self.frequencies[c] for c, a in table.forward_table.items() if a == aa)
            if not 0.95 <= total <= 1.05:
                raise ValueError(f"relative synonymous frequencies for {aa} must sum to one")
        if any(
            len(pair) != 6 or set(pair) - set("ACGT") or not math.isfinite(score)
            for pair, score in self.pair_scores.items()
        ):
            raise ValueError("codon pair scores require six-base keys and finite values")
        return self


class CodonMetrics(Model):
    region: str
    cai: float | None
    gc3_percent: float
    distinct_codons: int
    rare_codon_positions: tuple[int, ...]
    mean_pair_score: float | None
    scored_pairs: int
    total_pairs: int


class OptimizeSettings(Model):
    solver: Literal["greedy", "global"] = "greedy"
    search_budget: int = Field(default=10000, gt=0, le=10000000)
    #: "strict_subset" requires new_CpG_positions ⊆ original_CpG_positions, so no
    #: substitution may relocate a CpG to a coordinate that did not already carry
    #: one. "count_only" relaxes this to monotonic total reduction, which permits
    #: trading several CpGs for fewer at new coordinates; the audit then lists the
    #: introduced positions explicitly.
    cpg_policy: Literal["strict_subset", "count_only"] = "strict_subset"
    protected_regions: tuple[ProtectedRegion, ...] = ()
    protected_feature_types: tuple[str, ...] = ("regulatory", "promoter", "repeat_region")
    gc3_min: float = Field(default=0, ge=0, le=100)
    gc3_max: float = Field(default=100, ge=0, le=100)
    min_distinct_codons: int = Field(default=0, ge=0, le=64)
    rare_codon_threshold: float = Field(default=0.1, ge=0, le=1)
    max_rare_codons: int | None = Field(default=None, ge=0)
    prevent_screen_regressions: bool = True

    @model_validator(mode="after")
    def gc_bounds(self):
        if self.gc3_min > self.gc3_max:
            raise ValueError("GC3 minimum exceeds maximum")
        return self


class SearchResult(Model):
    solver: Literal["greedy", "global"]
    status: Literal[
        "locally_optimal", "proven_optimal", "best_found", "infeasible", "no_feasible_found"
    ]
    explored_states: int = Field(ge=0)
    budget: int | None = None
    objective: str = "CpG count, changed bases, descending codon preference, sequence"
    proof_scope: str


class ComponentAlternative(Model):
    name: str = Field(min_length=1)
    length: int = Field(gt=0)
    source: str = Field(min_length=1)
    citation: str = Field(min_length=1)


class VectorComponent(ProtectedRegion):
    kind: Literal["itr", "cds", "promoter", "polya", "linker", "regulatory", "other"]
    protected: bool = True
    alternatives: tuple[ComponentAlternative, ...] = ()


class VectorArchitecture(Model):
    genome_type: Literal["ssAAV", "scAAV"] = "ssAAV"
    representation: Literal["itr_to_itr_unit", "cassette"]
    # Explicit lengths for omitted ITRs; never silently assume they were included.
    omitted_itr_lengths: tuple[int, int] | None = None
    limit: int | None = Field(default=None, gt=0)
    components: tuple[VectorComponent, ...] = ()

    @model_validator(mode="after")
    def check_components(self):
        if self.representation == "cassette":
            if self.omitted_itr_lengths is None or any(n <= 0 for n in self.omitted_itr_lengths):
                raise ValueError("cassette input requires two positive omitted_itr_lengths")
        elif self.omitted_itr_lengths is not None:
            raise ValueError("ITR-to-ITR input cannot also add omitted ITRs")
        ordered = sorted(self.components, key=lambda c: c.start)
        if len({c.name for c in ordered}) != len(ordered):
            raise ValueError("component names must be unique")
        if any(a.end >= b.start for a, b in zip(ordered, ordered[1:], strict=False)):
            raise ValueError("budget components must not overlap")
        return self


class LintConfig(Model):
    max_length: int = Field(default=4700, gt=0)
    gc_window: int = Field(default=100, gt=0)
    gc_step: int = Field(default=10, gt=0)
    gc_min_percent: float = Field(default=30.0, ge=0, le=100)
    gc_max_percent: float = Field(default=70.0, ge=0, le=100)
    cpg_window: int = Field(default=100, ge=2)
    cpg_step: int = Field(default=10, gt=0)
    cpg_hotspot_threshold: int = Field(default=6, gt=0)
    forbidden_motifs: tuple[MotifDefinition, ...] = ()
    splice_scan: bool = False
    rnafold: bool = False
    rnafold_mode: Literal["off", "auto", "seqfold", "required"] = "off"
    structure_window: int = Field(default=150, ge=20)
    structure_step: int = Field(default=25, gt=0)
    structure_mfe_threshold: float = -15.0
    structure_temperature: float = Field(default=37.0, ge=0, le=100, allow_inf_nan=False)
    structure_ensemble: bool = False
    screens: ScreenConfig = Field(default_factory=ScreenConfig)
    vector: VectorArchitecture | None = None

    @model_validator(mode="after")
    def check_thresholds(self) -> "LintConfig":
        if self.gc_min_percent >= self.gc_max_percent:
            raise ValueError("minimum GC percentage must be less than maximum GC percentage")
        if self.gc_step > self.gc_window:
            raise ValueError("GC step cannot exceed GC window size")
        if self.cpg_step > self.cpg_window:
            raise ValueError("CpG step cannot exceed CpG window size")
        if self.structure_step > self.structure_window:
            raise ValueError("structure step cannot exceed structure window size")
        names = [motif.name.casefold() for motif in self.forbidden_motifs]
        if len(names) != len(set(names)):
            raise ValueError("forbidden motif names must be unique")
        return self


class WindowMetric(Model):
    # User-facing coordinates are one-based and inclusive.
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    gc_percent: float | None = None
    cpg_count: int = Field(ge=0)
    ambiguous_bases: int = Field(ge=0)


class MotifMatch(Model):
    name: str
    sequence: str
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    strand: Literal["+", "-", "both"]


class SpliceSignal(Model):
    kind: Literal["donor", "acceptor"]
    sequence: str
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    strand: Literal["+", "-"]


class FoldingWindow(Model):
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    gc_percent: float
    mfe_kcal_per_mol: float
    dot_bracket: str
    flagged: bool
    ensemble_energy: float | None = None
    mfe_probability: float | None = None
    ensemble_diversity: float | None = None


class LintMetrics(Model):
    length: int = Field(gt=0)
    within_length_limit: bool
    length_overage: int = Field(ge=0)
    gc_percent: float | None
    cpg_count: int = Field(ge=0)
    cpg_positions: tuple[int, ...]
    cpg_observed_expected: float | None
    gc_windows: tuple[WindowMetric, ...]
    cpg_hotspots: tuple[WindowMetric, ...]
    motif_matches: tuple[MotifMatch, ...]
    splice_signals: tuple[SpliceSignal, ...] = ()
    folding_tool_version: str | None = None
    folding_windows: tuple[FoldingWindow, ...] = ()
    folding_windows_skipped: int = Field(default=0, ge=0)
    screen_hits: tuple[ScreenHit, ...] = ()
    vector_budget: dict | None = None


class LintReport(Model):
    schema_version: str = "1.0"
    tool_version: str
    identifier: str
    sequence_sha256: str
    configuration: LintConfig
    validation: ValidationReport
    metrics: LintMetrics
    findings: tuple[Finding, ...]

    @property
    def valid(self) -> bool:
        return self.validation.valid and not any(item.severity == "error" for item in self.findings)


class BaseChange(Model):
    position: int = Field(gt=0)
    before: str = Field(min_length=1, max_length=1)
    after: str = Field(min_length=1, max_length=1)

    @model_validator(mode="after")
    def check_bases(self) -> "BaseChange":
        if self.before not in "ACGT" or self.after not in "ACGT":
            raise ValueError("base changes must use A/C/G/T")
        if self.before == self.after:
            raise ValueError("base change must alter the nucleotide")
        return self


class Mutation(Model):
    index: int = Field(gt=0)
    region: str
    amino_acid: str = Field(min_length=1, max_length=1)
    start: int = Field(gt=0)
    end: int = Field(gt=0)
    strand: Literal[1, -1]
    codon_before: str = Field(min_length=3, max_length=3)
    codon_after: str = Field(min_length=3, max_length=3)
    host_codon_frequency_before: float
    host_codon_frequency_after: float
    cpg_count_before: int = Field(ge=0)
    cpg_count_after: int = Field(ge=0)
    reason: str
    base_changes: tuple[BaseChange, ...]

    @model_validator(mode="after")
    def check_mutation(self) -> "Mutation":
        if set(self.codon_before + self.codon_after) - set("ACGT"):
            raise ValueError("mutation codons must contain only A/C/G/T")
        if self.codon_before == self.codon_after:
            raise ValueError("mutation codons must differ")
        if not self.base_changes:
            raise ValueError("mutation must contain at least one base change")
        positions = [change.position for change in self.base_changes]
        if len(positions) != len(set(positions)):
            raise ValueError("mutation base-change positions must be unique")
        return self


class RuntimeMetadata(Model):
    python: str
    platform: str
    biopython: str
    pydantic: str
    seqfold: str
    typer: str


class OptimizeReport(Model):
    schema_version: Literal["2.0"] = "2.0"
    tool_version: str
    identifier: str
    source_sha256: str
    optimized_sha256: str
    host: str
    host_profile: str
    max_length: int = Field(gt=0)
    cpg_depletion_requested: bool
    protein_translation_preserved: bool
    cpg_count_before: int = Field(ge=0)
    cpg_count_after: int = Field(ge=0)
    cpg_reduction: int = Field(ge=0)
    coding_regions: tuple[CodingRegion, ...]
    mutations: tuple[Mutation, ...]
    residual_cpg_positions: tuple[int, ...]
    #: Positions carrying a CpG after optimization that carried none before. Always
    #: empty under the default "strict_subset" policy; populated and reported as a
    #: warning under "count_only".
    introduced_cpg_positions: tuple[int, ...] = ()
    before_lint: LintReport
    after_lint: LintReport
    safety_checks_passed: bool
    safety_regressions: tuple[str, ...] = ()
    runtime: RuntimeMetadata
    warnings: tuple[str, ...] = ()
    edit_semantics: Literal["sequential", "joint"] = "sequential"
    settings: OptimizeSettings = Field(default_factory=OptimizeSettings)
    profile: CodonProfile
    profile_sha256: str
    search: SearchResult
    codon_metrics_before: tuple[CodonMetrics, ...] = ()
    codon_metrics_after: tuple[CodonMetrics, ...] = ()

    @model_validator(mode="after")
    def check_consistency(self) -> "OptimizeReport":
        sha256 = re.compile(r"[0-9a-f]{64}")
        if not sha256.fullmatch(self.source_sha256) or not sha256.fullmatch(
            self.optimized_sha256
        ):
            raise ValueError("sequence SHA-256 values must be 64 lowercase hexadecimal characters")
        if self.cpg_count_after > self.cpg_count_before:
            raise ValueError("optimized CpG count cannot exceed the source count")
        introduced = tuple(
            sorted(
                set(self.after_lint.metrics.cpg_positions)
                - set(self.before_lint.metrics.cpg_positions)
            )
        )
        if self.settings.cpg_policy == "strict_subset" and introduced:
            raise ValueError("optimized sequence introduces a CpG at a new position")
        if introduced != self.introduced_cpg_positions:
            raise ValueError("introduced CpG positions must match the before/after lint reports")
        if self.cpg_reduction != self.cpg_count_before - self.cpg_count_after:
            raise ValueError("CpG reduction must equal before minus after")
        if len(self.residual_cpg_positions) != self.cpg_count_after:
            raise ValueError("residual CpG positions must match the optimized CpG count")
        if tuple(sorted(set(self.residual_cpg_positions))) != self.residual_cpg_positions:
            raise ValueError("residual CpG positions must be sorted and unique")
        if self.identifier not in (self.before_lint.identifier, self.after_lint.identifier):
            raise ValueError("optimization and lint identifiers do not match")
        if self.before_lint.identifier != self.after_lint.identifier:
            raise ValueError("before/after lint identifiers do not match")
        if self.before_lint.configuration != self.after_lint.configuration:
            raise ValueError("before/after lint configurations do not match")
        if self.before_lint.sequence_sha256 != self.source_sha256:
            raise ValueError("source SHA-256 does not match the before-lint report")
        if self.after_lint.sequence_sha256 != self.optimized_sha256:
            raise ValueError("optimized SHA-256 does not match the after-lint report")
        if self.before_lint.metrics.cpg_count != self.cpg_count_before:
            raise ValueError("source CpG count does not match the before-lint report")
        if self.after_lint.metrics.cpg_count != self.cpg_count_after:
            raise ValueError("optimized CpG count does not match the after-lint report")
        if self.max_length != self.after_lint.configuration.max_length:
            raise ValueError("maximum length does not match the lint configuration")
        if self.safety_checks_passed == bool(self.safety_regressions):
            raise ValueError("safety status and regression list disagree")
        expected_indexes = tuple(range(1, len(self.mutations) + 1))
        if tuple(mutation.index for mutation in self.mutations) != expected_indexes:
            raise ValueError("mutation indexes must be contiguous and one-based")
        expected_cpg = self.cpg_count_before
        for mutation in self.mutations:
            if self.edit_semantics == "joint":
                if (mutation.cpg_count_before, mutation.cpg_count_after) != (
                    self.cpg_count_before, self.cpg_count_after
                ):
                    raise ValueError("joint edits must record whole-transaction CpG counts")
                expected_cpg = self.cpg_count_after
                continue
            if mutation.cpg_count_after >= mutation.cpg_count_before:
                raise ValueError("each sequential mutation must reduce the whole-payload CpG count")
            if mutation.cpg_count_before != expected_cpg:
                raise ValueError("mutation CpG counts do not form a continuous audit trail")
            expected_cpg = mutation.cpg_count_after
        if expected_cpg != self.cpg_count_after:
            raise ValueError("mutation audit trail does not end at the reported CpG count")
        return self

    @property
    def optimized(self) -> bool:
        return (
            self.cpg_depletion_requested
            and self.protein_translation_preserved
            and self.safety_checks_passed
            and self.after_lint.valid
            and self.search.status not in ("infeasible", "no_feasible_found")
        )
