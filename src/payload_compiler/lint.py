"""Deterministic, read-only payload metrics and lint rules."""

import re

from Bio.Seq import Seq

from payload_compiler import __version__
from payload_compiler.analysis import (
    fold_local_windows,
    fold_local_windows_seqfold,
    rnafold_available,
    scan_splice_signals,
)
from payload_compiler.design import vector_budget
from payload_compiler.exceptions import InputError
from payload_compiler.models import (
    Finding,
    LintConfig,
    LintMetrics,
    LintReport,
    MotifDefinition,
    MotifMatch,
    Payload,
    WindowMetric,
)
from payload_compiler.screening import screen_sequence
from payload_compiler.validation import DNA_BASES, validate_payload

IUPAC_PATTERN = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def parse_motif(spec: str) -> MotifDefinition:
    """Parse NAME=IUPAC_SEQUENCE from the CLI."""
    if "=" not in spec:
        raise InputError(f"Invalid motif '{spec}'; expected NAME=IUPAC_SEQUENCE")
    name, sequence = (part.strip() for part in spec.split("=", 1))
    sequence = sequence.upper()
    if not name or not sequence:
        raise InputError(f"Invalid motif '{spec}'; name and sequence are required")
    invalid = sorted(set(sequence) - IUPAC_PATTERN.keys())
    if invalid:
        raise InputError(f"Invalid motif '{name}'; unsupported DNA symbols: {', '.join(invalid)}")
    return MotifDefinition(name=name, sequence=sequence)


def _window_ranges(length: int, size: int, step: int) -> tuple[tuple[int, int], ...]:
    """Return zero-based half-open windows, including a final end-anchored window."""
    if length <= size:
        return ((0, length),)
    last = length - size
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    return tuple((start, start + size) for start in starts)


def _gc_percent(sequence: str) -> tuple[float | None, int]:
    resolved = sum(base in DNA_BASES for base in sequence)
    ambiguous = len(sequence) - resolved
    if not resolved:
        return None, ambiguous
    value = 100 * (sequence.count("G") + sequence.count("C")) / resolved
    return round(value, 4), ambiguous


def _cpg_positions(sequence: str, *, offset: int = 0) -> tuple[int, ...]:
    return tuple(offset + match.start() + 1 for match in re.finditer(r"(?=CG)", sequence))


def _windows(sequence: str, *, size: int, step: int) -> tuple[WindowMetric, ...]:
    windows = []
    for start, end in _window_ranges(len(sequence), size, step):
        fragment = sequence[start:end]
        gc_percent, ambiguous = _gc_percent(fragment)
        windows.append(
            WindowMetric(
                start=start + 1,
                end=end,
                gc_percent=gc_percent,
                cpg_count=len(_cpg_positions(fragment)),
                ambiguous_bases=ambiguous,
            )
        )
    return tuple(windows)


def _motif_regex(sequence: str) -> re.Pattern[str]:
    expanded = "".join(IUPAC_PATTERN[base] for base in sequence)
    return re.compile(f"(?=({expanded}))")


def _motif_matches(sequence: str, motifs: tuple[MotifDefinition, ...]) -> tuple[MotifMatch, ...]:
    matches: list[MotifMatch] = []
    for motif in motifs:
        reverse = str(Seq(motif.sequence).reverse_complement())
        orientations = (
            ((motif.sequence, "both"),)
            if reverse == motif.sequence
            else (
                (motif.sequence, "+"),
                (reverse, "-"),
            )
        )
        for pattern, strand in orientations:
            for match in _motif_regex(pattern).finditer(sequence):
                start = match.start() + 1
                matches.append(
                    MotifMatch(
                        name=motif.name,
                        sequence=motif.sequence,
                        start=start,
                        end=start + len(motif.sequence) - 1,
                        strand=strand,
                    )
                )
    return tuple(sorted(matches, key=lambda item: (item.start, item.end, item.name, item.strand)))


def lint_payload(payload: Payload, config: LintConfig | None = None) -> LintReport:
    config = config or LintConfig()
    validation = validate_payload(payload)
    sequence = payload.sequence
    findings: list[Finding] = []
    budget = vector_budget(payload, config.vector) if config.vector else None
    if budget and not budget["within_limit"]:
        findings.append(Finding(
            code="VECTOR_BUDGET_EXCEEDED", severity="error",
            message=f"Accounted vector unit exceeds its {budget['limit']}-nt budget "
                    f"by {-budget['headroom']} nt (including declared ITR overhead)",
        ))
    screen_hits = screen_sequence(payload, config.screens)
    for kind in sorted({h.kind for h in screen_hits}):
        selected = [h for h in screen_hits if h.kind == kind]
        findings.append(Finding(
            code=f"SCREEN_{kind.upper()}", severity="warning",
            message=(f"{len(selected)} {kind} screening hit(s); heuristic, "
                     "not biological validation"),
            positions=tuple(h.start for h in selected),
        ))

    if len(sequence) > config.max_length:
        findings.append(
            Finding(
                code="LENGTH_LIMIT_EXCEEDED",
                severity="error",
                message=(
                    f"Payload is {len(sequence) - config.max_length} nt over the "
                    f"{config.max_length} nt limit"
                ),
            )
        )

    gc_windows = _windows(sequence, size=config.gc_window, step=config.gc_step)
    low_gc = tuple(
        window.start
        for window in gc_windows
        if window.gc_percent is not None and window.gc_percent < config.gc_min_percent
    )
    high_gc = tuple(
        window.start
        for window in gc_windows
        if window.gc_percent is not None and window.gc_percent > config.gc_max_percent
    )
    if low_gc:
        findings.append(
            Finding(
                code="LOW_GC_WINDOW",
                severity="warning",
                message=f"{len(low_gc)} GC window(s) fall below {config.gc_min_percent:g}%",
                positions=low_gc,
            )
        )
    if high_gc:
        findings.append(
            Finding(
                code="HIGH_GC_WINDOW",
                severity="warning",
                message=f"{len(high_gc)} GC window(s) exceed {config.gc_max_percent:g}%",
                positions=high_gc,
            )
        )

    cpg_windows = _windows(sequence, size=config.cpg_window, step=config.cpg_step)
    hotspots = tuple(
        window for window in cpg_windows if window.cpg_count >= config.cpg_hotspot_threshold
    )
    if hotspots:
        findings.append(
            Finding(
                code="CPG_HOTSPOT",
                severity="warning",
                message=(
                    f"{len(hotspots)} window(s) contain at least "
                    f"{config.cpg_hotspot_threshold} CpG sites"
                ),
                positions=tuple(window.start for window in hotspots),
            )
        )

    motif_matches = _motif_matches(sequence, config.forbidden_motifs)
    for motif in config.forbidden_motifs:
        positions = tuple(match.start for match in motif_matches if match.name == motif.name)
        if positions:
            findings.append(
                Finding(
                    code="FORBIDDEN_MOTIF",
                    severity="error",
                    message=f"Forbidden motif '{motif.name}' occurs {len(positions)} time(s)",
                    region=motif.name,
                    positions=positions,
                )
            )

    splice_signals = scan_splice_signals(sequence) if config.splice_scan else ()
    if splice_signals:
        donor_count = sum(signal.kind == "donor" for signal in splice_signals)
        acceptor_count = sum(signal.kind == "acceptor" for signal in splice_signals)
        findings.append(
            Finding(
                code="SPLICE_CONSENSUS",
                severity="warning",
                message=(
                    f"Found {donor_count} donor- and {acceptor_count} acceptor-consensus "
                    "site(s); motif hits are not splice predictions"
                ),
                positions=tuple(signal.start for signal in splice_signals),
            )
        )

    folding_tool_version = None
    folding_windows = ()
    folding_windows_skipped = 0
    rnafold_requested = config.rnafold or config.rnafold_mode != "off"
    if rnafold_requested:
        use_native_rnafold = config.rnafold_mode == "required" or (
            (config.rnafold_mode == "auto" or config.rnafold) and rnafold_available()
        )
        folding_function = fold_local_windows if use_native_rnafold else fold_local_windows_seqfold
        folding_tool_version, folding_windows, folding_windows_skipped = folding_function(
            sequence,
            window_size=config.structure_window,
            step=config.structure_step,
            mfe_threshold=config.structure_mfe_threshold,
            temperature=config.structure_temperature,
            ensemble=config.structure_ensemble,
        )
        flagged_folds = tuple(fold for fold in folding_windows if fold.flagged)
        if flagged_folds:
            findings.append(
                Finding(
                    code="STRONG_LOCAL_FOLD",
                    severity="warning",
                    message=(
                        f"{len(flagged_folds)} local window(s) have folding MFE at or below "
                        f"{config.structure_mfe_threshold:g} kcal/mol"
                    ),
                    positions=tuple(fold.start for fold in flagged_folds),
                )
            )
        if folding_windows_skipped:
            findings.append(
                Finding(
                    code="RNAFOLD_AMBIGUOUS_WINDOWS_SKIPPED",
                    severity="warning",
                    message=(
                        f"Skipped {folding_windows_skipped} local window(s) containing "
                        "ambiguous bases"
                    ),
                )
            )

    gc_percent, _ = _gc_percent(sequence)
    cpg_positions = _cpg_positions(sequence)
    c_count, g_count = sequence.count("C"), sequence.count("G")
    resolved_length = sum(base in DNA_BASES for base in sequence)
    cpg_oe = None
    if c_count and g_count and resolved_length:
        cpg_oe = round(len(cpg_positions) * resolved_length / (c_count * g_count), 6)

    return LintReport(
        tool_version=__version__,
        identifier=payload.identifier,
        sequence_sha256=validation.sequence_sha256,
        configuration=config,
        validation=validation,
        metrics=LintMetrics(
            length=len(sequence),
            within_length_limit=len(sequence) <= config.max_length,
            length_overage=max(0, len(sequence) - config.max_length),
            gc_percent=gc_percent,
            cpg_count=len(cpg_positions),
            cpg_positions=cpg_positions,
            cpg_observed_expected=cpg_oe,
            gc_windows=gc_windows,
            cpg_hotspots=hotspots,
            motif_matches=motif_matches,
            splice_signals=splice_signals,
            folding_tool_version=folding_tool_version,
            folding_windows=folding_windows,
            folding_windows_skipped=folding_windows_skipped,
            screen_hits=screen_hits,
            vector_budget=budget,
        ),
        findings=tuple(findings),
    )
