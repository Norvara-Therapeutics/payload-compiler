"""Optional splice-consensus and ViennaRNA secondary-structure analyses."""

import re
import shutil
import subprocess
from importlib.metadata import version

from Bio.Seq import Seq
from seqfold import dg as seqfold_dg
from seqfold import dot_bracket as seqfold_dot_bracket
from seqfold import fold as seqfold_fold

from payload_compiler.exceptions import InputError
from payload_compiler.models import FoldingWindow, SpliceSignal

_DONOR = re.compile(r"(?=(GT[AG]AGT))")
_ACCEPTOR = re.compile(r"(?=([CT]{5,12}AGG))")
_FOLDING_LINE = re.compile(r"^([().]+)\s+\(\s*(-?\d+(?:\.\d+)?)\s*\)$")


def rnafold_available() -> bool:
    """Whether the ViennaRNA folding executable is on PATH."""
    return shutil.which("RNAfold") is not None


def scan_splice_signals(sequence: str) -> tuple[SpliceSignal, ...]:
    """Find canonical splice-consensus signals on either strand.

    This is a motif screen only; it does not predict spliceosome recognition.
    Donors use a common six-base U1 consensus and acceptors require a short
    polypyrimidine tract followed by the canonical AG|G boundary.
    """
    signals: list[SpliceSignal] = []
    length = len(sequence)
    for strand, oriented in (("+", sequence), ("-", str(Seq(sequence).reverse_complement()))):
        for kind, pattern in (("donor", _DONOR), ("acceptor", _ACCEPTOR)):
            for match in pattern.finditer(oriented):
                motif = match.group(1)
                offset = match.start(1)
                if strand == "+":
                    start, end = offset + 1, offset + len(motif)
                else:
                    start = length - offset - len(motif) + 1
                    end = length - offset
                signals.append(
                    SpliceSignal(
                        kind=kind,
                        sequence=motif,
                        start=start,
                        end=end,
                        strand=strand,
                    )
                )
    return tuple(sorted(signals, key=lambda item: (item.start, item.end, item.kind, item.strand)))


def _windows(length: int, size: int, step: int) -> tuple[tuple[int, int], ...]:
    if length <= size:
        return ((0, length),)
    last = length - size
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    return tuple((start, start + size) for start in starts)


def _gc_percent(sequence: str) -> float:
    return round(100 * sum(base in "GC" for base in sequence) / len(sequence), 4)


def fold_local_windows(
    sequence: str,
    *,
    window_size: int,
    step: int,
    mfe_threshold: float,
    timeout_seconds: float = 60.0,
    temperature: float = 37.0,
    ensemble: bool = False,
) -> tuple[str, tuple[FoldingWindow, ...], int]:
    """Fold resolved local windows with RNAfold, returning version, folds, skipped count."""
    executable = shutil.which("RNAfold")
    if executable is None:
        raise InputError(
            "RNAfold analysis was requested but the ViennaRNA RNAfold executable was not found"
        )

    window_ranges = _windows(len(sequence), window_size, step)
    inputs: list[tuple[int, int, str]] = []
    skipped = 0
    for start, end in window_ranges:
        fragment = sequence[start:end]
        if set(fragment) - set("ACGT"):
            skipped += 1
            continue
        inputs.append((start, end, fragment.replace("T", "U")))
    if not inputs:
        raise InputError("RNAfold analysis found no windows containing only A, C, G, and T")

    try:
        version_result = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=timeout_seconds
        )
        version_result.check_returncode()
        version = (version_result.stdout or version_result.stderr).strip().splitlines()[0]
        fasta = "".join(
            f">window_{index}\n{fragment}\n"
            for index, (_, _, fragment) in enumerate(inputs, start=1)
        )
        command = [executable, "--noPS", f"--temp={temperature:g}"]
        if ensemble:
            command.extend(["-p", "--noDP"])
        result = subprocess.run(
            command,
            input=fasta,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
        result.check_returncode()
    except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        raise InputError(f"RNAfold analysis failed: {exc}") from exc

    parsed: list[tuple[str, float]] = []
    ensemble_data: list[tuple[float, float, float]] = []
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    index = 0
    while index < len(lines):
        if lines[index].startswith(">"):
            if index + 2 >= len(lines):
                break
            structure_match = _FOLDING_LINE.fullmatch(lines[index + 2])
            if structure_match is None:
                raise InputError("Could not parse RNAfold output structure and energy")
            parsed.append((structure_match.group(1), float(structure_match.group(2))))
            if ensemble:
                stop = next((j for j in range(index + 1, len(lines))
                             if lines[j].startswith(">")), len(lines))
                block = "\n".join(lines[index + 3:stop])
                energy = re.search(r"\[\s*(-?\d+(?:\.\d+)?)\s*\]", block)
                stats = re.search(
                    r"frequency of mfe structure in ensemble\s+([\d.eE+\-]+);\s*"
                    r"ensemble diversity\s+([\d.eE+\-]+)", block
                )
                if energy is None or stats is None:
                    raise InputError("RNAfold did not return requested ensemble statistics")
                ensemble_data.append((float(energy[1]), float(stats[1]), float(stats[2])))
            index += 3
        else:
            index += 1
    if len(parsed) != len(inputs):
        raise InputError(
            f"RNAfold returned {len(parsed)} folds for {len(inputs)} requested windows"
        )

    if any(len(structure) != end - start
           for (start, end, _), (structure, _) in zip(inputs, parsed, strict=True)):
        raise InputError("RNAfold structure length does not match its input window")
    folds = tuple(
        FoldingWindow(
            start=start + 1,
            end=end,
            gc_percent=_gc_percent(sequence[start:end]),
            mfe_kcal_per_mol=mfe,
            dot_bracket=structure,
            flagged=mfe <= mfe_threshold,
            ensemble_energy=ensemble_data[i][0] if ensemble else None,
            mfe_probability=ensemble_data[i][1] if ensemble else None,
            ensemble_diversity=ensemble_data[i][2] if ensemble else None,
        )
        for i, ((start, end, _), (structure, mfe)) in enumerate(zip(inputs, parsed, strict=True))
    )
    return version, folds, skipped


def fold_local_windows_seqfold(
    sequence: str,
    *,
    window_size: int,
    step: int,
    mfe_threshold: float,
    temperature: float = 37.0,
    ensemble: bool = False,
) -> tuple[str, tuple[FoldingWindow, ...], int]:
    """Fold local RNA windows with the pip-installed seqfold backend."""
    if ensemble:
        raise InputError("Ensemble statistics require native RNAfold; seqfold provides MFE only")
    inputs: list[tuple[int, int, str]] = []
    skipped = 0
    for start, end in _windows(len(sequence), window_size, step):
        fragment = sequence[start:end]
        if set(fragment) - set("ACGT"):
            skipped += 1
            continue
        inputs.append((start, end, fragment.replace("T", "U")))
    if not inputs:
        raise InputError("Secondary-structure analysis found no unambiguous DNA windows")

    folds = []
    try:
        for start, end, fragment in inputs:
            structures = seqfold_fold(fragment, temp=temperature)
            mfe = float(seqfold_dg(fragment, temp=temperature))
            # seqfold uses a large positive sentinel when no favorable fold
            # exists. The physically relevant reference is the unfolded state
            # at 0 kcal/mol.
            dot_bracket = (
                seqfold_dot_bracket(fragment, structures) if mfe < 0 else "." * len(fragment)
            )
            mfe = min(0.0, mfe)
            folds.append(
                FoldingWindow(
                    start=start + 1,
                    end=end,
                    gc_percent=_gc_percent(sequence[start:end]),
                    mfe_kcal_per_mol=mfe,
                    dot_bracket=dot_bracket,
                    flagged=mfe <= mfe_threshold,
                )
            )
    except (RuntimeError, ValueError) as exc:
        raise InputError(f"seqfold analysis failed: {exc}") from exc
    return f"seqfold {version('seqfold')}", tuple(folds), skipped
