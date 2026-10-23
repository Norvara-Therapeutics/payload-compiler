"""Command-line entry points."""

import json
import os
import tempfile
from enum import StrEnum
from io import StringIO
from pathlib import Path
from typing import Annotated, Literal

import typer
from Bio import SeqIO
from Bio.Seq import Seq
from pydantic import ValidationError

from payload_compiler import __version__
from payload_compiler.audit import verify_optimization_audit
from payload_compiler.design import read_json_model
from payload_compiler.exceptions import InputError
from payload_compiler.lint import lint_payload, parse_motif
from payload_compiler.models import (
    CodonProfile,
    LintConfig,
    LintReport,
    OptimizeReport,
    OptimizeSettings,
    ProtectedRegion,
    ScoringMatrix,
    ScreenConfig,
    ValidationReport,
    VectorArchitecture,
)
from payload_compiler.optimize import optimize_payload
from payload_compiler.parser import SUFFIXES, read_payload
from payload_compiler.profiles import HOST_PROFILE_VERSIONS, HOST_PROFILES
from payload_compiler.reporting import render_lint_html, render_optimize_html
from payload_compiler.validation import validate_payload

app = typer.Typer(
    no_args_is_help=True,
    invoke_without_command=True,
    help="Optimize, validate, lint, report on DNA payloads, and inspect host profiles.",
)


class Format(StrEnum):
    fasta = "fasta"
    genbank = "genbank"
    embl = "embl"


class Solver(StrEnum):
    greedy = "greedy"
    global_ = "global"


class CpgPolicy(StrEnum):
    strict_subset = "strict-subset"
    count_only = "count-only"


class AssemblyStandard(StrEnum):
    golden_gate_bsai = "golden-gate-bsai"
    golden_gate_bsmbi = "golden-gate-bsmbi"
    biobrick = "biobrick"


@app.command()
def profiles(
    json_output: Annotated[bool, typer.Option("--json", help="Write profile data as JSON")] = False,
):
    """List bundled host codon-frequency profiles."""
    if json_output:
        data = {
            name: {"version": HOST_PROFILE_VERSIONS[name], "codons": HOST_PROFILES[name]}
            for name in sorted(HOST_PROFILES)
        }
        typer.echo(json.dumps(data, sort_keys=True, indent=2))
        return
    for name in sorted(HOST_PROFILES):
        typer.echo(f"{name}: {HOST_PROFILE_VERSIONS[name]}")


@app.callback()
def main(version: Annotated[bool, typer.Option("--version", is_eager=True)] = False):
    if version:
        typer.echo(f"payload-compiler {__version__}")
        raise typer.Exit()


def _json_report(report: ValidationReport | LintReport | OptimizeReport) -> str:
    data = report.model_dump(mode="json")
    if isinstance(report, OptimizeReport):
        data["optimized"] = report.optimized
    else:
        data["valid"] = report.valid
    return json.dumps(data, sort_keys=True, indent=2)


def _write_outputs_atomically(
    outputs: dict[Path, str], *, input_path: Path, remove: tuple[Path, ...] = ()
) -> None:
    """Stage a related output set and roll every path back if a commit fails."""
    targets = tuple(outputs) + remove
    _validate_output_paths(input_path, *targets)
    staged: dict[Path, Path] = {}
    backups: dict[Path, Path] = {}
    committed: list[Path] = []
    try:
        for path, content in outputs.items():
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="\n",
                dir=path.parent,
                prefix=f".{path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary:
                temporary.write(content)
                staged[path] = Path(temporary.name)
        for path in targets:
            if path.exists():
                with tempfile.NamedTemporaryFile(
                    dir=path.parent,
                    prefix=f".{path.name}.",
                    suffix=".backup",
                    delete=False,
                ) as placeholder:
                    backup = Path(placeholder.name)
                backup.unlink()
                os.replace(path, backup)
                backups[path] = backup
        for path, temporary in staged.items():
            os.replace(temporary, path)
            committed.append(path)
    except OSError as exc:
        for path in reversed(committed):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        for path, backup in backups.items():
            try:
                os.replace(backup, path)
            except OSError:
                pass
        raise InputError(f"Cannot commit output set: {exc}") from exc
    finally:
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)
        for backup in backups.values():
            backup.unlink(missing_ok=True)


def _write_output(path: Path, content: str, *, input_path: Path) -> None:
    _write_outputs_atomically({path: content}, input_path=input_path)


def _sequence_output_format(output: Path, explicit: Format | None) -> Format:
    if explicit is not None:
        return explicit
    inferred = SUFFIXES.get(output.suffix.lower())
    return Format(inferred or "fasta")


def _render_sequence_output(
    input_path: Path, payload, output_format: Format
) -> str:
    if output_format == Format.fasta:
        description = payload.description
        header = (
            description
            if description == payload.identifier
            or description.startswith(f"{payload.identifier} ")
            else f"{payload.identifier} {description}".rstrip()
        )
        return f">{header}\n" + "\n".join(
            payload.sequence[offset : offset + 60]
            for offset in range(0, len(payload.sequence), 60)
        ) + "\n"
    if payload.source_format != output_format.value:
        raise InputError(
            f"{output_format.value} output requires a matching annotated input; "
            "use FASTA output for format conversion"
        )
    try:
        with input_path.open(encoding="utf-8-sig") as handle:
            record = SeqIO.read(handle, output_format.value)
        if len(record.seq) != len(payload.sequence):
            raise InputError("Annotated output requires an unchanged sequence length")
        record.seq = Seq(payload.sequence)
        rendered = StringIO()
        SeqIO.write(record, rendered, output_format.value)
        return rendered.getvalue()
    except InputError:
        raise
    except (OSError, ValueError) as exc:
        raise InputError(f"Cannot render annotated {output_format.value} output: {exc}") from exc


def _validate_output_paths(input_path: Path, *outputs: Path | None) -> None:
    selected = [path for path in outputs if path is not None]
    resolved = [path.resolve() for path in selected]
    if input_path.resolve() in resolved:
        raise InputError("An output file cannot overwrite the input sequence")
    if len(resolved) != len(set(resolved)):
        raise InputError("Output files must use distinct paths")
    for path in selected:
        if not path.parent.is_dir():
            raise InputError(f"Report directory does not exist: {path.parent}")


def _emit_input_error(exc: InputError, *, json_output: bool) -> None:
    if json_output:
        typer.echo(
            json.dumps(
                {
                    "schema_version": "1.0",
                    "valid": False,
                    "error": {"code": "INPUT_ERROR", "message": str(exc)},
                },
                sort_keys=True,
                indent=2,
            )
        )
    else:
        typer.echo(f"Input error: {exc}", err=True)


def _screen_config(
    *,
    enable: bool,
    homopolymer_min: int,
    repeat_min_copies: int,
    repeat_max_unit: int,
    orf_min_codons: int,
    matrix_paths: tuple[Path, ...],
    screen_motifs: tuple[str, ...],
    assembly_standards: tuple[str, ...],
) -> ScreenConfig:
    """Build the heuristic screen configuration from CLI inputs."""
    matrices = tuple(read_json_model(path, ScoringMatrix) for path in matrix_paths)
    return ScreenConfig(
        enabled=enable,
        homopolymer_min=homopolymer_min,
        tandem_repeat_min_copies=repeat_min_copies,
        tandem_repeat_max_unit=repeat_max_unit,
        orf_min_codons=orf_min_codons,
        matrices=matrices,
        motifs=tuple(parse_motif(spec) for spec in screen_motifs),
        assembly_standards=tuple(assembly_standards),
    )


def _vector_architecture(path: Path | None) -> VectorArchitecture | None:
    return None if path is None else read_json_model(path, VectorArchitecture)


def _parse_protected(spec: str) -> ProtectedRegion:
    """Parse START:END[:NAME] one-based inclusive immutable intervals."""
    parts = spec.split(":")
    if len(parts) not in (2, 3):
        raise InputError(f"Invalid --protect '{spec}'; expected START:END[:NAME]")
    try:
        start, end = int(parts[0]), int(parts[1])
    except ValueError as exc:
        raise InputError(f"Invalid --protect '{spec}'; coordinates must be integers") from exc
    name = parts[2] if len(parts) == 3 and parts[2] else "protected"
    try:
        return ProtectedRegion(start=start, end=end, name=name)
    except ValidationError as exc:
        message = "; ".join(error["msg"] for error in exc.errors())
        raise InputError(f"Invalid --protect '{spec}': {message}") from exc


@app.command()
def validate(
    input_path: Annotated[
        Path, typer.Option("--input", "-i", help="Single-record FASTA, GenBank, or EMBL")
    ],
    cds: Annotated[
        list[str] | None,
        typer.Option(
            "--cds",
            help="FASTA CDS: START:END[:+|-], one-based inclusive; repeat for multiple CDSs",
        ),
    ] = None,
    input_format: Annotated[Format | None, typer.Option("--format")] = None,
    require_cds: Annotated[bool, typer.Option("--require-cds")] = False,
    json_output: Annotated[
        bool, typer.Option("--json", help="Write deterministic JSON to stdout")
    ] = False,
):
    """Check DNA alphabet, CDS coordinates, reading frames, and translations."""
    try:
        payload = read_payload(
            input_path,
            input_format=input_format.value if input_format else None,
            cds=tuple(cds or ()),
        )
        report = validate_payload(payload, require_cds=require_cds)
    except InputError as exc:
        _emit_input_error(exc, json_output=json_output)
        raise typer.Exit(2) from exc
    if json_output:
        typer.echo(_json_report(report))
    else:
        typer.echo(
            f"{'PASS' if report.valid else 'FAIL'} {report.identifier}: "
            f"{report.length} nt, {report.cds_count} CDS"
        )
        for finding in report.findings:
            context = f" [{finding.region}]" if finding.region else ""
            typer.echo(f"{finding.severity.upper()} {finding.code}{context}: {finding.message}")
    raise typer.Exit(0 if report.valid else 1)


@app.command()
def lint(
    input_path: Annotated[
        Path, typer.Option("--input", "-i", help="Single-record FASTA, GenBank, or EMBL")
    ],
    cds: Annotated[
        list[str] | None,
        typer.Option(
            "--cds",
            help="FASTA CDS: START:END[:+|-], one-based inclusive; repeat for multiple CDSs",
        ),
    ] = None,
    input_format: Annotated[Format | None, typer.Option("--format")] = None,
    max_length: Annotated[int, typer.Option("--max-length")] = 4700,
    gc_window: Annotated[int, typer.Option("--gc-window")] = 100,
    gc_step: Annotated[int, typer.Option("--gc-step")] = 10,
    gc_min: Annotated[float, typer.Option("--gc-min")] = 30.0,
    gc_max: Annotated[float, typer.Option("--gc-max")] = 70.0,
    cpg_window: Annotated[int, typer.Option("--cpg-window")] = 100,
    cpg_step: Annotated[int, typer.Option("--cpg-step")] = 10,
    cpg_hotspot: Annotated[int, typer.Option("--cpg-hotspot")] = 6,
    splice_scan: Annotated[
        bool, typer.Option("--splice-scan", help="Scan both strands for splice consensus motifs")
    ] = False,
    rnafold: Annotated[
        bool,
        typer.Option(
            "--rnafold",
            help="Fold local windows; uses native RNAfold when available, otherwise seqfold",
        ),
    ] = False,
    structure_window: Annotated[int, typer.Option("--structure-window")] = 150,
    structure_step: Annotated[int, typer.Option("--structure-step")] = 25,
    structure_mfe: Annotated[
        float, typer.Option("--structure-mfe", help="Flag folds at or below this kcal/mol value")
    ] = -15.0,
    motif: Annotated[
        list[str] | None,
        typer.Option(
            "--forbid-motif",
            help="Forbidden motif NAME=IUPAC_SEQUENCE; scans both strands; repeat as needed",
        ),
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Write deterministic JSON to stdout")
    ] = False,
    json_file: Annotated[
        Path | None, typer.Option("--json-output", help="Write deterministic JSON to a file")
    ] = None,
    html_report: Annotated[
        Path | None, typer.Option("--report", help="Write a standalone HTML report")
    ] = None,
    vector_path: Annotated[
        Path | None,
        typer.Option("--vector", help="VectorArchitecture JSON for ITR-aware size accounting"),
    ] = None,
    screen: Annotated[
        bool, typer.Option("--screen/--no-screen", help="Run heuristic sequence screens")
    ] = False,
    homopolymer_min: Annotated[int, typer.Option("--homopolymer-min")] = 8,
    repeat_min_copies: Annotated[int, typer.Option("--repeat-min-copies")] = 4,
    repeat_max_unit: Annotated[int, typer.Option("--repeat-max-unit")] = 6,
    orf_min_codons: Annotated[int, typer.Option("--orf-min-codons")] = 30,
    matrix: Annotated[
        list[Path] | None,
        typer.Option("--matrix", help="ScoringMatrix JSON file; repeat as needed"),
    ] = None,
    screen_motif: Annotated[
        list[str] | None,
        typer.Option("--screen-motif", help="Report NAME=IUPAC_SEQUENCE as a screen hit"),
    ] = None,
    assembly_standard: Annotated[
        list[AssemblyStandard] | None,
        typer.Option("--assembly-standard", help="Screen a restriction-site library"),
    ] = None,
    structure_ensemble: Annotated[
        bool,
        typer.Option("--structure-ensemble", help="Report ensemble statistics, not only the MFE"),
    ] = False,
):
    """Measure sequence features and optional splice/secondary-structure signals."""
    try:
        payload = read_payload(
            input_path,
            input_format=input_format.value if input_format else None,
            cds=tuple(cds or ()),
        )
        try:
            config = LintConfig(
                max_length=max_length,
                gc_window=gc_window,
                gc_step=gc_step,
                gc_min_percent=gc_min,
                gc_max_percent=gc_max,
                cpg_window=cpg_window,
                cpg_step=cpg_step,
                cpg_hotspot_threshold=cpg_hotspot,
                forbidden_motifs=tuple(parse_motif(spec) for spec in motif or ()),
                splice_scan=splice_scan,
                rnafold=rnafold,
                structure_window=structure_window,
                structure_step=structure_step,
                structure_mfe_threshold=structure_mfe,
                structure_ensemble=structure_ensemble,
                screens=_screen_config(
                    enable=screen,
                    homopolymer_min=homopolymer_min,
                    repeat_min_copies=repeat_min_copies,
                    repeat_max_unit=repeat_max_unit,
                    orf_min_codons=orf_min_codons,
                    matrix_paths=tuple(matrix or ()),
                    screen_motifs=tuple(screen_motif or ()),
                    assembly_standards=tuple(
                        standard.value for standard in assembly_standard or ()
                    ),
                ),
                vector=_vector_architecture(vector_path),
            )
        except ValidationError as exc:
            message = "; ".join(error["msg"] for error in exc.errors())
            raise InputError(f"Invalid lint configuration: {message}") from exc
        report = lint_payload(payload, config)
        rendered_json = _json_report(report)
        _validate_output_paths(input_path, json_file, html_report)
        lint_outputs: dict[Path, str] = {}
        if json_file is not None:
            lint_outputs[json_file] = rendered_json + "\n"
        if html_report is not None:
            lint_outputs[html_report] = render_lint_html(report)
        if lint_outputs:
            _write_outputs_atomically(lint_outputs, input_path=input_path)
    except InputError as exc:
        _emit_input_error(exc, json_output=json_output)
        raise typer.Exit(2) from exc

    if json_output:
        typer.echo(rendered_json)
    else:
        gc = "n/a" if report.metrics.gc_percent is None else f"{report.metrics.gc_percent:g}%"
        typer.echo(
            f"{'PASS' if report.valid else 'FAIL'} {report.identifier}: "
            f"{report.metrics.length}/{report.configuration.max_length} nt, "
            f"{report.metrics.cpg_count} CpG, {gc} GC"
        )
        budget = report.metrics.vector_budget
        if budget is not None:
            typer.echo(
                f"Vector: {budget['genome_type']} {budget['representation']}, "
                f"unit {budget['accounted_unit_length']} nt vs {budget['limit']} nt "
                f"limit (headroom {budget['headroom']} nt)"
            )
            for component in budget["components"]:
                shorter = component["shorter_alternatives"]
                detail = (
                    f"; {len(shorter)} shorter alternative(s)" if shorter else ""
                )
                typer.echo(
                    f"  {component['kind']:<10} {component['name']:<16} "
                    f"{component['length']:>5} nt "
                    f"{'protected' if component['protected'] else 'replaceable'}{detail}"
                )
        for finding in report.validation.findings + report.findings:
            context = f" [{finding.region}]" if finding.region else ""
            typer.echo(f"{finding.severity.upper()} {finding.code}{context}: {finding.message}")
        if json_file is not None:
            typer.echo(f"JSON report: {json_file}")
        if html_report is not None:
            typer.echo(f"HTML report: {html_report}")
    raise typer.Exit(0 if report.valid else 1)


@app.command()
def optimize(
    input_path: Annotated[
        Path, typer.Option("--input", "-i", help="Single-record FASTA, GenBank, or EMBL")
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Optimized FASTA output")],
    cds: Annotated[
        list[str] | None,
        typer.Option(
            "--cds",
            help="FASTA CDS: START:END[:+|-], one-based inclusive; repeat for multiple CDSs",
        ),
    ] = None,
    input_format: Annotated[Format | None, typer.Option("--format")] = None,
    output_format: Annotated[
        Format | None,
        typer.Option(
            "--output-format",
            help="Output format; otherwise inferred from the output suffix (default FASTA)",
        ),
    ] = None,
    host: Annotated[str, typer.Option("--host", help="Bundled codon preference profile")] = "human",
    max_length: Annotated[int, typer.Option("--max-length")] = 4700,
    deplete_cpg: Annotated[
        bool,
        typer.Option("--deplete-cpg", help="Enable synonymous CpG depletion"),
    ] = True,
    rnafold_mode: Annotated[
        Literal["off", "auto", "seqfold", "required"],
        typer.Option(
            "--rnafold",
            help="Folding mode: off, auto, seqfold, or require native RNAfold",
        ),
    ] = "off",
    splice_scan: Annotated[
        bool,
        typer.Option(
            "--splice-scan/--no-splice-scan",
            help="Prevent creation of new splice-consensus motifs",
        ),
    ] = True,
    structure_window: Annotated[int, typer.Option("--structure-window")] = 150,
    structure_step: Annotated[int, typer.Option("--structure-step")] = 25,
    structure_mfe: Annotated[
        float, typer.Option("--structure-mfe", help="Flag folds at or below this kcal/mol value")
    ] = -15.0,
    motif: Annotated[
        list[str] | None,
        typer.Option(
            "--forbid-motif",
            help="Avoid introducing NAME=IUPAC_SEQUENCE on either strand; repeat as needed",
        ),
    ] = None,
    report_file: Annotated[
        Path | None,
        typer.Option("--report", help="JSON audit (.json) or HTML visual diff (.html/.htm)"),
    ] = None,
    json_report_file: Annotated[
        Path | None, typer.Option("--json-report", help="JSON mutation audit report path")
    ] = None,
    diff_report: Annotated[
        Path | None, typer.Option("--diff-report", help="Standalone HTML visual mutation diff")
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Write the audit report to stdout")
    ] = False,
    solver: Annotated[
        Solver,
        typer.Option(
            "--solver",
            help="greedy: best single-codon edit per step. global: bounded exhaustive "
            "search that can prove optimality over the declared constraints",
        ),
    ] = Solver.greedy,
    search_budget: Annotated[
        int,
        typer.Option(
            "--search-budget",
            help="Maximum codon assignments the global solver may attempt before "
            "returning its best feasible incumbent",
        ),
    ] = 10000,
    cpg_policy: Annotated[
        CpgPolicy,
        typer.Option(
            "--cpg-policy",
            help="strict-subset: never place a CpG at a coordinate that lacked one. "
            "count-only: allow relocation if the total count still falls",
        ),
    ] = CpgPolicy.strict_subset,
    profile_path: Annotated[
        Path | None,
        typer.Option("--profile", help="User-supplied versioned CodonProfile JSON file"),
    ] = None,
    vector_path: Annotated[
        Path | None,
        typer.Option("--vector", help="VectorArchitecture JSON for ITR-aware size accounting"),
    ] = None,
    protect: Annotated[
        list[str] | None,
        typer.Option("--protect", help="Immutable interval START:END[:NAME]; repeat as needed"),
    ] = None,
    gc3_min: Annotated[float, typer.Option("--gc3-min", help="Minimum per-CDS GC3 percent")] = 0.0,
    gc3_max: Annotated[
        float, typer.Option("--gc3-max", help="Maximum per-CDS GC3 percent")
    ] = 100.0,
    min_distinct_codons: Annotated[
        int, typer.Option("--min-distinct-codons", help="Minimum distinct sense codons per CDS")
    ] = 0,
    rare_codon_threshold: Annotated[
        float,
        typer.Option(
            "--rare-codon-threshold",
            help="Relative frequency below which a codon is rare",
        ),
    ] = 0.1,
    max_rare_codons: Annotated[
        int | None, typer.Option("--max-rare-codons", help="Maximum rare codons permitted per CDS")
    ] = None,
    screen: Annotated[
        bool, typer.Option("--screen/--no-screen", help="Run heuristic sequence screens")
    ] = False,
    homopolymer_min: Annotated[int, typer.Option("--homopolymer-min")] = 8,
    repeat_min_copies: Annotated[int, typer.Option("--repeat-min-copies")] = 4,
    repeat_max_unit: Annotated[int, typer.Option("--repeat-max-unit")] = 6,
    orf_min_codons: Annotated[int, typer.Option("--orf-min-codons")] = 30,
    matrix: Annotated[
        list[Path] | None,
        typer.Option("--matrix", help="ScoringMatrix JSON file; repeat as needed"),
    ] = None,
    screen_motif: Annotated[
        list[str] | None,
        typer.Option("--screen-motif", help="Report NAME=IUPAC_SEQUENCE as a screen hit"),
    ] = None,
    assembly_standard: Annotated[
        list[AssemblyStandard] | None,
        typer.Option("--assembly-standard", help="Screen a restriction-site library"),
    ] = None,
    structure_ensemble: Annotated[
        bool,
        typer.Option("--structure-ensemble", help="Report ensemble statistics, not only the MFE"),
    ] = False,
):
    """Deplete CpGs and write a sequence plus JSON and HTML reports."""
    report_suffix = report_file.suffix.lower() if report_file is not None else None
    report_html_path = report_file if report_suffix in (".html", ".htm") else None
    report_json_path = report_file if report_suffix == ".json" else None
    audit_path = json_report_file or report_json_path or Path(f"{output}.report.json")
    diff_path = diff_report or report_html_path or Path(f"{output}.diff.html")
    try:
        if report_file is not None and report_suffix not in (".json", ".html", ".htm"):
            raise InputError("--report must end in .json, .html, or .htm")
        if report_json_path is not None and json_report_file is not None:
            raise InputError("Use only one of --report .json and --json-report")
        if report_html_path is not None and diff_report is not None:
            raise InputError("Use only one of --report .html and --diff-report")
        if not deplete_cpg:
            raise InputError("CpG depletion is required by the optimize command")
        _validate_output_paths(input_path, output, audit_path, diff_path)
        payload = read_payload(
            input_path,
            input_format=input_format.value if input_format else None,
            cds=tuple(cds or ()),
        )
        custom_profile = (
            None if profile_path is None else read_json_model(profile_path, CodonProfile)
        )
        try:
            config = LintConfig(
                max_length=max_length,
                forbidden_motifs=tuple(parse_motif(spec) for spec in motif or ()),
                splice_scan=splice_scan,
                rnafold_mode=rnafold_mode,
                structure_window=structure_window,
                structure_step=structure_step,
                structure_mfe_threshold=structure_mfe,
                structure_ensemble=structure_ensemble,
                screens=_screen_config(
                    enable=screen,
                    homopolymer_min=homopolymer_min,
                    repeat_min_copies=repeat_min_copies,
                    repeat_max_unit=repeat_max_unit,
                    orf_min_codons=orf_min_codons,
                    matrix_paths=tuple(matrix or ()),
                    screen_motifs=tuple(screen_motif or ()),
                    assembly_standards=tuple(
                        standard.value for standard in assembly_standard or ()
                    ),
                ),
                vector=_vector_architecture(vector_path),
            )
            settings = OptimizeSettings(
                solver="global" if solver is Solver.global_ else "greedy",
                search_budget=search_budget,
                cpg_policy=(
                    "count_only" if cpg_policy is CpgPolicy.count_only else "strict_subset"
                ),
                protected_regions=tuple(_parse_protected(spec) for spec in protect or ()),
                gc3_min=gc3_min,
                gc3_max=gc3_max,
                min_distinct_codons=min_distinct_codons,
                rare_codon_threshold=rare_codon_threshold,
                max_rare_codons=max_rare_codons,
            )
        except ValidationError as exc:
            message = "; ".join(error["msg"] for error in exc.errors())
            raise InputError(f"Invalid optimization configuration: {message}") from exc
        optimized_payload, report = optimize_payload(
            payload,
            host=host,
            lint_config=config,
            settings=settings,
            profile=custom_profile,
        )
        rendered_json = _json_report(report)
        rendered_diff = render_optimize_html(payload, optimized_payload, report)
        selected_output_format = _sequence_output_format(output, output_format)
        sequence_output = _render_sequence_output(
            input_path, optimized_payload, selected_output_format
        )
        output_set = {
            audit_path: rendered_json + "\n",
            diff_path: rendered_diff,
        }
        if report.optimized:
            output_set[output] = sequence_output
            _write_outputs_atomically(output_set, input_path=input_path)
        else:
            _write_outputs_atomically(output_set, input_path=input_path, remove=(output,))
    except InputError as exc:
        _emit_input_error(exc, json_output=json_output)
        raise typer.Exit(2) from exc

    if json_output:
        typer.echo(rendered_json)
    else:
        typer.echo(
            f"{'PASS' if report.optimized else 'FAIL'} {report.identifier}: "
            f"CpG {report.cpg_count_before} → {report.cpg_count_after} "
            f"({report.cpg_reduction} removed), {len(report.mutations)} mutation(s)"
        )
        typer.echo(
            f"Solver: {report.search.solver} [{report.search.status}] "
            f"over {report.search.explored_states} state(s); "
            f"CpG policy {report.settings.cpg_policy}"
        )
        typer.echo(f"Profile: {report.profile.name} {report.profile.version}")
        for metric in report.codon_metrics_after:
            cai = "n/a" if metric.cai is None else f"{metric.cai:.3f}"
            pair = (
                "n/a" if metric.mean_pair_score is None else f"{metric.mean_pair_score:.3f}"
            )
            typer.echo(
                f"  CDS {metric.region}: CAI {cai}, GC3 {metric.gc3_percent:.1f}%, "
                f"{metric.distinct_codons} distinct codons, "
                f"{len(metric.rare_codon_positions)} rare, pair-bias {pair}"
            )
        budget = report.after_lint.metrics.vector_budget
        if budget is not None:
            typer.echo(
                f"Vector: {budget['genome_type']} {budget['representation']}, "
                f"unit {budget['accounted_unit_length']} nt vs {budget['limit']} nt "
                f"limit (headroom {budget['headroom']} nt)"
            )
        screen_hits = report.after_lint.metrics.screen_hits
        if screen_hits:
            kinds = sorted({hit.kind for hit in screen_hits})
            typer.echo(
                f"Screens: {len(screen_hits)} hit(s) across {', '.join(kinds)} "
                "(heuristic, not biological validation)"
            )
        for warning in report.warnings:
            typer.echo(f"WARNING {warning}")
        for regression in report.safety_regressions:
            typer.echo(f"ERROR {regression}")
        if report.optimized:
            typer.echo(f"Sequence ({selected_output_format.value}): {output}")
        else:
            typer.echo(f"Optimization failed; no sequence output is present at {output}")
        typer.echo(f"Audit: {audit_path}")
        typer.echo(f"Visual diff: {diff_path}")
    raise typer.Exit(0 if report.optimized else 1)


@app.command()
def report(
    input_path: Annotated[
        Path, typer.Option("--input", "-i", help="Original FASTA, GenBank, or EMBL sequence")
    ],
    optimized_path: Annotated[Path, typer.Option("--optimized", help="Optimized sequence file")],
    audit_path: Annotated[Path, typer.Option("--audit", help="Optimization JSON audit")],
    output: Annotated[Path, typer.Option("--output", "-o", help="HTML report output")],
    cds: Annotated[
        list[str] | None,
        typer.Option(
            "--cds",
            help="FASTA CDS used to independently verify translation; repeat as needed",
        ),
    ] = None,
    input_format: Annotated[Format | None, typer.Option("--format")] = None,
    optimized_format: Annotated[Format | None, typer.Option("--optimized-format")] = None,
):
    """Render a visual diff from matching source, sequence, and audit files."""
    try:
        _validate_output_paths(input_path, output)
        protected_paths = {input_path.resolve(), optimized_path.resolve(), audit_path.resolve()}
        if output.resolve() in protected_paths:
            raise InputError(
                "HTML report output cannot overwrite an input, sequence, or audit file"
            )
        if input_path.resolve() == optimized_path.resolve():
            raise InputError("Original and optimized sequence paths must be different")
        source = read_payload(
            input_path,
            input_format=input_format.value if input_format else None,
            cds=tuple(cds or ()),
        )
        optimized = read_payload(
            optimized_path,
            input_format=optimized_format.value if optimized_format else None,
        )
        try:
            audit_data = json.loads(audit_path.read_text(encoding="utf-8"))
        except OSError as exc:
            raise InputError(f"Cannot read audit {audit_path}: {exc}") from exc
        except json.JSONDecodeError as exc:
            raise InputError(f"Invalid JSON audit {audit_path}: {exc}") from exc
        declared_status = (
            audit_data.pop("optimized", None) if isinstance(audit_data, dict) else None
        )
        try:
            audit = OptimizeReport.model_validate(audit_data)
        except ValidationError as exc:
            raise InputError(f"Invalid optimization audit: {exc}") from exc
        if declared_status is not None and declared_status != audit.optimized:
            raise InputError("Audit optimized status does not match its checks")
        verify_optimization_audit(source, optimized, audit)
        rendered = render_optimize_html(source, optimized, audit)
        _write_output(output, rendered, input_path=input_path)
    except InputError as exc:
        _emit_input_error(exc, json_output=False)
        raise typer.Exit(2) from exc
    typer.echo(f"HTML report: {output}")
