import json

import pytest
from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.exceptions import InputError
from payload_compiler.lint import lint_payload, parse_motif
from payload_compiler.models import LintConfig, MotifDefinition, Payload
from payload_compiler.reporting import render_lint_html

runner = CliRunner()


def payload(sequence: str, identifier: str = "test") -> Payload:
    return Payload(identifier=identifier, sequence=sequence, source_format="fasta")


def test_length_cpg_and_gc_metrics():
    report = lint_payload(
        payload("ACGCGTATAT"),
        LintConfig(gc_window=4, gc_step=3, cpg_window=5, cpg_step=5),
    )
    assert report.valid
    assert report.metrics.length == 10
    assert report.metrics.cpg_count == 2
    assert report.metrics.cpg_positions == (2, 4)
    assert report.metrics.cpg_observed_expected == 5.0
    assert report.metrics.gc_percent == 40.0
    assert [(w.start, w.end, w.gc_percent) for w in report.metrics.gc_windows] == [
        (1, 4, 75.0),
        (4, 7, 50.0),
        (7, 10, 0.0),
    ]


def test_final_window_is_end_anchored():
    report = lint_payload(
        payload("ACGTACGTACG"),
        LintConfig(gc_window=4, gc_step=4),
    )
    assert [(window.start, window.end) for window in report.metrics.gc_windows] == [
        (1, 4),
        (5, 8),
        (8, 11),
    ]


def test_short_sequence_has_one_window():
    report = lint_payload(payload("ACG"), LintConfig(gc_window=100, cpg_window=100))
    assert [(window.start, window.end) for window in report.metrics.gc_windows] == [(1, 3)]


def test_ambiguous_bases_are_excluded_from_gc_denominator():
    report = lint_payload(payload("GCNN"), LintConfig(gc_window=4, gc_step=4))
    assert report.metrics.gc_percent == 100.0
    assert report.metrics.gc_windows[0].ambiguous_bases == 2


def test_all_ambiguous_gc_is_unavailable():
    report = lint_payload(payload("NNNN"), LintConfig(gc_window=4, gc_step=4))
    assert report.metrics.gc_percent is None
    assert not any(f.code in {"LOW_GC_WINDOW", "HIGH_GC_WINDOW"} for f in report.findings)


def test_length_limit_is_an_error():
    report = lint_payload(payload("ACGT"), LintConfig(max_length=3))
    assert not report.valid
    assert report.metrics.length_overage == 1
    assert report.findings[0].code == "LENGTH_LIMIT_EXCEEDED"


def test_gc_and_cpg_findings_are_warnings():
    report = lint_payload(
        payload("CGCGCGCG"),
        LintConfig(
            gc_window=4,
            gc_step=4,
            gc_max_percent=60,
            cpg_window=4,
            cpg_step=4,
            cpg_hotspot_threshold=2,
        ),
    )
    assert report.valid
    assert [finding.code for finding in report.findings] == [
        "HIGH_GC_WINDOW",
        "CPG_HOTSPOT",
    ]
    assert [(window.start, window.cpg_count) for window in report.metrics.cpg_hotspots] == [
        (1, 2),
        (5, 2),
    ]


def test_motifs_scan_both_strands_and_palindromes_once():
    config = LintConfig(forbidden_motifs=(parse_motif("BsaI=GGTCTC"), parse_motif("EcoRI=GAATTC")))
    report = lint_payload(payload("AAGGTCTCTTGAGACCCAAGAATTCT"), config)
    matches = [(m.name, m.start, m.strand) for m in report.metrics.motif_matches]
    assert matches == [("BsaI", 3, "+"), ("BsaI", 11, "-"), ("EcoRI", 20, "both")]
    assert not report.valid
    assert [finding.region for finding in report.findings] == ["BsaI", "EcoRI"]


def test_overlapping_iupac_motifs():
    report = lint_payload(payload("AAAA"), LintConfig(forbidden_motifs=(parse_motif("run=AAN"),)))
    assert [(m.start, m.end) for m in report.metrics.motif_matches] == [(1, 3), (2, 4)]


@pytest.mark.parametrize("spec", ["BsaI", "=GGTCTC", "BsaI=", "bad=ACX"])
def test_invalid_motif_specs(spec):
    with pytest.raises(InputError):
        parse_motif(spec)


def test_direct_motif_model_normalizes_and_validates():
    assert MotifDefinition(name=" site ", sequence=" acn ").model_dump() == {
        "name": "site",
        "sequence": "ACN",
    }
    with pytest.raises(ValueError):
        MotifDefinition(name="site", sequence="ACX")


def test_html_report_content_and_escaping():
    report = lint_payload(
        payload("ACGCGT", identifier="<script>alert(1)</script>"),
        LintConfig(forbidden_motifs=(parse_motif("site=CG"),)),
    )
    first = render_lint_html(report)
    assert first == render_lint_html(report)
    assert first.startswith("<!doctype html>")
    assert "<script>alert(1)</script>" not in first
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in first
    assert "https://" not in first
    assert "Sequence overview" in first
    assert "Hotspot windows" in first


def test_cli_writes_deterministic_json_and_html(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">sample\nACGCGT\n")
    json_path = tmp_path / "report.json"
    html_path = tmp_path / "report.html"
    args = [
        "lint",
        "--input",
        str(source),
        "--max-length",
        "20",
        "--forbid-motif",
        "site=TTTT",
        "--json-output",
        str(json_path),
        "--report",
        str(html_path),
        "--json",
    ]
    first = runner.invoke(app, args)
    assert first.exit_code == 0, first.output
    stdout_report = json.loads(first.stdout)
    assert stdout_report["metrics"]["cpg_positions"] == [2, 4]
    assert stdout_report["valid"]
    first_json, first_html = json_path.read_bytes(), html_path.read_bytes()
    second = runner.invoke(app, args)
    assert second.exit_code == 0
    assert json_path.read_bytes() == first_json
    assert html_path.read_bytes() == first_html


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--max-length", "0"],
        ["--gc-min", "80", "--gc-max", "20"],
        ["--gc-window", "10", "--gc-step", "11"],
        ["--cpg-window", "1"],
        ["--cpg-window", "10", "--cpg-step", "11"],
        ["--forbid-motif", "bad=X"],
        ["--forbid-motif", "x=ATG", "--forbid-motif", "X=CGT"],
    ],
)
def test_cli_invalid_configuration(tmp_path, extra_args):
    source = tmp_path / "input.fa"
    source.write_text(">sample\nACGT\n")
    result = runner.invoke(app, ["lint", "--input", str(source), "--json", *extra_args])
    assert result.exit_code == 2
    assert json.loads(result.stdout)["error"]["code"] == "INPUT_ERROR"


def test_cli_protects_input_and_distinct_outputs(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">sample\nACGT\n")
    before = source.read_bytes()
    same_as_input = runner.invoke(app, ["lint", "--input", str(source), "--report", str(source)])
    assert same_as_input.exit_code == 2
    assert source.read_bytes() == before
    output = tmp_path / "same.out"
    same_outputs = runner.invoke(
        app,
        [
            "lint",
            "--input",
            str(source),
            "--json-output",
            str(output),
            "--report",
            str(output),
        ],
    )
    assert same_outputs.exit_code == 2
    assert not output.exists()
    json_output = tmp_path / "should-not-exist.json"
    mixed_outputs = runner.invoke(
        app,
        [
            "lint",
            "--input",
            str(source),
            "--json-output",
            str(json_output),
            "--report",
            str(source),
        ],
    )
    assert mixed_outputs.exit_code == 2
    assert not json_output.exists()


def test_cli_lint_failure_exit_code(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">sample\nACGT\n")
    result = runner.invoke(app, ["lint", "--input", str(source), "--max-length", "3", "--json"])
    assert result.exit_code == 1
    assert not json.loads(result.stdout)["valid"]
