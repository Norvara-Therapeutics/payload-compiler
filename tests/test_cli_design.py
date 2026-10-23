import json

import pytest
from Bio.Seq import Seq
from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.exceptions import InputError
from payload_compiler.models import CodingRegion, Payload
from payload_compiler.optimize import optimize_payload
from payload_compiler.parser import read_payload
from payload_compiler.validation import extract_cds

runner = CliRunner()


def _write_cpg_source(tmp_path):
    path = tmp_path / "source.fasta"
    path.write_text(">demo CpG fixture\nATGCGTTAA\n", encoding="utf-8")
    return path


def test_profiles_lists_only_substantiated_profiles():
    result = runner.invoke(app, ["profiles"])
    assert result.exit_code == 0
    assert "human:" in result.stdout
    # The withdrawn alias claimed a cell line the bundled table never measured.
    assert "human-hek293" not in result.stdout


def test_retired_cell_line_alias_is_an_actionable_error():
    source = Payload(
        identifier="retired", sequence="ATGCGTTAA", source_format="fasta",
        coding_regions=(CodingRegion(name="main", start=0, end=9),),
    )
    with pytest.raises(InputError) as excinfo:
        optimize_payload(source, host="human-hek293")
    message = str(excinfo.value)
    assert "never HEK293-derived" in message
    assert "--profile" in message


def test_optimizer_invariants_and_determinism(tmp_path):
    source_path = _write_cpg_source(tmp_path)
    source = read_payload(source_path, cds=("1:9",))
    first, first_report = optimize_payload(source, host="human")
    second, second_report = optimize_payload(source, host="human")

    assert first.sequence == second.sequence
    assert first_report.model_dump(mode="json") == second_report.model_dump(mode="json")
    assert first_report.cpg_count_after <= first_report.cpg_count_before
    assert first_report.protein_translation_preserved
    assert (
        Seq(extract_cds(source, source.coding_regions[0])).translate()
        == Seq(extract_cds(first, first.coding_regions[0])).translate()
    )
    assert [mutation.index for mutation in first_report.mutations] == list(
        range(1, len(first_report.mutations) + 1)
    )


def test_optimize_html_report_and_standalone_regeneration(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "optimized.fasta"
    html_report = tmp_path / "optimized.html"
    audit = tmp_path / "optimized.fasta.report.json"
    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(source),
            "--cds",
            "1:9",
            "--host",
            "human",
            "--deplete-cpg",
            "--output",
            str(optimized),
            "--report",
            str(html_report),
        ],
    )
    assert result.exit_code == 0, result.output
    assert optimized.exists() and audit.exists() and html_report.exists()
    audit_data = json.loads(audit.read_text(encoding="utf-8"))
    assert audit_data["optimized"]
    assert audit_data["host"] == "human"
    first_html = html_report.read_bytes()

    rebuilt = tmp_path / "rebuilt.html"
    report_result = runner.invoke(
        app,
        [
            "report",
            "--input",
            str(source),
            "--optimized",
            str(optimized),
            "--audit",
            str(audit),
            "--output",
            str(rebuilt),
        ],
    )
    assert report_result.exit_code == 0, report_result.output
    assert rebuilt.read_bytes() == first_html


def test_optimize_json_report_option_remains_supported(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "legacy.fasta"
    audit = tmp_path / "legacy.audit.json"
    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(source),
            "--cds",
            "1:9",
            "--output",
            str(optimized),
            "--report",
            str(audit),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(audit.read_text(encoding="utf-8"))["optimized"]
    assert (tmp_path / "legacy.fasta.diff.html").exists()


def test_report_rejects_mismatched_sequence_audit(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "optimized.fasta"
    result = runner.invoke(
        app,
        ["optimize", "--input", str(source), "--cds", "1:9", "--output", str(optimized)],
    )
    assert result.exit_code == 0, result.output
    optimized.write_text(">demo\nATGGCTTAA\n", encoding="utf-8")

    html_report = tmp_path / "mismatch.html"
    report_result = runner.invoke(
        app,
        [
            "report",
            "--input",
            str(source),
            "--optimized",
            str(optimized),
            "--audit",
            str(tmp_path / "optimized.fasta.report.json"),
            "--output",
            str(html_report),
        ],
    )
    assert report_result.exit_code == 2
    assert "SHA-256" in report_result.output
    assert not html_report.exists()


def test_report_rejects_internally_inconsistent_audit(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "optimized.fasta"
    result = runner.invoke(
        app,
        ["optimize", "--input", str(source), "--cds", "1:9", "--output", str(optimized)],
    )
    assert result.exit_code == 0, result.output
    audit = tmp_path / "optimized.fasta.report.json"
    data = json.loads(audit.read_text(encoding="utf-8"))
    data["cpg_reduction"] = 999
    audit.write_text(json.dumps(data), encoding="utf-8")

    html_report = tmp_path / "tampered.html"
    report_result = runner.invoke(
        app,
        [
            "report",
            "--input",
            str(source),
            "--optimized",
            str(optimized),
            "--audit",
            str(audit),
            "--output",
            str(html_report),
        ],
    )
    assert report_result.exit_code == 2
    assert "CpG reduction" in report_result.output
    assert not html_report.exists()


def test_failed_optimization_does_not_write_sequence(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "must-not-exist.fasta"
    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(source),
            "--cds",
            "1:9",
            "--forbid-motif",
            "terminal_stop=TAA",
            "--output",
            str(optimized),
        ],
    )
    assert result.exit_code == 1, result.output
    assert not optimized.exists()
    assert (tmp_path / "must-not-exist.fasta.report.json").exists()
    assert (tmp_path / "must-not-exist.fasta.diff.html").exists()


def test_failed_optimization_removes_stale_sequence_output(tmp_path):
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "stale.fasta"
    optimized.write_text(">old-result\nAAAA\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(source),
            "--cds",
            "1:9",
            "--forbid-motif",
            "terminal_stop=TAA",
            "--output",
            str(optimized),
        ],
    )

    assert result.exit_code == 1, result.output
    assert not optimized.exists()
    assert "no sequence output is present" in result.output


def test_atomic_output_set_rolls_back_on_commit_failure(tmp_path, monkeypatch):
    from payload_compiler import cli

    input_path = tmp_path / "input.fasta"
    first = tmp_path / "first.txt"
    second = tmp_path / "second.txt"
    input_path.write_text(">input\nACGT\n", encoding="utf-8")
    first.write_text("old-first", encoding="utf-8")
    second.write_text("old-second", encoding="utf-8")
    real_replace = cli.os.replace
    commit_count = 0

    def fail_second_commit(source, destination):
        nonlocal commit_count
        if str(source).endswith(".tmp"):
            commit_count += 1
            if commit_count == 2:
                raise OSError("simulated commit failure")
        return real_replace(source, destination)

    monkeypatch.setattr(cli.os, "replace", fail_second_commit)
    with pytest.raises(InputError, match="Cannot commit output set"):
        cli._write_outputs_atomically(
            {first: "new-first", second: "new-second"}, input_path=input_path
        )

    assert first.read_text(encoding="utf-8") == "old-first"
    assert second.read_text(encoding="utf-8") == "old-second"


def test_folding_auto_uses_seqfold_without_native_executable(tmp_path, monkeypatch):
    from payload_compiler import analysis

    monkeypatch.setattr(analysis.shutil, "which", lambda _: None)
    source = _write_cpg_source(tmp_path)
    optimized = tmp_path / "auto.fasta"
    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(source),
            "--cds",
            "1:9",
            "--rnafold",
            "auto",
            "--output",
            str(optimized),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert optimized.exists()
    audit = json.loads(result.stdout)
    assert audit["optimized"]
    assert audit["after_lint"]["metrics"]["folding_tool_version"].startswith("seqfold ")
    assert audit["after_lint"]["metrics"]["folding_windows"]
