"""Golden artifacts lock down serialized output across supported formats."""

import hashlib
import json
from pathlib import Path

from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.optimize import optimize_payload
from payload_compiler.parser import read_payload

FIXTURES = Path(__file__).parent / "fixtures"
runner = CliRunner()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_fasta_genbank_json_and_html_golden_artifacts(tmp_path):
    expected_hashes = dict(
        line.split() for line in (FIXTURES / "golden_artifacts.sha256").read_text().splitlines()
    )
    fasta_input = FIXTURES / "golden_input.fasta"
    fasta_payload = read_payload(fasta_input, cds=("1:12",))
    genbank_payload = read_payload(FIXTURES / "golden_input.gb")

    assert genbank_payload.sequence == fasta_payload.sequence
    genbank_locations = [
        (region.start, region.end, region.strand)
        for region in genbank_payload.coding_regions
    ]
    assert genbank_locations == [(0, 12, 1)]
    fasta_optimized, fasta_report = optimize_payload(fasta_payload, host="human")
    genbank_optimized, genbank_report = optimize_payload(genbank_payload, host="human")
    assert fasta_optimized.sequence == genbank_optimized.sequence == "ATGGCTGCTTAA"
    assert fasta_report.cpg_count_after == genbank_report.cpg_count_after == 0

    output = tmp_path / "golden.fasta"
    result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(fasta_input),
            "--cds",
            "1:12",
            "--host",
            "human",
            "--output",
            str(output),
            "--report",
            str(tmp_path / "golden.html"),
        ],
    )
    assert result.exit_code == 0, result.output

    artifacts = {
        "fasta": output,
        "json": tmp_path / "golden.fasta.report.json",
        "html": tmp_path / "golden.html",
    }
    assert output.read_bytes() == (FIXTURES / "golden_optimized.fasta").read_bytes()
    assert _sha256(artifacts["fasta"]) == expected_hashes["fasta"]
    # The audit and visual diff embed RuntimeMetadata (interpreter, platform and
    # package versions), so their digests are machine-specific and deliberately
    # absent from the fixture. Determinism is asserted by repeating the run here
    # instead of pinning a hash that could only ever match one machine.
    assert set(expected_hashes) == {"fasta"}
    repeat_dir = tmp_path / "repeat"
    repeat_dir.mkdir()
    repeat_output = repeat_dir / "golden.fasta"
    repeat = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(fasta_input),
            "--cds",
            "1:12",
            "--host",
            "human",
            "--output",
            str(repeat_output),
            "--report",
            str(repeat_dir / "golden.html"),
        ],
    )
    assert repeat.exit_code == 0, repeat.output
    for name, produced in (
        ("fasta", repeat_output),
        ("json", repeat_dir / "golden.fasta.report.json"),
        ("html", repeat_dir / "golden.html"),
    ):
        assert produced.read_bytes() == artifacts[name].read_bytes(), name
    audit = json.loads(artifacts["json"].read_text(encoding="utf-8"))
    assert audit["cpg_count_before"] == 1
    assert audit["cpg_count_after"] == 0
    assert audit["mutations"][0]["base_changes"] == [
        {"position": 6, "before": "C", "after": "T"}
    ]
    assert audit["schema_version"] == "2.0"
    assert set(audit["runtime"]) == {
        "python",
        "platform",
        "biopython",
        "pydantic",
        "seqfold",
        "typer",
    }

    repeated = tmp_path / "golden-repeated.fasta"
    repeated_result = runner.invoke(
        app,
        [
            "optimize",
            "--input",
            str(fasta_input),
            "--cds",
            "1:12",
            "--host",
            "human",
            "--output",
            str(repeated),
            "--report",
            str(tmp_path / "golden-repeated.html"),
        ],
    )
    assert repeated_result.exit_code == 0, repeated_result.output
    assert artifacts["json"].read_bytes() == (
        tmp_path / "golden-repeated.fasta.report.json"
    ).read_bytes()
    assert artifacts["html"].read_bytes() == (tmp_path / "golden-repeated.html").read_bytes()


def test_genbank_output_preserves_annotations(tmp_path):
    source = FIXTURES / "golden_input.gb"
    output = tmp_path / "optimized.gb"

    result = runner.invoke(
        app,
        ["optimize", "--input", str(source), "--output", str(output)],
    )

    assert result.exit_code == 0, result.output
    optimized = read_payload(output)
    original = read_payload(source)
    assert optimized.sequence == "ATGGCTGCTTAA"
    assert optimized.features == original.features
    assert optimized.coding_regions == original.coding_regions
