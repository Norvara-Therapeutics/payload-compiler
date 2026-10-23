"""Pins the phiX174 example quoted in the README quickstart.

Gene D of Escherichia phage phiX174 (NC_001422.1) is annotated at 390..848: a
459 nt in-frame CDS carrying 36 CpG sites. The quickstart depletes all of them
without altering the encoded protein or any base outside the CDS.
"""

import json
from pathlib import Path

from Bio.Seq import Seq
from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.models import LintConfig
from payload_compiler.optimize import optimize_payload
from payload_compiler.parser import read_payload

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "phix174.fasta"
GENE_D = (390, 848)

runner = CliRunner()


def test_example_genome_is_the_published_record():
    payload = read_payload(EXAMPLE)
    assert payload.identifier == "NC_001422.1"
    assert len(payload.sequence) == 5386
    assert set(payload.sequence) <= set("ACGT")


def test_gene_d_is_a_clean_in_frame_cds():
    sequence = read_payload(EXAMPLE).sequence
    start, end = GENE_D
    coding = sequence[start - 1 : end]
    assert len(coding) == 459
    protein = str(Seq(coding).translate(table=1))
    assert protein.endswith("*")
    assert "*" not in protein[:-1]
    assert coding.count("CG") == 36


def test_quickstart_depletes_every_cpg_in_the_cds_and_touches_nothing_else():
    source = read_payload(EXAMPLE, cds=(f"{GENE_D[0]}:{GENE_D[1]}",))
    optimized, report = optimize_payload(
        source, lint_config=LintConfig(max_length=len(source.sequence))
    )

    start, end = GENE_D
    before, after = source.sequence, optimized.sequence
    assert len(before) == len(after)
    assert before[: start - 1] == after[: start - 1]
    assert before[end:] == after[end:]
    assert after[start - 1 : end].count("CG") == 0

    assert Seq(before[start - 1 : end]).translate(table=1) == Seq(
        after[start - 1 : end]
    ).translate(table=1)

    # Whole-payload counts quoted in the README.
    assert report.cpg_count_before == 267
    assert report.cpg_count_after == 231
    assert report.cpg_reduction == 36
    assert len(report.mutations) == 35
    assert report.introduced_cpg_positions == ()
    assert report.optimized


def test_quickstart_command_line_reproduces_the_quoted_summary(tmp_path):
    output = tmp_path / "geneD.fasta"
    result = runner.invoke(
        app,
        [
            "optimize", "--input", str(EXAMPLE),
            "--cds", f"{GENE_D[0]}:{GENE_D[1]}",
            "--max-length", "5386",
            "--output", str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "PASS NC_001422.1: CpG 267 → 231 (36 removed), 35 mutation(s)" in result.stdout
    assert "Solver: greedy [locally_optimal] over 36 state(s); CpG policy strict_subset" in (
        result.stdout
    )
    assert output.exists()
    audit = json.loads((tmp_path / "geneD.fasta.report.json").read_text(encoding="utf-8"))
    assert audit["protein_translation_preserved"] is True
    assert audit["safety_checks_passed"] is True


def test_whole_genome_exceeds_the_default_payload_budget():
    """5,386 nt is over the 4,700 nt default, so lint reports a length error."""
    result = runner.invoke(app, ["lint", "--input", str(EXAMPLE)])
    assert result.exit_code == 1
    assert "LENGTH_LIMIT_EXCEEDED" in result.stdout
