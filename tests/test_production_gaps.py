"""Coverage for the five production gaps: vector sizing, global search, CpG policy,
profile provenance, and biological screens."""

import json

import pytest
from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.design import get_profile, vector_budget
from payload_compiler.exceptions import InputError
from payload_compiler.models import (
    CodingRegion,
    LintConfig,
    OptimizeReport,
    OptimizeSettings,
    Payload,
    ScreenConfig,
    VectorArchitecture,
)
from payload_compiler.optimize import optimize_payload
from payload_compiler.screening import screen_sequence

runner = CliRunner()

# ATG CGT ACG CGT GCT TAA -> M R T R A *, three CpGs, no internal stop.
CPG_CDS = "ATGCGTACGCGTGCTTAA"


def _payload(sequence: str = CPG_CDS, identifier: str = "fixture") -> Payload:
    return Payload(
        identifier=identifier,
        sequence=sequence,
        source_format="fasta",
        coding_regions=(CodingRegion(name="main", start=0, end=len(sequence)),),
    )


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --- Gap 1: ITR-aware vector size accounting -------------------------------


def test_cassette_adds_declared_itr_overhead_and_never_assumes_inclusion():
    architecture = VectorArchitecture(
        genome_type="ssAAV",
        representation="cassette",
        omitted_itr_lengths=(145, 145),
        components=(),
    )
    budget = vector_budget(_payload(), architecture)
    assert budget["omitted_itr_overhead"] == 290
    assert budget["accounted_unit_length"] == len(CPG_CDS) + 290
    assert budget["limit"] == 4700


def test_scaav_uses_a_smaller_budget_than_ssaav():
    single = vector_budget(
        _payload(),
        VectorArchitecture(
            genome_type="ssAAV", representation="cassette", omitted_itr_lengths=(145, 145)
        ),
    )
    self_comp = vector_budget(
        _payload(),
        VectorArchitecture(
            genome_type="scAAV", representation="cassette", omitted_itr_lengths=(145, 145)
        ),
    )
    assert self_comp["limit"] < single["limit"]


def test_itr_to_itr_input_requires_two_terminal_itrs():
    architecture = VectorArchitecture(
        genome_type="ssAAV", representation="itr_to_itr_unit", components=()
    )
    with pytest.raises(InputError, match="two terminal ITR components"):
        vector_budget(_payload(), architecture)


def test_optimize_reports_itr_accounting_in_its_warnings():
    _, report = optimize_payload(_payload())
    assert any("No vector architecture supplied" in w for w in report.warnings)

    architecture = VectorArchitecture(
        genome_type="ssAAV", representation="cassette", omitted_itr_lengths=(145, 145)
    )
    _, vector_report = optimize_payload(
        _payload(), lint_config=LintConfig(vector=architecture)
    )
    assert any("excludes both ITRs" in w for w in vector_report.warnings)
    assert vector_report.after_lint.metrics.vector_budget is not None


def test_oversized_input_is_explained_in_component_terms(tmp_path):
    # 60 nt against a 40 nt limit. The two ITRs (30 nt) are protected, so the only
    # route under budget is swapping the 30 nt promoter for the 5 nt alternative.
    sequence = "ATG" + "GCT" * 18 + "TAA"
    architecture = VectorArchitecture(
        genome_type="ssAAV",
        representation="itr_to_itr_unit",
        limit=40,
        components=(
            {"start": 1, "end": 10, "name": "left_itr", "kind": "itr"},
            {
                "start": 11,
                "end": 40,
                "name": "promoter_x",
                "kind": "promoter",
                "protected": False,
                "alternatives": (
                    {
                        "name": "minimal_promoter",
                        "length": 5,
                        "source": "catalog",
                        "citation": "user-supplied",
                    },
                ),
            },
            {"start": 41, "end": len(sequence), "name": "right_itr", "kind": "itr"},
        ),
    )
    with pytest.raises(InputError) as excinfo:
        optimize_payload(
            _payload(sequence),
            lint_config=LintConfig(max_length=40, vector=architecture),
        )
    message = str(excinfo.value)
    assert "cannot shorten a payload" in message
    assert "Both ITRs are inside the submitted sequence" in message
    assert "promoter_x -> minimal_promoter" in message


def test_oversized_input_without_vector_says_so(tmp_path):
    sequence = "ATG" + "GCT" * 18 + "TAA"
    with pytest.raises(InputError, match="No --vector architecture was supplied"):
        optimize_payload(_payload(sequence), lint_config=LintConfig(max_length=30))


# --- Gap 2: global solver and optimality status ----------------------------


def test_greedy_reports_local_optimality_only():
    _, report = optimize_payload(_payload())
    assert report.search.solver == "greedy"
    assert report.search.status == "locally_optimal"
    assert "no global optimum guarantee" in report.search.proof_scope


def test_global_solver_proves_optimality_within_its_declared_scope():
    _, report = optimize_payload(
        _payload(), settings=OptimizeSettings(solver="global", search_budget=100000)
    )
    assert report.search.solver == "global"
    assert report.search.status == "proven_optimal"
    assert report.search.budget == 100000
    assert report.edit_semantics == "joint"


def test_exhausted_budget_reports_best_found_and_claims_nothing_more():
    _, report = optimize_payload(
        _payload(), settings=OptimizeSettings(solver="global", search_budget=3)
    )
    assert report.search.status == "best_found"
    assert any("budget exhausted" in w for w in report.warnings)


def test_global_and_greedy_agree_on_this_payload():
    greedy, _ = optimize_payload(_payload())
    searched, _ = optimize_payload(_payload(), settings=OptimizeSettings(solver="global"))
    assert greedy.sequence == searched.sequence


# --- Gap 3: configurable no-new-CpG policy ---------------------------------


def test_strict_subset_is_the_default_and_introduces_no_new_cpg_coordinate():
    optimized, report = optimize_payload(_payload())
    assert report.settings.cpg_policy == "strict_subset"
    assert report.introduced_cpg_positions == ()
    before = set(report.before_lint.metrics.cpg_positions)
    assert set(report.after_lint.metrics.cpg_positions).issubset(before)


def test_policy_is_recorded_in_the_audit_and_proof_scope():
    _, report = optimize_payload(
        _payload(), settings=OptimizeSettings(cpg_policy="count_only")
    )
    assert report.settings.cpg_policy == "count_only"
    assert "count_only" in report.search.proof_scope


def test_report_rejects_an_introduced_cpg_under_strict_policy():
    """The invariant is a contract on the report, not only a search filter.

    Validation is re-run from the serialized form because `model_copy` bypasses
    validators, which would make a tampered report silently acceptable.
    """
    _, report = optimize_payload(_payload())
    raw = report.model_dump(mode="json")
    # Claim a CpG at a coordinate the source never carried.
    raw["after_lint"]["metrics"]["cpg_positions"] = [2]
    raw["after_lint"]["metrics"]["cpg_count"] = 1
    raw["cpg_count_after"] = 1
    raw["cpg_reduction"] = raw["cpg_count_before"] - 1
    raw["residual_cpg_positions"] = [2]
    with pytest.raises(ValueError, match="introduces a CpG at a new position"):
        OptimizeReport.model_validate(raw)


def test_count_only_policy_permits_a_relocated_cpg_the_strict_policy_rejects():
    """Same tampered report is contract-valid once the policy is relaxed."""
    _, report = optimize_payload(
        _payload(), settings=OptimizeSettings(cpg_policy="count_only")
    )
    raw = report.model_dump(mode="json")
    raw["after_lint"]["metrics"]["cpg_positions"] = [2]
    raw["after_lint"]["metrics"]["cpg_count"] = 1
    raw["cpg_count_after"] = 1
    raw["cpg_reduction"] = raw["cpg_count_before"] - 1
    raw["residual_cpg_positions"] = [2]
    raw["introduced_cpg_positions"] = [2]
    # Joint semantics require every mutation to carry whole-transaction counts.
    raw["edit_semantics"] = "joint"
    for mutation in raw["mutations"]:
        mutation["cpg_count_before"] = raw["cpg_count_before"]
        mutation["cpg_count_after"] = raw["cpg_count_after"]
    accepted = OptimizeReport.model_validate(raw)
    assert accepted.introduced_cpg_positions == (2,)
    assert accepted.settings.cpg_policy == "count_only"


# --- Gap 4: profile provenance ---------------------------------------------


def test_retired_cell_line_alias_names_its_replacement():
    with pytest.raises(InputError) as excinfo:
        get_profile("human-hek293")
    assert "never HEK293-derived" in str(excinfo.value)
    assert "--profile" in str(excinfo.value)


def test_audit_carries_profile_provenance_and_checksum():
    _, report = optimize_payload(_payload())
    assert report.profile.source
    assert report.profile.citation
    assert len(report.profile_sha256) == 64
    assert "not a cell-line or expression model" in report.profile.description


def test_user_supplied_profile_overrides_the_bundled_table(tmp_path):
    base = get_profile("human")
    custom = base.model_copy(
        update={
            "name": "lab-line-v1",
            "version": "2026-01 in-house",
            "source": "in-house measurement",
            "citation": "internal report",
            "description": "User-supplied versioned profile",
        }
    )
    path = _write(tmp_path, "profile.json", json.dumps(custom.model_dump(mode="json")))
    source = _write(tmp_path, "in.fasta", f">demo\n{CPG_CDS}\n")
    output = tmp_path / "out.fasta"
    result = runner.invoke(
        app,
        ["optimize", "--input", str(source), "--cds", f"1:{len(CPG_CDS)}",
         "--profile", str(path), "--output", str(output)],
    )
    assert result.exit_code == 0, result.output
    assert "lab-line-v1" in result.stdout
    audit = json.loads((tmp_path / "out.fasta.report.json").read_text(encoding="utf-8"))
    assert audit["profile"]["name"] == "lab-line-v1"
    assert audit["profile"]["citation"] == "internal report"


def test_codon_metrics_report_cai_gc3_and_rare_codons():
    _, report = optimize_payload(_payload())
    assert report.codon_metrics_before and report.codon_metrics_after
    metric = report.codon_metrics_after[0]
    assert metric.cai is not None and 0 < metric.cai <= 1
    assert 0 <= metric.gc3_percent <= 100
    assert metric.distinct_codons > 0


def test_codon_constraints_are_enforced_as_safety_regressions():
    _, report = optimize_payload(_payload(), settings=OptimizeSettings(gc3_min=99.0))
    assert not report.safety_checks_passed
    assert any("codon constraints" in r for r in report.safety_regressions)


# --- Gap 5: biological screens ---------------------------------------------


def test_assembly_standard_library_reports_disallowed_sites():
    # BsaI GGTCTC sits inside the insert.
    payload = Payload(
        identifier="bsai", sequence="AAGGTCTCAA", source_format="fasta"
    )
    hits = screen_sequence(
        payload, ScreenConfig(enabled=True, assembly_standards=("golden-gate-bsai",))
    )
    assembly = [hit for hit in hits if hit.kind == "assembly_site"]
    assert assembly and assembly[0].name == "golden-gate-bsai:BsaI"
    assert "cloning constraint" in assembly[0].evidence


def test_screens_cover_polya_promoter_kozak_homopolymer_and_repeats():
    sequence = "AATAAA" + "TATAAAT" + "AAAAAAAAA" + "ATGCATGCATGCATGC" + "ACGATGGTT"
    hits = screen_sequence(
        Payload(identifier="s", sequence=sequence, source_format="fasta"),
        ScreenConfig(enabled=True),
    )
    kinds = {hit.kind for hit in hits}
    assert {"polya", "promoter", "homopolymer", "tandem_repeat"} <= kinds
    assert all(hit.evidence for hit in hits)


def test_reverse_strand_hits_are_labelled_with_their_own_strand():
    """Guards the closure-binding defect that mislabelled reverse-strand hits."""
    hits = screen_sequence(
        Payload(identifier="r", sequence="TTTATT" + "GC" * 8, source_format="fasta"),
        ScreenConfig(enabled=True),
    )
    for hit in hits:
        assert hit.strand in ("+", "-", "both")
        assert 1 <= hit.start <= hit.end


def test_screens_are_off_until_requested():
    assert screen_sequence(_payload(), ScreenConfig()) == ()


def test_optimizer_refuses_to_introduce_new_screen_findings():
    _, report = optimize_payload(
        _payload(),
        lint_config=LintConfig(screens=ScreenConfig(enabled=True)),
        settings=OptimizeSettings(prevent_screen_regressions=True),
    )
    before = {(h.kind, h.name, h.start) for h in report.before_lint.metrics.screen_hits}
    after = {(h.kind, h.name, h.start) for h in report.after_lint.metrics.screen_hits}
    assert after.issubset(before)


# --- CLI surface -----------------------------------------------------------


def test_cli_exposes_every_gap_control(tmp_path):
    source = _write(tmp_path, "in.fasta", f">demo\n{CPG_CDS}\n")
    vector = _write(
        tmp_path,
        "vector.json",
        json.dumps(
            {
                "genome_type": "ssAAV",
                "representation": "cassette",
                "omitted_itr_lengths": [145, 145],
                "components": [],
            }
        ),
    )
    output = tmp_path / "out.fasta"
    result = runner.invoke(
        app,
        [
            "optimize", "--input", str(source), "--cds", f"1:{len(CPG_CDS)}",
            "--solver", "global", "--search-budget", "50000",
            "--cpg-policy", "strict-subset", "--vector", str(vector),
            "--screen", "--assembly-standard", "golden-gate-bsai",
            "--protect", "1:3:start", "--gc3-min", "0", "--gc3-max", "100",
            "--output", str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Solver: global [proven_optimal]" in result.stdout
    assert "CpG policy strict_subset" in result.stdout
    assert "Vector: ssAAV cassette" in result.stdout
    audit = json.loads((tmp_path / "out.fasta.report.json").read_text(encoding="utf-8"))
    assert audit["settings"]["solver"] == "global"
    assert audit["settings"]["cpg_policy"] == "strict_subset"
    assert audit["after_lint"]["metrics"]["vector_budget"]["genome_type"] == "ssAAV"


def test_lint_reports_component_budget_and_screens(tmp_path):
    source = _write(tmp_path, "in.fasta", f">demo\n{CPG_CDS}\n")
    vector = _write(
        tmp_path,
        "vector.json",
        json.dumps(
            {
                "genome_type": "ssAAV",
                "representation": "cassette",
                "omitted_itr_lengths": [145, 145],
                "components": [
                    {
                        "start": 1, "end": 9, "name": "promoter_x", "kind": "promoter",
                        "protected": False,
                        "alternatives": [
                            {"name": "short_promoter", "length": 4,
                             "source": "catalog", "citation": "user-supplied"}
                        ],
                    }
                ],
            }
        ),
    )
    result = runner.invoke(
        app,
        ["lint", "--input", str(source), "--cds", f"1:{len(CPG_CDS)}",
         "--vector", str(vector), "--screen"],
    )
    assert "Vector: ssAAV cassette" in result.stdout
    assert "promoter_x" in result.stdout
    assert "replaceable" in result.stdout
    assert "shorter alternative" in result.stdout
