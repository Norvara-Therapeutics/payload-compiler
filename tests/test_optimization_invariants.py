"""Regression and generated-case verification for synonymous optimization."""

import itertools
import os
import shutil

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from hypothesis import given, settings
from hypothesis import strategies as st

from payload_compiler.analysis import fold_local_windows, scan_splice_signals
from payload_compiler.exceptions import InputError
from payload_compiler.models import CodingRegion, LintConfig, Payload
from payload_compiler.optimize import optimize_payload
from payload_compiler.validation import extract_cds

SENSE_CODONS = tuple(sorted(CodonTable.unambiguous_dna_by_id[1].forward_table))


def _payload(
    sequence: str,
    regions: tuple[CodingRegion, ...],
    *,
    identifier: str = "invariant-fixture",
) -> Payload:
    return Payload(
        identifier=identifier,
        sequence=sequence,
        source_format="fasta",
        coding_regions=regions,
    )


def _translation(payload: Payload, region: CodingRegion) -> str:
    return str(Seq(extract_cds(payload, region)).translate(table=region.translation_table))


def test_cpg_spanning_adjacent_codons_is_removed_without_changing_translation():
    # GCC|GCT creates a CpG across the codon boundary (positions 6–7).
    source = _payload("ATGGCCGCTTAA", (CodingRegion(name="main", start=0, end=12),))

    optimized, report = optimize_payload(source)

    assert source.sequence[5:7] == "CG"
    assert optimized.sequence[5:7] != "CG"
    assert _translation(optimized, optimized.coding_regions[0]) == _translation(
        source, source.coding_regions[0]
    )
    assert report.cpg_count_after < report.cpg_count_before


def test_reverse_strand_cds_is_optimized_in_transcript_orientation():
    coding_sequence = "ATGCGTTAA"
    genomic_sequence = str(Seq(coding_sequence).reverse_complement())
    source = _payload(
        genomic_sequence,
        (CodingRegion(name="reverse", start=0, end=len(genomic_sequence), strand=-1),),
    )

    optimized, report = optimize_payload(source)

    assert optimized.sequence != source.sequence
    assert report.cpg_count_after <= report.cpg_count_before
    assert _translation(optimized, optimized.coding_regions[0]) == _translation(
        source, source.coding_regions[0]
    )


def test_overlapping_cds_annotations_all_keep_their_translation():
    sequence = "ATGCGTTAA"
    regions = (
        CodingRegion(name="primary", start=0, end=9),
        CodingRegion(name="overlap", start=0, end=9),
    )
    source = _payload(sequence, regions)

    optimized, report = optimize_payload(source)

    assert report.protein_translation_preserved
    assert [
        _translation(optimized, region) for region in optimized.coding_regions
    ] == [_translation(source, region) for region in source.coding_regions]


def test_different_frame_overlapping_cds_keep_both_translations():
    sequence = "AAAAAAAACG"
    regions = (
        CodingRegion(name="frame_0", start=0, end=9),
        CodingRegion(name="frame_1", start=1, end=10),
    )
    source = _payload(sequence, regions)

    optimized, report = optimize_payload(source)

    assert report.cpg_count_after == 0
    assert [_translation(optimized, region) for region in regions] == [
        _translation(source, region) for region in regions
    ]


def test_non_coding_flanks_and_ambiguous_bases_remain_unchanged():
    sequence = "NNNATGCGTTAAACNN"
    region = CodingRegion(name="main", start=3, end=12)
    source = _payload(sequence, (region,))

    optimized, report = optimize_payload(source)

    assert optimized.sequence[: region.start] == source.sequence[: region.start]
    assert optimized.sequence[region.end :] == source.sequence[region.end :]
    assert optimized.sequence.count("N") == source.sequence.count("N")
    assert report.cpg_count_after <= report.cpg_count_before


def test_joined_cds_is_optimized_without_changing_intronic_bases():
    sequence = "ATGNNNCGTTAA"
    region = CodingRegion(
        name="joined",
        start=0,
        end=12,
        parts=((0, 3), (6, 12)),
    )
    source = _payload(sequence, (region,))

    optimized, report = optimize_payload(source)

    assert optimized.sequence[3:6] == "NNN"
    assert _translation(optimized, region) == _translation(source, region)
    assert report.cpg_count_after < report.cpg_count_before


def test_ambiguous_coding_sequence_fails_before_any_optimization():
    source = _payload("ATGNNNTAA", (CodingRegion(name="main", start=0, end=9),))

    with pytest.raises(InputError, match="CDS_UNRESOLVED_BASE"):
        optimize_payload(source)


def test_already_optimal_payload_is_an_exact_deterministic_no_op():
    source = _payload("ATGGCTTAA", (CodingRegion(name="main", start=0, end=9),))

    first, first_report = optimize_payload(source)
    second, second_report = optimize_payload(source)

    assert first.sequence == source.sequence
    assert first_report.mutations == ()
    assert first_report.model_dump(mode="json") == second_report.model_dump(mode="json")
    assert first_report.cpg_count_after == first_report.cpg_count_before == 0


def test_optimizer_does_not_introduce_new_splice_consensus_signals():
    source = _payload(
        "CCCGGCGAG",
        (CodingRegion(name="main", start=0, end=9),),
    )
    assert not scan_splice_signals(source.sequence)

    optimized, report = optimize_payload(source, lint_config=LintConfig(splice_scan=True))

    assert not scan_splice_signals(optimized.sequence)
    assert report.safety_checks_passed
    assert report.optimized


def test_generated_synonymous_cases_preserve_translation_and_never_add_cpgs():
    """A compact deterministic property sweep over synonymous codon contexts."""
    table = CodonTable.unambiguous_dna_by_id[1]
    sense_codons = sorted(table.forward_table)
    # Cover a stable, broad selection of two-codon contexts without a runtime
    # dependency on a property-testing plugin or randomized test order.
    contexts = list(
        itertools.islice(
            (
                (left, right)
                for left, right in itertools.product(sense_codons, repeat=2)
                if "CG" in left + right
            ),
            128,
        )
    )
    for index, (left, right) in enumerate(contexts):
        coding = "ATG" + left + right + "TAA"
        source = _payload(
            coding,
            (CodingRegion(name="generated", start=0, end=len(coding)),),
            identifier=f"generated-{index}",
        )
        optimized, report = optimize_payload(source)

        assert _translation(optimized, optimized.coding_regions[0]) == _translation(
            source, source.coding_regions[0]
        )
        assert report.cpg_count_after <= report.cpg_count_before
        assert report.optimized_sha256
        assert optimize_payload(source)[0].sequence == optimized.sequence


@settings(max_examples=150, deadline=None)
@given(
    internal_codons=st.lists(st.sampled_from(SENSE_CODONS), min_size=1, max_size=12),
    strand=st.sampled_from((1, -1)),
    left_flank=st.text(alphabet="ACGTN", min_size=0, max_size=8),
    right_flank=st.text(alphabet="ACGTN", min_size=0, max_size=8),
)
def test_property_translation_noncoding_cpg_and_determinism(
    internal_codons, strand, left_flank, right_flank
):
    coding = "ATG" + "".join(internal_codons) + "TAA"
    genomic_coding = coding if strand == 1 else str(Seq(coding).reverse_complement())
    sequence = left_flank + genomic_coding + right_flank
    region = CodingRegion(
        name="generated",
        start=len(left_flank),
        end=len(left_flank) + len(genomic_coding),
        strand=strand,
    )
    source = _payload(sequence, (region,))

    first, first_report = optimize_payload(source)
    second, second_report = optimize_payload(source)

    assert _translation(first, region) == _translation(source, region)
    assert first.sequence[: region.start] == left_flank
    assert first.sequence[region.end :] == right_flank
    assert first_report.cpg_count_after <= first_report.cpg_count_before
    assert first.sequence == second.sequence
    assert first_report.model_dump(mode="json") == second_report.model_dump(mode="json")


def test_rnafold_present_path_parses_windows_and_skips_ambiguous_windows(tmp_path, monkeypatch):
    executable = tmp_path / "RNAfold"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if sys.argv[1:] == ['--version']:\n"
        "    print('RNAfold test-double 1.0')\n"
        "    raise SystemExit(0)\n"
        "records = [line.strip() for line in sys.stdin if line.strip()]\n"
        "for index in range(0, len(records), 2):\n"
        "    name, sequence = records[index:index + 2]\n"
        "    print(name)\n"
        "    print(sequence)\n"
        "    print('.' * len(sequence) + ' ( -1.20)')\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    monkeypatch.setenv(
        "PATH", os.pathsep.join((str(tmp_path), os.environ.get("PATH", "")))
    )

    version, folds, skipped = fold_local_windows(
        "ACGTNNNNACGT",
        window_size=4,
        step=4,
        mfe_threshold=-2.0,
    )

    assert version == "RNAfold test-double 1.0"
    assert [(window.start, window.end) for window in folds] == [(1, 4), (9, 12)]
    assert all(window.dot_bracket == "...." and not window.flagged for window in folds)
    assert skipped == 1


def test_newly_flagged_rnafold_window_blocks_sequence_publication(tmp_path, monkeypatch):
    executable = tmp_path / "RNAfold"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if sys.argv[1:] == ['--version']:\n"
        "    print('RNAfold safety-double 1.0')\n"
        "    raise SystemExit(0)\n"
        "records = [line.strip() for line in sys.stdin if line.strip()]\n"
        "for index in range(0, len(records), 2):\n"
        "    name, sequence = records[index:index + 2]\n"
        "    energy = -1.0 if 'CG' in sequence else -20.0\n"
        "    print(name)\n"
        "    print(sequence)\n"
        "    print('.' * len(sequence) + f' ( {energy:.2f})')\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    monkeypatch.setenv(
        "PATH", os.pathsep.join((str(tmp_path), os.environ.get("PATH", "")))
    )
    source = _payload("ATGCGTTAA", (CodingRegion(name="main", start=0, end=9),))

    optimized, report = optimize_payload(
        source,
        lint_config=LintConfig(
            rnafold_mode="required",
            structure_window=20,
            structure_step=10,
            structure_mfe_threshold=-15.0,
        ),
    )

    assert optimized.sequence != source.sequence
    assert not report.safety_checks_passed
    assert not report.optimized
    assert report.safety_regressions == (
        "optimization introduced 1 newly flagged RNAfold window(s)",
    )


def test_rnafold_required_path_fails_cleanly_when_unavailable(monkeypatch):
    from payload_compiler import analysis

    monkeypatch.setattr(analysis.shutil, "which", lambda _: None)

    with pytest.raises(InputError, match="not found"):
        fold_local_windows("ACGT", window_size=20, step=10, mfe_threshold=-15.0)


@pytest.mark.integration
def test_real_rnafold_executable_when_required_by_environment():
    required = os.environ.get("PAYLOAD_COMPILER_REQUIRE_REAL_RNAFOLD") == "1"
    if shutil.which("RNAfold") is None:
        if required:
            pytest.fail("Docker/CI environment requires the real RNAfold executable")
        pytest.skip("real RNAfold is not installed")

    version, folds, skipped = fold_local_windows(
        "GCGCGCAAAAGCGCGCAAAA",
        window_size=20,
        step=10,
        mfe_threshold=-2.0,
    )

    assert "RNAfold" in version
    assert len(folds) == 1
    assert len(folds[0].dot_bracket) == 20
    assert skipped == 0
