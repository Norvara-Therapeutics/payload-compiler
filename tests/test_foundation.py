import json

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError
from typer.testing import CliRunner

from payload_compiler.cli import app
from payload_compiler.exceptions import InputError
from payload_compiler.models import CodingRegion, Payload
from payload_compiler.parser import parse_cds, read_payload
from payload_compiler.validation import extract_cds, validate_payload

runner = CliRunner()


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "toy.fasta"
    path.write_text(">toy synthetic fixture\natggcttaa\n", encoding="utf-8")
    return path


def genbank(tmp_path, *, sequence="ATGGCTTAA", location=None, qualifiers=None):
    record = SeqRecord(Seq(sequence), id="toy", description="synthetic fixture")
    record.annotations = {"molecule_type": "DNA", "topology": "linear"}
    record.features = [
        SeqFeature(
            location if location is not None else SimpleLocation(0, len(sequence), strand=1),
            type="CDS",
            qualifiers=qualifiers if qualifiers is not None else {"translation": ["MA"]},
        )
    ]
    path = tmp_path / "toy.gb"
    SeqIO.write(record, path, "genbank")
    return path


def test_fasta_normalization_coordinates_and_no_mutation(fasta):
    before = fasta.read_bytes()
    payload = read_payload(fasta, cds=("1:9",))
    assert payload.sequence == "ATGGCTTAA"
    assert (payload.coding_regions[0].start, payload.coding_regions[0].end) == (0, 9)
    assert validate_payload(payload).valid
    assert fasta.read_bytes() == before


@pytest.mark.parametrize("spec", ["0:9", "9:1", "1:9:x", "1-9", "-1:9", "1:9:1"])
def test_invalid_coordinate_syntax(spec):
    with pytest.raises(InputError):
        parse_cds(spec)


def test_reverse_strand_and_flanking_bases(tmp_path):
    path = tmp_path / "reverse.fa"
    path.write_text(">reverse\nCCTTAAGCCATGG\n")
    payload = read_payload(path, cds=("3:11:-",))
    assert extract_cds(payload, payload.coding_regions[0]) == "ATGGCTTAA"
    assert validate_payload(payload).valid


@pytest.mark.parametrize(
    ("sequence", "spec", "code"),
    [
        ("ATGGCTTAA", "1:12", "CDS_OUT_OF_BOUNDS"),
        ("ATGGCTTAA", "1:8", "CDS_FRAME"),
        ("ATGTAGGCT", "1:9", "CDS_INTERNAL_STOP"),
        ("ATGNCGTAA", "1:9", "CDS_UNRESOLVED_BASE"),
        ("ATGUCTTAA", "1:9", "INVALID_BASE"),
    ],
)
def test_sequence_validation_failures(sequence, spec, code):
    payload = Payload(
        identifier="test",
        sequence=sequence,
        source_format="fasta",
        coding_regions=(parse_cds(spec),),
    )
    report = validate_payload(payload)
    assert not report.valid
    assert code in {finding.code for finding in report.findings}


def test_ambiguous_flank_and_missing_cds_policy():
    payload = Payload(identifier="test", sequence="NATGGCTTAA", source_format="fasta")
    report = validate_payload(payload)
    assert report.valid
    assert report.findings[0].positions == (1,)
    assert not validate_payload(payload, require_cds=True).valid
    annotated = payload.model_copy(update={"coding_regions": (parse_cds("2:10"),)})
    assert validate_payload(annotated, require_cds=True).valid


def test_overlap_is_reported(fasta):
    payload = read_payload(fasta, cds=("1:9", "4:9"))
    report = validate_payload(payload)
    assert report.valid
    assert [f.code for f in report.findings] == ["OVERLAPPING_CDS"]


def test_genbank_metadata_and_translation(tmp_path):
    payload = read_payload(genbank(tmp_path))
    assert payload.features[0].qualifiers["translation"] == ("MA",)
    assert payload.coding_regions[0].annotated_translation == "MA"
    assert validate_payload(payload).valid


def test_genbank_reverse_strand(tmp_path):
    payload = read_payload(
        genbank(tmp_path, sequence="TTAAGCCAT", location=SimpleLocation(0, 9, strand=-1))
    )
    assert extract_cds(payload, payload.coding_regions[0]) == "ATGGCTTAA"
    assert validate_payload(payload).valid


def test_translation_mismatch(tmp_path):
    payload = read_payload(genbank(tmp_path, qualifiers={"translation": ["MM"]}))
    report = validate_payload(payload)
    assert not report.valid
    assert report.findings[0].code == "TRANSLATION_MISMATCH"


def test_genbank_cds_uses_biological_identifier(tmp_path):
    payload = read_payload(
        genbank(
            tmp_path,
            qualifiers={"locus_tag": ["AAV_PAYLOAD"], "translation": ["MA"]},
        )
    )

    assert payload.coding_regions[0].name == "AAV_PAYLOAD"


def test_nonstandard_table_and_alternative_start(tmp_path):
    path = genbank(
        tmp_path, sequence="GTGGCTTAA", qualifiers={"transl_table": ["11"], "translation": ["MA"]}
    )
    assert validate_payload(read_payload(path)).valid


@pytest.mark.parametrize(
    "qualifiers",
    [
        {"codon_start": ["2"]},
        {"pseudo": [""]},
        {"transl_except": ["(pos:1..3,aa:Sec)"]},
        {"transl_table": ["999"]},
        {"transl_table": ["bad"]},
        {"translation": ["MA", "MM"]},
    ],
)
def test_unsupported_genbank_qualifiers(tmp_path, qualifiers):
    with pytest.raises(InputError):
        read_payload(genbank(tmp_path, qualifiers=qualifiers))


@pytest.mark.parametrize(
    "location",
    [
        SimpleLocation(BeforePosition(0), 9, strand=1),
        SimpleLocation(0, 9, strand=1, ref="REMOTE"),
    ],
)
def test_unsupported_locations(tmp_path, location):
    with pytest.raises(InputError):
        read_payload(genbank(tmp_path, location=location))


def test_exact_compound_cds_is_supported(tmp_path):
    location = CompoundLocation(
        [SimpleLocation(0, 3, strand=1), SimpleLocation(6, 9, strand=1)]
    )
    payload = read_payload(genbank(tmp_path, location=location, qualifiers={"translation": ["M*"]}))

    assert payload.coding_regions[0].parts == ((0, 3), (6, 9))
    assert extract_cds(payload, payload.coding_regions[0]) == "ATGTAA"
    assert validate_payload(payload).valid


def test_genbank_cannot_override_cds(tmp_path):
    with pytest.raises(InputError, match="FASTA only"):
        read_payload(genbank(tmp_path), cds=("1:9",))


@pytest.mark.parametrize("contents", ["", ">empty\n", ">\nATG\n", "ATG\n", ">a\nATG\n>b\nATG\n"])
def test_bad_fasta_inputs(tmp_path, contents):
    path = tmp_path / "bad.fa"
    path.write_text(contents)
    with pytest.raises(InputError):
        read_payload(path)


def test_unknown_extension_override(fasta):
    path = fasta.rename(fasta.with_suffix(".txt"))
    with pytest.raises(InputError, match="infer input format"):
        read_payload(path)
    assert read_payload(path, input_format="fasta").identifier == "toy"


def test_genbank_length_mismatch(tmp_path):
    path = genbank(tmp_path)
    path.write_text(path.read_text().replace("9 bp", "8 bp"))
    with pytest.raises(InputError):
        read_payload(path)


def test_model_constraints():
    with pytest.raises(ValidationError):
        CodingRegion(name="bad", start=9, end=3)
    with pytest.raises(ValidationError):
        CodingRegion(name="bad", start=0, end=3, strand=0)


def test_cli_json_determinism_and_input_preservation(fasta):
    before = fasta.read_bytes()
    args = ["validate", "--input", str(fasta), "--cds", "1:9", "--json"]
    first, second = runner.invoke(app, args), runner.invoke(app, args)
    assert first.exit_code == second.exit_code == 0
    assert first.stdout == second.stdout
    report = json.loads(first.stdout)
    assert report["valid"] and report["cds_count"] == 1
    assert report["length"] == 9
    assert len(report["sequence_sha256"]) == 64
    assert str(fasta.parent) not in first.stdout
    assert fasta.read_bytes() == before


def test_cli_validation_error(fasta):
    result = runner.invoke(app, ["validate", "--input", str(fasta), "--cds", "1:8", "--json"])
    assert result.exit_code == 1
    assert not json.loads(result.stdout)["valid"]


def test_cli_missing_file(tmp_path):
    result = runner.invoke(app, ["validate", "--input", str(tmp_path / "missing.fa"), "--json"])
    assert result.exit_code == 2
    assert json.loads(result.stdout)["error"]["code"] == "INPUT_ERROR"
    assert "Traceback" not in result.output


def test_cli_require_cds(fasta):
    result = runner.invoke(app, ["validate", "--input", str(fasta), "--require-cds"])
    assert result.exit_code == 1
    assert "NO_CDS" in result.stdout


@pytest.mark.parametrize("args", [["--help"], ["validate", "--help"], ["--version"]])
def test_cli_discovery(args):
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output


def test_context_dependent_translation_is_error():
    payload = Payload(
        identifier="test",
        sequence="ATGGCTTAA",
        source_format="fasta",
        coding_regions=(CodingRegion(name="test", start=0, end=9, translation_table=27),),
    )
    assert not validate_payload(payload).valid


def test_extract_rejects_out_of_bounds(fasta):
    with pytest.raises(ValueError):
        extract_cds(read_payload(fasta), parse_cds("1:12"))
