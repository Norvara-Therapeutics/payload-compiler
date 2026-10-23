"""Parse exactly one DNA record without guessing its reading frame."""

import re
import warnings
from pathlib import Path
from typing import Literal

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqFeature import CompoundLocation, ExactPosition, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError

from payload_compiler.exceptions import InputError
from payload_compiler.models import CodingRegion, FeatureAnnotation, Payload

InputFormat = Literal["fasta", "genbank", "embl"]
SUFFIXES: dict[str, InputFormat] = {
    ".fa": "fasta",
    ".fna": "fasta",
    ".fasta": "fasta",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".embl": "embl",
    ".emb": "embl",
}


def parse_cds(spec: str, index: int = 1) -> CodingRegion:
    """Convert one-based inclusive START:END[:+|-] into internal coordinates."""
    match = re.fullmatch(r"([1-9][0-9]*):([1-9][0-9]*)(?::([+-]))?", spec)
    if match is None:
        raise InputError(f"Invalid CDS '{spec}'; expected START:END[:+|-], one-based inclusive")
    start, end = int(match[1]), int(match[2])
    if end < start:
        raise InputError("CDS end must be at least its start, including on the reverse strand")
    return CodingRegion(
        name=f"cds_{index}", start=start - 1, end=end, strand=-1 if match[3] == "-" else 1
    )


def _annotated_cds(record: SeqRecord) -> tuple[CodingRegion, ...]:
    regions = []
    used_labels: set[str] = set()
    for feature in record.features:
        if feature.type != "CDS":
            continue
        qualifiers = feature.qualifiers
        preferred_label = next(
            (
                str(qualifiers[key][0]).strip()
                for key in ("locus_tag", "gene", "protein_id", "label")
                if qualifiers.get(key) and str(qualifiers[key][0]).strip()
            ),
            f"cds_{len(regions) + 1}",
        )
        label = preferred_label
        suffix = 2
        while label in used_labels:
            label = f"{preferred_label}_{suffix}"
            suffix += 1
        used_labels.add(label)
        loc = feature.location
        if not isinstance(loc, (SimpleLocation, CompoundLocation)):
            raise InputError(f"{label}: missing CDS location")
        parts = tuple(loc.parts) if isinstance(loc, CompoundLocation) else (loc,)
        if any(
            part.ref
            or part.ref_db
            or not all(type(bound) is ExactPosition for bound in (part.start, part.end))
            for part in parts
        ):
            raise InputError(f"{label}: remote or fuzzy CDS locations are not supported yet")
        strands = {part.strand for part in parts}
        if len(strands) != 1 or next(iter(strands)) not in (1, -1):
            raise InputError(f"{label}: CDS strand must be explicit")
        strand = next(iter(strands))
        for key in ("codon_start", "transl_table", "translation"):
            if key in qualifiers and len(qualifiers[key]) != 1:
                raise InputError(f"{label}: expected exactly one /{key} value")
        if qualifiers.get("codon_start", ["1"]) != ["1"]:
            raise InputError(f"{label}: partial CDS with codon_start != 1 is not supported yet")
        if any(key in qualifiers for key in ("transl_except", "exception", "pseudo", "pseudogene")):
            raise InputError(f"{label}: exceptional or pseudo CDS translation is not supported yet")
        try:
            regions.append(
                CodingRegion(
                    name=label,
                    start=min(int(part.start) for part in parts),
                    end=max(int(part.end) for part in parts),
                    strand=strand,
                    translation_table=int(qualifiers.get("transl_table", ["1"])[0]),
                    annotated_translation=(
                        "".join(qualifiers["translation"][0].split())
                        if "translation" in qualifiers
                        else None
                    ),
                    parts=(
                        tuple((int(part.start), int(part.end)) for part in parts)
                        if len(parts) > 1
                        else ()
                    ),
                )
            )
        except (ValueError, ValidationError) as exc:
            raise InputError(f"{label}: invalid CDS annotation: {exc}") from exc
    return tuple(regions)


def read_payload(
    path: Path, *, input_format: InputFormat | None = None, cds: tuple[str, ...] = ()
) -> Payload:
    """Read one FASTA, GenBank, or EMBL record. Explicit CDSs are FASTA-only."""
    fmt = input_format or SUFFIXES.get(path.suffix.lower())
    if fmt not in ("fasta", "genbank", "embl"):
        raise InputError("Cannot infer input format; use --format fasta, genbank, or embl")
    if fmt != "fasta" and cds:
        raise InputError("--cds is for FASTA only; GenBank and EMBL use annotated CDS features")
    try:
        with path.open(encoding="utf-8-sig") as handle:
            # Reject preambles that a permissive FASTA parser could silently discard.
            first = next((line.strip() for line in handle if line.strip()), "")
            if not first:
                raise InputError("Input file is empty")
            if fmt == "fasta" and not first.startswith(">"):
                raise InputError("FASTA must begin with a '>' header")
            if fmt == "genbank" and not first.startswith("LOCUS "):
                raise InputError("GenBank must begin with a LOCUS line")
            if fmt == "embl" and not first.startswith("ID "):
                raise InputError("EMBL must begin with an ID line")
            handle.seek(0)
            with warnings.catch_warnings():
                warnings.simplefilter("error", BiopythonParserWarning)
                record = SeqIO.read(handle, fmt)
        if not record.id or record.id == "<unknown id>":
            raise InputError("Sequence record must have a nonempty identifier")
        sequence = str(record.seq).upper()
        if not sequence:
            raise InputError("Sequence record is empty")
        molecule = record.annotations.get("molecule_type") or "DNA"
        if "DNA" not in molecule.upper():
            raise InputError(f"Expected DNA, found molecule type '{molecule}'")
        return Payload(
            identifier=record.id,
            description=record.description,
            sequence=sequence,
            source_format=fmt,
            topology=record.annotations.get("topology", "linear"),
            coding_regions=(
                _annotated_cds(record)
                if fmt != "fasta"
                else tuple(parse_cds(spec, i) for i, spec in enumerate(cds, 1))
            ),
            features=tuple(
                FeatureAnnotation(
                    kind=f.type,
                    location=str(f.location) if f.location is not None else None,
                    qualifiers={key: tuple(values) for key, values in f.qualifiers.items()},
                )
                for f in record.features
            ),
        )
    except InputError:
        raise
    except (OSError, UnicodeError, ValueError, BiopythonParserWarning) as exc:
        raise InputError(f"Cannot read {path.name}: {exc}") from exc
