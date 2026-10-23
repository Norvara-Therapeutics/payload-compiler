# payload-compiler

CpG refactoring, payload linting and AAV vector size accounting for DNA payloads
from the command line. Every edit is a synonymous codon substitution that
preserves each annotated translation exactly, recorded in a reviewable audit
trail and a self-contained HTML diff that opens from local disk.

## Install

Install the prerelease from PyPI:

```bash
python -m pip install --pre norvara-payload-compiler
```

The Python package and command keep the shorter project name:

```python
import payload_compiler
```

```bash
payload-compiler --version
```

## Quickstart

Deplete the CpG sites in phiX174 gene D (annotated at 390..848) without changing
the protein it encodes:

```bash
payload-compiler optimize --input examples/phix174.fasta --cds 390:848 \
  --max-length 5386 --output geneD.fasta
```

```
PASS NC_001422.1: CpG 267 → 231 (36 removed), 35 mutation(s)
Solver: greedy [locally_optimal] over 36 state(s); CpG policy strict_subset
Profile: human Homo sapiens CUTG, GenBank Release 160.0, 2007-06-15
  CDS cds_1: CAI 0.712, GC3 25.7%, 41 distinct codons, 2 rare, pair-bias n/a
WARNING 231 CpG site(s) remain outside editable CDS contexts or were retained under the selected constraints and search budget
WARNING No vector architecture supplied: length was checked against a flat 5386-nt limit with no ITR or component accounting
Sequence (fasta): geneD.fasta
Audit: geneD.fasta.report.json
Visual diff: geneD.fasta.diff.html
```

All 36 CpGs inside the CDS are gone, the 231 that remain lie outside it, and
nothing outside 390..848 is touched. The HTML diff aligns both sequences by
coordinate and lists every codon substitution.

## Commands

| Command | Purpose |
| --- | --- |
| `validate` | parse and check one record; report CDS findings |
| `lint` | read-only metrics, motif, splice, structure and screen findings |
| `optimize` | rewrite coding regions, with audit and visual diff |
| `report` | rebuild the HTML diff from an audit and verify it replays |
| `profiles` | list the bundled codon-use profiles |

## Options

```
-i, --input PATH         single-record FASTA, GenBank or EMBL
-o, --output PATH        optimized sequence output
    --cds START:END[:+|-]  FASTA reading frame, 1-based inclusive; repeatable
    --max-length N       payload length limit (default: 4700)
    --host NAME          bundled codon profile (default: human)
    --profile FILE       user-supplied versioned codon profile
    --solver {greedy,global}   greedy edits, or bounded exhaustive search
    --search-budget N    assignments the global solver may attempt
    --cpg-policy {strict-subset,count-only}
    --vector FILE        vector architecture for ITR-aware size accounting
    --protect START:END[:NAME]   intervals that must never change
    --forbid-motif NAME=IUPAC    motif never to introduce; repeatable
    --splice-scan        screen splice-consensus motifs
    --screen             enable the heuristic sequence screens
    --assembly-standard {golden-gate-bsai,golden-gate-bsmbi,biobrick}
    --rnafold {off,auto,seqfold,required}   folding backend
                         (on `lint` this is a plain on/off flag)
    --gc3-min N / --gc3-max N / --max-rare-codons N
```

`payload-compiler optimize --help` lists the complete set.

## Outputs

| Flag | Output |
| --- | --- |
| `--output PATH` | the optimized sequence; FASTA unless the suffix or `--output-format` selects the annotated source format |
| `--json-report PATH` | the JSON mutation audit (default `<output>.report.json`) |
| `--diff-report PATH` | the self-contained HTML diff (default `<output>.diff.html`) |
| `--json` | the audit on stdout |

Outputs are written as one rollback-capable set: if the run is unsafe, no stale
sequence is left at the output path. Repeating a run with identical inputs,
settings and software versions reproduces every byte.

## Reading the output

- Coordinates are 1-based and inclusive. `--cds 390:848` names the forward-strand
  frame; append `:-` for the reverse strand.
- Every accepted substitution is synonymous under the CDS translation table, and
  translation preservation is re-checked before anything is written. Start and
  stop codon identity is never altered.
- `--cpg-policy strict-subset` (the default) requires that the optimized sequence
  carry CpGs only at coordinates that already had one, so a CpG is never
  relocated. `count-only` relaxes this to a falling total and lists any new
  coordinates in `introduced_cpg_positions`.
- The solver line states what was established: `locally_optimal` means no
  improving single-codon edit remains; `proven_optimal` means the declared
  constraint space was exhausted; `best_found` means `--search-budget` was
  reached first. Proof is over the declared constraints only, never a claim of
  biological optimality.
- Codon-frequency profiles break ties and record provenance; `--profile` carries
  a name, version, source and citation, and its SHA-256 appears in every audit.
  A profile is not an expression model.
- `report` independently replays each base edit and recomputes hashes, CpG
  positions, codons and translations, so a mismatched sequence and audit cannot
  be presented as one result.

## Native ViennaRNA

Folding works out of the box through the bundled seqfold backend, which reports
a minimum free energy only. Installing ViennaRNA so that `RNAfold` is on `PATH`
additionally enables `--rnafold required` and the ensemble statistics behind
`--structure-ensemble`. The `Dockerfile` pins ViennaRNA 2.5.1 if you would
rather not install it on the host:

```bash
docker build -t payload-compiler .
docker run --rm -v "$PWD:/data" -w /data payload-compiler \
  payload-compiler lint --input examples/phix174.fasta --rnafold --structure-ensemble
```

With no command, the image runs the test suite against that pinned RNAfold,
which is how the native folding path is checked.

## Limitations

- Synonymous substitution preserves length, so an oversized payload cannot be
  shortened. With `--vector` the size problem is reported per component, with the
  shortest declared swap plan; the substitution itself is not performed.
- Vector component identities and lengths are user-declared. ITR sequence
  integrity is not verified and packaging limits are screening thresholds, not
  yield predictions.
- Splice, poly(A), promoter, Kozak, ORF, homopolymer and repeat screens are motif
  and threshold heuristics. They flag candidates for review and predict neither
  splicing, polyadenylation, promoter activity nor expression.
- Secondary structure is scored on local windows, which does not model the
  complete transcript.
- Sequences are treated as linear, and one record is processed per invocation.

## Development

```bash
pip install -e ".[dev]"
pytest -q
ruff check src tests
```

## License

Licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0).
See [LICENSE](LICENSE) for the complete terms.
