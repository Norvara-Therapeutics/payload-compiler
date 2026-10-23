"""Render a deterministic, standalone visual mutation report."""

# The embedded stylesheet and markup intentionally retain long lines.
# ruff: noqa: E501

from html import escape

from payload_compiler.models import OptimizeReport, Payload


def _overview(length: int, changed_positions: tuple[int, ...]) -> str:
    marks = "".join(
        f'<line x1="{20 + 960 * (position - 1) / max(1, length - 1):.2f}" y1="34" y2="61" />'
        for position in changed_positions
    )
    return f"""
    <svg class="map" viewBox="0 0 1000 82" role="img"
         aria-label="Payload overview with changed nucleotide positions">
      <text x="20" y="18">1</text><text x="980" y="18" text-anchor="end">{length} nt</text>
      <line class="axis" x1="20" y1="30" x2="980" y2="30" />
      <g class="change-markers">{marks}</g>
      <text x="20" y="76">Changed bases ({len(changed_positions)})</text>
    </svg>"""


def _sequence_rows(before: str, after: str, changed: set[int]) -> str:
    rows = []
    for offset in range(0, len(before), 60):
        end = min(offset + 60, len(before))
        tracks = []
        for label, sequence in (("Before", before), ("After", after)):
            bases = "".join(
                f'<span class="changed">{escape(base)}</span>'
                if position in changed
                else escape(base)
                for position, base in enumerate(sequence[offset:end], start=offset + 1)
            )
            tracks.append(
                f'<div class="track"><span class="track-label">{label}</span>'
                f"<code>{bases}</code></div>"
            )
        rows.append(
            f'<section class="sequence-chunk"><div class="coordinate">{offset + 1}–{end}</div>'
            + "".join(tracks)
            + "</section>"
        )
    return "".join(rows)


def _codon_rows(report: OptimizeReport) -> str:
    if not report.mutations:
        return '<tr><td colspan="9" class="muted">No silent substitutions were required</td></tr>'
    rows = []
    for mutation in report.mutations:
        changes = "; ".join(
            f"{change.position}: {change.before} → {change.after}"
            for change in mutation.base_changes
        )
        rows.append(
            "<tr>"
            f"<td>{mutation.index}</td><td>{escape(mutation.region)}</td>"
            f"<td>{mutation.start}–{mutation.end} ({'+' if mutation.strand == 1 else '−'})</td>"
            f"<td>{escape(mutation.amino_acid)}</td>"
            f"<td><code>{escape(mutation.codon_before)}</code></td>"
            f"<td><code>{escape(mutation.codon_after)}</code></td>"
            f"<td>{escape(changes)}</td>"
            f"<td>{mutation.cpg_count_before} → {mutation.cpg_count_after}</td>"
            f"<td>{mutation.host_codon_frequency_before:g} → "
            f"{mutation.host_codon_frequency_after:g}</td></tr>"
        )
    return "".join(rows)


def _warning_items(report: OptimizeReport) -> str:
    messages = [f"Safety regression: {message}" for message in report.safety_regressions]
    messages.extend(report.warnings)
    messages.extend(
        finding.message for finding in report.after_lint.findings if finding.severity == "warning"
    )
    if not messages:
        return '<p class="muted">No optimization warnings.</p>'
    return "<ul>" + "".join(f"<li>{escape(message)}</li>" for message in messages) + "</ul>"


def render_optimize_html(source: Payload, optimized: Payload, report: OptimizeReport) -> str:
    """Return a self-contained visual diff aligned to genomic coordinates."""
    if len(source.sequence) != len(optimized.sequence):
        raise ValueError("Visual diff requires input and output sequences of equal length")
    changed_positions = tuple(
        index + 1
        for index, (before, after) in enumerate(
            zip(source.sequence, optimized.sequence, strict=True)
        )
        if before != after
    )
    changed = set(changed_positions)
    status = "PASS" if report.optimized else "REVIEW"
    status_class = "pass" if report.optimized else "review"
    source_id = escape(source.identifier)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>payload-compiler mutation diff — {source_id}</title>
<style>
:root{{--ink:#17202a;--muted:#667085;--line:#d9e0e7;--panel:#f7f9fb;--accent:#3455db;--change:#b42318;--good:#18794e;--warn:#9a6700}}
*{{box-sizing:border-box}}body{{margin:0;background:#eef2f6;color:var(--ink);font:15px/1.5 system-ui,sans-serif}}
main{{max-width:1180px;margin:32px auto;padding:32px;background:white;border:1px solid var(--line);border-radius:14px}}
h1{{font-size:28px;margin:0}}h2{{font-size:19px;margin:30px 0 10px}}p{{margin:8px 0}}.sub,.muted{{color:var(--muted)}}
.status{{display:inline-block;margin-left:10px;padding:3px 9px;border-radius:999px;font-weight:700}}
.pass{{background:#e7f6ed;color:var(--good)}}.review{{background:#fff2cc;color:var(--warn)}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:24px 0}}
.card{{padding:16px;background:var(--panel);border:1px solid var(--line);border-radius:10px}}
.card b{{display:block;font-size:24px}}.map{{display:block;width:100%;height:auto;background:var(--panel);border-radius:10px}}
.map text{{font-size:11px;fill:var(--muted)}}.axis{{stroke:var(--line);stroke-width:4;stroke-linecap:round}}
.change-markers line{{stroke:var(--change);stroke-width:2}}.legend{{color:var(--change);font-weight:700}}
table{{width:100%;border-collapse:collapse}}th,td{{padding:9px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}}
th{{font-size:12px;text-transform:uppercase;color:var(--muted)}}code{{font-family:ui-monospace,monospace;font-size:13px}}
.sequence{{display:grid;gap:8px}}.sequence-chunk{{padding:10px;background:var(--panel);border:1px solid var(--line);border-radius:8px;overflow-x:auto}}
.coordinate{{font-size:11px;color:var(--muted);margin-bottom:4px}}.track{{display:grid;grid-template-columns:62px 1fr;gap:8px;white-space:pre;line-height:1.65}}
.track-label{{font-size:11px;color:var(--muted);font-weight:700}}.changed{{background:#feeceb;color:var(--change);font-weight:800;border-radius:2px}}
footer{{margin-top:28px;color:var(--muted);font-size:12px;word-break:break-all}}
@media(max-width:760px){{main{{margin:0;padding:18px;border-radius:0}}.cards{{grid-template-columns:1fr 1fr}}table{{display:block;overflow-x:auto}}}}
</style></head><body><main>
<header><h1>{source_id} <span class="status {status_class}">{status}</span></h1>
<div class="sub">Deterministic synonymous CpG-depletion visual diff · host {escape(report.host)} · profile {escape(report.host_profile)}</div></header>
<section class="cards">
 <div class="card"><span>Payload length</span><b>{len(source.sequence)} nt</b><small>limit {report.max_length} nt</small></div>
 <div class="card"><span>CpG sites</span><b>{report.cpg_count_before} → {report.cpg_count_after}</b><small>{report.cpg_reduction} removed</small></div>
 <div class="card"><span>Silent substitutions</span><b>{len(report.mutations)}</b><small>{len(changed_positions)} net base changes</small></div>
 <div class="card"><span>Protein translation</span><b>{"Preserved" if report.protein_translation_preserved else "Changed"}</b><small>{report.after_lint.valid and "Lint passed" or "Lint needs review"}</small></div>
</section>
<h2>Mutation map</h2><p class="legend">Red marks indicate changed nucleotide positions; coordinates are one-based.</p>
{_overview(len(source.sequence), changed_positions)}
<h2>Aligned sequence</h2><p class="sub">Input and optimized tracks share genomic coordinates. Changed positions are highlighted in red.</p>
<div class="sequence">{_sequence_rows(source.sequence, optimized.sequence, changed)}</div>
<h2>Mutation audit</h2><div class="sub">Each row records one synonymous codon replacement in optimization order.</div>
<table><thead><tr><th>#</th><th>CDS</th><th>Coordinates / strand</th><th>Amino acid</th><th>Codon before</th><th>Codon after</th><th>Genomic base changes</th><th>CpGs</th><th>Host codon frequency</th></tr></thead>
<tbody>{_codon_rows(report)}</tbody></table>
<h2>Review notes</h2>{_warning_items(report)}
    <footer>payload-compiler {escape(report.tool_version)} · schema {escape(report.schema_version)} · Python {escape(report.runtime.python)} · Biopython {escape(report.runtime.biopython)} · Pydantic {escape(report.runtime.pydantic)} · seqfold {escape(report.runtime.seqfold)} · Typer {escape(report.runtime.typer)} · platform {escape(report.runtime.platform)} · source SHA-256 {escape(report.source_sha256)} · optimized SHA-256 {escape(report.optimized_sha256)}</footer>
</main></body></html>\n"""
