"""Dependency-free HTML rendering for lint reports."""

# The embedded stylesheet and markup intentionally retain long lines.
# ruff: noqa: E501

from html import escape

from payload_compiler.models import Finding, LintReport


def _finding_rows(findings: tuple[Finding, ...]) -> str:
    if not findings:
        return '<tr><td colspan="4" class="muted">No findings</td></tr>'
    rows = []
    for finding in findings:
        positions = ", ".join(map(str, finding.positions)) or "—"
        rows.append(
            f'<tr><td><span class="pill {finding.severity}">{finding.severity}</span></td>'
            f"<td><code>{escape(finding.code)}</code></td>"
            f"<td>{escape(finding.message)}</td><td>{escape(positions)}</td></tr>"
        )
    return "".join(rows)


def _motif_rows(report: LintReport) -> str:
    if not report.metrics.motif_matches:
        return '<tr><td colspan="5" class="muted">No forbidden motifs found</td></tr>'
    return "".join(
        f"<tr><td>{escape(match.name)}</td><td><code>{escape(match.sequence)}</code></td>"
        f"<td>{match.start}</td><td>{match.end}</td><td>{match.strand}</td></tr>"
        for match in report.metrics.motif_matches
    )


def _splice_rows(report: LintReport) -> str:
    signals = report.metrics.splice_signals
    if not signals:
        message = (
            "No consensus sites found" if report.configuration.splice_scan else "Not requested"
        )
        return f'<tr><td colspan="5" class="muted">{message}</td></tr>'
    return "".join(
        f"<tr><td>{signal.kind.title()}</td><td><code>{escape(signal.sequence)}</code></td>"
        f"<td>{signal.start}</td><td>{signal.end}</td><td>{signal.strand}</td></tr>"
        for signal in signals
    )


def _folding_rows(report: LintReport) -> str:
    metrics = report.metrics
    if not metrics.folding_windows:
        message = (
            "No resolved windows available" if report.configuration.rnafold else "Not requested"
        )
        return f'<tr><td colspan="6" class="muted">{message}</td></tr>'
    return "".join(
        f"<tr><td>{window.start}</td><td>{window.end}</td><td>{window.gc_percent:g}%</td>"
        f"<td>{window.mfe_kcal_per_mol:g}</td><td>{'yes' if window.flagged else 'no'}</td>"
        f"<td><code>{escape(window.dot_bracket)}</code></td></tr>"
        for window in metrics.folding_windows
    )


def _window_rows(report: LintReport) -> str:
    rows = []
    for window in report.metrics.gc_windows:
        if window.gc_percent is None:
            continue
        if window.gc_percent < report.configuration.gc_min_percent:
            rows.append(("Low GC", window.start, window.end, f"{window.gc_percent:g}%"))
        elif window.gc_percent > report.configuration.gc_max_percent:
            rows.append(("High GC", window.start, window.end, f"{window.gc_percent:g}%"))
    rows.extend(
        ("CpG hotspot", window.start, window.end, str(window.cpg_count))
        for window in report.metrics.cpg_hotspots
    )
    if not rows:
        return '<tr><td colspan="4" class="muted">No hotspot windows found</td></tr>'
    return "".join(
        f"<tr><td>{kind}</td><td>{start}</td><td>{end}</td><td>{value}</td></tr>"
        for kind, start, end, value in rows
    )


def _overview_svg(report: LintReport) -> str:
    length = report.metrics.length
    cpg_marks = "".join(
        f'<line x1="{20 + 860 * (position - 1) / max(1, length - 1):.2f}" y1="39" '
        f'x2="{20 + 860 * (position - 1) / max(1, length - 1):.2f}" y2="61" />'
        for position in report.metrics.cpg_positions
    )
    motif_marks = "".join(
        f'<rect x="{20 + 860 * (match.start - 1) / length:.2f}" y="76" '
        f'width="{max(2, 860 * (match.end - match.start + 1) / length):.2f}" height="12" />'
        for match in report.metrics.motif_matches
    )
    gc_points = []
    for window in report.metrics.gc_windows:
        if window.gc_percent is None:
            continue
        midpoint = ((window.start - 1) + window.end) / 2
        x = 20 + 860 * midpoint / length
        y = 130 - window.gc_percent * 0.35
        gc_points.append(f"{x:.2f},{y:.2f}")
    if len(gc_points) == 1:
        x, y = gc_points[0].split(",")
        gc_plot = f'<circle cx="{x}" cy="{y}" r="3" />'
    else:
        gc_plot = f'<polyline points="{" ".join(gc_points)}" />' if gc_points else ""
    return f"""
    <svg class="overview" viewBox="0 0 900 155" role="img"
         aria-label="Sequence overview showing CpG sites, motif matches, and local GC percentage">
      <text x="20" y="22">1</text><text x="880" y="22" text-anchor="end">{length} nt</text>
      <line class="axis" x1="20" y1="32" x2="880" y2="32" />
      <text x="20" y="57">CpG</text><g class="cpg">{cpg_marks}</g>
      <text x="20" y="86">Motifs</text><g class="motif">{motif_marks}</g>
      <text x="20" y="120">GC%</text><line class="axis" x1="20" y1="130" x2="880" y2="130" />
      <g class="gc">{gc_plot}</g>
    </svg>"""


def render_lint_html(report: LintReport) -> str:
    """Return a self-contained report with no timestamps or external resources."""
    metrics = report.metrics
    all_findings = report.validation.findings + report.findings
    gc = "unavailable" if metrics.gc_percent is None else f"{metrics.gc_percent:g}%"
    cpg_oe = (
        "unavailable"
        if metrics.cpg_observed_expected is None
        else f"{metrics.cpg_observed_expected:g}"
    )
    status = "PASS" if report.valid else "FAIL"
    status_class = "pass" if report.valid else "fail"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>payload-compiler lint report — {escape(report.identifier)}</title>
<style>
:root{{--ink:#17202a;--muted:#667085;--line:#d9e0e7;--panel:#f7f9fb;--accent:#3455db;
--good:#18794e;--bad:#b42318;--warn:#9a6700}}*{{box-sizing:border-box}}body{{margin:0;
font:15px/1.5 system-ui,sans-serif;color:var(--ink);background:#eef2f6}}main{{max-width:1040px;
margin:32px auto;padding:32px;background:white;border:1px solid var(--line);border-radius:14px}}
h1{{font-size:28px;margin:0}}h2{{font-size:18px;margin:30px 0 10px}}.sub{{color:var(--muted)}}
.status{{display:inline-block;margin-left:10px;padding:3px 9px;border-radius:999px;font-weight:700}}
.pass{{background:#e7f6ed;color:var(--good)}}.fail{{background:#feeceb;color:var(--bad)}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:24px 0}}.card{{padding:16px;
background:var(--panel);border:1px solid var(--line);border-radius:10px}}.card b{{display:block;
font-size:24px}}table{{width:100%;border-collapse:collapse}}th,td{{padding:9px 10px;text-align:left;
border-bottom:1px solid var(--line);vertical-align:top}}th{{font-size:12px;text-transform:uppercase;
color:var(--muted)}}code{{font-size:13px}}.pill{{padding:2px 7px;border-radius:999px;font-size:12px}}
.pill.error{{background:#feeceb;color:var(--bad)}}.pill.warning{{background:#fff2cc;color:var(--warn)}}
.muted{{color:var(--muted)}}.overview{{width:100%;height:auto;background:var(--panel);border-radius:10px}}
.overview text{{font-size:11px;fill:var(--muted)}}.axis{{stroke:var(--line)}}.cpg line{{stroke:var(--accent);
stroke-width:1.5}}.motif rect{{fill:var(--bad);opacity:.75}}.gc polyline{{fill:none;stroke:var(--good);
stroke-width:2}}.gc circle{{fill:var(--good)}}footer{{margin-top:28px;color:var(--muted);font-size:12px;word-break:break-all}}
@media(max-width:700px){{main{{margin:0;padding:20px;border-radius:0}}.cards{{grid-template-columns:1fr 1fr}}}}
</style></head><body><main>
<header><h1>{escape(report.identifier)} <span class="status {status_class}">{status}</span></h1>
<div class="sub">Read-only payload lint report</div></header>
<section class="cards">
 <div class="card"><span>Length</span><b>{metrics.length} nt</b><small>limit {report.configuration.max_length}</small></div>
 <div class="card"><span>CpG sites</span><b>{metrics.cpg_count}</b><small>O/E {cpg_oe}</small></div>
 <div class="card"><span>Overall GC</span><b>{gc}</b><small>{len(metrics.gc_windows)} windows</small></div>
 <div class="card"><span>Forbidden motifs</span><b>{len(metrics.motif_matches)}</b><small>{len(report.configuration.forbidden_motifs)} rules</small></div>
</section>
<h2>Sequence overview</h2>{_overview_svg(report)}
<h2>Findings</h2><table><thead><tr><th>Level</th><th>Code</th><th>Finding</th><th>Positions</th></tr></thead>
<tbody>{_finding_rows(all_findings)}</tbody></table>
<h2>Hotspot windows</h2><table><thead><tr><th>Type</th><th>Start</th><th>End</th><th>Value</th></tr></thead>
<tbody>{_window_rows(report)}</tbody></table>
<h2>Forbidden motif matches</h2><table><thead><tr><th>Name</th><th>Pattern</th><th>Start</th><th>End</th><th>Strand</th></tr></thead>
<tbody>{_motif_rows(report)}</tbody></table>
<h2>Splice consensus scan</h2><table><thead><tr><th>Type</th><th>Sequence</th><th>Start</th><th>End</th><th>Strand</th></tr></thead>
<tbody>{_splice_rows(report)}</tbody></table>
<h2>Local RNAfold windows</h2><p class="sub">Tool: {escape(metrics.folding_tool_version or "not requested")} · 
{metrics.folding_windows_skipped} ambiguous window(s) skipped</p>
<table><thead><tr><th>Start</th><th>End</th><th>GC</th><th>MFE (kcal/mol)</th><th>Flagged</th><th>Dot-bracket structure</th></tr></thead>
<tbody>{_folding_rows(report)}</tbody></table>
<footer>payload-compiler {escape(report.tool_version)} · schema {escape(report.schema_version)} · SHA-256 {escape(report.sequence_sha256)}</footer>
</main></body></html>\n"""
