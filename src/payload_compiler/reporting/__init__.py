"""Render deterministic, standalone reports."""

from payload_compiler.reporting.diff import render_optimize_html
from payload_compiler.reporting.html import render_lint_html

__all__ = ["render_lint_html", "render_optimize_html"]
