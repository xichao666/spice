"""Create SVG comparison charts from the CSV files without external packages."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape


COLORS = {
    "direct_newton": "#2563eb",
    "source_stepping": "#059669",
    "sequential_source_stepping": "#7c3aed",
    "gmin_stepping": "#ea580c",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read one UTF-8 CSV file into dictionaries."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def svg_header(width: int, height: int) -> list[str]:
    """Start a self-contained SVG document."""
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#1f2937}.axis{stroke:#6b7280}.grid{stroke:#d1d5db}.small{font-size:12px}.title{font-size:20px;font-weight:600}.legend{font-size:13px}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
    ]


def bar_panel(parts: list[str], rows: list[dict[str, str]], key: str,
              title: str, top: int, width: int) -> None:
    """Draw one labelled bar panel for a numeric summary-column."""
    left, right, height = 95, width - 35, 250
    values = [float(row[key]) for row in rows]
    maximum = max(values) if max(values) > 0.0 else 1.0
    bar_width = (right - left) / max(len(rows), 1) * 0.68

    parts.append(f'<text class="title" x="{left}" y="{top - 22}">{escape(title)}</text>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top + height}" x2="{right}" y2="{top + height}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + height}"/>')
    for tick in range(5):
        value = maximum * tick / 4.0
        y = top + height - height * tick / 4.0
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>')
        parts.append(f'<text class="small" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{value:.3g}</text>')
    for index, row in enumerate(rows):
        x = left + (index + 0.5) * (right - left) / len(rows) - bar_width / 2.0
        bar_height = height * float(row[key]) / maximum
        color = COLORS.get(row["algorithm"], "#64748b")
        if row["converged"] != "1":
            color = "#dc2626"
        y = top + height - bar_height
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="{color}"/>')
        label = f'{row["circuit"].replace("_netlist", "")} / {row["algorithm"]}'
        parts.append(f'<text class="small" transform="translate({x + bar_width / 2:.1f},{top + height + 16}) rotate(45)" text-anchor="start">{escape(label)}</text>')


def write_summary_chart(rows: list[dict[str, str]], output_directory: Path) -> None:
    """Write one SVG containing total-iteration and runtime bar charts."""
    width, height = 1450, 780
    parts = svg_header(width, height)
    bar_panel(parts, rows, "total_newton_iterations", "Total Newton iterations (red = failed path)", 70, width)
    bar_panel(parts, rows, "total_runtime_ms", "Measured runtime (milliseconds)", 430, width)
    parts.append('</svg>')
    (output_directory / "dc_algorithm_comparison.svg").write_text("\n".join(parts), encoding="utf-8")


def write_step_chart(output_directory: Path) -> None:
    """Write a multi-panel SVG of Newton iterations at accepted continuation points."""
    grouped: dict[str, list[tuple[str, list[dict[str, str]]]]] = defaultdict(list)
    for step_path in output_directory.glob("*/*/dc_steps.csv"):
        rows = read_csv(step_path)
        if rows:
            grouped[step_path.parents[1].name].append((step_path.parent.name, rows))
    if not grouped:
        return

    width, panel_height = 1250, 300
    parts = svg_header(width, 70 + panel_height * len(grouped))
    for panel_index, (circuit, series) in enumerate(sorted(grouped.items())):
        left, right = 95, width - 35
        top, height = 55 + panel_index * panel_height, 180
        maximum = max(int(row["newton_iterations"]) for _, rows in series for row in rows)
        maximum = max(maximum, 1)
        parts.append(f'<text class="title" x="{left}" y="{top - 18}">{escape(circuit)}: Newton iterations at accepted points</text>')
        parts.append(f'<line class="axis" x1="{left}" y1="{top + height}" x2="{right}" y2="{top + height}"/>')
        parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + height}"/>')
        for tick in range(5):
            value, y = maximum * tick / 4.0, top + height - height * tick / 4.0
            parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>')
            parts.append(f'<text class="small" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{value:.3g}</text>')
        for series_index, (algorithm, rows) in enumerate(series):
            maximum_step = max(int(row["step"]) for row in rows)
            points = []
            for row in rows:
                x = left + (right - left) * int(row["step"]) / max(maximum_step, 1)
                y = top + height - height * int(row["newton_iterations"]) / maximum
                points.append(f"{x:.1f},{y:.1f}")
            color = COLORS.get(algorithm, "#64748b")
            parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{" ".join(points)}"/>')
            parts.append(f'<text class="legend" x="{left + 230 * series_index}" y="{top + height + 48}" fill="{color}">{escape(algorithm)}</text>')
    parts.append('</svg>')
    (output_directory / "dc_step_iterations.svg").write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    """Validate the report directory and create both SVG figures."""
    if len(sys.argv) != 2:
        print(f"Usage: {Path(sys.argv[0]).name} experiment-output-directory")
        return 1
    output_directory = Path(sys.argv[1])
    summary_path = output_directory / "dc_algorithm_comparison.csv"
    if not summary_path.is_file():
        print(f"Missing comparison CSV: {summary_path}")
        return 1
    rows = read_csv(summary_path)
    if not rows:
        print("The comparison CSV has no data rows.")
        return 1
    write_summary_chart(rows, output_directory)
    write_step_chart(output_directory)
    print(f"Figures written to: {output_directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
