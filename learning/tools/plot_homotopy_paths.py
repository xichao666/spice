"""将固定点 Homotopy 实验导出的 CSV 绘制为“未知量—lambda”路径图。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


# 输出为白底科学图；SVG 保留可编辑文字，PNG 用于快速查看。
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "font.size": 9,
})


def read_csv(path: Path) -> list[dict[str, float]]:
    """读取 C 程序输出的数值 CSV，不改变或筛除任何路径点。"""
    with path.open(newline="", encoding="utf-8") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def draw_case(data_directory: Path, output_directory: Path, case_name: str) -> None:
    """每个未知量占一个子图，避免不同量纲掩盖微小的分支变化。"""
    path_rows = read_csv(data_directory / f"{case_name}_path.csv")
    solution_rows = read_csv(data_directory / f"{case_name}_solutions.csv")
    variables = [key for key in path_rows[0] if key.startswith("x")]
    columns = 2 if len(variables) <= 8 else 3
    rows = (len(variables) + columns - 1) // columns
    figure, axes = plt.subplots(rows, columns, figsize=(4.7 * columns, 2.45 * rows),
                                sharex=True, constrained_layout=True)
    axes_list = list(axes.flat) if hasattr(axes, "flat") else [axes]

    lambdas = [row["lambda"] for row in path_rows]
    for index, variable in enumerate(variables):
        axis = axes_list[index]
        values = [row[variable] for row in path_rows]
        axis.plot(lambdas, values, color="#2166ac", linewidth=1.15)
        axis.axvline(1.0, color="#b2182b", linestyle="--", linewidth=0.9)
        # 工作点由原始 F(x)=0 的 Newton 复核而来，故以实心点单独标识。
        axis.scatter([1.0] * len(solution_rows), [row[variable] for row in solution_rows],
                     color="#b2182b", s=16, zorder=3,
                     label="DC operating points" if index == 0 else None)
        axis.set_title(variable, loc="left", fontweight="bold")
        axis.set_ylabel("value")
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.55)
        axis.margins(x=0.02)
    for axis in axes_list[len(variables):]:
        axis.remove()
    for axis in axes_list[:len(variables)]:
        axis.set_xlabel("lambda")
    axes_list[0].legend(loc="best")
    figure.suptitle(f"{case_name}: fixed-point Homotopy paths", fontweight="bold", y=1.01)
    output_directory.mkdir(parents=True, exist_ok=True)
    base = output_directory / f"{case_name.lower()}_homotopy_paths"
    figure.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    figure.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="绘制 Homotopy 路径图")
    parser.add_argument("data_directory", type=Path, help="C 实验程序输出的 CSV 目录")
    parser.add_argument("output_directory", type=Path, help="图片输出目录")
    arguments = parser.parse_args()
    for case_name in ("Schmitt1", "Schmitt2", "Chua"):
        draw_case(arguments.data_directory, arguments.output_directory, case_name)


if __name__ == "__main__":
    main()
