#!/usr/bin/env python3
"""Visualize spatial cell-type weights produced by RCTD."""

from __future__ import annotations

import argparse
import math
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


BASE_GRID_SPOTS = 50
BASE_PLOT_SIZE = 4.37
LEGEND_COLUMN_WIDTH = 1.8
LEGEND_LABEL_WRAP_WIDTH = 32
LEGEND_LABEL_INCHES_PER_CHAR = 0.055
LEGEND_LABEL_EXTRA_WIDTH = 0.75
COLORBAR_WIDTH = 0.45
COLORBAR_GAP = 0.55
VERTICAL_MARGIN = 0.73
PLOT_EDGE_PAD = 0.9
SPOT_SIDE_LENGTH = 0.82
TOP_DOMINANT_TYPES = 15
OTHER_CELL_TYPE = "Other"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot spatial cell-type weights from RCTD."
    )
    parser.add_argument("--weights-file", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--x-spots-number", required=True, type=int)
    parser.add_argument("--y-spots-number", required=True, type=int)
    parser.add_argument(
        "--orientation",
        choices=("normal", "horizontal", "vertical", "rotate"),
        default="normal",
    )
    parser.add_argument("--swap-xy", action="store_true")
    parser.add_argument("--rotate", type=int, choices=(0, 90, 180, 270), default=0)
    return parser.parse_args()


def load_weights(path: Path) -> tuple[pd.DataFrame, list[str]]:
    if not path.is_file():
        raise SystemExit(f"Weights file not found: {path}")

    table = pd.read_csv(path)
    metadata_columns = ["barcode", "x", "y"]
    missing = [column for column in metadata_columns if column not in table.columns]
    if missing:
        raise SystemExit(
            "Weights file is missing required column(s): " + ", ".join(missing)
        )
    cell_types = [column for column in table.columns if column not in metadata_columns]
    if not cell_types:
        raise SystemExit("Weights file does not contain any cell-type columns")

    numeric_columns = ["x", "y", *cell_types]
    try:
        table[numeric_columns] = table[numeric_columns].apply(
            pd.to_numeric, errors="raise"
        )
    except (TypeError, ValueError) as error:
        raise SystemExit(
            f"Coordinates and cell-type weights must be numeric: {error}"
        ) from error

    values = table[cell_types].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise SystemExit("Cell-type weights must be finite and non-negative")
    row_sums = values.sum(axis=1)
    if (row_sums <= 0).any():
        raise SystemExit("Every spot must have a positive total cell-type weight")
    table[cell_types] = values / row_sums[:, None]
    return table, cell_types


def transform_coordinates(
    frame: pd.DataFrame,
    x_spots: int,
    y_spots: int,
    orientation: str,
    swap_xy: bool,
    rotate: int,
) -> tuple[pd.DataFrame, int, int]:
    transformed = frame.copy()
    if orientation == "horizontal":
        transformed["x"] = x_spots - 1 - transformed["x"]
    elif orientation == "vertical":
        transformed["y"] = y_spots - 1 - transformed["y"]
    elif orientation == "rotate":
        transformed["x"] = x_spots - 1 - transformed["x"]
        transformed["y"] = y_spots - 1 - transformed["y"]
    if swap_xy:
        transformed["x"], transformed["y"] = (
            transformed["y"].copy(),
            transformed["x"].copy(),
        )
        x_spots, y_spots = y_spots, x_spots
    if rotate == 90:
        old_x = transformed["x"].copy()
        transformed["x"] = y_spots - 1 - transformed["y"]
        transformed["y"] = old_x
        x_spots, y_spots = y_spots, x_spots
    elif rotate == 180:
        transformed["x"] = x_spots - 1 - transformed["x"]
        transformed["y"] = y_spots - 1 - transformed["y"]
    elif rotate == 270:
        old_x = transformed["x"].copy()
        transformed["x"] = transformed["y"]
        transformed["y"] = x_spots - 1 - old_x
        x_spots, y_spots = y_spots, x_spots
    transformed["plot_x"] = transformed["x"]
    transformed["plot_y"] = transformed["y"]
    return transformed, x_spots, y_spots


def figure_layout(
    x_spots: int, y_spots: int, side_panel_width: float
) -> tuple[tuple[float, float], list[float]]:
    inches_per_spot = BASE_PLOT_SIZE / BASE_GRID_SPOTS
    plot_width = max(BASE_PLOT_SIZE, x_spots * inches_per_spot)
    plot_height = max(BASE_PLOT_SIZE, y_spots * inches_per_spot)
    return (
        (plot_width + side_panel_width, plot_height + VERTICAL_MARGIN),
        [plot_width, side_panel_width],
    )


def grid_cells(frame: pd.DataFrame) -> list[Rectangle]:
    half_side = SPOT_SIDE_LENGTH / 2
    return [
        Rectangle(
            (row.plot_x - half_side, row.plot_y - half_side),
            SPOT_SIDE_LENGTH,
            SPOT_SIDE_LENGTH,
        )
        for row in frame.itertuples(index=False)
    ]


def configure_spatial_axis(axis, x_spots: int, y_spots: int, title: str) -> None:
    axis.set_xlim(-PLOT_EDGE_PAD, x_spots - 1 + PLOT_EDGE_PAD)
    axis.set_ylim(-PLOT_EDGE_PAD, y_spots - 1 + PLOT_EDGE_PAD)
    axis.set_aspect("equal")
    axis.invert_yaxis()
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_title(title)


def save_continuous_plot(
    frame: pd.DataFrame,
    values: np.ndarray,
    output_path: Path,
    title: str,
    colorbar_label: str,
    x_spots: int,
    y_spots: int,
    cmap: str = "Reds",
) -> None:
    figure_size, _ = figure_layout(
        x_spots, y_spots, COLORBAR_GAP + COLORBAR_WIDTH
    )
    figure, (axis, gap_axis, colorbar_axis) = plt.subplots(
        ncols=3,
        figsize=figure_size,
        gridspec_kw={
            "width_ratios": [
                figure_size[0] - COLORBAR_GAP - COLORBAR_WIDTH,
                COLORBAR_GAP,
                COLORBAR_WIDTH,
            ]
        },
    )
    gap_axis.axis("off")
    normalization = mcolors.Normalize(vmin=0.0, vmax=1.0)
    collection = PatchCollection(
        grid_cells(frame),
        cmap=cmap,
        norm=normalization,
        edgecolors="none",
    )
    collection.set_array(np.asarray(values, dtype=float))
    axis.add_collection(collection)
    configure_spatial_axis(axis, x_spots, y_spots, title)
    colorbar = figure.colorbar(
        collection, cax=colorbar_axis, ticks=np.linspace(0.0, 1.0, 5)
    )
    colorbar.set_label(colorbar_label)
    colorbar.ax.yaxis.set_label_position("left")
    colorbar.ax.tick_params(labelsize=8)
    figure.subplots_adjust(
        left=0.028, right=0.90, top=0.92, bottom=0.04, wspace=0.0
    )
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def categorical_colors(cell_types: list[str]) -> dict[str, str]:
    if len(cell_types) <= TOP_DOMINANT_TYPES:
        palette = [
            color
            for index, color in enumerate(plt.get_cmap("tab20").colors)
            if index not in (14, 15)
        ]
        colors = palette[: len(cell_types)]
    else:
        palette = plt.get_cmap("gist_ncar", len(cell_types))
        colors = [palette(index) for index in range(len(cell_types))]
    return {
        cell_type: mcolors.to_hex(colors[index])
        for index, cell_type in enumerate(cell_types)
    }


def wrap_legend_label(label: str) -> str:
    return textwrap.fill(
        label,
        width=LEGEND_LABEL_WRAP_WIDTH,
        break_long_words=False,
        break_on_hyphens=False,
    )


def legend_column_width(labels: list[str]) -> float:
    if not labels:
        return LEGEND_COLUMN_WIDTH
    wrapped_line_lengths = [
        len(line)
        for label in labels
        for line in wrap_legend_label(label).splitlines()
    ]
    longest_line = max(wrapped_line_lengths, default=0)
    return max(
        LEGEND_COLUMN_WIDTH,
        longest_line * LEGEND_LABEL_INCHES_PER_CHAR + LEGEND_LABEL_EXTRA_WIDTH,
    )


def save_dominant_plot(
    frame: pd.DataFrame,
    cell_types: list[str],
    output_path: Path,
    x_spots: int,
    y_spots: int,
) -> None:
    dominant_set = set(frame["dominant_cell_type"])
    dominant_counts = frame["dominant_cell_type"].value_counts()
    cell_type_order = {cell_type: index for index, cell_type in enumerate(cell_types)}
    ranked_types = sorted(
        dominant_set,
        key=lambda cell_type: (
            -int(dominant_counts[cell_type]),
            cell_type_order[cell_type],
        ),
    )
    dominant_types = ranked_types[:TOP_DOMINANT_TYPES]
    has_other = len(ranked_types) > TOP_DOMINANT_TYPES
    displayed_types = [
        *dominant_types,
        *([OTHER_CELL_TYPE] if has_other else []),
    ]
    legend_columns = max(1, math.ceil(len(displayed_types) / 25))
    legend_width = legend_column_width(displayed_types) * legend_columns
    figure_size, width_ratios = figure_layout(x_spots, y_spots, legend_width)
    figure, (axis, legend_axis) = plt.subplots(
        ncols=2,
        figsize=figure_size,
        gridspec_kw={"width_ratios": width_ratios},
    )
    color_map = categorical_colors(dominant_types)
    if has_other:
        color_map[OTHER_CELL_TYPE] = "#bdbdbd"
    display_values = frame["dominant_cell_type"].where(
        frame["dominant_cell_type"].isin(dominant_types),
        OTHER_CELL_TYPE,
    )
    cell_colors = [color_map[value] for value in display_values]
    axis.add_collection(
        PatchCollection(
            grid_cells(frame),
            facecolors=cell_colors,
            edgecolors="none",
        )
    )
    configure_spatial_axis(
        axis,
        x_spots,
        y_spots,
        f"RCTD dominant cell type (top {TOP_DOMINANT_TYPES} + Other)",
    )
    legend_axis.axis("off")
    legend_axis.legend(
        handles=[
            Patch(facecolor=color_map[cell_type], label=wrap_legend_label(cell_type))
            for cell_type in displayed_types
        ],
        title="cell type",
        loc="upper left",
        frameon=False,
        fontsize=7,
        title_fontsize=9,
        borderaxespad=0,
        handletextpad=0.45,
        labelspacing=0.4,
        ncol=legend_columns,
        columnspacing=0.8,
    )
    figure.subplots_adjust(
        left=0.028, right=0.995, top=0.92, bottom=0.04, wspace=0.02
    )
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def safe_file_component(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned or "cell_type"


def main() -> None:
    args = parse_args()
    if args.x_spots_number < 1 or args.y_spots_number < 1:
        raise SystemExit("x/y spots number must be positive integers")

    frame, cell_types = load_weights(args.weights_file)
    coordinates = frame[["x", "y"]].to_numpy(dtype=float)
    if not np.allclose(coordinates, np.rint(coordinates)):
        raise SystemExit("Spatial plots require integer x/y coordinates")
    frame[["x", "y"]] = np.rint(coordinates).astype(int)
    outside = (
        (frame["x"] < 0)
        | (frame["x"] >= args.x_spots_number)
        | (frame["y"] < 0)
        | (frame["y"] >= args.y_spots_number)
    )
    if outside.any():
        invalid = frame.loc[outside, ["x", "y"]].iloc[0]
        raise SystemExit(
            f"Coordinate ({invalid['x']}, {invalid['y']}) is outside the configured grid"
        )

    weights = frame[cell_types].to_numpy(dtype=float)
    dominant_indices = np.argmax(weights, axis=1)
    frame["dominant_cell_type"] = [cell_types[index] for index in dominant_indices]
    frame["max_weight"] = weights[np.arange(len(frame)), dominant_indices]
    if len(cell_types) == 1:
        frame["entropy"] = 0.0
    else:
        entropy_terms = np.zeros_like(weights)
        positive = weights > 0
        entropy_terms[positive] = weights[positive] * np.log(weights[positive])
        frame["entropy"] = -entropy_terms.sum(axis=1) / np.log(len(cell_types))

    frame, plot_x_spots, plot_y_spots = transform_coordinates(
        frame,
        x_spots=args.x_spots_number,
        y_spots=args.y_spots_number,
        orientation=args.orientation,
        swap_xy=args.swap_xy,
        rotate=args.rotate,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell_type_dir = args.output_dir / "cell_types"
    cell_type_dir.mkdir(parents=True, exist_ok=True)
    for old_plot in cell_type_dir.glob("*.png"):
        old_plot.unlink()

    save_dominant_plot(
        frame,
        cell_types,
        args.output_dir / "dominant_cell_type.png",
        plot_x_spots,
        plot_y_spots,
    )
    save_continuous_plot(
        frame,
        frame["max_weight"].to_numpy(),
        args.output_dir / "max_weight.png",
        "RCTD maximum cell-type weight",
        "maximum weight",
        plot_x_spots,
        plot_y_spots,
    )
    save_continuous_plot(
        frame,
        frame["entropy"].to_numpy(),
        args.output_dir / "entropy.png",
        "RCTD normalized cell-type entropy",
        "normalized entropy",
        plot_x_spots,
        plot_y_spots,
    )
    for index, cell_type in enumerate(cell_types, start=1):
        filename = f"{index:03d}_{safe_file_component(cell_type)}.png"
        save_continuous_plot(
            frame,
            frame[cell_type].to_numpy(),
            cell_type_dir / filename,
            cell_type,
            "normalized weight",
            plot_x_spots,
            plot_y_spots,
        )

    print(
        f"Saved RCTD visualizations for {len(cell_types)} cell types to "
        f"{args.output_dir}"
    )


if __name__ == "__main__":
    main()
