#!/usr/bin/env python3
"""Cluster RCTD weights, optionally with expression, into spatial domains."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
from banksy.cluster_methods import run_Leiden_partition
from banksy.embed_banksy import generate_banksy_matrix
from banksy.initialize_banksy import initialize_banksy
from banksy_utils.umap_pca import pca_umap
import matplotlib
import scanpy as sc

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse


BASE_GRID_SPOTS = 50
BASE_PLOT_SIZE = 4.37
LEGEND_WIDTH = 0.63
VERTICAL_MARGIN = 0.73
SPOT_SIDE_LENGTH = 0.82
PLOT_EDGE_PAD = 0.9
DOWNSAMPLE_FACTOR = 10
MIN_OUTPUT_DIMENSION = 1500
CELL_TYPE_DIFFERENCE_THRESHOLD = 0.05
PEARSON_RESIDUAL_THETA = 100


def spatial_figure_layout(
    x_spots: int, y_spots: int, cluster_count: int
) -> tuple[tuple[float, float], list[float], dict]:
    """Scale the spatial panel and make enough room for large legends."""
    inches_per_spot = BASE_PLOT_SIZE / BASE_GRID_SPOTS
    plot_width = max(BASE_PLOT_SIZE, x_spots * inches_per_spot)
    plot_height = max(BASE_PLOT_SIZE, y_spots * inches_per_spot)
    figure_height = plot_height + VERTICAL_MARGIN
    legend_font_size = 8 if cluster_count <= 30 else 7 if cluster_count <= 60 else 6
    row_height = 0.22 if legend_font_size >= 8 else 0.19 if legend_font_size == 7 else 0.16
    max_legend_rows = max(8, int((figure_height - 0.45) / row_height))
    legend_columns = max(1, (cluster_count + max_legend_rows - 1) // max_legend_rows)
    legend_width = max(LEGEND_WIDTH, 0.64 * legend_columns + 0.1)
    return (
        (plot_width + legend_width, figure_height),
        [plot_width, legend_width],
        {
            "ncols": legend_columns,
            "fontsize": legend_font_size,
            "title_fontsize": legend_font_size + 1,
        },
    )


def draw_spot_grid(axis, frame: pd.DataFrame, colors) -> None:
    """Draw equal-sized grid cells with a uniform gap in data coordinates."""
    half_side = SPOT_SIDE_LENGTH / 2
    cells = [
        Rectangle(
            (row.plot_x - half_side, row.plot_y - half_side),
            SPOT_SIDE_LENGTH,
            SPOT_SIDE_LENGTH,
        )
        for row in frame.itertuples(index=False)
    ]
    collection = PatchCollection(
        cells,
        facecolors=list(colors),
        edgecolors="none",
    )
    axis.add_collection(collection)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster the cell_type_weights.csv output from spacexr.R "
            "into spatial domains with BANKSY and Leiden."
        )
    )
    parser.add_argument(
        "-i",
        "--weights-file",
        required=True,
        type=Path,
        help="RCTD cell_type_weights.csv file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "New or empty output directory (default: banksy_output next to "
            "the input file)."
        ),
    )
    parser.add_argument(
        "--expression-h5ad",
        type=Path,
        help=(
            "Optional raw-count expression H5AD. Existing metadata and embeddings "
            "are discarded before expression PCA is recomputed."
        ),
    )
    parser.add_argument(
        "--expression-hvg-count",
        type=int,
        default=3000,
        help="Number of expression highly variable genes (default: 3000).",
    )
    parser.add_argument(
        "--expression-components",
        type=int,
        default=20,
        help="Number of expression embedding components to use (default: 20).",
    )
    parser.add_argument(
        "--expression-weight",
        type=float,
        default=0.3,
        help=(
            "Target expression contribution in the joint feature variance; "
            "RCTD receives 1 minus this value (default: 0.3)."
        ),
    )
    parser.add_argument(
        "--lambda-param",
        type=float,
        default=0.8,
        help="Contribution of spatial-neighbour features (default: 0.8).",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Leiden clustering resolution (default: 1.0).",
    )
    parser.add_argument(
        "--spatial-neighbors",
        type=int,
        default=30,
        help="Number of spatial neighbours used by BANKSY (default: 30).",
    )
    parser.add_argument(
        "--cluster-neighbors",
        type=int,
        default=50,
        help="Number of neighbours used for Leiden clustering (default: 50).",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        default=20,
        help="Number of BANKSY principal components (default: 20).",
    )
    parser.add_argument(
        "--max-m",
        type=int,
        choices=(0, 1),
        default=1,
        help="Maximum azimuthal transform order (default: 1).",
    )
    parser.add_argument(
        "--neighbor-decay",
        choices=("scaled_gaussian", "reciprocal"),
        default="scaled_gaussian",
        help="Spatial-neighbour weight decay (default: scaled_gaussian).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Leiden partition seed (default: 42).",
    )
    parser.add_argument(
        "--x-spots-number",
        type=int,
        required=True,
        help="Number of spots along the x axis in the complete DBiT grid.",
    )
    parser.add_argument(
        "--y-spots-number",
        type=int,
        required=True,
        help="Number of spots along the y axis in the complete DBiT grid.",
    )
    parser.add_argument(
        "--length-spot",
        type=int,
        required=True,
        help="Spot square size in pixels before resizing (default: 20).",
    )
    parser.add_argument(
        "--interval",
        type=int,
        required=True,
        help="Gap between neighbouring spots in pixels before resizing (default: 20).",
    )
    parser.add_argument(
        "--pixel-length",
        type=float,
        required=True,
        help="Pixel scaling factor used by frame_filtered.py (default: 1.0).",
    )
    parser.add_argument(
        "--orientation",
        choices=("normal", "horizontal", "vertical", "rotate"),
        help="Grid origin orientation: normal, horizontal, vertical, or rotate.",
    )
    parser.add_argument(
        "--swap_xy",
        action="store_true",
        help="Swap x and y axes after applying orientation.",
    )
    parser.add_argument(
        "--rotate",
        type=int,
        choices=(0, 90, 180, 270),
        default=0,
        help="Rotate the plotted spatial grid clockwise (default: 0).",
    )
    parser.add_argument(
        "--no-normalize-weights",
        action="store_true",
        help="Do not normalize each spot's RCTD weights to sum to one.",
    )
    parser.add_argument(
        "--subcluster",
        action="store_true",
        help="Recluster large BANKSY domains once to produce sub-domain labels.",
    )
    parser.add_argument(
        "--subcluster-min-parent-spots",
        type=int,
        default=40,
        help="Minimum parent-domain spot count required for subclustering (default: 40).",
    )
    parser.add_argument(
        "--subcluster-min-spots",
        type=int,
        default=10,
        help="Minimum accepted subcluster spot count (default: 10).",
    )
    parser.add_argument(
        "--subcluster-max-depth",
        type=int,
        choices=(1,),
        default=1,
        help="Maximum subclustering depth. This first implementation supports 1.",
    )
    parser.add_argument(
        "--subcluster-resolution",
        type=float,
        default=0.5,
        help="Leiden resolution for within-domain subclustering (default: 0.5).",
    )
    parser.add_argument(
        "--subcluster-spatial-neighbors",
        type=int,
        default=8,
        help="BANKSY spatial neighbours for within-domain subclustering (default: 8).",
    )
    parser.add_argument(
        "--subcluster-cluster-neighbors",
        type=int,
        default=12,
        help="Leiden neighbours for within-domain subclustering (default: 12).",
    )
    parser.add_argument(
        "--subcluster-lambda-param",
        type=float,
        help="Spatial-neighbour contribution for subclustering (default: --lambda-param).",
    )
    parser.add_argument(
        "--subcluster-pca-components",
        type=int,
        help="PCA components for subclustering (default: --pca-components).",
    )
    parser.add_argument(
        "--subcluster-max-dominant-fraction",
        type=float,
        default=0.9,
        help=(
            "Reject subclustering if the largest subcluster keeps more than this "
            "fraction of the parent domain (default: 0.9)."
        ),
    )
    parser.add_argument(
        "--subcluster-min-differential-cell-types",
        type=int,
        default=2,
        help=(
            "Merge two subclusters when fewer than this many cell types differ "
            f"by at least {CELL_TYPE_DIFFERENCE_THRESHOLD:g} in mean RCTD weight "
            "(default: 2)."
        ),
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.weights_file.is_file():
        raise SystemExit(f"Weights file not found: {args.weights_file}")
    if args.output_dir is None:
        args.output_dir = args.weights_file.parent / "banksy_output"
    if args.expression_h5ad is not None and not args.expression_h5ad.is_file():
        raise SystemExit(f"Expression H5AD not found: {args.expression_h5ad}")
    if args.expression_hvg_count < 1:
        raise SystemExit("--expression-hvg-count must be a positive integer")
    if args.expression_components < 1:
        raise SystemExit("--expression-components must be a positive integer")
    if not 0.0 < args.expression_weight < 1.0:
        raise SystemExit("--expression-weight must be between 0 and 1")
    if not 0.0 <= args.lambda_param <= 1.0:
        raise SystemExit("--lambda-param must be between 0 and 1")
    if args.resolution <= 0:
        raise SystemExit("--resolution must be greater than zero")
    if args.spatial_neighbors < 1 or args.cluster_neighbors < 1:
        raise SystemExit("Neighbour counts must be positive integers")
    if args.pca_components < 1:
        raise SystemExit("--pca-components must be a positive integer")
    if args.x_spots_number < 1 or args.y_spots_number < 1:
        raise SystemExit("x/y spots number must be positive integers")
    if args.length_spot < 1:
        raise SystemExit("--length-spot must be a positive integer")
    if args.interval < 0:
        raise SystemExit("--interval must be a non-negative integer")
    if args.pixel_length <= 0:
        raise SystemExit("--pixel-length must be greater than zero")
    if args.subcluster_min_parent_spots < 3:
        raise SystemExit("--subcluster-min-parent-spots must be at least 3")
    if args.subcluster_min_spots < 1:
        raise SystemExit("--subcluster-min-spots must be a positive integer")
    if args.subcluster_resolution <= 0:
        raise SystemExit("--subcluster-resolution must be greater than zero")
    if args.subcluster_spatial_neighbors < 1 or args.subcluster_cluster_neighbors < 1:
        raise SystemExit("Subclustering neighbour counts must be positive integers")
    if args.subcluster_lambda_param is not None and not 0.0 <= args.subcluster_lambda_param <= 1.0:
        raise SystemExit("--subcluster-lambda-param must be between 0 and 1")
    if args.subcluster_pca_components is not None and args.subcluster_pca_components < 1:
        raise SystemExit("--subcluster-pca-components must be a positive integer")
    if not 0.0 < args.subcluster_max_dominant_fraction <= 1.0:
        raise SystemExit("--subcluster-max-dominant-fraction must be between 0 and 1")
    if args.subcluster_min_differential_cell_types < 0:
        raise SystemExit("--subcluster-min-differential-cell-types must be non-negative")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def load_weights(path: Path, normalize: bool):
    table = pd.read_csv(path)
    metadata_columns = ["barcode", "x", "y"]
    missing = [column for column in metadata_columns if column not in table.columns]
    if missing:
        raise SystemExit(
            "Weights file is missing required column(s): " + ", ".join(missing)
        )

    feature_columns = [column for column in table.columns if column not in metadata_columns]
    if not feature_columns:
        raise SystemExit("Weights file does not contain any cell-type columns")
    if table["barcode"].isna().any() or table["barcode"].astype(str).duplicated().any():
        raise SystemExit("barcode values must be non-empty and unique")

    numeric_columns = ["x", "y", *feature_columns]
    try:
        numeric = table[numeric_columns].apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise SystemExit(f"Coordinates and cell-type weights must be numeric: {error}") from error
    values = numeric[feature_columns].to_numpy(dtype=np.float64)
    coordinates = numeric[["x", "y"]].to_numpy(dtype=np.float64)
    if not np.isfinite(coordinates).all():
        raise SystemExit("Spatial coordinates contain missing or non-finite values")
    if not np.isfinite(values).all() or (values < 0).any():
        raise SystemExit("Cell-type weights must be finite and non-negative")

    row_sums = values.sum(axis=1)
    if (row_sums <= 0).any():
        bad_barcodes = table.loc[row_sums <= 0, "barcode"].astype(str).tolist()
        preview = ", ".join(bad_barcodes[:5])
        raise SystemExit(f"Cell-type weights sum to zero for spot(s): {preview}")
    if normalize:
        values = values / row_sums[:, None]

    barcodes = table["barcode"].astype(str).to_numpy()
    obs = pd.DataFrame(
        {"x": coordinates[:, 0], "y": coordinates[:, 1]}, index=barcodes
    )
    obs.index.name = "barcode"
    var = pd.DataFrame(
        {"modality": "rctd", "feature_weight": 1.0},
        index=pd.Index(feature_columns, name="feature"),
    )
    adata = ad.AnnData(X=values, obs=obs, var=var)
    adata.obsm["spatial"] = coordinates
    return adata, table, feature_columns


def add_expression_pca(
    adata,
    path: Path,
    requested_hvgs: int,
    requested_components: int,
    seed: int,
):
    source = ad.read_h5ad(path, backed="r")
    try:
        if not source.obs_names.is_unique:
            raise SystemExit("Expression H5AD obs names must be unique")
        indexer = source.obs_names.get_indexer(adata.obs_names)
        if (indexer < 0).any():
            missing = adata.obs_names[indexer < 0].tolist()
            preview = ", ".join(map(str, missing[:5]))
            raise SystemExit(
                f"Expression H5AD is missing RCTD barcode(s): {preview}"
            )

        if "spatial" in source.obsm:
            expression_coordinates = np.asarray(
                source.obsm["spatial"], dtype=np.float64
            )[indexer]
        elif {"x", "y"}.issubset(source.obs.columns):
            expression_coordinates = source.obs.iloc[indexer][
                ["x", "y"]
            ].to_numpy(dtype=np.float64)
        else:
            raise SystemExit(
                "Expression H5AD must contain obsm['spatial'] or obs x/y coordinates"
            )
        if expression_coordinates.shape != adata.obsm["spatial"].shape or not np.allclose(
            expression_coordinates,
            adata.obsm["spatial"],
            rtol=0.0,
            atol=1e-8,
        ):
            raise SystemExit(
                "Expression and RCTD spatial coordinates differ after barcode alignment"
            )

        aligned = source[indexer, :].to_memory()
    finally:
        source.file.close()

    # Keep only raw counts and identifiers. All original obs/var metadata,
    # embeddings, graphs, layers, and unstructured annotations are discarded.
    expression = ad.AnnData(
        X=aligned.X.copy(),
        obs=pd.DataFrame(
            index=pd.Index(adata.obs_names.astype(str), name="barcode")
        ),
        var=pd.DataFrame(
            index=pd.Index(aligned.var_names.astype(str), name="gene")
        ),
    )
    del aligned
    expression.var_names_make_unique()
    count_values = (
        expression.X.data
        if sparse.issparse(expression.X)
        else np.asarray(expression.X).ravel()
    )
    if not np.isfinite(count_values).all() or (count_values < 0).any():
        raise SystemExit("Expression X must contain finite, non-negative raw counts")
    if not np.allclose(count_values, np.rint(count_values)):
        raise SystemExit("Expression X must contain raw integer counts for PCA recomputation")

    hvg_count = min(requested_hvgs, expression.n_vars)
    sc.experimental.pp.highly_variable_genes(
        expression,
        flavor="pearson_residuals",
        n_top_genes=hvg_count,
        theta=PEARSON_RESIDUAL_THETA,
        inplace=True,
    )
    hvg_mask = expression.var["highly_variable"].to_numpy()
    hvg_totals = np.asarray(expression[:, hvg_mask].X.sum(axis=1)).ravel()
    if (hvg_totals <= 0).any():
        bad_barcodes = expression.obs_names[hvg_totals <= 0].tolist()
        preview = ", ".join(map(str, bad_barcodes[:5]))
        raise SystemExit(
            "Expression PCA cannot be computed for spot(s) with zero counts "
            f"across selected HVGs: {preview}"
        )

    components = min(
        requested_components,
        expression.n_obs - 1,
        hvg_count - 1,
    )
    if components < 1:
        raise SystemExit("Not enough expression spots or genes to recompute PCA")
    expression_pca = expression[:, hvg_mask].copy()
    sc.experimental.pp.normalize_pearson_residuals(
        expression_pca,
        theta=PEARSON_RESIDUAL_THETA,
        inplace=True,
    )
    sc.pp.pca(
        expression_pca,
        n_comps=components,
        random_state=seed,
    )
    aligned_embedding = np.asarray(
        expression_pca.obsm["X_pca"][:, :components], dtype=np.float64
    )
    if not np.isfinite(aligned_embedding).all():
        raise SystemExit("Recomputed expression PCA contains non-finite values")
    explained_variance_ratio = np.asarray(
        expression_pca.uns["pca"]["variance_ratio"][:components],
        dtype=np.float64,
    )
    if (
        not np.isfinite(explained_variance_ratio).all()
        or (explained_variance_ratio < 0).any()
        or explained_variance_ratio.sum() <= 0
    ):
        raise SystemExit("Expression PCA explained variance ratios are invalid")
    expression_feature_weights = (
        explained_variance_ratio / explained_variance_ratio.sum()
    )

    expression_names = [
        f"expression_pc_{component + 1:02d}" for component in range(components)
    ]
    joint_values = np.concatenate(
        [dense_values(adata.X), aligned_embedding], axis=1
    )
    joint_var = pd.concat(
        [
            adata.var.copy(),
            pd.DataFrame(
                {
                    "modality": "expression",
                    "feature_weight": expression_feature_weights,
                },
                index=pd.Index(expression_names, name=adata.var_names.name),
            ),
        ]
    )
    joint = ad.AnnData(X=joint_values, obs=adata.obs.copy(), var=joint_var)
    joint.obsm["spatial"] = adata.obsm["spatial"].copy()
    metadata = {
        "h5ad": str(path.resolve()),
        "input_metadata_discarded": True,
        "preprocessing": "pearson_residuals",
        "pearson_residual_theta": PEARSON_RESIDUAL_THETA,
        "input_spots": int(expression.n_obs),
        "input_genes": int(expression.n_vars),
        "requested_hvgs": requested_hvgs,
        "used_hvgs": int(expression_pca.n_vars),
        "requested_components": requested_components,
        "used_components": components,
        "expression_pca_recomputed": True,
        "pca_variance_weighted": True,
        "selected_pca_explained_variance_ratio": explained_variance_ratio.tolist(),
        "normalized_pca_feature_weights": expression_feature_weights.tolist(),
    }
    return joint, metadata


def transform_coordinates(
    frame,
    x_spots: int,
    y_spots: int,
    orientation: str,
    swap_xy: bool,
    rotate: int,
):
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


def spot_frame_size(
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
) -> tuple[int, int]:
    width = int(x_spots * length_spot + (x_spots - 1) * interval)
    height = int(y_spots * length_spot + (y_spots - 1) * interval)
    return width, height


def render_spot_frame(
    frame: pd.DataFrame,
    color_map: dict,
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
) -> np.ndarray:
    width, height = spot_frame_size(x_spots, y_spots, length_spot, interval)
    canvas = np.zeros((height, width, 4), dtype=np.uint8)
    for row in frame.itertuples(index=False):
        x_idx = int(row.x)
        y_idx = int(row.y)
        x_start = x_idx * (length_spot + interval)
        y_start = y_idx * (length_spot + interval)
        x_end = x_start + length_spot
        y_end = y_start + length_spot
        canvas[y_start:y_end, x_start:x_end, :] = color_map[row.cluster]
    return canvas


def resize_spot_image(
    image: Image.Image,
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
    pixel_length: float,
) -> Image.Image:
    width, height = spot_frame_size(x_spots, y_spots, length_spot, interval)
    resized_width = max(
        MIN_OUTPUT_DIMENSION,
        int(width / pixel_length) // DOWNSAMPLE_FACTOR,
    )
    resized_height = max(
        MIN_OUTPUT_DIMENSION,
        int(height / pixel_length) // DOWNSAMPLE_FACTOR,
    )
    return image.resize((resized_width, resized_height), resample=Image.NEAREST)


def run_banksy_partition(
    adata,
    lambda_param: float,
    resolution: float,
    spatial_neighbors_arg: int,
    cluster_neighbors_arg: int,
    pca_components_arg: int,
    max_m: int,
    neighbor_decay: str,
    seed: int,
    expression_weight: float,
):
    if adata.n_obs < 3:
        raise SystemExit("BANKSY clustering requires at least three spots")

    spatial_neighbors = min(spatial_neighbors_arg, adata.n_obs - 1)
    banksy_dict = initialize_banksy(
        adata,
        coord_keys=("x", "y", "spatial"),
        num_neighbours=spatial_neighbors,
        nbr_weight_decay=neighbor_decay,
        max_m=max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
    )
    banksy_dict, banksy_adata = generate_banksy_matrix(
        adata,
        banksy_dict,
        lambda_list=[float(lambda_param)],
        max_m=max_m,
        plot_std=False,
        save_matrix=False,
        verbose=False,
    )
    modalities = banksy_adata.var["modality"].astype(str)
    if "expression" in set(modalities):
        original_features = ~banksy_adata.var["is_nbr"].astype(bool)
        feature_weights = banksy_adata.var["feature_weight"].to_numpy(
            dtype=np.float64
        )
        column_scale_squared = np.zeros(banksy_adata.n_vars, dtype=np.float64)
        modality_targets = {
            "rctd": 1.0 - expression_weight,
            "expression": expression_weight,
        }
        for modality, target_weight in modality_targets.items():
            modality_mask = modalities.eq(modality).to_numpy()
            original_mask = original_features.to_numpy() & modality_mask
            original_weight_sum = feature_weights[original_mask].sum()
            if original_mask.sum() < 1 or original_weight_sum <= 0:
                raise RuntimeError(
                    f"Joint BANKSY input has invalid {modality} feature weights"
                )
            column_scale_squared[modality_mask] = (
                target_weight
                * feature_weights[modality_mask]
                / original_weight_sum
            )
        column_scales = np.sqrt(column_scale_squared)
        banksy_adata.X = np.asarray(banksy_adata.X) * column_scales[None, :]

    max_components = min(banksy_adata.n_obs - 1, banksy_adata.n_vars)
    pca_components = min(pca_components_arg, max_components)
    pca_umap(
        banksy_dict,
        pca_dims=[pca_components],
        plt_remaining_var=False,
        add_umap=False,
    )
    cluster_neighbors = min(cluster_neighbors_arg, adata.n_obs - 1)
    results, _ = run_Leiden_partition(
        banksy_dict=banksy_dict,
        resolutions=[float(resolution)],
        num_nn=cluster_neighbors,
        num_iterations=-1,
        partition_seed=seed,
        match_labels=False,
        annotations=None,
    )
    if len(results) != 1:
        raise RuntimeError(f"Expected one BANKSY result, received {len(results)}")

    result = results.iloc[0]
    labels = result["labels"]
    labels = labels.dense if hasattr(labels, "dense") else labels
    result_adata = result["adata"]
    result_adata.obsm["spatial"] = result_adata.obs[["x", "y"]].to_numpy()
    result_adata.obs["banksy_cluster"] = [str(label) for label in labels]
    result_adata.obs["banksy_cluster"] = result_adata.obs["banksy_cluster"].astype(
        "category"
    )
    return result_adata, labels, spatial_neighbors, cluster_neighbors, pca_components


def natural_cluster_key(label) -> tuple:
    parts = str(label).split(".")
    key = []
    for part in parts:
        key.append(int(part) if part.isdigit() else part)
    return tuple(key)


def safe_cluster_label(label) -> str:
    return str(label).replace(".", "_")


def dense_values(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def merge_similar_subclusters(sub_adata, sub_labels, min_differential_cell_types: int):
    if min_differential_cell_types <= 0:
        return sub_labels, []

    labels = np.asarray([str(label) for label in sub_labels], dtype=object)
    unique_labels = sorted(pd.unique(labels), key=natural_cluster_key)
    if len(unique_labels) < 2:
        return labels, []

    rctd_mask = sub_adata.var["modality"].astype(str).eq("rctd").to_numpy()
    values = dense_values(sub_adata.X)[:, rctd_mask]
    means = {
        label: values[labels == label].mean(axis=0)
        for label in unique_labels
    }
    roots = {label: label for label in unique_labels}

    def find(label):
        while roots[label] != label:
            roots[label] = roots[roots[label]]
            label = roots[label]
        return label

    def union(left, right):
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return left_root
        new_root = min(left_root, right_root, key=natural_cluster_key)
        old_root = right_root if new_root == left_root else left_root
        roots[old_root] = new_root
        return new_root

    merged_pairs = []
    for left_index, left in enumerate(unique_labels):
        for right in unique_labels[left_index + 1:]:
            differential_count = int(
                np.count_nonzero(
                    np.abs(means[left] - means[right])
                    >= CELL_TYPE_DIFFERENCE_THRESHOLD
                )
            )
            if differential_count < min_differential_cell_types:
                union(left, right)
                merged_pairs.append(
                    {
                        "left": left,
                        "right": right,
                        "differential_cell_types": differential_count,
                    }
                )

    if not merged_pairs:
        return labels, []

    merged_labels = np.asarray([find(label) for label in labels], dtype=object)
    return merged_labels, merged_pairs


def subcluster_domains(adata, parent_labels, args: argparse.Namespace):
    final_labels = np.asarray([str(label) for label in parent_labels], dtype=object)
    parent_labels = np.asarray([str(label) for label in parent_labels], dtype=object)
    summary = []

    sub_lambda = (
        args.lambda_param
        if args.subcluster_lambda_param is None
        else args.subcluster_lambda_param
    )
    sub_pca_components = (
        args.pca_components
        if args.subcluster_pca_components is None
        else args.subcluster_pca_components
    )

    for parent in sorted(pd.unique(parent_labels), key=natural_cluster_key):
        parent_mask = parent_labels == parent
        parent_size = int(parent_mask.sum())
        parent_summary = {
            "parent": parent,
            "parent_spots": parent_size,
            "accepted": False,
            "reason": "",
            "subcluster_counts": {},
        }
        if parent_size < args.subcluster_min_parent_spots:
            parent_summary["reason"] = "parent_below_min_spots"
            summary.append(parent_summary)
            continue

        sub_adata = adata[parent_mask].copy()
        try:
            _, sub_labels, sub_spatial_neighbors, sub_cluster_neighbors, sub_pca = (
                run_banksy_partition(
                    adata=sub_adata,
                    lambda_param=sub_lambda,
                    resolution=args.subcluster_resolution,
                    spatial_neighbors_arg=args.subcluster_spatial_neighbors,
                    cluster_neighbors_arg=args.subcluster_cluster_neighbors,
                    pca_components_arg=sub_pca_components,
                    max_m=args.max_m,
                    neighbor_decay=args.neighbor_decay,
                    seed=args.seed,
                    expression_weight=args.expression_weight,
                )
            )
        except Exception as error:
            parent_summary["reason"] = f"subclustering_failed: {error}"
            summary.append(parent_summary)
            continue

        sub_labels = np.asarray([str(label) for label in sub_labels], dtype=object)
        raw_counts = pd.Series(sub_labels).value_counts().sort_index()
        parent_summary["raw_subcluster_counts"] = {
            str(label): int(count) for label, count in raw_counts.items()
        }
        sub_labels, merged_pairs = merge_similar_subclusters(
            sub_adata,
            sub_labels,
            min_differential_cell_types=args.subcluster_min_differential_cell_types,
        )
        counts = pd.Series(sub_labels).value_counts().sort_index()
        parent_summary["subcluster_counts"] = {
            str(label): int(count) for label, count in counts.items()
        }
        parent_summary["merged_subcluster_pairs"] = merged_pairs
        parent_summary["min_differential_cell_types"] = (
            args.subcluster_min_differential_cell_types
        )
        parent_summary["cell_type_difference_threshold"] = (
            CELL_TYPE_DIFFERENCE_THRESHOLD
        )
        parent_summary["spatial_neighbors"] = int(sub_spatial_neighbors)
        parent_summary["cluster_neighbors"] = int(sub_cluster_neighbors)
        parent_summary["pca_components"] = int(sub_pca)

        if len(counts) < 2:
            parent_summary["reason"] = "fewer_than_two_subclusters"
            summary.append(parent_summary)
            continue
        if int(counts.min()) < args.subcluster_min_spots:
            parent_summary["reason"] = "subcluster_below_min_spots"
            summary.append(parent_summary)
            continue
        if float(counts.max() / parent_size) > args.subcluster_max_dominant_fraction:
            parent_summary["reason"] = "dominant_subcluster_too_large"
            summary.append(parent_summary)
            continue

        sub_label_map = {
            label: index
            for index, label in enumerate(sorted(counts.index, key=natural_cluster_key))
        }
        final_labels[parent_mask] = [
            f"{parent}.{sub_label_map[label]}" for label in sub_labels
        ]
        parent_summary["accepted"] = True
        parent_summary["reason"] = "accepted"
        summary.append(parent_summary)

    return final_labels, summary


def save_outputs(
    adata,
    source_table,
    feature_columns: list[str],
    parent_labels,
    final_labels,
    args: argparse.Namespace,
    spatial_neighbors: int,
    cluster_neighbors: int,
    pca_components: int,
    subcluster_summary: list[dict] | None = None,
) -> None:
    parent_labels = np.asarray([str(label) for label in parent_labels], dtype=object)
    final_labels = np.asarray([str(label) for label in final_labels], dtype=object)
    clusters = source_table.copy()
    clusters["banksy_cluster"] = parent_labels
    if args.subcluster:
        clusters["banksy_subcluster"] = final_labels
    adata.obs["banksy_cluster"] = parent_labels
    adata.obs["banksy_cluster"] = adata.obs["banksy_cluster"].astype("category")
    if args.subcluster:
        adata.obs["banksy_subcluster"] = final_labels
        adata.obs["banksy_subcluster"] = adata.obs["banksy_subcluster"].astype(
            "category"
        )
    clusters.to_csv(args.output_dir / "banksy_clusters.csv", index=False)
    adata.write_h5ad(args.output_dir / "banksy_result.h5ad")

    plot_table = pd.DataFrame(
        {
            "x": adata.obs["x"].to_numpy(),
            "y": adata.obs["y"].to_numpy(),
            "cluster": final_labels,
        }
    )
    coordinates = plot_table[["x", "y"]].to_numpy(dtype=float)
    if not np.allclose(coordinates, np.rint(coordinates)):
        raise SystemExit("Square grid output requires integer x/y coordinates")
    plot_table[["x", "y"]] = np.rint(coordinates).astype(int)
    outside_grid = (
        (plot_table["x"] < 0)
        | (plot_table["x"] >= args.x_spots_number)
        | (plot_table["y"] < 0)
        | (plot_table["y"] >= args.y_spots_number)
    )
    if outside_grid.any():
        invalid = plot_table.loc[outside_grid, ["x", "y"]].iloc[0]
        raise SystemExit(
            f"Coordinate ({invalid['x']}, {invalid['y']}) is outside the configured "
            f"{args.x_spots_number}x{args.y_spots_number} grid"
        )
    grid_table = plot_table.copy()
    plot_table, plot_x_spots, plot_y_spots = transform_coordinates(
        plot_table,
        x_spots=args.x_spots_number,
        y_spots=args.y_spots_number,
        orientation=args.orientation,
        swap_xy=args.swap_xy,
        rotate=args.rotate,
    )

    cluster_ids = sorted(pd.unique(plot_table["cluster"]), key=natural_cluster_key)
    cmap = plt.get_cmap("tab20", len(cluster_ids))
    cluster_colors = {
        cluster_id: np.rint(np.asarray(cmap(index)) * 255).astype(np.uint8)
        for index, cluster_id in enumerate(cluster_ids)
    }
    legend_colors = {
        cluster_id: mcolors.to_hex(cluster_colors[cluster_id] / 255.0)
        for cluster_id in cluster_ids
    }
    point_colors = plot_table["cluster"].map(legend_colors)

    figure_size, width_ratios, legend_options = spatial_figure_layout(
        plot_x_spots,
        plot_y_spots,
        cluster_count=len(cluster_ids),
    )
    figure, (axis, legend_axis) = plt.subplots(
        ncols=2,
        figsize=figure_size,
        gridspec_kw={"width_ratios": width_ratios},
    )
    draw_spot_grid(axis, plot_table, point_colors)
    axis.set_xlim(-PLOT_EDGE_PAD, plot_x_spots - 1 + PLOT_EDGE_PAD)
    axis.set_ylim(-PLOT_EDGE_PAD, plot_y_spots - 1 + PLOT_EDGE_PAD)
    axis.set_aspect("equal")
    axis.invert_yaxis()
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_title("BANKSY spatial domains")
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="None",
            markerfacecolor=legend_colors[cluster_id],
            markeredgewidth=0,
            label=str(cluster_id),
        )
        for cluster_id in cluster_ids
    ]
    legend_axis.axis("off")
    legend_axis.legend(
        handles=legend_handles,
        title="cluster",
        loc="upper left",
        frameon=False,
        fontsize=legend_options["fontsize"],
        title_fontsize=legend_options["title_fontsize"],
        ncols=legend_options["ncols"],
        borderaxespad=0,
        handletextpad=0.45,
        columnspacing=0.8,
        labelspacing=0.45,
        handlelength=1.15,
        handleheight=1.15,
        bbox_to_anchor=(0.0, 1.01),
    )
    figure.subplots_adjust(left=0.028, right=0.998, top=0.92, bottom=0.04, wspace=0.0)
    figure.savefig(args.output_dir / "banksy_clusters.png", dpi=300, bbox_inches="tight")
    plt.close(figure)

    grid_dir = args.output_dir / "cluster_grids"
    grid_dir.mkdir(parents=True, exist_ok=True)
    for cluster_id in cluster_ids:
        members = grid_table[grid_table["cluster"] == cluster_id]
        grid = render_spot_frame(
            members[["x", "y", "cluster"]],
            color_map=cluster_colors,
            x_spots=args.x_spots_number,
            y_spots=args.y_spots_number,
            length_spot=args.length_spot,
            interval=args.interval,
        )
        resize_spot_image(
            Image.fromarray(grid, mode="RGBA"),
            x_spots=args.x_spots_number,
            y_spots=args.y_spots_number,
            length_spot=args.length_spot,
            interval=args.interval,
            pixel_length=args.pixel_length,
        ).save(grid_dir / f"cluster_{safe_cluster_label(cluster_id)}.png")

    parameters = {
        "weights_file": str(args.weights_file.resolve()),
        "normalized_weights": not args.no_normalize_weights,
        "cell_type_columns": feature_columns,
        "expression": getattr(args, "expression_metadata", None),
        "expression_weight": (
            args.expression_weight if args.expression_h5ad is not None else None
        ),
        "rctd_weight": (
            1.0 - args.expression_weight
            if args.expression_h5ad is not None
            else None
        ),
        "joint_pca_recomputed": args.expression_h5ad is not None,
        "modality_balancing": (
            "scale each BANKSY feature by sqrt(modality_weight * "
            "normalized_original_feature_weight); RCTD features are equal and "
            "expression PCs follow explained variance"
            if args.expression_h5ad is not None
            else None
        ),
        "lambda_param": args.lambda_param,
        "resolution": args.resolution,
        "spatial_neighbors": spatial_neighbors,
        "cluster_neighbors": cluster_neighbors,
        "pca_components": pca_components,
        "max_m": args.max_m,
        "neighbor_decay": args.neighbor_decay,
        "seed": args.seed,
        "x_spots_number": args.x_spots_number,
        "y_spots_number": args.y_spots_number,
        "length_spot": args.length_spot,
        "interval": args.interval,
        "pixel_length": args.pixel_length,
        "orientation": args.orientation,
        "swap_xy": args.swap_xy,
        "rotate": args.rotate,
        "subcluster": args.subcluster,
    }
    if args.subcluster:
        parameters["subcluster_parameters"] = {
            "min_parent_spots": args.subcluster_min_parent_spots,
            "min_spots": args.subcluster_min_spots,
            "max_depth": args.subcluster_max_depth,
            "lambda_param": (
                args.lambda_param
                if args.subcluster_lambda_param is None
                else args.subcluster_lambda_param
            ),
            "resolution": args.subcluster_resolution,
            "spatial_neighbors": args.subcluster_spatial_neighbors,
            "cluster_neighbors": args.subcluster_cluster_neighbors,
            "pca_components": (
                args.pca_components
                if args.subcluster_pca_components is None
                else args.subcluster_pca_components
            ),
            "max_dominant_fraction": args.subcluster_max_dominant_fraction,
            "min_differential_cell_types": args.subcluster_min_differential_cell_types,
            "cell_type_difference_threshold": CELL_TYPE_DIFFERENCE_THRESHOLD,
        }
        parameters["subcluster_summary"] = subcluster_summary or []
    with (args.output_dir / "run_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(parameters, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def main() -> None:
    args = parse_args()
    validate_args(args)
    adata, source_table, feature_columns = load_weights(
        args.weights_file, normalize=not args.no_normalize_weights
    )
    args.expression_metadata = None
    if args.expression_h5ad is not None:
        adata, args.expression_metadata = add_expression_pca(
            adata,
            path=args.expression_h5ad,
            requested_hvgs=args.expression_hvg_count,
            requested_components=args.expression_components,
            seed=args.seed,
        )
    result = run_banksy_partition(
        adata=adata,
        lambda_param=args.lambda_param,
        resolution=args.resolution,
        spatial_neighbors_arg=args.spatial_neighbors,
        cluster_neighbors_arg=args.cluster_neighbors,
        pca_components_arg=args.pca_components,
        max_m=args.max_m,
        neighbor_decay=args.neighbor_decay,
        seed=args.seed,
        expression_weight=args.expression_weight,
    )
    final_labels = np.asarray([str(label) for label in result[1]], dtype=object)
    subcluster_summary = []
    if args.subcluster:
        final_labels, subcluster_summary = subcluster_domains(adata, result[1], args)
    save_outputs(
        result[0],
        source_table,
        feature_columns,
        result[1],
        final_labels,
        args,
        result[2],
        result[3],
        result[4],
        subcluster_summary,
    )
    print(f"BANKSY clustering completed: {args.output_dir}")


if __name__ == "__main__":
    main()
