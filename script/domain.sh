#!/bin/bash
set -o pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)} || exit 1
R_SCRIPT="$SCRIPT_DIR/R/spacexr.R"
PYTHON_SCRIPT="$SCRIPT_DIR/python/banksy_cluster.py"
RCTD_PLOT_SCRIPT="$SCRIPT_DIR/python/rctd_visualize.py"
MERGE_SCRIPT="$SCRIPT_DIR/python/merge_on_gray.py"
CHIP_FILE="$REPO_DIR/config/chip.sh"

if [[ ! -f "$CHIP_FILE" ]]; then
    echo "Error: chip preset file not found: $CHIP_FILE" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CHIP_FILE"

if [[ $# -ne 2 ]]; then
    echo "Error: expected DBiT config and DBiT-spatial-Brain config arguments" >&2
    exit 1
fi
config_file=$1
brain_config_file=$2
for required_config in "$config_file" "$brain_config_file"; do
    if [[ ! -f "$required_config" ]]; then
        echo "Error: config file not found: $required_config" >&2
        exit 1
    fi
done

# shellcheck source=/dev/null
source "$config_file"
# shellcheck source=/dev/null
source "$brain_config_file"

normalize_dir_path() {
    local path=$1
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path=${path%/}
    done
    printf '%s\n' "$path"
}

run_id=${SLURM_JOB_ID:-domain_$$}
scratch_root=""

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ -n ${scratch_root:-} && -d $scratch_root ]]; then
        if (( status != 0 )) && [[ -n ${scratch_deconv:-} && -d $scratch_deconv && -n ${deconv_output:-} ]]; then
            echo "Recovering domain-analysis scratch outputs after exit status $status: $deconv_output" >&2
            mkdir -p "$deconv_output" && cp -a "$scratch_deconv/." "$deconv_output/" || \
                echo "Warning: failed to recover scratch outputs: $scratch_deconv" >&2
        fi
        rm -rf -- "$scratch_root" || echo "Warning: failed to clean scratch directory: $scratch_root" >&2
    fi
    exit "$status"
}

enable_cleanup() {
    trap cleanup_scratch EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP
}

run_pixi() {
    local environment=$1
    shift
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$environment" "$@"
    )
}

rctd_complete() {
    local directory=$1
    local file
    for file in cell_type_weights.csv rctd_result.rds run_info.csv session_info.txt; do
        [[ -e "$directory/$file" ]] || return 1
    done
    return 0
}

copy_rctd_outputs() {
    local source_dir=$1
    local destination_dir=$2
    local file
    mkdir -p "$destination_dir" || return 1
    for file in cell_type_weights.csv rctd_result.rds run_info.csv session_info.txt; do
        cp -a "$source_dir/$file" "$destination_dir/$file" || return 1
    done
}

pixi_env_dir=${pixi_env_dir:-$REPO_DIR}
spatial_h5ad="${mrna_dir:-}/raw/clustered.tissuefiltered.h5ad"
deconv_output="${mrna_output_path:-}/deconv"
banksy_output="$deconv_output/banksy_output"
orientation=${orientation}
swap_xy=${swap_xy}
display_rotate=${rotate}
if [[ -z ${chip:-} ]]; then
    echo "Error: chip must be set in the config. Valid chips: $(chip_preset_names_csv)." >&2
    exit 1
fi
apply_chip_preset "$chip" || {
    echo "Error: unsupported chip '$chip'." >&2
    echo "Valid chips: $(chip_preset_names_csv)." >&2
    exit 1
}

if [[ ! -d "$pixi_env_dir" ]]; then
    echo "Error: pixi environment directory not found: $pixi_env_dir" >&2
    exit 1
fi
if [[ ! -f "$R_SCRIPT" || ! -f "$PYTHON_SCRIPT" || ! -f "$RCTD_PLOT_SCRIPT" ]]; then
    echo "Error: Domain Analysis worker script is missing." >&2
    exit 1
fi
if [[ ! -f "$MERGE_SCRIPT" ]]; then
    echo "Error: merge_on_gray.py not found: $MERGE_SCRIPT" >&2
    exit 1
fi
if [[ -z ${mrna_dir:-} || -z ${mrna_output_path:-} ]]; then
    echo "Error: mrna_dir and mrna_output_path must be set by the mRNA step." >&2
    exit 1
fi

case "$orientation" in normal|horizontal|vertical|rotate) ;; *)
    echo "Error: invalid orientation: $orientation" >&2; exit 1 ;;
esac
case "$swap_xy" in
    True|true|TRUE|1|yes|YES) swap_xy=True ;;
    False|false|FALSE|0|no|NO|"") swap_xy=False ;;
    *) echo "Error: swap_xy must be True or False." >&2; exit 1 ;;
esac
case "$banksy_subcluster" in
    True|true|TRUE|1|yes|YES) banksy_subcluster=True ;;
    False|false|FALSE|0|no|NO|"") banksy_subcluster=False ;;
    *) echo "Error: banksy_subcluster must be True or False." >&2; exit 1 ;;
esac

if [[ -z "$display_rotate" ]]; then
    echo "Error: rotate must be set in the config." >&2
    exit 1
fi
case "$display_rotate" in 0|90|180|270) ;; *)
    echo "Error: rotate must be 0, 90, 180, or 270." >&2; exit 1 ;;
esac

use_scratch=false
scratch_deconv=""

if [[ -n ${scratch:-} ]]; then
    scratch=$(normalize_dir_path "$scratch")
    scratch_root="$scratch/dbit-brain/$run_id"
    scratch_deconv="$scratch_root/deconv"
    use_scratch=true
    enable_cleanup
fi

if $use_scratch; then
    run_spatial_h5ad="$scratch_root/$(basename "$spatial_h5ad")"
else
    run_spatial_h5ad="$spatial_h5ad"
fi

if rctd_complete "$deconv_output"; then
    echo "Step1 RCTD already complete, skipping."
else
    reference_dir=${rctd_reference_dir:-}
    reference_dir=$(normalize_dir_path "$reference_dir")
    if $use_scratch; then
        # Only clean our own run_id-scoped scratch directory (safe for concurrency)
        if [[ -d "$scratch_deconv" ]]; then
            echo "Warning: removing incomplete scratch RCTD output: $scratch_deconv" >&2
            rm -rf -- "$scratch_deconv" || exit 1
        fi
        mkdir -p "$scratch_root" || exit 1
        cp "$spatial_h5ad" "$run_spatial_h5ad" || exit 1
        run_deconv_output=$scratch_deconv
    else
        # Non-scratch mode: clean shared deconv_output (not concurrency-safe; use scratch for parallel runs)
        if [[ -d "$deconv_output" ]]; then
            echo "Warning: removing incomplete RCTD output: $deconv_output" >&2
            rm -rf -- "$deconv_output" || exit 1
        fi
        mkdir -p "$(dirname "$deconv_output")" || exit 1
        run_deconv_output=$deconv_output
    fi

    rctd_args=(
        Rscript "$R_SCRIPT"
        --reference-dir "$reference_dir"
        --spatial-h5ad "$run_spatial_h5ad"
        --output-dir "$run_deconv_output"
        --reference-barcode-column "${rctd_reference_barcode_column}"
        --reference-numi-column "${rctd_reference_numi_column}"
        --cell-type-column "${rctd_cell_type_column}"
        --dataset-column "${rctd_reference_dataset_column:-Dataset}"
        --reference-gene-column "${rctd_reference_gene_column}"
        --max-cells-per-type "${rctd_max_cells_per_type}"
        --cores "${rctd_cores}"
        --doublet-mode "${rctd_mode}"
        --reference-min-umi "${rctd_ref_min_umi}"
        --spatial-min-umi "${rctd_spa_min_umi}"
        --seed "${seed}"
    )
    if [[ -n ${rctd_reference_cache:-} ]]; then
        rctd_args+=(--reference-cache "$rctd_reference_cache")
    fi

    echo "Running Step1 RCTD: $run_deconv_output"
    run_pixi rctd "${rctd_args[@]}" || exit 1
    if ! rctd_complete "$run_deconv_output"; then
        echo "Error: RCTD finished without all required outputs." >&2
        exit 1
    fi
    if $use_scratch; then
        copy_rctd_outputs "$scratch_deconv" "$deconv_output" || exit 1
    fi
fi

plot_args=(
    python -B "$RCTD_PLOT_SCRIPT"
    --weights-file "$deconv_output/cell_type_weights.csv"
    --output-dir "$deconv_output/rctd_plots"
    --x-spots-number "$x_spots_number"
    --y-spots-number "$y_spots_number"
    --orientation "$orientation"
    --rotate "$display_rotate"
)
[[ "$swap_xy" == True ]] && plot_args+=(--swap-xy)

echo "Plotting Step1 RCTD cell-type weights: $deconv_output/rctd_plots"
run_pixi default "${plot_args[@]}" || exit 1

run_banksy_output=$banksy_output

banksy_args=(
    python -B "$PYTHON_SCRIPT"
    --weights-file "$deconv_output/cell_type_weights.csv"
    --output-dir "$run_banksy_output"
    --lambda-param "${banksy_lambda}"
    --resolution "${banksy_resolution}"
    --spatial-neighbors "${banksy_spatial_neighbors}"
    --cluster-neighbors "${banksy_cluster_neighbors}"
    --pca-components "${banksy_pca_components}"
    --max-m "${banksy_max_m}"
    --neighbor-decay "${banksy_neighbor_decay}"
    --seed "${seed}"
    --x-spots-number "$x_spots_number"
    --y-spots-number "$y_spots_number"
    --length-spot "$length_spot"
    --interval "$interval"
    --pixel-length "$pixel_length"
    --orientation "$orientation"
    --rotate "$display_rotate"
    --expression-h5ad "$spatial_h5ad"
    --expression-hvg-count "${banksy_expression_hvg_count}"
    --expression-components "${banksy_expression_components}"
    --expression-weight "${banksy_expression_weight}"
)
[[ "$swap_xy" == True ]] && banksy_args+=(--swap_xy)
if [[ "$banksy_subcluster" == True ]]; then
    banksy_args+=(
        --subcluster
        --subcluster-min-parent-spots "$banksy_subcluster_min_parent_spots"
        --subcluster-min-spots "$banksy_subcluster_min_spots"
        --subcluster-max-depth "$banksy_subcluster_max_depth"
        --subcluster-resolution "$banksy_subcluster_resolution"
        --subcluster-spatial-neighbors "$banksy_subcluster_spatial_neighbors"
        --subcluster-cluster-neighbors "$banksy_subcluster_cluster_neighbors"
        --subcluster-max-dominant-fraction "$banksy_subcluster_max_dominant_fraction"
        --subcluster-min-differential-cell-types "$banksy_subcluster_min_differential_cell_types"
    )
    [[ -n "$banksy_subcluster_lambda" ]] && banksy_args+=(--subcluster-lambda-param "$banksy_subcluster_lambda")
    [[ -n "$banksy_subcluster_pca_components" ]] && banksy_args+=(--subcluster-pca-components "$banksy_subcluster_pca_components")
fi

echo "Running Step2 BANKSY: $run_banksy_output"
run_pixi default "${banksy_args[@]}" || exit 1

grid_dir="$banksy_output/cluster_grids"
if [[ -n ${gray_path:-} && -f "$gray_path" && -d "$grid_dir" ]]; then
    merge_args=(python "$MERGE_SCRIPT" --gray "$gray_path" --orientation "$orientation")
    if [[ "$swap_xy" == True ]]; then
        merge_args+=(--swap_xy)
    fi
    shopt -s nullglob
    grid_frames=("$grid_dir"/cluster_*.png)
    shopt -u nullglob
    if [[ ${#grid_frames[@]} -gt 0 ]]; then
        for frame in "${grid_frames[@]}"; do
            merge_args+=(--frame "$frame")
        done
        echo "Merging BANKSY cluster grids onto gray image: $grid_dir"
        run_pixi default "${merge_args[@]}" || exit 1
    fi
elif [[ -n ${gray_path:-} && ! -f "$gray_path" ]]; then
    echo "Warning: gray image not found, skipping cluster-grid merge: $gray_path" >&2
fi

echo "Domain Analysis completed: $deconv_output"
