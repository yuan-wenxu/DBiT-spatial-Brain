#!/bin/bash
set -o pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
START_DIR=$(pwd -P) || exit 1
PROGRAM_NAME=$(basename "$0")
DBIT_BRAIN_CONFIG_FILE=${DBIT_BRAIN_CONFIG_FILE:-$REPO_DIR/config/dbit.brain.config.sh}

show_help() {
    cat <<EOF
Usage: $PROGRAM_NAME <step> [--config <file>]

Steps:
  domain        Run RCTD deconvolution and BANKSY domain analysis

Options:
  --config FILE Configuration file (default: ./dbit.config.sh)
  -h, --help    Show this help message and exit

EOF
}

require_option_value() {
    if [[ $# -lt 2 || $2 == --* ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 1
    fi
}

if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then
    show_help
    exit 0
fi

step=$1
shift
case "$step" in
    domain) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: domain." >&2
        exit 1
        ;;
esac

config_file=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            require_option_value "$@"
            config_file=$2
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: unknown option or argument '$1'." >&2
            exit 1
            ;;
    esac
done

if [[ -z "$config_file" ]]; then
    config_file="$START_DIR/dbit.config.sh"
fi
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    echo "Pass --config <file>, or create dbit.config.sh in the current directory." >&2
    exit 1
fi
config_abs=$(realpath "$config_file") || exit 1
echo "Using config: $config_abs"

if [[ ! -f "$DBIT_BRAIN_CONFIG_FILE" ]]; then
    echo "Error: DBiT-spatial-Brain config file not found: $DBIT_BRAIN_CONFIG_FILE" >&2
    echo "Copy config/dbit.brain.config.example.sh to config/dbit.brain.config.sh and edit it." >&2
    exit 1
fi
DBIT_BRAIN_CONFIG_FILE=$(realpath "$DBIT_BRAIN_CONFIG_FILE") || exit 1

# shellcheck source=/dev/null
source "$config_abs"
# shellcheck source=/dev/null
source "$DBIT_BRAIN_CONFIG_FILE"

case "${execution_mode:-}" in
    local|hpc) ;;
    *)
        echo "Error: execution_mode must be local or hpc." >&2
        exit 1
        ;;
esac

case "$step" in
    domain)
        script="$SCRIPT_DIR/domain.sh"
        cpus=${sbatch_domain_cpus:-}
        partition=${sbatch_domain_partition:-}
        memory=${sbatch_domain_mem:-}
        walltime=${sbatch_domain_time:-}
        ;;
esac

if [[ ! -f "$script" ]]; then
    echo "Error: worker script not found: $script" >&2
    exit 1
fi

if [[ "$execution_mode" == local ]]; then
    echo "Running $step locally"
    exec bash "$script" "$config_abs" "$DBIT_BRAIN_CONFIG_FILE"
fi

if [[ ! "$cpus" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: sbatch_${step}_cpus must be a positive integer." >&2
    exit 1
fi
if [[ -z "$partition" || -z "$memory" || -z "$walltime" ]]; then
    echo "Error: sbatch_${step}_partition, sbatch_${step}_mem, and sbatch_${step}_time are required in HPC mode." >&2
    exit 1
fi
if [[ -z ${sbatch_job_name_prefix:-} || -z ${sbatch_output:-} || -z ${sbatch_error:-} ]]; then
    echo "Error: sbatch_job_name_prefix, sbatch_output, and sbatch_error are required in HPC mode." >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: sbatch executable not found." >&2
    exit 1
fi

sbatch_args=(
    -J "${sbatch_job_name_prefix}_${step}"
    -c "$cpus"
    -p "$partition"
    --mem="$memory"
    --time="$walltime"
    -o "$sbatch_output"
    -e "$sbatch_error"
)
case "${sbatch_requeue:-false}" in
    True|true|TRUE|1|yes|YES) sbatch_args+=(--requeue) ;;
    False|false|FALSE|0|no|NO|"") ;;
    *) echo "Error: sbatch_requeue must be True or False." >&2; exit 1 ;;
esac

printf -v wrapped_command 'exec bash %q %q %q' \
    "$script" "$config_abs" "$DBIT_BRAIN_CONFIG_FILE"
echo "Submitting $step"
sbatch "${sbatch_args[@]}" --wrap="$wrapped_command"
