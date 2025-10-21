#!/usr/bin/env bash
# Ultra-minimal CPLEX env setup for COBRApy/optlang (no Python API required)

set -euo pipefail

log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }
log_error() { echo "[ERROR] $*" >&2; }

# Find candidate CPLEX install roots
find_cplex_roots() {
    local candidates=()
    local patterns=(
        "$HOME/ilog/cplex*"
        "/opt/ibm/ILOG/CPLEX_Studio*"
        "/opt/cplex*"
        "/usr/local/cplex*"
        "${CONDA_PREFIX:-}/ilog/cplex*"
        "${CONDA_PREFIX:-}/../ilog/cplex*"
    )

    for pattern in "${patterns[@]}"; do
        for dir in $(eval echo "$pattern" 2>/dev/null || true); do
            if [[ -d "$dir/cplex/bin/x86-64_linux" ]]; then
                candidates+=("$dir")
            fi
        done
    done

    printf '%s\n' "${candidates[@]}" | sort -u -V -r
}

# Set environment vars based on detected install
set_cplex_env() {
    local cplex_root="$1"
    local cplex_home="$cplex_root/cplex"
    local lib_path="$cplex_home/bin/x86-64_linux"

    export CPLEX_HOME="$cplex_home"
    export LD_LIBRARY_PATH="$lib_path:${LD_LIBRARY_PATH:-}"

    log_info "✓ CPLEX_HOME set to: $CPLEX_HOME"
    log_info "✓ LD_LIBRARY_PATH updated with: $lib_path"
}

# Write activation/deactivation scripts
write_conda_hooks() {
    [[ -z "${CONDA_PREFIX:-}" ]] && return 1

    local act_dir="$CONDA_PREFIX/etc/conda/activate.d"
    local deact_dir="$CONDA_PREFIX/etc/conda/deactivate.d"

    mkdir -p "$act_dir" "$deact_dir"

    cat > "$act_dir/cplex.sh" <<EOF
#!/bin/bash
export CPLEX_HOME="$CPLEX_HOME"
export LD_LIBRARY_PATH="$CPLEX_HOME/bin/x86-64_linux:\$LD_LIBRARY_PATH"
EOF

    cat > "$deact_dir/cplex.sh" <<EOF
#!/bin/bash
unset CPLEX_HOME
# Leave LD_LIBRARY_PATH intact to avoid breaking other tools
EOF

    log_info "✓ Conda activation hook saved: $act_dir/cplex.sh"
    log_info "✓ Conda deactivation hook saved: $deact_dir/cplex.sh"
}

# Main routine
main() {
    log_info "Scanning for CPLEX installations..."

    local roots
    roots=$(find_cplex_roots)

    if [[ -z "$roots" ]]; then
        log_error "❌ No valid CPLEX installations found"
        return 1
    fi

    while IFS= read -r cplex_root; do
        log_info "Trying: $cplex_root"
        set_cplex_env "$cplex_root"
        write_conda_hooks || log_warn "Conda hooks not written (not in conda env)"
        log_info "✓ CPLEX configuration complete"
        return 0
    done <<< "$roots"

    log_error "❌ Failed to configure any CPLEX environment"
    return 1
}

main "$@"
