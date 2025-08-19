#!/usr/bin/env bash
# CPLEX Discovery and Environment Setup Script for CarveWe
# This script finds CPLEX installations and configures environment variables
# for optlang/COBRApy compatibility

set -euo pipefail

# Function to log messages with timestamps
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

# Function to check if CPLEX is already configured
is_cplex_configured() {
    if [[ -n "${CPLEX_STUDIO_DIR:-}" ]] && [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
        if python -c "import cplex" 2>/dev/null; then
            return 0
        fi
    fi
    return 1
}

# Function to find CPLEX installations
find_cplex_installations() {
    local cplex_paths=()
    
    # Common CPLEX installation directories
    local search_dirs=(
        "/opt/ibm/ILOG/CPLEX_Studio*"
        "/Applications/CPLEX_Studio*"
        "$HOME/CPLEX_Studio*"
        "$HOME/Applications/CPLEX_Studio*"
        "/usr/local/CPLEX_Studio*"
        "/opt/CPLEX_Studio*"
        "${CONDA_PREFIX:-}/cplex*"
        "${CONDA_PREFIX:-}/../cplex*"
    )
    
    # Also check environment variables that might point to CPLEX
    if [[ -n "${CPLEX_STUDIO_DIR:-}" ]] && [[ -d "$CPLEX_STUDIO_DIR" ]]; then
        cplex_paths+=("$CPLEX_STUDIO_DIR")
    fi
    
    if [[ -n "${ILOG_CPLEX_PATH:-}" ]] && [[ -d "$ILOG_CPLEX_PATH" ]]; then
        cplex_paths+=("$ILOG_CPLEX_PATH")
    fi
    
    # Search for CPLEX installations
    for pattern in "${search_dirs[@]}"; do
        # Use shell expansion to find matching directories
        for dir in $pattern; do
            if [[ -d "$dir" ]]; then
                # Verify this looks like a real CPLEX installation
                if [[ -d "$dir/cplex" ]] && [[ -d "$dir/cplex/bin" ]]; then
                    cplex_paths+=("$dir")
                fi
            fi
        done 2>/dev/null || true
    done
    
    # Remove duplicates and sort by version (newest first)
    printf '%s\n' "${cplex_paths[@]}" | sort -u -V -r
}

# Function to validate CPLEX installation
validate_cplex_installation() {
    local cplex_dir="$1"
    
    # Check for required directories and files
    if [[ ! -d "$cplex_dir/cplex" ]]; then
        return 1
    fi
    
    if [[ ! -d "$cplex_dir/cplex/bin" ]]; then
        return 1
    fi
    
    # Look for Python API
    local python_api_found=false
    local python_dirs=(
        "$cplex_dir/cplex/python"
        "$cplex_dir/python"
    )
    
    for python_dir in "${python_dirs[@]}"; do
        if [[ -d "$python_dir" ]]; then
            # Check for setup.py or cplex module
            if find "$python_dir" -name "setup.py" -o -name "cplex" -type d | grep -q .; then
                python_api_found=true
                break
            fi
        fi
    done
    
    if [[ "$python_api_found" != true ]]; then
        log_warn "CPLEX installation at $cplex_dir appears to lack Python API"
        return 1
    fi
    
    return 0
}

# Function to get Python version info
get_python_version() {
    python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown"
}

# Function to find the best CPLEX Python API path
find_cplex_python_path() {
    local cplex_dir="$1"
    local python_version="$2"
    
    # Common Python API locations
    local python_paths=(
        "$cplex_dir/cplex/python/${python_version}/x86-64_linux"
        "$cplex_dir/cplex/python/${python_version}/x86-64_osx"
        "$cplex_dir/cplex/python/${python_version}/win64"
        "$cplex_dir/cplex/python/${python_version}"
        "$cplex_dir/python/${python_version}/x86-64_linux"
        "$cplex_dir/python/${python_version}/x86-64_osx"
        "$cplex_dir/python/${python_version}/win64"
        "$cplex_dir/python/${python_version}"
        "$cplex_dir/cplex/python"
        "$cplex_dir/python"
    )
    
    for path in "${python_paths[@]}"; do
        if [[ -d "$path" ]] && [[ -f "$path/setup.py" || -d "$path/cplex" ]]; then
            echo "$path"
            return 0
        fi
    done
    
    return 1
}

# Function to find CPLEX library paths
find_cplex_lib_paths() {
    local cplex_dir="$1"
    local lib_paths=()
    
    # Common library locations
    local lib_dirs=(
        "$cplex_dir/cplex/bin/x86-64_linux"
        "$cplex_dir/cplex/bin/x86-64_osx" 
        "$cplex_dir/cplex/bin/win64"
        "$cplex_dir/cplex/lib/x86-64_linux/static_pic"
        "$cplex_dir/cplex/lib/x86-64_osx/static_pic"
        "$cplex_dir/cplex/lib/win64/static_pic"
        "$cplex_dir/bin"
        "$cplex_dir/lib"
    )
    
    for lib_dir in "${lib_dirs[@]}"; do
        if [[ -d "$lib_dir" ]]; then
            # Check if it contains CPLEX libraries
            if find "$lib_dir" -name "*cplex*" -o -name "*ilocplex*" | grep -q .; then
                lib_paths+=("$lib_dir")
            fi
        fi
    done
    
    printf '%s\n' "${lib_paths[@]}"
}

# Function to configure CPLEX environment
configure_cplex_environment() {
    local cplex_dir="$1"
    local python_version="$2"
    
    log_info "Configuring CPLEX environment for $cplex_dir"
    
    # Set CPLEX_STUDIO_DIR
    export CPLEX_STUDIO_DIR="$cplex_dir"
    
    # Find and set Python path
    local python_path
    if python_path=$(find_cplex_python_path "$cplex_dir" "$python_version"); then
        log_info "Found CPLEX Python API at: $python_path"
        
        # Add to PYTHONPATH
        if [[ -n "${PYTHONPATH:-}" ]]; then
            export PYTHONPATH="$python_path:$PYTHONPATH"
        else
            export PYTHONPATH="$python_path"
        fi
    else
        log_warn "Could not find CPLEX Python API path"
        return 1
    fi
    
    # Find and set library paths
    local lib_paths
    lib_paths=$(find_cplex_lib_paths "$cplex_dir")
    
    if [[ -n "$lib_paths" ]]; then
        while IFS= read -r lib_path; do
            log_info "Adding CPLEX library path: $lib_path"
            
            # Add to LD_LIBRARY_PATH (Linux)
            if [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
                export LD_LIBRARY_PATH="$lib_path:$LD_LIBRARY_PATH"
            else
                export LD_LIBRARY_PATH="$lib_path"
            fi
            
            # Add to DYLD_LIBRARY_PATH (macOS)
            if [[ -n "${DYLD_LIBRARY_PATH:-}" ]]; then
                export DYLD_LIBRARY_PATH="$lib_path:$DYLD_LIBRARY_PATH"
            else
                export DYLD_LIBRARY_PATH="$lib_path"
            fi
            
        done <<< "$lib_paths"
    else
        log_warn "Could not find CPLEX library paths"
        return 1
    fi
    
    return 0
}

# Function to save environment variables to conda activation script
save_conda_env_vars() {
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        log_warn "Not in a conda environment, cannot save persistent environment variables"
        return 1
    fi
    
    local activate_dir="$CONDA_PREFIX/etc/conda/activate.d"
    local deactivate_dir="$CONDA_PREFIX/etc/conda/deactivate.d"
    local activate_script="$activate_dir/cplex_setup.sh"
    local deactivate_script="$deactivate_dir/cplex_setup.sh"
    
    # Create directories if they don't exist
    mkdir -p "$activate_dir" "$deactivate_dir"
    
    # Write activation script
    cat > "$activate_script" << EOF
#!/usr/bin/env bash
# Auto-generated CPLEX environment setup for CarveWe
# Generated on $(date)

export CPLEX_STUDIO_DIR="${CPLEX_STUDIO_DIR:-}"
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}"
EOF
    
    # Write deactivation script
    cat > "$deactivate_script" << EOF
#!/usr/bin/env bash
# Auto-generated CPLEX environment cleanup for CarveWe

unset CPLEX_STUDIO_DIR
# Note: We don't unset PYTHONPATH, LD_LIBRARY_PATH, DYLD_LIBRARY_PATH
# as they might contain other important paths
EOF
    
    chmod +x "$activate_script" "$deactivate_script"
    
    log_info "Saved conda environment activation scripts:"
    log_info "  - Activate: $activate_script"
    log_info "  - Deactivate: $deactivate_script"
}

# Function to test CPLEX availability
test_cplex_availability() {
    log_info "Testing CPLEX availability..."
    
    # Test Python import
    if python -c "import cplex; print(f'CPLEX version: {cplex.__version__}')" 2>/dev/null; then
        log_info "✓ CPLEX Python API is accessible"
    else
        log_error "✗ CPLEX Python API is not accessible"
        return 1
    fi
    
    # Test optlang integration
    if python -c "import optlang.cplex_interface; print('✓ optlang CPLEX interface is available')" 2>/dev/null; then
        log_info "✓ optlang CPLEX interface is available"
    else
        log_error "✗ optlang CPLEX interface is not available"
        return 1
    fi
    
    return 0
}

# Main execution
main() {
    log_info "Starting CPLEX discovery and configuration..."
    
    # Check if CPLEX is already configured and working
    if is_cplex_configured; then
        log_info "CPLEX appears to already be configured and working"
        if test_cplex_availability; then
            log_info "CPLEX configuration is valid, nothing to do"
            return 0
        else
            log_warn "CPLEX seems configured but not working, reconfiguring..."
        fi
    fi
    
    # Get Python version
    local python_version
    python_version=$(get_python_version)
    log_info "Python version: $python_version"
    
    # Find CPLEX installations
    log_info "Searching for CPLEX installations..."
    local cplex_installations
    cplex_installations=$(find_cplex_installations)
    
    if [[ -z "$cplex_installations" ]]; then
        log_error "No CPLEX installations found!"
        log_error "Please ensure CPLEX is installed and try again."
        log_error "Common installation locations:"
        log_error "  - /opt/ibm/ILOG/CPLEX_Studio*"
        log_error "  - /Applications/CPLEX_Studio* (macOS)"
        log_error "  - \$HOME/CPLEX_Studio*"
        return 1
    fi
    
    # Try each installation until one works
    local configured=false
    while IFS= read -r cplex_dir; do
        log_info "Checking CPLEX installation: $cplex_dir"
        
        if validate_cplex_installation "$cplex_dir"; then
            log_info "Valid CPLEX installation found: $cplex_dir"
            
            if configure_cplex_environment "$cplex_dir" "$python_version"; then
                if test_cplex_availability; then
                    log_info "Successfully configured CPLEX from: $cplex_dir"
                    configured=true
                    
                    # Save to conda environment if possible
                    if save_conda_env_vars; then
                        log_info "Environment variables saved for future activations"
                    fi
                    
                    break
                else
                    log_warn "CPLEX configuration failed for: $cplex_dir"
                fi
            else
                log_warn "Could not configure environment for: $cplex_dir"
            fi
        else
            log_warn "Invalid CPLEX installation: $cplex_dir"
        fi
    done <<< "$cplex_installations"
    
    if [[ "$configured" != true ]]; then
        log_error "Failed to configure any CPLEX installation"
        return 1
    fi
    
    log_info "CPLEX configuration completed successfully!"
    return 0
}

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi