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
        # Standard IBM CPLEX installations
        "/opt/ibm/ILOG/CPLEX_Studio*"
        "/Applications/CPLEX_Studio*"
        "$HOME/CPLEX_Studio*"
        "$HOME/Applications/CPLEX_Studio*"
        "/usr/local/CPLEX_Studio*"
        "/opt/CPLEX_Studio*"
        
        # Conda-based installations (various patterns)
        "${CONDA_PREFIX:-}/cplex*"
        "${CONDA_PREFIX:-}/../cplex*"
        "${CONDA_PREFIX:-}/ilog/cplex*"
        "${CONDA_PREFIX:-}/../ilog/cplex*"
        
        # Alternative naming patterns
        "/opt/cplex*"
        "$HOME/cplex*"
        "/usr/local/cplex*"
        
        # Search in all conda environments for ilog folders
        "${CONDA_PREFIX:-}/../*/ilog/cplex*"
        "${CONDA_PREFIX:-}/../../*/ilog/cplex*"
    )
    
    # Also search for "ilog" directories that might contain CPLEX
    local ilog_search_dirs=(
        "/opt/ilog"
        "$HOME/ilog" 
        "/usr/local/ilog"
        "${CONDA_PREFIX:-}/ilog"
        "${CONDA_PREFIX:-}/../ilog"
        "${CONDA_PREFIX:-}/../../*/ilog"  # Search other conda envs
    )
    
    # Search for ilog directories and look for CPLEX inside them
    for ilog_pattern in "${ilog_search_dirs[@]}"; do
        for ilog_dir in $ilog_pattern; do
            if [[ -d "$ilog_dir" ]]; then
                log_info "Found ilog directory: $ilog_dir"
                # Look for CPLEX installations inside ilog directory
                for cplex_subdir in "$ilog_dir"/cplex* "$ilog_dir"/CPLEX*; do
                    if [[ -d "$cplex_subdir" ]]; then
                        log_info "Found potential CPLEX in ilog: $cplex_subdir"
                        cplex_paths+=("$cplex_subdir")
                    fi
                done
            fi
        done 2>/dev/null || true
    done
    
    # Also check environment variables that might point to CPLEX
    if [[ -n "${CPLEX_STUDIO_DIR:-}" ]] && [[ -d "$CPLEX_STUDIO_DIR" ]]; then
        cplex_paths+=("$CPLEX_STUDIO_DIR")
    fi
    
    if [[ -n "${ILOG_CPLEX_PATH:-}" ]] && [[ -d "$ILOG_CPLEX_PATH" ]]; then
        cplex_paths+=("$ILOG_CPLEX_PATH")
    fi
    
    # Search for CPLEX installations using the main patterns
    for pattern in "${search_dirs[@]}"; do
        # Use shell expansion to find matching directories
        for dir in $pattern; do
            if [[ -d "$dir" ]]; then
                log_info "Found potential CPLEX directory: $dir"
                cplex_paths+=("$dir")
            fi
        done 2>/dev/null || true
    done
    
    # Remove duplicates and sort by version (newest first)
    printf '%s\n' "${cplex_paths[@]}" | sort -u -V -r
}

# Function to validate CPLEX installation
validate_cplex_installation() {
    local cplex_dir="$1"
    
    log_info "Validating CPLEX installation: $cplex_dir"
    
    # Check for various CPLEX directory structures
    local valid_structure=false
    
    # Structure 1: Traditional CPLEX_Studio layout
    if [[ -d "$cplex_dir/cplex" ]] && [[ -d "$cplex_dir/cplex/bin" ]]; then
        log_info "Found traditional CPLEX Studio structure"
        valid_structure=true
    fi
    
    # Structure 2: Direct CPLEX installation (like cplex_22.1.2)
    if [[ -d "$cplex_dir/bin" ]] && [[ -d "$cplex_dir/lib" ]]; then
        # Look for CPLEX-specific files
        if find "$cplex_dir" -name "*cplex*" -o -name "*ilocplex*" | grep -q .; then
            log_info "Found direct CPLEX installation structure"
            valid_structure=true
        fi
    fi
    
    # Structure 3: Look for Python API indicators
    local python_dirs=(
        "$cplex_dir/python"
        "$cplex_dir/cplex/python"
        "$cplex_dir/lib/python"
    )
    
    for python_dir in "${python_dirs[@]}"; do
        if [[ -d "$python_dir" ]]; then
            if find "$python_dir" -name "setup.py" -o -name "cplex" -type d | grep -q .; then
                log_info "Found Python API in: $python_dir"
                valid_structure=true
                break
            fi
        fi
    done
    
    if [[ "$valid_structure" != true ]]; then
        log_warn "CPLEX installation at $cplex_dir does not match expected structure"
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
    
    # Common Python API locations for different installation types
    local python_paths=(
        # Traditional CPLEX Studio structure
        "$cplex_dir/cplex/python/${python_version}/x86-64_linux"
        "$cplex_dir/cplex/python/${python_version}/x86-64_osx"
        "$cplex_dir/cplex/python/${python_version}/win64"
        "$cplex_dir/cplex/python/${python_version}"
        "$cplex_dir/cplex/python"
        
        # Direct CPLEX installation structure  
        "$cplex_dir/python/${python_version}/x86-64_linux"
        "$cplex_dir/python/${python_version}/x86-64_osx"
        "$cplex_dir/python/${python_version}/win64"
        "$cplex_dir/python/${python_version}"
        "$cplex_dir/python"
        
        # Alternative locations
        "$cplex_dir/lib/python/${python_version}"
        "$cplex_dir/lib/python"
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
    
    # Common library locations for different installation types
    local lib_dirs=(
        # Traditional CPLEX Studio structure
        "$cplex_dir/cplex/bin/x86-64_linux"
        "$cplex_dir/cplex/bin/x86-64_osx" 
        "$cplex_dir/cplex/bin/win64"
        "$cplex_dir/cplex/lib/x86-64_linux/static_pic"
        "$cplex_dir/cplex/lib/x86-64_osx/static_pic"
        "$cplex_dir/cplex/lib/win64/static_pic"
        
        # Direct CPLEX installation structure
        "$cplex_dir/bin/x86-64_linux"
        "$cplex_dir/bin/x86-64_osx"
        "$cplex_dir/bin/win64"
        "$cplex_dir/bin"
        "$cplex_dir/lib/x86-64_linux"
        "$cplex_dir/lib/x86-64_osx"
        "$cplex_dir/lib/win64"
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

# Function to detect CPLEX license type and validate sufficiency
detect_cplex_license_type() {
    local cplex_dir="$1"
    
    log_info "Detecting CPLEX license type and capabilities..."
    
    # Try to detect if this is the limited Community Edition
    # Community Edition is typically distributed as "CPLEX Studio" with specific version patterns
    local is_community_edition=false
    
    # Check for Community Edition indicators
    if [[ "$cplex_dir" == *"CPLEX_Studio"* ]]; then
        # Look for Community Edition specific files or version indicators
        if [[ -f "$cplex_dir/cplex/include/ilcplex/cplex.h" ]]; then
            # Try to extract version info from header file
            local version_info
            version_info=$(grep -i "community\|student\|academic\|limited" "$cplex_dir/cplex/include/ilcplex/cplex.h" 2>/dev/null || true)
            if [[ -n "$version_info" ]]; then
                log_warn "Detected Community/Student/Limited edition indicators"
                is_community_edition=true
            fi
        fi
        
        # Check version number patterns - Community editions often have specific version ranges
        local dir_name=$(basename "$cplex_dir")
        if [[ "$dir_name" =~ CPLEX_Studio([0-9]+) ]]; then
            local version_num="${BASH_REMATCH[1]}"
            log_info "Detected CPLEX Studio version: $version_num"
            
            # Community editions are typically distributed as standalone "Studio" packages
            # while full versions are part of larger ILOG/IBM packages
            if [[ ! -d "$(dirname "$cplex_dir")/opl" ]] && [[ ! -d "$(dirname "$cplex_dir")/concert" ]]; then
                log_warn "Installation appears to be standalone CPLEX Studio (possibly Community Edition)"
                is_community_edition=true
            fi
        fi
    fi
    
    # Try runtime detection by attempting to solve a test problem
    if command -v python >/dev/null 2>&1; then
        log_info "Performing runtime license validation..."
        
        local test_result
        test_result=$(python -c "
import sys
try:
    import cplex
    # Create a test problem with more than 1000 variables (Community Edition limit)
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)
    
    # Add 1001 variables to test the limit
    var_names = [f'x{i}' for i in range(1001)]
    prob.variables.add(names=var_names, lb=[0.0]*1001, ub=[1.0]*1001)
    
    # Try to solve - this will fail on Community Edition
    prob.solve()
    print('FULL_LICENSE')
except Exception as e:
    error_msg = str(e).lower()
    if 'size limit' in error_msg or 'license' in error_msg or 'limit' in error_msg:
        print('COMMUNITY_EDITION')
    else:
        print('UNKNOWN_ERROR')
except ImportError:
    print('IMPORT_ERROR')
" 2>/dev/null || echo "PYTHON_ERROR")
        
        case "$test_result" in
            "COMMUNITY_EDITION")
                log_error "CPLEX Community Edition detected (1000 variable limit)"
                log_error "CarveWe requires a full CPLEX license for metabolic network reconstruction"
                log_error "Please obtain an academic or commercial CPLEX license"
                return 1
                ;;
            "FULL_LICENSE")
                log_info "✓ Full CPLEX license detected - no size limitations"
                return 0
                ;;
            "UNKNOWN_ERROR"|"PYTHON_ERROR"|"IMPORT_ERROR")
                log_warn "Could not determine CPLEX license type through runtime testing"
                if [[ "$is_community_edition" == true ]]; then
                    log_error "Installation appears to be Community Edition based on directory structure"
                    log_error "CarveWe requires a full CPLEX license for metabolic network reconstruction"
                    return 1
                else
                    log_warn "Proceeding with configuration, but license validation failed"
                    return 0
                fi
                ;;
        esac
    else
        log_warn "Python not available for license testing"
        if [[ "$is_community_edition" == true ]]; then
            log_error "Installation appears to be Community Edition, and runtime validation unavailable"
            log_error "CarveWe requires a full CPLEX license for metabolic network reconstruction"
            return 1
        fi
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
            
            # Check if license is sufficient for CarveWe
            if ! detect_cplex_license_type "$cplex_dir"; then
                log_warn "CPLEX installation has insufficient license for CarveWe: $cplex_dir"
                continue  # Try next installation
            fi
            
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
        log_error "Failed to configure any suitable CPLEX installation"
        log_error ""
        log_error "CarveWe requires a full CPLEX license (Academic or Commercial)"
        log_error "The free Community Edition (1000 variable limit) is insufficient"
        log_error ""
        log_error "To obtain CPLEX:"
        log_error "  Academic users: https://www.ibm.com/academic/home"  
        log_error "  Commercial users: Contact IBM for licensing"
        log_error "  Alternative: Use other solvers (SCIP, Gurobi) if supported"
        return 1
    fi
    
    log_info "CPLEX configuration completed successfully!"
    return 0
}

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi