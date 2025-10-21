#!/usr/bin/env bash
set -euxo pipefail

mkdir -p "$PREFIX/bin"

# 1) Install entry point scripts to a private location (not directly in user's PATH)
mkdir -p "$PREFIX/libexec/carvewe"
if [ -d "bin" ]; then
    find bin -maxdepth 1 -type f -exec install -m 0755 {} "$PREFIX/libexec/carvewe/" \;
fi

# 2) Optional: install reference data
if [ -d "ref_data" ]; then
    mkdir -p "$PREFIX/share/carvewe"
    cp -R ref_data/* "$PREFIX/share/carvewe"/
fi

# >>> NEW: Install CPLEX activation script for conda environment integration
# Create conda activation directory structure
mkdir -p "$PREFIX/etc/conda/activate.d"
mkdir -p "$PREFIX/etc/conda/deactivate.d"

# Install the CPLEX discovery and activation script
if [ -f "scripts/detect_cplex.sh" ]; then
    # Install the main script as a standalone command users can run manually
    install -m 0755 scripts/detect_cplex.sh "$PREFIX/bin/setup-cplex"
    
    # Create conda activation hook that runs the setup
    cat > "$PREFIX/etc/conda/activate.d/carvewe_cplex.sh" << 'EOF'
#!/usr/bin/env bash
# CarveWe CPLEX auto-configuration
# This runs the CPLEX discovery script on first activation if not already configured

# Only run if CPLEX isn't already working and we're in interactive mode
if [[ -z "${CPLEX_STUDIO_DIR:-}" ]] && [[ $- == *i* ]] 2>/dev/null; then
    # Check if CPLEX is needed (COBRApy with optlang available)
    if python -c "import cobra, optlang" 2>/dev/null; then
        # Check if CPLEX is already accessible
        if ! python -c "import optlang.cplex_interface" 2>/dev/null; then
            echo "CarveWe: CPLEX not detected, attempting auto-configuration..."
            if command -v setup-cplex >/dev/null 2>&1; then
                setup-cplex || echo "CarveWe: CPLEX auto-configuration failed. Run 'setup-cplex' manually if needed."
            fi
        fi
    fi
fi
EOF
    chmod 0755 "$PREFIX/etc/conda/activate.d/carvewe_cplex.sh"
    
    # Create minimal deactivation script
    cat > "$PREFIX/etc/conda/deactivate.d/carvewe_cplex.sh" << 'EOF'
#!/usr/bin/env bash
# CarveWe CPLEX deactivation
# Note: We don't unset environment variables here as they're managed
# by the cplex_setup.sh scripts created by the activation process
EOF
    chmod 0755 "$PREFIX/etc/conda/deactivate.d/carvewe_cplex.sh"
    
    echo "CPLEX activation script installed successfully"
else
    echo "WARNING: CPLEX activation script not found at scripts/detect_cplex.sh"
    echo "Users will need to manually configure CPLEX integration"
fi
# <<< END NEW CPLEX INTEGRATION

# 3) Create user-facing CLI commands (clean names without .sh extensions)
make_launcher () {
    local target="$1"  # e.g., CarveWe.sh
    local link="$2"    # e.g., carvewe

    (
        cd "$PREFIX/bin"
        # Create a shim that calls the actual script in libexec
        cat > "$link" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
LIBEXEC_DIR="$(dirname "$HERE")/libexec/carvewe"
exec "$LIBEXEC_DIR/__TARGET__" "$@"
EOF
        # Replace placeholder with the actual target filename
        sed -i.bak "s|__TARGET__|$target|g" "$link" || \
            perl -0777 -pe "s/__TARGET__/$target/g" -i "$link"
        rm -f "$link.bak"
        chmod 0755 "$link"
    )
}

# Ensure the actual entry scripts exist in libexec (installed above)
test -x "$PREFIX/libexec/carvewe/CarveWe.sh"
test -x "$PREFIX/libexec/carvewe/genome-aligner.sh"

# Create user-facing command names (clean names without .sh extensions)
make_launcher "CarveWe.sh" "carvewe"
make_launcher "genome-aligner.sh" "genome-aligner"

# 4) Do NOT install publication repro scripts
echo "Skipping reproduce_publication/ (not part of end-user package)."