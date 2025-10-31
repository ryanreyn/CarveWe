#!/usr/bin/env bash
set -euxo pipefail

mkdir -p "$PREFIX/bin"

# 1) Install entry point scripts to a private location
mkdir -p "$PREFIX/libexec/carvewe"
if [ -d "bin" ]; then
    find bin -maxdepth 1 -type f -exec install -m 0755 {} "$PREFIX/libexec/carvewe/" \;
fi

# 2) Install reference data
if [ -d "ref_data" ]; then
    mkdir -p "$PREFIX/share/carvewe"
    cp -R ref_data/* "$PREFIX/share/carvewe"/
fi

# Install documentation
mkdir -p "$PREFIX/share/carvewe/docs"
if [ -f "docs/CPLEX_SETUP.md" ]; then
    cp docs/CPLEX_SETUP.md "$PREFIX/share/carvewe/docs/"
fi

# 3) Install CPLEX environment configuration
mkdir -p "$PREFIX/etc/conda/activate.d"
mkdir -p "$PREFIX/etc/conda/deactivate.d"

# Create activation script for CPLEX
cat > "$PREFIX/etc/conda/activate.d/carvewe_cplex.sh" << 'EOF'
#!/usr/bin/env bash
# CarveWe CPLEX Configuration
# This script sets up CPLEX environment variables if CPLEX is installed

# Only proceed if ilog directory exists
if [ -d "$CONDA_PREFIX/ilog" ]; then
    # Find CPLEX installation (handle multiple versions, use newest)
    CPLEX_VERSION_DIR=$(ls -1 "$CONDA_PREFIX/ilog" | grep "^cplex_" | sort -V | tail -n 1)
    
    if [ -n "$CPLEX_VERSION_DIR" ]; then
        export CPLEX_HOME="$CONDA_PREFIX/ilog/$CPLEX_VERSION_DIR/cplex"
        
        # Add CPLEX libraries to library path
        if [ -d "$CPLEX_HOME/bin/x86-64_linux" ]; then
            export LD_LIBRARY_PATH="$CPLEX_HOME/bin/x86-64_linux:${LD_LIBRARY_PATH:-}"
        fi
        
        # Set CPLEX_STUDIO_DIR for other tools that may need it
        export CPLEX_STUDIO_DIR="$CONDA_PREFIX/ilog/$CPLEX_VERSION_DIR"
        
        # Mark that we configured CPLEX (for deactivation)
        export _CARVEWE_CPLEX_CONFIGURED=1
    fi
fi
EOF
chmod 0755 "$PREFIX/etc/conda/activate.d/carvewe_cplex.sh"

# Create deactivation script
cat > "$PREFIX/etc/conda/deactivate.d/carvewe_cplex.sh" << 'EOF'
#!/usr/bin/env bash
# CarveWe CPLEX Deactivation
# Clean up CPLEX environment variables

if [ -n "${_CARVEWE_CPLEX_CONFIGURED:-}" ]; then
    unset CPLEX_HOME
    unset CPLEX_STUDIO_DIR
    unset _CARVEWE_CPLEX_CONFIGURED
    # Note: We don't unset LD_LIBRARY_PATH as other tools may use it
    # It will be reset when the environment is deactivated
fi
EOF
chmod 0755 "$PREFIX/etc/conda/deactivate.d/carvewe_cplex.sh"

# 4) Install CPLEX diagnostic tool
cat > "$PREFIX/bin/carvewe-check-cplex" << 'EOF'
#!/usr/bin/env python
"""Check CPLEX installation and configuration for CarveWe"""

import sys
import os
from pathlib import Path

def check_cplex():
    print("=" * 60)
    print("CarveWe CPLEX Configuration Check")
    print("=" * 60)
    
    # Check environment variables
    conda_prefix = os.environ.get('CONDA_PREFIX')
    cplex_home = os.environ.get('CPLEX_HOME')
    cplex_studio = os.environ.get('CPLEX_STUDIO_DIR')
    ld_library_path = os.environ.get('LD_LIBRARY_PATH', '')
    
    print(f"\n1. Environment Variables:")
    print(f"   CONDA_PREFIX: {conda_prefix or 'NOT SET'}")
    print(f"   CPLEX_HOME: {cplex_home or 'NOT SET'}")
    print(f"   CPLEX_STUDIO_DIR: {cplex_studio or 'NOT SET'}")
    
    # Check for ilog installation
    print(f"\n2. CPLEX Binary Installation:")
    if conda_prefix:
        ilog_path = Path(conda_prefix) / 'ilog'
        if ilog_path.exists():
            cplex_versions = list(ilog_path.glob('cplex_*'))
            if cplex_versions:
                print(f"   ✓ Found CPLEX at: {ilog_path}")
                for v in cplex_versions:
                    print(f"     - {v.name}")
            else:
                print(f"   ✗ ilog directory exists but no CPLEX version found")
        else:
            print(f"   ✗ No ilog directory at: {ilog_path}")
            print(f"      Install CPLEX binaries with:")
            print(f"      ./cplex_studio_installer.bin -i console -DUSER_INSTALL_DIR=$CONDA_PREFIX/ilog")
    
    # Check pip cplex package
    print(f"\n3. CPLEX Python API:")
    try:
        import cplex
        print(f"   ✓ cplex module found: {cplex.__file__}")
        
        # Test if it's full or community edition
        try:
            c = cplex.Cplex()
            c.variables.add(names=[f'x{i}' for i in range(2000)])
            c.solve()
            print(f"   ✓ FULL CPLEX (2000 variable test passed)")
            print(f"   ✓ Academic license active!")
        except Exception as e:
            if 'size' in str(e).lower() or 'limit' in str(e).lower():
                print(f"   ✗ COMMUNITY EDITION detected (limited to 1000 vars)")
                print(f"      This will NOT work for CarveWe models!")
            else:
                print(f"   ? Test failed with: {e}")
    except ImportError:
        print(f"   ✗ cplex module not found")
        print(f"      Install with: pip install cplex")
    
    # Check library linking
    print(f"\n4. Library Linking:")
    if 'cplex' in sys.modules:
        try:
            import cplex
            # Find the shared library
            cplex_internal = Path(cplex.__file__).parent / '_internal'
            if cplex_internal.exists():
                so_files = list(cplex_internal.glob('*.so'))
                if so_files:
                    print(f"   ✓ CPLEX shared libraries found")
                    # Check if LD_LIBRARY_PATH includes CPLEX
                    if cplex_home and cplex_home in ld_library_path:
                        print(f"   ✓ LD_LIBRARY_PATH includes CPLEX libraries")
                    elif cplex_home:
                        print(f"   ⚠ LD_LIBRARY_PATH may not include CPLEX")
                        print(f"      Try: conda deactivate && conda activate {Path(conda_prefix).name}")
        except Exception as e:
            print(f"   ? Could not verify linking: {e}")
    
    # Check COBRApy integration
    print(f"\n5. COBRApy Integration:")
    try:
        import cobra
        print(f"   ✓ COBRApy installed")
        try:
            from cobra.core import Configuration
            solvers = Configuration().solvers
            print(f"   Available solvers: {', '.join(solvers)}")
            if 'cplex' in solvers:
                print(f"   ✓ CPLEX available as solver")
            else:
                print(f"   ✗ CPLEX not available as solver")
        except:
            print(f"   ? Could not check solver configuration")
    except ImportError:
        print(f"   ✗ COBRApy not installed")
    
    print(f"\n" + "=" * 60)
    print("Next Steps:")
    if not conda_prefix:
        print("  - Activate your conda environment")
    elif not (conda_prefix and Path(conda_prefix, 'ilog').exists()):
        print("  1. Install CPLEX binaries from IBM Academic Initiative")
        print("  2. Run: ./cplex_studio_installer.bin -i console \\")
        print(f"          -DUSER_INSTALL_DIR=$CONDA_PREFIX/ilog")
        print("  3. Run: pip install cplex")
        print("  4. Reactivate environment: conda deactivate && conda activate <env>")
    elif 'cplex' not in sys.modules:
        print("  1. Install: pip install cplex")
        print("  2. Reactivate environment")
    else:
        print("  ✓ CPLEX appears to be configured correctly!")
    print("=" * 60)

    # Point to local docs
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        doc_path = Path(conda_prefix) / "share" / "carvewe" / "docs" / "CPLEX_SETUP.md"
        if doc_path.exists():
            print(f"  Local guide: {doc_path}")
    
    # Always show GitHub link
    print(f"  Online guide: https://github.com/ryanreyn/CarveWe/blob/main/docs/CPLEX_SETUP.md")
    print("=" * 60)

if __name__ == '__main__':
    check_cplex()
EOF
chmod 0755 "$PREFIX/bin/carvewe-check-cplex"

# Create a command to view the cplex installation doc
cat > "$PREFIX/bin/carvewe-docs" << 'EOF'
#!/usr/bin/env bash
# Show CarveWe documentation

DOCS_DIR="$(dirname "$(dirname "$0")")/share/carvewe/docs"

if [ -f "$DOCS_DIR/CPLEX_SETUP.md" ]; then
    if command -v less >/dev/null 2>&1; then
        less "$DOCS_DIR/CPLEX_SETUP.md"
    elif command -v more >/dev/null 2>&1; then
        more "$DOCS_DIR/CPLEX_SETUP.md"
    else
        cat "$DOCS_DIR/CPLEX_SETUP.md"
    fi
else
    echo "Documentation not found locally."
    echo "Visit: https://github.com/ryanreyn/CarveWe/blob/main/docs/CPLEX_SETUP.md"
fi
EOF
chmod 0755 "$PREFIX/bin/carvewe-docs"

# 5) Create user-facing CLI commands
make_launcher () {
    local target="$1"
    local link="$2"
    (
        cd "$PREFIX/bin"
        cat > "$link" <<'LAUNCHER_EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
LIBEXEC_DIR="$(dirname "$HERE")/libexec/carvewe"
exec "$LIBEXEC_DIR/__TARGET__" "$@"
LAUNCHER_EOF
        sed -i.bak "s|__TARGET__|$target|g" "$link" || \
            perl -0777 -pe "s/__TARGET__/$target/g" -i "$link"
        rm -f "$link.bak"
        chmod 0755 "$link"
    )
}

test -x "$PREFIX/libexec/carvewe/CarveWe.sh"
test -x "$PREFIX/libexec/carvewe/genome-aligner.sh"

make_launcher "CarveWe.sh" "carvewe"
make_launcher "genome-aligner.sh" "genome-aligner"

echo "CarveWe installation complete!"
echo "To configure CPLEX, run: carvewe-check-cplex"