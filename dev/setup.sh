#!/usr/bin/env bash
# Setup script for CarveWe development environment

# Get the directory where this script lives and navigate to project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Set up development environment variables
export CARVEWE_INSTALL_DIR="${PROJECT_ROOT}"
export CARVEWE_SCRIPTS_DIR="${PROJECT_ROOT}/bin"
export CARVEWE_DATA_DIR="${PROJECT_ROOT}/ref_data"

# Add scripts to PATH for development
export PATH="${CARVEWE_SCRIPTS_DIR}:${PATH}"

echo "CarveWe development environment set up:"
echo "CARVEWE_INSTALL_DIR=${CARVEWE_INSTALL_DIR}"
echo "CARVEWE_SCRIPTS_DIR=${CARVEWE_SCRIPTS_DIR}"
echo "CARVEWE_DATA_DIR=${CARVEWE_DATA_DIR}"
echo ""
echo "To use:"
echo "source dev-setup.sh"