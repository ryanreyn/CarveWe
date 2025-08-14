#!/usr/bin/env bash
set -euxo pipefail

# Entrypoints are available?
command -v carvewe
command -v genome-aligner

# Help should succeed quickly
CarveWe.sh --help >/dev/null
genome-aligner.sh --help >/dev/null

# Data dir (optional)
test -d "${CARVEWE_SHARE:-$PREFIX/share/carvewe}" || true

echo "Smoke tests passed."
