#!/usr/bin/env bash
set -euxo pipefail

mkdir -p "$PREFIX/bin"

# 1) Install all scripts from bin/ preserving exec bits
if [ -d "bin" ]; then
  find bin -maxdepth 1 -type f -exec install -m 0755 {} "$PREFIX/bin/" \;
fi

# 2) Optional: install reference data
if [ -d "ref_data" ]; then
  mkdir -p "$PREFIX/share/carvewe"
  cp -R ref_data/* "$PREFIX/share/carvewe"/
fi

# 3) Create friendly CLI names (symlink, with portable shim fallback)
make_launcher () {
  local target="$1"   # e.g., CarveWe.sh
  local link="$2"     # e.g., carvewe
  (
    cd "$PREFIX/bin"
    # Prefer a relative symlink for relocatability
    if ln -s "$target" "$link" 2>/dev/null; then
      :
    else
      # Fallback: write a tiny shim that execs the real script next to it
      cat > "$link" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
exec "$HERE/__TARGET__" "$@"
EOF
      # Replace placeholder with the actual target filename
      sed -i.bak "s|__TARGET__|$target|g" "$link" || \
      perl -0777 -pe "s/__TARGET__/$target/g" -i "$link"
      rm -f "$link.bak"
      chmod 0755 "$link"
    fi
  )
}

# Ensure the actual entry scripts exist (installed above)
test -x "$PREFIX/bin/CarveWe.sh"
test -x "$PREFIX/bin/genome-aligner.sh"

# Create user-facing command names (no .sh)
make_launcher "CarveWe.sh" "carvewe"
make_launcher "genome-aligner.sh" "genome-aligner"

# 4) Do NOT install publication repro scripts
echo "Skipping reproduce_publication/ (not part of end-user package)."
