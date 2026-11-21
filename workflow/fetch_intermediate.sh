#!/bin/bash
# download and unpack the archived intermediate outputs of steps 1 to 9 (per-group fits, per-cell records, contrasts,
# chromatin and spatial summaries) into $1; the archive location is in config/intermediate.txt
set -euo pipefail
dest=${1:-results}; mkdir -p "$dest"
url=$(grep -v '^#' config/intermediate.txt | head -1)
[ -s "$dest/.intermediate_ok" ] && { echo "intermediate outputs already in $dest"; exit 0; }
tmp=$(mktemp -d)
if [[ "$url" == gh:* ]]; then gh release download "${url#gh:}" -p 'loyscope_intermediate.tar.gz' -D "$tmp"
else curl -sSL --retry 5 -o "$tmp/loyscope_intermediate.tar.gz" "$url"; fi
( cd "$tmp" && sha256sum -c "$OLDPWD/config/intermediate.sha256" 2>/dev/null || shasum -a 256 -c "$OLDPWD/config/intermediate.sha256" )
tar xzf "$tmp/loyscope_intermediate.tar.gz" -C "$dest"; rm -rf "$tmp"; date > "$dest/.intermediate_ok"
