#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sha256sum -c SHA256SUMS >/dev/null
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 -B check_addendum.py > "$T/normal.json"
python3 -B -O check_addendum.py > "$T/optimized.json"
cmp "$T/normal.json" normal.json
cmp "$T/optimized.json" optimized.json
cmp "$T/normal.json" "$T/optimized.json"
printf 'PASS: bounded two-object cycle profile and pinned incoming E10 isomorphism.\n'
