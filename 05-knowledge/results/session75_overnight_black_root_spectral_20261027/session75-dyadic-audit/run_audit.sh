#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT
python3 independent_verify.py --producer producer_revision2 > "$temporary/normal.json"
python3 -O independent_verify.py --producer producer_revision2 > "$temporary/optimized.json"
cmp "$temporary/normal.json" "$temporary/optimized.json"
cmp "$temporary/normal.json" receipt.compared.json
echo 'PASS: independent normal/-O exact reconstruction, frozen revision2 comparison, byte-identical receipt'
