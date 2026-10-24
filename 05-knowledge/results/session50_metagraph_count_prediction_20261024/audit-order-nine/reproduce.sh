#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
normal_bin=$(mktemp /tmp/order-nine-audit-normal.XXXXXX)
optimized_bin=$(mktemp /tmp/order-nine-audit-optimized.XXXXXX)
trap 'rm -f "$normal_bin" "$optimized_bin"' EXIT
g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit.cpp -o "$normal_bin"
"$normal_bin" archive/merged_nodes.tsv.txt compiled-normal > compiled-normal.out
g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -pedantic audit.cpp -o "$optimized_bin"
"$optimized_bin" archive/merged_nodes.tsv.txt compiled-optimized > compiled-optimized.out
cmp compiled-normal.json compiled-optimized.json
cmp compiled-normal-parents.tsv compiled-optimized-parents.tsv
cmp compiled-normal-regular.tsv compiled-optimized-regular.tsv
PYTHONDONTWRITEBYTECODE=1 python3 check.py > normal.out
PYTHONDONTWRITEBYTECODE=1 python3 -O check.py --output optimized.json > optimized.out
cmp normal.json optimized.json
python3 verify_manifest.py
