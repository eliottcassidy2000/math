#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_fifteen.py > normal.json
python3 -O verify_fifteen.py > optimized.json
cmp normal.json optimized.json
printf 'C2 independent fifteen-move audit PASS.\n'
