#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
g++ -std=c++17 -O2 -Wall -Wextra forest_audit.cpp -o forest_audit
./forest_audit 9 > forest_audit.out
python3 independent_check.py > independent-normal.json
python3 -O independent_check.py > independent-optimized.json
cmp independent-normal.json independent-optimized.json
printf 'Independent audit PASS; normal and optimized JSON agree.\n'
