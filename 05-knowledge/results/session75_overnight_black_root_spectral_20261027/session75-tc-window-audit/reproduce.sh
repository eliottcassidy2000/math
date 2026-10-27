#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sha256sum -c SOURCE_SHA256SUMS
python3 window_configs.py jobs-normal > configs-normal.out
python3 -O window_configs.py jobs-optimized > configs-optimized.out
cmp configs-normal.out configs-optimized.out
python3 window_configs_n12.py jobs-n12 > configs-n12.out
python3 -O window_configs_n12.py jobs-n12-opt > configs-n12-opt.out
cmp configs-n12.out configs-n12-opt.out
python3 audit_jobs.py > audit-jobs-normal.json
python3 -O audit_jobs.py > audit-jobs-optimized.json
cmp audit-jobs-normal.json audit-jobs-optimized.json
python3 make_portable.py
g++ -O3 -std=c++17 -o tc_window_scan_portable tc_window_scan_portable.cpp
g++ -O3 -std=c++17 -o test_aut test_aut.cpp
./test_aut > test-aut.out
python3 run_scan.py normal > replay-normal.progress
python3 run_scan.py n12 > replay-n12.progress
python3 summarize.py > SUMMARY.json
python3 -O summarize.py > SUMMARY-optimized.json
cmp SUMMARY.json SUMMARY-optimized.json
cat SUMMARY.json
