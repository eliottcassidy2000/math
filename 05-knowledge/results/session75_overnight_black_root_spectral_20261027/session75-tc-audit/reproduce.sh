#!/bin/sh
set -eu
cd "$(dirname "$0")"
python analytic_audit.py > analytic.normal.json
python -O analytic_audit.py > analytic.optimized.json
cmp analytic.normal.json analytic.optimized.json
python group_checks.py > group_checks.normal.json
python -O group_checks.py > group_checks.optimized.json
cmp group_checks.normal.json group_checks.optimized.json
g++ -std=c++17 -O3 -Wall -Wextra -Werror cayley_independent.cpp -o cayley_independent
./cayley_independent > cayley_independent.out
./cayley_independent input < vt15_independent_inputs.txt > vt15_independent.out
g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Werror cayley_independent.cpp -o cayley_independent.ndebug
./cayley_independent.ndebug > cayley_independent.ndebug.out
./cayley_independent.ndebug input < vt15_independent_inputs.txt > vt15_independent.ndebug.out
cmp cayley_independent.out cayley_independent.ndebug.out
cmp vt15_independent.out vt15_independent.ndebug.out
python 04-computation/experiments/tc_fibres_n12_20261027/vt15_noncayley.py > vt15.normal.out
python -O 04-computation/experiments/tc_fibres_n12_20261027/vt15_noncayley.py > vt15.optimized.out
cmp vt15.normal.out vt15.optimized.out
python verify_receipts.py
