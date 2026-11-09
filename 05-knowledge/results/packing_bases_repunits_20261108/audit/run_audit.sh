#!/usr/bin/env bash
# run_audit.sh - reproduces every audit computation for opus-2026-11-08-S86 (THM-4940 / HYP-9530).
# Usage: bash run_audit.sh WORKDIR      (WORKDIR: any scratch directory; about 2 minutes in total)
# Needs gcc and Python 3 (standard library only).  Uses kind-pasteur's classgen.c (independent prior tool) for the
# class lists of orders <= 9; everything else is the audit's own code.  No session script is imported or copied.
set -euo pipefail
W=${1:?usage: run_audit.sh WORKDIR}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W/cls"
gcc -O2 -o "$W/classgen" "$HERE/../../first_creators_genealogy_20261107/classgen.c"
gcc -O2 -o "$W/aud_fewcycle" "$HERE/aud_fewcycle.c"
gcc -O2 -o "$W/aud_twins_search" "$HERE/aud_twins_search.c"
"$W/classgen" 9 "$W/cls/cls"                                  # A000568 classes of orders 1..9
for n in 1 2 3 4 5 6 7 8 9; do                                # own recomputation of every class (OCF, strong, SC)
  "$W/aud_fewcycle" 0 2 "$n" "$W/cls/cls_$n.txt" > "$W/census_$n.txt"
done 2> "$HERE/aud_census.log"
"$W/aud_fewcycle" 80 0 13 > "$W/fc80.txt" 2> "$HERE/aud_fewcycle_80.log"            # strong, alpha_1 <= 80
"$W/aud_fewcycle" 1023 1 12 > "$W/ic1023.txt" 2> "$HERE/aud_intercyclic_1023.log"   # strong intercyclic, < 1024
"$W/aud_fewcycle" 8191 1 0 > "$W/ic8191.txt" 2> "$HERE/aud_intercyclic_8191.log"    # strong intercyclic, < 8192
"$W/aud_fewcycle" 0 3 20 "$W/fc80.txt" > "$HERE/aud_invariance.out"                 # canonical-form invariance
"$W/aud_fewcycle" 0 3 20 "$W/ic1023.txt" >> "$HERE/aud_invariance.out"
"$W/aud_fewcycle" 0 4 8 "$W" > "$HERE/aud_fibre_le8.out"                            # same-H move components
python "$HERE/aud_spectra.py" "$W/fc80.txt" "$W/ic1023.txt" > "$HERE/aud_spectra.out"
python "$HERE/aud_compare.py" "$W" "$HERE/.." "$HERE/../../session224_first_occurrence_20261107/census/birth_fibers.csv" > "$HERE/aud_compare.out"
python "$HERE/aud_single_apex.py" > "$HERE/aud_single_apex.out"
# carried twins for the 38 carry-free singletons in [1024, 8192) (randomised search; witnesses re-verified)
python "$HERE/aud_singles_hi.py" "$W/ic8191.txt" 1024 7845 > "$W/singles_hi.txt"   # 37 values at order 10
"$W/aud_twins_search" 10 300 < "$W/singles_hi.txt" > "$HERE/aud_twins_hi.out"
python "$HERE/aud_singles_hi.py" "$W/ic8191.txt" 7845 8192 > "$W/singles_top.txt"   # a = 7845: H = 15691 needs order 11
"$W/aud_twins_search" 11 300 < "$W/singles_top.txt" >> "$HERE/aud_twins_hi.out"
# follow-up audit of the post-checkpoint commits (PREREGISTRATION-2, coverage_products, ternary_attainment)
python "$HERE/aud_followup.py" "$W" "$HERE/.." > "$HERE/aud_followup.out"
echo "done"
