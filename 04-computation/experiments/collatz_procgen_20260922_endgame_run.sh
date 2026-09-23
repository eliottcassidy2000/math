#!/usr/bin/env bash
# Reproduce 05-knowledge/results/collatz_procgen_20260922_endgame.out  (Q2 endgame lane, session collatz-procgen-20260922)
#
#   bash 04-computation/experiments/collatz_procgen_20260922_endgame_run.sh > 05-knowledge/results/collatz_procgen_20260922_endgame.out
#   FULL=1 bash ...   additionally: scan of the 1/2 class to 10^11 (3.5 min) and the exhaustive Psi-DFS to 10^18
#                     (about 30 min, 1 MB); FULL=2 also runs the DFS to 10^19 (several hours).
#   DUMP=<path to the dimension lane's bwd_bad_r41.txt>  recomputes the 98-point census from the dump.
# Quick mode: about 1.5 minutes on one core, peak memory < 300 MB (the exact prover on the census is the slow part).
set -euo pipefail
cd "$(dirname "$0")"
P=collatz_procgen_20260922
T=$(mktemp -d)
clang -O3 -march=native -o "$T/psi" ${P}_endgame_psi.c -lm

echo "### 1. exact chain calculus: Psi, budget constant, explicit escapes, frontier recursion, thread table"
python3 ${P}_endgame_chain.py 80 ${DUMP:-}
echo
echo "### 2. exact hostility prover: census, Theorem F_b families, T^-1(2) points, Psi-candidates"
python3 -u ${P}_endgame_hostile.py 22 24 GHIJ
echo
echo "### 3. Diophantine inputs (LTE, Yu, Senge-Straus/Stewart, digit-exhausted family) and renewal structure"
python3 ${P}_endgame_diophantine.py
echo
echo "### 4. Psi-orbit numerics (128-bit)"
python3 -c "
import sys; sys.path.insert(0,'.')
from ${P}_endgame_chain import census_points
for x in census_points():
    d=x.denominator; e=d.bit_length()-1; print(x.numerator, e)
" > "$T/points.txt"
echo "--- 4a. every m = 14 (mod 27) up to 10^9"
"$T/psi" scan 1000000000
echo "--- 4b. random m = (p + 3^j W)/2^e <= 10^18 near each of the 98 census threads, 2000 per (thread, j)"
"$T/psi" sample "$T/points.txt" 1000000000000000000 2000 20260922 | tail -8
echo "--- 4c. exhaustive DFS over Psi-alive classes, m <= 10^13 (quick)"
"$T/psi" dfs 27 10000000000000
echo "--- 4d. extreme orbits"
for m in 6082250 60091390742 19296859109 909072835151420597 829812238225934675 962019445183081187 5728828425613736 4847486; do "$T/psi" orbit $m; done
if [ "${FULL:-0}" != "0" ]; then
  echo "--- 4e. FULL: every m = 14 (mod 27) up to 10^11"
  "$T/psi" scan 100000000000
  echo "--- 4f. FULL: exhaustive DFS over Psi-alive classes mod 3^38, m <= 10^18"
  "$T/psi" dfs 38 1000000000000000000
fi
if [ "${FULL:-0}" = "2" ]; then
  echo "--- 4g. FULL=2: exhaustive DFS over Psi-alive classes mod 3^40, m <= 10^19"
  "$T/psi" dfs 40 10000000000000000000
fi
rm -rf "$T"
