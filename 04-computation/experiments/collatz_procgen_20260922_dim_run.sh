#!/usr/bin/env bash
# Reproduce collatz_procgen_20260922_dim.out (exceptional-set dimension lane).
#   QUICK (default, ~10 min, <0.3 GB): validation of both new thread programs against the old counts,
#     forward threads to m=50 (P0=26), backward threads to r=33 (P0=16), alpha table, Theorem-P closure,
#     Z[1/3] census to j=18, analysis.
#   FULL=1 (hours; forward 0.8 GB, backward 0.4 GB; run at most two at once): forward P0=30 to m=64,
#     backward W128/P0=17 to r=41, census to j=25 on the m=61 dump, descent routes for |x|>3/2.
set -euo pipefail
cd "$(dirname "$0")"
P=collatz_procgen_20260922
T=${OUTDIR:-$(mktemp -d)}; mkdir -p "$T"
CC="clang -O3 -march=native"
$CC -o $T/fwd    ${P}_dim_forward.c -lm
$CC -o $T/fwdtt  ${P}_dim_forward_tt.c -lm
$CC -o $T/bwd    ${P}_dim_backward.c -lm
$CC -DW128 -DP0C=17 -DP0C_M0=387420489ULL -o $T/bwd17w ${P}_dim_backward.c -lm
$CC -o $T/alpha  ${P}_dim_alpha.c -lm
$CC -o $T/olddp  ${P}_e_forward_dp_full.c -lm
$CC -o $T/oldbdp ${P}_e_reverse_dp.c -lm

echo "### 1. validation"
echo "--- old float DP (forward) m<=24 vs new exact integer DP (P0=24): counts must agree"
$T/olddp 24 | sed -E 's/^m= *([0-9]+) .*bad=([0-9]+).*/\1 \2/' > $T/old_fwd.txt
$T/fwd 24 24 | grep '^m=' | sed -E 's/^m=([0-9]+) .*bad=([0-9]+).*/\1 \2/' > $T/new_fwd.txt
diff -q $T/old_fwd.txt $T/new_fwd.txt >/dev/null && echo "CHECK: forward counts m=1..24 identical (old float DP vs new integer DP): $(tail -1 $T/new_fwd.txt)" || { echo "MISMATCH forward DP"; exit 1; }
echo "--- combined-lift DFS vs single-lift DFS (P0=20, m<=44): bad lists must be identical"
mkdir -p $T/c $T/s
$T/fwd 20 44 $T/c > $T/c.log; SINGLE_LIFT=1 $T/fwd 20 44 $T/s > $T/s.log
for m in 27 32 35 40 43 44; do cmp -s <(sort -n $T/c/fwd_bad_m$m.txt) <(sort -n $T/s/fwd_bad_m$m.txt) && echo "CHECK m=$m: identical ($(wc -l < $T/c/fwd_bad_m$m.txt) classes)" || { echo "MISMATCH m=$m"; exit 1; }; done
echo "--- transposition-table version (independent search code) must give the same counts"
$T/fwdtt 20 44 22 | grep -E '^m=(27|32|35|40|43|44) ' | awk '{print $1,$4}'
echo "--- known forward counts (choice-ladder note): m=16:124 21:391 24:369 27:255 32:561 35:454 40:1030"
grep -E '^m=(16|21|24|27|32|35|40) ' $T/c.log | awk '{print $1,$4}' | tr '\n' ' '; echo
echo "--- backward: old float DP r<=14 vs new; known r=19:67 24:94 27:134 30:190 31:157, minrep r=30: 7167614823536, r=31: 20230651018121"
$T/oldbdp 14 | grep '^r=' | sed -E 's/^r= *([0-9]+) .*bad=([0-9]+).*/\1 \2/' > $T/old_bwd.txt
$T/bwd 12 31 > $T/bwd_val.log
grep '^r=' $T/bwd_val.log | sed -E 's/^r=([0-9]+) .*bad=([0-9]+).*/\1 \2/' | head -14 > $T/new_bwd.txt
diff -q $T/old_bwd.txt $T/new_bwd.txt >/dev/null && echo "CHECK: backward counts r=1..14 identical (old float DP vs new integer DP)" || { echo "MISMATCH backward DP"; exit 1; }
grep -E '^r=(19|24|27|30|31) ' $T/bwd_val.log | awk '{print $1,$3,$5}' | tr '\n' ' '; echo

echo "### 2. forward threads (exact counts)"
if [ "${FULL:-0}" = 1 ]; then mkdir -p $T/fwd30; $T/fwd 30 64 $T/fwd30 | tee $T/fwd.log; FD=$T/fwd30; FM=61;
else mkdir -p $T/fwd26; $T/fwd 26 50 $T/fwd26 | tee $T/fwd.log; FD=$T/fwd26; FM=50; fi
echo "### 3. backward threads (exact counts)"
if [ "${FULL:-0}" = 1 ]; then mkdir -p $T/bwd17; $T/bwd17w 17 41 0 $T/bwd17 | tee $T/bwd.log; BD=$T/bwd17;
else mkdir -p $T/bwd16; $T/bwd 16 33 0 $T/bwd16 | tee $T/bwd.log; BD=$T/bwd16; fi
echo "### 4. alpha(i) = min x3-moves of an E-path from -1 at its i-th halving (Lemma A / Theorem F)"
$T/alpha 26 56 -1 1
echo "### 5. closure of {-1} under the perturbation lemma (Theorem P): generations, PROVED hostile points"
(cd $T && python3 "$OLDPWD/${P}_dim_generations.py" $T/alpha 26 56 $T/gen_points_56.txt)
echo "### 6. exact prover: examples (hostile certificates; descents incl. slow descenders beyond |x|=3/2)"
python3 ${P}_dim_hostile_prover.py -1 -13/9 -97/81 -355/243 -3211/2187 -87209/59049 -781553/531441 -11/9 -37/27 -7585/6561 -90025/59049 -9905/6561
echo "### 7. census of Z[1/3] points in exceptional classes, classified by the exact prover"
if [ "${FULL:-0}" = 1 ]; then python3 ${P}_dim_z13scan.py $FD/fwd_bad_m$FM.txt $FM 25 30000000;
else python3 ${P}_dim_z13scan.py $FD/fwd_bad_m$FM.txt $FM 18 3000000; fi
echo "### 7b. tree structure of the certified census (parents, relative costs, depth)"
CEN=$(ls $FD/hostile_z13_M${FM}_J*.txt | head -1)
python3 ${P}_dim_census_tree.py $CEN $T/alpha 26 $T/gen_points_56.txt
echo "### 8. growth analysis"
python3 ${P}_dim_analysis.py $T/fwd.log $FD $T/bwd.log $BD $T/gen_points_56.txt
echo "outputs in $T"
