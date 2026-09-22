#!/usr/bin/env bash
# Reproduce collatz_procgen_20260922_choice_ladder.out (about 3-5 minutes; DFS deepening dominates).
set -euo pipefail
cd "$(dirname "$0")"
P=collatz_procgen_20260922
T=$(mktemp -d)
for f in e_forward_dp e_forward_dp_full e_forward_dfs e_reverse_dp e_reverse_dump deterministic_counts loops_through_one; do
  clang -O3 -o $T/$f ${P}_$f.c -lm
done
echo "### 1. deterministic counts (Collatz T mod 2^m; greedy G mod 3^(J+1))"; $T/deterministic_counts
echo "### 2. E-graph forward (Q1) full choice, DP mod 2^m"; (cd $T && ./e_forward_dp_full 26 1)
echo "### 3. partial choice E_S (excursion allowed iff residue in S)"
for spec in "1" "1 0" "2 2" "2 0" "3 4" "3 0" "3 2" "3 6" "4 8" "4 0"; do echo "--- J,S = $spec"; $T/e_forward_dp 22 $spec | tail -3; done
echo "### 4. DFS deepening of Q1 bad threads (independent code path), levels 27..32"
(cd $T && ./e_forward_dfs 26 32 q1bad_dump.txt && sort -n q1deep_m27.txt > dfs27.txt)
(cd $T && ./e_forward_dp_full 27 1 >/dev/null && awk '{print $1}' q1bad_dump.txt | sort -n > dp27.txt && if diff -q dp27.txt dfs27.txt >/dev/null; then echo "CHECK: DP(m=27) == DFS(m=27) bad sets identical ($(wc -l < dp27.txt) classes)"; else echo "CHECK FAILED"; exit 1; fi)
echo "### 5. E-graph reverse (Q2) DP mod 3^(r+1)"; $T/e_reverse_dp 14
echo "### 6. independent Python DFS check of Q2 bad lists"; for r in 6 8 10; do python3 ${P}_e_reverse_dfs_check.py $r; done
echo "### 7. reverse E-loops through 1 (minimal K per length s)"; $T/loops_through_one 40
echo "### 8. minimal escape price from the 2-adic point -1"; python3 ${P}_minus_one_escape_price.py
rm -rf $T
