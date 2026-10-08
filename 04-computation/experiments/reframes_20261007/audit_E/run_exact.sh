#!/bin/bash
# Sequential rigorous shortest-merge / exact-absorption runs (one heavy process at a time). Output: e_shortest_exact.out
D=/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_E
cd "$D"
: > e_shortest_exact.out
python3 e_shortest_exact.py reset >> e_shortest_exact.out 2>&1
# confirm the author's first-merge times (L = claimed value captures every absorption at times <= L)
for pl in "3 3" "5 11" "7 4" "9 11" "11 18" "13 23" "15 5" "17 13" "21 8" "23 18" "31 6"; do
  set -- $pl
  python3 e_shortest_exact.py $1 $2 3000000 >> e_shortest_exact.out 2>&1
done
# p where the author found nothing to depth 26 (pruning there was |f| > 1e4, not rigorous)
for p in 19 25 27 29; do
  python3 e_shortest_exact.py $p 26 3000000 >> e_shortest_exact.out 2>&1
done
# exact P(absorbed by L) for larger L (rigorous lower bounds on q_p)
python3 e_shortest_exact.py 5 22 3000000 >> e_shortest_exact.out 2>&1
python3 e_shortest_exact.py 7 22 3000000 >> e_shortest_exact.out 2>&1
python3 e_shortest_exact.py 9 24 3000000 >> e_shortest_exact.out 2>&1
python3 e_shortest_exact.py 11 28 3000000 >> e_shortest_exact.out 2>&1
python3 e_shortest_exact.py 13 30 3000000 >> e_shortest_exact.out 2>&1
python3 e_shortest_exact.py 19 30 3000000 >> e_shortest_exact.out 2>&1
echo DONE >> e_shortest_exact.out
