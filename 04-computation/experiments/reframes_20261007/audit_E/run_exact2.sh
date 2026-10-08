#!/bin/bash
# Deeper exact enumerations (rigorous lower bounds on q_11, q_13). Output: e_shortest_exact2.out
D=/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_E
cd "$D"
: > e_shortest_exact2.out
python3 e_shortest_exact.py 11 32 6000000 >> e_shortest_exact2.out 2>&1
python3 e_shortest_exact.py 13 33 6000000 >> e_shortest_exact2.out 2>&1
echo DONE >> e_shortest_exact2.out
