#!/bin/bash
# Fresh-seed confirmation batches for q_5 (one heavy process at a time). Output: e_p5_confirm.out
D=/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_E
cd "$D"
: > e_p5_confirm.out
python3 e_chain_tail.py 5 400000 50000 40 2201 >> e_p5_confirm.out 2>&1
python3 e_direct_orbits.py 5 400000 400 1024 505 >> e_p5_confirm.out 2>&1
echo DONE >> e_p5_confirm.out
