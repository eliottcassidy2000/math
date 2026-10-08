#!/bin/bash
# Sequential long-horizon chain runs (own chain), after run_direct.sh and run_exact.sh. Output: e_chain_tail.out
D=/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_E
cd "$D"
: > e_chain_tail.out
python3 e_chain_tail.py 5 400000 50000 40 2101 >> e_chain_tail.out 2>&1
python3 e_chain_tail.py 7 400000 20000 40 2102 >> e_chain_tail.out 2>&1
python3 e_chain_tail.py 9 400000 20000 40 2103 >> e_chain_tail.out 2>&1
echo DONE >> e_chain_tail.out
