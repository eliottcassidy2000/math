#!/bin/bash
# Sequential direct-orbit runs (one heavy process at a time). Output: e_direct_orbits.out
D=/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_E
cd "$D"
: > e_direct_orbits.out
python3 e_direct_orbits.py 5 600000 1000 1536 501 >> e_direct_orbits.out 2>&1
python3 e_direct_orbits.py 7 600000 400 1024 701 >> e_direct_orbits.out 2>&1
python3 e_direct_orbits.py 9 600000 400 1024 901 >> e_direct_orbits.out 2>&1
python3 e_direct_orbits.py 11 4000000 200 512 1101 >> e_direct_orbits.out 2>&1
python3 e_direct_orbits.py 13 4000000 200 512 1301 >> e_direct_orbits.out 2>&1
python3 e_direct_orbits.py 5 100000 400 3000 503 >> e_direct_orbits.out 2>&1
echo DONE >> e_direct_orbits.out
