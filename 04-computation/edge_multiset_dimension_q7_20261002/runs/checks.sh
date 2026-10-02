#!/bin/bash
# Cheap independent checks (single core, nice).  Outputs in logs/.
cd "$(dirname "$0")/.."
N="nice -n 10"
$N python3 py/pycheck.py selftest > logs/pycheck_selftest.out 2>&1
$N python3 py/xcheck.py 400 7 > logs/xcheck.out 2>&1
# orbit representative files: every listed set is lex-min in its orbit (explicit group), no duplicates
: > logs/repcheck.out
# negative control: {0,1,3} (not lex-min), {0,1,2} twice (duplicate), {1,2,4} (no 0): expect noncanonical=2 duplicates=1
$N ./src/repcheck 5 3 data/neg_n5_a3.bin > logs/repcheck_negative.out
$N ./src/repcheck 5 3 data/neg_n5_a3.bin full >> logs/repcheck_negative.out
for a in $(seq 0 32); do $N ./src/repcheck 5 $a data/q5/reps_n5_a$a.bin full >> logs/repcheck.out; done
for a in $(seq 0 8); do $N ./src/repcheck 6 $a data/q6/reps_n6_a$a.bin full >> logs/repcheck.out; done
for a in 9 10; do $N ./src/repcheck 6 $a data/q6/reps_n6_a$a.bin >> logs/repcheck.out; done
# normal-form reduction, empirically (random sets -> constructed normal form is in the search domain and in the same orbit)
: > logs/nf_test.out
$N python3 py/nf_test.py 6 14 400 1 >> logs/nf_test.out 2>&1
$N python3 py/nf_test.py 6 15 400 2 >> logs/nf_test.out 2>&1
for k in 11 12 13; do $N python3 py/nf_test.py 7 $k 150 $k >> logs/nf_test.out 2>&1; done
echo CHECKS-DONE >> logs/nf_test.out
