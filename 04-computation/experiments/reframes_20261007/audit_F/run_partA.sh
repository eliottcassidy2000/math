#!/bin/bash
# audit_F part A: independent big-integer simulations, run sequentially (one heavy process at a time)
cd /Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_F
OUT=partA_bigint.out
: > $OUT
python3 bigint_pairs.py x+1      2000  1024  110 >> $OUT
python3 bigint_pairs.py 3x+1     4000 65536  101 >> $OUT
python3 bigint_pairs.py 5x+1     5000  4096  104 >> $OUT
python3 bigint_pairs.py Z3_124   2000 16384  106 >> $OUT
python3 bigint_pairs.py Z3_125   2000 65536  102 >> $OUT
python3 bigint_pairs.py Z3_1416  1500  4096  105 >> $OUT
python3 bigint_pairs.py Z3_157   1000  4096  109 >> $OUT
python3 bigint_pairs.py Z5_12311 1500 16384  107 >> $OUT
python3 bigint_pairs.py Z5_12471 1500 16384  108 >> $OUT
python3 bigint_pairs.py Z5_12371 2000 65536  103 >> $OUT
echo DONE >> $OUT
