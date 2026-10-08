#!/bin/bash
# audit_F part C: counterexample hunt (sequential; one heavy process at a time)
cd /Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_F
OUT=partC_maps.out
: > $OUT
echo "=== C1 obstruction maps (proof + simulation)" >> $OUT
python3 obstruction_check.py 400 4096 >> $OUT
echo "=== C2 further contracting maps, direct big-integer orbits" >> $OUT
# rank 1 with tied multipliers
python3 bigint_pairs.py '5:1,0;2,3;2,1;2,4;2,2'            1000 16384 201 >> $OUT
# rank 2
python3 bigint_pairs.py '3:1,0;2,1;7,1'                     1000 16384 202 >> $OUT
python3 bigint_pairs.py '3:1,0;4,-1;5,-1'                   1000 16384 203 >> $OUT
python3 bigint_pairs.py '3:1,0;2,1;11,2'                    1000 16384 204 >> $OUT
python3 bigint_pairs.py '5:1,0;3,2;9,2;2,4;1,1'             1000 16384 205 >> $OUT
python3 bigint_pairs.py '4:1,0;3,1;5,2;1,1'                 1000 16384 206 >> $OUT
python3 bigint_pairs.py '7:1,0;2,5;3,1;4,2;6,4;1,2;1,1'     1000 16384 207 >> $OUT
python3 bigint_pairs.py '5:1,0;6,4;11,3;1,2;1,1'            1000 16384 208 >> $OUT
# rank 2, other sheets of the HYP's Z3 (1,2,5)
python3 bigint_pairs.py '3:1,0;2,4;5,2'                     1000 16384 209 >> $OUT
python3 bigint_pairs.py '3:1,0;2,1;5,-1'                    1000 16384 210 >> $OUT
python3 bigint_pairs.py '3:1,0;2,-2;5,5'                    1000 16384 211 >> $OUT
# rank 3
python3 bigint_pairs.py '4:1,0;3,1;5,2;7,3'                 1000 16384 212 >> $OUT
python3 bigint_pairs.py '5:1,0;8,-3;3,4;7,4;12,2'           1000 16384 213 >> $OUT
python3 bigint_pairs.py '5:1,0;2,3;3,4;7,4;6,1'             1000 16384 214 >> $OUT
python3 bigint_pairs.py '7:1,0;2,5;3,1;5,6;1,3;1,2;1,1'     1000 16384 215 >> $OUT
python3 bigint_pairs.py '5:1,0;6,4;11,3;16,2;1,1'           1000 16384 216 >> $OUT
# rank 4
python3 bigint_pairs.py '5:1,0;2,3;3,4;7,4;11,1'            1000 16384 217 >> $OUT
# long runs: weak rank-2 and weak rank-3
python3 bigint_pairs.py '3:1,0;2,1;11,2'                    600 65536 218 >> $OUT
python3 bigint_pairs.py '5:1,0;8,-3;3,4;7,4;12,2'           600 65536 219 >> $OUT
python3 bigint_pairs.py '4:1,0;3,1;5,2;7,3'                 600 65536 220 >> $OUT
echo DONE >> $OUT
