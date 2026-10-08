#!/bin/bash
# audit_F: anisotropy test (C3) and battery checks (E), sequential
cd /Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/reframes_20261007/audit_F
OUT0=partC4_extra.out
: > $OUT0
echo "=== 3x+5 with offset 5 (conjugate to 3x+1 with offset 1): should merge" >> $OUT0
python3 bigint_pairs.py 3x+5 2000 4096 230 --offset 5 >> $OUT0
echo "=== product-like Z_6 map m_i = alpha(i mod 2) beta(i mod 3), alpha=(1,5), beta=(1,1,7): rank 2" >> $OUT0
python3 bigint_pairs.py '6:1,0;5,1;7,4;5,3;1,2;35,5' 1000 16384 231 >> $OUT0
echo DONE >> $OUT0
OUT=partC3_anisotropy.out
: > $OUT
python3 anisotropy.py Z3_125   600 16384 301 6 >> $OUT
python3 anisotropy.py Z5_12311 600 16384 302 6 >> $OUT
python3 anisotropy.py Z5_12371 600 16384 303 6 >> $OUT
python3 anisotropy.py '4:1,0;3,1;5,2;7,3' 600 16384 304 6 >> $OUT
python3 anisotropy.py '3:1,0;2,1;11,2'    400 16384 305 6 >> $OUT
echo DONE >> $OUT
OUT2=partE_battery.out
: > $OUT2
python3 battery_check.py mersenne 12800 >> $OUT2
python3 battery_check.py epsgen 1000000 >> $OUT2
python3 battery_check.py landing 4499 >> $OUT2
echo DONE >> $OUT2
python3 anchor_check.py > partE_anchors.out
echo DONE >> partE_anchors.out
