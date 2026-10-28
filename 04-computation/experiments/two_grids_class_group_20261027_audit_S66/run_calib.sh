#!/bin/sh
# Calibration: representation counts (deep-mode multiset-split solutions) at B = 200 for every n <= 200 in the residue
# classes mod 16 of the persistent misses (6, 54: 6 mod 16; 60, 124, 140: 12 mod 16; 75: 11 mod 16; 78: 14 mod 16).
cd "$(dirname "$0")"
for n in 6 22 38 54 70 86 102 118 134 150 166 182 198  12 28 44 60 76 92 108 124 140 156 172 188  11 27 43 59 75 91 107 123 139 155 171 187  14 30 46 62 78 94 110 126 142 158 174 190; do
  ./a3_search2.exe deep $n 200 0 | grep -E "^deep n = .*processed" >> a3_calib200.out
done
echo DONECALIB >> a3_calib200.out
