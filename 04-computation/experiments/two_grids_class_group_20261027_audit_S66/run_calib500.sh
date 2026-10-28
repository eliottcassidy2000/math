#!/bin/sh
# growth of counts from B = 200 to B = 500 in the residue class of 6 (6 mod 16), then B = 800 for 6 and 54
cd "$(dirname "$0")"
for n in 22 38 86 134 166; do
  ./a3_search2.exe deep $n 500 0 | grep -E "^deep n = .*processed" >> a3_calib500.out
done
echo DONECALIB500 >> a3_calib500.out
./run_deep800.sh
