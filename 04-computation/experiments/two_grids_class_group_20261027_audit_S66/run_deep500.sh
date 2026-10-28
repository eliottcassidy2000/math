#!/bin/sh
# B = 500 for the seven persistent misses, then calibration counts at B = 300 for their residue classes mod 16
cd "$(dirname "$0")"
for n in 6 54 60 124 140 78 75; do
  ./a3_search2.exe deep $n 500 5 | grep -E "^deep n = .*processed|hit:" >> a3_deep500.out
done
echo DONE500 >> a3_deep500.out
