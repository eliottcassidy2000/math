#!/bin/sh
# targeted order-12 sampling for the K1 junction window just above the order-10 window (H filtered to [1500, 6500])
i=0
for q in 0.91 0.915 0.92 0.925 0.93 0.935 0.94 0.945; do
  i=$((i+1)); ./tiltsample 12 $q 2500000 $((3000+i)) 1500 6500 > j12_$i.txt &
done
wait
