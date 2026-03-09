#!/bin/bash
cd /Users/e/Documents/GitHub/math/04-computation
for n in $(seq 36 60); do
    val=$(./k_uniform_fast_enum 6 $n 8 2>/dev/null | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
    echo "$n $val" >> b309860.txt
    echo "n=$n done ($(date +%H:%M:%S))"
done
echo "Done extending b309860.txt"
