#!/bin/bash
cd /Users/e/Documents/GitHub/math/04-computation
for n in $(seq 41 70); do
    val=$(./a051249_gmp_enum $n 8 2>/dev/null | grep "^a(" | sed 's/a([0-9]*) = //')
    echo "$n $val" >> b051249.txt
    echo "n=$n done ($(date +%H:%M:%S))"
done
echo "Done extending b051249.txt"
