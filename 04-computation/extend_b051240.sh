#!/bin/bash
# Extend b051240.txt from n=51 to n=80
cd /Users/e/Documents/GitHub/math/04-computation
for n in $(seq 51 80); do
    val=$(./a051240_gmp_enum $n 8 2>/dev/null | grep "^a(" | sed 's/a([0-9]*) = //')
    echo "$n $val" >> b051240.txt
    echo "n=$n done ($(date +%H:%M:%S))"
done
echo "Done extending b051240.txt"
