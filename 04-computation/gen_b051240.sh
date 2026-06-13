#!/bin/bash
echo "0 1" > b051240.txt
echo "1 1" >> b051240.txt
echo "2 1" >> b051240.txt
echo "3 1" >> b051240.txt
for n in $(seq 4 50); do
    val=$(./k_uniform_gmp_enum 4 $n 8 | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
    echo "$n $val" >> b051240.txt
    time_info=$(./k_uniform_gmp_enum 4 $n 8 2>/dev/null | grep "Time:")
    echo "n=$n done ($time_info)"
done
