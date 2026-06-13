#!/bin/bash
# A309860: 6-uniform hypergraphs
for n in $(seq 0 35); do
    if [ $n -lt 6 ]; then
        echo "$n 1" >> b309860.txt
    else
        val=$(./k_uniform_gmp_enum 6 $n 8 | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
        echo "$n $val" >> b309860.txt
        echo "k=6 n=$n done"
    fi
done
