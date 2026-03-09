#!/bin/bash
# A051249: 5-uniform hypergraphs
for n in $(seq 0 40); do
    if [ $n -lt 5 ]; then
        echo "$n 1" >> b051249.txt
    else
        val=$(./k_uniform_gmp_enum 5 $n 8 | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
        echo "$n $val" >> b051249.txt
        echo "k=5 n=$n done"
    fi
done
