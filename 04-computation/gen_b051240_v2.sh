#!/bin/bash
# Generate extended b-file for A051240 using dedicated closed-form enumerator
outfile="b051240_v2.txt"
> $outfile
for n in $(seq 0 80); do
    if [ $n -lt 4 ]; then
        echo "$n 1" >> $outfile
    else
        val=$(./a051240_gmp_enum $n 8 | grep "^a(" | sed 's/a([0-9]*) = //')
        echo "$n $val" >> $outfile
        echo "n=$n done"
    fi
done
echo "Done: $outfile"
