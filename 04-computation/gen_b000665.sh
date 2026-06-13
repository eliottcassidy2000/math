#!/bin/bash
echo "0 1" > b000665.txt
echo "1 1" >> b000665.txt
echo "2 1" >> b000665.txt
for n in $(seq 3 80); do
    val=$(./a000665_gmp_enum $n 8 | grep "^a(" | sed 's/a([0-9]*) = //')
    echo "$n $val" >> b000665.txt
    echo "n=$n done"
done
