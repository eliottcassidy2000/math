#!/bin/bash
# Generate b-files for all k-uniform hypergraph sequences
# A309861 (k=7), A309862 (k=8), A309863 (k=9), A309864 (k=10)
# Also extends A309865 triangle

cd /Users/e/Documents/GitHub/math/04-computation

for k in 7 8 9 10; do
    case $k in
        7) seqid="309861" ;;
        8) seqid="309862" ;;
        9) seqid="309863" ;;
        10) seqid="309864" ;;
    esac

    outfile="b${seqid}.txt"
    echo "=== Generating $outfile (k=$k uniform) ==="
    > $outfile

    max_n=$((35 + k))  # higher k → smaller values → can go further
    if [ $k -ge 9 ]; then max_n=30; fi

    for n in $(seq 0 $max_n); do
        if [ $n -lt $k ]; then
            echo "$n 1" >> $outfile
        else
            val=$(./k_uniform_gmp_enum $k $n 8 | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
            echo "$n $val" >> $outfile
        fi
    done
    echo "Done: $outfile (n=0..$max_n)"
done

# Also generate A309865 triangle rows 15..25
echo "=== Generating A309865 triangle extension ==="
> b309865_ext.txt
for n in $(seq 0 25); do
    for k in $(seq 0 $n); do
        if [ $k -eq 0 ] || [ $k -eq $n ]; then
            val=2
        elif [ $n -lt $k ]; then
            val=1
        else
            val=$(./k_uniform_gmp_enum $k $n 8 | grep "^a(" | sed 's/a([0-9]*,[0-9]*) = //')
        fi
        echo "T($n,$k) = $val" >> b309865_ext.txt
    done
    echo "Row n=$n done"
done
echo "Done: b309865_ext.txt"
