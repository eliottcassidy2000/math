#!/bin/bash
# Compute a(151)..a(200) to extend b000568.txt
# Uses a000568_gmp_enum_v4 with 8 threads

DIR="$(cd "$(dirname "$0")" && pwd)"
BINARY="$DIR/a000568_gmp_enum_v4"
BFILE="$DIR/b000568.txt"
OUTDIR="$DIR"

echo "Computing a(151)..a(200) with 8 threads"

for n in $(seq 151 200); do
    echo "=== n=$n ==="
    result=$("$BINARY" $n 8 2>&1)
    # Extract the value (first line, after "a(N) = ")
    value=$(echo "$result" | head -1 | sed "s/a($n) = //")
    timing=$(echo "$result" | tail -1)
    echo "$n $value" >> "$BFILE"
    echo "$timing"
done

echo "Done. b-file extended to n=200."
