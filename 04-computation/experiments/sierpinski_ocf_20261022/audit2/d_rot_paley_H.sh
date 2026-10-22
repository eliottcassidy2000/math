#!/bin/bash
# Audit2 (D): exact H(R_n) (odd n <= 27) and H(P_19), H(P_23) by inclusion-exclusion with rotation orbits (hcount.c).
H=$1
for n in 3 5 7 9 11 13 15 17 19 21 23 25 27; do
  s=""; for j in $(seq 1 $(( (n-1)/2 ))); do s="$s $j"; done
  t0=$(date +%s); v=$($H iecirc $n $s); t1=$(date +%s)
  echo "R_$n H=$v ($((t1-t0))s)"
done
for q in 19 23; do
  s=$(python3 -c "q=$q; print(' '.join(str(x) for x in sorted({(t*t)%q for t in range(1,q)})))")
  t0=$(date +%s); v=$($H iecirc $q $s); t1=$(date +%s)
  echo "P_$q H=$v ($((t1-t0))s)"
done
