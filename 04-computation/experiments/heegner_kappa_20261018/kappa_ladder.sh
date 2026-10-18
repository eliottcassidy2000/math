#!/bin/zsh
# For each q, find the least k (from KSTART) at which tabucol_gq finds a proper colouring of G_q = Cay(F_(q^2), mu_(q+1)).
# Usage: ./kappa_ladder.sh KSTART TIMEOUT ITERS q1 q2 ...
KSTART=$1; T=$2; IT=$3; shift 3
for q in "$@"; do
  k=$KSTART
  while [ $k -le 16 ]; do
    out=$(timeout $T ./tabucol_gq $q $k $IT 13 2>&1)
    echo "$out"
    if echo "$out" | grep -q "final_conflicts(verified)=0"; then break; fi
    k=$((k+1))
  done
done
