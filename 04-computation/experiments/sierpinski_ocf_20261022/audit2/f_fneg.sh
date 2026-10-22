#!/bin/bash
# Audit2 (F): F(R_n,-1) for odd n <= 21 with hcount.c ("fneg" mode).  Usage: f_fneg.sh HCOUNT_BINARY
# Outputs f_fneg_small.out (n <= 15) and f_fneg.out (n = 17, 19, 21).
HC=$1
for n in 3 5 7 9 11 13 15; do s=""; for j in $(seq 1 $(( (n-1)/2 ))); do s="$s $j"; done; $HC fneg $n $s; done > f_fneg_small.out
for n in 17 19 21; do s=""; for j in $(seq 1 $(( (n-1)/2 ))); do s="$s $j"; done; $HC fneg $n $s; done > f_fneg.out
