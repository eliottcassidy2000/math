#!/bin/sh
# k*(n): least k with f(n+k,k,k) = g(n), the catalytic limit taken from lattice_W2500.txt (lattice_debt.c);
# f from merge_debt.c with K = 4n + 8.  If f(n+K,K,K) > g(n) the threshold lies beyond K and is reported as "> K".
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; do
  K=$((4 * n + 8)); g=$(grep "^g($n) = " lattice_W2500.txt | sed "s/^g([0-9]*) = \([0-9]*\).*/\1/")
  ./merge_debt $n $K | awk -v n=$n -v g=$g -v K=$K '{ split($0, A, ")="); split(A[2], B, " "); k = NR - 1; if (B[1] == g && ks == "") ks = k; if (k == 0) f0 = B[1]; last = B[1] }
    END { if (ks == "") ks = "> " K " (f at K = " last ")"; print "n=" n " f(n,0,0)=" f0 " g(n)=" g " k*=" ks }'
done
