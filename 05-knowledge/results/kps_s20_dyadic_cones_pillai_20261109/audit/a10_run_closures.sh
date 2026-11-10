#!/bin/sh
# audit closures (HYP-9656 P1, P2, P3, P5/P6, prefix table) with a10_closure.exe; matrices from a09 (amat_le*.txt)
S="$1"
./a10_closure.exe "$S/amat_le7.txt" 10000000000 1 0 22 "$S/p1_vals_exact.txt"
./a10_closure.exe "$S/amat_le6.txt" 100000000 63 0 23
./a10_closure.exe "$S/amat_le6.txt" 1000000000 63 0 23
./a10_closure.exe "$S/amat_le6.txt" 1000000000 63 1 23
for r in 3 7 15 31 63 127 255 511 1023; do ./a10_closure.exe "$S/amat_le5.txt" 1000000000 $r 0 25; done
./a10_closure.exe "$S/amat_le5.txt" 1000000000 1023 1 25
./a10_closure.exe "$S/amat_le5.txt" 4000000000 2047 0 26 "$S/p5_vals_exact.txt"
./a10_closure.exe "$S/amat_le5.txt" 4000000000 2047 1 26 "$S/p5_vals_auth.txt"
