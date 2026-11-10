#!/bin/sh
# tail probes of HYP-9656 (A) to HMAX = 10^13 (three-word-key variant a10b_closure)
S="$1"
./a10b_closure.exe "$S/amat_le6.txt" 10000000000000 255 0 25 "$S/tail_le6_r255.txt"
./a10b_closure.exe "$S/amat_le6.txt" 10000000000000 63 0 23 "$S/tail_le6_r63.txt"
./a10b_closure.exe "$S/amat_le5.txt" 10000000000000 1023 0 26 "$S/tail_le5_r1023.txt"
./a10b_closure.exe "$S/amat_le7.txt" 10000000000000 1 0 22 "$S/tail_le7_r1.txt"
