#!/bin/sh
# out-of-sample stress tests of HYP-9656 (A) (finiteness of the value set) beyond the packet's runs
S="$1"
./a10_closure.exe "$S/amat_le5.txt" 50000000000 4095 0 27
./a10_closure.exe "$S/amat_le6.txt" 50000000000 255 0 26
./a10_closure.exe "$S/amat_le7.txt" 50000000000 3 0 24
