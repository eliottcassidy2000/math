#!/usr/bin/env python3
"""Audit of THM-4848's Z_13 row: the <3>-cyclotomic class (H = {1,3,9}), k = 2, two of its labellings."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_runobs import circulant, observer2_fast, blocks_data, bareiss_det, detB1_coeffs, sigma, norm_Qzeta, partial_factor
def say(*a): print(" ".join(str(x) for x in a), flush=True)
for S in ({1, 2, 3, 5, 6, 9}, {4, 7, 8, 10, 11, 12}):
    H = [m for m in range(1, 13) if {(m * x) % 13 for x in S} == S]
    M, tup, idx = observer2_fast(13, circulant(S, 13))
    R, Ms = blocks_data(M, tup, idx, 13)
    B0 = [[sum(x) for x in row] for row in Ms]
    dB0 = bareiss_det(B0)
    c, used = detB1_coeffs(Ms, 13)
    inv = all(sigma(c, m, 13) == c for m in H)
    NB = norm_Qzeta(c, 13)
    pf = partial_factor(NB, 10 ** 7)
    say(f"Z_13 S={sorted(S)}: H={H}; det B_0 = {partial_factor(dB0)}; det B_1 fixed by H: {inv}; "
        f"N(det B_1) partial factorisation: {pf}; residues mod 13 of its primes: {{p: p % 13 for p in pf if isinstance(p, int)}}")
