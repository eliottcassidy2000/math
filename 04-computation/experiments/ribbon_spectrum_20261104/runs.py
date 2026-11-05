#!/usr/bin/env python3
"""runs.py -- closed-form probes: H of constant-letter ribbon words (runs) and spectral data of the letters."""
import sys
from fractions import Fraction
from ribbon import transfer_system, letters, count_with_matrices, word_str
import numpy as np

def run_seq(w, a, nmax, sysw):
    return [count_with_matrices(w, [a] * n, sysw) for n in range(1, nmax + 1)]

for w in (1, 2, 3):
    sysw = transfer_system(w)
    order, states, mats, final = sysw
    d = len(order)
    print(f"=== w={w}: {d} states")
    for a in letters(w):
        if not a:
            continue
        seq = run_seq(w, a, 14, sysw)
        # dense matrix and spectral radius
        M = np.zeros((d, d))
        for (j, i), m in mats[a].items():
            M[j, i] = m
        ev = np.linalg.eigvals(M)
        rho = max(abs(ev))
        # ratio guess
        rat = [seq[i + 1] / seq[i] for i in range(len(seq) - 1)]
        print(f"letter {sorted(a)}: H(run_n) n=1..14: {seq}")
        print(f"    spectral radius {rho:.5f}; last ratios {[round(r,4) for r in rat[-3:]]}")
        # test 2^(n-2)+1
        print("    2^(n-2)+1 match:", [seq[n-1] == 2**(n-2) + 1 for n in range(3, 15)])
