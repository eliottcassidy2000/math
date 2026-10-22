#!/usr/bin/env python3
"""Brute-force O_5 and O_6 (ordered mutually orthogonal root 5-/6-frames) for several irreducible systems,
compared with the EGF recursion used for the true l_5; also R1/R2 O_5 via components."""
import sys; sys.path.insert(0, '.')
import numpy as np
from math import factorial
from fractions import Fraction as Fr
from rootsys import roots, gram, count_tuples
src = open('b_paper2_formulas.py').read()
ns = {}
exec(src.split("# brute-force check of O_m")[0], ns)   # defines Z, logcoeffs, series_mul, M, R1, R2 ...
Z = ns['Z']
for t, k in [('A', 9), ('D', 6), ('D', 8), ('E', 6), ('E', 7), ('E', 8)]:
    R, sc = roots(t, k); IP = gram(R, sc)
    for m in (5, 6):
        bf = R.shape[0]*count_tuples(IP, 2*np.eye(m, dtype=np.int64), first_candidates=[0])
        rec = Z(t, k)[m]*factorial(m)
        print(f"{t}{k}: O_{m} brute force = {bf}, recursion = {rec}, agree: {bf == rec}")
