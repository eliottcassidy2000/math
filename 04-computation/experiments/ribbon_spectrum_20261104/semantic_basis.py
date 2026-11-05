#!/usr/bin/env python3
"""semantic_basis.py -- print the reduced integer matrices A_a of span 2 and find a 'semantic' Q-basis of
the observable space made of short-suffix functionals u -> H(u v), |v| <= 2 (greedy by length)."""
import itertools, numpy as np
from fractions import Fraction
from intrep import build, vecmat, in_module_coeffs
from ribbon import letters, word_str
w = 2
basis, A, y0, hcoef, order, states, mats, final = build(w, verbose=False)
r = len(basis); d = len(order); L = letters(w)
names = {frozenset(): '-', frozenset([1]): '1', frozenset([2]): '2', frozenset([1, 2]): '12'}
for a in L:
    print(f"A_{{{names[a]}}} (rows = output coordinates):")
    for k in range(r):
        print("   ", " ".join(f"{x:3d}" for x in A[a][k]))
phi = [0] * d
for j in final: phi[j] = 1
# functionals u -> H(u v): phi M_v, as rows in the reduced coordinates
funcs = []
for n in range(0, 3):
    for v in itertools.product(L, repeat=n):
        f = phi
        for a in v: f = vecmat(f, mats[a], d)   # careful: H(u v) = phi M_{v_n} ... M_{v_1} x(u): apply letters in order
        coeffs = in_module_coeffs(basis, f)
        funcs.append((n, v, coeffs))
# note: vecmat(f, M) computes f.M, so after the loop f = phi M_{v_1} M_{v_2} ...; the correct order is
# phi M_{v_n} ... M_{v_1}, so rebuild with reversed order
funcs = []
for n in range(0, 3):
    for v in itertools.product(L, repeat=n):
        f = phi
        for a in reversed(v): f = vecmat(f, mats[a], d)
        funcs.append((n, v, in_module_coeffs(basis, f)))
chosen = []
mat = []
for n, v, c in funcs:
    cand = mat + [c]
    if np.linalg.matrix_rank(np.array(cand, dtype=float)) > len(mat):
        mat = cand; chosen.append((n, v))
print("a Q-basis of the 14 observables by short suffixes v (u -> H(u v)):", [''.join(names[a] + '.' for a in v).rstrip('.') or 'eps' for n, v in chosen])
print("number of suffix functionals of length <= 2:", len(funcs), "; rank:", np.linalg.matrix_rank(np.array([c for n, v, c in funcs], dtype=float)))
