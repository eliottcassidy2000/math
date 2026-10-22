"""Auditor: t0 for the 13 doubly regular 23-tournaments from mac-mini's forward-step coefficients f_k
(sierpinski_ocf_20261022/drt23_forest_exact.json), via sum_L (n-L)! N_L t^L = sum_k f_k (1+t)^k (1-t)^(n-1-k)
(THM-4757 section 5).  The conversion is first checked against my own forest numbers on random small tournaments
(F(T,x) by brute force over permutations).  Usage: python drt23_check.py"""
import os, sys, json, itertools, math, random
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn
HERE = os.path.dirname(os.path.abspath(__file__))
t = sp.symbols('t')

def c_from_f(f, n):
    poly = sp.expand(sum(int(fk) * (1 + t) ** k * (1 - t) ** (n - 1 - k) for k, fk in enumerate(f)))
    return [int(poly.coeff(t, L)) for L in range(n)]

rng = random.Random(3)
for n in (5, 6, 7):
    code = ''.join(rng.choice('01') for _ in range(n * (n - 1) // 2))
    A = tourn.from_code(code)
    f = [0] * n
    for s in itertools.permutations(range(n)):
        f[sum(A[a, b] for a, b in zip(s, s[1:]))] += 1
    (H, c, _), = tourn.forest_numbers([code], n)
    print(f'n={n}: conversion check {c_from_f(f, n) == c} (H = {H} = f_(n-1) = {f[n - 1]})')
D = json.load(open(os.path.join(HERE, '..', '..', 'sierpinski_ocf_20261022', 'drt23_forest_exact.json')))
n = 23
res = []
for d in D:
    c = c_from_f(d['f'], n)
    N = tourn.N_from_c(c, n)
    assert all(N[L] == 0 for L in range(1, n, 2))
    t0, w, kind = tourn.smallest_zero(N, n)
    res.append(t0)
    print(f"{d['label']}: H = {d['H']}, t0 = {t0:.6f}, w = {w.real:.6f}{'+' if w.imag >= 0 else '-'}{abs(w.imag):.1e}i ({kind})")
print(f'min {min(res):.5f} median {np.median(res):.5f} max {max(res):.5f}')
