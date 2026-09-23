#!/usr/bin/env python3
"""Orchestrator's independent fairness check of the HYP-9128 super-blocks (N = 16, 64).
Uses the lane's construction only to obtain the integers f_{i,r}; then, independently of Lemma S and the class-polynomial
algebra, rebuilds Phi_block(p) = sum_{L=N}^{4N-1} p^L q sum_z e_{L,z} p^z q^{R_L - z}  (q = 1-p) as an exact integer
polynomial in p and checks: (i) Phi_block(p) == Phi_block(1-p) identically (fairness of the block);
(ii) Phi_block == w^N S_N(w), w = p(1-p); (iii) realizability |e| <= C(R,z), e == C(R,z) mod 2; (iv) T(L) = L+1+R_L <= ceil(1.59 L)."""
import sys, importlib.util
from math import comb
from fractions import Fraction
spec = importlib.util.spec_from_file_location("fin", "04-computation/experiments/amm12592_procgen_20260923_hyp9128_finite.py")
fin = importlib.util.module_from_spec(spec); spec.loader.exec_module(fin)
def padd(a, b):
    n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: out[i+j] += x*y
    return out
def ppow(a, k):
    r = [1]
    for _ in range(k): r = pmul(r, a)
    return r
def compose_one_minus(a):   # a(1-p)
    r = [0]; one_minus = [1, -1]; pw = [1]
    for c in a:
        r = padd(r, [c*v for v in pw]); pw = pmul(pw, one_minus)
    return r
def trim(a):
    while len(a) > 1 and a[-1] == 0: a = a[:-1]
    return a
for N in (16, 64):
    out, data = fin.margins_certificate(N)
    ok, err, sites, f = fin.lemma_R_round(N, data)
    t, R, a = fin.profile(N)
    n = 3*N - 1
    q = [1, -1]; p = [0, 1]
    Phi = [0]
    box_ok = True; dl_ok = True
    for i in range(n+1):
        L = N + i
        if not (L + 1 + R[i] <= fin.ceil_frac(Fraction(159,100)*L) and L + 1 + R[i] <= 4*N): dl_ok = False
        E = [0]
        for r in range(R[i]+1):
            e = (-1)**(i+r) * f[i][r]
            if abs(e) > comb(R[i], r) or (e - comb(R[i], r)) % 2: box_ok = False
            if e:
                E = padd(E, [e*v for v in pmul(ppow(p, r), ppow(q, R[i]-r))])
        term = pmul(pmul(ppow(p, L), q), E)
        Phi = padd(Phi, term)
    Phi = trim(Phi)
    sym = trim(compose_one_minus(Phi)) == Phi
    S = fin.state_S(N)
    w = [0, 1, -1]
    target = [0]
    for k, s in enumerate(S):
        if s: target = padd(target, [s*v for v in ppow(w, N + k)])
    eq = trim(target) == Phi
    print(f"N={N}: levels {N}..{4*N-1}; realizable (box+parity) {box_ok}; deadlines <= ceil(1.59L) and <= 4N {dl_ok}; "
          f"Phi(p)==Phi(1-p) {sym}; Phi == w^N S_N(w) {eq}; deg Phi = {len(Phi)-1}")
