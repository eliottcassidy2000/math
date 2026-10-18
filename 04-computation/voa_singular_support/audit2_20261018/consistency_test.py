#!/usr/bin/env python3
"""Module-consistency test: X(Y v) - (-1)^{|X||Y|} Y(X v) == [X,Y} v for random modes and random PBW states.
Also checks that the computed Gram matrices are (super)symmetric.
usage: python consistency_test.py n2 c  |  python consistency_test.py w3 c"""
import sys, random
from fractions import Fraction
from aud_algebras import N2, W3

random.seed(12345)
kind, c = sys.argv[1], Fraction(sys.argv[2])
va = N2(c) if kind == 'n2' else W3(c)
P = va.P
ntests, fails = 0, 0
for trial in range(600):
    W2 = random.randint(0, 14)
    B = va.basis(W2)
    if not B:
        continue
    v = random.choice(B)
    t1 = random.randrange(len(va.gens)); t2 = random.randrange(len(va.gens))
    def rmode(t):
        h2 = va.gens[t][1]
        par_int = (h2 % 2 == 0)
        while True:
            n = random.randint(-8, 8)
            n2 = 2 * n if par_int else 2 * n + 1
            return (t, n2)
    X, Y = rmode(t1), rmode(t2)
    lhs = {}
    for m, cc in va.apply_vec(X, va.apply(Y, v)).items():
        lhs[m] = (lhs.get(m, 0) + cc) % P
    sgn = P - 1 if (va.par(X) and va.par(Y)) else 1
    for m, cc in va.apply_vec(Y, va.apply(X, v)).items():
        lhs[m] = (lhs.get(m, 0) - sgn * cc) % P
    lhs = {m: cc for m, cc in lhs.items() if cc}
    rhs = va.bracket_apply(X, Y, v)
    ntests += 1
    if lhs != rhs:
        fails += 1
        if fails < 5:
            print('FAIL', X, Y, v, lhs, rhs)
print(f'{kind} c={c}: {ntests} consistency tests, {fails} failures')
# Gram symmetry
G = va.gram_all(16 if kind == 'n2' else 14)
asym = 0
for key, (Bb, Gm) in G.items():
    if Gm is None:
        continue
    n = len(Bb)
    for i in range(n):
        for j in range(n):
            if int(Gm[i, j]) != int(Gm[j, i]):
                asym += 1
print('Gram blocks:', len(G), 'asymmetric entries:', asym)
