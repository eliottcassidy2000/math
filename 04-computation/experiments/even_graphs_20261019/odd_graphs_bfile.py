#!/usr/bin/env python3
"""Odd graphs (graphs with an automorphism reversing an odd number of edges of an orientation) up to isomorphism:
odd(n) = A000088(n) - A000568(n) = sum over cycle types with an even part of 2^{c_2(lambda)}/z_lambda (Royle et al. Thm 4.4).
Also labelled-even and other small sequences from the census for the record."""
from math import gcd, factorial
from fractions import Fraction
from collections import Counter
def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for r in partitions(n - k, k): yield (k,) + r
def zlam(p):
    c = Counter(p); z = 1
    for L, m in c.items(): z *= L ** m * factorial(m)
    return z
def c2(lam):
    return sum(L // 2 for L in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
N = 60
odd = []; even = []
for n in range(1, N + 1):
    o = Fraction(0); e = Fraction(0)
    for lam in partitions(n):
        w = Fraction(2 ** c2(lam), zlam(lam))
        if any(L % 2 == 0 for L in lam): o += w
        else: e += w
    assert o.denominator == 1 and e.denominator == 1
    odd.append(int(o)); even.append(int(e))
print("even graphs (= A000568):", even[:12])
print("odd graphs (A000088 - A000568), n = 1..12:", odd[:12])
with open('odd_graphs_bfile.txt', 'w') as f:
    for n in range(1, N + 1): f.write(f"{n} {odd[n-1]}\n")
print("odd graphs b-file written to n =", N)
