#!/usr/bin/env python3
"""Finite AP probe accompanying ``LRCMod13Blocking.lean``.

Corrected scope (codex-2026-07-18-S75): the denominator-300 enumeration gives
lower bounds for the true supremum on eleven displayed APs.  It does not prove
an all-AP or all-family equality classification.  Independently, the witness
``t=1/(2a+11d)`` proves every row with ``a>d`` has margin above ``1/13``.
"""
from math import gcd
from fractions import Fraction as Fr

def Mstar_bounded(V, QMAX=300):
    b = Fr(0)
    for q in range(2, QMAX+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(min((a*v) % q, q - ((a*v) % q)) for v in V)
            if Fr(m, q) > b: b = Fr(m, q)
    return b

print('finite AP probe: max over reduced denominators q<=300 (a lower bound for M)')
for a, d in [(1,1),(2,2),(3,3),(5,5),(2,1),(3,1),(1,2),(1,3),(2,3),(7,1),(1,7)]:
    C = [a + d*k for k in range(12)]
    lower = Mstar_bounded(C)
    elementary = Fr(a, 2*a + 11*d)
    if a > d:
        assert elementary > Fr(1, 13)
        assert lower >= elementary
    suffix = '  (known homogeneous exact row)' if a == d else ''
    print(f'  a={a},d={d}: M_q<=300={lower!s:>7}  '
          f't=1/(2a+11d) gives {elementary!s:>7}  '
          f'a==d(mod13)={((a-d)%13==0)}{suffix}')
print()
print('PROVED here: a>d => M>1/13; all-nonzero AP residues mod13 <=> a==d (mod13).')
print('NOT proved here: the a<d branch, a=d as an integer identity, or the general n=12 tight locus.')
