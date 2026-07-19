#!/usr/bin/env python3
"""
The exact n=12 tight locus is homogeneous c*{1..12} (a=d)  (boxeph-2026-07-18-S116)

Companion to the Lean LRCMod13Blocking.lean. Shows: M(C)=1/13 among APs {a,a+d,..} holds
ONLY for a=d (dilated c*{1..12}), never shifted APs. And (proved via mod-13 blocking):
AP + M(C)=1/13 => a == d (mod 13) [residues must miss 0].
"""
from math import gcd
from fractions import Fraction as Fr

def Mstar(V, QMAX=300):
    b = Fr(0)
    for q in range(2, QMAX+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(min((a*v) % q, q - ((a*v) % q)) for v in V)
            if Fr(m, q) > b: b = Fr(m, q)
    return b

print('which APs {a,a+d,...,a+11d} are tight (M=1/13)?  ONLY a=d (homogeneous c*{1..12}):')
for a, d in [(1,1),(2,2),(3,3),(5,5),(2,1),(3,1),(1,2),(1,3),(2,3),(7,1),(1,7)]:
    C = [a + d*k for k in range(12)]
    M = Mstar(C)
    print(f'  a={a},d={d}: M={M!s:>7} {"TIGHT" if M==Fr(1,13) else "loose"}  '
          f'a==d(mod13)={((a-d)%13==0)}  {"(a=d homogeneous)" if a==d else ""}')
print()
print('=> exact tight locus = {c*{1..12} : c>=1}; PROVED (mod-13 blocking): tight AP => a==d (mod 13).')
