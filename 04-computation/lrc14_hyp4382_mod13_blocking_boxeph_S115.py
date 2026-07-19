#!/usr/bin/env python3
"""
Sharpening HYP-4382: mod-13 pair-blocking is PROVED but necessary-not-sufficient (boxeph-S115)

PROVED: M(C)=1/13, 13 not| c_i  =>  {+-c_i mod 13} = {1..12} (mod-13 pair-blocking).
  (t=b/13 gives M >= min_i ||c_i b/13||, so min_i |c_i b mod 13| <= 1 => c_i b == 0,+-1 mod 13.)
NECESSARY not SUFFICIENT: non-AP complete-mod-13 families reach only M=1/12, not 1/13.
So HYP-4382's AP rigidity is the offset-vanishing (mod-val coord), orthogonal to mod-13; = the crux.
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

def is_AP(C):
    C = sorted(C); d = C[1]-C[0]
    return all(C[i+1]-C[i] == d for i in range(len(C)-1))

def complete_mod13(C): return set(c % 13 for c in C) == set(range(1, 13))

print('PROVED necessary condition + necessary-not-sufficient demonstration:')
for name, C in [('AP {1..12}', list(range(1,13))),
                ('dilated 5*{1..12}', [5*i for i in range(1,13)]),
                ('non-AP {1..11,25} (complete mod13)', list(range(1,12))+[25]),
                ('non-AP {2..12,14} (complete mod13)', list(range(2,13))+[14]),
                ('AP-diff-13 {1,14,..144} (NOT blocker)', [1+13*k for k in range(12)])]:
    M = Mstar(C)
    print(f'  {name:<40} complete13={complete_mod13(C)!s:>5} AP={is_AP(C)!s:>5} M={M!s:>7} (=1/13? {M==Fr(1,13)})')
print()
print('Only the (dilated) AP is tight; complete-mod-13 non-APs give 1/12 (beaten at another q).')
print('=> mod-13 blocking PROVED but necessary-not-sufficient; AP rigidity = offset-vanishing = the crux.')
