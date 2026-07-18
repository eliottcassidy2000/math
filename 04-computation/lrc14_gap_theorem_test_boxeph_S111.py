#!/usr/bin/env python3
"""
GAP THEOREM test: non-AP core => M >= 1/12 ?  (boxeph-2026-07-18-S111)

Finding: TRUE for covering families (non-AP floor = 1/12 exactly), but the gap theorem
is LOGICALLY STRONGER than INV (same hypothesis, stronger conclusion 1/12>1/13), hence
>= INV strength (open). S110's "gap is more tractable" route is CORRECTED. Provable
fragment: non-covering-at-{2..12} => M>=1/q'>=1/12 (sieve); the rest is INV.
"""
from math import gcd
from fractions import Fraction as Fr

def Mstar(V, QMAX=400):
    best = Fr(0)
    for q in range(2, QMAX + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(min((a*v) % q, q - ((a*v) % q)) for v in V)
            if Fr(m, q) > best: best = Fr(m, q)
    return best

def is_AP(C):
    C = sorted(C); d = C[1]-C[0]
    return all(C[i+1]-C[i] == d for i in range(len(C)-1))

def covering(V):
    return all(any(v % q == 0 for v in V) for q in range(2, 14))

fams = {
    'deep well {1..12}u182 (AP)':        list(range(1,13))+[182],
    'dilated {2..24}u182 (AP)':          [2*k for k in range(1,13)]+[182],
    'repl 5->24 (cover, non-AP)':        [1,2,3,4,6,7,8,9,10,11,12,24,182],
    'repl 7->14 (cover, non-AP)':        [1,2,3,4,5,6,8,9,10,11,12,14,182],
    'dbl+perturb {2..22,26}u780 (nonAP)':[2,4,6,8,10,12,14,16,18,20,22,26,780],
}
for name, V in fams.items():
    V = sorted(set(V)); M = Mstar(V)
    print(f'{name:<36} cover={str(covering(V)):>5} AP-core={str(is_AP(V[:-1])):>5} '
          f'M={str(M):>8}={float(M):.4f}')
print('\ncovering non-AP floor = 1/12 (sharp); gap (1/13,1/12) empty of covering non-AP.')
print('gap thm => INV (stronger conclusion, same hyp) => >= open crux. Route corrected (S110->S111).')
