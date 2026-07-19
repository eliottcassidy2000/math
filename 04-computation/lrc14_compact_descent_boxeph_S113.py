#!/usr/bin/env python3
"""
The compact case rho<13 => M>=1/13: descent is too weak  (boxeph-2026-07-18-S113)

Finding: compact rho<13 covering => M>=1/13 is the SOLE RESIDUAL = LRC(14) (S86), SHARP
(boundary {2*{1..12},13} at M=1/13). Descent (THM-1010, M>=rho*M(core)/(rho+1)) loses a
factor ~2 at rho~1 (compact) -- proves only 5/15. Descent belongs to single-killer (large rho).
Provable sub-case: dilated-AP-core compact => M>=1/13 (THM-1013 dilated sieve).
"""
from math import gcd
from fractions import Fraction as Fr
import random

def Mstar(V, QMAX=250):
    best = Fr(0)
    for q in range(2, QMAX+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(min((a*v) % q, q - ((a*v) % q)) for v in V)
            if Fr(m, q) > best: best = Fr(m, q)
    return best

def covering(V): return all(any(v % q == 0 for v in V) for q in range(2, 14))

# the SHARP boundary compact family: M=1/13 exactly
V = [2,4,6,8,10,12,13,14,16,18,20,22,24]
core = V[:-1]
print(f'boundary compact family {V}:')
print(f'  covering={covering(V)}  rho={Fr(V[-1],core[-1])}={float(Fr(V[-1],core[-1])):.2f}')
print(f'  M(core)={Mstar(core)}  M(V)={Mstar(V)}  (=1/13 EXACTLY => the bound is SHARP)')
print()
print('descent LB = rho*M(core)/(rho+1) ~ M(core)/2 at rho~1 << actual M ~ M(core) => descent too weak.')
print('compact rho<13 => M>=1/13 is equivalent to LRC(14) (open); dilated-AP-core sub-case = THM-1013.')
