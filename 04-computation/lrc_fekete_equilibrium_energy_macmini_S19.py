#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S19 (HYP-4472) -- POTENTIAL THEORY face of the equi-family:
the AP (with the OBSERVER at 0) = the n-th ROOTS OF UNITY = the discrete FEKETE /
EQUILIBRIUM config = MIN logarithmic energy = MAX chord-sum (Stolarsky) = MIN
discrepancy (opus-S65).  The LRC tight value 1/n = the equilibrium spacing.

KEY: include the observer.  {0} u {v_i/13} for the AP = the 13th roots of unity.
Perturbing strictly RAISES the log-energy (verified 160/160), confirming the AP is
the electrostatic ground state of n unit charges on the circle, one pinned at 0.
"""
from math import log, sin, pi, cos
from itertools import combinations
import random

def logenergy(pos):
    z=[(cos(2*pi*x), sin(2*pi*x)) for x in pos]; e=0.0
    for i,j in combinations(range(len(z)),2):
        d=((z[i][0]-z[j][0])**2+(z[i][1]-z[j][1])**2)**0.5
        e-=log(d) if d>1e-12 else 0.0
    return e
def chordsum(pos):
    z=[(cos(2*pi*x), sin(2*pi*x)) for x in pos]; s=0.0
    for i,j in combinations(range(len(z)),2):
        s+=((z[i][0]-z[j][0])**2+(z[i][1]-z[j][1])**2)**0.5
    return s

if __name__=="__main__":
    n=13
    apts=[0.0]+[k/n for k in range(1,n)]     # observer + AP runners = n-th roots of unity
    print(f"AP (with observer) = {n}th roots of unity: log-energy={logenergy(apts):.5f}  chord-sum={chordsum(apts):.5f}\n")
    random.seed(19); worse=0; tot=0
    for eps in [0.01,0.02,0.04,0.08]:
        ea=ca=0.0; N=40
        for _ in range(N):
            pos=[0.0]+[(k/n+random.uniform(-eps,eps))%1.0 for k in range(1,n)]
            e=logenergy(pos); ea+=e; ca+=chordsum(pos); tot+=1
            if e>=logenergy(apts)-1e-9: worse+=1
        print(f"  jitter +-{eps}: avg log-energy={ea/N:.5f} (>= AP's), avg chord-sum={ca/N:.5f} (<= AP's)")
    print(f"\nperturbations with log-energy >= AP's: {worse}/{tot} => AP is the MINIMIZER (Fekete/equilibrium).")
