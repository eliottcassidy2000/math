"""
mac-mini-2026-07-08-S58 -- DIRECT rho_K for structured (2-block / AP) clusters: does a good
PERIOD actually exist, or is the soft discrepancy bound just loose?

rho_K(E,Vmax) = (1/Vmax) #{ j=0..Vmax-1 : maxgap{ frac(e_i j / Vmax) } > 1/7 }.
For structured E the good set has #arcs ~ beta*spread with beta>rho*, so the discrepancy bound
#{good} >= Vmax*rho* - #arcs is VACUOUS. Compute rho_K EXACTLY on the integer grid to see if
good periods exist anyway (equidistribution of the cyclic orbit {j*(e_i)/Vmax}).

If rho_K > 0 for all Vmax >= spread: good periods exist (soft bound was just loose).
If rho_K = 0 for some Vmax: a real arithmetic obstruction (needs quantitative equidistribution).
"""
import numpy as np
from math import gcd
from functools import reduce
from fractions import Fraction as F

def rho_K_exact(E, Vmax):
    """exact count of good periods j: maxgap{frac(e_i j/Vmax)}>1/7, using integer residues."""
    E = sorted(set(e % Vmax for e in E))
    ngood = 0
    for j in range(Vmax):
        ph = sorted({(e*j) % Vmax for e in E})   # phases * Vmax (integers mod Vmax)
        m = len(ph)
        # circular gaps (in units of 1/Vmax); maxgap>1/7 <=> maxgap_int > Vmax/7
        mg = 0
        for i in range(m):
            g = (ph[(i+1) % m] - ph[i]) % Vmax
            if i == m-1: g = (ph[0] + Vmax - ph[m-1])
            mg = max(mg, g)
        # also handle single distinct phase
        if m == 1: mg = Vmax
        if mg * 7 > Vmax:     # maxgap > 1/7
            ngood += 1
    return ngood, ngood/Vmax
def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("DIRECT rho_K on the integer grid -- do good periods exist for structured clusters?\n")
# 2-block co-offset clusters (the worst beta) at spread s, swept over Vmax >= s
def twoblock(k, s):
    a = k//2; b = k - a
    return sorted(set(list(range(a)) + [s-b+1+i for i in range(b)]))
def apshape(k, s):
    d = s//(k-1); return sorted(set([d*i for i in range(k-1)] + [s]))

for k in (11, 13):
    for name, mk in [('2-block', twoblock), ('AP', apshape)]:
        for s in [30, 60]:
            E = mk(k, s)
            if len(E) != k or not primitive(E):
                # tweak to primitive
                E = sorted(set(E) | {1}); E = E[:k]
            mn = 1.0; argmin = None; nzero = 0
            for Vmax in range(s, s+140):
                if Vmax < max(E): continue
                ng, rk = rho_K_exact(E, Vmax)
                if rk < mn: mn = rk; argmin = Vmax
                if ng == 0: nzero += 1
            print(f"k={k} {name:8s} s={s:3d} E={E}")
            print(f"     over Vmax in [{max(s,max(E))}, {s+139}]: min rho_K = {mn:.4f} at Vmax={argmin}; "
                  f"#(Vmax with rho_K=0) = {nzero}  {'*** ZERO EXISTS ***' if nzero else '(all have a good period)'}")
