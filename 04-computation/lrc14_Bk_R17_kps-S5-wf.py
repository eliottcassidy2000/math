#!/usr/bin/env python3
"""Claim (4): R_1/7 = rho*/(meas(G_P)*mu_1/7) over consecutive E, k=8..13. EXACT."""
from fractions import Fraction as F
from itertools import combinations
import sys
src=open('04-computation/lrc14_Bk_verify_gpintersectionun_kps-S5-wf.py').read().split('if __name__')[0]
ns={}; exec(src,ns)
mu_exact=ns['mu_exact']; safe_GP=ns['safe_GP']; meas=ns['meas']
good_E=ns['good_E']; intersect_intervals=ns['intersect_intervals']
th17=F(1,7)
Rmin=None; Rarg=None
for k in range(8,14):
    sz=13-k; E=list(range(k)); muc=mu_exact(E,th17); GE=good_E(E,th17)
    kmin=None; karg=None
    for Pset in combinations(range(1,14), sz):
        GP=safe_GP(list(Pset)); gp=meas(GP)
        rho=meas(intersect_intervals(GP,GE))
        denom=gp*muc
        if denom>0:
            R=rho/denom
            if kmin is None or R<kmin: kmin=R; karg=Pset
            if Rmin is None or R<Rmin: Rmin=R; Rarg=(k,Pset)
    print(f"k={k}: min R={kmin}={float(kmin):.6f} at P={karg}", flush=True)
print(f"OVERALL min R_1/7={Rmin}={float(Rmin):.6f} at {Rarg}", flush=True)
print(f"claimed >= 67053/84241 = {float(F(67053,84241)):.6f}", flush=True)
print(f"Rmin >= 67053/84241 ? {Rmin>=F(67053,84241)}", flush=True)
print(f"Rmin == 67053/84241 (exact) ? {Rmin==F(67053,84241)}", flush=True)
print("DONE", flush=True)
