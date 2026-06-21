#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exhaustive genuine-wide DOUBLET window certificate (kind-pasteur-S28).
For ALL primitive bounded bases B (size k-2, span<=14, 0 in B), g in {1,2}, M in [15, M*(k)]:
verify p0(B u {M,M+g}) < cap_k. The doublet analogue of THM-563's 12805-base check.
Tail M>M* is the rigorous Mordell-Tornheim bound (HYP-2817, T=12zeta3, N<=15)."""
import importlib.util, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
spec=importlib.util.spec_from_file_location("wrs","04-computation/lrc14_wide_resonance_sup_kpswf7.py")
wrs=importlib.util.module_from_spec(spec); spec.loader.exec_module(wrs)
p0=wrs.p0_fast; CAP=wrs.CAP
def prim(E): return reduce(gcd,[e for e in E if e],0)==1
MSTAR={9:42,10:44,11:38,12:30}
for k in (9,10,11,12):
    cap=float(CAP[k]); nb=k-2; Mstar=MSTAR[k]
    viol=0; worst=0; argw=None; nchk=0; nbases=0
    for sub in itertools.combinations(range(1,15), nb-1):
        B=(0,)+sub
        if not prim(B): continue
        nbases+=1
        for g in (1,2):
            for M in range(15,Mstar+1):
                E=tuple(sorted(set(list(B)+[M,M+g])))
                if len(E)!=k or not prim(E): continue
                pv=float(p0(E)); nchk+=1
                if pv>=cap: viol+=1; print(f"  VIOLATION k={k} B={B} g={g} M={M} p0={pv}")
                if pv>worst: worst=pv; argw=(B,g,M)
    print(f"k={k}: {nbases} bases, {nchk} configs, window[15,{Mstar}], cap={cap:.4f} | violations={viol} | worst p0={worst:.5f} margin={cap-worst:.5f} at {argw}",flush=True)
print("DONE: exhaustive doublet window certificate.")
