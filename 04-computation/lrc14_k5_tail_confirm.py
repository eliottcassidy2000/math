#!/usr/bin/env python3
"""Final confirm: k=5 gcd=1 tail shapes (spread 17-18) + scaling-tail limit, line-buffered."""
import sys, itertools
from fractions import Fraction as F
from math import gcd; from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
import importlib.util
spec=importlib.util.spec_from_file_location('m','/Users/e/Documents/GitHub/math/04-computation/lrc14_k5_exact_closure_macmini.py')
mm=importlib.util.module_from_spec(spec); sys.argv=['x']; spec.loader.exec_module(mm)
Ps=[list(p) for p in itertools.combinations(range(1,14),8)]
def minP_rho(E):
    g=mm.good_set_exact(list(E)); mn=F(2); argp=None
    for P in Ps:
        r=mm.meas(mm.intersect(g,mm.safe_set(P)))
        if r<mn: mn=r; argp=tuple(P)
    return mn,argp
floor=F(95,2548)
print("FLOOR (bounded-spread, exhaustive s<=16) = 95/2548 =",float(floor),flush=True)
print("\nk=5 gcd=1 tail shapes, spread 17-18 (the genuinely-new, non-scaling shapes):",flush=True)
gmin=(F(2),None,None); below=0; nshapes=0
for s in (17,18):
    smin=(F(2),None)
    for T in itertools.combinations(range(1,s),3):
        E=tuple([0]+list(T)+[s])
        if reduce(gcd,E)!=1: continue
        nshapes+=1
        r,P=minP_rho(E)
        if r<smin[0]: smin=(r,E)
        if r<floor: below+=1
        if r<gmin[0]: gmin=(r,E,P)
    print(f"  s={s}: min rho* over gcd=1 shapes = {float(smin[0]):.6f} at {smin[1]}",flush=True)
print(f"  gcd=1 shapes checked: {nshapes}; # below 95/2548: {below}",flush=True)
print(f"  global min (s17-18 gcd1): {gmin[0]} = {float(gmin[0]):.6f} at E={gmin[1]} P={gmin[2]}",flush=True)
