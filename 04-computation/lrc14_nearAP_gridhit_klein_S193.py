#!/usr/bin/env python3
"""
klein-2026-07-08-S193: STRESS-TEST the large-spread pigeonhole on ADVERSARIAL
clusters -- the near-dilated-AP family e = d*{0..9} u {p} (the extremal LOW-rho*
shape), which last session's RANDOM gridhit did not cover.

Concern: for a block, #arcs/spread = (k+1)/(k-1) ~ 1.2 > 1, and for near-AP
rho* (= mu_{1/7} * meas(G_P)) may be only ~0.5, so the pigeonhole
#{good j} >= rho*Vmax - #arcs could go NEGATIVE (c > rho*).  Does a good ruler
period STILL exist (the LOWER bound being vacuous does not mean #good=0)?

Report, at worst-case Vmax = spread + 14 (Vmin=14): #arcs, rho*, the pigeonhole
lower bound rho*Vmax-#arcs, and the ACTUAL #good_j.  k=11, P={1,2}.
"""
import numpy as np
from math import gcd

INV7=1.0/7.0

def maxgap_mask_and_arcs(e, N):
    e=np.asarray(e,float); x=(np.arange(N)+0.5)/N
    Ph=np.mod(np.outer(x,e),1.0); Ph.sort(axis=1)
    g=np.diff(Ph,axis=1); g=np.concatenate([g,(1.0-Ph[:,-1]+Ph[:,0])[:,None]],axis=1)
    good=(g.max(axis=1)>INV7+1e-12)
    gi=good.astype(int); edges=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((edges==1).sum());
    if gi.all(): nc=1
    return good, nc

def gridhit(e, Vmax, P):
    j=np.arange(Vmax); x=(j+0.5)/Vmax
    inGP=np.ones(Vmax,bool)
    for p in P:
        dd=np.abs(((p*x+0.5)%1.0)-0.5); inGP&=(dd>=1.0/14-1e-12)
    e=np.asarray(e,float); Ph=np.mod(np.outer(x,e),1.0); Ph.sort(axis=1)
    g=np.diff(Ph,axis=1); g=np.concatenate([g,(1.0-Ph[:,-1]+Ph[:,0])[:,None]],axis=1)
    good=inGP&(g.max(axis=1)>INV7+1e-12)
    return int(good.sum()), float(inGP.mean())

P=[1,2]
print("Adversarial near-dilated-AP  e = d*{0..9} u {p},  k=11, P={1,2}, Vmax=spread+14")
print(f"{'d':>4} {'p':>4} {'spread':>7} {'Vmax':>7} {'#arcs':>6} {'rho*':>6} {'measGP':>7} "
      f"{'rho*Vmax-#arcs':>14} {'ACTUAL #good':>12}")
worst=10**9
for d in (5,10,20,40,80,150,300):
    for p in (1,2,7,11):
        if gcd(p,d)!=1: continue
        e=list(d*np.arange(10))+[p]
        e=sorted(set(e))
        if len(e)<11: continue
        spread=max(e)-min(e)
        Vmax=spread+14
        good_mask,nc=maxgap_mask_and_arcs(e, 200000)
        rho_cluster=good_mask.mean()
        ng,measGP=gridhit(e,Vmax,P)
        rho_star=rho_cluster*measGP    # approx (independence heuristic; ng is exact)
        pig=rho_star*Vmax-nc
        worst=min(worst,ng)
        flag="  <-- FAIL" if ng==0 else ("  <-- pigeonhole vacuous" if pig<1 else "")
        print(f"{d:>4} {p:>4} {spread:>7} {Vmax:>7} {nc:>6} {rho_cluster:>6.3f} {measGP:>7.3f} "
              f"{pig:>14.1f} {ng:>12}{flag}")

print(f"\nWORST actual #good over all near-AP tested = {worst}")
print("If worst>=1 despite pigeonhole vacuous, the naive #arcs<rho*Vmax bound is")
print("TOO WEAK for near-AP -- the good periods exist for a finer (resonance) reason.")
print("If worst=0, the large-spread half genuinely FAILS on near-AP at Vmax=spread+14.")
