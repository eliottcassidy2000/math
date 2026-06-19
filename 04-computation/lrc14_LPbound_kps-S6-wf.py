#!/usr/bin/env python3
"""
lrc14_LPbound_kps-S6-wf.py
RIGOROUS uniform upper bound on net(E)=meas{x: 8 pts frac(ex) form a 1/7-net}
via a test-function / averaging LP.

FRAMEWORK.  Let W:[0,1)->R_{>=0} (circle).  Suppose for EVERY configuration of
k points p_1..p_k on the circle with covering radius <= 1/14 (max gap <= 1/7),
we have  sum_i W(p_i) >= T.   (T = guaranteed lower bound = min over 1/7-nets.)
Then for x in N(E):  sum_{e in E} W(frac(ex)) >= T.
Integrate over x in [0,1):
   T * meas(N(E)) <= int_{N(E)} sum_e W(frac ex) dx <= int_0^1 sum_e W(frac ex) dx
                   = W(0) + (k-1) * int_0^1 W   [since for e!=0, int W(frac ex)dx = int W;
                                                  for e=0, W(frac 0)=W(0) constant].
So   meas(N(E)) <= ( W(0) + (k-1) * Iavg ) / T,   Iavg = int_0^1 W.
This holds for ALL E with 0 in E, |E|=k.  UNIFORM. No spread assumption.

We need:  ( W(0) + 7*Iavg ) / T  <  cap_8 = 2243/5880 = 0.381463.

We OPTIMIZE W over piecewise-constant nonneg functions on a fine grid, computing
T = min over 1/7-nets of sum W(p_i) exactly-ish (T is a covering-LP itself).  Since
the maximizer is consecutive (equispaced-ish), we expect T achieved near equispaced
configs.  We do a numeric optimization to find a GOOD W, then we will VERIFY the
resulting bound with exact rationals and an exact/rigorous lower bound on T.

This script: numeric exploration to find whether such W EXISTS (feasibility of the
LP bound < cap_8).  If yes, we then build an exact certificate.
"""
import sys
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

# We need T = min over k-point 1/7-nets of sum_i W(p_i).  A 1/7-net: k points,
# all cyclic gaps <= 1/7.  Minimizing a sum of a function over such configs.
# We approximate by sampling many random 1/7-nets AND by structured ones, taking min.
# For a *valid* upper bound on net(E) we need a *valid lower bound* on T; for the
# EXPLORATION we use the min over a rich sample (gives an optimistic T; we re-verify).

K=8
RAD=1.0/7.0   # max gap
NG=70         # grid resolution for W (W constant on each of NG cells)

def gen_nets(num):
    """Random k-point configs with all gaps<=1/7, returned as point arrays."""
    nets=[]
    rng=np.random.default_rng(0)
    # gaps: k nonneg summing to 1, each <= 1/7. Sample via Dirichlet + rejection/repair.
    for _ in range(num):
        # sample gaps in simplex with cap 1/7
        while True:
            g=rng.dirichlet(np.ones(K))
            if g.max()<=RAD+1e-12: break
        # random rotation and the config always includes a point; place first at 0..
        start=rng.random()
        pts=(start+np.concatenate([[0],np.cumsum(g)[:-1]]))%1.0
        nets.append(np.sort(pts))
    # add structured: equispaced rotated (the suspected minimizer of net but here for T)
    for s in np.linspace(0,1,200,endpoint=False):
        pts=(s+np.arange(K)/K)%1.0
        nets.append(np.sort(pts))
    # near-equispaced with one merged
    for s in np.linspace(0,1,100,endpoint=False):
        g=np.full(K,1.0/K);
        pts=(s+np.concatenate([[0],np.cumsum(g)[:-1]]))%1.0
        nets.append(np.sort(pts))
    return nets

def cell(u):
    return int(u*NG)%NG

def Tval(Wvec, nets):
    T=1e9
    for pts in nets:
        s=sum(Wvec[cell(p)] for p in pts)
        if s<T: T=s
    return T

def bound(Wvec, nets):
    Iavg=Wvec.mean()         # int_0^1 W = mean of cell values
    W0=Wvec[cell(0.0)]       # W(0)
    T=Tval(Wvec,nets)
    if T<=0: return 1e9,T,Iavg,W0
    return (W0+(K-1)*Iavg)/T, T,Iavg,W0

if __name__=="__main__":
    nets=gen_nets(4000)
    print(f"K={K}, max-gap={RAD:.4f}, NG={NG} cells, {len(nets)} sample 1/7-nets")
    print(f"cap_8 = 0.381463  (target: find W with bound < cap_8)")

    rng=np.random.default_rng(1)
    # Start: W = indicator-ish that is LARGE away from where equispaced pts cluster?
    # Heuristic init: W(u) peaked, then optimize by coordinate descent to minimize bound.
    best=None; bestW=None
    for init in range(6):
        if init==0:
            W=np.ones(NG)
        elif init==1:
            # peak near 0 (since 0 always a point, W(0) costs us; avoid) -> peak near 1/14?
            W=np.array([1.0 if (0.05<((i+0.5)/NG)%1<0.45 or 0.55<((i+0.5)/NG)%1<0.95) else 0.2 for i in range(NG)])
        else:
            W=rng.random(NG)+0.1
        cur,T,Ia,W0=bound(W,nets)
        # coordinate descent
        for sweep in range(40):
            improved=False
            order=rng.permutation(NG)
            for i in order:
                old=W[i]
                for cand in (old*0.5, old*1.5, old+0.1, max(0.0,old-0.1), 0.0, old*0.8, old*1.2):
                    if cand<0: continue
                    W[i]=cand
                    nb,_,_,_=bound(W,nets)
                    if nb<cur-1e-9:
                        cur=nb; improved=True
                    else:
                        W[i]=old
            if not improved: break
        cur,T,Ia,W0=bound(W,nets)
        if best is None or cur<best:
            best=cur; bestW=W.copy()
        print(f"  init {init}: bound ~ {cur:.4f}  (T={T:.4f}, Iavg={Ia:.4f}, W0={W0:.4f})")
    print(f"\nBEST numeric bound ~ {best:.4f}  vs cap_8=0.381463  "
          f"=> {'FEASIBLE (< cap_8) -- LP route can work!' if best<0.381463 else 'NOT below cap_8 with this setup'}")
    if bestW is not None:
        print("  best W (cell values, rounded):")
        print("   ", np.round(bestW,2).tolist())
