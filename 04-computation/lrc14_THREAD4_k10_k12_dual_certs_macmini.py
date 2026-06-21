#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_k10_k12_dual_certs_macmini.py  (THREAD 4, mac-mini, 2026-06-21)

EXACT-RATIONAL Delsarte dual certificates for the BINDING ROWS k=10 and k=12,
ready for Lean formalization (mirror gK8/gK9/gK11 in LRCFactorialAtom.lean).

The Delsarte leg requires:  max_E L_y(E) <= cap_k  with a Krawtchouk-nonneg dual
g(t) = sum_{r<=t} y_r C(t,r), g(0)>=1, g(t)>=0 (t=1..6), so that for every shape E,
    measS7(E) = q_0(E) <= sum_t g(t) q_t(E) = L_y(E).

ALREADY FORMALIZED (LRCFactorialAtom.lean):
  k=8        : yK8 *10 = (10,-10,10,-9,6,0,0),  g=(10,0,0,1,0,0,10)
  k=9,10     : yK9 *18 = (18,-13,8,-3,0,0,0),   g=(18,5,0,0,2,3,0)
  k=11,12,13 : yK11*6  = (6,-3,1,0,0,0,0),      g=(6,3,1,0,0,1,3)

The task wants STANDALONE k=10 and k=12 certificates (their own gK10, gK12 lemmas),
with the optimal (largest-margin) dual found by LP, snapped to exact rationals and
VERIFIED exactly, then the integer (scaled) covector for Lean.

For each k we:
  (1) run the exact LP over span<=BANK shapes to find the MIN-max dual at degrees 2,3,4;
  (2) snap the LP optimum to exact rationals (the LP basis is rational);
  (3) VERIFY exactly: g(0)>=scale, g(t)>=0, max_E L_y <= cap_k, with exact margin;
  (4) print the integer (scaled) yK, gK covectors + readout, ready for Lean.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass
import numpy as np
from scipy.optimize import linprog

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}
BANK = {10:14, 12:15}     # bounded-span finite-check region

def atom_q(E):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(bps); q=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*abs(e)*xm)%7 for e in Enz)
        occ=len([s for s in hit if s!=0]); missed=6-occ
        q[missed]+=L
    return q

def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1

def collect_shapes(k,span):
    return [(0,)+rest for rest in itertools.combinations(range(1,span+1),k-1)
            if primitive((0,)+rest)]

def Sr_of_q(q): return [sum(comb(t,r)*q[t] for t in range(7)) for r in range(7)]

def snap(x, dens=(1,2,3,4,5,6,7,8,9,10,12,14,15,18,20,21,24,28,30,36,40,42,45,56,60,63,70,72,84,90,98,105,126,140,180,210,252,280,360,420,504,840,2940,5880)):
    """Snap a float to the simplest rational with small denominator near it."""
    best=None; berr=1e9
    for d in dens:
        n=round(x*d)
        f=F(n,d); e=abs(float(f)-x)
        if e<berr-1e-15 or (abs(e-berr)<1e-15 and (best is None or d<best.denominator)):
            berr=e; best=f
    return best

def run_k(k):
    cap=CAP[k]; span=BANK[k]
    shapes=collect_shapes(k,span)
    Q=[atom_q(E) for E in shapes]
    Srs=[Sr_of_q(q) for q in Q]
    print("\n"+"#"*86)
    print(f"# BINDING ROW k={k}:  cap_{k} = {cap} = {float(cap):.8f}   ({len(shapes)} span<={span} shapes)")
    print("#"*86)
    # find consec argmax q0 (the measS7 maximizer floor)
    q0max=max(q[0] for q in Q)
    print(f"  max_E measS7 (q0) = {q0max} = {float(q0max):.8f}  (the unbeatable g>=0 floor)")
    results={}
    for R in (2,3,4):
        # vars y_1..y_R (y_0=1) + M. minimize M.
        nv=R+1
        c=np.zeros(nv); c[R]=1.0
        A=[]; b=[]
        # g(t)>=0, t=1..6:  -sum y_r C(t,r) <= 1
        for t in range(1,7):
            row=np.zeros(nv)
            for r in range(1,R+1): row[r-1]=-comb(t,r)
            A.append(row); b.append(1.0)
        # 1 + sum y_r S_r(E) <= M
        for Sr in Srs:
            row=np.zeros(nv)
            for r in range(1,R+1): row[r-1]=float(Sr[r])
            row[R]=-1.0
            A.append(row); b.append(-1.0)
        bounds=[(None,None)]*R+[(0,None)]
        res=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method="highs")
        if not res.success:
            print(f"  deg {R}: LP failed"); continue
        # snap y to exact
        y=[F(1)]+[snap(v) for v in res.x[:R]]
        # build exact g and verify exactly
        g=[sum(y[r]*comb(t,r) for r in range(R+1)) for t in range(7)]
        feas = g[0]>=1 and all(g[t]>=0 for t in range(1,7))
        # exact max_E L_y
        best=F(-10); bestE=None
        for E,q in zip(shapes,Q):
            ly=sum(g[t]*q[t] for t in range(7))
            if ly>best: best,bestE=ly,E
        margin=cap-best
        results[R]=(y,g,feas,best,bestE,margin)
        ok = feas and margin>0
        print(f"\n  deg {R}: y = {[str(v) for v in y]}")
        print(f"          g = {[str(v) for v in g]}   feasible(g(0)>=1,g>=0)={feas}")
        print(f"          max_E L_y = {best} = {float(best):.8f}   argmax={list(bestE)}")
        print(f"          margin cap_{k} - max L_y = {margin} = {float(margin):+.8f}  [{'CLEARS' if ok else 'FAILS'}]")
    # pick best feasible by margin
    feas_results={R:v for R,v in results.items() if v[2] and v[5]>0}
    if feas_results:
        bestR=max(feas_results, key=lambda R: feas_results[R][5])
        y,g,feas,best,bestE,margin=feas_results[bestR]
        # integer scale
        scale=reduce(lambda a,b: a*b//gcd(a,b), [v.denominator for v in g], 1)
        gint=[v*scale for v in g]; yint=[v*scale for v in y]
        gint=[int(v) for v in gint]; yint=[int(v) for v in yint]
        print(f"\n  >>> CHOSEN k={k} dual: degree {bestR}, margin {float(margin):+.8f} (largest)")
        print(f"      EXACT y = {[str(v) for v in y]}")
        print(f"      EXACT g = {[str(v) for v in g]}")
        print(f"      INTEGER (x{scale}) yK{k} = {yint}")
        print(f"      INTEGER (x{scale}) gK{k} = {gint}")
        # readout in q form: L_y * scale = sum_t gint[t] q_t
        terms=[f"{gint[t]}*q{t}" for t in range(7) if gint[t]!=0]
        print(f"      Lean readout: L_yK{k}(q) = {' + '.join(terms)}")
        print(f"      Lean bound:   {gint[0]}*q0 <= L_yK{k}(q)   (since gint[t]>=0)")
        print(f"      finite check: max_E L_yK{k} = {best*scale} <= {scale}*cap_{k} = {scale*cap}")
        return (k,bestR,y,g,yint,gint,scale,best,margin)
    return None

def main():
    print("EXACT k=10 and k=12 Delsarte dual certificates for Lean (THREAD 4)")
    out=[]
    for k in (10,12):
        r=run_k(k)
        if r: out.append(r)
    print("\n"+"="*86)
    print("LEAN-READY SUMMARY")
    print("="*86)
    for (k,R,y,g,yint,gint,scale,best,margin) in out:
        print(f"  k={k}: deg {R}, scale x{scale}")
        print(f"     yK{k} = {yint}   gK{k} = {gint}   margin = {margin} (+{float(margin):.6f})")

if __name__=="__main__":
    main()
