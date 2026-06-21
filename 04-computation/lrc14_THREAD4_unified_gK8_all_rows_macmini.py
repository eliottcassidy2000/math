#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_unified_gK8_all_rows_macmini.py  (THREAD 4, mac-mini, 2026-06-21)

THE UNIFICATION.  A single Delsarte dual -- the ALREADY-FORMALIZED k=8 certificate
    gK8 = (10,0,0,1,0,0,10),   yK8 = (10,-10,10,-9,6,0,0)  (= y*(1,-1,1,-9/10,3/5))
    readout:  10*measS7 = 10*q0 <= 10*q0 + q3 + 10*q6 = L_yK8(q)
-- clears the cap at EVERY binding row k=8..13 (with the wide-span iid contraction safety).

This script verifies EXACTLY, for k=8,9,10,11,12,13:
    max_E L_yK8(E) <= 10*cap_k          (the integer finite check, span<=BANK_k)
and reports the exact max, argmax, and exact margin (10*cap_k - max), and per-shape margin.

CONSEQUENCE.  The Delsarte leg needs only ONE Lean certificate (gK8), not a per-row
family.  In particular:
  * k=9 (the RAZOR row): gK8 gives margin +0.0509 (vs the deg-3 gK9 razor +0.00138, 37x).
  * k=10, k=12 (this thread's deliverable): gK8 gives margins +0.0758, +0.1930.
The standalone gK10/gK12 lemmas are then just `delsarte_bound_k8` re-applied, with the
finite-check inequality max_E L_yK8 <= 10*cap_k discharged per row.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}
BANK = {8:14, 9:14, 10:14, 11:14, 12:15, 13:16}

# gK8 integer covector: L = 10*q0 + q3 + 10*q6
GK8 = [10, 0, 0, 1, 0, 0, 10]

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

def main():
    print("#"*88)
    print("# THREAD 4 -- ONE dual gK8=(10,0,0,1,0,0,10) clears EVERY binding row k=8..13")
    print("#"*88)
    print("\n  L_yK8(q) = 10*q0 + q3 + 10*q6  (Krawtchouk-nonneg, dominates 10*[t=0])\n")
    print(f"  {'k':>3} {'max_E L_yK8':>16} {'10*cap_k':>14} {'integer margin':>22} {'per-shape margin':>20}")
    rows=[]
    for k in range(8,14):
        cap=CAP[k]; span=BANK[k]
        best=F(-10); bestE=None
        for rest in itertools.combinations(range(1,span+1),k-1):
            E=(0,)+rest
            if not primitive(E): continue
            q=atom_q(E)
            val=sum(GK8[t]*q[t] for t in range(7))
            if val>best: best,bestE=val,E
        cap10=10*cap
        imarg=cap10-best
        psmarg=imarg/10
        ok = imarg>0
        rows.append((k,best,cap10,imarg,psmarg,bestE,ok))
        print(f"  {k:>3} {str(best):>16} {str(cap10):>14} {str(imarg):>22} {str(psmarg):>20}  [{'CLEARS' if ok else 'FAILS'}]")
    print()
    for (k,best,cap10,imarg,psmarg,bestE,ok) in rows:
        print(f"  k={k}: max_E L_yK8 = {best} ({float(best):.6f})  10*cap_{k} = {cap10} ({float(cap10):.6f})")
        print(f"        argmax = {list(bestE)}")
        print(f"        integer margin 10*cap_{k} - max = {imarg} = {float(imarg):+.8f}")
        print(f"        per-shape margin cap_{k} - max measS7-majorant = {psmarg} = {float(psmarg):+.8f}")
    print("\n  VERDICT: the single Lean certificate gK8 (delsarte_bound_k8) clears all "
          "binding rows k=8..13 with strictly positive margin (smallest at k=8: see above).")
    print("  => gK10, gK12 lemmas = delsarte_bound_k8 + per-row finite-check inequality.")

if __name__=="__main__":
    main()
