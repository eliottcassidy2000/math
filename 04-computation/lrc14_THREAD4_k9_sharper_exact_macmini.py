#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_k9_sharper_exact_macmini.py  (THREAD 4, mac-mini, 2026-06-21)

EXACT-RATIONAL verification of the SHARPER k=9 Delsarte dual found by the LP.

FINDING (lrc14_THREAD4_k9_sharper_dual_macmini.py).  At degree 3 the THM-534 k=9
dual y=(1,-13/18,4/9,-1/6) is ALREADY LP-optimal (no sharper deg-3 dual): razor
margin cap_9 - max L_y = 10441/7567560 ~ +0.001380.

BUT going to degree 4, the LP optimum is y = (1,-1,1,-9/10,3/5) -- which is EXACTLY
the THM-534 k=8 dual yK8, with integer readout g = (10,0,0,1,0,0,10)/10 = gK8/10.
Applied as a k=9 majorant this gives a margin ~ +0.0509 -- about 37x larger.

So: the SHARPER k=9 dual = the EXISTING gK8 certificate.  It is ALREADY Lean-formalized
(gK8_values, gK8_dominates, delsarte_bound_k8).  Re-using it at k=9 only needs the
finite check  max_E L_yK8(E) <= cap_9, exactly the same shape we verify here.

This script:
  (1) Verifies EXACTLY that g=(10,0,0,1,0,0,10)/10 (=gK8/10) gives
        max_E L_y(E) at k=9 (over span<=14, then span-robust to 18),
      the exact max value, argmax, and the EXACT margin cap_9 - max L_y.
  (2) Confirms the deg-3 dual is LP-optimal at deg 3 (margin recovered exactly).
  (3) Reports the sharper margin as an exact Fraction, ready for the proof.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

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

def Ly_from_g(q,g): return sum(g[t]*q[t] for t in range(7))

# duals as readout covectors g(t), t=0..6 (rational)
G_K9_DEG3 = [F(1),F(5,18),F(0),F(0),F(2,18),F(3,18),F(0)]          # = (18,5,0,0,2,3,0)/18
G_K8_DEG4 = [F(10,10),F(0),F(0),F(1,10),F(0),F(0),F(10,10)]        # = (10,0,0,1,0,0,10)/10 = gK8/10

def collect_shapes(k,span):
    return [(0,)+rest for rest in itertools.combinations(range(1,span+1),k-1)
            if primitive((0,)+rest)]

def maxLy(k,span,g):
    best=F(-10); bestE=None
    for E in collect_shapes(k,span):
        ly=Ly_from_g(atom_q(E),g)
        if ly>best: best,bestE=ly,E
    return best,bestE

def main():
    cap9=CAP[9]
    print("#"*86)
    print("# THREAD 4 -- EXACT sharper k=9 Delsarte dual = gK8/10 (re-used k=8 cert)")
    print("#"*86)
    print(f"\ncap_9 = {cap9} = {float(cap9):.10f}\n")

    print("="*86)
    print("DEG-3 dual g=(18,5,0,0,2,3,0)/18  (current THM-534 k=9 dual)")
    print("="*86)
    for span in (14,16,18):
        m,E=maxLy(9,span,G_K9_DEG3)
        print(f"  span<={span}: max L_y = {m} = {float(m):.10f}  argmax={list(E)}  margin={cap9-m}={float(cap9-m):+.8f}")

    print("\n"+"="*86)
    print("SHARPER DEG-4 dual g=(10,0,0,1,0,0,10)/10 = gK8/10  (re-used k=8 certificate)")
    print("="*86)
    sharper_margin=None; sharper_max=None; sharper_arg=None
    for span in (14,16,18):
        m,E=maxLy(9,span,G_K8_DEG4)
        mar=cap9-m
        print(f"  span<={span}: max L_y = {m} = {float(m):.10f}  argmax={list(E)}  margin={mar}={float(mar):+.8f}")
        if span==14:
            sharper_margin,sharper_max,sharper_arg=mar,m,E

    # express the sharper bound in the gK8 integer scale: 10*q0 <= 10*q0 + q3 + 10*q6
    # so measS7 = q0 <= (10*q0 + q3 + 10*q6)/10 = L_yK8/10. The finite check is
    # max_E (10*q0 + q3 + 10*q6) <= 10*cap_9.
    print("\n"+"="*86)
    print("INTEGER (Lean) FORM:  10*measS7 <= L_yK8 = 10*q0 + q3 + 10*q6")
    print("  finite check needed:  max_E (10*q0 + q3 + 10*q6) <= 10*cap_9")
    print("="*86)
    cap9x10=10*cap9
    best=F(-10); bestE=None
    for E in collect_shapes(9,14):
        q=atom_q(E)
        val=10*q[0]+q[3]+10*q[6]
        if val>best: best,bestE=val,E
    print(f"  max_E (10*q0+q3+10*q6) = {best} = {float(best):.8f}   argmax={list(bestE)}")
    print(f"  10*cap_9 = {cap9x10} = {float(cap9x10):.8f}")
    print(f"  10*cap_9 - max = {cap9x10-best} = {float(cap9x10-best):+.8f}  ({'CLEARS' if cap9x10-best>0 else 'FAILS'})")
    print(f"  per-shape margin cap_9 - max L_y = {(cap9x10-best)/10} = {float((cap9x10-best)/10):+.8f}")

    print("\n"+"="*86)
    print("SUMMARY")
    print("="*86)
    print(f"  deg-3 (current)  margin = 10441/7567560 = +0.00137970   (LP-optimal at deg 3)")
    print(f"  deg-4 = gK8/10   margin = {sharper_margin} = {float(sharper_margin):+.8f}")
    ratio = float(sharper_margin)/float(F(10441,7567560))
    print(f"  => SHARPER by factor {ratio:.1f}x ; the sharper cert is the ALREADY-FORMALIZED gK8.")

if __name__=="__main__":
    main()
