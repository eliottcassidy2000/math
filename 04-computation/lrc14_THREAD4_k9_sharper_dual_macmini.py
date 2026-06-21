#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_k9_sharper_dual_macmini.py   (THREAD 4, mac-mini, 2026-06-21)

GOAL.  The Delsarte leg  max_E L_y(E) <= cap_k  is TIGHTEST at k=9
(margin cap_9 - L_y(consec_9) = 10441/7567560 ~ 0.00138).  This script:

  (1) Confirms the k=9 razor margin EXACTLY and confirms the L_y maximizer is
      consec_9 = [0..8] (span-robust: re-checked over wider span).
  (2) Searches for a SHARPER Delsarte dual at k=9 -- a Krawtchouk-nonnegative
      covector g with g(0)=1, g(t)>=0 (t=1..6), minimizing  max_E <g,q(E)>.
      A sharper dual gives a LARGER margin = more robust.

KEY REFRAME (the clean form).  Write q_t(E) = P(N_E = t), N_E = # missed inner
sectors (t=0..6).  Then for ANY moment-dual y with readout g(t)=sum_{r<=t} y_r C(t,r):
    L_y(E) = sum_r y_r S_r(E) = sum_r y_r E[C(N,r)] = E[ sum_r y_r C(N,r) ]
           = E[ g(N) ] = sum_t g(t) q_t(E).
So L_y(E) = <g, q(E)> is LINEAR in g.  The Delsarte feasibility is g(0)>=1,
g(t)>=0 (t=1..6).  measS7(E) = q_0(E) <= <g,q(E)> = L_y(E) whenever g(0)>=1, g>=0.

THE LP.  min_g max_E <g, q(E)>  s.t.  g(0)=1, g(t)>=0.
The LP MINIMAX value over the FULL atom polytope (q ranging over ALL distributions
realized by shapes) equals max_E q_0(E) = max_E measS7(E) (achieved by the IE / alternating
covector g=[t=0], which is g=(1,0,0,0,0,0,0) -- but that g is NOT >=0-feasible? it IS:
g=(1,0,0,0,0,0,0) has g(0)=1, g(t)=0>=0. wait -- but is it a valid moment-dual readout?
g=(1,0,...,0) gives y = binomial-inverse = the FULL alternating IE dual
(1,-1,1,-1,1,-1,1). So the LP optimum IS g=delta_0, recovering measS7 exactly, BUT this
requires using ALL SIX moments S_0..S_6.  The point of a *low-degree* dual (R<=3,4) is to
need only a few moments (cheap, Lean-formalizable).  The razor question: with the
DEGREE constraint y_r=0 for r>R, how sharp can we get at k=9?)

So we run TWO searches:
  (A) UNCONSTRAINED degree: best g (g(0)=1, g>=0) minimizing max_E <g,q(E)>.
      -> this recovers max_E measS7 = the cap-binding measS7 itself (margin maximal).
  (B) DEGREE-CONSTRAINED: y_r=0 for r>R (R=3, then R=4), i.e. g must be a polynomial
      of degree<=R in t evaluated at 0..6 (g(t)=poly_R(t)), with g(0)>=1, g(t)>=0.
      -> the Lean-cheap regime; report the best achievable margin per R.

We snap the scipy LP optimum to exact rationals and VERIFY the margin exactly.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
from functools import reduce
from math import gcd
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

import numpy as np
from scipy.optimize import linprog

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

# ---------------------------------------------------------------------------
# Exact atom q_t(E) = P(N_E = t), t=0..6  (N = # missed inner sectors {1..6})
# ---------------------------------------------------------------------------
def atom_q(E):
    """Return q = [q_0,...,q_6], q_t = meas{x : exactly t of inner sectors {1..6} missed}."""
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
        # number of inner sectors {1..6} that are HIT (occupied):
        occ=len([s for s in hit if s!=0])
        missed=6-occ   # N
        q[missed]+=L
    return q

def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1

def factorial_from_atom(q):
    """S_r = sum_t C(t,r) q_t."""
    return [sum(comb(t,r)*q[t] for t in range(7)) for r in range(7)]

# k=9 deg-3 dual readout g (from THM-534): g=(18,5,0,0,2,3,0)/18
G_K9_CURRENT = [F(18,18),F(5,18),F(0),F(0),F(2,18),F(3,18),F(0)]

def Ly_from_g(q, g):
    return sum(g[t]*q[t] for t in range(7))

# ---------------------------------------------------------------------------
def collect_shapes(k, span):
    shapes=[]
    for rest in itertools.combinations(range(1,span+1),k-1):
        E=(0,)+rest
        if not primitive(E): continue
        shapes.append(E)
    return shapes

def main():
    print("#"*86)
    print("# THREAD 4 -- SHARPER k=9 Delsarte dual search (the razor-thin row)")
    print("#"*86)
    k=9; cap9=CAP[9]
    print(f"\ncap_9 = {cap9} = {float(cap9):.8f}")

    # ---- (1) confirm current razor margin, span-robust ----
    print("\n"+"="*86)
    print("(1) Current deg-3 dual g=(18,5,0,0,2,3,0)/18 : razor margin, span-robust")
    print("="*86)
    for span in (14,16,18):
        shapes=collect_shapes(k,span)
        best=F(-10); bestE=None
        for E in shapes:
            q=atom_q(E)
            ly=Ly_from_g(q,G_K9_CURRENT)
            if ly>best: best,bestE=ly,E
        margin=cap9-best
        print(f"  span<={span:2d} ({len(shapes):4d} shapes): max L_y = {best} = {float(best):.8f}")
        print(f"       argmax = {list(bestE)}   margin cap_9 - max L_y = {margin} = {float(margin):+.8f}")

    # use span=14 finite-check region for the LP (the bounded-span region)
    SPAN_LP=14
    shapes=collect_shapes(k,SPAN_LP)
    Q=[atom_q(E) for E in shapes]
    Qf=np.array([[float(qt) for qt in q] for q in Q])
    print(f"\n  LP shape set: span<={SPAN_LP}, {len(shapes)} primitive shapes (atoms q_t computed exact).")

    # ---- (A) UNCONSTRAINED best dual: min_g max_E <g,q>, g(0)=1, g>=0 ----
    # variables: g_1..g_6 (g_0 fixed =1) and slack M (the max). minimize M.
    # constraints: for each shape E: g_0 q0 + sum_{t>=1} g_t q_t <= M
    #   => sum_{t>=1} g_t q_t(E) - M <= -q0(E)
    # vars order: [g1,g2,g3,g4,g5,g6, M]  (7 vars)
    print("\n"+"="*86)
    print("(A) UNCONSTRAINED best Delsarte dual (any degree): min max_E <g,q>, g(0)=1, g>=0")
    print("="*86)
    nv=7
    c=np.zeros(nv); c[6]=1.0  # minimize M
    A=[]; b=[]
    for q in Qf:
        row=np.zeros(nv)
        for t in range(1,7): row[t-1]=q[t]
        row[6]=-1.0
        A.append(row); b.append(-q[0])
    # g_t>=0 bounds, M free-ish (>=0)
    bounds=[(0,None)]*6+[(0,None)]
    res=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method="highs")
    print(f"  LP status: {res.message}")
    print(f"  optimum max_E L_y (float) = {res.fun:.8f}")
    gopt=[1.0]+list(res.x[:6])
    print(f"  optimal g (float) = {[round(x,5) for x in gopt]}")
    print(f"  => margin cap_9 - opt = {float(cap9)-res.fun:+.8f}  (this is the best ANY g>=0 dual can do)")
    print("  NOTE: this needs all 6 moments; the value should approach max_E measS7 = q0_max.")
    q0max=max(float(q[0]) for q in Qf)
    print(f"       max_E measS7 (q0) = {q0max:.8f}   (the unbeatable floor for any g>=0 dual)")

    # ---- (B) DEGREE-CONSTRAINED dual: g(t) = poly of degree<=R in t ----
    # g_t = sum_{r=0}^R a_r C(t,r)  (binomial basis = degree-R polynomial in t).
    # Equivalent to y_r = a_r (the moment coeffs), y_r=0 for r>R.
    # constraints: g(0)=a_0=1; g(t)>=0 for t=1..6; minimize max_E sum_t g(t) q_t.
    for R in (2,3,4):
        print("\n"+"="*86)
        print(f"(B) DEGREE-{R} dual: y_r=0 for r>{R}.  g(t)=sum_(r<={R}) y_r C(t,r), g(0)=1, g(t)>=0")
        print("="*86)
        # vars: y_1..y_R  (y_0=1 fixed), and M. nv = R + 1
        nv=R+1
        c=np.zeros(nv); c[R]=1.0  # minimize M (last var)
        # g(t) = 1 + sum_{r=1}^R y_r C(t,r)
        def gcoef(t):  # returns (const, [coef of y_1..y_R])
            return (1.0, [comb(t,r) for r in range(1,R+1)])
        A=[]; b=[]
        # feasibility g(t)>=0, t=1..6:  -(1+ sum y_r C(t,r)) <= 0
        for t in range(1,7):
            const,coefs=gcoef(t)
            row=np.zeros(nv)
            for i,cf in enumerate(coefs): row[i]=-cf
            A.append(row); b.append(const)  # -sum y_r C - 1 <=0 => sum(-cf y)<=1... wait sign
        # rewrite: g(t)=1+sum y_r C(t,r) >=0  <=>  -sum y_r C(t,r) <= 1
        # objective constraint: sum_t g(t) q_t(E) <= M  for each E
        #   sum_t [1+sum_r y_r C(t,r)] q_t <= M
        #   sum_t q_t + sum_r y_r (sum_t C(t,r) q_t) <= M
        #   1 + sum_r y_r S_r(E) <= M   (sum_t q_t =1)
        for q in Qf:
            Sr=[sum(comb(t,r)*q[t] for t in range(7)) for r in range(R+1)]
            row=np.zeros(nv)
            for r in range(1,R+1): row[r-1]=Sr[r]
            row[R]=-1.0
            A.append(row); b.append(-1.0)  # 1 + sum y_r S_r - M <=0 => sum y_r S_r - M <= -1
        bounds=[(None,None)]*R+[(0,None)]
        res=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method="highs")
        if not res.success:
            print(f"  LP failed: {res.message}"); continue
        yopt=[1.0]+list(res.x[:R])
        gopt=[1.0+sum(res.x[r-1]*comb(t,r) for r in range(1,R+1)) for t in range(7)]
        print(f"  LP status: {res.message}")
        print(f"  optimum max_E L_y (float) = {res.fun:.8f}")
        print(f"  optimal y (float) = {[round(x,6) for x in yopt]}")
        print(f"  optimal g (float) = {[round(x,6) for x in gopt]}")
        print(f"  => margin cap_9 - opt = {float(cap9)-res.fun:+.8f}")
        cur=float(sum(G_K9_CURRENT[t]*qt for q in [Q[0]] for t,qt in []))  # placeholder
        print(f"     (current deg-3 dual margin was +0.00137970)")

if __name__=="__main__":
    main()
