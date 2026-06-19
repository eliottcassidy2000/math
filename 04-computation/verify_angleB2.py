#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL part 2: verify the EXPONENT (k-1), the J(A) Fourier baseline,
and that L_y converges to L_inf for dissociated sets, using the INDEPENDENT
breakpoint engine.
"""
import sys, itertools
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb

# ---- independent re-implementation of the breakpoint engine ----
def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        p[N_at(E, mid)] += (hi-lo)
    return p

def g_poly(k):
    g=[]
    for t in range(7):
        if k==8: g.append(F((t-1)*(t-2)*(t-4)*(t-5),40))
        elif k in (9,10): g.append(F(-(t-2)*(t-3)*(t-6),36))
        else: g.append(F((t-3)*(t-4),12))
    return g

def L_y_exact(E,k):
    p=dist_p(E); g=g_poly(k)
    return sum(p[t]*g[t] for t in range(7)), p

def y_coeffs(k):
    g=g_poly(k); y=[]
    for r in range(7):
        y.append(sum((-1)**(r-j)*comb(r,j)*g[j] for j in range(r+1)))
    return y

def L_inf(k):
    y=y_coeffs(k)
    return sum(y[r]*comb(6,r)*F(7-r,7)**(k-1) for r in range(7))

# ---------------------------------------------------------------
# DIRECT exponent test via factorial moments S_r and the J(A) baseline.
# Compute S_r = E[C(N,r)] exactly from dist_p, and compare the moment-based
# L_inf prediction. Also directly compute one J(A) for a strongly dissociated
# set and compare to (1-r/7)^(k-1) vs (1-r/7)^k.
# ---------------------------------------------------------------

def S_moments(E):
    p=dist_p(E)
    return [sum(p[t]*comb(t,r) for t in range(7)) for r in range(7)]

# J(A): measure that every generator avoids ALL sectors in set A (A subset of {1..6})
def J_exact(E, A):
    E=sorted(set(E)); Aset=set(A)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        # avoid A iff no generator lands in any sector of A
        ok=True
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator)
            s=(v.numerator*7)//v.denominator
            if s in Aset: ok=False; break
        if ok: tot+=(hi-lo)
    return tot

print("=== EXPONENT test: J(A) for a strongly dissociated set vs (1-r/7)^? ===")
# Build a dissociated set: powers-like / large gaps to kill small integer relations.
# k=8 example: 0 plus 7 generators with big multiplicative gaps.
for k, gens in [(8,[0,1,11,113,1129,3457,7919,9973])]:
    print(f"k={k}, E={gens}")
    for r in [1,2,3]:
        A=tuple(range(1,r+1))  # any A with |A|=r
        J=J_exact(gens,A)
        base_km1=F(7-r,7)**(k-1)
        base_k  =F(7-r,7)**k
        print(f"  r={r}, A={A}: J={float(J):.6f}  (1-r/7)^(k-1)={float(base_km1):.6f}  (1-r/7)^k={float(base_k):.6f}")
        print(f"        |J-(k-1)form|={float(abs(J-base_km1)):.6f}   |J-(k)form|={float(abs(J-base_k)):.6f}")
print()

print("=== L_y convergence to L_inf along increasingly dissociated sets ===")
for k in [8,9]:
    Li=L_inf(k)
    print(f"k={k}: L_inf={float(Li):.6f}")
    seqs = {
        'consec': list(range(k)),
        'mild':  [0]+[2**i+1 for i in range(k-1)],
    }
    # geometric dissociated families with growing base (capped to keep engine fast)
    for base in [3,5]:
        E=[0]+[base**i for i in range(k-1)]
        if max(E) <= 80000:
            seqs[f'geo{base}']=E
    for name,E in seqs.items():
        try:
            Ly,_=L_y_exact(E,k)
            print(f"   {name:8s} E={E}  L_y={float(Ly):.6f}  diff_to_inf={float(Ly-Li):+.6f}")
        except Exception as ex:
            print(f"   {name}: skipped ({ex})")
    print()
