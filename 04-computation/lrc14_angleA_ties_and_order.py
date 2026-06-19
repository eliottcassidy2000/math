#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE A part 2: characterize ties (who equals consec) and find the working order.

Findings from part 1: consec maximizes p_0 AND L_y exactly (0 violations), but
worst gap = 0 => there are TIES. Hypotheses:
  - ties on L_y are exactly DILATIONS c*consec (scale invariance THM-531).
  - but is p_0 also tied only by dilations? Or do non-dilation sets tie p_0?

Also: since usual SD fails (dilations break it), test whether the right object is
NOT N_E but the *factorial moments* S_r, and whether consec maximizes EACH S_r
that g uses. g(t)=sum_r y_r C(t,r), L_y = sum_r y_r S_r. So L_y<=L_y(consec) would
follow from S_r(E) <= S_r(consec) for every r with y_r>0 -- IF y_r>0 for all used r.
Check the SIGN of y_r (the moment weights) for k=8,9.

kind-pasteur-2026-06-19 ANGLE-A.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

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
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            val = Fraction((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = Fraction(-(t-2)*(t-3)*(t-6), 36)
        else:
            val = Fraction((t-3)*(t-4), 12)
        g.append(val)
    return g

def binom(n,r):
    from math import comb
    return comb(n,r)

def moments_from_p(p):
    """S_r = E[C(N,r)] = sum_t p_t C(t,r), r=0..6."""
    return [sum(p[t]*binom(t,r) for t in range(7)) for r in range(7)]

def y_from_g(g):
    """recover y_r so that g(t)=sum_r y_r C(t,r). Newton forward differences:
       y_r = (Delta^r g)(0) = sum_{j=0}^r (-1)^{r-j} C(r,j) g(j)."""
    y=[]
    for r in range(7):
        val=Fraction(0)
        for j in range(r+1):
            val += (-1)**(r-j)*binom(r,j)*g[j]
        y.append(val)
    return y

def L_y(p,k):
    g=g_poly(k); return sum(p[t]*g[t] for t in range(7))

def gen_competitors(k, maxspread):
    out=[]
    for combo in itertools.combinations(range(1,maxspread+1), k-1):
        E=[0]+list(combo)
        if reduce(gcd,E)!=1: continue
        out.append(E)
    return out

if __name__=="__main__":
    for k in [8,9]:
        print(f"\n{'='*70}\nk={k}\n{'='*70}")
        g=g_poly(k)
        y=y_from_g(g)
        print(f"g={[str(x) for x in g]}")
        print(f"y_r (moment weights, L_y=sum y_r S_r): {[str(x) for x in y]}")
        print(f"  signs of y_r: {[ '+' if v>0 else ('-' if v<0 else '0') for v in y]}")
        used=[r for r in range(7) if y[r]!=0]
        print(f"  used moments r in {used}")
        C=list(range(k))
        pc=dist_p(C); Sc=moments_from_p(pc); Lc=L_y(pc,k)
        print(f"  consec moments S_r = {[str(x) for x in Sc]}")
        # verify L_y = sum y_r S_r
        chk=sum(y[r]*Sc[r] for r in range(7))
        print(f"  check sum y_r S_r = {chk} vs L_y = {Lc}  {'OK' if chk==Lc else 'MISMATCH'}")

        comps=gen_competitors(k, {8:14,9:13}[k])
        # for each used moment, does consec maximize it?
        mom_max = {r:True for r in used}
        mom_beats = {r:0 for r in used}
        # ties on L_y: are they all dilations?
        nondil_ties=[]
        Cset=set(C)
        def is_dilation(E):
            E=sorted(E)
            if E[0]!=0: return False
            g0=reduce(gcd,E)
            return [e//g0 for e in E]==C
        for E in comps:
            if E==C: continue
            p=dist_p(E); S=moments_from_p(p); L=L_y(p,k)
            for r in used:
                if S[r]>Sc[r]:
                    mom_max[r]=False; mom_beats[r]+=1
            if L==Lc and not is_dilation(E):
                nondil_ties.append(E)
        print(f"  --- per-moment: does consec maximize S_r? ---")
        for r in used:
            print(f"    S_{r}: y_{r}={'+' if y[r]>0 else '-'}  consec-max={'YES' if mom_max[r] else f'NO ({mom_beats[r]} beat)'}")
        print(f"  non-dilation ties on L_y: {len(nondil_ties)}")
        for E in nondil_ties[:10]:
            print(f"     {E}")
