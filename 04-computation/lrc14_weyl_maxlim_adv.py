#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL part 2: pin down max_{B0} LIM(B0) and the worst finite-w single-far meas.
Key questions:
 (1) Reproduce maxLIM_k exactly. Is the maximizing base consec(k-1)?
 (2) Does the bounded near-AP base [0..k-2, k] give a HIGHER LIM (the prompt's worst non-AP)?
 (3) Most important: is sup over single-far sets (any base in a window, any large w)
     STRICTLY below M_k?  Search finite-w meas over wide window & many w.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def missed(E, x):
    hit = set()
    for e in E:
        v = e*x; v = v - (v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return [j for j in range(1,7) if j not in hit]

def dist(E):
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        t = len(missed(E, mid))
        p[t] += (hi-lo)
    return p

def meas(E): return dist(E)[0]
def primitive(E):
    nz=[e for e in E if e!=0]; return reduce(gcd,nz)==1 if nz else False
def LIM(B0):
    p=dist(B0); return p[0]+p[1]/7

if __name__=="__main__":
    for k in (8,9):
        Mk = meas(list(range(k)))
        print(f"########## k={k}  M_k = {Mk} = {float(Mk):.6f} ##########")

        # (1)+(2): max LIM over all primitive (k-1)-bases in window {0..W}
        W = 14
        rest = list(range(1, W+1))
        bestlim = Fraction(0); bestB = None
        # also track consec-base and near-AP base explicitly
        consec_base = list(range(k-1))
        nearAP_base = list(range(k-1)) ; nearAP_base[-1]=k-1  # [0..k-3,k-1]
        for combo in itertools.combinations(rest, k-2):
            B0 = (0,)+combo
            # base need not be primitive (far w fixes), but LIM is base-only;
            # for a valid k-set the base+w can always be made primitive, so consider all.
            l = LIM(list(B0))
            if l > bestlim:
                bestlim = l; bestB = list(B0)
        gap = Mk - bestlim
        print(f"  maxLIM over bases in {{0..{W}}} = {float(bestlim):.6f} at {bestB}")
        print(f"  LIM(consec base {consec_base}) = {float(LIM(consec_base)):.6f}")
        print(f"  M_k - maxLIM = {float(gap):+.6f}   (need > 0; reduction target)")
        assert gap > 0, "REDUCTION FAILS: a base LIM meets/exceeds M_k!"

        # (3): worst finite-w single-far meas, wide window & many w
        worst = Fraction(0); worstE=None; viol=[]
        wlist = [37,53,101,211,307,503,809,1009,2003]
        base_W = 11
        for combo in itertools.combinations(range(1,base_W+1), k-2):
            B0 = [0]+list(combo)
            for w in wlist:
                if w <= max(B0): continue
                E = B0+[w]
                if not primitive(E): continue
                m = meas(E)
                if m > worst: worst=m; worstE=list(E)
                if m >= Mk: viol.append((float(m),list(E)))
        print(f"  worst single-far meas (window<= {base_W}, w<= {max(wlist)}) = {float(worst):.6f} at {worstE}")
        print(f"    margin to M_k = {float(Mk-worst):+.6f}")
        if viol:
            print(f"    *** {len(viol)} single-far sets MEET/EXCEED M_k: {viol[:3]}")
        else:
            print(f"    NO single-far set reached M_k. reduction holds on grid.")
        print()
