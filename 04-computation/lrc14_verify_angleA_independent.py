#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of the 'convex-order / consec maximizes L_y' claim.
Own breakpoint engine (Fractions), own moment computation, own searches.
"""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd

# ---------- independent engine ----------
def dist_p(E):
    """exact distribution p_t (t=0..6) of N(x)= # of inner sectors 1..6 missed by orbit."""
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
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
        hit = set()
        for e in E:
            v = e*mid
            v = v - (v.numerator//v.denominator)
            s = (v.numerator*7)//v.denominator
            hit.add(s)
        N = sum(1 for j in range(1,7) if j not in hit)
        p[N] += (hi-lo)
    return p

def C_binom(n, r):
    from math import comb
    return comb(n, r) if 0 <= r <= n else 0

def moments_from_p(p):
    """S_r = E[C(N,r)], r=0..6."""
    return [sum(p[t]*C_binom(t, r) for t in range(7)) for r in range(7)]

def g_poly(k):
    g=[]
    for t in range(7):
        if k==8: val=Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
        elif k in (9,10): val=Fraction(-(t-2)*(t-3)*(t-6),36)
        else: val=Fraction((t-3)*(t-4),12)
        g.append(val)
    return g

def L_y(E,k):
    p=dist_p(E); g=g_poly(k)
    return sum(p[t]*g[t] for t in range(7)), p

def p0_inclusion_exclusion(p):
    """check p_0 = sum_r (-1)^r S_r."""
    S=moments_from_p(p)
    return sum((-1)**r * S[r] for r in range(7))

def consec(k): return list(range(k))

# y_r weights claimed
Y = {
 8:[Fraction(1),Fraction(-1),Fraction(1),Fraction(-9,10),Fraction(3,5),Fraction(0),Fraction(0)],
 9:[Fraction(1),Fraction(-13,18),Fraction(4,9),Fraction(-1,6),Fraction(0),Fraction(0),Fraction(0)],
}

def gen_sets(k, maxval):
    """all E with 0 in E, |E|=k, max element <= maxval, sorted distinct."""
    for rest in itertools.combinations(range(1, maxval+1), k-1):
        yield [0]+list(rest)

if __name__=="__main__":
    import sys
    sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

    # ---- claim (2): exact p0 and L_y values for consec ----
    print("=== exact consec values ===")
    for k in [8,9]:
        C=consec(k); L,p=L_y(C,k)
        S=moments_from_p(p)
        print(f"k={k}: consec={C}")
        print(f"   p0=meas(S7)={p[0]} = {float(p[0]):.6f}")
        print(f"   L_y={L} = {float(L):.6f}")
        # verify inclusion-exclusion p0 = sum (-1)^r S_r
        p0ie=p0_inclusion_exclusion(p)
        print(f"   incl-excl p0 reconstruct: {p0ie}  match={p0ie==p[0]}")
        # verify y_r weights reproduce L_y: L_y = sum y_r S_r
        if k in Y:
            Ly_from_S=sum(Y[k][r]*S[r] for r in range(7))
            print(f"   sum y_r S_r = {Ly_from_S}  match L_y={Ly_from_S==L}")
            # verify g(t)=sum y_r C(t,r)
            g=g_poly(k); ok=all(g[t]==sum(Y[k][r]*C_binom(t,r) for r in range(7)) for t in range(7))
            print(f"   g(t)=sum y_r C(t,r) for all t: {ok}")

    # ---- claims (1)(2): consec strict max of L_y AND of p0 over bounded spread ----
    print("\n=== exhaustive search: consec max of L_y and p0? ===")
    cfg={8:14, 9:13}
    for k,mx in cfg.items():
        C=consec(k); Lc,pc=L_y(C,k); p0c=pc[0]
        nsets=0; beatL=[]; tieL=[]; beatP0=[]; tieP0=[]
        for E in gen_sets(k,mx):
            if reduce(gcd,E)!=1:  # only primitive; dilations covered by scale-inv
                continue
            nsets+=1
            L,p=L_y(E,k)
            if L>Lc: beatL.append((E,float(L)))
            elif L==Lc and E!=C: tieL.append(E)
            if p[0]>p0c: beatP0.append((E,float(p[0])))
            elif p[0]==p0c and E!=C: tieP0.append(E)
        print(f"k={k} maxval<={mx}: {nsets} primitive sets")
        print(f"   L_y beats={len(beatL)} ties(non-consec)={len(tieL)}")
        if beatL: print("      BEAT L:",beatL[:5])
        if tieL: print("      TIE L:",tieL[:5])
        print(f"   p0  beats={len(beatP0)} ties(non-consec)={len(tieP0)}")
        if beatP0: print("      BEAT p0:",beatP0[:5])
        if tieP0: print("      TIE p0:",tieP0[:5])
