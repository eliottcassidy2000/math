#!/usr/bin/env python3
"""
LRC <-> additive combinatorics (mac-mini-2026-07-03-S36): explore BDFKK / Goncalves-Radchenko connections.
The safe measure mu(S)=meas{t: all ||v_i t||>=1/14} = sum_{m in L(S)} prod_i c(m_i), where
 L(S)={m in Z^13: sum m_i v_i = 0} (the RELATION LATTICE), c(0)=6/7, c(m)=-hat{1_D}(m), hat{1_D}(m)=sin(pi m/7)/(pi m).
 => TIGHT (mu=0) <=> the relation-lattice theta-sum EXACTLY cancels the main term (6/7)^13.
Tests:
 (1) additive energy E(S)=#{(a,b,c,d) in S^4: a+b=c+d}, and small-relation counts, for AP/GW/deep-well/covering.
 (2) the theta-sum mu(S) truncated to |m_i|<=K (verify mu(AP)~0, mu(GW)~0, mu(covering)>0), and how fast it converges.
 (3) does 'tight <=> maximal additive energy' hold? (the inverse-theorem framing of GAP-A).
"""
from fractions import Fraction as F
from math import sin, pi, gcd
from functools import reduce
from itertools import product
import numpy as np

def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
def additive_energy(S):
    S = list(S); from collections import Counter
    cnt = Counter(a+b for a in S for b in S)
    return sum(c*c for c in cnt.values())
def rep_relations(S, K):
    """count m in {-K..K}^len(S), m!=0, sum m_i S_i = 0 (small relations)."""
    # too big to brute for 13 dims at K>=2; instead count via meet-in-the-middle on the multiset of sums
    # count solutions to sum m_i v_i=0 with |m_i|<=K: use generating/convolution over a bounded range.
    from collections import Counter
    half = len(S)//2
    A = S[:half]; B = S[half:]
    ca = Counter()
    for m in product(range(-K,K+1), repeat=len(A)):
        ca[sum(mi*ai for mi,ai in zip(m,A))] += 1
    total = 0
    for m in product(range(-K,K+1), repeat=len(B)):
        total += ca[-sum(mi*bi for mi,bi in zip(m,B))]
    return total - 1  # subtract the all-zero
def mu_theta(S, K):
    """truncated theta-sum mu ~ sum_{|m_i|<=K, sum m_i S_i=0} prod c(m_i), c(0)=6/7, c(m)=-sin(pi m/7)/(pi m)."""
    def c(m):
        if m == 0: return 6.0/7.0
        if m % 7 == 0: return 0.0
        return -sin(pi*m/7)/(pi*m)
    from collections import defaultdict
    half = len(S)//2
    A = S[:half]; B = S[half:]
    # accumulate partial products by partial sum
    da = defaultdict(float)
    for m in product(range(-K,K+1), repeat=len(A)):
        p = 1.0
        for mi in m: p *= c(mi)
        if p != 0.0: da[sum(mi*ai for mi,ai in zip(m,A))] += p
    tot = 0.0
    for m in product(range(-K,K+1), repeat=len(B)):
        p = 1.0
        for mi in m: p *= c(mi)
        if p == 0.0: continue
        s = sum(mi*bi for mi,bi in zip(m,B))
        tot += p * da.get(-s, 0.0)
    return tot

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    AP = list(range(1,14))
    GW = [1,2,3,4,5,6,7,8,9,10,11,13,24]
    DW = list(range(1,13))+[182]
    COV = [2,4,6,8,10,12,14,3,9,15,5,11,13]  # a covering-ish primitive family
    fams = [("AP", AP), ("GW", GW), ("deep well", DW), ("2..14", list(range(2,15))), ("cov-ish", COV)]
    out("(1) additive energy E(S) and small-relation counts (K=1: {-1,0,1}^13):")
    out(f"{'family':>11} {'M':>8} {'tight?':>6} {'add.energy':>11} {'#rel |m|<=1':>12}")
    for name, S in fams:
        M = M_exact(S)
        E = additive_energy(S)
        r1 = rep_relations(S, 1)
        out(f"{name:>11} {float(M):>8.5f} {str(M==F(1,14)):>6} {E:>11} {r1:>12}")
    out("\n(2) truncated theta-sum mu(S) at K=1,2 (tight => mu->0; loose => mu>0):")
    out(f"{'family':>11} {'muK1':>10} {'muK2':>10} {'(6/7)^13':>10}")
    base = (6.0/7.0)**13
    for name, S in fams:
        m1 = mu_theta(S, 1); m2 = mu_theta(S, 2)
        out(f"{name:>11} {m1:>10.5f} {m2:>10.5f} {base:>10.5f}")
    out("\n=> if tight (AP,GW,DW) have HIGH additive energy + theta-sum -> 0, and loose have lower energy + mu>0:")
    out("   GAP-A = inverse theorem (mu-extremal => AP); GAP-B = additive-energy/relation-lattice bound (BDFKK/GR).")
