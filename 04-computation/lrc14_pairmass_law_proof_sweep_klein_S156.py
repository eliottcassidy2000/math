#!/usr/bin/env python3
"""
klein-2026-07-07-S156 -- THM-638 verification sweep: the SIGNED PAIR-MASS LAW (general
rational threshold), the triple-mass law probe, and the Bonferroni-3 endpoint floor at k=8.

THE LAW (derivation in THM-638; this script is the exhaustive verification):
  For coprime q1,q2 >= 1 and theta = c/s in lowest terms (0 < c/s < 1), with
  r_i := c*q_i mod s in {0..s-1}:
   SAME-SIGN:  meas{frac(q1 x), frac(q2 x) in (0, theta]} = theta^2 + G+/(s^2 q1 q2),
               G+ = s*min(r1,r2) - r1 r2 = min(r)*(s - max(r))  >= 0 ALWAYS.
   MIXED-SIGN: meas{frac(q1 x) in (0,theta], frac(q2 x) in [1-theta,1)} =
               theta^2 + G-/(s^2 q1 q2),  G- = s*max(0, r1+r2-s) - r1 r2  <= 0 ALWAYS.
  Non-coprime (d1,d2) reduce to (q1,q2)=(d1,d2)/g with the SAME m (x -> gx measure-pres).
  Proof method: Bezout offset enumeration + grid-count c(t) = a1 + 1[frac(Nt) < r1/s]
  + integrate the indicator over the second window (two floor identities). Half a page.

CONSEQUENCES CHECKED HERE:
  - G+ = 0 iff s | q1 c or s | q2 c (threshold-resonant directions: the 'apex invisibility')
  - G- = -r1 r2 when r1 + r2 <= s: mixed masses VANISH iff q1 q2 = r1 r2 (small exact cases)
  - k=8 Hunter floor arithmetic: 6*theta^2 = 6/49 (7 same-sign events, 6-edge tree)
  - k=9 CORRECTION (MISTAKE-122): 8 events -> 7-edge tree -> bare floor 1 - 8/7 + 7/49 = 0
    EXACTLY (my S155 '1/49 > 0' used 8 edges -- wrong; trees on n events have n-1 edges).

TRIPLE-MASS LAW PROBE: is m123 = theta^3 + theta*(G12/(q1q2)+G13/(q1q3)+G23/(q2q3)) + H/(q1q2q3)
with H constant on residue classes (r1,r2,r3)?  (exact triple engine = interval lists)

BONFERRONI-3 ENDPOINT FLOOR at k=8 (exact per shape, RIGOROUS unconditionally):
  W_end >= 1 - S1 + S2 - S3 = S2 - S3 (S1 = 1 exactly at k=8);  adversarial min over
  8-sets vs the R-route bar 0.197 and the Hunter floor 6/49 = 0.1224.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def intervals_A(q, theta: F, sign):
    """the hit set for signed difference sign*q: list of (lo,hi] subintervals of [0,1)."""
    out = []
    for j in range(q):
        if sign > 0:
            lo, hi = F(j, q), F(j, q) + theta/q
        else:
            lo, hi = (F(j+1, q) - theta/q), F(j+1, q)
        out.append((lo, hi))
    return out

def m_pair_exact(d1, d2, theta: F):
    """exact pair mass for signed d1,d2 (nonzero ints)."""
    g = gcd(abs(d1), abs(d2)); q1, q2 = abs(d1)//g, abs(d2)//g
    A = intervals_A(q1, theta, 1 if d1 > 0 else -1)
    B = intervals_A(q2, theta, 1 if d2 > 0 else -1)
    tot = F(0)
    for (a0, a1) in A:
        for (b0, b1) in B:
            lo, hi = max(a0, b0), min(a1, b1)
            if hi > lo: tot += hi - lo
    return tot

def law_pair(d1, d2, theta: F):
    """the THM-638 prediction."""
    c, s = theta.numerator, theta.denominator
    g = gcd(abs(d1), abs(d2)); q1, q2 = abs(d1)//g, abs(d2)//g
    r1, r2 = (c*q1) % s, (c*q2) % s
    if (d1 > 0) == (d2 > 0):
        G = s*min(r1, r2) - r1*r2
    else:
        G = s*max(0, r1 + r2 - s) - r1*r2
    return theta*theta + F(G, s*s*q1*q2)

def m_triple_exact(d1, d2, d3, theta: F):
    """exact triple mass (all same sign assumed positive here) via interval lists."""
    A = intervals_A(d1, theta, 1)
    # intersect A with B
    B = intervals_A(d2, theta, 1)
    AB = []
    for (a0, a1) in A:
        for (b0, b1) in B:
            lo, hi = max(a0, b0), min(a1, b1)
            if hi > lo: AB.append((lo, hi))
    C = intervals_A(d3, theta, 1)
    tot = F(0)
    for (a0, a1) in AB:
        for (c0, c1) in C:
            lo, hi = max(a0, c0), min(a1, c1)
            if hi > lo: tot += hi - lo
    return tot

if __name__ == "__main__":
    rng = np.random.default_rng(41560)

    print("=== 1. SAME-SIGN law sweep: all coprime pairs q1<q2<=60, theta=1/7 (exhaustive exact) ===")
    bad = 0; tested = 0
    TH = F(1, 7)
    for q1 in range(1, 61):
        for q2 in range(q1+1, 61):
            if gcd(q1, q2) != 1: continue
            tested += 1
            if m_pair_exact(q1, q2, TH) != law_pair(q1, q2, TH): bad += 1
    print(f"  tested {tested} coprime pairs: {bad} violations")

    print("\n=== 2. MIXED-SIGN law sweep: q1<=40, q2<=40 coprime, theta=1/7 ===")
    bad = 0; tested = 0
    for q1 in range(1, 41):
        for q2 in range(1, 41):
            if gcd(q1, q2) != 1: continue
            tested += 1
            if m_pair_exact(q1, -q2, TH) != law_pair(q1, -q2, TH): bad += 1
    print(f"  tested {tested}: {bad} violations")

    print("\n=== 3. GENERAL theta sweep: theta in {1/4, 1/5, 1/6, 1/9, 2/7, 3/8}, pairs <= 30, both signs ===")
    for theta in (F(1,4), F(1,5), F(1,6), F(1,9), F(2,7), F(3,8)):
        bad = 0; tested = 0
        for q1 in range(1, 31):
            for q2 in range(1, 31):
                if q1 == q2 and q1 > 1: continue
                if gcd(q1, q2) != 1: continue
                tested += 1
                if m_pair_exact(q1, q2, theta) != law_pair(q1, q2, theta): bad += 1
                if m_pair_exact(q1, -q2, theta) != law_pair(q1, -q2, theta): bad += 1
        print(f"  theta={theta}: tested {tested} pairs x 2 signs: {bad} violations")

    print("\n=== 4. corollary checks ===")
    # G+ >= 0, G- <= 0, vanishing conditions
    okpos = all(law_pair(q1, q2, TH) >= TH*TH for q1 in range(1, 30) for q2 in range(q1+1, 30) if gcd(q1,q2)==1)
    okneg = all(law_pair(q1, -q2, TH) <= TH*TH for q1 in range(1, 30) for q2 in range(1, 30) if gcd(q1,q2)==1)
    print(f"  same-sign m >= theta^2 on sweep: {okpos};  mixed-sign m <= theta^2 on sweep: {okneg}")
    print(f"  k=8 Hunter floor: 1 - 7*(1/7) + 6*(1/49) = {1 - 7*F(1,7) + 6*F(1,49)} = 6/49 ~ {float(6*F(1,49)):.4f}")
    print(f"  k=9 CORRECTION (MISTAKE-122): 1 - 8*(1/7) + 7*(1/49) = {1 - 8*F(1,7) + 7*F(1,49)}  (bare floor exactly 0, NOT 1/49)")

    print("\n=== 5. TRIPLE-MASS LAW probe: H := (m123 - theta^3 - theta*sum G_ij/(qq)) * q1q2q3 ===")
    from collections import defaultdict
    tab = defaultdict(set)
    trips = []
    for q1 in range(1, 12):
        for q2 in range(q1+1, 16):
            for q3 in range(q2+1, 20):
                if gcd(q1,q2)!=1 or gcd(q1,q3)!=1 or gcd(q2,q3)!=1: continue
                trips.append((q1,q2,q3))
    rng.shuffle(trips)
    for (q1,q2,q3) in trips[:120]:
        m3 = m_triple_exact(q1, q2, q3, TH)
        pred2 = TH**3 + TH*((law_pair(q1,q2,TH)-TH*TH) + (law_pair(q1,q3,TH)-TH*TH) + (law_pair(q2,q3,TH)-TH*TH))
        H = (m3 - pred2) * q1*q2*q3
        tab[(q1%7, q2%7, q3%7)].add(H)
    const_classes = sum(1 for v in tab.values() if len(v) == 1)
    multi = [(k, sorted(map(float, v))) for k, v in tab.items() if len(v) > 1]
    print(f"  residue classes seen: {len(tab)}; constant: {const_classes}; varying: {len(multi)}")
    for k, v in multi[:6]:
        print(f"    class {k}: H values {v[:4]}{'...' if len(v)>4 else ''}")
    allH = [h for v in tab.values() for h in v]
    print(f"  H range: [{float(min(allH)):.4f}, {float(max(allH)):.4f}]  (law form holds iff classes constant)")

    print("\n=== 6. BONFERRONI-3 ENDPOINT FLOOR at k=8 (exact S2, S3 on the 7 endpoint differences) ===")
    def bonf3_floor(E):
        e_top = max(E)
        D = sorted(e_top - e for e in E if e != e_top)   # 7 positive differences
        S2 = sum(m_pair_exact(D[a], D[b], TH) for a in range(7) for b in range(a+1, 7))
        S3 = F(0)
        for a in range(7):
            for b in range(a+1, 7):
                for c in range(b+1, 7):
                    g12 = gcd(D[a], D[b]); g = gcd(g12, D[c])
                    S3 += m_triple_exact(D[a]//g, D[b]//g, D[c]//g, TH) if g > 1 else m_triple_exact(D[a], D[b], D[c], TH)
        return S2 - S3, S2, S3
    bank = {
        "AP {1..8}": [1,2,3,4,5,6,7,8],
        "spread": [0,5,11,17,26,33,41,50],
        "two-cluster": [0,1,2,3,100,101,102,103],
        "S155-adv": [1,11,12,15,25,28,32,33],
    }
    x = (np.arange(40009)+0.5)/40009
    for nm, E in bank.items():
        fl, S2, S3 = bonf3_floor(E)
        # true W_end numeric
        e_top = max(E); ok = np.ones_like(x, bool)
        for e in E:
            if e == e_top: continue
            v = ((e_top-e)*x) % 1.0
            ok &= ~((v > 0) & (v <= 1/7))
        print(f"  {nm:>14}: S2 = {float(S2):.4f}  S3 = {float(S3):.4f}  Bonf3 floor = {float(fl):.4f}   true W_end = {ok.mean():.4f}")
    # adversarial min of the Bonf3 floor (numeric search, exact confirm)
    def bonf3_numeric(E, xg):
        e_top = max(E)
        D = [e_top - e for e in E if e != e_top]
        hits = np.stack([(((d*xg) % 1.0) > 0) & (((d*xg) % 1.0) <= 1/7) for d in D])
        S1 = hits.mean(axis=1).sum()
        S2 = 0.0; S3 = 0.0
        for a in range(7):
            for b in range(a+1, 7):
                pab = (hits[a] & hits[b])
                S2 += pab.mean()
                for c in range(b+1, 7):
                    S3 += (pab & hits[c]).mean()
        return 1 - S1 + S2 - S3
    xs = (np.arange(6007)+0.5)/6007
    gmin = (2.0, None)
    for trial in range(24):
        Hh = int(rng.choice([9, 14, 22, 40, 70]))
        E = sorted(rng.choice(np.arange(0, Hh+1), size=8, replace=False).tolist())
        cur = bonf3_numeric(E, xs)
        for step in range(35):
            i = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 30, 60, 90]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 8: continue
            c = bonf3_numeric(cand, xs)
            if c < cur - 1e-4: E, cur = cand, c
        v = bonf3_numeric(E, x)
        if v < gmin[0]: gmin = (v, tuple(E))
    print(f"  ADVERSARIAL MIN Bonf3 endpoint floor = {gmin[0]:.4f} at {gmin[1]}")
    fl, S2, S3 = bonf3_floor(list(gmin[1]))
    print(f"    exact confirm: floor = {fl} = {float(fl):.5f}  (R-route bar 0.197; Hunter floor 6/49 = {float(F(6,49)):.4f})")
