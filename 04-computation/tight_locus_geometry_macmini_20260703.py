#!/usr/bin/env python3
"""
GEOMETRY OF THE TIGHT LOCUS (mac-mini-2026-07-03-S31) -- the open core (rigidity M=1/14 => {AP,GW}).
Building on THM-610 Lemma 2 (tight => 14|q*, runners on dilated 14th-root configs). Explore:
 (1) the known tight families AP={1..13}, GW={1..11,13,24}: verify M=1/14, q*, #tight points, binding pairs,
     residues mod 14, do they cover the +-units {1,3,5,9,11,13}?
 (2) CONJECTURE 'primitive tight => q*=14': search for ANY primitive tight family (M=1/14) with q*>14.
 (3) the danger-arc COVERING reformulation: N(t)=#{i:||v_i t||<1/14}; tight <=> N>=1 a.e. (safe set = points).
     integral N = 13/7; measure where N=0 (the loose/safe measure L). Check L(AP)=L(GW)=0.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0); pts = []
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            t = F(a, q); m = min(nd(v*t) for v in sp)
            if m > best: best = m; pts = [t]
            elif m == best: pts.append(t)
    return best, sorted(set(pts))
def tight_points(sp, M):
    """all t=a/q (q over breakpoints) achieving min-dist == M."""
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    pts = []
    for q in sorted(Q):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            t = F(a, q)
            if min(nd(v*t) for v in sp) == M: pts.append(t)
    return sorted(set(pts))
def binding_pairs(sp, t, M):
    """runners at distance exactly M, and which are at +M vs -M (mod 1)."""
    plus = [v for v in sp if (v*t) % 1 == M]          # phase = M (just above 0)
    minus = [v for v in sp if (v*t) % 1 == 1-M]        # phase = 1-M (just below 1)
    return plus, minus

units14 = {1,3,5,9,11,13}
if __name__ == "__main__":
    print("(1) KNOWN TIGHT FAMILIES: AP={1..13}, GW={1..11,13,24}")
    print("="*82)
    for name, S in [("AP", list(range(1,14))), ("GW", [1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        M, _ = M_exact(S)
        tp = tight_points(S, M)
        qs = sorted(set(t.denominator for t in tp))
        res = sorted(set(v % 14 for v in S))
        covu = units14.issubset(set(res))
        # binding pair at the first tight point
        t0 = tp[0]; plus, minus = binding_pairs(S, t0, M)
        print(f"  {name:>3} {S}")
        print(f"      M={M}={float(M):.5f}  #tight pts={len(tp)}  q* denominators={qs}  primitive={gcd_all(S)==1}")
        print(f"      residues mod14={res}  covers +-units{{1,3,5,9,11,13}}? {covu}  (missing: {sorted(set(range(1,14))-set(res))})")
        print(f"      at t*={t0}: binding +M runners={plus}, -M runners={minus}, sum pairs sum to q*={t0.denominator}?")
        for p in plus:
            for mns in minus:
                if (p+mns) % t0.denominator == 0:
                    print(f"         pair ({p},{mns}): sum={p+mns} = {(p+mns)//t0.denominator} x q*({t0.denominator})")

    print("\n(2) CONJECTURE 'primitive tight => q*=14': search perturbations of AP for M=1/14 with q*>14")
    print("="*82)
    target = F(1,14); found_gt14 = []; found_tight = 0
    base = list(range(1,14))
    # perturb: replace 1-2 speeds with a value congruent mod 14 (keeps residues) or nearby, small
    cnt = 0
    for drop in range(13):
        for add in range(15, 200):
            S = base[:drop] + base[drop+1:] + [add]
            if len(set(S)) != 13: continue
            if gcd_all(S) != 1: continue
            cnt += 1
            M, _ = M_exact(S)
            if M == target:
                found_tight += 1
                tp = tight_points(S, M); qs = sorted(set(t.denominator for t in tp))
                if any(q > 14 for q in qs):
                    found_gt14.append((S, qs))
    # also two-speed perturbations congruent mod 14 (the GW mechanism: replace k by k+14)
    for drop1, drop2 in combinations(range(13), 2):
        for a1 in [base[drop1]+14, base[drop1]+28]:
            for a2 in [base[drop2]+14, base[drop2]+28]:
                S = [base[i] for i in range(13) if i not in (drop1,drop2)] + [a1, a2]
                if len(set(S)) != 13 or gcd_all(S) != 1: continue
                M, _ = M_exact(S)
                if M == target:
                    found_tight += 1
                    tp = tight_points(S, M); qs = sorted(set(t.denominator for t in tp))
                    if any(q > 14 for q in qs): found_gt14.append((S, qs))
    print(f"  tested ~{cnt}+ perturbations; found {found_tight} tight (M=1/14) families")
    print(f"  primitive tight families with q*>14: {len(found_gt14)}")
    for S, qs in found_gt14[:8]:
        print(f"     {sorted(S)}  q*={qs}")
    if not found_gt14:
        print("  => NO primitive tight family with q*>14 found (supports 'primitive tight => q*=14')")

    print("\n(3) COVERING reformulation: L(S)=safe measure (where N(t)=0). tight <=> L=0 (safe set = points).")
    for name, S in [("AP", list(range(1,14))), ("GW", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
                    ("deep well", list(range(1,13))+[182]), ("loose {2..14}", list(range(2,15)))]:
        # L = measure of {t: all ||v_i t||>=1/14}. Compute via fine exact scan over lcm-ish grid.
        D = 14 * reduce(lambda a,b: a*b//gcd(a,b), S[:6], 1)  # coarse but exact-ish denom
        D = min(D, 200000)
        safe = sum(1 for a in range(D) if all(nd(v*F(a,D))>=F(1,14) for v in S))
        M,_ = M_exact(S)
        print(f"  {name:>12}: M={float(M):.5f}  safe-fraction≈{safe/D:.5f}  (tight<=>0; loose>0)")
    print("\n=> tight locus = L=0 boundary; Lemma 2 (14|q*) + binding pair + primitive=>q*=14 (if (2) holds)")
    print("   would confine primitive tight families to denom 14 => a FINITE mod-14 residue problem.")
