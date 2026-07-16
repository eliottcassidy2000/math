#!/usr/bin/env python3
"""death-star-2026-07-16-S17: the 6617 identity — exact computation of BOTH sides.
flat(F) = lambda{x : S2(x) = 6/7} (clock theorem THM-878 flat set, tight AP {1..13})
corridor G = good_{1/14}({1..12}) (THM-853 corridor value at lam = 1/14).
Claim in canon (kps cont.28 / mac-mini backlog): lambda(F) = 6617/97020 = 2 x 6617/194040
= 2 m(G). This script: exact interval lists for F and G, the measure identity, and the
component-level correspondence (is there a 2:1 map? which one?).
S2(x) = sum_{d=1}^{12} (13-d) * max(0, 1/7 - ||d x||)   [pair overlaps of the 13 arcs]
Slope balance note: term d has slope ±d(13-d), symmetric under d <-> 13-d — the flat
directions pair difference d with 13-d.  Naive half-turn F = G u (G+1/2) is FALSE
(S2 has one-sided kink at 4/7: slopes +10/+16 — checked below); the true map is found here.
"""
from fractions import Fraction as F_
from math import gcd

SEV = F_(1, 7)

def dist(x):
    f = x - x.__floor__()
    return min(f, 1 - f)

def S2(x):
    tot = F_(0)
    for d in range(1, 13):
        dd = dist(d * x)
        if dd < SEV:
            tot += (13 - d) * (SEV - dd)
    return tot

def flat_set():
    """Exact flat set of S2 on [0,1]: kinks at a/(7d) and a/d, d<=12."""
    kinks = set()
    for d in range(1, 13):
        for a in range(0, 7 * d + 1):
            kinks.add(F_(a, 7 * d))
        for a in range(0, d + 1):
            kinks.add(F_(a, d))
    ks = sorted(k for k in kinks if 0 <= k <= 1)
    target = F_(6, 7)
    segs = []
    for i in range(len(ks) - 1):
        a, b = ks[i], ks[i + 1]
        if S2(a) == target and S2(b) == target:
            # piecewise linear between adjacent kinks: flat on [a,b]
            segs.append((a, b))
    # merge adjacent
    merged = []
    for a, b in segs:
        if merged and merged[-1][1] == a:
            merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]

def good_set_14(speeds):
    """Exact closed good set {x: ||vx|| >= 1/14 for v in speeds} on [0,1]."""
    G = [(F_(0), F_(1))]
    for w in sorted(speeds, reverse=True):
        arcs = []
        for c in range(0, w + 1):
            a = F_(14 * c + 1, 14 * w)
            b = F_(14 * c + 13, 14 * w)
            if b < 0 or a > 1:
                continue
            arcs.append((max(a, F_(0)), min(b, F_(1))))
        out = []
        i = j = 0
        while i < len(G) and j < len(arcs):
            a1, a2 = G[i]; b1, b2 = arcs[j]
            lo, hi = max(a1, b1), min(a2, b2)
            if lo <= hi:
                out.append((lo, hi))
            if a2 < b2: i += 1
            else: j += 1
        G = out
    return G

if __name__ == "__main__":
    Fset = flat_set()
    lamF = sum(b - a for a, b in Fset)
    print(f"FLAT SET F of S2 (tight AP): {len(Fset)} components, lambda = {lamF} ≈ {float(lamF):.6f}")
    print(f"  THM-878 value 6617/97020 ≈ {float(F_(6617,97020)):.6f}  match: {lamF == F_(6617,97020)}")
    for a, b in Fset:
        print(f"    [{a}, {b}]  len {b-a}")
    G = good_set_14(range(1, 13))
    mG = sum(b - a for a, b in G)
    print(f"\nCORRIDOR good set G = good_1/14({{1..12}}): {len(G)} components, m = {mG} ≈ {float(mG):.6f}")
    print(f"  corridor value 6617/194040 ≈ {float(F_(6617,194040)):.6f}  match: {mG == F_(6617,194040)}")
    for a, b in G:
        print(f"    [{a}, {b}]  len {b-a}")
    print(f"\nIDENTITY lambda(F) = 2 m(G): {lamF == 2*mG}")
    # correspondence hunt
    print("\nCORRESPONDENCE:")
    # (i) is G a subset of F?
    def inside(x, ivs):
        return any(a <= x <= b for a, b in ivs)
    sub = all(inside(a, Fset) and inside(b, Fset) for a, b in G)
    print(f"  G ⊆ F (endpoints test): {sub}")
    # (ii) half-turn twin
    Gh = sorted(((a + F_(1,2)) % 1, (b + F_(1,2)) % 1) for a, b in G)
    subh = all(inside(a, Fset) and inside(b, Fset) for a, b in Gh)
    print(f"  G + 1/2 ⊆ F: {subh}")
    # (iii) component length multiset comparison: F lengths vs doubled G lengths vs G u G+1/2
    lF = sorted(b - a for a, b in Fset)
    lG = sorted(b - a for a, b in G)
    print(f"  F lengths:      {[str(x) for x in lF]}")
    print(f"  G lengths:      {[str(x) for x in lG]}")
    print(f"  2x G lengths:   {[str(2*x) for x in lG]}")
    # (iv) is F = union of G and reflected/shifted copies? test x -> x+1/2, x -> 1-x, x -> 1/2-x
    def norm(ivs):
        return sorted((a, b) for a, b in ivs)
    cands = {
        "G ∪ (G+1/2)": None,
        "G ∪ (1/2−G)": None,
        "G ∪ (1−G)": None,
    }
    def shift_half(ivs): return norm(((a+F_(1,2))%1, (b+F_(1,2))%1) for a,b in ivs)
    def refl_half(ivs):  return norm(sorted(((F_(1,2)-b)%1, (F_(1,2)-a)%1)) for a,b in ivs)
    def refl_one(ivs):   return norm(sorted((1-b, 1-a)) for a,b in ivs)
    def union(A, B):
        pts = sorted(set(A) | set(B))
        # naive union merge
        allints = sorted(list(A)+list(B))
        out=[]
        for a,b in allints:
            if out and a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a,b])
        return [(a,b) for a,b in out]
    tests = {
        "G ∪ (G+1/2)": union(norm(G), shift_half(G)),
        "G ∪ (1/2−G)": union(norm(G), refl_half(G)),
        "G ∪ (1−G)": union(norm(G), refl_one(G)),
    }
    for name, U in tests.items():
        print(f"  F == {name}: {norm(Fset) == norm(U)}   (meas {sum(b-a for a,b in U)})")
    # (v) THE DOUBLING MAP: F == 2·G (mod 1)?
    G_pos = [(a, b) for a, b in G if b > a]          # positive-length components
    dbl = []
    for a, b in G_pos:
        a2, b2 = 2 * a, 2 * b
        # each component is short (< 1/2), so image is one interval mod 1
        sh = a2.__floor__()
        dbl.append((a2 - sh, b2 - sh))
    dbl = norm(dbl)
    print(f"\n  THE DOUBLING TEST: F == 2·G mod 1 (component-wise): {norm(Fset) == dbl}")
    print(f"  #G positive components: {len(G_pos)} -> {len(dbl)} images; injective (no overlaps): "
          f"{all(dbl[i][1] < dbl[i+1][0] for i in range(len(dbl)-1))}")
    # no half-pair check: G contains no {y, y+1/2} pair
    halfpair = False
    for a, b in G_pos:
        for c, d in G_pos:
            lo = max(a + F_(1,2), c); hi = min(b + F_(1,2), d)
            if lo < hi:
                halfpair = True
    print(f"  G ∩ (G+1/2) has positive measure: {halfpair}  (must be False for Jacobian-2)")
    # isolated points map too?
    iso_G = [a for a, b in G if a == b]
    print(f"  G isolated points: {[str(x) for x in iso_G]} -> doubles "
          f"{[str((2*x) - (2*x).__floor__()) for x in iso_G]} (isolated flat points of S2: "
          f"{[str(S2((2*x) - (2*x).__floor__())) for x in iso_G]})")
