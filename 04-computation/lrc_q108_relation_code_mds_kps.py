#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_relation_code_mds_kps.py   (kind-pasteur 2026-06-21, HYP-2723)

The MDS / arc / coding-theory lens on OPEN-Q-108 (LRC(14) wide-cover crux).

Frontier (mac-mini HYP-2719): the carrier error
    corr(E) = measS7(E) - iid_k = Sum_{0 != n in Lambda(E)} K(n),
    Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }  (the OFFSET RELATION LATTICE
                = LRC twin of the cycle space; support-size seam).
This script frames Lambda(E) as a LINEAR CODE [k, k-1, d], d = MINIMUM SUPPORT
of a nonzero relation, and tests:
  (a) corr anti-tracks d: AP/consec MINIMIZE d (anti-MDS, hardest); Sidon/arc
      MAXIMIZE d (MDS / general position, easiest, corr ~ 0).
  (b) the weight enumerator A_s(E) = #(support-s primitive relations, |coef|<=B)
      -- A_3 (3-APs + Schur triples = additive energy) is the leading binding term.
  (c) the '56' / tournament bijection hint.

EXACT arithmetic for measS7; integer search for the relation code.
"""
import itertools
from fractions import Fraction as Fr
from math import comb, gcd
from functools import reduce

# ---------------------------------------------------------------- measS7 (exact)
def stirling2(n, k):
    # S(n,k)
    if k == 0:
        return 1 if n == 0 else 0
    S = [[0]*(k+1) for _ in range(n+1)]
    S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, min(i, k)+1):
            S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][k]

def iid_k(k, p=7):
    # P(all p sectors hit) for k independent uniform phases = p! S(k,p)/p^k
    from math import factorial
    return Fr(factorial(p) * stirling2(k, p), p**k)

def measS7(E, p=7):
    """Exact measure of { x in [0,1) : {frac(e x): e in E} hits all p sectors }.
       E: list of nonnegative ints. Sector(y)=floor(p*frac(y))."""
    E = [int(e) for e in E]
    bps = set([Fr(0), Fr(1)])
    for e in E:
        if e == 0:
            continue
        for t in range(0, p*e):            # x = t/(p e), boundaries of frac(e x)
            bps.add(Fr(t, p*e))
    pts = sorted(bps)
    total = Fr(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        sectors = set()
        for e in E:
            y = (e*mid) % 1                # Fraction mod 1
            s = int(p*y)                   # floor(p*frac)
            sectors.add(s)
            if len(sectors) == p:
                break
        if len(sectors) == p:
            total += (b - a)
    return total

def corr(E, p=7):
    return measS7(E, p) - iid_k(len(E), p)

# ---------------------------------------------------------- relation code Lambda(E)
def support_spectrum(E, B=2, max_support=5):
    """Count PRIMITIVE nonzero relations sum n_i e_i = 0, |n_i|<=B, by support.
       Returns dict support -> count, plus min support (= min distance d)."""
    E = [int(e) for e in E]
    k = len(E)
    counts = {s: 0 for s in range(2, max_support+1)}
    seen = set()
    idxs = range(k)
    for s in range(2, max_support+1):
        for combo in itertools.combinations(idxs, s):
            # coeffs in [-B,B]\{0} on these s positions, others 0
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs):
                    continue
                if sum(c*E[i] for c, i in zip(coefs, combo)) != 0:
                    continue
                g = reduce(gcd, [abs(c) for c in coefs])
                prim = tuple(c//g for c in coefs)
                # canonical sign (first nonzero positive)
                if prim[0] < 0:
                    prim = tuple(-c for c in prim)
                key = (combo, prim)
                if key in seen:
                    continue
                # also skip if this relation actually has smaller support (some prim coef 0 -- impossible here)
                seen.add(key)
                counts[s] += 1
    nz = [s for s in counts if counts[s] > 0]
    dmin = min(nz) if nz else None
    return counts, dmin

# --------------------------------------------------------------------- batteries
def is_sidon(E):
    sums = {}
    E = list(E)
    for i in range(len(E)):
        for j in range(i, len(E)):
            s = E[i]+E[j]
            if s in sums:
                return False
            sums[s] = 1
    return True

def main():
    p = 7
    print("="*78)
    print("HYP-2723: the relation code Lambda(E) and the LRC carrier error corr(E)")
    print("="*78)
    print(f"iid_k (k=8,9,10) = {[float(iid_k(k)) for k in (8,9,10)]}")
    print()

    # battery of k=8 sets (offsets; 0 included as the pinned observer-phase)
    battery = {
        "consec/AP {0..7}":      [0,1,2,3,4,5,6,7],
        "AP step2 {0,2..14}":    [0,2,4,6,8,10,12,14],
        "dyadic {0,1,2,4,8,..}": [0,1,2,4,8,16,32,64],
        "Sidon/arc (Mian-Chowla)":[0,1,3,7,12,20,30,44],
        "wide near-consec":      [0,1,2,3,4,5,6,40],
        "two-block":             [0,1,2,3,40,41,42,43],
        "random wide":           [0,5,9,14,22,33,41,50],
    }
    print(f"{'set':<26}{'measS7':>10}{'corr':>10}{'dmin':>6}{'A2':>5}{'A3':>5}{'A4':>5}{'Sidon?':>8}")
    rows = []
    for name, E in battery.items():
        c = corr(E, p)
        m = measS7(E, p)
        counts, dmin = support_spectrum(E, B=2, max_support=4)
        sid = is_sidon([e for e in E if e != 0])
        rows.append((name, E, float(m), float(c), dmin, counts, sid))
        print(f"{name:<26}{float(m):>10.4f}{float(c):>10.4f}{str(dmin):>6}"
              f"{counts.get(2,0):>5}{counts.get(3,0):>5}{counts.get(4,0):>5}{str(sid):>8}")

    print()
    print("READING: hardness ~ corr (want corr small for the wide bound). HYP-2723 predicts")
    print("  AP/consec: small dmin, large A3 (additive energy), LARGEST corr (hardest).")
    print("  Sidon/arc: large dmin, A3=0, corr ~ 0 (easiest).")
    # quick correlation A3 vs corr
    import statistics
    A3 = [r[5].get(3,0) for r in rows]
    cr = [r[3] for r in rows]
    if len(set(A3)) > 1:
        # Pearson
        n=len(A3); mA=sum(A3)/n; mC=sum(cr)/n
        num=sum((a-mA)*(c-mC) for a,c in zip(A3,cr))
        den=(sum((a-mA)**2 for a in A3)*sum((c-mC)**2 for c in cr))**0.5
        print(f"\n  Pearson(A3 support-3 count, corr) = {num/den:+.3f}  (predict positive: more 3-relations => larger error)")

    # ---- the '56' probe: support-3 relation hypergraph types at small k ----
    print()
    print("="*78)
    print("'56' PROBE: count distinct support-3 relation-HYPERGRAPHS over k-sets")
    print("(the user's hint: 3 -> 56 challenger shapes = A000568(6)=56 tournaments on 6)")
    print("="*78)
    A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456}
    for k in range(4, 9):
        # enumerate the *types* of the support-3 relation hypergraph among consec-like
        # k-sets in a small window; record how many distinct (as a first coarse probe)
        # Here: over all k-subsets of {0..k+2}, the multiset of support-3 relation patterns.
        from collections import Counter
        types = set()
        base = range(0, k+3)
        cnt = 0
        for E in itertools.combinations(base, k):
            counts, dmin = support_spectrum(list(E), B=2, max_support=3)
            # signature: support-3 count + dmin (coarse 'shape')
            types.add((counts.get(2,0), counts.get(3,0), dmin))
            cnt += 1
        a568 = A000568.get(k, '?')
        print(f"  k={k}: {cnt} k-subsets of [0,{k+2}] -> {len(types)} distinct (A2,A3,dmin) shapes  "
              f"[A000568({k})={a568}; A000568(6)=56]")
    print("\nDONE.")

if __name__ == "__main__":
    main()
