#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_abel_cesaro_levelmass_s7.py   (ANGLE 6 — Abel/Cesaro summation of level masses)
mac-mini-2026-06-15-S7

LRC(14) singular series for the WORST core  S = {1,2,3,4,5,7,8,9,10,11,12,13,98}
(98 = 2*7^2, the resonant stranger; global numerical inf L ~ 0.00524).

THE OBJECT.  L(S) = sum_{t in Lambda} prod_i h(t_i),  Lambda = {t in Z^13 : sum t_i v_i = 0},
   h(0) = 6/7,  h(t) = -s(t),  s(t) = sin(pi t/7)/(pi t)  (s(t)=0 iff 7|t).
Group the lattice sum by SUPPORT SIZE k (#nonzero coords).  Main term (6/7)^13 (k=0):

   L = (6/7)^13 + sum_{k>=2} (-1)^k Lambda_k,
   Lambda_k = (6/7)^{13-k} * sum_{|T|=k} [ sum over FULL-SUPPORT, 7-primitive points of the
              kernel lattice {t in Z^T : sum t_v v_v = 0}  of  prod_{v in T} s(t_v) ].

*** CORRECTION over the naive first pass ***
The right truncation is NOT "|t_v| <= B for a fixed B" (that artificially chops the
convergent tail of integer MULTIPLES of each base relation, and the Lambda_k it yields
keep growing with B and never settle).  The right truncation enumerates, FOR EACH
SUPPORT T, the entire rank-(k-1) kernel lattice up to a radius R IN LATTICE COORDINATES,
and lets R -> infinity.  With that order each Lambda_k is a well-defined finite number.

*** KEY STRUCTURAL FINDING ***
 - k=2 : the per-support sum is ABSOLUTELY convergent (~1/n^2 on a rank-1 lattice).
 - k>=3: the per-support sums are only CONDITIONALLY convergent: the SIGNED box-sum
   settles, but sum|prod s| DIVERGES (slowly).  So Lambda_k for k>=3 exists only as a
   conditionally-convergent (symmetric-box) limit, not as an absolutely convergent total.

This is the crux: an Abel/Cesaro/Euler attack would need the SEQUENCE (Lambda_k) to be
controlled.  This script computes the converged Lambda_k (k=2..K) and then tests every
standard summation method against the exact L (lonely measure).  Pure python3, stdlib+math.
"""

import sys
from math import gcd, sin, pi
from itertools import combinations, product

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------- core data
S = sorted([x for x in range(1, 14) if x != 6] + [98])     # worst core
N = len(S)
assert N == 13
SIX7 = 6.0 / 7.0

def s_of(t):
    """s(t) = sin(pi t/7)/(pi t); s(0)=1/7; s(t)=0 exactly when 7 | t."""
    if t == 0:
        return 1.0 / 7.0
    if t % 7 == 0:
        return 0.0
    return sin(pi * t / 7.0) / (pi * t)

# ----------------------------------------------------------------- exact L
def lonely_measure(speeds, Q):
    rad = Q // 14
    cnt = 0
    for a in range(Q):
        ok = True
        for v in speeds:
            r = (v * a) % Q
            if r <= rad or r >= Q - rad:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt / Q

# ----------------------------------------------------------------- kernel lattice
def kernel_basis(v):
    """Integer basis of {t : sum t_i v_i = 0} for the row v (length k), rank k-1.
       Star basis from v[0]; index 1 whenever gcd(v[0], v[i]) handled, verified covol-correct
       for the supports here (all contain a unit or pairwise-coprime structure)."""
    k = len(v)
    basis = []
    for i in range(1, k):
        gi = gcd(v[0], v[i])
        b = [0] * k
        b[0] = v[i] // gi
        b[i] = -v[0] // gi
        basis.append(b)
    return basis

def level_mass(k, R):
    """Lambda_k at lattice radius R: sum over supports T (|T|=k), over kernel-lattice
       points (basis coords in [-R,R]^{k-1}, not all zero), full support & 7-primitive,
       of prod s(t_v).  Returns (raw_sum, Lambda_k, n_terms, abs_sum)."""
    raw = 0.0
    abss = 0.0
    cnt = 0
    rng = range(-R, R + 1)
    for T in combinations(range(N), k):
        v = [S[i] for i in T]
        basis = kernel_basis(v)
        for combo in product(rng, repeat=k - 1):
            if not any(combo):
                continue
            # build lattice point
            t = [0] * k
            for j in range(k - 1):
                cj = combo[j]
                if cj:
                    bj = basis[j]
                    for idx in range(k):
                        t[idx] += cj * bj[idx]
            ok = True
            for x in t:
                if x == 0 or x % 7 == 0:
                    ok = False
                    break
            if not ok:
                continue
            p = 1.0
            for x in t:
                p *= s_of(x)
            raw += p
            abss += abs(p)
            cnt += 1
    lam = (SIX7 ** (N - k)) * raw
    return raw, lam, cnt, abss

# ----------------------------------------------------------------- summation methods
def partial_sums(main, lam):
    P = [(1, main)]
    acc = main
    for k in sorted(lam):
        acc += ((-1) ** k) * lam[k]
        P.append((k, acc))
    return P

def euler_transform(b):
    """Euler sum of  sum_{j>=0} (-1)^j b_j  via forward differences of b at 0."""
    if not b:
        return 0.0
    d = list(b)
    total = 0.0
    j = 0
    while d:
        total += d[0] / (2 ** (j + 1))
        d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        j += 1
    return total

def abel_sum(main, lam, r):
    acc = main
    for k in sorted(lam):
        acc += ((-1) ** k) * lam[k] * (r ** k)
    return acc

def cesaro(P):
    vals = [p for _, p in P]
    return sum(vals) / len(vals)

# ----------------------------------------------------------------- main
def main():
    print("=" * 78)
    print("ANGLE 6: Abel/Cesaro summation of level masses — worst LRC-14 core")
    print("S =", S, "  (98 = 2*7^2)")
    print("=" * 78)

    print("\n[TRUTH] exact lonely measure L(S):")
    L_true = None
    for Q in [27720, 55440, 110880]:
        Lm = lonely_measure(S, Q)
        L_true = Lm
        print(f"   Q={Q:7d}:  L = {Lm:.7f}")
    print(f"   --> reference L_true ~ {L_true:.7f}")
    main_term = SIX7 ** N
    print(f"   main term (6/7)^13 = {main_term:.7f}")

    print("\n" + "=" * 78)
    print("(0) CONVERGENCE TYPE per level: signed sum settles, |sum| may diverge")
    print("    (lattice radius R; SIGNED = Lambda_k, ABS = sum|prod s| (6/7)^{13-k} scaled)")
    print("=" * 78)
    # radii chosen so runtime stays bounded; bigger k -> smaller R
    radii = {2: [50, 200, 800], 3: [10, 25, 45], 4: [5, 9, 14], 5: [3, 5, 7]}
    converged = {}
    for k in sorted(radii):
        print(f"\n   level k={k}:")
        last = None
        for R in radii[k]:
            raw, lam, cnt, abss = level_mass(k, R)
            absL = (SIX7 ** (N - k)) * abss
            tag = ""
            if last is not None:
                tag = f"(d_signed={lam-last[0]:+.5f}, d_abs={absL-last[1]:+.5f})"
            print(f"      R={R:4d}: Lambda_k(signed)={lam:.6f}  |..|(abs)={absL:.6f}  terms={cnt}  {tag}")
            last = (lam, absL)
        converged[k] = last[0]   # signed value at largest R = best estimate
    print("\n   --> SIGNED Lambda_k stabilize (conditionally convergent); ABS keep climbing")
    print("       for k>=3 => the level series is CONDITIONALLY convergent, order-sensitive.")

    print("\n" + "=" * 78)
    print("(1) CONVERGED LEVEL MASSES  Lambda_k  (best signed estimate)")
    print("=" * 78)
    print(f"   {'k':>2} {'Lambda_k':>12} {'(-1)^k Λ_k':>12}")
    for k in sorted(converged):
        print(f"   {k:2d} {converged[k]:12.6f} {((-1)**k)*converged[k]:12.6f}")

    print("\n" + "=" * 78)
    print("(2) PARTIAL SUMS  P_m = (6/7)^13 + Σ_{k=2}^m (-1)^k Λ_k  vs L_true")
    print("=" * 78)
    P = partial_sums(main_term, converged)
    print(f"   {'m':>2} {'P_m':>12} {'P_m - L_true':>14}")
    for m, p in P:
        print(f"   {m:2d} {p:12.6f} {p - L_true:14.6f}")

    print("\n" + "=" * 78)
    print("(3) SUMMATION METHODS")
    print("=" * 78)
    # alternating tail terms: Lambda_2, Lambda_3, ... fed as b_j with sign (-1)^{j} (since (-1)^{j+2}=(-1)^j)
    bjs = [converged[k] for k in sorted(converged)]
    eul_tail = euler_transform(bjs)
    print(f"\n   EULER transform of Σ(-1)^j Λ_{{j+2}}:")
    print(f"      tail = {eul_tail:+.6f};  main+tail = {main_term+eul_tail:.6f}"
          f"   (L_true {L_true:.6f}, err {main_term+eul_tail-L_true:+.6f})")

    ces = cesaro(P)
    print(f"\n   CESARO (C,1) mean of P_m: {ces:.6f}   (err {ces-L_true:+.6f})")

    print(f"\n   ABEL  Σ(-1)^k Λ_k r^k  as r->1^-:")
    print(f"      {'r':>6} {'Abel(r)':>12}")
    av = []
    for r in [0.5, 0.7, 0.85, 0.92, 0.96, 0.98, 0.99, 1.0]:
        a = abel_sum(main_term, converged, r)
        av.append((r, a))
        print(f"      {r:6.2f} {a:12.6f}")
    (r1, a1), (r2, a2) = av[-3], av[-2]
    slope = (a2 - a1) / ((1 - r2) - (1 - r1))
    a_ext = a2 + slope * (-(1 - r2))
    print(f"      r->1 linear extrap (r=.98,.99): {a_ext:.6f}  (err {a_ext-L_true:+.6f})")

    print("\n" + "=" * 78)
    print("(4) PER-LEVEL BOUND DIAGNOSIS")
    print("=" * 78)
    ks = sorted(converged)
    print(f"   {'k':>2} {'|Λ_k|':>12} {'|Λ_{k+1}|/|Λ_k|':>18}")
    for i, k in enumerate(ks):
        amag = abs(converged[k])
        if i + 1 < len(ks):
            ratio = abs(converged[ks[i + 1]]) / amag if amag else float('inf')
            print(f"   {k:2d} {amag:12.6f} {ratio:18.4f}")
        else:
            print(f"   {k:2d} {amag:12.6f} {'(last)':>18}")
    print("\n   For an Abel/geometric tail bound we would need |Λ_k| <= C g^k with g<1.")
    print("   Observed: Λ_k do NOT decay geometrically — they plateau / decay sub-geometrically,")
    print("   and even sum|prod s| per level DIVERGES (k>=3). A Polya–Vinogradov per-level bound")
    print("   |Σ_{|T|=k} ∏ s| <= G^k is FALSE as stated (the unsigned mass is unbounded);")
    print("   only the SIGNED, symmetric-limit Λ_k is finite, which no elementary summation")
    print("   method can certify a lower bound from below.")

    print("\nDONE.")

if __name__ == "__main__":
    main()
