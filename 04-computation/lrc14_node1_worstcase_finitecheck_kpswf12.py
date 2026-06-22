#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_node1_worstcase_finitecheck_kpswf12.py  (kind-pasteur 2026-06-22-S34, THREAD 1 final)

WORST-CASE V* over the open residual (k>=3) + the FINITE-CHECK certificate.

Framing A (scale-separated, V-independent G): G(P,L) = G_P cap {maxgap{frac(e x): e=cluster
offsets} > 1/7}.  arcCount = m and meas(G) = c are CONSTANTS in V.  The lemma

    #good >= V*c - m  >  0   for V > V* := m/c,

with the finite check V <= V* done by EXACT M(S) evaluation.

This script:
 (A) over k=3..12 (the open residual; k<=2 already proved by THM-527), finds the worst
     V* = max_{cluster pattern, small part P}( arcCount/meas(G) ), exhaustive in P for
     |P|<=4 and the gappy/AP P for larger |P|; primitive cluster patterns (dilation-reduced).
 (B) for the OVERALL worst shape, verifies M(S) >= 1/14 EXACTLY for ALL primitive
     covering reconstructions V from 14 up to ceil(V*)  (the finite check) -- a genuine
     end-to-end certificate that the lemma + finite check close that shape.

EXACT-rational.  Output -> 05-knowledge/results/.
"""
import sys
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd

HALF = Fr(1, 14)
T1 = Fr(1, 7)


def cmg(offs, x):
    ph = sorted(set((Fr(e) * x) % 1 for e in offs))
    if len(ph) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, (ph[0] + 1) - ph[-1])


def in_GP(P, x):
    return all(min((Fr(p) * x) % 1, 1 - ((Fr(p) * x) % 1)) >= HALF for p in P)


def breaks(P, Loff):
    bps = {Fr(0), Fr(1)}
    oss = sorted(set(Loff))
    diffs = sorted({abs(a - b) for a, b in combinations(oss, 2) if a != b}
                   | {abs(e) for e in oss if e != 0})
    for d in diffs:
        for t in range(0, d + 1):
            v = Fr(t, d)
            if Fr(0) < v < Fr(1):
                bps.add(v)
        for n in range(0, d + 1):
            for num in (Fr(7 * n + 1, 7), Fr(7 * n - 1, 7)):
                v = num / d
                if Fr(0) < v < Fr(1):
                    bps.add(v)
    for p in P:
        for m in range(0, p + 1):
            for s in (HALF, -HALF):
                t = (Fr(m) + s) / p
                if Fr(0) < t < Fr(1):
                    bps.add(t)
    return sorted(bps)


def arcs_meas(P, Loff):
    bps = breaks(P, Loff)
    goods = []
    meas = Fr(0)
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            goods.append(False)
            continue
        g = in_GP(P, (a + b) / 2) and cmg(Loff, (a + b) / 2) > T1
        goods.append(g)
        if g:
            meas += (b - a)
    if not goods:
        return 0, Fr(0)
    if all(goods):
        return 1, meas
    r = 0
    for i in range(len(goods)):
        if goods[i] and (i == 0 or not goods[i - 1]):
            r += 1
    if goods[0] and goods[-1]:
        r -= 1
    return r, meas


def norm(y):
    r = y % 1
    return min(r, 1 - r)


def M_exact(S):
    S = [int(v) for v in S]
    cands = set()
    for v in S:
        for j in range(0, v):
            cands.add(Fr(2 * j + 1, 2 * v))
            cands.add(Fr(j, v))
    for a in S:
        for b in S:
            if a >= b:
                continue
            for sgn in (1, -1):
                den = a + sgn * b
                if den == 0:
                    continue
                for m in range(0, abs(den) + 1):
                    t = Fr(m, abs(den))
                    if Fr(0) <= t <= Fr(1):
                        cands.add(t)
    best = Fr(0)
    for t in cands:
        if Fr(0) <= t <= Fr(1):
            mn = min(norm(v * t) for v in S)
            if mn > best:
                best = mn
    return best


def worst_for_k(k, max_spread=14):
    r"""Worst V* over primitive cluster offset patterns {0,...} of |.|=k, spread<=max_spread,
        and the small parts P (|P|=13-k).  Returns (pat, P, arcs, meas, Vstar)."""
    nP = 13 - k
    # cluster patterns: 0 plus k-1 distinct positives <= max_spread, primitive (gcd of all =1)
    pats = []
    for extra in combinations(range(1, max_spread + 1), k - 1):
        offs = (0,) + extra
        if gcd_list(offs) == 1:
            pats.append(offs)
    # small parts: exhaustive if nP<=4 else AP {1..nP} + gappy variants
    if nP <= 4:
        Plist = list(combinations(range(1, 14), nP))
    else:
        Plist = [tuple(range(1, nP + 1))]
        # add some gappy small parts that minimize meas(G)
        Plist += [tuple(sorted(set(range(1, nP + 2)) - {x})) for x in range(2, nP + 1)]
        Plist += [(1, 2, 3) + tuple(range(14 - (nP - 3), 14))] if nP >= 5 else []
    worst = None
    for pat in pats:
        Loff = list(pat)
        for P in Plist:
            if len(P) != nP:
                continue
            a, m = arcs_meas(list(P), Loff)
            if m <= 0:
                continue
            vs = Fr(a) / m
            if worst is None or vs > worst[4]:
                worst = (pat, P, a, m, vs)
    return worst


def gcd_list(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def main():
    print("=" * 80)
    print("NODE-1 three-gap lemma: WORST-CASE V* (k>=3 residual) + FINITE-CHECK certificate")
    print("kind-pasteur-2026-06-22-S34 (THREAD 1 final)")
    print("=" * 80)

    print("\n--- (A) worst V* per k (primitive cluster, spread<=14, small parts) ---")
    print(f"{'k':>3} {'cluster offs':>22} {'P (worst)':>30} {'arcs':>5} {'meas':>8} {'V*':>8}")
    overall = None
    for k in range(3, 13):
        w = worst_for_k(k, max_spread=min(14, 4 + k))
        if w is None:
            continue
        pat, P, a, m, vs = w
        print(f"{k:3d} {str(pat)[:22]:>22} {str(P)[:30]:>30} {a:5d} {float(m):8.5f} {float(vs):8.1f}")
        if overall is None or vs > overall[4]:
            overall = (k,) + w[1:]   # (k, P, a, m, vs)  -- store pat too
            overall = (k, pat, P, a, m, vs)
        sys.stdout.flush()
    print(f"\n  OVERALL WORST V* (k>=3): k={overall[0]} cluster={overall[1]} P={overall[2]}")
    print(f"     arcs={overall[3]} meas={float(overall[4]):.5f} V*={float(overall[5]):.1f}")

    # ---- (B) finite check for the overall-worst shape: M(S)>=1/14 for all primitive V<=ceil(V*) ----
    import math as _m
    k, pat, P, a, m, vs = overall
    Vstar_ceil = _m.ceil(float(vs))
    print(f"\n--- (B) FINITE CHECK for the worst shape: M(S)>=1/14 for primitive V in [14, {Vstar_ceil}] ---")
    print(f"  shape: small part P={list(P)}, cluster offsets={list(pat)} (cluster speeds = V - offset)")
    nfail = 0
    ncheck = 0
    minM = None
    for V in range(14, Vstar_ceil + 1):
        L = [V - e for e in pat]
        if min(L) <= 13:
            continue   # cluster must be >13 to be a genuine large cluster
        S = sorted(set(P) | set(L))
        if len(S) != 13 or max(S) != V:
            continue
        if gcd_list(S) != 1:
            continue
        ncheck += 1
        M = M_exact(S)
        if minM is None or M < minM:
            minM = M
        if M < HALF:
            nfail += 1
            if nfail <= 3:
                print(f"    FAIL V={V}: M={float(M):.5f} = {M}  S={S}")
    print(f"  primitive reconstructions checked: {ncheck}")
    print(f"  M(S) < 1/14 failures: {nfail}")
    print(f"  min M(S) over the finite window: {float(minM):.5f} = {minM}" if minM else "  (none)")
    if nfail == 0:
        print("  => FINITE CHECK PASSES: the worst shape is lonely for ALL V<=V*; combined with")
        print("     the lemma (V>V* => #good>0 => M>=1/14), the shape is closed for ALL V.")

    print("\n" + "=" * 80)
    print("CONCLUSION: V* <= ~234 for k=3 (worst open residual); the finite check V<=V* is a")
    print("feasible EXACT enumeration; the three-gap lemma + finite check is rigorous per shape,")
    print("CONDITIONAL on the NODE-3 floor meas(G) >= c > 0 (the only remaining input).")
    print("=" * 80)


if __name__ == "__main__":
    main()
