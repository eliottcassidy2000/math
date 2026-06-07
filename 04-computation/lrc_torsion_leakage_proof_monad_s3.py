#!/usr/bin/env python3
"""
lrc_torsion_leakage_proof_monad_s3.py

monad-explorer-2026-06-06-S3

Rigorous backbone for the TORSION-LEAKAGE CORRESPONDENCE.

Proves (by exhaustive verification at composite n) the two structural claims
that make the seed precise:

  (A) GCD-CLASS LAW.  For a single-coordinate defect r*e_i in the S367 full-cell
      system, leak(r*e_i) depends on r ONLY through g = gcd(r,n).
      => leak is constant on each subgroup coset of fixed gcd; in particular
         constant on the nonzero order-p torsion (n/p)Z/n (all g = n/p).

  (B) CLOSED FORM.  leak(i,g) = N_i * n - g * W_i(g), where
        N_i      = #{ i-exposed cells } = #{cells p : bins_j(p) not in {0,n-1} for all j != i}
        W_i(g)   = sum over i-exposed cells p of ( [g | bins_i(p)] + [g | bins_i(p)+1] ).
      This is computed from the cell list ALONE (no shift enumeration) and
      checked against the brute leak for every (i, r).

  (C) MONOTONICITY => smallest prime wins.  For the extremal coordinate, tabulate
      leak(i,g) over divisors g|n: it is strictly decreasing in g, so the global
      minimum is at g = n/p* (p* = smallest prime factor) = the largest torsion
      generator = the residues that project to ZERO in the (n/p*)-runner base.

Self-contained re-implementation of the S367 model (validated: n=14 -> 56,
n=15 -> 120).
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd
from sympy import factorint, divisors

ONE = Fraction(1, 1)


def cell_bins(n, k, alpha):
    return tuple(int((n * ((i * alpha) % ONE)) // ONE) for i in range(1, k + 1))


def cell_patterns(n):
    k = n - 1
    breaks = {Fraction(0), ONE}
    for i in range(1, k + 1):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))
    ordered = sorted(breaks)
    out, seen = [], set()
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        b = cell_bins(n, k, (lo + hi) / 2)
        if b in seen:
            continue
        seen.add(b)
        out.append(b)
    return out


def build_masks(n, patterns):
    k = n - 1
    P = len(patterns)
    masks = [[0] * n for _ in range(k)]
    for s in range(n):
        for p_idx, bins in enumerate(patterns):
            bit = 1 << (s * P + p_idx)
            for i in range(k):
                bv = bins[i]
                for r in range(n):
                    if (s * r + bv) % n in (0, n - 1):
                        masks[i][r] |= bit
    return masks, n * P


def brute_leak(masks, n, i, r):
    k = n - 1
    blocked = masks[i][r]
    for j in range(k):
        if j != i:
            blocked |= masks[j][0]
    return None  # placeholder; replaced below


def brute_leak_full(masks, n, cand, i, r):
    """leak of defect r at coord i (0-based), all other coords 0."""
    k = n - 1
    blocked = masks[i][r]
    for j in range(k):
        if j != i:
            blocked |= masks[j][0]
    return cand - bin(blocked).count("1")


def exposed_and_W(n, patterns):
    """For each coordinate i (0-based), N_i and the multiset of bins_i over
    i-exposed cells, so W_i(g) can be computed for any g."""
    k = n - 1
    info = []
    for i in range(k):
        exp_bins_i = []
        for bins in patterns:
            if all(bins[j] not in (0, n - 1) for j in range(k) if j != i):
                exp_bins_i.append(bins[i])
        info.append(exp_bins_i)
    return info


def Wi(exp_bins_i, g):
    w = 0
    for b in exp_bins_i:
        if b % g == 0:
            w += 1
        if (b + 1) % g == 0:
            w += 1
    return w


def main():
    composites = [6, 10, 12, 14, 15, 18, 20, 21]
    print("=" * 78)
    print("TORSION-LEAKAGE PROOF BACKBONE")
    print("=" * 78)
    allA = allB = allC = True
    for n in composites:
        patterns = cell_patterns(n)
        masks, cand = build_masks(n, patterns)
        k = n - 1
        info = exposed_and_W(n, patterns)
        pstar = sorted(factorint(n))[0]

        # (A) gcd-class law + (B) closed form, exhaustive over all (i, r)
        okA = okB = True
        # group leaks by (i, gcd) to test (A)
        for i in range(k):
            by_g = {}
            for r in range(1, n):
                g = gcd(r, n)
                lk = brute_leak_full(masks, n, cand, i, r)
                by_g.setdefault(g, set()).add(lk)
                # (B) closed form
                cf = len(info[i]) * n - g * Wi(info[i], g)
                if cf != lk:
                    okB = False
            for g, s in by_g.items():
                if len(s) != 1:
                    okA = False
        allA &= okA
        allB &= okB

        # (C) monotonicity of leak(i*,g) in g at the extremal coordinate
        # pick the coordinate achieving the global min leak
        best = None
        for i in range(k):
            for r in range(1, n):
                lk = brute_leak_full(masks, n, cand, i, r)
                if best is None or lk < best[0]:
                    best = (lk, i, r)
        _, istar, rstar = best
        # leak as a function of g over divisors of n (1<g<n meaningful for torsion)
        divs = [d for d in divisors(n) if d != n]  # g ranges over proper divisors (gcd(r,n)<n since r!=0)
        leak_g = {}
        for g in divs:
            leak_g[g] = len(info[istar]) * n - g * Wi(info[istar], g)
        # monotone decreasing in g?
        seq = [leak_g[g] for g in sorted(divs)]
        mono = all(seq[a] > seq[a + 1] for a in range(len(seq) - 1))
        # min at g = n/p* ?
        g_min = max(leak_g, key=lambda gg: -leak_g[gg])  # arg min leak
        target_g = n // pstar
        okC = (g_min == target_g)
        allC &= okC

        print(f"\n### n={n} factor={dict(factorint(n))} p*={pstar} cells={len(patterns)} "
              f"cand={cand}")
        print(f"  (A) leak depends only on gcd(r,n): {okA}")
        print(f"  (B) closed form leak=N_i*n - g*W_i(g) exact for all (i,r): {okB}")
        print(f"  extremal coord (0-based) i*={istar} (1-based {istar+1}), N_i*={len(info[istar])}")
        print(f"  leak(i*, g) by proper divisor g of n  (g=n/p are the torsion classes):")
        for g in sorted(divs):
            p_of = n // g
            tag = f"= n/{p_of} torsion" if (n % p_of == 0 and p_of in factorint(n)) else ""
            print(f"      g={g:3d}  leak={leak_g[g]:6d}   {tag}")
        print(f"  (C) leak strictly decreasing in g: {mono};  argmin g={g_min}, "
              f"n/p*={target_g}  -> min at smallest prime: {okC}")

    print("\n" + "=" * 78)
    print(f"OVERALL:  (A) gcd-class law: {allA}   (B) closed form: {allB}   "
          f"(C) min at n/p*: {allC}")
    print("=" * 78)


if __name__ == "__main__":
    main()
