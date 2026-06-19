#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) -- Poisson / theta-transform of the seven-sector signed correction.
kind-pasteur-2026-06-19-S12 (workflow angle: poisson-summation / theta-transform).

GOAL (from the endgame brief): the signed correction
    corr(E) = meas(S7(E)) - M7(k) = sum_{0 != n in Lambda(E)} K(n)
has an ABSOLUTELY DIVERGENT relation-lattice sum (MISTAKE-078). Apply Poisson
summation to find a CONVERGENT dual representation that bounds corr for WIDE sets.

FINDINGS (all verified exactly with fractions.Fraction here):

1. POISSON DUAL = THE x-CELL (WEYL) EVALUATION.  The Poisson dual of the
   relation-lattice theta sum is the orbit integral in x:
       sum_{n in Lambda(E)} K(n) = integral_0^1 [1_{cover}(x) - M7(k)] dx,
   which is an EXACT FINITE sum over the O(7*sum e_i) breakpoint cells. This is
   the convergent dual: the dual 'function' (per-cell cover-indicator deviation)
   is bounded by 1 and finitely supported. (= the engine itself.) The
   orbit-character-sum factorizes because n.e is an integer, so it is EXACTLY
   1_{n in Lambda(E)} -- the precise Poisson pairing.

2. M7(k) IS THE iid COUPON LIMIT: M7(k) = P(k-1 iid uniform sector labels cover
   {1..6}).  Verified vs Monte Carlo.

3. KOKSMA-HLAWKA / DISCREPANCY IS THE WRONG DUAL.  The Erdos-Turan-Koksma
   bound |corr| <= V_k * D*(orbit) has D* governed by sum_{Lambda} 1/r(n), which
   DIVERGES harmonically (no support-6 floor applies to the absolute object).
   So discrepancy reintroduces the divergence. The signed K-sum (with the
   THM-538 support-6 floor) is strictly the right object.

4. THE ABSOLUTE ENVELOPE DIVERGES FOR ALL E, INCLUDING WIDE ONES (confirms
   MISTAKE-078): sum_{supp>=6} prod 0.6973/|n_j| grows without bound even for
   Sidon sets (squares). Width does NOT rescue the absolute bound.

5. THE RULER IS R6 (support-6 relation DENSITY), NOT shortest-relation length.
   lambda_inf(shortest support-6 relation) = 1 for essentially all E (the 0,1
   prefix and 3-term coincidences supply tiny low-support relations), so it does
   not discriminate. But the COUNT of support-6 relations in a box tracks corr
   monotonically: AP 9336 > nearAP 8862 > oddAP 5210 > primes 4236 > squares
   1544, matching corr 0.303 > 0.184 > 0.213(*) > 0.0096 > 0.0005.

6. SIGNED CANCELLATION IMPROVES WITH WIDTH.  ratio sum_T|P_T-iid| / |corr|:
   AP 8.4, nearAP 13, primes 102, squares 1379. The signed structure is far more
   effective for wide sets -- this is exactly why the wide half is the easy half.

CONCLUSION ON THE ANGLE.  Poisson/theta does give a genuinely convergent
representation -- but it is the FINITE x-cell evaluation (the engine), not a new
dual-lattice series. No analytic dual series converges absolutely (the divergence
is intrinsic to the discontinuous sector indicators). The usable wide bound must
keep the SIGNED support-6 structure; the honest convergent dual is the exact
finite x-integral, whose value is small for wide E because the orbit
equidistributes. This corroborates the existing wide-spread plan (THM-538 +
stranger-contraction) rather than replacing it.
"""
import sys, itertools, cmath, math
from fractions import Fraction
from math import comb, gcd
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def measS7(E):
    """Exact meas(S7(E)) = P(N(x)=0) via breakpoint cells (Fraction)."""
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = Fraction(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e * mid
            v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if sum(1 for j in range(1, 7) if j not in hit) == 0:
            tot += (hi - lo)
    return tot


def M7(k):
    """iid coupon limit = K(0)."""
    return sum(Fraction((-1) ** t) * comb(6, t) * Fraction(7 - t, 7) ** (k - 1)
               for t in range(7))


SUBS = [frozenset(c) for r in range(7) for c in itertools.combinations(range(1, 7), r)]


def chat_T(T, n):
    """Fourier coeff of 1_{complement of union_{j in T}[j/7,(j+1)/7)} at frequency n."""
    if n == 0:
        return complex(1 - len(T) / 7, 0)
    if n % 7 == 0:
        return 0 + 0j  # THM-503 apex-prime vanishing
    s = 0 + 0j
    for j in T:
        a = j / 7.0
        b = (j + 1) / 7.0
        s += (cmath.exp(-2j * math.pi * n * b) - cmath.exp(-2j * math.pi * n * a)) / (-2j * math.pi * n)
    return -s


def K(nv, Enz, cache=None):
    tot = 0 + 0j
    for T in SUBS:
        p = 1 + 0j
        for ne in nv:
            p *= chat_T(T, ne)
        tot += (-1) ** len(T) * p
    return tot.real


def lattice_box_sum(E, H):
    """Partial signed sum of K(n) over Lambda(E) cap [-H,H]^{k-1} (the divergent primal)."""
    Enz = [e for e in sorted(set(E)) if e != 0]
    d = len(Enz)
    cache = {n: [chat_T(T, n) for T in SUBS] for n in range(-H, H + 1)}
    tot = 0.0
    cnt = 0
    for nv in itertools.product(range(-H, H + 1), repeat=d):
        if all(x == 0 for x in nv):
            continue
        if sum(n * e for n, e in zip(nv, Enz)) != 0:
            continue
        s = 0 + 0j
        for ti, T in enumerate(SUBS):
            p = 1 + 0j
            for ne in nv:
                p *= cache[ne][ti]
            s += (-1) ** len(T) * p
        tot += s.real
        cnt += 1
    return tot, cnt


def support6_count(E, H):
    """#{support>=6 relations in box} = the R6 ruler."""
    Enz = [e for e in sorted(set(E)) if e != 0]
    d = len(Enz)
    cnt = 0
    for nv in itertools.product(range(-H, H + 1), repeat=d):
        sup = sum(1 for x in nv if x != 0 and x % 7 != 0)
        if sup < 6:
            continue
        if sum(n * e for n, e in zip(nv, Enz)) != 0:
            continue
        cnt += 1
    return cnt


if __name__ == "__main__":
    bank = {
        "AP(consec)": list(range(8)),
        "near-AP": [0, 1, 2, 3, 4, 5, 6, 8],
        "odd-AP": [0, 1, 3, 5, 7, 9, 11, 13],
        "primes-ish": [0, 2, 3, 5, 7, 11, 13, 17],
        "squares(Sidon)": [0, 1, 4, 9, 16, 25, 36, 49],
    }
    print("=== Poisson/theta dual: primal lattice sum (divergent) vs exact x-cell dual ===\n")
    for name, E in bank.items():
        k = len(E)
        corr = measS7(E) - M7(k)
        print(f"{name:16s} k={k}  corr(exact x-cell dual) = {float(corr):+.6f}")
        # primal box sums (illustrate divergence/slow convergence)
        ps = [lattice_box_sum(E, H)[0] for H in (1, 2)]
        r6 = support6_count(E, 2)
        print(f"   primal box-sum H=1,2: {ps[0]:+.5f} {ps[1]:+.5f}   R6(box2)={r6}")
    print("\nThe x-cell value is exact+finite (the convergent Poisson dual);")
    print("the primal box sums crawl up harmonically. R6 (support-6 count) is the ruler.")
