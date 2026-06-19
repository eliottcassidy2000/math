#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFY of the LRC(n) runner-recursion angle.
Completely independent reimplementation (fractions.Fraction, no copied logic).

Three claims to check:
(1) Factorization table for C=2n-1, n=8..20; n=14 -> C=27=3^3 unique non-trivial prime power.
(2) AP-to-cap margin on uniform binding row |P|=5; n=14 tightest (+0.054, only one < 0.10).
(3) Closed form: p_all-missed(consec_k) = 1/(s*(k-1)), s=n/2, EXACTLY for all even n, all k.

Also: decoupling table E={0,1,2,3}+{0,d} meas_S << cap.

Plus independent counterexample hunt on the recursion (3) and on (2)'s extremality.
"""
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
import itertools

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')


def factor(m):
    f = {}
    d = 2
    x = m
    while d * d <= x:
        while x % d == 0:
            f[d] = f.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        f[x] = f.get(x, 0) + 1
    return f


def frac_part(v):
    """fractional part of Fraction v in [0,1)."""
    return v - (v.numerator // v.denominator)


def sectors_hit(E, x, s):
    """set of sectors 0..s-1 hit by {frac(e*x): e in E}."""
    hit = set()
    for e in E:
        v = frac_part(e * x)
        hit.add((v.numerator * s) // v.denominator)
    return hit


def dist_N(E, s):
    """
    Exact distribution of N(x) = #{ j in 1..s-1 : sector j missed }, over x in [0,1).
    Returns list p[0..s-1] of Fractions summing to 1.
    Breakpoints: x where some e*x crosses a wall a/s, i.e. x = a/(s*e).
    Independent: I sample the OPEN midpoint of each elementary interval.
    """
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        # walls of e*x at integer+j/s ; x = (i + j/s)/e for the relevant range
        # equivalently x = a/(s*e) for a = 0..s*e
        for a in range(0, s * e + 1):
            bps.add(Fraction(a, s * e))
    bps = sorted(b for b in bps if Fraction(0) <= b <= Fraction(1))
    p = [Fraction(0)] * s
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        hit = sectors_hit(E, mid, s)
        N = sum(1 for j in range(1, s) if j not in hit)
        p[N] += (hi - lo)
    return p


def meas_S(E, s):
    """measure of x where ALL sectors 1..s-1 are hit (N=0)."""
    return dist_N(E, s)[0]


def p_all_missed(E, s):
    """measure of x where ALL s-1 sectors are missed (N = s-1)."""
    return dist_N(E, s)[s - 1]


# -------- cap_k as canonical min over complement sets P, |P| = n-1-k --------
# The lonely-runner cover: speeds are 1..n-1 (n-1 nonzero), our E sits inside.
# cap_k = min over |P|=n-1-k of meas(G_P) where G_P is the "all sectors hit" set
# for the set P (plus the runner himself at speed... ). Per thread, cap is min meas_S
# over complementary configs of size n-1-k. We reproduce THM-530 canonical:
# for n=14, cap_8 = min over |P|=5 of meas_S(P, s=7). Let's just compute meas_S over
# all size-(n-1-k) subsets of speeds {1..n-1} and take the min.

def cap_k(n, k):
    s = n // 2
    speeds = list(range(1, n))   # 1..n-1
    psize = (n - 1) - k
    if psize < 0:
        return None
    best = None
    bestP = None
    for P in itertools.combinations(speeds, psize):
        # include 0 (the runner reference) as in E having 0
        E = (0,) + P
        m = meas_S(E, s)
        if best is None or m < best:
            best = m
            bestP = P
    return best, bestP


def main():
    ns = [8, 10, 12, 14, 16, 18, 20]

    print("=" * 70)
    print("CLAIM 1: factorization of C=2n-1, n=14 unique non-trivial prime power")
    print("=" * 70)
    print(f"{'n':>3} {'C=2n-1':>7} {'factor':>12} {'#div':>5} {'maxexp':>7} {'pp?':>5}")
    for n in ns:
        C = 2 * n - 1
        f = factor(C)
        facstr = "*".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(f.items()))
        ndiv = 1
        for e in f.values():
            ndiv *= (e + 1)
        maxexp = max(f.values())
        is_pp = len(f) == 1
        is_ntpp = is_pp and maxexp >= 2
        tag = "NTPP" if is_ntpp else ("pp" if is_pp else "")
        print(f"{n:>3} {C:>7} {facstr:>12} {ndiv:>5} {maxexp:>7} {tag:>5}")

    print("\n" + "=" * 70)
    print("CLAIM 3: p_all-missed(consec_k) = 1/(s*(k-1)) EXACTLY (s=n/2)")
    print("=" * 70)
    mismatches = 0
    total = 0
    for n in ns:
        s = n // 2
        # test a range of k
        for k in range(2, n):
            E = list(range(k))   # consec, includes 0
            am = p_all_missed(E, s)
            pred = Fraction(1, s * (k - 1))
            total += 1
            ok = (am == pred)
            if not ok:
                mismatches += 1
                print(f"  MISMATCH n={n} k={k}: actual={am}  pred={pred}")
        print(f"  n={n} (s={s}): k=2..{n-1} all checked")
    print(f"\n  TOTAL checked = {total}, mismatches = {mismatches}")

    print("\n" + "=" * 70)
    print("CLAIM 2: AP-to-cap margin on uniform binding row |P|=5 (k = n-1-5 = n-6)")
    print("Wait: |P|=5 means k = n-1-5. Check the claim's framing precisely.")
    print("=" * 70)
    # The claim says: uniform binding row |P|=5, AP-to-cap margin per n.
    # |P|=5 => k = (n-1) - 5 = n-6.  For n=14: k=8. matches "dangerous k=8".
    # margin = cap_k - meas_S(consec_k)  (AP=consec is below cap; margin = cap - AP)
    print(f"{'n':>3} {'k=n-6':>6} {'meas_S(AP)':>14} {'cap':>14} {'margin':>10}")
    margins = {}
    for n in ns:
        s = n // 2
        k = n - 6
        if k < 1:
            print(f"{n:>3}  k={k} (too small, skip)")
            continue
        ap = meas_S(list(range(k)), s)
        cap, bestP = cap_k(n, k)
        margin = cap - ap
        margins[n] = float(margin)
        flag = "  <-- TIGHTEST" if False else ""
        print(f"{n:>3} {k:>6} {str(ap):>14} {str(cap):>14} {float(margin):>10.4f}  AP={float(ap):.4f} cap={float(cap):.4f}")
    print(f"\n  margins: {margins}")
    if margins:
        tn = min(margins, key=margins.get)
        print(f"  TIGHTEST margin at n={tn} ({margins[tn]:.4f})")
        below10 = [n for n, m in margins.items() if m < 0.10]
        print(f"  margins < 0.10 : {below10}")

    print("\n" + "=" * 70)
    print("DECOUPLING: E={0,1,2,3}+{0,d}  meas_S(s=7) << cap_8 ~ 0.381")
    print("=" * 70)
    s = 7
    for d in [5, 6, 7, 8, 12, 100]:
        A = [0, 1, 2, 3]
        E = sorted(set(a + b for a in A for b in [0, d]))
        ms = meas_S(E, s)
        print(f"  d={d:>4}: E={E}  meas_S={float(ms):.5f}")

    print("\n" + "=" * 70)
    print("VALIDATION against canon: cap_8(n=14)=2243/5880, meas_S(consec_8,s=7)=481/1470")
    print("=" * 70)
    cap8, P8 = cap_k(14, 8)
    print(f"  cap_8(n=14) = {cap8} = {float(cap8):.6f}  (canon 2243/5880={float(Fraction(2243,5880)):.6f})  bestP={P8}")
    print(f"    matches 2243/5880? {cap8 == Fraction(2243,5880)}")
    ms8 = meas_S(list(range(8)), 7)
    print(f"  meas_S(consec_8,s=7) = {ms8} = {float(ms8):.6f}  (canon 481/1470={float(Fraction(481,1470)):.6f})")
    print(f"    matches 481/1470? {ms8 == Fraction(481,1470)}")
    pam = p_all_missed(list(range(8)), 7)
    print(f"  p_all-missed(consec_8,s=7) = {pam}  (canon 1/49={float(Fraction(1,49)):.6f})")
    print(f"    matches 1/49? {pam == Fraction(1,49)}")


if __name__ == "__main__":
    main()
