#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_witness_2q_closedform_kpswf11.py  (kind-pasteur 2026-06-22, THREAD 1 part 5)

Closed-form / structural analysis of the witness-route floors m_P^{(2q)} and the
q-uniform pattern.  Companion to lrc_witness_route_2q_v2_kpswf11.py.

Three quantities, per q (n=2q):
  (A) m_P^{(n)}  = admissible floor = min meas(G_P) over |P|=n-4 (k=3 cluster).
  (B) phi_q      = floor-case binding value = min meas(G_P) over |P|=q-1
                   (= min G2 at k=q, since consec cluster makes GOOD=[0,1) a.e.).
  (C) the simplest analytic anchor: meas(G_{1}) = 1 - 1/q = (q-1)/q  (single small
      runner), the trivial pigeonhole lower bound for any nonempty P.

We tabulate (A),(B) for q=2..8 (incl composite q for the pattern), factor the
denominators, and test candidate closed forms.  HONEST verdict on closed form.
"""
import itertools
import sys
from fractions import Fraction as Fr


def meas_GP(P, n):
    thr = Fr(1, n)
    bps = {Fr(0), Fr(1)}
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, n - 1):
                v = Fr(n * m + r, n * p)
                if 0 < v < 1:
                    bps.add(v)
    pts = sorted(bps)
    tot = Fr(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        if all(min((Fr(p) * mid) % 1, 1 - ((Fr(p) * mid) % 1)) >= thr for p in P):
            tot += b - a
    return tot


def min_measGP(n, npart):
    best = None
    bP = None
    for P in itertools.combinations(range(1, n), npart):
        m = meas_GP(list(P), n)
        if best is None or m < best:
            best = m
            bP = P
    return best, bP


def factor(num):
    if num == 0:
        return "0"
    m = num
    f = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(f.items()))


def main():
    print("=" * 84)
    print("THREAD 1 (5): closed-form / structural analysis of m_P^{(2q)}")
    print("=" * 84)

    qs = list(range(2, 9))  # include composite q to see the pattern is number-theoretic
    print(f"\n--- (A) admissible floor m_P = min meas(G_P) over |P|=n-4 (k=3) ---")
    print(f"    (and (B) floor-case binding phi_q = min meas(G_P) over |P|=q-1) ---")
    print(f"\n  {'q':>3} {'n':>4} {'prime?':>7} "
          f"{'m_P (|P|=n-4)':>18} {'=dec':>9}  {'phi_q (|P|=q-1)':>18} {'=dec':>9}")
    print("  " + "-" * 78)
    A = {}
    B = {}
    for q in qs:
        n = 2 * q
        isp = all(q % d for d in range(2, q)) and q > 1
        np_A = n - 4
        np_B = q - 1
        mA, _ = min_measGP(n, np_A) if np_A >= 1 else (Fr(0), ())
        mB, _ = min_measGP(n, np_B) if np_B >= 1 else (Fr(1), ())
        A[q] = mA
        B[q] = mB
        print(f"  {q:>3} {n:>4} {str(isp):>7} {str(mA):>18} {float(mA):>9.6f}  "
              f"{str(mB):>18} {float(mB):>9.6f}")
        sys.stdout.flush()

    print(f"\n--- denominators (number-theoretic structure) ---")
    for q in qs:
        n = 2 * q
        print(f"  q={q} (n={n}): m_P denom = {A[q].denominator} = {factor(A[q].denominator)}")
    print()
    for q in qs:
        n = 2 * q
        print(f"  q={q} (n={n}): phi_q denom = {B[q].denominator} = {factor(B[q].denominator)}")

    print(f"\n--- candidate closed forms (HONEST: test, don't assume) ---")
    print("  The denominators carry q (or q^2) and factors n-1, n-3 (=2q-1, 2q-3),")
    print("  i.e. lcm-like products from the grid t/(n p).  Test simple guesses:")
    for q in (3, 5, 7):
        n = 2 * q
        mP = A[q]
        # guesses
        g1 = Fr(1, 2 * q - 1)             # 1/(n-1)
        g2 = Fr(q - 1, q) ** 3            # pigeonhole-product (independent 3 runners)
        print(f"  q={q}: m_P={float(mP):.5f}  vs 1/(n-1)={float(g1):.5f}  "
              f"((q-1)/q)^3={float(g2):.5f}")
    print("\n  VERDICT: no single elementary closed form fits all q for the MIN floor")
    print("  (it is a discrete optimum over P-subsets, number-theoretic in q).")
    print("  The q-UNIFORM content is the MECHANISM + monotone slack, not a formula.")

    # the clean analytic lower bound that IS uniform:
    print(f"\n--- the UNIFORM analytic lower bound (this is the q-uniform anchor) ---")
    print("  For ANY nonempty P:  meas(G_P) >= 1 - sum_{p in P} 2/n = 1 - 2|P|/n  (union bd),")
    print("  and the floor-case |P| = n-1-k <= q-1 (k>=q), so")
    print("     meas(G_P) >= 1 - 2(q-1)/(2q) = 1 - (q-1)/q = 1/q  > 0   UNIFORMLY.")
    print("  Sharper (the binding floor case k=q, |P|=q-1):")
    for q in (3, 5, 7):
        n = 2 * q
        ub = 1 - Fr(2 * (q - 1), n)  # = 1/q
        print(f"     q={q}: crude union LB on meas(G_P), |P|=q-1: 1-2(q-1)/n = {ub} = "
              f"{float(ub):.5f};  actual min = {float(B[q]):.5f}  "
              f"(slack {float(B[q]/ub):.2f}x)")
    print("\n  => the witness floor in the floor cases is BOUNDED BELOW by 1/q UNIFORMLY")
    print("     (crude union bound), and the true min is several-fold larger.")
    print("\nDONE.")


if __name__ == "__main__":
    main()
