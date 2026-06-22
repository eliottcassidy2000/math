#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_witness_2q_synthesis_kpswf11.py  (kind-pasteur 2026-06-22, THREAD 1 deliverable)

THE DELIVERABLE: per-q floors, the full slack table, and the q-uniform witness-floor
pattern for LRC(n=2q), q in {3,5,7}.

Self-contained (recomputes everything exactly).  Uses the canonical-worst-P shortlist
+ consec cluster, which the companion full-enumeration scripts (v2 + robust checks)
confirmed is the binding shape.  This script is fast (consec + worst-P only) and
reproduces the exact floors and slacks.

OUTPUT (the three Thread-1 deliverables):
  (I)   per-q admissible floors m_P^{(2q)}  and floor-case binding phi_q
  (II)  the full slack table  min-G2(k) / m_P  for k=q..2q-1
  (III) the q-uniform pattern + the uniform analytic lower bound 1/q
"""
import itertools
import sys
from fractions import Fraction as Fr


# ----- exact arc primitives -----
def cmg(E, x):
    ph = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(ph) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(ph, ph[1:]):
        if b - a > g:
            g = b - a
    return max(g, (ph[0] + 1) - ph[-1])


def garcs(E, q):
    gt = Fr(1, q)
    bps = {Fr(0), Fr(1)}
    diffs = set()
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d:
                diffs.add(d)
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        for m in range(0, q * d + 1):
            for s in (1, -1):
                v = Fr(q * m + s, q * d)
                if 0 < v < 1:
                    bps.add(v)
    bps = sorted(bps)
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        if cmg(E, (x0 + x1) / 2) > gt:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def parcs(P, n):
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
    bps = sorted(bps)
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if all(min((Fr(p) * mid) % 1, 1 - ((Fr(p) * mid) % 1)) >= thr for p in P):
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def im(A, B):
    i = j = 0
    t = Fr(0)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            t += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return t


def am(a):
    return sum((b - x for x, b in a), Fr(0))


def meas_GP(P, n):
    return am(parcs(list(P), n))


def min_measGP(n, npart):
    best = None
    bP = None
    for P in itertools.combinations(range(1, n), npart):
        m = meas_GP(P, n)
        if best is None or m < best:
            best = m
            bP = P
    return best, bP


def minG2_consec(n, q, k):
    """min over worst P (this |P|) of G2(P, consec_k).  Consec is the binding shape
       (verified by full-enumeration companions)."""
    npart = (n - 1) - k
    E = list(range(k))
    ga = garcs(E, q)
    if npart == 0:
        return am(ga), ()
    best = None
    bP = None
    for P in itertools.combinations(range(1, n), npart):
        r = im(ga, parcs(list(P), n))
        if best is None or r < best:
            best = r
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
    print("=" * 86)
    print("THREAD 1 DELIVERABLE: WITNESS ROUTE for LRC(2q), q in {3,5,7}  (kpswf11)")
    print("=" * 86)

    floors = {}
    tables = {}
    for q in (3, 5, 7):
        n = 2 * q
        mP, argP = min_measGP(n, n - 4)        # admissible floor (|P|=n-4, k=3)
        phi, argphi = min_measGP(n, q - 1)     # floor-case binding (|P|=q-1, k=q)
        floors[q] = (mP, argP, phi, argphi, n)
        rows = []
        for k in range(q, 2 * q):
            g2, P = minG2_consec(n, q, k)
            rows.append((k, (n - 1) - k, g2, g2 / mP))
        tables[q] = rows
        sys.stdout.flush()

    # (I) floors
    print("\n(I) PER-q ADMISSIBLE FLOORS")
    print(f"  {'q':>3} {'n':>4} {'LRC(n)':>8} {'m_P (|P|=n-4)':>18} {'=dec':>9} "
          f"{'phi_q (|P|=q-1)':>18} {'=dec':>9}")
    print("  " + "-" * 80)
    for q in (3, 5, 7):
        mP, _, phi, _, n = floors[q]
        st = "PROVED" if n < 14 else "OPEN"
        print(f"  {q:>3} {n:>4} {st:>8} {str(mP):>18} {float(mP):>9.6f} "
              f"{str(phi):>18} {float(phi):>9.6f}")
    print("\n  m_P denominators (number-theoretic, NOT a clean closed form):")
    for q in (3, 5, 7):
        mP, _, _, _, n = floors[q]
        print(f"    q={q} (n={n}): {mP.denominator} = {factor(mP.denominator)}")

    # (II) slack table
    print("\n(II) THE SLACK TABLE  min-G2(k) / m_P  over floor cases k=q..2q-1")
    for q in (3, 5, 7):
        mP, _, _, _, n = floors[q]
        print(f"\n  q={q} (n={n}, m_P={float(mP):.6f}):")
        print(f"    {'k':>3} {'|P|':>4} {'min G2':>16} {'=dec':>9} {'G2/m_P':>8} "
              f"{'regime':>10}")
        for (k, npart, g2, ratio) in tables[q]:
            reg = "BINDING" if k == q else "dense"
            print(f"    {k:>3} {npart:>4} {str(g2):>16} {float(g2):>9.6f} "
                  f"{float(ratio):>8.3f} {reg:>10}")
        dense = [r for r in tables[q] if r[0] >= q + 1]
        if dense:
            print(f"    => DENSE-k (k>=q+1) slack in [{float(min(r[3] for r in dense)):.2f}, "
                  f"{float(max(r[3] for r in dense)):.2f}]x")

    # (III) q-uniform pattern
    print("\n(III) THE q-UNIFORM WITNESS-FLOOR PATTERN")
    print("""
  * PIGEONHOLE closes k < q UNIFORMLY: k points => maxgap >= 1/k > 1/q for every x
    => GOOD(E) = [0,1) => G2 = meas(G_P) >= m_P.  Elementary, no E-structure.
  * The k=q BOUNDARY: consec_q = {0,..,q-1} has maxgap > 1/q for a.e. x (the
    equality set x=a/q has measure 0).  So GOOD(consec_q)=[0,1) a.e. and the binding
    floor-case value is min G2 = phi_q = min_{|P|=q-1} meas(G_P).
  * phi_q >= m_P for all q (the floor holds at the binding case) with WIDENING margin:""")
    print(f"      {'q':>3} {'phi_q/m_P':>10}")
    for q in (3, 5, 7):
        mP, _, phi, _, n = floors[q]
        print(f"      {q:>3} {float(phi / mP):>10.3f}")
    print("""  * DENSE k>q: |P|=n-1-k < q-1 shrinks => meas(G_P) grows => slack grows in k
    AND in q (q=3: ~1.1x, q=5: ~3x, q=7: ~5-8x).
  * UNIFORM analytic lower bound (union bound): for the floor cases |P| <= q-1,
        meas(G_P) >= 1 - 2|P|/n >= 1 - 2(q-1)/(2q) = 1/q > 0,
    so the witness floor is bounded below by 1/q UNIFORMLY in q.  The exact min is
    several-fold larger (phi_q ~ 1/4 as q grows, vs 1/q crude bound).
  * Since LRC(6) and LRC(10) are PROVED, q=3 and q=5 are a CONSISTENCY CHECK: the
    witness floor MUST be positive there -- and it is (tight 1.0x at q=3, 2.7x at
    q=5).  The SAME mechanism gives the q=7 floor (4.97x at the binding k=7), which
    is the LRC(14) witness floor.""")
    print("\nDONE.")


if __name__ == "__main__":
    main()
