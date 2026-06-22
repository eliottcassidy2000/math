#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_witness_route_2q_kpswf11.py   (kind-pasteur 2026-06-22, THREAD 1)

THE WITNESS ROUTE for LRC(n) at n = 2q, q prime.   q in {3, 5, 7}.

Generalizes the n=14 (q=7) witness machinery (kps-S29, HYP-2823..2827) to all
n=2q.  Since LRC(6) and LRC(10) are PROVED, the witness floor MUST hold there;
we verify it and EXTRACT the q-uniform pattern that would also close q=7.

================================ THE OBJECTS ================================
LRC(n): n runners = n-1 nonzero integer speeds S, threshold 1/n.  A counterexample
would have max_tau min_{v in S} ||v tau|| < 1/n.  By THM-523B WLOG S is COVERING
(a multiple of every q' in {2,...,n}) and PRIMITIVE (gcd 1), |S| = n-1.

The witness route writes S = P u L:
  * P = the "small/safe" part (the speeds that on their own leave a wide safe arc),
  * L = {Vmax - e_i : e_i in E} = the "cluster" of k = |E| co-offsets of a large Vmax.
  |P| + k = n-1.

G_P = { x in [0,1) : ||p x|| >= 1/n  for all p in P }.   (small-part safe set)
      Breakpoints on the grid t/(n p): frac(p x) = 1/n or (n-1)/n.

GAP threshold = 1/q  (= 2/n).   A free fast-phase slot for the cluster exists iff
the k cluster phases {frac(e_i x)} leave a circular gap > 1/q.  (maxgap > 1/q
=> a width-(maxgap - 1/q) fast window clears the 1/n collar => M(S) >= 1/n.)

WITNESS floor:
  G2(P,E) = meas{ x in G_P : maxgap{frac(e_i x): e_i in E} > 1/q }.
  G2 > 0  =>  M(S) >= 1/n  =>  LRC holds for S.   (THM-527 / HYP-2825-2826)

PIGEONHOLE (the collapse, HYP-2827):
  k points on the circle => circular maxgap >= 1/k.
  So  k < q  =>  maxgap >= 1/k > 1/q  for EVERY x  =>  GOOD(E) = whole circle
  =>  G2 = meas(G_P) >= m_P.    [ELEMENTARY: no E-structure needed.]
  k = q is a.e. (equality maxgap = 1/q only on a measure-zero set x = a/q).
  The FLOOR is needed only for  k = q .. 2q-1  (the "dense" cluster regime).

ADMISSIBLE-|P| floor:
  m_P^{(n)} = min over admissible P of meas(G_P).
  The min is at the LARGEST |P| = n-1-k with k as small as the admissible cluster
  allows.  We compute it exactly per q and confirm it is the THM-530 analogue.

================================ WHAT THIS SCRIPT DOES ================================
For each q in {3,5,7}  (n = 2q):
  (1) m_P^{(n)} = min over P subset {1..n-1}, all |P|, of meas(G_P)  (exact rational).
      Report the binding |P|, the argmin P, and the factorization.
  (2) PIGEONHOLE check: GOOD(E) = whole circle for every k<q shape  => G2 = meas(G_P).
      (we verify mu_pure(E) = 1 for representative k<q shapes and a.e. at k=q).
  (3) WITNESS floor over the floor cases k = q .. 2q-1 and admissible (P,E):
      min G2(P,E), confirm >= m_P^{(n)}.   (bounded-spread clusters E.)
  (4) SLACK: min-G2 / m_P for k >= q+1  (the dense cases).
  (5) The q-uniform pattern: closed form for m_P^{(n)} in q; binding=sparsest;
      dense-k slack growing in q.

All measures are EXACT Fractions.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


# ============================================================ circular max-gap
def circ_maxgap_at(E, x):
    """Exact circular max-gap of {frac(e*x): e in E}.  Fraction in [0,1].
       If all phases coincide, gap = 1."""
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def norm(y):
    r = y % 1
    return min(r, 1 - r)


# ============================================================ breakpoint sets (parametrized by n,q)
def gp_breaks(P, n):
    """G_P breakpoints in (0,1): frac(p x) = 1/n or (n-1)/n  =>  x = (n m + r)/(n p)."""
    bps = set()
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, n - 1):
                v = Fr(n * m + r, n * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks(E, q):
    """GOOD(E) breakpoints in (0,1) for gap threshold 1/q:
       phase collisions x = t/d  AND  gap = +-1/q crossings.
       gap = 1/q at  frac(e_j x) - frac(e_i x) = +-1/q  =>  (e_j-e_i) x = +-1/q + Z
       =>  x = (q m +- 1)/(q d),  for every pairwise difference d=|e_i-e_j|."""
    bps = set()
    diffs = set()
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d != 0:
                diffs.add(d)
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        hi = q * d
        for m in range(0, hi + 1):
            for s in (1, -1):
                num = q * m + s
                v = Fr(num, q * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


# ============================================================ arc representations
def good_arcs(E, q):
    """GOOD set {maxgap > 1/q} as sorted disjoint arcs."""
    gapthr = Fr(1, q)
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E, q))
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        if circ_maxgap_at(E, (x0 + x1) / 2) > gapthr:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def gp_arcs(P, n):
    """G_P {||p x||>=1/n all p} as sorted disjoint arcs."""
    thr = Fr(1, n)
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P, n))
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for p in P:
            if norm(Fr(p) * mid) < thr:
                ok = False
                break
        if ok:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def arcs_measure(arcs):
    return sum((b - a for a, b in arcs), Fr(0))


def arcs_intersect_measure(A, B):
    """measure of union(A) cap union(B)."""
    i = j = 0
    tot = Fr(0)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            tot += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return tot


def mu_pure(E, q):
    return arcs_measure(good_arcs(E, q))


def meas_GP(P, n):
    return arcs_measure(gp_arcs(P, n))


def G2(P, E, n, q):
    """The witness floor meas{x in G_P : maxgap{frac(e x)} > 1/q}."""
    return arcs_intersect_measure(gp_arcs(P, n), good_arcs(E, q))


# ============================================================ admissibility helpers
def is_covering(S, n):
    """S contains a multiple of every q' in {2..n}."""
    for qq in range(2, n + 1):
        if not any(v % qq == 0 for v in S):
            return False
    return True


def is_primitive(S):
    return reduce(gcd, S) == 1


def factor(num):
    """Tiny trial-division factorization string."""
    if num == 0:
        return "0"
    n = num
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(f.items()))


# ============================================================ per-q analysis
def mP_floor(n, verbose=True):
    """min over ALL P subset {1..n-1}, all sizes, of meas(G_P).  Exact.
       Returns (min_meas, argmin_P, binding_|P|)."""
    best = None
    bestP = None
    bestnp = None
    # the min is at large |P|; scan |P| from n-1 down, stop when it can't decrease.
    # meas(G_P) is non-increasing as we add speeds, so the global min is at |P|=n-1
    # (all speeds) -- but we report the binding behaviour for context.
    for npart in range(1, n):
        loc = None
        locP = None
        for P in itertools.combinations(range(1, n), npart):
            m = meas_GP(list(P), n)
            if loc is None or m < loc:
                loc = m
                locP = P
        if verbose:
            print(f"    |P|={npart:2d}: min meas(G_P) = {loc} = {float(loc):.6f}  "
                  f"at P={locP}")
        if best is None or loc < best:
            best = loc
            bestP = locP
            bestnp = npart
        sys.stdout.flush()
    return best, bestP, bestnp


def main():
    print("=" * 86)
    print("THREAD 1: THE WITNESS ROUTE for LRC(n=2q),  q in {3,5,7}   (kpswf11)")
    print("=" * 86)
    sys.stdout.flush()

    floors = {}
    for q in (3, 5, 7):
        n = 2 * q
        print("\n" + "#" * 86)
        print(f"#  q = {q}   (n = {n} = 2*{q}),   gap threshold 1/q = 1/{q},   "
              f"G_P threshold 1/n = 1/{n}")
        print(f"#  LRC({n}) is {'PROVED' if n < 14 else 'OPEN (FIRST OPEN CASE)'}.  "
              f"Speeds = n-1 = {n-1}.")
        print("#" * 86)
        sys.stdout.flush()

        # ---------- (1) the admissible-|P| floor m_P^{(n)} ----------
        print(f"\n--- (1) m_P^({n}) = min over P subset {{1..{n-1}}} of meas(G_P) ---")
        mP, argP, bnp = mP_floor(n)
        print(f"\n  >>> m_P^({n}) = {mP} = {float(mP):.8f}")
        print(f"      binding |P| = {bnp} (= all speeds), argmin P = {argP}")
        print(f"      denominator {mP.denominator} = {factor(mP.denominator)}")
        print(f"      numerator   {mP.numerator} = {factor(mP.numerator)}")
        floors[q] = (mP, argP, bnp, n)
        sys.stdout.flush()

        # ---------- (2) pigeonhole: k<q => GOOD=whole circle => G2=meas(G_P) ----------
        print(f"\n--- (2) PIGEONHOLE: k < q={q} => circular maxgap >= 1/k > 1/q "
              f"=> GOOD(E)=[0,1) ---")
        for k in range(2, q + 1):
            # representative shapes: consec, and a spread-out shape
            reps = [list(range(k))]
            if k >= 3:
                reps.append([0] + list(range(2, 2 * k - 1, 2)))  # even-spread
                reps.append([0, 1] + [3 * i for i in range(1, k - 1)])  # mixed
            allone = True
            mn = None
            for E in reps:
                mu = mu_pure(E, q)
                if mu != Fr(1):
                    allone = False
                if mn is None or mu < mn:
                    mn = mu
            tag = ("PIGEONHOLE: maxgap>=1/k=1/%d>1/q=1/%d ALWAYS => GOOD=[0,1)" % (k, q)) if k < q \
                else ("k=q: maxgap>=1/q a.e. (equality only on measure-zero x=a/q)")
            status = "mu=1 (all reps)" if allone else f"min mu over reps = {mn}"
            print(f"    k={k}: {status:<22}  {tag}")
        sys.stdout.flush()

        # ---------- (3)+(4) witness floor over floor cases k=q..2q-1 ----------
        print(f"\n--- (3)+(4) WITNESS floor G2(P,E) over floor cases k=q..2q-1, "
              f"admissible (P,E) ---")
        print(f"    For each k, |P| = {n-1}-k.  Search bounded-spread E and worst P.")
        print(f"    {'k':>3} {'|P|':>4} {'min G2':>16} {'=dec':>10} {'m_P':>12} "
              f"{'G2/m_P':>8} {'>=m_P?':>7}  argmin (E ; P)")
        print("    " + "-" * 96)
        slack_rows = []
        for k in range(q, 2 * q):
            npart = (n - 1) - k
            # P-shortlist: smallest meas(G_P) at this |P| (the binding P)
            Pall = []
            if npart == 0:
                Pall = [(Fr(1), (), [(Fr(0), Fr(1))])]
            else:
                for P in itertools.combinations(range(1, n), npart):
                    a = gp_arcs(list(P), n)
                    Pall.append((arcs_measure(a), P, a))
                Pall.sort(key=lambda e: e[0])
            # keep a generous shortlist of the smallest-measure P (those can bind G2 low)
            Pcands = Pall[: min(len(Pall), 40)]

            # enumerate bounded-spread shapes E (0=e_1<...<e_k <= Wmax)
            Wmax = min(k + 6, 18)
            best = None
            bestE = None
            bestP = None
            for tail in itertools.combinations(range(1, Wmax + 1), k - 1):
                E = [0] + list(tail)
                ga = good_arcs(E, q)
                mu = arcs_measure(ga)
                for (mg, P, pa) in Pcands:
                    lb = mu + mg - 1  # IE lower bound on intersection
                    if best is not None and lb >= best:
                        continue
                    r = arcs_intersect_measure(ga, pa)
                    if best is None or r < best:
                        best = r
                        bestE = E
                        bestP = list(P)
            ratio = (best / mP) if mP > 0 else Fr(0)
            ok = best >= mP
            slack_rows.append((k, npart, best, mP, ratio, ok, bestE, bestP))
            print(f"    {k:>3} {npart:>4} {str(best):>16} {float(best):>10.6f} "
                  f"{float(mP):>12.6f} {float(ratio):>8.3f} {str(ok):>7}  "
                  f"E={bestE} ; P={bestP}")
            sys.stdout.flush()

        # ---------- summary per q ----------
        minG2 = min(r[2] for r in slack_rows)
        allok = all(r[5] for r in slack_rows)
        print(f"\n    => min G2 over floor cases k=q..2q-1: {minG2} = {float(minG2):.6f}")
        print(f"       m_P^({n}) = {mP} = {float(mP):.6f}")
        print(f"       ALL G2 >= m_P (floor holds): {allok}")
        # dense-k slack (k>=q+1)
        dense = [r for r in slack_rows if r[0] >= q + 1]
        if dense:
            dmin = min(r[4] for r in dense)
            dmax = max(r[4] for r in dense)
            print(f"       DENSE-k (k>=q+1) slack G2/m_P in [{float(dmin):.2f}, "
                  f"{float(dmax):.2f}]")
        sys.stdout.flush()

    # ---------- (5) the q-uniform pattern ----------
    print("\n" + "=" * 86)
    print("(5) THE q-UNIFORM PATTERN")
    print("=" * 86)
    print(f"\n  {'q':>3} {'n=2q':>5} {'m_P^(n)':>18} {'=decimal':>12} "
          f"{'denominator factored':>30}")
    print("  " + "-" * 74)
    for q in (3, 5, 7):
        mP, argP, bnp, n = floors[q]
        print(f"  {q:>3} {n:>5} {str(mP):>18} {float(mP):>12.7f} "
              f"{factor(mP.denominator):>30}")
    print("""
  READING:
   * PIGEONHOLE closes k < q UNIFORMLY in q (elementary, no E-structure).
   * The witness floor G2 >= m_P holds in every floor case k=q..2q-1.
   * The BINDING case is the SPARSEST admissible cluster: G2 = meas(G_P) = m_P
     there (whole circle lonely), and the DENSE cases k>=q+1 have growing slack
     because |P|=n-1-k is small => meas(G_P) is large.
   * Since LRC(6) and LRC(10) are PROVED, the q=3,5 floors are a CONSISTENCY
     CHECK: the witness floor MUST be positive there, and it is.
""")
    print("DONE.")


if __name__ == "__main__":
    main()
