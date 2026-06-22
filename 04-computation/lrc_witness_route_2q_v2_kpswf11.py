#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_witness_route_2q_v2_kpswf11.py   (kind-pasteur 2026-06-22, THREAD 1, corrected)

THE WITNESS ROUTE for LRC(n=2q),  q in {3,5,7}.   [corrected m_P + faster]

Fixes over v1:
  * m_P^{(n)} = the ADMISSIBLE floor = min meas(G_P) over P COMPLETABLE to an
    admissible (covering+primitive) S = P u L with a k>=3 large cluster.
    Binding |P| = n-4 (k=3, the sparsest coordinated-growth core).  v1 wrongly
    let |P| run to n-1 (all speeds), which is degenerate (no cluster) and gives 0.
    For n=14 this reproduces the canon m_P = 14249/252252 at |P|=10.
  * Faster floor-case search (worst-P shortlist; bounded shapes), and it now runs
    ALL floor cases k=q..2q-1 for every q.

Objects (parametrized by n=2q):
  G_P     = {x : ||p x|| >= 1/n  all p in P}            (G_P threshold 1/n)
  GOOD(E) = {x : circ maxgap{frac(e x): e in E} > 1/q}  (gap threshold 1/q)
  G2(P,E) = meas(G_P cap GOOD(E))  >= rho*   [G2>0 => M(S)>=1/n => LRC holds]

Pigeonhole: k < q => maxgap >= 1/k > 1/q for every x => GOOD=[0,1) => G2=meas(G_P).
Floor needed only for k = q .. 2q-1.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


# ---------------- core measure primitives (parametrized) ----------------
def circ_maxgap_at(E, x):
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


def gp_breaks(P, n):
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
                v = Fr(q * m + s, q * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_arcs(E, q):
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
    thr = Fr(1, n)
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P, n))
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if all(norm(Fr(p) * mid) >= thr for p in P):
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def arcs_measure(arcs):
    return sum((b - a for a, b in arcs), Fr(0))


def arcs_intersect_measure(A, B):
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


def meas_GP(P, n):
    return arcs_measure(gp_arcs(P, n))


def mu_pure(E, q):
    return arcs_measure(good_arcs(E, q))


# ---------------- admissibility ----------------
def is_covering(S, n):
    return all(any(v % qq == 0 for v in S) for qq in range(2, n + 1))


def is_primitive(S):
    return reduce(gcd, S) == 1


def completable_k3(P, n, vmax_lim=400, spread_lim=None):
    """Is P completable to a covering+primitive admissible S = P u L with a
       k=3 large cluster L={Vmax-e_i}, all L>n-1?  Returns True/False (existence)."""
    P = list(P)
    if spread_lim is None:
        spread_lim = n + 2
    for Vmax in range(n, vmax_lim):
        for E in itertools.combinations(range(0, spread_lim + 1), 3):
            L = [Vmax - e for e in E]
            if min(L) <= n - 1:
                continue
            if len(set(L)) != 3:
                continue
            S = sorted(set(P) | set(L))
            if len(S) != len(P) + 3:
                continue
            if is_primitive(S) and is_covering(S, n):
                return True
    return False


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


# ---------------- the admissible floor m_P ----------------
def mP_admissible(n, verbose=True):
    """m_P^{(n)} = min meas(G_P) over admissible P.
       Binding at |P| = n-4 (k=3 cluster).  We compute the raw min meas(G_P) at
       |P|=n-4 AND verify the argmin is completable (covering+primitive).  We also
       report |P|=n-3 (k=2) for context."""
    results = {}
    # binding family: |P| = n-4 (k=3). Also show n-3 (k=2) and a couple above.
    for npart in (n - 3, n - 4):
        loc = None
        locP = None
        for P in itertools.combinations(range(1, n), npart):
            m = meas_GP(list(P), n)
            if loc is None or m < loc:
                loc = m
                locP = P
        results[npart] = (loc, locP)
        if verbose:
            kcl = (n - 1) - npart
            comp = completable_k3(locP, n) if npart == n - 4 else None
            ctag = "" if comp is None else f"  completable(k=3)={comp}"
            print(f"    |P|={npart:2d} (cluster k={kcl}): min meas(G_P) = {loc} = "
                  f"{float(loc):.7f}  at P={locP}{ctag}")
        sys.stdout.flush()
    mP, argP = results[n - 4]
    return mP, argP, n - 4, results


def worst_P_shortlist(n, npart, topk=30):
    """The smallest-meas(G_P) P's at this |P| (the binding ones for G2-min)."""
    if npart == 0:
        return [(Fr(1), (), [(Fr(0), Fr(1))])]
    entries = []
    for P in itertools.combinations(range(1, n), npart):
        a = gp_arcs(list(P), n)
        entries.append((arcs_measure(a), P, a))
    entries.sort(key=lambda e: e[0])
    return entries[:topk]


def main():
    print("=" * 90)
    print("THREAD 1: THE WITNESS ROUTE for LRC(n=2q),  q in {3,5,7}   (kpswf11 v2, corrected)")
    print("=" * 90)
    sys.stdout.flush()

    floors = {}
    slack_tables = {}
    for q in (3, 5, 7):
        n = 2 * q
        proved = "PROVED" if n < 14 else "OPEN (first open case)"
        print("\n" + "#" * 90)
        print(f"#  q={q}  n={n}=2*{q}   gap thr 1/q=1/{q}   G_P thr 1/n=1/{n}   "
              f"speeds={n-1}   LRC({n}) {proved}")
        print("#" * 90)
        sys.stdout.flush()

        # (1) admissible floor m_P
        print(f"\n--- (1) admissible floor m_P^({n}) (binding |P|=n-4=k=3 cluster) ---")
        mP, argP, bnp, byp = mP_admissible(n)
        print(f"\n  >>> m_P^({n}) = {mP} = {float(mP):.8f}   at |P|={bnp}, P={argP}")
        print(f"      denominator {mP.denominator} = {factor(mP.denominator)}")
        floors[q] = (mP, argP, bnp, n, byp)
        sys.stdout.flush()

        # (2) pigeonhole
        print(f"\n--- (2) PIGEONHOLE: k<q={q} => maxgap>=1/k>1/q => GOOD=[0,1) => G2=meas(G_P) ---")
        for k in range(2, q + 1):
            reps = [list(range(k))]
            if k >= 3:
                reps.append([0] + list(range(2, 2 * k - 1, 2)))
                reps.append([0, 1] + [3 * i for i in range(1, k - 1)])
            mus = [mu_pure(E, q) for E in reps]
            allone = all(m == Fr(1) for m in mus)
            if k < q:
                tag = f"maxgap>=1/{k}>1/{q} ALWAYS => GOOD=[0,1)  (mu=1 confirmed: {allone})"
            else:
                tag = f"k=q: GOOD=[0,1) a.e. (mu=1 confirmed: {allone}; equality x=a/q meas 0)"
            print(f"    k={k}: {tag}")
        sys.stdout.flush()

        # (3)+(4) witness floor over floor cases k=q..2q-1
        print(f"\n--- (3)+(4) WITNESS floor G2 over floor cases k=q..2q-1 ---")
        print(f"    {'k':>3} {'|P|':>4} {'min G2':>16} {'=dec':>10} "
              f"{'m_P':>11} {'G2/m_P':>8} {'>=m_P?':>7}  argmin (E ; P)")
        print("    " + "-" * 92)
        rows = []
        for k in range(q, 2 * q):
            npart = (n - 1) - k
            Pcands = worst_P_shortlist(n, npart, topk=30)
            # shapes: consec is the canonical worst (most-uniform cluster); also a
            # handful of perforated/spread shapes. Keep bounded.
            Wmax = min(k + 5, 17)
            best = None
            bestE = None
            bestP = None
            for tail in itertools.combinations(range(1, Wmax + 1), k - 1):
                E = [0] + list(tail)
                ga = good_arcs(E, q)
                mu = arcs_measure(ga)
                for (mg, P, pa) in Pcands:
                    lb = mu + mg - 1
                    if best is not None and lb >= best:
                        continue
                    r = arcs_intersect_measure(ga, pa)
                    if best is None or r < best:
                        best = r
                        bestE = E
                        bestP = list(P)
            ratio = best / mP if mP > 0 else Fr(0)
            ok = best >= mP
            rows.append((k, npart, best, ratio, ok, bestE, bestP))
            print(f"    {k:>3} {npart:>4} {str(best):>16} {float(best):>10.6f} "
                  f"{float(mP):>11.6f} {float(ratio):>8.3f} {str(ok):>7}  "
                  f"E={bestE} ; P={bestP}")
            sys.stdout.flush()
        slack_tables[q] = rows
        minG2 = min(r[2] for r in rows)
        allok = all(r[4] for r in rows)
        print(f"\n    => min G2 over k=q..2q-1: {minG2} = {float(minG2):.6f};  "
              f"m_P = {float(mP):.6f};  ALL G2>=m_P: {allok}")
        dense = [r for r in rows if r[0] >= q + 1]
        if dense:
            print(f"       DENSE-k (k>=q+1) slack G2/m_P in "
                  f"[{float(min(r[3] for r in dense)):.2f}, "
                  f"{float(max(r[3] for r in dense)):.2f}]")
        sys.stdout.flush()

    # (5) q-uniform pattern
    print("\n" + "=" * 90)
    print("(5) THE q-UNIFORM PATTERN -- per-q floors and the slack")
    print("=" * 90)
    print(f"\n  {'q':>3} {'n':>4} {'m_P^(n)':>20} {'=decimal':>11}  "
          f"{'denominator = factored':>34}")
    print("  " + "-" * 78)
    for q in (3, 5, 7):
        mP, argP, bnp, n, byp = floors[q]
        print(f"  {q:>3} {n:>4} {str(mP):>20} {float(mP):>11.7f}  "
              f"{factor(mP.denominator):>34}")
    print("\n  binding |P|=n-4 (k=3) floor AND the |P|=n-3 (k=2) floor, per q:")
    for q in (3, 5, 7):
        mP, argP, bnp, n, byp = floors[q]
        k2 = byp[n - 3]
        k3 = byp[n - 4]
        print(f"    q={q}: |P|={n-3} (k=2): {k2[0]}={float(k2[0]):.6f}   "
              f"|P|={n-4} (k=3): {k3[0]}={float(k3[0]):.6f}")
    print("\nDONE.")


if __name__ == "__main__":
    main()
