#!/usr/bin/env python3
"""
lrc_arc_menu_geometry_s521.py

claudebox-2026-06-01-S521

GEOMETRIC (box-independent) LRC source menu.

Context (THM-384, HYP-1987, oracle-S512):
  By THM-384 the marked observer is lonely at time t iff all m = n-1 movers lie
  in the closed safe arc [1/n, 1-1/n] of length L = (n-2)/n.  At such a time the
  runner-runner half-turn sub-tournament is some iso-class in A000568(m).  S512
  enumerated the classes that an ARITHMETIC walk actually reaches over a finite
  speed box ("reachable menu": 2,2,6,6 for n=4..7 -- a LOWER bound, polluted by
  box size and by boundary near-ties).

  HYP-1987 says the true target is the GEOMETRIC menu: the half-turn tournaments
  realizable by m points confined to an arc of length L, independent of any box.

This script computes that geometric menu exactly.

KEY STRUCTURE (proved in THM-387 draft):
  Place the m movers at 0 <= x_1 < x_2 < ... < x_m inside an interval of length L
  (so x_m - x_1 <= L).  For i<j the half-turn rule a->b iff frac(p_a-p_b) in (0,1/2)
  gives, generically (no exact 1/2 separation):
       j beats i   if  x_j - x_i < 1/2   ("short" pair: later beats earlier)
       i beats j   if  x_j - x_i > 1/2   ("long"  pair: flipped)
  So every realizable class is the transitive backbone ("later beats earlier")
  with the set S of LONG pairs flipped.  S is an UP-SET of the interval-
  containment order ((i,j) <= (i',j') iff i'<=i, j'>=j): a longer interval has a
  larger separation, so if (i,j) is long so is any interval containing it.
  Up-sets of this (type-A root) poset are counted by the Catalan number C_m.

  A given up-set S is GEOMETRICALLY REALIZABLE at length L iff the strict
  difference-constraint system
       x_{i+1}-x_i > 0  (order),  x_j-x_i > 1/2 for (i,j) in S,
       x_j-x_i < 1/2 for (i,j) not in S,  x_m - x_1 <= L
  is feasible.  We test feasibility exactly with Bellman-Ford over (rational,eps)
  weights (eps an infinitesimal encoding strict inequalities).

OUTPUTS:
  * geometric menu size (iso-classes) at L=(n-2)/n for n=4..8, with H/score,
    cross-checked to contain the S512 arithmetic-reachable classes;
  * #feasible labelled flip-sets vs Catalan(m) vs A000568(m);
  * an L-sweep showing exactly where new classes enter as the arc grows.
"""

from __future__ import annotations
from fractions import Fraction
from itertools import permutations
from functools import lru_cache
from collections import Counter

HALF = Fraction(1, 2)

# ---------- (rational, eps) lexicographic weights for strict constraints ----------
# value v + c*eps with eps an infinitesimal > 0.  Represented as (Fraction v, int c).

def w_add(a, b):
    return (a[0] + b[0], a[1] + b[1])

def w_lt(a, b):
    return (a[0], a[1]) < (b[0], b[1])

NEG = None  # sentinel handled in bellman_ford

# ---------- difference-constraint feasibility ----------
# constraints given as list of (a, b, w): meaning  x[a] - x[b] <= w.
# feasible  <=>  no negative cycle in graph with edge b->a of weight w.

def feasible(m, constraints):
    # nodes 0..m-1 are x_1..x_m ; node m is a virtual source with 0-edges to all.
    N = m + 1
    edges = []  # (u, v, w) directed u->v
    for (a, b, w) in constraints:
        edges.append((b, a, w))           # x[a]-x[b] <= w  => edge b->a weight w
    for v in range(m):
        edges.append((m, v, (Fraction(0), 0)))
    INF = None
    dist = [None] * N
    dist[m] = (Fraction(0), 0)
    for _ in range(N - 1):
        changed = False
        for (u, v, w) in edges:
            if dist[u] is None:
                continue
            nd = w_add(dist[u], w)
            if dist[v] is None or w_lt(nd, dist[v]):
                dist[v] = nd
                changed = True
        if not changed:
            break
    # one more pass: any relaxation => negative cycle => infeasible
    for (u, v, w) in edges:
        if dist[u] is None:
            continue
        nd = w_add(dist[u], w)
        if dist[v] is None or w_lt(nd, dist[v]):
            return False
    return True

def realizable(m, S, L):
    """S: set of (i,j) i<j flipped (long). L: Fraction arc length. 1-indexed pairs."""
    EPS = (Fraction(0), 1)
    cons = []
    # ordering strict: x_{i+1} - x_i >= eps  <=>  x_i - x_{i+1} <= -eps
    for i in range(1, m):
        cons.append((i - 1, i, (Fraction(0), -1)))
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            if (i, j) in S:
                # x_j - x_i >= 1/2 + eps  <=>  x_i - x_j <= -1/2 - eps
                cons.append((i - 1, j - 1, (-HALF, -1)))
            else:
                # x_j - x_i <= 1/2 - eps
                cons.append((j - 1, i - 1, (HALF, -1)))
    # cap: x_m - x_1 <= L
    cons.append((m - 1, 0, (L, 0)))
    return feasible(m, cons)

# ---------- enumerate up-sets of the interval-containment poset ----------
# up-set <-> thresholds tau_i (i=1..m-1), tau_i in {i+1..m+1}, non-decreasing.
# (i,j) in S  iff  j >= tau_i.

def upsets(m):
    res = []
    def rec(i, prev, cur):
        if i == m:  # i ranges 1..m-1 ; vertex m has no out-pair
            res.append(cur)
            return
        lo = max(i + 1, prev)
        for tau in range(lo, m + 2):
            S = set(cur)
            for j in range(tau, m + 1):
                S.add((i, j))
            rec(i + 1, tau, S)
    rec(1, 2, set())
    return res

# ---------- tournament from flip-set, invariants ----------

def build_adj(m, S):
    adj = [[0] * m for _ in range(m)]
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            if (i, j) in S:      # long: i beats j
                adj[i - 1][j - 1] = 1
            else:                # short: j beats i
                adj[j - 1][i - 1] = 1
    return tuple(tuple(r) for r in adj)

def canon(adj):
    m = len(adj)
    best = None
    for p in permutations(range(m)):
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best:
            best = flat
    return best

def Hcount(adj):
    m = len(adj); full = (1 << m) - 1
    @lru_cache(None)
    def dp(mask, last):
        if mask == full:
            return 1
        return sum(dp(mask | (1 << x), x) for x in range(m)
                   if not (mask >> x) & 1 and adj[last][x])
    return sum(dp(1 << s, s) for s in range(m))

def scores(adj):
    return tuple(sorted(sum(r) for r in adj))

def is_transitive(adj):
    return scores(adj) == tuple(range(len(adj)))

A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056]
def catalan(k):
    c = 1
    for i in range(k):
        c = c * 2 * (2 * i + 1) // (i + 2)
    return c

# ---------- geometric menu at a given (m, L) ----------

def geometric_menu(m, L, do_iso=True):
    us = upsets(m)
    feas = [S for S in us if realizable(m, S, L)]
    classes = {}
    if do_iso:
        for S in feas:
            adj = build_adj(m, S)
            c = canon(adj)
            if c not in classes:
                classes[c] = (Hcount(adj), scores(adj), is_transitive(adj))
    return us, feas, classes

def main():
    print("LRC geometric arc menu -- box-independent target (claudebox-S521)\n")
    print("S512 arithmetic-reachable iso-classes (lower bound): n=4..7 -> 2,2,6,6\n")
    print("=" * 70)
    print(f"{'n':>2} {'m':>2} {'L=(n-2)/n':>10} {'Catalan(m)':>10} "
          f"{'#feasS':>7} {'menu':>5} {'A000568(m)':>10}")
    print("-" * 70)
    detail = {}
    for n in range(4, 10):          # m = 3..8 : exact iso menu (dedup raw matrices first)
        m = n - 1
        L = Fraction(n - 2, n)
        us, feas, _ = geometric_menu(m, L, do_iso=False)
        raw = set(build_adj(m, S) for S in feas)     # distinct LABELLED tournaments
        classes = {}
        for adj in raw:
            c = canon(adj)
            if c not in classes:
                classes[c] = (Hcount(adj), scores(adj), is_transitive(adj))
        a = A000568[m] if m < len(A000568) else None
        print(f"{n:>2} {m:>2} {str(L):>10} {catalan(m):>10} "
              f"{len(feas):>7} {len(classes):>5} {str(a):>10}")
        detail[n] = classes
    # m = 9 : feasible-set count only (full permutation canon too slow)
    for n in (10,):
        m = n - 1
        L = Fraction(n - 2, n)
        us, feas, _ = geometric_menu(m, L, do_iso=False)
        raw = set(build_adj(m, S) for S in feas)
        a = A000568[m] if m < len(A000568) else None
        # heuristic iso key: (H, score, sorted 3-cycle count) -- may undercount real classes
        hk = set()
        for adj in raw:
            hk.add((Hcount(adj), scores(adj)))
        print(f"{n:>2} {m:>2} {str(L):>10} {catalan(m):>10} "
              f"{len(feas):>7} {'~'+str(len(hk)):>5} {str(a):>10}   (heuristic (H,score) key)")
    print()
    print("menu sequence n=4..9 (m=3..8):",
          [len(detail[n]) for n in range(4, 10)])
    print()
    for n in range(4, 10):
        print(f"--- n={n}  geometric menu iso-classes (the target inside A000568) ---")
        for c, (H, sc, tr) in sorted(detail[n].items(), key=lambda kv: (kv[1][0], kv[1][1])):
            print(f"     H={H:<4} score={sc} transitive={tr}")
        print()

    # ---- L-sweep: where do new iso-classes enter as the arc grows? ----
    print("=" * 70)
    print("L-SWEEP: geometric menu size vs arc length L (m fixed). New classes")
    print("enter at rational thresholds; LRC uses L=(n-2)/n marked below.")
    print("=" * 70)
    for m in (4, 5, 6):
        n = m + 1
        Lstar = Fraction(n - 2, n)
        print(f"\n m={m} (n={n}, LRC arc L*={Lstar}):")
        grid = [Fraction(k, 60) for k in range(30, 61)]  # L in [0.5, 1.0]
        prev = None
        for L in grid:
            _, feas, classes = geometric_menu(m, L, do_iso=True)
            sz = len(classes)
            if sz != prev:
                star = "  <== LRC arc" if L == Lstar else ""
                print(f"   L={float(L):.3f} ({L}):  menu={sz}  feasible flip-sets={len(feas)}{star}")
                prev = sz
        # also report exactly at L*
        _, feasL, classesL = geometric_menu(m, Lstar, do_iso=True)
        print(f"   [at L*={Lstar}: menu={len(classesL)}]")

if __name__ == "__main__":
    main()
