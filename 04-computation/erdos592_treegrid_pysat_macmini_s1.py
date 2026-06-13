#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, part 4 — Q(n,t) decided exactly with a real SAT solver (pysat/Glucose3).
mac-mini-2026-06-09-S1  (T768, HYP-2346, THM-453)

Q(n,t): does a triangle-free graph exist on the leaves [t]^n of the complete t-ary
tree of height n hitting EVERY binary subgrid (= no independent binary subgrid)?

THE BRIDGE (corrected, compactness-tight — see draft doc §3.2):
  R(n,2) := least T such that EVERY triangle-free graph on the T-grid has an
            independent binary subgrid  ( = least T with Q(n,T) UNSAT ).
  * If Q(n,T) is SAT for ALL T, then by König/diagonalization there is a triangle-free
    graph on the full omega-grid with no independent binary subgrid; since every
    subset of type omega^n contains a full subgrid (hence a binary one), that graph
    witnesses omega^n -/-> (omega^n, 3).
  * Hence omega^2 -> (omega^2,3) (Specker) FORCES R(2,2) < infinity: Q(2,t) must die.
  * Q(3,t) SAT for all t (+ a uniform rule extracted from witnesses) would re-prove
    Specker's negative omega^3 -/-> (omega^3,3) with an explicit strong witness.
  Calibration: R(1,2) = 3 (Ramsey shadow; verified by DPLL in part 3).

CEGAR: subgrid clauses added lazily (find an independent binary subgrid in the current
model by backtracking; add its clause; repeat) — scales past full clause generation.
"""

import itertools, sys, time
from pysat.solvers import Glucose3

def leaves(n, t):
    return list(itertools.product(range(t), repeat=n))

def find_independent_binary_subgrid(n, t, adj, L, idx):
    """Backtracking: find a binary subgrid (2^n leaves) that is independent in adj.
    Returns the leaf tuple list or None.  Subgrid = choose 2 children per node."""
    # Represent partial subgrid construction level by level (DFS over the binary tree).
    # We build the set of chosen leaves; constraint: pairwise non-adjacent.
    chosen = []

    def build(prefix):
        """Complete the binary subtree below 'prefix' (a tuple in [t]^k)."""
        k = len(prefix)
        if k == n:
            leaf = prefix
            for c in chosen:
                if adj[idx[leaf]][idx[c]]:
                    return False
            chosen.append(leaf)
            return True
        for c1 in range(t):
            ok1 = build(prefix + (c1,))
            if not ok1:
                continue
            n1 = 2 ** (n - k - 1)
            for c2 in range(c1 + 1, t):
                if build(prefix + (c2,)):
                    return True
                # keep c1 subtree, try next c2
            del chosen[-n1:]
        return False

    if build(()):
        return list(chosen)
    return None

def solve_Q(n, t, verbose=True, max_rounds=100000):
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    evar = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            evar[(i, j)] = cnt

    sol = Glucose3()
    ntri = 0
    for i, j, k in itertools.combinations(range(N), 3):
        sol.add_clause([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])
        ntri += 1

    rounds = 0
    added = 0
    t0 = time.time()
    while True:
        rounds += 1
        if rounds > max_rounds:
            return None, None, rounds, added
        if not sol.solve():
            if verbose:
                print(f"   UNSAT after {added} subgrid clauses, {rounds} rounds, "
                      f"{time.time()-t0:.1f}s")
            return False, None, rounds, added
        model = set(l for l in sol.get_model() if l > 0)
        adj = [[False] * N for _ in range(N)]
        edges = set()
        for (i, j), v in evar.items():
            if v in model:
                adj[i][j] = adj[j][i] = True
                edges.add((i, j))
        bad = find_independent_binary_subgrid(n, t, adj, L, idx)
        if bad is None:
            if verbose:
                print(f"   SAT: witness with {len(edges)} edges "
                      f"({added} lazy clauses, {rounds} rounds, {time.time()-t0:.1f}s)")
            return True, (edges, L), rounds, added
        cl = []
        for a, b in itertools.combinations(bad, 2):
            i, j = idx[a], idx[b]
            if i > j: i, j = j, i
            cl.append(evar[(i, j)])
        sol.add_clause(cl)
        added += 1

def meet_level(a, b):
    m = 0
    while m < len(a) and a[m] == b[m]:
        m += 1
    return m

def analyse(edges, L, n, t, label):
    from collections import Counter
    hist = Counter()
    detail = Counter()
    for (i, j) in edges:
        a, b = L[i], L[j]
        m = meet_level(a, b)
        hist[m] += 1
        if a > b:
            a, b = b, a
        # coarse within-split feature: child gap at split + suffix comparison profile
        if m < n:
            gap = b[m] - a[m]
            sufcmp = tuple(("<" if a[k] < b[k] else (">" if a[k] > b[k] else "="))
                           for k in range(m + 1, n))
            detail[(m, gap, sufcmp)] += 1
    print(f"   [{label}] meet-level histogram: {dict(sorted(hist.items()))}")
    top = detail.most_common(12)
    print(f"   [{label}] top (meetlevel, childgap, suffix-cmp) classes:")
    for cls, c in top:
        print(f"        {cls}: {c}")

def main():
    print("=" * 78)
    print("Q(n,t) decided with Glucose3 + CEGAR  (R(n,2) hunt)")
    print("=" * 78)

    print("\n--- n=2: Specker positive forces a finite cutoff R(2,2) ---")
    for t in range(2, 9):
        print(f"n=2 t={t}: leaves={t*t}", flush=True)
        res, wit, rounds, added = solve_Q(2, t)
        if res is True and wit:
            analyse(wit[0], wit[1], 2, t, f"Q(2,{t}) witness")
        if res is False:
            print(f"   *** R(2,2) = {t} : every triangle-free graph on the {t}-grid "
                  f"has an independent binary subgrid ***")
            break

    print("\n--- n=3: hunting strong-witness shadows (expect SAT to persist) ---")
    for t in range(2, 6):
        print(f"n=3 t={t}: leaves={t**3}", flush=True)
        res, wit, rounds, added = solve_Q(3, t)
        if res is True and wit:
            analyse(wit[0], wit[1], 3, t, f"Q(3,{t}) witness")
            if t <= 3:
                print("   full witness edge list:")
                for (i, j) in sorted(wit[0]):
                    print("      ", wit[1][i], "--", wit[1][j])
        if res is False:
            print(f"   *** Q(3,{t}) UNSAT — R(3,2) = {t} (finite!) — would be evidence "
                  f"AGAINST a strong omega^3 witness ***")
            break

if __name__ == "__main__":
    main()
