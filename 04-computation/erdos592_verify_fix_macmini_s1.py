#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, part 6 — FIXED exhaustive subgrid verifier + full rerun of Q(n,t)
and invariant Q(n,t).   mac-mini-2026-06-09-S1

BUG FOUND (would-be MISTAKE-066): find_independent_binary_subgrid in part 4 committed
greedily to the FIRST consistent subtree under each chosen child and never explored
alternative subtrees for the same child; an incomplete search can return "no independent
subgrid" when one exists, invalidating CEGAR's final SAT certificates. (Spotted because
the t=11 invariant witness visibly left rectangles {2,3}x{4,5} unhit at gap 1.)

This file: exhaustive generator-based search (sound and complete), then re-run of all
headline computations. Free-graph UNSATs and constructive grid-existence results from
parts 1-2 are unaffected (their positive witnesses are sound); every SAT certificate
from part 4/5 is re-derived here.
"""

import itertools, time
from pysat.solvers import Glucose3

def leaves(n, t):
    return list(itertools.product(range(t), repeat=n))

def find_independent_binary_subgrid_exhaustive(n, t, adj, L, idx):
    """Complete backtracking over ALL binary subgrids; returns one independent subgrid
    (list of 2^n leaves) or None.  chosen = current leaf stack."""
    chosen = []

    def gen(prefix):
        """Yield once for every completion of the binary subtree below prefix that is
        pairwise-independent with everything in `chosen` (state pushed on chosen)."""
        k = len(prefix)
        if k == n:
            i = idx[prefix]
            for c in chosen:
                if adj[i][idx[c]]:
                    return
            chosen.append(prefix)
            yield
            chosen.pop()
            return
        for c1 in range(t):
            for _ in gen(prefix + (c1,)):
                for c2 in range(c1 + 1, t):
                    yield from gen(prefix + (c2,))

    for _ in gen(()):
        return list(chosen)
    return None

def solve_Q(n, t, invariant=False, verbose=True, tlimit=600):
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    qvar = {}
    cnt = [0]

    def q(x, y):
        if x > y:
            x, y = y, x
        if invariant:
            a, a2 = x[0], y[0]
            key = (a2 - a, x[1:], y[1:]) if a2 != a else (0,) + tuple(sorted((x[1:], y[1:])))
        else:
            key = (x, y)
        if key not in qvar:
            cnt[0] += 1
            qvar[key] = cnt[0]
        return qvar[key]

    sol = Glucose3()
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if c not in seen:
            seen.add(c)
            sol.add_clause(list(c))

    rounds, added = 0, 0
    t0 = time.time()
    while True:
        rounds += 1
        if time.time() - t0 > tlimit:
            if verbose:
                print(f"   TIMEOUT after {added} clauses, {time.time()-t0:.0f}s")
            return None, None
        if not sol.solve():
            if verbose:
                print(f"   UNSAT ({'invariant ' if invariant else ''}exact; {added} lazy "
                      f"clauses, {time.time()-t0:.1f}s)")
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        adj = [[False] * N for _ in range(N)]
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    adj[i][j] = adj[j][i] = True
                    edges.add((i, j))
        bad = find_independent_binary_subgrid_exhaustive(n, t, adj, L, idx)
        if bad is None:
            if verbose:
                print(f"   SAT (verified by COMPLETE search): {len(edges)} edges "
                      f"({added} lazy clauses, {time.time()-t0:.1f}s)")
            return True, (edges, L, model, qvar)
        sol.add_clause(sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2))))
        added += 1

def main():
    print("=" * 78)
    print("RERUN with complete verifier — Q(n,t) free and invariant")
    print("=" * 78)
    results = {}
    print("\n--- free n=2 ---")
    for t in range(2, 15):
        print(f"Q(2,{t}):", flush=True)
        res, wit = solve_Q(2, t)
        results[("free", 2, t)] = res
        if res is False:
            print(f"   *** R(2,2) = {t} (exact, complete verifier) ***")
            break
        if res is None:
            break
    print("\n--- free n=3 ---")
    for t in range(2, 7):
        print(f"Q(3,{t}):", flush=True)
        res, wit = solve_Q(3, t, tlimit=900)
        results[("free", 3, t)] = res
        if res is False:
            print(f"   *** R(3,2) = {t} (exact) ***")
            break
        if res is None:
            break
    print("\n--- invariant n=2 ---")
    for t in range(2, 15):
        print(f"invQ(2,{t}):", flush=True)
        res, wit = solve_Q(2, t, invariant=True)
        results[("inv", 2, t)] = res
        if res is False:
            print(f"   *** invariant cutoff at t = {t} ***")
            break
        if res is None:
            break
    print("\n--- invariant n=3 ---")
    for t in range(2, 7):
        print(f"invQ(3,{t}):", flush=True)
        res, wit = solve_Q(3, t, invariant=True, tlimit=900)
        results[("inv", 3, t)] = res
        if res is False:
            print(f"   *** invariant n=3 cutoff at t = {t} ***")
            break
        if res is None:
            break
    print("\nSUMMARY:", results)

if __name__ == "__main__":
    main()
