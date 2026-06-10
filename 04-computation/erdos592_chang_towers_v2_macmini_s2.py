#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, session 2 part 2b — Chang tower sweep, FAST verifier.
mac-mini-2026-06-09-S2  (THM-460 D; supersedes the v1 sweep loop, same definitions)

Speedups over v1 (same complete semantics, cross-validated again):
  * all binary 2-grids precomputed once per ambient (4-subsets passing the
    THM-460 B3(ii) split-hierarchy check);
  * per CEGAR round: adjacency bitmasks; a tower exists iff SOME independent
    2-grid Γ has ≥ 2 below-elements non-adjacent to Γ with an independent pair
    among them.  Since the model is triangle-free (clauses upfront), any 3 such
    candidates contain an independent pair — only the count-2 case needs a check.
  * M=3 variant: an independent (pair < 2-grid < 3-grid) stack; 3-grids
    precomputed structurally (8-subsets are too many: built from compatible
    2-grid pairs).
"""

import itertools, time
from pysat.solvers import Glucose3
from erdos592_chang_towers_macmini_s2 import ambient, split, is_binary_grid, brute_independent_tower

def all_2grids(L):
    """All binary 2-grids as index 4-tuples (sorted)."""
    N = len(L)
    out = []
    for c in itertools.combinations(range(N), 4):
        if is_binary_grid(tuple(L[i] for i in c)):
            out.append(c)
    return out

def all_3grids(L, grids2):
    """Binary 3-grids = two 2-grids A < B, cross-split common & above internal."""
    N = len(L)
    out = []
    by_min = {}
    for g in grids2:
        by_min.setdefault(g[0], []).append(g)
    for A in grids2:
        for mn in range(A[-1] + 1, N):
            for B in by_min.get(mn, []):
                S = A + B
                if is_binary_grid(tuple(L[i] for i in S)):
                    out.append(S)
    return out

def solve_chang2(M, s, C, tlimit=2400, verbose=True):
    L = ambient(s, C)
    N = len(L)
    grids2 = all_2grids(L)
    grids3 = all_3grids(L, grids2) if M >= 3 else []
    if verbose:
        print(f"   [ambient {N} vertices, {len(grids2)} 2-grids"
              + (f", {len(grids3)} 3-grids]" if M >= 3 else "]"), flush=True)
    evar = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            evar[(i, j)] = cnt
    sol = Glucose3()
    for i, j, k in itertools.combinations(range(N), 3):
        sol.add_clause([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])

    def epair(i, j):
        return evar[(min(i, j), max(i, j))]

    t0 = time.time()
    added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT M={M} s={s} C={C} ({added} cls, {time.time()-t0:.0f}s)", flush=True)
            return None
        if not sol.solve():
            if verbose: print(f"   UNSAT  M={M} s={s} C={C}  (lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        nb = [0] * N            # neighbor bitmask
        for (i, j), v in evar.items():
            if v in model:
                nb[i] |= 1 << j
                nb[j] |= 1 << i
        full = (1 << N) - 1

        def indep(ids):
            return all(not (nb[a] >> b) & 1 for a, b in itertools.combinations(ids, 2))

        def below_cands(g):
            mask = (1 << g[0]) - 1   # indices < min(g)  (L sorted = index order)
            for x in g:
                mask &= ~nb[x]
            for x in g:
                mask &= ~(1 << x)
            return mask

        def pick_pair(mask):
            ids = []
            m = mask
            while m:
                b = m & -m
                ids.append(b.bit_length() - 1)
                m ^= b
            if len(ids) < 2:
                return None
            for a, b in itertools.combinations(ids, 2):
                if not (nb[a] >> b) & 1:
                    return (a, b)
            return None  # all pairs adjacent: only possible if len(ids)==2

        bad = None
        if M == 2:
            for g in grids2:
                if not indep(g):
                    continue
                pr = pick_pair(below_cands(g))
                if pr:
                    bad = list(pr) + list(g)
                    break
        else:  # M == 3: pair < 2-grid < 3-grid, jointly independent
            for g3 in grids3:
                if not indep(g3):
                    continue
                m3 = below_cands(g3)
                # find an independent 2-grid among candidates below, then a pair below that
                for g2 in grids2:
                    if g2[-1] >= g3[0]:
                        continue
                    if any(not (m3 >> x) & 1 for x in g2):
                        continue
                    if not indep(g2):
                        continue
                    mask = below_cands(g2) & m3
                    pr = pick_pair(mask)
                    if pr:
                        bad = list(pr) + list(g2) + list(g3)
                        break
                if bad:
                    break
        if bad is None:
            edges = sum(1 for v in evar.values() if v in model)
            if verbose: print(f"   SAT    M={M} s={s} C={C}  ({edges} edges, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return True
        sol.add_clause([epair(a, b) for a, b in itertools.combinations(bad, 2)])
        added += 1

def crossvalidate():
    import random
    rng = random.Random(5)
    print("=== v2 cross-validation vs brute (M=2) ===", flush=True)
    for (s, C) in ((2, 3), (3, 2), (2, 4)):
        L = ambient(s, C); N = len(L)
        idx = {v: i for i, v in enumerate(L)}
        grids2 = all_2grids(L)
        bad = 0
        for tr in range(30):
            adj = [[False] * N for _ in range(N)]
            nb = [0] * N
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < 0.35 * rng.random():
                        adj[i][j] = adj[j][i] = True
                        nb[i] |= 1 << j; nb[j] |= 1 << i
            # v2-style search
            found = None
            for g in grids2:
                if any((nb[a] >> b) & 1 for a, b in itertools.combinations(g, 2)):
                    continue
                mask = (1 << g[0]) - 1
                for x in g:
                    mask &= ~nb[x]; mask &= ~(1 << x)
                ids = [b for b in range(N) if (mask >> b) & 1]
                pr = None
                for a, b in itertools.combinations(ids, 2):
                    if not adj[a][b]:
                        pr = (a, b); break
                if pr:
                    found = pr + g
                    break
            r2 = brute_independent_tower(s, C, 2, L, idx, adj)
            if (found is None) != (r2 is None):
                bad += 1
                print(f"   DISAGREE s={s} C={C} trial {tr}", flush=True)
        print(f"   s={s} C={C}: 30 trials {'OK' if bad == 0 else str(bad) + ' BAD'}", flush=True)

def main():
    crossvalidate()
    print("\n=== M=2 Chang sweep (fast verifier) ===", flush=True)
    first_unsat = None
    for (s, C) in ((2, 3), (3, 2), (2, 4), (3, 3), (2, 5), (4, 2), (2, 6), (3, 4),
                   (2, 7), (4, 3), (2, 8), (5, 2), (3, 5)):
        if C ** s > 130:
            print(f"   skip s={s} C={C} (size {C**s})", flush=True)
            continue
        r = solve_chang2(2, s, C)
        if r is False and first_unsat is None:
            first_unsat = (s, C)
    print(f"\nfirst UNSAT in sweep: {first_unsat}", flush=True)
    print("\n=== M=3 spot checks ===", flush=True)
    for (s, C) in ((3, 3), (2, 5), (4, 2), (2, 6)):
        solve_chang2(3, s, C, tlimit=1800)

if __name__ == "__main__":
    main()
