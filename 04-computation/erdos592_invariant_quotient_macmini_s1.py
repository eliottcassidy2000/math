#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, part 5 — TRANSLATION-INVARIANT Q(n,t): the graded-relation quotient.
mac-mini-2026-06-09-S1  (T768, HYP-2346, THM-453 part F)

Same question as Q(n,t) (triangle-free graph on [t]^n leaves hitting every binary
subgrid) but the graph is forced TRANSLATION-INVARIANT in the top-level coordinate:
the edge variable for {x,y} depends only on (a'-a, (suffix of x), (suffix of y)) where
a,a' are the first coordinates — i.e. exactly the data (R, {B_g}) of THM-453 part F1:
  within-row graph R on [t]^{n-1} (a'=a), and gap-graded cross relations
  B_g ⊆ [t]^{n-1} × [t]^{n-1} (g = a'-a >= 1).
Triangle-freeness on the grid then encodes the graded composition-freeness
  B_{g1} ∘ B_{g2} ∩ B_{g1+g2} = ∅   (+ mixed R conditions, + R triangle-free),
and subgrid-hitting is the rectangle/subgrid density floor.

PREDICTION (THM-453 F2/F3): invariant Q(2,t) dies EARLY (additive Schur-type
obstruction: density vs composition-freeness over (N,+)); invariant Q(3,t) survives
longer (columns are hierarchical, density demand drops to column-subgrids only).

CEGAR as in part 4; exact (Glucose3).
"""

import itertools, sys, time
from pysat.solvers import Glucose3
from erdos592_treegrid_pysat_macmini_s1 import leaves, find_independent_binary_subgrid

def solve_Q_invariant(n, t, verbose=True):
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)

    # quotient variable for a pair {x,y}
    qvar = {}
    cnt = [0]
    def q(x, y):
        if x > y:
            x, y = y, x
        a, a2 = x[0], y[0]
        key = (a2 - a, x[1:], y[1:]) if a2 != a else (0,) + tuple(sorted((x[1:], y[1:])))
        if key not in qvar:
            cnt[0] += 1
            qvar[key] = cnt[0]
        return qvar[key]

    sol = Glucose3()
    # triangle clauses on the grid, mapped through the quotient (dedup)
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k]))))
        if len(set(map(abs, c))) < 3:
            # a triangle whose two sides map to the SAME quotient var: the clause
            # (-v or -v or -w) = (-v or -w); keep reduced clause
            c = tuple(sorted(set(c)))
        if c not in seen:
            seen.add(c)
            sol.add_clause(list(c))
    ntri = len(seen)

    rounds, added = 0, 0
    t0 = time.time()
    while True:
        rounds += 1
        if not sol.solve():
            if verbose:
                print(f"   INVARIANT-UNSAT  (quotient vars={cnt[0]}, tri-clauses={ntri}, "
                      f"{added} lazy subgrid clauses, {time.time()-t0:.1f}s)")
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        adj = [[False] * N for _ in range(N)]
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    adj[i][j] = adj[j][i] = True
                    edges.add((i, j))
        bad = find_independent_binary_subgrid(n, t, adj, L, idx)
        if bad is None:
            if verbose:
                print(f"   INVARIANT-SAT: witness with {len(edges)} grid edges "
                      f"(quotient vars={cnt[0]}, {added} lazy clauses, {time.time()-t0:.1f}s)")
            return True, (edges, L, model, qvar)
        cl = sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1

def describe_invariant_witness(model, qvarmap, n, t):
    """Print the (R, B_g) structure of an invariant witness."""
    from collections import defaultdict
    Bg = defaultdict(list)
    R = []
    for key, v in qvarmap.items():
        if v in model:
            if key[0] == 0:
                R.append((key[1], key[2]))
            else:
                Bg[key[0]].append((key[1], key[2]))
    print(f"      R (within-row): {len(R)} edges:", sorted(R)[:14])
    for g in sorted(Bg):
        print(f"      B_{g}: {len(Bg[g])} pairs:", sorted(Bg[g])[:12],
              "..." if len(Bg[g]) > 12 else "")

def main():
    print("=" * 78)
    print("TRANSLATION-INVARIANT Q(n,t): the (R,{B_g}) graded-relation quotient")
    print("=" * 78)
    print("\n--- n=2 invariant (prediction: dies early — Schur-type obstruction) ---")
    for t in range(2, 13):
        print(f"inv n=2 t={t}:", flush=True)
        res, wit = solve_Q_invariant(2, t)
        if res and wit:
            edges, L, model, qm = wit
            describe_invariant_witness(model, qm, 2, t)
        if res is False:
            print(f"   *** invariant R(2,2) = {t} ***")
            break

    print("\n--- n=3 invariant (prediction: survives — hierarchical columns) ---")
    for t in range(2, 7):
        print(f"inv n=3 t={t}:", flush=True)
        res, wit = solve_Q_invariant(3, t)
        if res and wit:
            edges, L, model, qm = wit
            if t <= 4:
                describe_invariant_witness(model, qm, 3, t)
        if res is False:
            print(f"   *** invariant Q(3,{t}) UNSAT at t={t} ***")
            break

if __name__ == "__main__":
    main()
