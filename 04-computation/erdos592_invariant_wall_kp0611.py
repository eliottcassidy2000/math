#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-470 lab, part 2: THE t=7 WALL (kind-pasteur-2026-06-11-S1, HYP-2396).

Trigger: the dyadic 1-jet algebra F2J is feature-UNSAT PER-SIZE at (3,7) — and
F2J refines F2, so the bi-dyadic algebra also dies per-size at t=7 (one size past
THM-465 A's frontier).  Gap-determined feature algebras are ANTITONE in t (a
witness at t+1 restricts to one at t: same classes, fewer demands), so each such
algebra's SAT region is an interval [3, t_dead).

Master question: the FINEST gap-determined algebra — the full translation-
invariant one, Finv(x,y) = the gap vector itself — at (3,7).
  * If feature-UNSAT: NO translation-invariant strong witness exists on [7]^3,
    subsuming every gap-determined rung (F2, F2J, F2X, leading-digit, ...) at all
    t >= 7; together with R(1,2)=3 and R(2,2)=5 (where invariant = free, THM-453 G)
    this is the first support for R(n,2) = 2n+1 — which would REFUTE HYP-2363
    (strong-witness frontier) and redirect the constructive-Specker program toward
    non-invariant or non-strong witnesses.
  * If SAT: the wall is algebra-specific; the invariant ladder continues above the
    jet rung and the table is a new witness one size past the known frontier.

Control first: Finv at (3,6) must be SAT (F2 SAT there already implies it).
"""

import itertools, time
from pysat.solvers import Glucose3
from erdos592_satverifier_frontier_macmini_s2 import leaves, SubgridVerifier, seed_subgrids
from erdos592_bidyadic_rule_macmini_s3 import solve_featured, independent_verify


def Finv(x, y):
    return tuple(y[i] - x[i] for i in range(len(x)))


def solve_featured_verbose(n, t, F, name, tlimit=900):
    """solve_featured clone with per-round instrumentation (progress every 100
    CEGAR rounds) — same encoding, same complete verifier."""
    L = leaves(n, t); N = len(L)
    qv = {}; cnt = [0]

    def q(x, y):
        if x > y:
            x, y = y, x
        key = F(x, y)
        if key not in qv:
            cnt[0] += 1; qv[key] = cnt[0]
        return qv[key]

    sol = Glucose3()
    seen = set()
    tb = time.time()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if c not in seen:
            seen.add(c); sol.add_clause(list(c))
    ntri = len(seen)
    for g in seed_subgrids(n, t):
        c = tuple(sorted(set(q(a, b) for a, b in itertools.combinations(g, 2))))
        if c not in seen:
            seen.add(c); sol.add_clause(list(c))
    print(f"   [{name} ({n},{t}): {cnt[0]} features, {ntri} tri + {len(seen)-ntri} seed clauses, "
          f"build {time.time()-tb:.1f}s]", flush=True)
    ver = SubgridVerifier(n, t)
    t0 = time.time(); added = 0
    while True:
        if time.time() - t0 > tlimit:
            print(f"   TIMEOUT {name} ({n},{t}) lazy={added}", flush=True)
            return None, None
        if not sol.solve():
            print(f"   feature-UNSAT {name} ({n},{t}) (lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    edges.add((i, j))
        bad = ver.find(edges)
        if bad is None:
            print(f"   SAT {name} ({n},{t}) ({len(edges)} edges, lazy={added}, "
                  f"{time.time()-t0:.1f}s)", flush=True)
            return True, (edges, L, model, qv)
        idx = {v: i for i, v in enumerate(L)}
        cl = sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1
        if added % 100 == 0:
            print(f"      ...{name} ({n},{t}) round {added}, {time.time()-t0:.1f}s", flush=True)


def main():
    t0 = time.time()
    print("=== plumbing control: invariant algebra at (3,4) — mac-mini S2 had "
          "invQ(3,4) SAT (317 edges, different code) ===", flush=True)
    r, w = solve_featured_verbose(3, 4, Finv, "INV", tlimit=1200)
    if r:
        ok = independent_verify(3, 4, w[0])
        print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)

    print("\n=== THE WALL: full invariant algebra at (3,7) ===", flush=True)
    r7, w7 = solve_featured_verbose(3, 7, Finv, "INV", tlimit=7200)
    if r7:
        ok = independent_verify(3, 7, w7[0])
        print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
        print("   -> wall is algebra-specific; invariant ladder lives at t=7", flush=True)
    elif r7 is False:
        print("   -> NO translation-invariant strong witness on [7]^3;", flush=True)
        print("      all gap-determined algebras dead for ALL t >= 7 (antitone);", flush=True)
        print("      R(n,2) = 2n+1 pattern: 3, 5, 7?  (HYP-2396)", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
