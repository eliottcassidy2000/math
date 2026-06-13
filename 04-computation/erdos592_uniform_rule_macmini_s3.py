#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 session 3, part 3 — A SINGLE t-UNIFORM BI-DYADIC RULE (THM-465 core).
mac-mini-2026-06-10-S1

(a) Freeze the 42-class table found at t=5; build its rule graph at t = 6, 7, 8 and
    verify (triangle-free + dominates all binary subgrids) with the complete verifier.
(b) SIMULTANEOUS-t SAT: one shared set of feature variables; triangle + subgrid
    constraints from t = 4,5,6,7 conjoined.  A model = ONE table valid at all four
    sizes — the uniform-rule candidate for the constructive strong Specker
    (THM-453 D1 route).  Verify the found table independently at every t including
    one size beyond (t=8).

Note: for t ≤ 8 all gaps are ≤ 7, so v2 ≤ 2 and the F2 feature space is the SAME
finite set across these sizes — a fixed table is well-defined.  (v2 = 3 first
appears at gap 8, i.e. t = 9; capping v2 at 2 extends the rule format beyond.)
"""

import itertools, time
from pysat.solvers import Glucose3
from erdos592_satverifier_frontier_macmini_s2 import SubgridVerifier, leaves, seed_subgrids
from erdos592_bidyadic_rule_macmini_s3 import F2, independent_verify

def rule_graph(n, t, blue_classes):
    L = leaves(n, t)
    edges = set()
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if F2(L[i], L[j]) in blue_classes:
                edges.add((i, j))
    return edges, L

def check_rule(n, t, blue_classes, label):
    edges, L = rule_graph(n, t, blue_classes)
    N = len(L)
    adj = [[False] * N for _ in range(N)]
    for (i, j) in edges:
        adj[i][j] = adj[j][i] = True
    for i, j, k in itertools.combinations(range(N), 3):
        if adj[i][j] and adj[i][k] and adj[j][k]:
            print(f"   [{label}] t={t}: TRIANGLE — rule fails triangle-freeness", flush=True)
            return False
    ver = SubgridVerifier(n, t)
    bad = ver.find(edges)
    if bad is None:
        print(f"   [{label}] t={t}: VERIFIED witness ({len(edges)} edges)", flush=True)
        return True
    print(f"   [{label}] t={t}: independent subgrid exists — rule fails domination", flush=True)
    return False

def simultaneous_uniform(n, ts, tlimit=3600):
    """One feature-variable set; constraints from every t in ts. CEGAR per t."""
    qvar = {}
    cnt = [0]
    def q(x, y):
        if x > y:
            x, y = y, x
        key = F2(x, y)
        if key not in qvar:
            cnt[0] += 1
            qvar[key] = cnt[0]
        return qvar[key]
    sol = Glucose3()
    seen = set()
    Ls = {}
    vers = {}
    for t in ts:
        L = leaves(n, t)
        Ls[t] = L
        N = len(L)
        for i, j, k in itertools.combinations(range(N), 3):
            c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
            if c not in seen:
                seen.add(c)
                sol.add_clause(list(c))
        for g in seed_subgrids(n, t):
            sol.add_clause(sorted(set(q(a, b) for a, b in itertools.combinations(g, 2))))
        vers[t] = SubgridVerifier(n, t)
    t0 = time.time()
    added = 0
    while True:
        if time.time() - t0 > tlimit:
            print("   simultaneous: TIMEOUT", flush=True)
            return None
        if not sol.solve():
            print(f"   simultaneous over t={ts}: feature-UNSAT (lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        blue = set(k for k, v in qvar.items() if v in model)
        violated = False
        for t in ts:
            L = Ls[t]
            N = len(L)
            edges = set()
            for i in range(N):
                for j in range(i + 1, N):
                    if F2(L[i], L[j]) in blue:
                        edges.add((i, j))
            bad = vers[t].find(edges)
            if bad is not None:
                sol.add_clause(sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2))))
                added += 1
                violated = True
                break
        if not violated:
            print(f"   simultaneous over t={ts}: SAT — ONE table works at all sizes "
                  f"({len(blue)} blue classes, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return blue

def main():
    # (a) frozen t=5 table re-derived, then tested at 6,7,8
    print("=== (a) re-derive a t=5 table, freeze, test at t=6,7,8 ===", flush=True)
    from erdos592_bidyadic_rule_macmini_s3 import solve_featured
    r, wit = solve_featured(3, 5, F2, "F2", verbose=False)
    assert r
    edges, L, model, qvar = wit
    T5 = set(k for k, v in qvar.items() if v in model)
    print(f"   frozen table from t=5: {len(T5)} blue classes", flush=True)
    for t in (5, 6, 7):
        check_rule(3, t, T5, "frozen-T5")

    # (b) simultaneous t=4..7
    print("\n=== (b) simultaneous-t SAT for one uniform table (t=4,5,6,7) ===", flush=True)
    blue = simultaneous_uniform(3, (4, 5, 6, 7))
    if blue:
        print("   UNIFORM TABLE:", flush=True)
        for k in sorted(blue):
            print("      ", k)
        print("\n=== (c) the uniform table tested ONE SIZE BEYOND its training: t=8 ===", flush=True)
        check_rule(3, 8, blue, "uniform")

if __name__ == "__main__":
    main()
