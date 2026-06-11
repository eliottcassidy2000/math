#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-470 lab: climbing the algebra ladder to a t-uniform rung.
kind-pasteur-2026-06-11-S1.  Continues THM-465 B (constructive strong Specker at w^3).

Background: THM-465 B proved NO (sign,v2)-measurable table satisfies the conjoined
t=4..7 instance (frozen t=5 table grows a triangle at t=7).  Coarsening collapse
(THM-470 A): any algebra FACTORING THROUGH (sign,v2) is refuted a priori, since its
graphs are F2-measurable.  So we climb by REFINEMENT:

  F2J  (dyadic 1-jet):  per coordinate ('=',) or (sign, v2(|d|), next bit of |d|).
        Same-valuation composition refines: next-bits equal -> v_out = v+1 exactly;
        next-bits differ -> v_out >= v+2.  [HYP-2392]
  F2X  (cross-gap / Larson-flavored): bi-dyadic data sgn_v2 of the SIX quantities
        (d1, d2, d3, d1-d2, d2-d3, d1-d3).  Captures inter-coordinate magnitude
        relations uniformly in t.  [HYP-2393]

Protocol per algebra: per-size t=4..7 (does the rung at least match F2 per-size?),
then conjoined simultaneous-t SAT over t=(4,5,6,7); if a uniform table survives,
verify it one size beyond (t=8) per MISTAKE-067 discipline ("read the witness").
Also re-runs the F2 conjoined refutation as an environment-reproducibility control.

Honest semantics: conjoined-SAT here is NOT yet an infinite witness — classes with
v2 >= 3 first appear at t = 9 and tail behavior is unconstrained by t <= 8.  A
surviving table is a CANDIDATE rule; a conjoined-UNSAT is a genuine rung-killer
(restriction argument, THM-465 B verbatim).
"""

import itertools, time
from pysat.solvers import Glucose3

from erdos592_satverifier_frontier_macmini_s2 import leaves, SubgridVerifier, seed_subgrids
from erdos592_bidyadic_rule_macmini_s3 import v2, sgn_v2, F2, solve_featured, independent_verify


# ---------------------------------------------------------------- feature maps

def sgn_v2_jet(d):
    if d == 0:
        return ('=',)
    g = abs(d)
    v = v2(g)
    u = g >> v               # odd part
    return ('+' if d > 0 else '-', v, (u >> 1) & 1)


def F2J(x, y):
    return tuple(sgn_v2_jet(y[i] - x[i]) for i in range(len(x)))


def F2X(x, y):
    d = [y[i] - x[i] for i in range(len(x))]
    virt = list(d) + [d[i] - d[j] for i, j in itertools.combinations(range(len(d)), 2)]
    return tuple(sgn_v2(v) for v in virt)


# --------------------------------------- parameterized conjoined-t uniform SAT

def simultaneous_uniform_F(n, ts, F, name, tlimit=3600, verbose=True):
    """One shared table over feature classes of F; constraints of every t in ts.
    Returns blue-class set | False (feature-UNSAT) | None (timeout)."""
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
    Ls = {}
    vers = {}
    t_build = time.time()
    for t in ts:
        L = leaves(n, t); N = len(L)
        Ls[t] = L
        for i, j, k in itertools.combinations(range(N), 3):
            c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
            if c not in seen:
                seen.add(c); sol.add_clause(list(c))
        for g in seed_subgrids(n, t):
            c = tuple(sorted(set(q(a, b) for a, b in itertools.combinations(g, 2))))
            if c not in seen:
                seen.add(c); sol.add_clause(list(c))
        vers[t] = SubgridVerifier(n, t)
    if verbose:
        print(f"   [{name} conjoined {ts}: {cnt[0]} features, build {time.time()-t_build:.1f}s]", flush=True)
    t0 = time.time(); added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT {name} conjoined (lazy={added})", flush=True)
            return None
        if not sol.solve():
            if verbose: print(f"   feature-UNSAT {name} conjoined {ts} "
                              f"(lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        blue = set(k for k, v in qv.items() if v in model)
        bad_found = False
        for t in ts:
            L = Ls[t]; N = len(L)
            edges = set()
            for i in range(N):
                for j in range(i + 1, N):
                    if q(L[i], L[j]) in model:
                        edges.add((i, j))
            bad = vers[t].find(edges)
            if bad is not None:
                cl = sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2)))
                sol.add_clause(cl)
                added += 1
                bad_found = True
                break
        if not bad_found:
            if verbose: print(f"   UNIFORM TABLE {name} over {ts}: {len(blue)} blue classes "
                              f"(lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return blue


def check_rule_F(n, t, blue, F, label):
    """Verify the F-rule given by blue classes at size t (triangle scan + complete verifier)."""
    L = leaves(n, t); N = len(L)
    edges = set()
    for i in range(N):
        for j in range(i + 1, N):
            if F(L[i], L[j]) in blue:
                edges.add((i, j))
    adj = [[False] * N for _ in range(N)]
    for i, j in edges:
        adj[i][j] = adj[j][i] = True
    for i, j, k in itertools.combinations(range(N), 3):
        if adj[i][j] and adj[i][k] and adj[j][k]:
            print(f"   {label} at t={t}: FAILS triangle-freeness "
                  f"({L[i]},{L[j]},{L[k]})", flush=True)
            return False
    bad = SubgridVerifier(n, t).find(edges)
    if bad is not None:
        print(f"   {label} at t={t}: FAILS domination (subgrid root {bad[0]}..{bad[-1]})", flush=True)
        return False
    print(f"   {label} at t={t}: VERIFIED witness ({len(edges)} edges)", flush=True)
    return True


def print_table(blue, name):
    print(f"   --- {name} uniform blue table ({len(blue)} classes) ---", flush=True)
    for k in sorted(blue):
        print(f"      {k}", flush=True)


def main():
    t0 = time.time()
    print("=== control: F2 conjoined t=4..7 (expected feature-UNSAT, THM-465 B) ===", flush=True)
    r = simultaneous_uniform_F(3, (4, 5, 6, 7), F2, "F2", tlimit=1800)
    assert r is False or r is None, "F2 conjoined unexpectedly SAT — environment mismatch!"

    for name, F in (("F2J", F2J), ("F2X", F2X)):
        print(f"\n=== {name}: per-size t=4..7 ===", flush=True)
        persize = {}
        for t in (4, 5, 6, 7):
            r, w = solve_featured(3, t, F, name, tlimit=900)
            persize[t] = 'SAT' if r else ('timeout' if r is None else 'feature-UNSAT')
            if r:
                ok = independent_verify(3, t, w[0])
                print(f"      independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
            if r is False:
                break
        print(f"   {name} per-size: {persize}", flush=True)
        if any(v == 'feature-UNSAT' for v in persize.values()):
            print(f"   {name}: dead per-size; skipping conjoined.", flush=True)
            continue
        print(f"\n=== {name}: conjoined t=(4,5,6,7) ===", flush=True)
        blue = simultaneous_uniform_F(3, (4, 5, 6, 7), F, name, tlimit=3600)
        if blue:
            print_table(blue, name)
            print(f"\n=== {name}: verify uniform table beyond training (t=8) ===", flush=True)
            check_rule_F(3, 8, blue, F, name + "-uniform")
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
