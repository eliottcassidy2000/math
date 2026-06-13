#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 session 3, part 1 — THE BI-DYADIC RULE SHORTCUT (HYP-2374, THM-465;
POKE Steering Task 2).   mac-mini-2026-06-10-S1

Idea: instead of deciding Q(3,5) by raw SAT (80k+ CEGAR clauses, timeout), quotient
the game by a small FEATURE ALGEBRA and decide the quotient — then VERIFY any witness
explicitly with the complete verifier (the easy direction).  A verified witness at
t=5 settles Q(3,5) = SAT with an explicit, closed-form-ish rule; feature-UNSAT only
means the algebra is too coarse (logged honestly).

Feature algebras for a pair x <lex y in [t]^3 (gaps d_i in coordinate i):
  F1 "signs"      : (per-coordinate sign pattern)                       [known-weak]
  F2 "bi-dyadic"  : sign + v2(|gap|) per coordinate                     [the candidate]
  F3 "bi-dyadic+" : F2 + v2-class of the meet level's child values     [fallback]

S1 background: pure order-patterns provably cannot witness at n=3; the S2 dyadic
result (B_g through v2(row-gap)) says the ROW grading factors through v2 — F2 asks
whether v2 works at ALL THREE coordinate levels simultaneously.
"""

import itertools, time
from pysat.solvers import Glucose3
from erdos592_satverifier_frontier_macmini_s2 import SubgridVerifier, leaves, seed_subgrids

def v2(g):
    k = 0
    while g % 2 == 0:
        g //= 2; k += 1
    return k

def sgn_v2(d):
    if d == 0:
        return ('=',)
    return ('+' if d > 0 else '-', v2(abs(d)))

def F1(x, y):
    return tuple(('=',) if y[i] == x[i] else ('+' if y[i] > x[i] else '-',) for i in range(len(x)))

def F2(x, y):
    return tuple(sgn_v2(y[i] - x[i]) for i in range(len(x)))

def F3(x, y):
    base = F2(x, y)
    m = 0
    while m < len(x) and x[m] == y[m]:
        m += 1
    return base + ((m, x[m] % 2, y[m] % 2),) if m < len(x) else base + (('eq',),)

def solve_featured(n, t, F, name, tlimit=900, verbose=True):
    L = leaves(n, t)
    N = len(L)
    qvar = {}
    cnt = [0]
    def q(x, y):
        if x > y:
            x, y = y, x
        key = F(x, y)
        if key not in qvar:
            cnt[0] += 1
            qvar[key] = cnt[0]
        return qvar[key]
    sol = Glucose3()
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if len(c) == 1:
            # a triangle whose three pair-features coincide: that feature is forced OFF
            c = c
        if c not in seen:
            seen.add(c)
            sol.add_clause(list(c))
    for g in seed_subgrids(n, t):
        sol.add_clause(sorted(set(q(a, b) for a, b in itertools.combinations(g, 2))))
    ver = SubgridVerifier(n, t)
    t0 = time.time(); added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose: print(f"   [{name}] TIMEOUT n={n} t={t}", flush=True)
            return None, None
        if not sol.solve():
            if verbose: print(f"   [{name}] feature-UNSAT n={n} t={t} (features={cnt[0]}, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    edges.add((i, j))
        bad = ver.find(edges)
        if bad is None:
            if verbose: print(f"   [{name}] SAT n={n} t={t} ({len(edges)} edges, features={cnt[0]}, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return True, (edges, L, model, qvar)
        sol.add_clause(sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2))))
        added += 1

def independent_verify(n, t, edges):
    """Fresh complete verification of an explicit graph (paranoia per MISTAKE-067):
    new verifier instance + triangle check."""
    L = leaves(n, t); idx = {v: i for i, v in enumerate(L)}; N = len(L)
    adj = [[False]*N for _ in range(N)]
    for (i, j) in edges:
        adj[i][j] = adj[j][i] = True
    for i, j, k in itertools.combinations(range(N), 3):
        assert not (adj[i][j] and adj[i][k] and adj[j][k]), "TRIANGLE in witness!"
    ver = SubgridVerifier(n, t)
    bad = ver.find(set(edges))
    return bad is None

def print_rule(model, qvar, name):
    on = sorted([k for k, v in qvar.items() if v in model])
    off = sorted([k for k, v in qvar.items() if v not in model])
    print(f"   [{name}] BLUE feature classes ({len(on)}):")
    for k in on:
        print("      ", k)
    print(f"   [{name}] ({len(off)} classes red)")

def main():
    print("=== stage 1: which algebras carry a witness at (3,4)? ===", flush=True)
    winners = []
    for F, name in ((F1, "F1-signs"), (F2, "F2-bidyadic"), (F3, "F3-bidyadic+meet")):
        r, wit = solve_featured(3, 4, F, name)
        if r:
            winners.append((F, name))
    print("\n=== stage 2: winners attempt (3,5) — the frontier ===", flush=True)
    for F, name in winners:
        r, wit = solve_featured(3, 5, F, name, tlimit=2400)
        if r:
            edges, L, model, qvar = wit
            ok = independent_verify(3, 5, edges)
            print(f"   [{name}] INDEPENDENT RE-VERIFICATION at t=5: {'PASS' if ok else 'FAIL'}", flush=True)
            if ok:
                print(f"   *** Q(3,5) = SAT, settled by explicit {name} witness ({len(edges)} edges) ***", flush=True)
                print_rule(model, qvar, name)
                print("\n=== stage 3: push the SAME rule family to t=6 ===", flush=True)
                r6, wit6 = solve_featured(3, 6, F, name, tlimit=3600)
                if r6:
                    e6, L6, m6, q6 = wit6
                    ok6 = independent_verify(3, 6, e6)
                    print(f"   t=6 independent verification: {'PASS' if ok6 else 'FAIL'}", flush=True)
                break

if __name__ == "__main__":
    main()
