#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 session 3, part 2 — THE SEAM-FOLLOWS-BRANCHING EXPERIMENT (HYP-2373,
THM-464).   mac-mini-2026-06-10-S1

The b-ary subgrid game Q_b(n,t): triangle-free graph on [t]^n leaves with no
independent b-ary subgrid (choose b children per node; b^n leaves).
  * Bridge (THM-453 D1 verbatim for any fixed b): Q_b SAT ∀t ⟹ ω^n ↛ (ω^n,3)².
  * Calibration row (PROVED one-liner): at n=1 the b-ary subgrids are ALL b-subsets,
    so a witness ⟺ triangle-free with independence number ≤ b−1, hence
    R_b(1) = R(3,b) — classical Ramsey numbers are the height-1 row of the table.
    (R_2(1)=3, R_3(1)=6, R_4(1)=9, R_5(1)=14.)

THE EXPERIMENT at n=2: cutoffs for {b=2, b=3} × {free, invariant, v2-graded,
v3-graded}.  S2 found binary: free = invariant = v2 = 5.  HYP-2373 predicts the
sufficient grading follows the branching (v_b costs nothing, v_p for p≠b costs);
the v3-in-binary and v2-in-ternary cells are the discriminating controls.

Verifier: persistent SAT instance generalized to branching b (exactly-b active
children per active node), relaxation-literal assumptions as in S2; cross-validated
against brute force on small cases below.
"""

import itertools, time, random
from pysat.solvers import Glucose3

def leaves(n, t):
    return list(itertools.product(range(t), repeat=n))

def prefixes(n, t):
    out = [()]
    for ln in range(1, n + 1):
        out += list(itertools.product(range(t), repeat=ln))
    return out

class BVerifier:
    """Find an independent b-ary subgrid of [t]^n w.r.t. assumption-supplied edges."""
    def __init__(self, n, t, b):
        self.n, self.t, self.b = n, t, b
        self.L = leaves(n, t)
        self.var = {}
        self.cnt = 0
        self.sol = Glucose3()
        P = prefixes(n, t)
        for p in P:
            self.cnt += 1
            self.var[('a', p)] = self.cnt
        a = lambda p: self.var[('a', p)]
        self.sol.add_clause([a(())])
        for p in P:
            if len(p) == n:
                continue
            kids = [p + (c,) for c in range(t)]
            for k in kids:
                self.sol.add_clause([-a(k), a(p)])
            # at most b active children
            for sub in itertools.combinations(kids, b + 1):
                self.sol.add_clause([-a(x) for x in sub])
            # active => at least b active children: every (t-b+1)-subset hits one
            for sub in itertools.combinations(kids, t - b + 1):
                self.sol.add_clause([-a(p)] + [a(x) for x in sub])
        self.rel = {}
        for i in range(len(self.L)):
            for j in range(i + 1, len(self.L)):
                self.cnt += 1
                self.rel[(i, j)] = self.cnt
                self.sol.add_clause([-a(self.L[i]), -a(self.L[j]), self.cnt])

    def find(self, edges):
        if not self.sol.solve(assumptions=[-self.rel[e] for e in edges]):
            return None
        model = set(l for l in self.sol.get_model() if l > 0)
        return [v for v in self.L if self.var[('a', v)] in model]

def brute_independent_bgrid(n, t, b, adj, L):
    idx = {v: i for i, v in enumerate(L)}
    chosen = []
    def gen(prefix):
        k = len(prefix)
        if k == n:
            i = idx[prefix]
            for c in chosen:
                if adj[i][idx[c]]:
                    return
            chosen.append(prefix); yield; chosen.pop()
            return
        for combo in itertools.combinations(range(t), b):
            def rec(ci):
                if ci == b:
                    yield
                    return
                for _ in gen(prefix + (combo[ci],)):
                    yield from rec(ci + 1)
            yield from rec(0)
    for _ in gen(()):
        return list(chosen)
    return None

def v2(g):
    k = 0
    while g % 2 == 0:
        g //= 2; k += 1
    return k

def v3(g):
    k = 0
    while g % 3 == 0:
        g //= 3; k += 1
    return k

def keyfun(mode):
    def q(x, y):
        a0, a2 = x[0], y[0]
        if a0 == a2:
            return (0,) + tuple(sorted((x[1:], y[1:])))
        g = a2 - a0
        if mode == "free":
            return (x, y)
        if mode == "inv":
            return ('g', g, x[1:], y[1:])
        if mode == "v2":
            return ('v2', v2(g), x[1:], y[1:])
        if mode == "v3":
            return ('v3', v3(g), x[1:], y[1:])
    return q

def solve_game(n, t, b, mode, tlimit=1200, verbose=True):
    L = leaves(n, t)
    N = len(L)
    qv = {}; cnt = [0]
    kf = keyfun(mode)
    def q(x, y):
        if x > y:
            x, y = y, x
        key = kf(x, y)
        if key not in qv:
            cnt[0] += 1; qv[key] = cnt[0]
        return qv[key]
    sol = Glucose3()
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if c not in seen:
            seen.add(c); sol.add_clause(list(c))
    ver = BVerifier(n, t, b)
    t0 = time.time(); added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT b={b} {mode} n={n} t={t}", flush=True)
            return None
        if not sol.solve():
            if verbose: print(f"   UNSAT b={b} {mode:4s} n={n} t={t} (lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    edges.add((i, j))
        bad = ver.find(edges)
        if bad is None:
            if verbose: print(f"   SAT   b={b} {mode:4s} n={n} t={t} ({len(edges)} edges, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return True
        idx = {v: i for i, v in enumerate(L)}
        cl = sorted(set(q(a, c) for a, c in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1

def crossvalidate():
    print("=== BVerifier cross-validation (b=3) vs brute ===", flush=True)
    rng = random.Random(17)
    for (n, t) in ((1, 4), (1, 5), (2, 3), (2, 4)):
        L = leaves(n, t); N = len(L)
        if 3 ** n > N:
            continue
        ver = BVerifier(n, t, 3)
        bad = 0
        for tr in range(25):
            adj = [[False] * N for _ in range(N)]
            edges = set()
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < 0.3 * rng.random():
                        adj[i][j] = adj[j][i] = True
                        edges.add((i, j))
            r1 = ver.find(edges)
            r2 = brute_independent_bgrid(n, t, 3, adj, L)
            if (r1 is None) != (r2 is None):
                bad += 1
                print(f"   DISAGREE n={n} t={t} tr={tr}", flush=True)
            if r1 is not None:
                idx = {v: i for i, v in enumerate(L)}
                assert len(r1) == 3 ** n
                for a, c in itertools.combinations(r1, 2):
                    assert not adj[idx[a]][idx[c]]
        print(f"   n={n} t={t}: 25 trials {'OK' if bad == 0 else f'{bad} BAD'}", flush=True)

def main():
    crossvalidate()
    print("\n=== calibration: R_3(1) should be R(3,3) = 6 ===", flush=True)
    for t in (4, 5, 6):
        solve_game(1, t, 3, "free")
    print("\n=== binary game n=2: the v3 CONTROL (free/inv/v2 known = 5) ===", flush=True)
    for t in (4, 5):
        solve_game(2, t, 2, "v3")
    print("\n=== ternary game n=2: free / inv / v2 / v3 sweep ===", flush=True)
    cut = {}
    for mode in ("free", "inv", "v3", "v2"):
        for t in range(3, 13):
            r = solve_game(2, t, 3, mode, tlimit=1800)
            if r is False:
                cut[mode] = t
                print(f"   --> ternary n=2 cutoff [{mode}] = {t}", flush=True)
                break
            if r is None:
                cut[mode] = f">{t-1} (timeout)"
                break
    print("\nSUMMARY ternary cutoffs:", cut, flush=True)

if __name__ == "__main__":
    main()
