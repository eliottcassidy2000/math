#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, session 2 part 1 — SAT-based COMPLETE verifier + the Q(3,5) frontier.
mac-mini-2026-06-09-S2  (T768, HYP-2363/2364, THM-453 follow-up)

Replaces the Python-DFS complete verifier (the Q(3,4) bottleneck: 692s) with a
persistent SAT instance:
  * activity vars a_p for every tree prefix p; leaves: s_v = a_v.
  * shape: root active; child active => parent active; active internal node has
    >= 2 active children (for each child c: a_p -> OR_{c'!=c} a_{c'}) and
    <= 2 active children (forbid triples).  Active leaves = selected subgrid leaves.
    With exactly-2 branching enforced at every active internal node, models are
    EXACTLY the binary subgrids.
  * independence on the CURRENT graph via relaxation literals: for every leaf pair
    {u,v} one clause  (-s_u OR -s_v OR r_uv)  added ONCE; per CEGAR round we solve
    under assumptions { -r_e : e an edge of the candidate model }.  Complete search,
    millisecond solves, nothing rebuilt.
SOUNDNESS CHECK: cross-validated against the part-6 exhaustive DFS verifier on all
small cases before use (printed below); both must agree (MISTAKE-067 paranoia).

Also: clause PRE-SEEDING (all subgrids with contiguous child pairs) to cut CEGAR
rounds; invariant (gap-graded quotient) mode included.

Targets: validate; decide Q(3,5); invariant cutoffs n=2,3; curiosity Q(4,3).
"""

import itertools, time, sys
from pysat.solvers import Glucose3

def leaves(n, t):
    return list(itertools.product(range(t), repeat=n))

def prefixes(n, t):
    out = [()]
    for ln in range(1, n + 1):
        out += list(itertools.product(range(t), repeat=ln))
    return out

# ---------------------------------------------------------------------------------
# the old exhaustive DFS verifier (from part 6) for cross-validation
def find_subgrid_dfs(n, t, adj, L, idx):
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
        for c1 in range(t):
            for _ in gen(prefix + (c1,)):
                for c2 in range(c1 + 1, t):
                    yield from gen(prefix + (c2,))
    for _ in gen(()):
        return list(chosen)
    return None

# ---------------------------------------------------------------------------------

class SubgridVerifier:
    """Persistent SAT verifier: find an independent binary subgrid of [t]^n w.r.t.
    an edge set given per-call via assumptions."""
    def __init__(self, n, t):
        self.n, self.t = n, t
        self.L = leaves(n, t)
        self.idx = {v: i for i, v in enumerate(self.L)}
        self.var = {}
        self.cnt = 0
        self.sol = Glucose3()
        P = prefixes(n, t)
        for p in P:
            self.cnt += 1
            self.var[('a', p)] = self.cnt
        a = lambda p: self.var[('a', p)]
        # root active
        self.sol.add_clause([a(())])
        for p in P:
            if len(p) == n:
                continue
            kids = [p + (c,) for c in range(self.t)]
            # child => parent
            for k in kids:
                self.sol.add_clause([-a(k), a(p)])
            # at most 2 active children
            for trip in itertools.combinations(kids, 3):
                self.sol.add_clause([-a(trip[0]), -a(trip[1]), -a(trip[2])])
            # active => at least 2 active children
            for c in range(self.t):
                self.sol.add_clause([-a(p)] + [a(p + (c2,)) for c2 in range(self.t) if c2 != c])
        # relaxation literals for every leaf pair
        self.rel = {}
        for i in range(len(self.L)):
            for j in range(i + 1, len(self.L)):
                self.cnt += 1
                self.rel[(i, j)] = self.cnt
                self.sol.add_clause([-a(self.L[i]), -a(self.L[j]), self.cnt])

    def find(self, edges):
        """edges: set of (i,j) with i<j. Returns independent subgrid leaves or None."""
        assum = [-self.rel[e] for e in edges]
        if not self.sol.solve(assumptions=assum):
            return None
        model = set(l for l in self.sol.get_model() if l > 0)
        sel = [v for v in self.L if self.var[('a', v)] in model]
        return sel

# ---------------------------------------------------------------------------------

def seed_subgrids(n, t, limit=60000):
    """Subgrids with contiguous child pairs (c, c+1) at every node — a structured seed."""
    out = []
    def build(level):
        if level == n:
            return [((),)]
        subs = build(level + 1)
        res = []
        for c in range(t - 1):
            for s1 in subs:
                for s2 in subs:
                    res.append(tuple((c,) + x for x in s1) + tuple((c + 1,) + x for x in s2))
                    if len(res) > limit:
                        return res
        return res
    return build(0)

def solve_Q2(n, t, invariant=False, tlimit=1800, seed=True, verbose=True):
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    qvar = {}
    cnt = [0]
    def q(x, y):
        if x > y:
            x, y = y, x
        if invariant:
            a0, a2 = x[0], y[0]
            key = (a2 - a0, x[1:], y[1:]) if a2 != a0 else (0,) + tuple(sorted((x[1:], y[1:])))
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
    nseed = 0
    if seed:
        for g in seed_subgrids(n, t):
            cl = sorted(set(q(a, b) for a, b in itertools.combinations(g, 2)))
            sol.add_clause(cl)
            nseed += 1

    ver = SubgridVerifier(n, t)
    rounds, added = 0, 0
    t0 = time.time()
    while True:
        rounds += 1
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT ({added} lazy clauses, {time.time()-t0:.0f}s)")
            return None, None
        if not sol.solve():
            if verbose:
                print(f"   UNSAT ({'inv ' if invariant else ''}n={n},t={t}; seed={nseed}, "
                      f"lazy={added}, {time.time()-t0:.1f}s)")
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    edges.add((i, j))
        bad = ver.find(edges)
        if bad is None:
            if verbose:
                print(f"   SAT  ({'inv ' if invariant else ''}n={n},t={t}; {len(edges)} edges, "
                      f"seed={nseed}, lazy={added}, {time.time()-t0:.1f}s)")
            return True, (edges, L, model, qvar)
        cl = sorted(set(q(a, b) for a, b in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1

def crossvalidate():
    print("=== cross-validation: SAT verifier vs exhaustive DFS verifier ===")
    import random
    rng = random.Random(7)
    ok = True
    for (n, t) in ((2, 3), (2, 4), (3, 2), (3, 3)):
        L = leaves(n, t); idx = {v: i for i, v in enumerate(L)}; N = len(L)
        ver = SubgridVerifier(n, t)
        for trial in range(40):
            edges = set()
            adj = [[False]*N for _ in range(N)]
            for i in range(N):
                for j in range(i+1, N):
                    if rng.random() < (0.08 + 0.2*rng.random()):
                        edges.add((i, j)); adj[i][j] = adj[j][i] = True
            r1 = ver.find(edges)
            r2 = find_subgrid_dfs(n, t, adj, L, idx)
            if (r1 is None) != (r2 is None):
                ok = False
                print(f"   DISAGREEMENT n={n} t={t} trial={trial}: SAT={r1} DFS={r2}")
            if r1 is not None:
                # check r1 is genuinely independent + a subgrid shape (leaf count)
                assert len(r1) == 2**n, (n, t, r1)
                for a, b in itertools.combinations(r1, 2):
                    i, j = idx[a], idx[b]
                    assert not adj[i][j], "verifier returned non-independent set!"
        print(f"   n={n} t={t}: 40 random graphs, verdicts agree, witnesses valid")
    print("   CROSS-VALIDATION:", "PASS" if ok else "FAIL")
    return ok

def main():
    if not crossvalidate():
        print("ABORT: verifier mismatch"); return
    print("\n=== sanity re-derive (must match part 6): ===")
    for (n, t, expect) in ((2, 4, True), (2, 5, False), (3, 3, True)):
        r, _ = solve_Q2(n, t)
        print(f"   Q({n},{t}) = {r}  (expected {expect})  {'OK' if r == expect else '*** MISMATCH ***'}")

    print("\n=== THE FRONTIER: Q(3,5) free ===", flush=True)
    r, wit = solve_Q2(3, 5, tlimit=3600)
    if r is True:
        print("   --> R(3,2) > 5: strong-witness shadow persists (HYP-2363 supported)")
    elif r is False:
        print("   --> R(3,2) = 5 !!! same cutoff as n=2 — strong witnesses die at ω³ too")

    print("\n=== invariant cutoffs ===", flush=True)
    for t in range(2, 7):
        r, wit = solve_Q2(2, t, invariant=True, tlimit=900)
        print(f"   invQ(2,{t}) = {r}")
        if r is False:
            print(f"   --> invariant n=2 cutoff = {t}")
            break
    for t in range(2, 7):
        r, wit = solve_Q2(3, t, invariant=True, tlimit=1800)
        print(f"   invQ(3,{t}) = {r}")
        if r is False:
            print(f"   --> invariant n=3 cutoff = {t}")
            break

    print("\n=== curiosity: n=4 ===", flush=True)
    for t in (2, 3):
        r, wit = solve_Q2(4, t, tlimit=900)
        print(f"   Q(4,{t}) = {r}")

if __name__ == "__main__":
    main()
