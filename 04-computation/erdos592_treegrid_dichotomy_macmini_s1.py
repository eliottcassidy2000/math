#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, part 3 — THE FINITE TREE-GRID DICHOTOMY (the abstract core).
mac-mini-2026-06-09-S1  (T768, HYP-2346, THM-453)

QUESTION Q(n,t): on the leaf set L = [t]^n of the complete t-ary tree of height n,
does there exist a TRIANGLE-FREE graph G such that EVERY binary subgrid contains an
edge of G?  (binary subgrid = choose 2 children at every step hierarchically =
the leaves of an embedded complete binary subtree of height n: 2^n leaves.)

WHY THIS IS THE FINITE SHADOW OF ERDŐS 592 AT omega^n:
  * X ⊆ omega^n has order type omega^n  <=>  X contains a full (omega-branching) grid.
  * alpha -/-> (alpha,3) witness = triangle-free graph killing all full grids.
  [CORRECTED BRIDGE — the original docstring had the implication BACKWARDS
   (see MISTAKES.md MISTAKE-066 and THM-453 part D for the proved version):
     * Q(n,t) SAT for ALL t  ==>  omega^n -/-> (omega^n,3)   (König/compactness:
       coherent restrictions glue to a triangle-free graph on the full grid with no
       independent binary subgrid; full-type sets contain binary subgrids).
     * Hence omega^n -> (omega^n,3) (e.g. Specker at n=2) ==> Q(n,t) UNSAT from some
       finite cutoff R(n,2).  A negative ordinal relation does NOT formally imply
       Q SAT for all t (witnesses need not avoid independent binary subgrids).]
  * n=1 calibration: subgrids = all pairs => G must be complete & triangle-free:
    Q(1,t) SAT iff t <= 2.  (Ramsey omega -> (omega,3) shadow.)
  * EXPECTATION (Specker dichotomy): Q(2,t) dies at small t (omega^2 positive);
    Q(3,t) lives for all t (omega^3 negative), and the SAT witnesses reveal the
    structure of Specker's coloring.

Exact method: DPLL with unit propagation on edge variables.
  clauses: triangle (¬e1 ∨ ¬e2 ∨ ¬e3) for every vertex triple (downward closed need),
           subgrid hitting (e_{p1} ∨ ... ∨ e_{pk}), k = C(2^n, 2) pairs of its leaves.
Plus a fast greedy/local-search to find witnesses quickly when SAT.
"""

import itertools, random, sys
sys.setrecursionlimit(100000)

# ----------------------------------------------------------------------------------

def leaves(n, t):
    return list(itertools.product(range(t), repeat=n))

def binary_subgrids(n, t):
    """All binary subgrids as tuples of leaves (2^n leaves each)."""
    grids = []
    def rec(prefix_choices):
        # represent a subgrid as: at each internal node reached, a sorted pair of children
        # build recursively: returns list of leaf-suffix-sets
        pass
    # explicit recursive construction
    def build(level):
        """Return list of (subgrid leaf-set) for a subtree at given level."""
        if level == n:
            return [((),)]  # single empty-suffix leaf
        subs = build(level + 1)
        out = []
        for c1, c2 in itertools.combinations(range(t), 2):
            for s1 in subs:
                for s2 in subs:
                    out.append(tuple((c1,) + leaf for leaf in s1)
                               + tuple((c2,) + leaf for leaf in s2))
        return out
    return build(0)

def count_binary_subgrids(n, t):
    c = 1
    for level in range(n):
        # number of internal nodes at this level in the binary subtree: 2^level
        c *= (t * (t - 1) // 2) ** (2 ** level)
    return c

# ----------------------------------------------------------------------------------
# DPLL SAT solver (small, exact)
# ----------------------------------------------------------------------------------

class DPLL:
    def __init__(self, nvars, clauses):
        self.nvars = nvars
        self.clauses = clauses  # list of lists of signed ints (1-indexed)

    def solve(self, node_budget=8_000_000):
        clauses = [tuple(c) for c in self.clauses]
        assign = {}
        self.nodes = 0
        self.budget = node_budget

        def unit_propagate(assign, clauses):
            changed = True
            while changed:
                changed = False
                for cl in clauses:
                    unassigned = []
                    sat = False
                    for lit in cl:
                        v = abs(lit)
                        if v in assign:
                            if (lit > 0) == assign[v]:
                                sat = True
                                break
                        else:
                            unassigned.append(lit)
                    if sat:
                        continue
                    if not unassigned:
                        return False  # conflict
                    if len(unassigned) == 1:
                        lit = unassigned[0]
                        assign[abs(lit)] = lit > 0
                        changed = True
            return True

        def rec(assign):
            self.nodes += 1
            if self.nodes > self.budget:
                raise TimeoutError("node budget exceeded")
            assign = dict(assign)
            if not unit_propagate(assign, clauses):
                return None
            # pick unassigned var appearing in most short unsatisfied clauses
            best, bestv = None, -1
            counts = {}
            done = True
            for cl in clauses:
                sat = False
                un = []
                for lit in cl:
                    v = abs(lit)
                    if v in assign:
                        if (lit > 0) == assign[v]:
                            sat = True
                            break
                    else:
                        un.append(v)
                if not sat and un:
                    done = False
                    w = 1.0 / (1 << min(len(un), 20))
                    for v in un:
                        counts[v] = counts.get(v, 0) + w
            if done:
                return assign
            best = max(counts, key=counts.get)
            for val in (True, False):
                assign2 = dict(assign)
                assign2[best] = val
                r = rec(assign2)
                if r is not None:
                    return r
            return None

        return rec(assign)

# ----------------------------------------------------------------------------------
# Greedy / local search witness finder (fast SAT direction)
# ----------------------------------------------------------------------------------

def greedy_witness(n, t, seed=0, iters=200000):
    """Try to build triangle-free graph hitting all binary subgrids."""
    rng = random.Random(seed)
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    grids = binary_subgrids(n, t)
    # subgrid -> list of pairs (i<j)
    gpairs = []
    for g in grids:
        ps = []
        for a, b in itertools.combinations(g, 2):
            i, j = idx[a], idx[b]
            if i > j: i, j = j, i
            ps.append((i, j))
        gpairs.append(ps)
    adj = [[False] * N for _ in range(N)]
    edges = set()

    def creates_triangle(i, j):
        for k in range(N):
            if adj[i][k] and adj[j][k]:
                return True
        return False

    # pair -> subgrids containing it
    from collections import defaultdict
    pair2g = defaultdict(list)
    for gi, ps in enumerate(gpairs):
        for p in ps:
            pair2g[p].append(gi)

    unhit = set(range(len(gpairs)))
    hitcount = [0] * len(gpairs)

    def add(i, j):
        adj[i][j] = adj[j][i] = True
        edges.add((i, j))
        for gi in pair2g[(i, j)]:
            hitcount[gi] += 1
            unhit.discard(gi)

    # greedy: repeatedly pick the unhit subgrid with fewest addable pairs; add the pair
    # covering most other unhit subgrids without creating a triangle
    for _ in range(iters):
        if not unhit:
            return edges, L
        gi = min(unhit, key=lambda g: sum(1 for p in gpairs[g]
                                          if p not in edges and not creates_triangle(*p)))
        cands = [p for p in gpairs[gi] if p not in edges and not creates_triangle(*p)]
        if not cands:
            return None, L  # stuck (greedy failure, not proof of UNSAT)
        best = max(cands, key=lambda p: sum(1 for g2 in pair2g[p] if g2 in unhit)
                   + rng.random() * 0.5)
        add(*best)
    return None, L

# ----------------------------------------------------------------------------------

def make_clauses(n, t):
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    evar = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            evar[(i, j)] = cnt
    clauses = []
    for i, j, k in itertools.combinations(range(N), 3):
        clauses.append([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])
    for g in binary_subgrids(n, t):
        cl = []
        for a, b in itertools.combinations(g, 2):
            i, j = idx[a], idx[b]
            if i > j: i, j = j, i
            cl.append(evar[(i, j)])
        clauses.append(cl)
    return cnt, clauses, evar, L

def analyse_witness(edges, L, n, t):
    """Describe witness structure: classify edges by tree-relation of endpoints."""
    from collections import Counter
    def meet_level(a, b):
        m = 0
        while m < n and a[m] == b[m]:
            m += 1
        return m  # 0 = split at root
    cnt = Counter()
    for (i, j) in edges:
        a, b = L[i], L[j]
        ml = meet_level(a, b)
        # after split: relative order of children at split, then suffix relation coarse
        cnt[ml] += 1
    return dict(sorted(cnt.items()))

def main():
    print("=" * 78)
    print("Q(n,t): triangle-free graph on [t]^n leaves hitting all binary subgrids?")
    print("=" * 78)

    # calibration n=1
    for t in (2, 3):
        nv, cls, evar, L = make_clauses(1, t)
        res = DPLL(nv, cls).solve()
        print(f"n=1 t={t}: #leaves={len(L)} SAT={res is not None}   (expect SAT iff t<=2; Ramsey shadow)")

    # n=2
    print("\n--- n=2 (omega^2: Specker POSITIVE => expect UNSAT from some small t) ---")
    for t in (2, 3, 4, 5):
        ng = count_binary_subgrids(2, t)
        print(f"n=2 t={t}: leaves={t*t} subgrids={ng}", flush=True)
        w, L = greedy_witness(2, t, seed=1)
        if w is None:
            for s in range(2, 6):
                w, L = greedy_witness(2, t, seed=s)
                if w is not None:
                    break
        if w is not None:
            print(f"   greedy WITNESS found: {len(w)} edges; meet-level histogram:",
                  analyse_witness(w, L, 2, t))
            continue
        # greedy failed -> try exact DPLL
        nv, cls, evar, L = make_clauses(2, t)
        print(f"   greedy failed; DPLL on {nv} vars, {len(cls)} clauses...", flush=True)
        try:
            res = DPLL(nv, cls).solve()
            if res is None:
                print(f"   *** Q(2,{t}) UNSAT (exact) ***")
            else:
                edges = {(i, j) for (i, j), v in evar.items() if res.get(v, False)}
                print(f"   DPLL witness: {len(edges)} edges; meet-level histogram:",
                      analyse_witness(edges, L, 2, t))
        except TimeoutError:
            print("   DPLL budget exceeded — UNKNOWN")

    # n=3
    print("\n--- n=3 (omega^3: Specker NEGATIVE => expect SAT all t; read structure) ---")
    for t in (2, 3):
        ng = count_binary_subgrids(3, t)
        print(f"n=3 t={t}: leaves={t**3} subgrids={ng}", flush=True)
        w, L = greedy_witness(3, t, seed=1)
        if w is None:
            for s in range(2, 8):
                w, L = greedy_witness(3, t, seed=s)
                if w is not None:
                    break
        if w is not None:
            print(f"   greedy WITNESS: {len(w)} edges; meet-level histogram:",
                  analyse_witness(w, L, 3, t))
            print("   sample edges (leaf pairs):")
            for (i, j) in sorted(w)[:12]:
                print("      ", L[i], "--", L[j])
        else:
            nv, cls, evar, L = make_clauses(3, t)
            print(f"   greedy failed; DPLL on {nv} vars, {len(cls)} clauses...", flush=True)
            try:
                res = DPLL(nv, cls).solve()
                if res is None:
                    print(f"   *** Q(3,{t}) UNSAT (exact) — would CONTRADICT expectation ***")
                else:
                    edges = {(i, j) for (i, j), v in evar.items() if res.get(v, False)}
                    print(f"   DPLL witness: {len(edges)} edges; histogram:",
                          analyse_witness(edges, L, 3, t))
            except TimeoutError:
                print("   DPLL budget exceeded — UNKNOWN")

if __name__ == "__main__":
    main()
