#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-471 lab, part 2: the XOR feature quotient for the c=2 shape games.
kind-pasteur-2026-06-11-S1.

At c = 2 the ambient A(m; s, 2) is {0,1}^(s^m) and the natural invariant feature
of a pair is its DISAGREEMENT SET D(x,y) ⊆ positions.  Triangle composition is
F2-LINEAR: D(x,z) = D(x,y) Δ D(y,z).  Hence an XOR-measurable rule T (blue iff
D ∈ T) is triangle-free iff T is SUM-FREE IN THE GROUP F2^(s^m) (no D1, D2 ∈ T
with D1ΔD2 ∈ T — note D1=D2 gives 0 ∉ T automatically) — the THM-469 sum-free
seam reappearing as a binary cap-set condition.  Hitting: CEGAR against the
complete ShapeFinder (MISTAKE-067 discipline), exactly the THM-465 shortcut
ported to the THM-460 miniatures.

feature-UNSAT semantics as always: the XOR algebra is too coarse, NOT game-UNSAT.
A SAT verdict + independent re-verification settles the cell.
"""

import itertools, time
from pysat.solvers import Glucose3

from erdos592_shape_miniatures_kp0611 import (
    ambient_points, bt_size, ShapeFinder, is_bt, BudgetExceeded)


def solve_shape_xor(m, M, conv, s, tlimit=2400, verbose=True):
    """XOR-quotient CEGAR for the c=2 shape game Q^(m)(M; s, 2)."""
    L = ambient_points(m, s, 2)
    N = len(L)
    npos = s ** m
    size = bt_size(m, M, conv)
    if size > N:
        print(f"   VACUOUS m={m} M={M} {conv} (s,2): |BT|={size} > N={N}", flush=True)
        return 'vacuous'
    # disagreement set as a bitmask over positions
    def D(i, j):
        x, y = L[i], L[j]
        d = 0
        for k in range(npos):
            if x[k] != y[k]:
                d |= 1 << k
        return d

    qv = {}
    cnt = [0]

    def q(d):
        if d not in qv:
            cnt[0] += 1
            qv[d] = cnt[0]
        return qv[d]

    # triangle clauses purely group-theoretically: for any D1 < D2 (nonzero,
    # distinct), the triple {D1, D2, D1^D2} is realizable in {0,1}^npos as
    # x=0, y=1_{D1}, z=1_{D1^D2}... realizability: ANY pair of distinct nonzero
    # D1, D2 is realized by x=0, y with disagreement D1, z = y ^ 1_{D2}; all
    # three points lie in the ambient (c=2: every 0/1 vector is a point). So
    # ALL group triples are realizable -- add them all.
    sol = Glucose3()
    seen = set()
    t_build = time.time()
    for d1 in range(1, 1 << npos):
        for d2 in range(d1 + 1, 1 << npos):
            d3 = d1 ^ d2
            if d3 < d2:
                continue  # canonical: d1 < d2 < d3 (each triple counted once)
            c = tuple(sorted(set((-q(d1), -q(d2), -q(d3)))))
            if c not in seen:
                seen.add(c)
                sol.add_clause(list(c))
    if verbose:
        print(f"   [xor quotient m={m} M={M} {conv} (s={s},c=2): N={N}, "
              f"{cnt[0]} classes, {len(seen)} triangle clauses, build {time.time()-t_build:.1f}s]",
              flush=True)
    finder = ShapeFinder(m, M, conv, L)
    t0 = time.time()
    added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose:
                print(f"   TIMEOUT xor m={m} M={M} {conv} (s={s}) lazy={added}", flush=True)
            return None
        if not sol.solve():
            if verbose:
                print(f"   feature-UNSAT xor m={m} M={M} {conv} (s={s},c=2) "
                      f"(lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        blue = set(d for d, v in qv.items() if v in model)
        nb = [0] * N
        edges = 0
        for i in range(N):
            for j in range(i + 1, N):
                if D(i, j) in blue:
                    nb[i] |= 1 << j
                    nb[j] |= 1 << i
                    edges += 1
        finder.set_graph(nb, budget=30_000_000, anchor0=True)
        # anchor0 is complete here: the rule graph is XOR-measurable, so the
        # independent-shape family is XOR-translation-invariant (THM-471 C)
        try:
            bad = finder.find_bt()
        except BudgetExceeded:
            if verbose:
                print(f"   TIMEOUT xor m={m} M={M} {conv} (s={s}) — finder node budget "
                      f"exceeded (lazy={added}, {time.time()-t0:.1f}s) [NOT a SAT certificate]",
                      flush=True)
            return None
        if bad is None:
            if verbose:
                print(f"   SAT xor m={m} M={M} {conv} (s={s},c=2) N={N} "
                      f"({edges} edges, |T|={len(blue)}, lazy={added}, {time.time()-t0:.1f}s)",
                      flush=True)
            return True, blue, L, nb
        cl = sorted(set(q(D(a, b)) for a, b in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1


def reverify(m, M, conv, s, blue, L):
    """Fresh independent verification: rebuild graph, O(N^3) triangle scan,
    fresh ShapeFinder."""
    N = len(L)
    npos = s ** m
    def D(i, j):
        x, y = L[i], L[j]
        d = 0
        for k in range(npos):
            if x[k] != y[k]:
                d |= 1 << k
        return d
    nb = [0] * N
    adj = [[False] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            if D(i, j) in blue:
                nb[i] |= 1 << j
                nb[j] |= 1 << i
                adj[i][j] = adj[j][i] = True
    for i, j, k in itertools.combinations(range(N), 3):
        assert not (adj[i][j] and adj[i][k] and adj[j][k]), "triangle!"
    f = ShapeFinder(m, M, conv, L)
    # anchor0 valid: same XOR-measurable graph; no budget — re-verification
    # must be complete or not return
    f.set_graph(nb, anchor0=True)
    return f.find_bt() is None


def main():
    t0 = time.time()
    print("=== calibration: m=2 cells already decided raw ===", flush=True)
    # m=2 M=1 (2,2) was raw-UNSAT; xor quotient must be feature-UNSAT or timeout
    r = solve_shape_xor(2, 1, 'j1', 2, tlimit=600)
    print(f"   (raw UNSAT at this cell; xor verdict above must NOT be SAT)", flush=True)

    print("\n=== m=2 M=2 j0 at s=2 (xor quotient) ===", flush=True)
    r22 = solve_shape_xor(2, 2, 'j0', 2, tlimit=1200)
    if isinstance(r22, tuple):
        ok = reverify(2, 2, 'j0', 2, r22[1], r22[2])
        print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)

    print("\n=== m=3 PROBES at s=2 (xor quotient; the open-case cells) ===", flush=True)
    for (M, conv) in ((1, 'j1'), (2, 'j0')):
        r3 = solve_shape_xor(3, M, conv, 2, tlimit=3000)
        if isinstance(r3, tuple):
            ok = reverify(3, M, conv, 2, r3[1], r3[2])
            print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
            print(f"   blue table (disagreement-set masks): {sorted(r3[1])}", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
