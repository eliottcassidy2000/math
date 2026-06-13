#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, session 2 part 3 — MINIMUM witness sizes (new integer sequences) +
the dyadic fingerprint of invariant witnesses.
mac-mini-2026-06-09-S2  (THM-453 E addendum, HYP-2364)

f_grid(t)  := min edges of a triangle-free graph on [t]^2 (binary tree leaves)
              dominating all binary subgrids (defined for t < R(2,2)=5: t=2,3,4).
f_tower(s,C):= same for the Chang M=2 tower game on [C]^s (THM-460), small cases.

Method: CEGAR as before + cardinality bound via pysat ITotalizer; decrease k until
UNSAT — each bound certified with the COMPLETE verifier in the loop.

Also: for the minimal t=4 grid witness, print the gap-graded structure
(B_1,B_2,B_3 between rows by gap, and v_2(g) grouping) to test the dyadic ansatz.
"""

import itertools, time
from pysat.solvers import Glucose3
from pysat.card import ITotalizer

from erdos592_satverifier_frontier_macmini_s2 import SubgridVerifier, leaves

def min_witness_grid(t, tlimit=600):
    n = 2
    L = leaves(n, t)
    idx = {v: i for i, v in enumerate(L)}
    N = len(L)
    evar = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            evar[(i, j)] = cnt
    best = None
    bestE = None
    # outer loop over decreasing cardinality bounds
    ub = N * N  # no bound first
    while True:
        sol = Glucose3()
        top = cnt
        for i, j, k in itertools.combinations(range(N), 3):
            sol.add_clause([-evar[(i, j)], -evar[(i, k)], -evar[(j, k)]])
        tot = None
        if best is not None:
            tot = ITotalizer(lits=list(evar.values()), ubound=best - 1, top_id=cnt + 1)
            for cl in tot.cnf.clauses:
                sol.add_clause(cl)
            sol.add_clause([-tot.rhs[best - 1]])  # at most best-1 edges
        ver = SubgridVerifier(n, t)
        t0 = time.time()
        while True:
            if time.time() - t0 > tlimit:
                return best, bestE, "TIMEOUT"
            if not sol.solve():
                return best, bestE, "OPT"
            model = set(l for l in sol.get_model() if l > 0)
            edges = set(e for e, v in evar.items() if v in model)
            bad = ver.find(edges)
            if bad is None:
                best = len(edges)
                bestE = edges
                break  # tighten bound and restart
            sol.add_clause([evar[(min(idx[a], idx[b]), max(idx[a], idx[b]))]
                            for a, b in itertools.combinations(bad, 2)])

def analyse_gaps(edges, L, t):
    from collections import defaultdict
    Bg = defaultdict(set)
    R = defaultdict(set)
    for (i, j) in edges:
        a, b = L[i], L[j]
        if a[0] == b[0]:
            R[a[0]].add((a[1], b[1]))
        else:
            g = abs(b[0] - a[0])
            lo, hi = (a, b) if a[0] < b[0] else (b, a)
            Bg[g].add((lo[1], hi[1]))
    print("      within-row edges by row:", {k: sorted(v) for k, v in sorted(R.items())})
    for g in sorted(Bg):
        print(f"      gap {g} (v2={bin(g)[::-1].index('1') if g else '-'}): {sorted(Bg[g])}")

def main():
    print("=== minimum grid-witness sizes f_grid(t), t=2..4 (R(2,2)=5) ===", flush=True)
    seq = []
    for t in (2, 3, 4):
        best, E, status = min_witness_grid(t)
        seq.append(best)
        print(f"   f_grid({t}) = {best}   [{status}]", flush=True)
        if t == 4 and E is not None:
            print("   minimal t=4 witness structure (gap-graded):")
            analyse_gaps(E, leaves(2, t), t)
    print("   sequence f_grid:", seq, " (t=2,3,4; not in OEIS to our knowledge — check)")

if __name__ == "__main__":
    main()
