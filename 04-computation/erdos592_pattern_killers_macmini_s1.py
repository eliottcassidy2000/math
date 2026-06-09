#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős 592 lab, part 2 — hunting GRID-KILLER pattern sets at n=3 (omega^3 level).
mac-mini-2026-06-09-S1  (T768, HYP-2339/2340, follows erdos592_shuffle_pattern_lab)

Goal: decide whether Specker's negative witness for omega^3 -/-> (omega^3,3) can be
PATTERN-MEASURABLE: a triangle-free set S of pair-patterns (all 31, shared values
allowed) such that every t-grid in [V] contains an S-pair (S "kills" grids).

Method:
 1. enumerate ALL maximal triangle-free pattern sets over the full 31-pattern algebra
    (exact; triple table exhaustive over [0,9));
 2. for each, exact backtracking search for S-free t-grids in [V];
 3. report killers (bounded t) vs avoidable (growing t).

If NO pattern set kills grids while Specker's theorem says a witness exists, then the
witness is NOT pattern-measurable (refutes the strong form of HYP-2339 at n=3) and must
use scheme/value arithmetic — itself a sharp finding.
"""

import itertools, sys
sys.setrecursionlimit(100000)

from erdos592_shuffle_pattern_lab_macmini_s1 import (
    pattern, all_patterns, triangle_table, is_triangle_free,
    exists_S_free_grid)

def maximal_triangle_free_sets_capped(pats, triples, cap=20000):
    pats = list(pats)
    results = set()
    count = [0]

    def extend(S, candidates):
        if count[0] > cap:
            return
        extended = False
        for i, p in enumerate(candidates):
            S2 = S | {p}
            if is_triangle_free(S2, triples):
                extended = True
                extend(S2, candidates[i + 1:])
        if not extended:
            if all(not is_triangle_free(S | {q}, triples) for q in pats if q not in S):
                if frozenset(S) not in results:
                    results.add(frozenset(S))
                    count[0] += 1

    extend(set(), pats)
    return sorted(results, key=lambda s: (-len(s), sorted(s)))

def main():
    n = 3
    pats = all_patterns(n)
    print(f"n={n}: {len(pats)} patterns total")
    tri = triangle_table(n)
    print(f"triple table: {len(tri)} realizable pattern-triples")

    mtf = maximal_triangle_free_sets_capped(pats, tri)
    print(f"\nmaximal triangle-free sets over ALL patterns: {len(mtf)}")
    sizes = {}
    for S in mtf:
        sizes[len(S)] = sizes.get(len(S), 0) + 1
    print("size distribution:", dict(sorted(sizes.items())))

    print("\n--- grid hunt (exact backtracking) ---")
    print("t-grid = full t-branching height-3 tree, values increasing along branches.")
    killers = []
    for idx, S in enumerate(mtf):
        caps = []
        for (V, tmax) in ((8, 3), (11, 3), (14, 4)):
            best = 0
            for t in range(1, tmax + 1):
                if exists_S_free_grid(n, t, V, S):
                    best = t
                else:
                    break
            caps.append(best)
        flag = ""
        if caps[-1] <= 1:
            flag = "  <<<< KILLER (no 2-grid free at V=14)"
            killers.append((S, caps))
        elif caps[-1] == 2:
            flag = "  <-- candidate (capped at t=2 thru V=14)"
            killers.append((S, caps))
        print(f"[{idx}] |S|={len(S)} caps(V=8,11,14)={caps}{flag}")
        for p in sorted(S):
            print("      ", p)

    print(f"\nKILLER/candidate sets: {len(killers)}")
    if killers:
        print("Pushing the best candidates to larger V to test boundedness...")
        for S, caps in killers[:6]:
            for (V, t) in ((18, 2), (18, 3), (24, 2), (24, 3)):
                r = exists_S_free_grid(n, t, V, S)
                print(f"   S(size {len(S)}) V={V} t={t}: S-free grid exists = {r}")
            print("   S =", sorted(S))

if __name__ == "__main__":
    main()
