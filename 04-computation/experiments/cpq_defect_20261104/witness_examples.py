#!/usr/bin/env python3
"""witness_examples.py -- print the witness configurations (sign-odd graph types with an N-even H-fixed representative)
at the first defect levels of the C_pq (p,q,1^k) families, with the odd automorphism found by nauty."""
import sys, pynauty
import cpq_defect_census as m
def describe(p, q, k):
    row, wit = m.census(p, q, k, verbose=False)
    print("(%d,%d,1^%d): defect=%d W=%d E_T=%d E_even=%d E_wit=%d" % (p, q, k, row["defect"], row["W"], row["E_T"], row["E_even"], row["E_wit"]))
    n = p + q + k
    for cert, (Cp, Cq, cross, adjF, col) in wit.items():
        adj = m.build(p, q, k, Cp, Cq, cross, adjF, col, False)
        edges = [(u, v) for u in adj for v in adj[u] if u < v]
        g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
        gens, s1, s2, orb, no = pynauty.autgrp(g)
        odd = [gg for gg in gens if m.reversal_parity(gg, edges)]
        deg = sorted(len(adj[v]) for v in range(n))
        print("  Cp=%s Cq=%s cross=%d F-graph=%s colours=%s |Aut|=%d degrees=%s odd generator=%s" % (
            sorted(Cp), sorted(Cq), cross, {u: adjF[u] for u in adjF if adjF[u]}, col, int(round(s1 * 10 ** s2)), deg, odd[0] if odd else None))
for (p, q, k) in [(3, 5, 2), (3, 5, 3), (5, 7, 2), (5, 7, 3), (3, 7, 3), (3, 7, 4)]:
    describe(p, q, k)
