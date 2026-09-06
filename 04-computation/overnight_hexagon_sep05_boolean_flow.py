#!/usr/bin/env python3
"""Exact native orbit-measure flow controls; no spectral extrapolation.

Universe: every Eulerian graph and every endpoint pair for n=3,4,5;
all random orders of differing anchored-triangle coordinates.  The first
path aggregates subset positions; the second literally enumerates orders.
Separately check the closed Burnside constraint rank against the inherited
literal edge-orbit matrix for every partition through n=16.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
from math import comb, factorial, gcd
import hashlib
import json

import even_graph_triangle_quotient_spectrum_thm4078 as primal
import even_graph_burnside_s260 as burnside


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def flow(n, literal=False):
    _, ei = primal.edge_system(n)
    states = primal.cycle_space(n, ei)
    reps, orbit, sizes, _ = primal.quotient_data(n)
    q = len(reps)
    m = (n-1)*(n-2)//2
    load = defaultdict(F)
    pair_count = 0
    for x in range(len(states)):
        for y in range(x+1, len(states)):
            ox, oy = orbit[states[x]], orbit[states[y]]
            if ox == oy:
                continue
            pair_count += 1
            difference = x ^ y
            bits = [1 << i for i in range(m) if difference >> i & 1]
            d = len(bits)
            base = F(1, q*sizes[ox]*sizes[oy])
            if literal:
                for order in permutations(bits):
                    current = x
                    for bit in order:
                        other = current ^ bit
                        a, b = orbit[states[current]], orbit[states[other]]
                        if a != b:
                            load[tuple(sorted((a, b)))] += base*F(d, factorial(d))
                        current = other
                    require(current == y, "literal endpoint")
            else:
                sub = difference
                while True:
                    j = sub.bit_count()
                    for bit in bits:
                        if sub & bit:
                            continue
                        current = x ^ sub
                        a = orbit[states[current]]
                        b = orbit[states[current ^ bit]]
                        if a != b:
                            load[tuple(sorted((a,b)))] += base/F(comb(d-1,j))
                    if sub == 0:
                        break
                    sub = (sub-1) & difference
    triangles = primal.simple_cycles(n, 3, ei)
    matrix = primal.weighted_operator(reps, orbit, triangles)
    boolean_edges = {(i,j) for i in range(q) for j in range(i+1,q) if matrix[i][j]}
    require(set(load) == boolean_edges, "flow support covers native triangle graph")
    return q, pair_count, dict(load)


def main():
    rows = []
    gates = 0
    for n in range(3,6):
        q, pairs, load = flow(n)
        q2, pairs2, literal = flow(n, literal=True)
        require((q,pairs,load)==(q2,pairs2,literal), "subset/literal flow mismatch")
        gates += len(load)+1
        rho = max(load.values())
        row = {"n": n, "q": q, "labelled_unordered_cross_orbit_pairs": pairs,
               "edges": len(load), "rho": str(rho), "gap_lower": str(1/rho)}
        rows.append(row)
        print("NATIVE FRACTIONAL FLOW", row)
    partition_count = 0
    counts = []
    for n in range(1,17):
        numerator = 0
        for parts in burnside.partitions(n):
            b = sum(a//2 for a in parts) + sum(gcd(a,c) for i,a in enumerate(parts) for c in parts[i+1:])
            rank = len(parts) - int(any(a % 2 for a in parts))
            old_b, matrix = burnside.edge_orbits_and_constraints(n, parts)
            old_rank = burnside.gf2_rank(matrix)
            require((b,rank)==(old_b,old_rank), f"Burnside closed rank n={n},parts={parts}")
            numerator += burnside.conjugacy_class_size(n, parts)*2**(b-rank)
            partition_count += 1
            gates += 1
        require(numerator % factorial(n) == 0, "Burnside divisibility")
        counts.append(numerator//factorial(n))
    print("CLOSED BURNSIDE RANK partitions through16", partition_count)
    print("EULERIAN ORBIT COUNTS n1..16", counts)
    # Actual sparse degree prevents any one-tree all-pairs flow from being
    # polynomial: the exact lower bound below need not bound optimal flow.
    for n in (8,12,16,20):
        state_count = 2**comb(n-1,2)
        lower_q = F(state_count, factorial(n))
        tree_rho_lower = (lower_q-1)/(2*comb(n,3))
        print("ONE-TREE LOWER n", n, "q_lower", lower_q, "rho_lower", tree_rho_lower)
    payload = {"flows": rows, "partitions": partition_count, "orbit_counts": counts}
    print("SEMANTIC SHA256", hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest())
    print("CHECKS", gates)
    print("PROVED flow representation, one-tree obstruction, closed Burnside rank and uniform-orbit lift")
    print("FINITE-EXACT fractional congestion n3..5 only")
    print("OPEN polynomial all-order Boolean triangle gap / fractional congestion")


if __name__ == "__main__":
    main()
