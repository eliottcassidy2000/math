#!/usr/bin/env python3
"""Kernel-design search over every strict 2x3 cyclic topology of paid cost 10.

The first-axis slot is active (three tied horizontal rows).  If v of the
four vertical slots are active, the search uses 7-v seeds, so every candidate
has g=10,u=6, score 5/3.  Only vertical subsets containing a square are used;
forest cases are already excluded by THM-3875.
"""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement
import random
import sys

import sympy as sp

sys.path.insert(0, "04-computation")
from ak_forcing_engine import AKInstance, force

N = 6
GROUND = 6
CROSS = ((0, 3), (1, 4), (2, 5))
VERTICAL = ((0, 1), (1, 2), (3, 4), (4, 5))


def incidence(edges):
    rows = []
    for u, v, _ in edges:
        row = [0]*N
        if u != GROUND:
            row[u] += 1
        if v != GROUND:
            row[v] -= 1
        rows.append(row)
    return sp.Matrix(rows)


def make_topology(vertical_indices, seed_positions):
    # Parameter 0 is the shared horizontal slope; active verticals then get
    # one parameter each; seed rows use all remaining parameters.  Total=8.
    edges = [(u, v, 0) for u, v in CROSS]
    next_param = 1
    for idx in vertical_indices:
        edges.append((*VERTICAL[idx], next_param))
        next_param += 1
    for pos in seed_positions:
        edges.append((pos, GROUND, next_param))
        next_param += 1
    if next_param != 8 or len(edges) != 10:
        raise RuntimeError("cost/parameter accounting mismatch")
    B = incidence(edges)
    if B.rank() != N:
        return None
    cycle_basis = sp.Matrix.hstack(*B.T.nullspace()).T
    if cycle_basis.shape != (4, 10):
        raise RuntimeError("cycle dimension mismatch")
    # W=C diag(rho) B, recorded by parameter coefficient.
    coeff = [[[Fraction(0) for _ in range(8)] for _ in range(N)] for _ in range(4)]
    for i in range(4):
        for e, (_, _, param) in enumerate(edges):
            ce = Fraction(cycle_basis[i, e])
            for v in range(N):
                coeff[i][v][param] += ce * Fraction(B[e, v])
    return (tuple(vertical_indices), tuple(seed_positions), tuple(edges), B, coeff)


def topologies():
    result = []
    for vcount in (2, 3, 4):
        for vertical_indices in combinations(range(4), vcount):
            if not ({0, 2}.issubset(vertical_indices) or {1, 3}.issubset(vertical_indices)):
                continue
            seed_count = 7-vcount
            for positions in combinations_with_replacement(range(N), seed_count):
                topology = make_topology(vertical_indices, positions)
                if topology is not None:
                    result.append(topology)
    return tuple(result)


def design_matrix(coeff, kernel):
    return sp.Matrix([
        [sum(coeff[i][v][p]*Fraction(k[v]) for v in range(N)) for p in range(8)]
        for i in range(4) for k in kernel
    ])


def rank_mod(rows, columns, prime):
    a = [[row[c] % prime for c in columns] for row in rows]
    rr = 0
    for col in range(len(columns)):
        pivot = next((i for i in range(rr, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rr], a[pivot] = a[pivot], a[rr]
        inv = pow(a[rr][col], -1, prime)
        a[rr] = [(x*inv) % prime for x in a[rr]]
        for i in range(len(a)):
            if i != rr and a[i][col]:
                q = a[i][col]
                a[i] = [(x-q*y) % prime for x, y in zip(a[i], a[rr])]
        rr += 1
        if rr == len(a):
            break
    return rr


def closure_mod(topology, params, prime=1000003):
    _, _, edges, _, _ = topology
    vals = []
    for x in params:
        x = Fraction(x)
        if x.denominator % prime == 0:
            return None
        vals.append(x.numerator % prime * pow(x.denominator % prime, -1, prime) % prime)
    rows = []
    for u, v, p in edges:
        row = [0]*(2*N)
        rho = vals[p]
        if u != GROUND:
            row[2*u], row[2*u+1] = 1, rho
        if v != GROUND:
            row[2*v], row[2*v+1] = -1, -rho
        rows.append(row)
    live = set(range(N))
    trace = []
    while live:
        cols = [2*v+j for v in sorted(live) for j in (0, 1)]
        r = rank_mod(rows, cols, prime)
        fired = tuple(v for v in sorted(live)
                      if rank_mod(rows, [c for c in cols if c != 2*v+1], prime) == r-1)
        if not fired:
            break
        live.difference_update(fired)
        trace.append(fired)
    return not live, tuple(trace)


def raw(rho):
    rho = Fraction(rho)
    return (rho.denominator+rho.numerator, rho.denominator-rho.numerator)


def build(topology, params):
    vertical_indices, seed_positions, _, _, _ = topology
    labels = [raw(x) for x in params]
    fs = [{(1,): labels[0]}, {}]
    for offset, idx in enumerate(vertical_indices, 1):
        key = ((1, 1), (1, 2), (2, 1), (2, 2))[idx]
        fs[1][key] = labels[offset]
    seed_start = 1+len(vertical_indices)
    coords = ((1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3))
    seeds = [(coords[pos], labels[seed_start+i]) for i, pos in enumerate(seed_positions)]
    X = [(0, 0)] + list(dict.fromkeys(labels))
    return AKInstance(X, [2, 3], fs, [], seeds)


def main():
    rng = random.Random(20260823)
    universe = topologies()
    print(f"topologies={len(universe)}")
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    gates = rich = modular_hits = 0
    trace_counts = {}
    fire_sets = tuple(f for size in (1, 2, 3, 4) for f in combinations(range(N), size))
    for trial in range(trials):
        topology = rng.choice(universe)
        fire = set(rng.choice(fire_sets))
        rest = [v for v in range(N) if v not in fire]
        kernel_dim = rng.choice(tuple(d for d in (2, 3, 4) if d <= len(rest))) if len(rest) >= 2 else None
        if kernel_dim is None:
            continue
        while True:
            kernel = [[0]*N for _ in range(kernel_dim)]
            for v in rest:
                for j in range(kernel_dim):
                    kernel[j][v] = rng.randint(-3, 3)
            if (sp.Matrix([[kernel[j][v] for j in range(kernel_dim)] for v in rest]).rank() == kernel_dim
                    and all(any(kernel[j][v] for j in range(kernel_dim)) for v in rest)):
                break
        M = design_matrix(topology[4], tuple(kernel))
        null = M.nullspace()
        gates += 1
        if len(null) <= 1:
            continue
        rich += 1
        uniform = tuple(Fraction(1) for _ in range(8))
        basis = []
        for vector in null:
            q = tuple(Fraction(x) for x in vector)
            basis.append(tuple(q[i]-q[0]*uniform[i] for i in range(8)))
        for _ in range(8):
            weights = [rng.randint(-3, 3) for _ in basis]
            params = tuple(sum(Fraction(weights[j])*basis[j][i] for j in range(len(basis))) for i in range(8))
            if len(set(params)) == 1:
                continue
            mod = closure_mod(topology, params)
            if mod is None:
                continue
            trace_counts[mod[1]] = trace_counts.get(mod[1], 0) + 1
            if not mod[0]:
                continue
            modular_hits += 1
            inst = build(topology, params)
            ok, _, trace = force(inst, "strict")
            if ok:
                print("status=STRICT_2x3_SCORE_5_OVER_3_FOUND")
                print(f"trial={trial};vertical_indices={topology[0]};seed_positions={topology[1]}")
                print(f"prescribed_fire={tuple(sorted(fire))};kernel={tuple(tuple(k) for k in kernel)}")
                print(f"params={params};trace={tuple(tuple(x) for x in trace)}")
                print(f"score={inst.score()};m={inst.m()};r={inst.r()};fs={inst.fs};R={inst.R}")
                return
    print("status=NO_SCORE_5_OVER_3_IN_ALL_TOPOLOGY_KERNEL_SEARCH")
    print(f"gates={gates};rich={rich};modular_hits={modular_hits}")
    print(f"top_traces={tuple(sorted(trace_counts.items(), key=lambda item: (-item[1], item[0]))[:12])}")


if __name__ == "__main__":
    main()
