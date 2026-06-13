#!/usr/bin/env python3
"""
CORRECT computation of I(Omega(T), x) for the n=9 counterexample.

CRITICAL FIX: Omega(T) vertices are DIRECTED odd cycles, not vertex sets.
A vertex set of size k can support multiple directed cycles.

For a 3-cycle on {a,b,c}: exactly 1 directed cycle (up to... wait).
Actually for a 3-cycle vertex set in a tournament, there's exactly ONE
directed 3-cycle (the tournament determines the direction).

For a 5-cycle vertex set: there can be 1 or 2 directed Hamiltonian cycles.
For k vertices: up to (k-1)!/2 directed cycles, but tournament constraints limit this.

Two directed cycles are adjacent in Omega iff they share a vertex.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
import numpy as np
from itertools import combinations, permutations
import time

T = [
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
]
n = 9

def find_directed_odd_cycles(T, n, max_len=None):
    """Find all directed odd cycles. Each cycle is a tuple representing
    the cyclic order, normalized (smallest vertex first, second < last)."""
    if max_len is None:
        max_len = n
    cycles = []

    for length in range(3, max_len + 1, 2):
        for subset in combinations(range(n), length):
            # Find all directed Hamiltonian cycles on this subset
            # Fix first vertex = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                is_cycle = True
                for i in range(length):
                    if not T[cycle[i]][cycle[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize: rotate so smallest is first,
                    # then ensure cycle goes in the canonical direction
                    # (second element < last element)
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    if len(rotated) > 2 and rotated[1] > rotated[-1]:
                        rotated = (rotated[0],) + tuple(reversed(rotated[1:]))
                    cycles.append(rotated)

    # Deduplicate
    return list(set(cycles))

print("=" * 70)
print("CORRECT OMEGA(T) COMPUTATION — DIRECTED CYCLES")
print("=" * 70)

t0 = time.time()
all_cycles = find_directed_odd_cycles(T, n, 9)
t1 = time.time()
print(f"\nTotal directed odd cycles: {len(all_cycles)} (in {t1-t0:.2f}s)")

# Count by length
from collections import Counter
lengths = Counter(len(c) for c in all_cycles)
for k in sorted(lengths.keys()):
    print(f"  Length {k}: {lengths[k]} directed cycles")

# But wait — the definition says "directed odd cycles". A directed cycle
# a->b->c->a and its reverse c->b->a->c are DIFFERENT directed cycles.
# Let me re-check: in a tournament, for {a,b,c}, if a->b->c->a is a cycle,
# then a->c->b->a is NOT (since c->b would require T[c][b]=1, but T[b][c]=1).
# So for 3-cycles: each vertex set gives exactly 1 directed cycle.
# For 5-cycles: each directed Hamiltonian cycle is distinct. Its reverse is NOT
# a valid cycle (arcs are reversed). So the count is correct if we don't
# canonicalize by direction.

# Let me recount WITHOUT the direction normalization
def find_directed_odd_cycles_v2(T, n, max_len=None):
    """Find all directed odd cycles as sequences (up to rotation)."""
    if max_len is None:
        max_len = n
    cycles = set()

    for length in range(3, max_len + 1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                is_cycle = True
                for i in range(length):
                    if not T[cycle[i]][cycle[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize by rotation only (not reversal)
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)

    return list(cycles)

all_cycles_v2 = find_directed_odd_cycles_v2(T, n, 9)
print(f"\nDirected odd cycles (rotation-normalized only): {len(all_cycles_v2)}")
lengths2 = Counter(len(c) for c in all_cycles_v2)
for k in sorted(lengths2.keys()):
    print(f"  Length {k}: {lengths2[k]} directed cycles")

# The vertex set of each cycle
cycle_vsets = [frozenset(c) for c in all_cycles_v2]

# Build adjacency (share at least one vertex)
m = len(all_cycles_v2)
adj = {}
for i in range(m):
    adj[i] = frozenset(j for j in range(m) if j != i and cycle_vsets[i] & cycle_vsets[j])

degrees = [len(adj[i]) for i in range(m)]
print(f"\nOmega graph: {m} vertices, avg degree {sum(degrees)/m:.1f}")

# Independence polynomial via deletion-contraction
def independence_poly_dc(adj, vertices=None):
    memo = {}
    if vertices is None:
        vertices = frozenset(adj.keys())

    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            memo[verts] = [1]
            return [1]

        v = max(verts, key=lambda u: len(adj[u] & verts))
        remaining = verts - {v}
        p1 = solve(remaining)
        closed_nbrs = (adj[v] & verts) | {v}
        p2 = solve(verts - closed_nbrs)

        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]

        memo[verts] = result
        return result

    return solve(vertices)

print("\nComputing independence polynomial...")
t0 = time.time()
coeffs = independence_poly_dc(adj)
t1 = time.time()
print(f"  Time: {t1-t0:.2f}s")

deg = len(coeffs) - 1
while deg > 0 and coeffs[deg] == 0:
    deg -= 1

print(f"  I(Omega, x) = {' + '.join(str(coeffs[i]) + '*x^' + str(i) for i in range(deg+1) if coeffs[i] > 0)}")
print(f"  I(Omega, 2) = {sum(coeffs[i] * 2**i for i in range(deg+1))}")

# Verify H(T)
h = 0
for perm in permutations(range(n)):
    if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
        h += 1
print(f"  H(T) = {h}")
print(f"  Match: {sum(coeffs[i] * 2**i for i in range(deg+1)) == h}")

# Real-rootedness
if deg > 0:
    poly_c = [coeffs[deg - i] for i in range(deg + 1)]
    roots = np.roots(poly_c)
    all_real = all(abs(r.imag) < 1e-8 for r in roots)
    print(f"\n  Roots: {roots}")
    print(f"  All real: {all_real}")

    # Newton's inequalities
    for k in range(1, deg):
        lhs = coeffs[k]**2
        rhs = coeffs[k-1] * coeffs[k+1] * (k+1) / k
        print(f"  Newton k={k}: {coeffs[k]}^2={lhs} vs {coeffs[k-1]}*{coeffs[k+1]}*{(k+1)/k:.2f}={rhs:.1f}: {'OK' if lhs >= rhs else 'FAIL'}")
