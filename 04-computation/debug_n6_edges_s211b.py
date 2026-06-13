#!/usr/bin/env python3
"""
debug_n6_edges_s211b.py -- opus-2026-03-23-S211b

Debug the n=6 edge count discrepancy (291 vs 290).
Compute edges via THREE independent methods:
1. Orbit-based (current S211 approach)
2. F-matrix adjacency (F[C][D] > 0)
3. Direct labeled enumeration

Author: opus-2026-03-23-S211b
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            c = dp[(mask,v)]
            if c == 0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += c
    return sum(dp[(full,v)] for v in range(n))

def canonical(A, n):
    """Canonical form using score-based pruning for n=6."""
    scores = [sum(A[i]) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n):
        score_groups[scores[v]].append(v)
    groups = [score_groups[s] for s in sorted(score_groups.keys())]

    best = None
    def gen_perms(grps):
        if not grps:
            yield []
            return
        for p in permutations(grps[0]):
            for rest in gen_perms(grps[1:]):
                yield list(p) + rest

    for perm in gen_perms(groups):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

n = 6
m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
print("n=%d, m=%d, 2^m=%d" % (n, m, 2**m))

# Build iso classes
cmap = {}; sizes = {}; reps = {}; H_vals = {}
tc = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canonical(A, n)
    if c not in cmap:
        cid = len(cmap)
        cmap[c] = cid
        sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
        H_vals[cid] = count_hp(A, n)
    sizes[cmap[c]] += 1
    tc[bits] = cmap[c]

nc = len(cmap)
print("nc = %d classes" % nc)

# METHOD 1: F-matrix adjacency
F = np.zeros((nc, nc), dtype=int)
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m):
        bits2 = bits ^ (1<<k)
        cj = tc[bits2]
        F[ci][cj] += 1

# Adjacency: F[C][D] > 0 means C and D are adjacent
adj = (F > 0).astype(int)
np.fill_diagonal(adj, 0)
E_adj = np.sum(adj) // 2  # undirected edges (excluding self-loops)
# Self-loops: F[C][C] > 0
self_loop_classes = sum(1 for c in range(nc) if F[c][c] > 0)
print("\nMETHOD 1 (F-matrix adjacency):")
print("  Undirected edges (no self-loops): %d" % E_adj)
print("  Classes with self-loops: %d" % self_loop_classes)

# METHOD 2: Direct labeled enumeration
# For each pair (bits1, bits2) differing in exactly 1 bit,
# they are in different classes or same class.
# Edge C-D exists iff some labeled tournament in C has an arc reversal in D.
edges_set = set()
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m):
        bits2 = bits ^ (1<<k)
        cj = tc[bits2]
        if ci != cj:
            a, b = min(ci, cj), max(ci, cj)
            edges_set.add((a, b))
E_direct = len(edges_set)
print("\nMETHOD 2 (direct labeled enumeration):")
print("  Undirected edges (no self-loops): %d" % E_direct)

# METHOD 3: Orbit-based
# For each class, find arc-orbits and their targets
print("\nMETHOD 3 (orbit-based):")
T_total = 0; D_total = 0
orbit_edges = set()
for cid in range(nc):
    A = reps[cid]

    # Find automorphisms
    auts = []
    for perm in permutations(range(n)):
        ok = all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if ok: auts.append(perm)

    # Arc orbits
    arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
    seen = set(); orbits = []
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            mapped = (perm[arc[0]], perm[arc[1]])
            orbit.add(mapped)
        orbits.append(orbit)
        seen.update(orbit)

    T_total += len(orbits)

    # Find target of each orbit
    targets = []
    for orbit in orbits:
        e = next(iter(orbit))
        A2 = [row[:] for row in A]
        A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
        target = cmap[canonical(A2, n)]
        targets.append(target)
        if cid != target:
            a, b = min(cid, target), max(cid, target)
            orbit_edges.add((a, b))

    mult = Counter(targets)
    D_total += sum(comb(v, 2) for v in mult.values())

    # Verify: m arcs total should equal sum of orbit sizes
    total_arcs = sum(len(o) for o in orbits)
    if total_arcs != m:
        print("  BUG: class %d has %d arcs in orbits, expected %d" % (cid, total_arcs, m))

E_orbit = len(orbit_edges)
E_formula = (T_total - D_total) // 2
print("  T = %d (should be %d)" % (T_total, 704))
print("  D = %d" % D_total)
print("  E = (T-D)/2 = %d" % E_formula)
print("  Distinct edges from orbits = %d" % E_orbit)

# Compare all three
print("\nCOMPARISON:")
print("  Method 1 (F-matrix adj): E = %d" % E_adj)
print("  Method 2 (direct enum):  E = %d" % E_direct)
print("  Method 3 (orbit formula): E = %d" % E_formula)
print("  Method 3 (orbit set):    E = %d" % E_orbit)

# Check: are methods 1 and 2 the same?
if E_adj == E_direct:
    print("\n  Methods 1,2 AGREE: E = %d" % E_adj)
    if E_formula != E_adj:
        print("  Method 3 DISAGREES by %d" % (E_formula - E_adj))
        # Find the problematic class(es)
        print("\n  Investigating orbit-formula discrepancy...")

        # For each class, check if mult computation gives correct edge count
        for cid in range(nc):
            A = reps[cid]
            auts = []
            for perm in permutations(range(n)):
                ok = all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
                if ok: auts.append(perm)

            arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
            seen = set(); orbits = []
            for arc in arcs:
                if arc in seen: continue
                orbit = set()
                for perm in auts:
                    mapped = (perm[arc[0]], perm[arc[1]])
                    orbit.add(mapped)
                orbits.append(orbit)
                seen.update(orbit)

            targets = []
            for orbit in orbits:
                e = next(iter(orbit))
                A2 = [row[:] for row in A]
                A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
                target = cmap[canonical(A2, n)]
                targets.append(target)

            mult = Counter(targets)
            n_orbits = len(orbits)
            s2_c = sum(v**2 for v in mult.values())
            d_c = sum(comb(v, 2) for v in mult.values())

            # From F matrix: arcs per tournament targeting D
            # F[cid][D] / sizes[cid] should = sum of orbit sizes targeting D / 1
            # But orbit sizes vary!
            for D_target, mult_val in mult.items():
                # Count arcs (not orbits) targeting D
                arcs_to_D = sum(len(orbits[i]) for i, t in enumerate(targets) if t == D_target)
                f_val = F[cid][D_target]
                f_per_T = f_val // sizes[cid]
                if f_per_T != arcs_to_D:
                    print("  Class %d -> %d: orbits=%d, arcs=%d, F/size=%d, F=%d, size=%d" %
                          (cid, D_target, mult_val, arcs_to_D, f_per_T, f_val, sizes[cid]))

print("\nDone.")
