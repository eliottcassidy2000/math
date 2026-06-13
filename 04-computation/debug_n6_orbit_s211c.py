#!/usr/bin/env python3
"""
debug_n6_orbit_s211c.py -- opus-2026-03-23-S211c

Find exactly which class has the D discrepancy at n=6.
"""

import sys
import numpy as np
from math import comb
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def canonical(A, n):
    scores = [sum(A[i]) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n): score_groups[scores[v]].append(v)
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

n = 6; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

cmap = {}; sizes = {}; reps = {}; tc = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canonical(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1
    tc[bits] = cmap[c]
nc = len(cmap)

# Build F matrix
F = np.zeros((nc, nc), dtype=int)
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m):
        F[ci][tc[bits^(1<<k)]] += 1

# For each class: compute D two ways
# Way A: from orbit counting
# Way B: from F-matrix multiplicities

for cid in range(nc):
    A = reps[cid]
    sz = sizes[cid]

    # Automorphisms
    auts = []
    for perm in permutations(range(n)):
        if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j):
            auts.append(perm)

    # Arc orbits
    arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
    seen = set(); orbits = []
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            orbit.add((perm[arc[0]], perm[arc[1]]))
        orbits.append(frozenset(orbit))
        seen.update(orbit)

    # Orbit targets
    targets_A = []
    for orbit in orbits:
        e = min(orbit)  # deterministic choice
        A2 = [row[:] for row in A]
        A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
        targets_A.append(cmap[canonical(A2, n)])

    mult_A = Counter(targets_A)
    D_A = sum(comb(v, 2) for v in mult_A.values())

    # From F: mult(C->D) should be orbits targeting D
    # F[C][D] counts (labeled tournament, arc) pairs going C->D
    # F[C][D] = sz * (arcs per T targeting D)
    # arcs per T targeting D = sum of orbit sizes for orbits targeting D
    # orbits targeting D = mult_A[D]
    # But we want to also verify using F

    # F-based orbit count (arcs_per_T / orbit_size doesn't work uniformly)
    # Instead: total distinct D targets from F = #{D : F[C][D] > 0}
    n_targets_F = sum(1 for D in range(nc) if F[cid][D] > 0)
    n_targets_A = len(mult_A)

    # Check consistency: T_class = # orbits = len(orbits)
    n_orbits = len(orbits)
    T_from_F = sum(F[cid][D] // sz for D in range(nc) if F[cid][D] > 0)
    # Wait: T_from_F should be m (arcs per tournament), not orbits.
    # Actually F[cid][D]/sz = arcs per T in class cid that go to D.
    # sum_D F[cid][D]/sz = m (total arcs per T).
    arcs_per_T = sum(F[cid][D] for D in range(nc)) // sz
    # Orbits: m / avg_orbit_size. Orbit sizes = |aut(C)| * (orbit type).
    # Actually just use Burnside: #orbits = (1/|Aut|) * sum_sigma #fixed_arcs(sigma)

    if D_A != 0:
        # Compare with F-based D: need to figure out orbit structure from F
        # From F: arcs_to_D = F[C][D]/sz for each D
        # If all orbits have size |Aut|, then orbits_to_D = arcs_to_D / |Aut|
        # But orbits can be smaller if arc is fixed by some automorphisms

        # Actually: orbits_to_D = mult_A[D] from our computation
        # F[C][D]/sz = sum_{orbits O targeting D} |O|
        # |O| = |Aut| / |Stab(e)| for any e in O

        # Let's just check if D_A matches what it should be
        # D_total_check from known E=290: D = T - 2E = 704 - 580 = 124
        pass

    if n_orbits != m // len(auts) * len(auts) + (m % len(auts) != 0):
        pass  # Orbits don't have to all be the same size

    # Print classes with any degeneracy
    if D_A > 0:
        print("Class %d (|Aut|=%d, size=%d, H=%d): %d orbits, D=%d" %
              (cid, len(auts), sz, 0, n_orbits, D_A))
        print("  mult: %s" % dict(mult_A))
        # Show orbit sizes
        osizes = [len(o) for o in orbits]
        print("  orbit sizes: %s" % Counter(osizes))

# Now compute total D
print("\nTotal D from orbit method:")
D_total = 0
for cid in range(nc):
    A = reps[cid]
    auts = []
    for perm in permutations(range(n)):
        if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j):
            auts.append(perm)
    arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
    seen = set(); orbits = []
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            orbit.add((perm[arc[0]], perm[arc[1]]))
        orbits.append(orbit)
        seen.update(orbit)
    targets = []
    for orbit in orbits:
        e = min(orbit)
        A2 = [row[:] for row in A]
        A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
        targets.append(cmap[canonical(A2, n)])
    mult = Counter(targets)
    d = sum(comb(v, 2) for v in mult.values())
    D_total += d

print("D_total = %d (should be 124, T-2E = 704-580 = 124)" % D_total)
print("E = (T-D)/2 = (704-%d)/2 = %d" % (D_total, (704-D_total)//2))

# Also check: what does the orbit computation give for arcs not
# mapping back to the right adjacency?
# Is there an orbit where the representative arc reversal gives wrong target?

# Let's verify: for classes with |Aut|>1, check each orbit
print("\nVerifying orbits for classes with |Aut| > 1:")
for cid in range(nc):
    A = reps[cid]
    auts = []
    for perm in permutations(range(n)):
        if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j):
            auts.append(perm)
    if len(auts) <= 1: continue

    arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
    seen = set(); orbits = []
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            orbit.add((perm[arc[0]], perm[arc[1]]))
        orbits.append(orbit)
        seen.update(orbit)

    # For each orbit, verify ALL elements give same target
    for oi, orbit in enumerate(orbits):
        target_set = set()
        for e in orbit:
            A2 = [row[:] for row in A]
            A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
            target_set.add(cmap[canonical(A2, n)])
        if len(target_set) > 1:
            print("  BUG: class %d orbit %d has %d distinct targets: %s" %
                  (cid, oi, len(target_set), target_set))

print("Verification complete.")
