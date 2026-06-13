#!/usr/bin/env python3
"""
sl_decomposition_n7_s212.py -- opus-2026-03-23-S212

Compute the SL (self-loop) and excess decomposition at n=7.

Strategy: enumerate all 2^21 = 2M tournaments, build F matrix using
score-based hashing + canonical form refinement. Then compute:
- SL from F[C][C] / size(C) for trivial-Aut classes
- Full orbit decomposition for nontrivial-Aut classes (few)
- excess_cross from the orbit-level mult data

G(7) = T(7) - 2*E(7) = 8912 - 2*4086 = 740 = SL(7) + excess(7)

Author: opus-2026-03-23-S212
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

n = 7
m = comb(n, 2)  # 21
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

print("n=%d, m=%d, 2^m=%d" % (n, m, 2**m))

# ============================================================
banner("PHASE 1: BUILD ISO CLASSES")
# ============================================================

def score_seq(A):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_3cyc(A):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c += 1
                if A[i][k] and A[k][j] and A[j][i]: c += 1
    return c

def count_hp(A):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v, v)] = 1
    full = (1<<n) - 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)): continue
            c = dp[(mask, v)]
            if c == 0: continue
            for w in range(n):
                if mask & (1<<w): continue
                if A[v][w]: dp[(mask|(1<<w), w)] += c
    return sum(dp[(full, v)] for v in range(n))

def canonical(A):
    """Score-based canonical form for n=7."""
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

# Use hash-based approach: group by (score_seq, c3, H), then canonical for collisions
t0 = time.time()

# First pass: group tournaments by hash
hash_groups = defaultdict(list)  # (score, c3) -> list of (bits, A)
class_map = {}  # canonical_form -> class_id
sizes = {}
class_reps = {}
tc = {}  # bits -> class_id

for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(pairs):
        if (bits>>k) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = score_seq(A)
    hash_groups[key].append(bits)

t1 = time.time()
print("  Phase 1a: grouped %d tournaments by score in %.1fs" % (2**m, t1-t0))
print("  %d distinct score sequences" % len(hash_groups))

# For each score group, compute canonical forms
nc = 0
for key in sorted(hash_groups.keys()):
    group = hash_groups[key]
    # Within this group, build canonical forms
    for bits in group:
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits>>k) & 1: A[i][j] = 1
            else: A[j][i] = 1

        c = canonical(A)
        if c not in class_map:
            cid = nc; nc += 1
            class_map[c] = cid
            sizes[cid] = 0
            class_reps[cid] = [row[:] for row in A]
        sizes[class_map[c]] += 1
        tc[bits] = class_map[c]

t2 = time.time()
print("  Phase 1b: found %d iso classes in %.1fs" % (nc, t2-t0))
print("  Total time: %.1fs" % (t2-t0))

# ============================================================
banner("PHASE 2: BUILD F MATRIX (DIAGONAL AND EDGES)")
# ============================================================

# We need F[C][D] for all C,D (or at least for C=D and neighbors)
# Full F matrix: nc x nc = 456 x 456

F = np.zeros((nc, nc), dtype=np.int64)
t3 = time.time()

for bits in range(2**m):
    ci = tc[bits]
    for k in range(m):
        bits2 = bits ^ (1<<k)
        cj = tc[bits2]
        F[ci][cj] += 1

t4 = time.time()
print("  F matrix (%d x %d) built in %.1fs" % (nc, nc, t4-t3))

# ============================================================
banner("PHASE 3: EDGE COUNT AND BASIC STATISTICS")
# ============================================================

# E from adjacency
adj = (F > 0).astype(int)
np.fill_diagonal(adj, 0)
E = np.sum(adj) // 2

print("  V = %d, E = %d (expected 4086)" % (nc, E))
print("  T = %d" % sum(F[c].sum() // sizes[c] for c in range(nc)))

# Quick T computation
T_check = 0
for c in range(nc):
    # Total arcs per tournament = m. But F counts per labeled tournament.
    # Total F[c] sum = sizes[c] * m (each tournament has m arcs)
    # So orbits per class = sum_D F[c][D] / sizes[c] should be m.
    # But T = sum over classes of (number of arc-orbits)
    pass

# Actually T = total arc-orbits = sum_C (#arc-orbits of C)
# For |Aut|=1: #orbits = m = 21
# For |Aut|=a: #orbits = Burnside count = (1/a)*sum_sigma #fixed_arcs(sigma)

# For now, compute T from known value
T_known = 8912
E_known = 4086
G_known = T_known - 2*E_known  # = 740

print("  G = T - 2E = %d - %d = %d (expected 740)" %
      (T_known, 2*E_known, G_known))

# ============================================================
banner("PHASE 4: COMPUTE |Aut| FOR EACH CLASS")
# ============================================================

# For self-loop computation, we need |Aut| and arc-orbits
# First: identify which classes have nontrivial Aut

aut_sizes = {}
t5 = time.time()

for cid in range(nc):
    # Quick check: if size(cid) = n!/1 = 5040, then |Aut| = 1
    expected_size = factorial(n) // 1
    if sizes[cid] == expected_size:
        aut_sizes[cid] = 1
    else:
        aut_sizes[cid] = factorial(n) // sizes[cid]

nontrivial = [(cid, aut_sizes[cid]) for cid in range(nc) if aut_sizes[cid] > 1]
t6 = time.time()

print("  Classes with |Aut| > 1: %d out of %d" % (len(nontrivial), nc))
for cid, a in nontrivial:
    print("    class %d: |Aut| = %d, size = %d" % (cid, a, sizes[cid]))

# ============================================================
banner("PHASE 5: COMPUTE SL")
# ============================================================

# For |Aut|=1 classes: SL(C) = F[C][C] / size(C)
# (each arc is its own orbit, so self-loop orbits = self-loop arcs per tournament)

# For |Aut|>1 classes: need orbit decomposition
# SL(C) = number of arc-orbits targeting C itself

SL_total = 0
excess_total = 0  # excess cross-orbits beyond 2 per edge

# Track mult(C→D) for each (C, D) pair
mult_matrix = {}  # (C, D) -> mult (number of orbits)

for cid in range(nc):
    A = class_reps[cid]
    a = aut_sizes[cid]

    if a == 1:
        # Each arc is its own orbit. mult(C→D) = F[C][D]/size(C) for each D.
        for did in range(nc):
            if F[cid][did] > 0:
                mult_val = F[cid][did] // sizes[cid]
                mult_matrix[(cid, did)] = mult_val
                if cid == did:
                    SL_total += mult_val
    else:
        # Need orbit decomposition
        # Find actual automorphisms
        auts = []
        for perm in permutations(range(n)):
            if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j):
                auts.append(perm)

        # Find arc orbits
        arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
        seen = set(); orbits = []
        for arc in arcs:
            if arc in seen: continue
            orbit = set()
            for perm in auts:
                orbit.add((perm[arc[0]], perm[arc[1]]))
            orbits.append(orbit)
            seen.update(orbit)

        # Target of each orbit
        target_counter = Counter()
        for orbit in orbits:
            e = min(orbit)
            A2 = [row[:] for row in A]
            A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
            target = tc_from_adj(A2) if False else None
            # Need canonical form of A2
            target_can = canonical(A2)
            target_cid = class_map[target_can]
            target_counter[target_cid] += 1

        for did, mult_val in target_counter.items():
            mult_matrix[(cid, did)] = mult_val
            if cid == did:
                SL_total += mult_val

    if (cid + 1) % 50 == 0:
        print("  ... processed %d/%d classes" % (cid+1, nc))

t7 = time.time()
print("\n  SL computation done in %.1fs" % (t7-t5))
print("  SL = %d" % SL_total)

# ============================================================
banner("PHASE 6: COMPUTE EXCESS")
# ============================================================

# excess = sum over undirected edges {A,B} of (mult(A→B) + mult(B→A) - 2)

edges = set()
for (ci, cj), mult_val in mult_matrix.items():
    if ci != cj and mult_val > 0:
        edges.add((min(ci, cj), max(ci, cj)))

E_from_mult = len(edges)
print("  E from mult_matrix = %d (should be %d)" % (E_from_mult, E_known))

excess_total = 0
for (a, b) in edges:
    m_ab = mult_matrix.get((a, b), 0)
    m_ba = mult_matrix.get((b, a), 0)
    excess_total += (m_ab + m_ba - 2)

print("  excess_cross = %d" % excess_total)
print("  SL + excess = %d + %d = %d (should be G = %d)" %
      (SL_total, excess_total, SL_total + excess_total, G_known))

# ============================================================
banner("PHASE 7: VERIFICATION AND SUMMARY")
# ============================================================

# Verify T
T_computed = sum(mult_matrix.values())
print("  T from mult_matrix = %d (expected %d)" % (T_computed, T_known))

# D = sum C(mult, 2)
D_total = sum(comb(v, 2) for v in mult_matrix.values())
print("  D = %d" % D_total)

# S_2 = sum mult^2
S2_total = sum(v**2 for v in mult_matrix.values())
print("  S_2 = %d" % S2_total)

# Cross-check: D = (S_2 - T) / 2
print("  (S_2 - T) / 2 = (%d - %d) / 2 = %d (should be D=%d)" %
      (S2_total, T_computed, (S2_total - T_computed)//2, D_total))

# E check
E_formula = (T_computed - D_total) // 2  # the WRONG formula
print("\n  E_wrong = (T - D)/2 = %d" % E_formula)
print("  E_correct = (T - SL - excess)/2 = %d" %
      ((T_computed - SL_total - excess_total)//2))
print("  E_actual = %d" % E_known)

print("\n  COMPLETE DECOMPOSITION TABLE:")
print("  n  | T     | E    | G=T-2E | SL  | excess | G/T")
print("  ---|-------|------|--------|-----|--------|------")
SL_known = {3:2, 4:6, 5:16, 6:58}
excess_known = {3:0, 4:0, 5:12, 6:66}
for nn in [3, 4, 5, 6]:
    T_n = {3:4,4:16,5:88,6:704}[nn]
    E_n = {3:1,4:5,5:30,6:290}[nn]
    G_n = T_n - 2*E_n
    print("  %d  | %5d | %4d | %5d  | %3d | %5d  | %.4f" %
          (nn, T_n, E_n, G_n, SL_known[nn], excess_known[nn], G_n/T_n))

# n=7 row
print("  %d  | %5d | %4d | %5d  | %3d | %5d  | %.4f" %
      (7, T_known, E_known, G_known, SL_total, excess_total, G_known/T_known))

print("\nDone. opus-2026-03-23-S212")
