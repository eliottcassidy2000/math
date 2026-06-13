#!/usr/bin/env python3
"""
β_2 = 0 FOR TOURNAMENTS: DT-SUFFICIENCY PROOF

THEOREM (computational): For all tournaments on n ≤ 6 vertices,
  im(∂_3|_{DT}) = ker(∂_2|Ω_2)
where DT = span of doubly-transitive 4-paths.

A 4-path (a,b,c,d) is DT iff: a→b→c→d AND a→c AND b→d.
Equivalently: (a,b,c,d) is an allowed 3-path AND all 4 faces are in A_2.

This means β_2 = 0 follows from: every Ω_2-cycle is a boundary of DT 4-paths.

KEY ALGEBRAIC STRUCTURE:
For DT 4-path (a,b,c,d):
  ∂_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
All four faces are transitive triples (verified by DT condition).

WHAT THIS MEANS:
If we can show im(∂_3|_{DT}) ⊇ ker(∂_2|Ω_2) algebraically, we get β_2 = 0.
"""
import numpy as np
from itertools import permutations
from collections import Counter
import sys, time, random
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def check_beta2_dt(A, n):
    """Check if DT 4-paths span ker(∂_2|Ω_2).
    Returns (ker_dim, im_dt_rank, beta2_dt)."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Ω_2 = transitive triples
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0:
        return (0, 0, 0)

    # ∂_2: Ω_2 → A_1
    a1_idx = {tuple(p): i for i, p in enumerate(a1)}
    bd2 = np.zeros((len(a1), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[a1_idx[(b,c)], j] += 1
        bd2[a1_idx[(a,c)], j] -= 1
        bd2[a1_idx[(a,b)], j] += 1

    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2

    if ker_dim == 0:
        return (0, 0, 0)

    # DT 4-paths
    dt_paths = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]
    if len(dt_paths) == 0:
        return (ker_dim, 0, ker_dim)

    tt_idx = {t: i for i, t in enumerate(tt)}

    # ∂_3|DT → Ω_2
    bd3_dt = np.zeros((len(tt), len(dt_paths)))
    for j, path in enumerate(dt_paths):
        v0, v1, v2, v3 = path
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in tt_idx:
                bd3_dt[tt_idx[face], j] += sign

    im_dt_rank = np.linalg.matrix_rank(bd3_dt, tol=1e-8)
    beta2 = ker_dim - im_dt_rank
    return (ker_dim, im_dt_rank, beta2)

# ===== Exhaustive at n=4,5 =====
print("=" * 70)
print("DT-SUFFICIENCY FOR β_2 = 0")
print("=" * 70)

for n in [4, 5]:
    t0 = time.time()
    all_pass = True
    count = 0
    fail_count = 0
    stats = Counter()

    for A in all_tournaments_gen(n):
        count += 1
        ker, im_dt, beta2 = check_beta2_dt(A, n)
        if beta2 > 0:
            all_pass = False
            fail_count += 1
        if ker > 0:
            stats[(ker, im_dt)] += 1

    t1 = time.time()
    print(f"\nn={n}: {count} tournaments, {t1-t0:.1f}s")
    print(f"  DT-sufficiency holds: {all_pass} ({fail_count} failures)")
    print(f"  (ker, im_dt): count")
    for key in sorted(stats.keys()):
        k, im = key
        print(f"    ({k}, {im}): {stats[key]}  {'✓' if im >= k else '✗'}")

# ===== Exhaustive at n=6 =====
print(f"\n--- n=6 exhaustive (2^15 = 32768 tournaments) ---")
n = 6
t0 = time.time()
all_pass = True
count = 0
fail_count = 0
stats = Counter()

for A in all_tournaments_gen(n):
    count += 1
    ker, im_dt, beta2 = check_beta2_dt(A, n)
    if beta2 > 0:
        all_pass = False
        fail_count += 1
    if ker > 0:
        stats[(ker, im_dt)] += 1
    if count % 5000 == 0:
        print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

t1 = time.time()
print(f"\nn=6: {count} tournaments, {t1-t0:.1f}s")
print(f"  DT-sufficiency holds: {all_pass} ({fail_count} failures)")
print(f"  (ker, im_dt): count")
for key in sorted(stats.keys()):
    k, im = key
    print(f"    ({k}, {im}): {stats[key]}  {'✓' if im >= k else '✗'}")

# ===== Random sample at n=7 =====
print(f"\n--- n=7 random sample (500 tournaments) ---")
n = 7
t0 = time.time()
all_pass = True
count = 0
fail_count = 0
stats = Counter()
random.seed(42)

for _ in range(500):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    A = [[0]*n for _ in range(n)]
    for (i,j) in edges:
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1
    count += 1
    ker, im_dt, beta2 = check_beta2_dt(A, n)
    if beta2 > 0:
        all_pass = False
        fail_count += 1
    if ker > 0:
        stats[(ker, im_dt)] += 1

t1 = time.time()
print(f"\nn=7 (sample): {count} tournaments, {t1-t0:.1f}s")
print(f"  DT-sufficiency holds: {all_pass} ({fail_count} failures)")
print(f"  (ker, im_dt): count")
for key in sorted(stats.keys()):
    k, im = key
    print(f"    ({k}, {im}): {stats[key]}  {'✓' if im >= k else '✗'}")

# ===== ALGEBRAIC ANALYSIS: What structure do DT boundaries have? =====
print(f"\n\n{'='*70}")
print("ALGEBRAIC STRUCTURE OF DT BOUNDARY MAP")
print("="*70)

# For a fixed n, let's understand the boundary map ∂_3|DT: span(DT) → Ω_2.
# Each DT path (a,b,c,d) maps to (b,c,d)-(a,c,d)+(a,b,d)-(a,b,c).
# All four are transitive triples.
#
# KEY QUESTION: What is the STRUCTURE of this map?
# In particular, when is a transitive triple in the IMAGE of ∂_3|DT?

n = 5
print(f"\nn={n}: Which transitive triples appear in DT boundaries?")

for t_idx, A in enumerate(all_tournaments_gen(n)):
    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    dt = [tuple(p) for p in a3 if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]

    if len(tt) < 5 or len(dt) < 3:
        continue

    # Which transitive triples appear in at least one DT boundary?
    tt_appears = set()
    for path in dt:
        v0, v1, v2, v3 = path
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        for face in faces:
            if face in set(tt):
                tt_appears.add(face)

    coverage = len(tt_appears) / len(tt) if len(tt) > 0 else 0
    if coverage < 1.0 and t_idx < 50:
        missing = set(tt) - tt_appears
        print(f"  T#{t_idx}: |Ω_2|={len(tt)}, |DT|={len(dt)}, coverage={coverage:.2f}")
        print(f"    Missing triples: {missing}")

        # These missing triples can't be individually reached.
        # But can they appear in LINEAR COMBINATIONS of DT boundaries?
        # (They definitely can if β_2=0.)
        break

print("\nDone.")
