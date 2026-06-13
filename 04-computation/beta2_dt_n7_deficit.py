#!/usr/bin/env python3
"""
beta2_dt_n7_deficit.py - DT filling deficit at n=7,8 (sampled)

Results so far:
  n=5: DT fills ALL Z_2 (deficit = 0 always)
  n=6: DT deficit <= 2 (0: 97.1%, 1: 2.2%, 2: 0.7%)

Questions:
  1. What is the max deficit at n=7?
  2. Does DT+cancellation still fill everything?
  3. What tournaments have the highest deficit?
  4. Is there a formula for |DT| in terms of score sequence?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def build_adj_rand(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def analyze_dt_filling(A, n):
    """Compute DT filling analysis."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if d_om2 == 0:
        return {'dim_Z2': 0, 'deficit': 0, 'cancel_ok': True}

    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))

    # Z_2
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    dim_Z2 = d_om2 - rk2

    if dim_Z2 == 0 or not a3:
        return {'dim_Z2': dim_Z2, 'deficit': 0, 'cancel_ok': True,
                'n_dt': 0, 'dim_Om3': 0, 'n_A3': 0}

    bd3 = build_full_boundary_matrix(a3, a2)

    # DT
    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

    if dt_idx:
        V_dt = np.zeros((len(a3), len(dt_idx)))
        for j, idx in enumerate(dt_idx):
            V_dt[idx, j] = 1
        bd3_dt = bd3 @ V_dt
        coords_dt, _, _, _ = np.linalg.lstsq(om2, bd3_dt, rcond=None)
        rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)
    else:
        rk_dt = 0

    deficit = dim_Z2 - rk_dt

    # Cancellation pairs
    cancel_ok = True
    if deficit > 0:
        bad02 = defaultdict(list)
        bad13 = defaultdict(list)
        for i, p in enumerate(a3):
            a_, b_, c_, d_ = p
            if not A[a_][c_]:
                bad02[(a_, c_, d_)].append(i)
            if not A[b_][d_]:
                bad13[(a_, b_, d_)].append(i)

        all_vecs = []
        for idx in dt_idx:
            v = np.zeros(len(a3))
            v[idx] = 1
            all_vecs.append(v)
        for group_key, indices in bad02.items():
            if len(indices) >= 2:
                for j in range(1, len(indices)):
                    v = np.zeros(len(a3))
                    v[indices[0]] = 1
                    v[indices[j]] = -1
                    all_vecs.append(v)
        for group_key, indices in bad13.items():
            if len(indices) >= 2:
                for j in range(1, len(indices)):
                    v = np.zeros(len(a3))
                    v[indices[0]] = 1
                    v[indices[j]] = -1
                    all_vecs.append(v)

        if all_vecs:
            V_all = np.column_stack(all_vecs)
            bd3_all = bd3 @ V_all
            coords_all, _, _, _ = np.linalg.lstsq(om2, bd3_all, rcond=None)
            rk_all = np.linalg.matrix_rank(coords_all, tol=1e-8)
            cancel_ok = rk_all >= dim_Z2
        else:
            cancel_ok = False

    # Omega_3
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    return {
        'dim_Z2': dim_Z2, 'deficit': deficit, 'cancel_ok': cancel_ok,
        'n_dt': len(dt_idx), 'dim_Om3': d_om3, 'n_A3': len(a3),
        'rk_dt': rk_dt,
    }


# ============================================================
# TEST 1: n=7 sampled
# ============================================================
print("=" * 70)
print("DT FILLING DEFICIT AT n=7 (sampled)")
print("=" * 70)

random.seed(42)
rng = random.Random(42)
n = 7
samples = 3000

deficit_dist = Counter()
cancel_fails = 0
max_deficit = 0
deficit_examples = []

t0 = time.time()
for s in range(samples):
    if s % 500 == 0 and s > 0:
        dt = time.time() - t0
        print(f"  {s}/{samples} ({dt:.0f}s), max_deficit={max_deficit}, cancel_fails={cancel_fails}")

    A = build_adj_rand(n, rng)
    res = analyze_dt_filling(A, n)

    deficit = res['deficit']
    deficit_dist[deficit] += 1

    if deficit > max_deficit:
        max_deficit = deficit
        # Save example
        scores = tuple(sorted([sum(row) for row in A]))
        deficit_examples.append((scores, deficit, res))

    if not res.get('cancel_ok', True):
        cancel_fails += 1

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")
print(f"DT deficit distribution: {dict(sorted(deficit_dist.items()))}")
print(f"Max deficit: {max_deficit}")
print(f"Cancel fails: {cancel_fails}")

if cancel_fails == 0:
    print("DT + cancellation fills ALL Z_2!")

# Show high-deficit examples
print("\nHigh-deficit examples:")
for scores, deficit, res in deficit_examples[-5:]:
    print(f"  scores={scores}, deficit={deficit}, Z2={res['dim_Z2']}, "
          f"DT={res['n_dt']}, rk_DT={res.get('rk_dt', '?')}, Om3={res['dim_Om3']}")


# ============================================================
# TEST 2: n=8 sampled (smaller sample)
# ============================================================
print(f"\n{'='*70}")
print("DT FILLING DEFICIT AT n=8 (sampled)")
print("=" * 70)

rng = random.Random(123)
n = 8
samples = 1000

deficit_dist_8 = Counter()
cancel_fails_8 = 0
max_deficit_8 = 0

t0 = time.time()
for s in range(samples):
    if s % 200 == 0 and s > 0:
        dt = time.time() - t0
        print(f"  {s}/{samples} ({dt:.0f}s), max_deficit={max_deficit_8}, cancel_fails={cancel_fails_8}")

    A = build_adj_rand(n, rng)
    res = analyze_dt_filling(A, n)

    deficit = res['deficit']
    deficit_dist_8[deficit] += 1
    if deficit > max_deficit_8:
        max_deficit_8 = deficit
    if not res.get('cancel_ok', True):
        cancel_fails_8 += 1

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")
print(f"DT deficit distribution: {dict(sorted(deficit_dist_8.items()))}")
print(f"Max deficit: {max_deficit_8}")
print(f"Cancel fails: {cancel_fails_8}")

if cancel_fails_8 == 0:
    print("DT + cancellation fills ALL Z_2 at n=8!")


# ============================================================
# ANALYSIS: DT count formula
# ============================================================
print(f"\n{'='*70}")
print("DT COUNT FORMULA")
print("=" * 70)

# A DT 3-path (a,b,c,d) requires:
# a->b, b->c, c->d (allowed 3-path)
# a->c (face 02 is transitive)
# b->d (face 13 is transitive)
#
# This means: in the subtournament on {a,b,c,d}:
# a,b,c form a transitive triple (a->b->c, a->c)
# b,c,d form a transitive triple (b->c->d, b->d)
# So: a->b->c->d with a->c, b->d
# The arc between a and d can go either way.
#
# Count of DT paths = sum over all 4-vertex sets of (DT paths in that set)
#
# For a 4-vertex set {w,x,y,z} that forms a transitive tournament w>x>y>z:
# The only DT path is (w,x,y,z). So contribution = 1.
#
# For a 4-vertex set with one 3-cycle:
# Score sequence is some permutation of (0,1,2,3)... no.
# Actually, 4-vertex tournaments have 3 types by c_3:
# c_3=0: transitive (score (0,1,2,3)), 1 DT path
# c_3=1: exactly one 3-cycle, multiple score types
# c_3=4: all 4 triples are 3-cycles (score (1,1,2,2) type? No...)
#
# Wait, c_3 counts 3-cycles among the C(4,3)=4 triples.
# For n=4: c_3 can be 0, 1, or 4.
# c_3=4 is impossible: if ALL 4 triples are cycles, the score must be (1,1,2,2)
# But (1,1,2,2) has c_3 = 1, not 4. Let me check.

# At n=4, score sequence determines the tournament type:
# (0,1,2,3): transitive, c_3 = 0
# (0,1,2,3) is the only way to have 0+1+2+3=6=C(4,2) total score
# (1,1,1,3): c_3 = ?
# (0,2,2,2): c_3 = ?
# (1,1,2,2): c_3 = ?

n4 = 4
total4 = 1 << (n4*(n4-1)//2)
c3_by_score = defaultdict(set)
dt_by_c3_4 = Counter()

for bits in range(total4):
    A = build_adj(n4, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = 0
    for i in range(n4):
        for j in range(i+1, n4):
            for k in range(j+1, n4):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1
    c3_by_score[scores].add(c3)

    # DT paths
    a3 = [(a,b,c,d) for a in range(n4) for b in range(n4) for c in range(n4) for d in range(n4)
           if len({a,b,c,d})==4 and A[a][b] and A[b][c] and A[c][d] and A[a][c] and A[b][d]]
    dt_by_c3_4[(scores, c3)] += len(a3)

print("n=4 tournaments: c_3 by score sequence")
for scores in sorted(c3_by_score.keys()):
    print(f"  {scores}: c_3 in {c3_by_score[scores]}")

print("\nn=4: DT paths by (score, c_3)")
# Count tournaments per type
n_by_type = Counter()
for bits in range(total4):
    A = build_adj(n4, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = 0
    for i in range(n4):
        for j in range(i+1, n4):
            for k in range(j+1, n4):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1
    n_by_type[(scores, c3)] += 1

for key in sorted(dt_by_c3_4.keys()):
    avg_dt = dt_by_c3_4[key] / n_by_type[key]
    print(f"  {key}: avg DT = {avg_dt:.1f} ({n_by_type[key]} tournaments)")


# ============================================================
# KEY INSIGHT: DT paths per 4-vertex set
# ============================================================
print(f"\n{'='*70}")
print("KEY: DT paths per 4-vertex sub-tournament")
print("=" * 70)

# For each 4-vertex sub-tournament type, count DT paths:
# Transitive (c_3=0): 1 DT path
# One 3-cycle (c_3=1): 0 DT paths? Let me check.

print("\nChecking DT count per 4-vertex sub-tournament:")
# Build all n=4 tournaments and count DT paths
dt_per_c3 = Counter()
n_per_c3 = Counter()

for bits in range(total4):
    A = build_adj(n4, bits)
    c3 = 0
    for i in range(n4):
        for j in range(i+1, n4):
            for k in range(j+1, n4):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1

    n_dt = 0
    for a in range(n4):
        for b in range(n4):
            if b == a or not A[a][b]:
                continue
            for c in range(n4):
                if c == a or c == b or not A[b][c]:
                    continue
                for d in range(n4):
                    if d == a or d == b or d == c or not A[c][d]:
                        continue
                    if A[a][c] and A[b][d]:
                        n_dt += 1

    dt_per_c3[c3] += n_dt
    n_per_c3[c3] += 1

for c3 in sorted(dt_per_c3.keys()):
    avg = dt_per_c3[c3] / n_per_c3[c3]
    print(f"  c_3={c3}: avg DT = {avg:.2f} ({n_per_c3[c3]} tournaments)")

# So the number of DT paths relates to the number of transitive 4-sets
# DT = sum over transitive 4-sets of 1 = (number of transitive 4-tuples)
# A transitive 4-tuple (a,b,c,d) has a->b->c->d, a->c, b->d
# This is precisely: a 4-chain in a transitive sub-tournament of size >= 4

# For the TOTAL number of DT paths in a tournament T on n vertices:
# |DT| = sum over 4-vertex subsets S of (number of DT paths in T|_S)
#       = (number of transitive 4-vertex subsets) * 1 + contributions from c_3=1 sets
# But we need to check if c_3=1 sets contribute DT paths.

print("\nDoes c_3 = 1 contribute DT paths?")
for bits in range(total4):
    A = build_adj(n4, bits)
    c3 = 0
    for i in range(n4):
        for j in range(i+1, n4):
            for k in range(j+1, n4):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1
    if c3 != 1:
        continue

    n_dt = 0
    for a in range(n4):
        for b in range(n4):
            if b == a or not A[a][b]:
                continue
            for c in range(n4):
                if c == a or c == b or not A[b][c]:
                    continue
                for d in range(n4):
                    if d == a or d == b or d == c or not A[c][d]:
                        continue
                    if A[a][c] and A[b][d]:
                        n_dt += 1
    if n_dt > 0:
        scores = [sum(row) for row in A]
        print(f"  YES: c_3=1, DT={n_dt}, scores={scores}")
        # Show the DT path
        for a in range(n4):
            for b in range(n4):
                if b == a or not A[a][b]:
                    continue
                for c in range(n4):
                    if c == a or c == b or not A[b][c]:
                        continue
                    for d in range(n4):
                        if d == a or d == b or d == c or not A[c][d]:
                            continue
                        if A[a][c] and A[b][d]:
                            print(f"    DT path: ({a},{b},{c},{d})")
        break


# ============================================================
# FORMULA: |DT| = #{transitive 4-tuples} = 4*t_4 where t_4 = # transitive 4-sets
# Wait, per 4-set: transitive gives 1 DT, c_3=1 might give some too
# ============================================================
print(f"\n{'='*70}")
print("FORMULA: |DT| as function of combinatorial invariants")
print("=" * 70)

# At n=5, let's verify
n = 5
total = 1 << (n*(n-1)//2)

from itertools import combinations

for bits in range(min(total, 200)):
    A = build_adj(n, bits)
    a3 = enumerate_allowed_paths(A, n, 3)
    n_dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    # Count transitive 4-sets
    t4 = 0
    for quad in combinations(range(n), 4):
        c3_sub = 0
        for triple in combinations(quad, 3):
            i, j, k = triple
            si = A[i][j] + A[i][k]
            sj = A[j][i] + A[j][k]
            sk = A[k][i] + A[k][j]
            if max(si, sj, sk) < 2:
                c3_sub += 1
        if c3_sub == 0:
            t4 += 1

    # Check formula: |DT| = t4 + (contribution from non-transitive 4-sets)
    # We need to compute non-transitive contribution
    dt_from_nontrans = 0
    for quad in combinations(range(n), 4):
        c3_sub = 0
        for triple in combinations(quad, 3):
            i, j, k = triple
            si = A[i][j] + A[i][k]
            sj = A[j][i] + A[j][k]
            sk = A[k][i] + A[k][j]
            if max(si, sj, sk) < 2:
                c3_sub += 1
        if c3_sub > 0:
            # Count DT paths within this 4-set
            q = list(quad)
            for a in q:
                for b in q:
                    if b == a or not A[a][b]:
                        continue
                    for c in q:
                        if c == a or c == b or not A[b][c]:
                            continue
                        for d in q:
                            if d == a or d == b or d == c or not A[c][d]:
                                continue
                            if A[a][c] and A[b][d]:
                                dt_from_nontrans += 1

    if bits < 5:
        print(f"  bits={bits}: |DT|={n_dt}, t4={t4}, non-trans contribution={dt_from_nontrans}")

# So we need to understand the structure of c_3=1 four-sets that still have DT paths
# This would be: a 4-set with exactly one 3-cycle but still having a transitive 4-chain
# This means the 3-cycle doesn't involve the "top" or "bottom" vertex of the 4-chain


print("\nDone.")
