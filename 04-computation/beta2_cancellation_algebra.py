#!/usr/bin/env python3
"""
beta2_cancellation_algebra.py — The cancellation chain mechanism for β₂=0

KEY INSIGHT FROM PREVIOUS ANALYSIS:
- Two 3-paths sharing the same bad face → their difference is in Ω₃
- Example: (a,b₁,c,d) and (a,b₂,c,d) with bad face (a,c,d) cancel
- This is the SAME mechanism as at dimension 2:
  (a,b₁,c) - (a,b₂,c) with bad face (a,c) cancel to give Ω₂ element

ALGEBRAIC STRUCTURE:
For a 3-path (a,b,c,d) ∈ A₃, the boundary ∂ has 4 faces:
  ∂(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

A face (x,y,z) is "bad" if it's NOT in A₂ (i.e., not all edges present).
In a tournament, (x,y,z) is in A₂ iff x→y AND y→z (both edges directed right).

For (a,b,c,d) ∈ A₃ (a→b→c→d), the faces:
  (b,c,d): b→c→d ✓ always in A₂
  (a,b,c): a→b→c ✓ always in A₂
  (a,c,d): a→c? depends. c→d ✓. So bad iff c→a.
  (a,b,d): a→b ✓. b→d? depends. So bad iff d→b.

So a 3-path has 0, 1, or 2 bad faces, depending on edges a-c and b-d.
- 0 bad faces: a→c AND b→d (= DT path) → path itself is in Ω₃
- 1 bad face: exactly one of c→a or d→b
- 2 bad faces: c→a AND d→b

For pairs of 3-paths sharing a bad face, their difference is in Ω₃.
The question: does this cancellation mechanism ALWAYS generate enough
Ω₃ elements to fill all Z₂?

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ======================================================================
print("="*70)
print("CANCELLATION ALGEBRA FOR β₂=0")
print("="*70)

n = 5
print(f"\nn = {n}")

# For each tournament, classify 3-paths by their bad face count
bad_face_stats = Counter()
cancellation_groups = defaultdict(list)  # bad_face -> list of 3-paths sharing it

total = 0
for A in all_tournaments(n):
    total += 1
    if total > 100: break  # sample

    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    bad_face_groups = defaultdict(list)  # bad_face tuple -> list of 3-paths

    for p in a3:
        a_, b_, c_, d_ = p
        bad_faces = []

        # Face (a,c,d): bad iff c→a
        if not A[a_][c_]:  # c→a
            bad_faces.append(('face_acd', (a_,c_,d_)))

        # Face (a,b,d): bad iff d→b
        if not A[b_][d_]:  # d→b
            bad_faces.append(('face_abd', (a_,b_,d_)))

        bad_face_stats[len(bad_faces)] += 1

        for bf_type, bf in bad_faces:
            bad_face_groups[bf].append(p)

    # How many bad faces have multiplicity > 1?
    for bf, paths in bad_face_groups.items():
        if len(paths) > 1:
            cancellation_groups[len(paths)].append((bf, paths))

print(f"\nBad face count distribution (sample of {total} tournaments):")
for k, v in sorted(bad_face_stats.items()):
    print(f"  {k} bad faces: {v} 3-paths")

print(f"\nCancellation group sizes:")
for size, groups in sorted(cancellation_groups.items()):
    print(f"  Multiplicity {size}: {len(groups)} groups")

# ======================================================================
# DETAILED: Show cancellation chains and their boundaries
print(f"\n\n{'='*70}")
print("DETAILED CANCELLATION ANALYSIS")
print("="*70)

# Take a specific tournament and trace the full mechanism
n = 5
A = None
for T in all_tournaments(n):
    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (T[i][j] and T[j][k] and T[k][i]) or (T[j][i] and T[i][k] and T[k][j]))
    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    if t3 == 4 and scores == (1,2,2,2,3):
        A = T
        break

if A is None:
    # Just use first tournament with cancellation chains
    for T in all_tournaments(n):
        A = T
        break

print(f"\nTournament: t3={sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n) if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))}")
print(f"Adjacency:")
for i in range(n):
    print(f"  {A[i]}")

a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

print(f"\n|A₁|={len(a1)}, |A₂|={len(a2)}, |A₃|={len(a3)}")

# Classify all 3-paths
print(f"\nAll 3-paths classified by bad faces:")
dt_paths = []
one_bad = defaultdict(list)
two_bad = []

for p in a3:
    a_, b_, c_, d_ = p
    bads = []
    if not A[a_][c_]:
        bads.append(('acd', a_, c_, d_))
    if not A[b_][d_]:
        bads.append(('abd', a_, b_, d_))

    if len(bads) == 0:
        dt_paths.append(p)
        print(f"  {p}: DT (0 bad faces)")
    elif len(bads) == 1:
        key = (bads[0][0], bads[0][1:])
        one_bad[key].append(p)
        print(f"  {p}: 1 bad face {bads[0]}")
    else:
        two_bad.append(p)
        print(f"  {p}: 2 bad faces {bads}")

print(f"\nDT: {len(dt_paths)}, 1-bad: {sum(len(v) for v in one_bad.values())}, 2-bad: {len(two_bad)}")

# Show cancellation pairs
print(f"\nCancellation pairs (sharing same bad face):")
for key, paths in one_bad.items():
    if len(paths) >= 2:
        print(f"  Bad face {key}: {len(paths)} paths: {paths}")
        # The differences p₁ - p₂ are in Ω₃
        # Show what ∂(p₁ - p₂) looks like
        p1, p2 = paths[0], paths[1]
        print(f"    ∂({p1}) - ∂({p2}):")

        # Compute boundaries
        for p, sign_prefix in [(p1, '+'), (p2, '-')]:
            faces = [(p[1:], +1), (p[:1]+p[2:], -1), (p[:2]+p[3:], +1), (p[:3], -1)]
            for face, sgn in faces:
                in_a2 = face in a2
                total_sign = '+' if sgn > 0 else '-'
                if sign_prefix == '-':
                    total_sign = '-' if total_sign == '+' else '+'
                status = "✓" if in_a2 else "BAD"
                print(f"      {sign_prefix} {total_sign}{face} [{status}]")

# ======================================================================
# BOUNDARY ANALYSIS: What does ∂(cancellation chain) look like in Z₂?
print(f"\n\n{'='*70}")
print("BOUNDARY OF CANCELLATION CHAINS IN Z₂")
print("="*70)

# For each cancellation pair (p₁, p₂) sharing bad face f:
# ∂(p₁ - p₂) consists of the NON-BAD faces of p₁ minus those of p₂
# The bad face terms cancel!

# At the level of A₂ coordinates:
# ∂(a,b₁,c,d) = +(b₁,c,d) -(a,c,d) +(a,b₁,d) -(a,b₁,c)
# ∂(a,b₂,c,d) = +(b₂,c,d) -(a,c,d) +(a,b₂,d) -(a,b₂,c)
# Difference: +(b₁,c,d)-(b₂,c,d) +(a,b₁,d)-(a,b₂,d) -(a,b₁,c)+(a,b₂,c)
# The -(a,c,d) terms CANCEL!

# This is exactly analogous to the Ω₂ cancellation:
# (a,b₁,c) - (a,b₂,c) has boundary +(b₁,c)-(b₂,c) +(a,b₁)-(a,b₂)
# The -(a,c) terms CANCEL!

print(f"\nPattern: Ω₃ cancellation chains → boundary in Ω₂")
print(f"  just as Ω₂ cancellation chains → boundary in Ω₁")
print(f"\nThis is a TELESCOPING structure across dimensions:")
print(f"  dim p cancellation = pairs of p-paths sharing same bad (p-1)-face")
print(f"  dim p+1 cancellation fills dim p cancellation cycles")

# ======================================================================
# KEY TEST: Does the cancellation mechanism explain ALL of Z₂?
print(f"\n\n{'='*70}")
print("DOES CANCELLATION + DT EXPLAIN ALL OF Z₂?")
print("="*70)

n = 5
dt_fills = 0
cancel_needed = 0
total_tested = 0

for A in all_tournaments(n):
    total_tested += 1

    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a2_set = set(a2)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if d_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2, a1)
    S = np.linalg.svd(bd2 @ om2, compute_uv=False)
    z2 = d_om2 - sum(s > 1e-8 for s in S)

    # DT paths
    dt = [p for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]]]

    # Cancellation pairs: paths sharing a bad face
    bad_face_groups = defaultdict(list)
    for p in a3:
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad_face_groups[(a_,c_,d_)].append(p)
        if not A[b_][d_]:
            bad_face_groups[(a_,b_,d_)].append(p)

    cancel_pairs = []
    for bf, paths in bad_face_groups.items():
        if len(paths) >= 2:
            # Add all pairwise differences
            for i in range(len(paths)):
                for j in range(i+1, len(paths)):
                    cancel_pairs.append((paths[i], paths[j]))

    # Build combined Ω₃ space: DT singles + cancellation differences
    # Each DT path is a basis vector
    # Each cancellation difference p₁ - p₂ is a vector

    vectors = []
    for p in dt:
        v = np.zeros(len(a3))
        v[a3.index(p)] = 1
        vectors.append(v)

    for p1, p2 in cancel_pairs:
        v = np.zeros(len(a3))
        v[a3.index(p1)] = 1
        v[a3.index(p2)] = -1
        vectors.append(v)

    if vectors:
        V = np.column_stack(vectors)
        # Compute im(∂₃ on this space)
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_V = bd3 @ V

        # Project onto Ω₂
        if om2.ndim == 2 and om2.shape[1] > 0:
            coords, _, _, _ = np.linalg.lstsq(om2, bd3_V, rcond=None)
            rk = np.linalg.matrix_rank(coords, tol=1e-8)
        else:
            rk = 0
    else:
        rk = 0

    if rk >= z2:
        dt_fills += 1
    else:
        cancel_needed += 1
        if cancel_needed <= 3:
            t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                     if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))
            print(f"  DT+cancel insufficient: rk={rk}, Z₂={z2}, t3={t3}")
            print(f"    DT: {len(dt)}, cancel pairs: {len(cancel_pairs)}")

print(f"\n  DT + single cancellation pairs fill Z₂: {dt_fills}/{total_tested}")
print(f"  Need higher-order cancellations: {cancel_needed}/{total_tested}")

# ======================================================================
# What about 2-bad-face paths? Can they contribute?
print(f"\n\n{'='*70}")
print("2-BAD-FACE PATHS AND HIGHER-ORDER CANCELLATIONS")
print("="*70)

# A 3-path with 2 bad faces: c→a AND d→b
# Can THREE such paths cancel both bad faces?
# e.g., (a,b₁,c,d₁) - (a,b₂,c,d₁) - (a,b₁,c,d₂) + (a,b₂,c,d₂)
# This cancels both bad face types if b₁,b₂ share the c→a bad face
# and d₁,d₂ share the b→d bad face

# But this is already captured by taking differences of cancellation pairs
# at the 1-bad-face level (which we did above)

# Let me check: do ALL Ω₃ elements arise from:
# 1. DT paths (0 bad faces)
# 2. Pairwise differences of 1-bad-face paths sharing same bad face
# 3. Linear combinations of 2-bad-face paths where BOTH bad faces cancel

n = 5
om3_decomposition = Counter()

for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a2_set = set(a2)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    for j in range(d_om3):
        col = om3[:, j]
        support = [i for i in range(len(col)) if abs(col[i]) > 1e-10]
        paths = [a3[i] for i in support]

        # Classify bad faces of each path in support
        max_bad = 0
        for p in paths:
            a_, b_, c_, d_ = p
            bads = 0
            if not A[a_][c_]: bads += 1
            if not A[b_][d_]: bads += 1
            max_bad = max(max_bad, bads)

        if len(support) == 1 and max_bad == 0:
            om3_decomposition['DT_single'] += 1
        elif max_bad <= 1:
            om3_decomposition['1bad_cancel'] += 1
        else:
            om3_decomposition['2bad_cancel'] += 1

print(f"\nΩ₃ element classification at n={n}:")
for key, cnt in sorted(om3_decomposition.items()):
    print(f"  {key}: {cnt}")

# ======================================================================
# FORMULA: dim(Ω₃) in terms of path structure
print(f"\n\n{'='*70}")
print("Ω₃ DIMENSION FORMULA")
print("="*70)

# dim(Ω₃) = |DT| + (cancellation chains)
# Cancellation chains come from bad face multiplicities

for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a2_set = set(a2)

    dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    # Count bad face multiplicities
    bad_face_mult = defaultdict(int)
    for p in a3:
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad_face_mult[(a_,c_,d_)] += 1
        if not A[b_][d_]:
            bad_face_mult[(a_,b_,d_)] += 1

    # cancellation dimension = Σ (mult - 1) for mult > 1
    # But this might overcount if paths have 2 bad faces
    cancel_dim = sum(max(m-1, 0) for m in bad_face_mult.values())

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))

    predicted = dt + cancel_dim
    match = "✓" if predicted == d_om3 else f"✗ (pred={predicted})"

    break  # just first

print(f"  DT={dt}, cancel_dim={cancel_dim}, pred={predicted}, actual Ω₃={d_om3} {match}")

# Do this for all tournaments
print(f"\n  Testing formula dim(Ω₃) = |DT| + Σ(mult-1) at n={n}:")
matches = 0
mismatches = 0
for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])
    bad_face_mult = defaultdict(int)
    for p in a3:
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad_face_mult[(a_,c_,d_)] += 1
        if not A[b_][d_]:
            bad_face_mult[(a_,b_,d_)] += 1
    cancel_dim = sum(max(m-1, 0) for m in bad_face_mult.values())

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dt + cancel_dim == d_om3:
        matches += 1
    else:
        mismatches += 1
        if mismatches <= 3:
            t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                     if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))
            print(f"    MISMATCH: DT={dt}, cancel={cancel_dim}, pred={dt+cancel_dim}, actual={d_om3}, t3={t3}")

print(f"\n  Matches: {matches}/{matches+mismatches}")
if mismatches == 0:
    print(f"  ✓ UNIVERSAL: dim(Ω₃) = |DT| + Σ(mult-1) for ALL n={n} tournaments!")
else:
    print(f"  Formula fails for {mismatches} tournaments")

# Now test at n=6
n = 6
print(f"\n  Testing formula at n={n} (sample):")
matches_6 = 0
mismatches_6 = 0
total_6 = 0

pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)
import random
random.seed(42)

for _ in range(2000):
    total_6 += 1
    mask = random.randint(0, (1 << m) - 1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])
    bad_face_mult = defaultdict(int)
    for p in a3:
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad_face_mult[(a_,c_,d_)] += 1
        if not A[b_][d_]:
            bad_face_mult[(a_,b_,d_)] += 1
    cancel_dim = sum(max(m_val-1, 0) for m_val in bad_face_mult.values())

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dt + cancel_dim == d_om3:
        matches_6 += 1
    else:
        mismatches_6 += 1
        if mismatches_6 <= 3:
            print(f"    MISMATCH: DT={dt}, cancel={cancel_dim}, pred={dt+cancel_dim}, actual={d_om3}")

    if total_6 % 500 == 0:
        print(f"    ... {total_6}", flush=True)

print(f"\n  n={n}: Matches: {matches_6}/{total_6}")

print(f"\n{'='*70}")
print("DONE")
print("="*70)
