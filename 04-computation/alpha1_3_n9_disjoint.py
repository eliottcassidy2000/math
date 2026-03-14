#!/usr/bin/env python3
"""
alpha1_3_n9_disjoint.py — opus-2026-03-14-S71g

CRITICAL CHECK: At n=9, can α₁=3 with α₂=0?

Construction: 3 disjoint directed 3-cycles covering all 9 vertices,
with transitive inter-group structure (no inter-group 3-cycles).

Groups: {0,1,2}, {3,4,5}, {6,7,8}
Within each: directed 3-cycle 0→1→2→0, 3→4→5→3, 6→7→8→6
Between: group 0 beats group 1, group 1 beats group 2, group 0 beats group 2
(transitively ordered)

Prediction: H = I(Ω, 2) = (1+2)^3 = 27 (3 isolated simplex bricks)
"""

from itertools import permutations, combinations

def count_hp(A, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def count_directed_odd_cycles(A, n, max_length=None):
    """Count all directed odd cycles (canonicalized by min vertex start)."""
    if max_length is None:
        max_length = n
    cycles_by_length = {}
    for k in range(3, min(n, max_length) + 1, 2):
        count = 0
        for verts in combinations(range(n), k):
            for p in permutations(verts[1:]):
                order = [verts[0]] + list(p)
                ok = True
                for idx in range(k):
                    if A[order[idx]][order[(idx+1) % k]] != 1:
                        ok = False
                        break
                if ok:
                    count += 1
        cycles_by_length[k] = count
    return cycles_by_length

# ============================================================
# Construct the n=9 tournament
# ============================================================

n = 9
A = [[0]*n for _ in range(n)]

# Groups: G0={0,1,2}, G1={3,4,5}, G2={6,7,8}
groups = [[0,1,2], [3,4,5], [6,7,8]]

# Within-group 3-cycles
for g in groups:
    a, b, c = g
    A[a][b] = 1  # a→b
    A[b][c] = 1  # b→c
    A[c][a] = 1  # c→a

# Between-group: G0 beats G1, G1 beats G2, G0 beats G2 (transitive)
for a in groups[0]:
    for b in groups[1]:
        A[a][b] = 1
    for c in groups[2]:
        A[a][c] = 1
for b in groups[1]:
    for c in groups[2]:
        A[b][c] = 1

# Verify: every pair has exactly one arc
print("=" * 70)
print("α₁=3 DISJOINT CONSTRUCTION at n=9")
print("=" * 70)

for i in range(n):
    for j in range(n):
        if i == j:
            assert A[i][j] == 0
        else:
            assert A[i][j] + A[j][i] == 1, f"Error at ({i},{j})"

# Score sequence
degs = [sum(A[i]) for i in range(n)]
print(f"\nOut-degrees: {degs}")
print(f"Score sequence: {tuple(sorted(degs))}")

# Count cycles by length
print(f"\nCounting directed odd cycles...")
cycles = count_directed_odd_cycles(A, n, max_length=9)
for k, count in sorted(cycles.items()):
    print(f"  d_{k} = {count}")

alpha1 = sum(cycles.values())
print(f"  α₁ = {alpha1}")

# Count HP
print(f"\nComputing H(T) via Held-Karp DP...")
H = count_hp(A, n)
print(f"  H = {H}")

# Check OCF prediction
print(f"\nOCF prediction:")
if alpha1 == 3 and cycles.get(3, 0) == 3:
    # Check if cycles are disjoint (they are by construction)
    print(f"  3 disjoint 3-cycles → Ω = 3 isolated vertices")
    print(f"  I(Ω, 2) = (1+2)^3 = 27")
    print(f"  H matches OCF: {H == 27}")
else:
    print(f"  Unexpected α₁={alpha1}, not exactly 3 disjoint 3-cycles")
    print(f"  Cannot apply simple formula")

# ============================================================
# Try alternate inter-group orderings
# ============================================================

print(f"\n{'='*70}")
print("ALTERNATE INTER-GROUP ORDERINGS")
print(f"{'='*70}")

# What if inter-group arcs form a 3-cycle instead of transitive?
# G0 beats G1, G1 beats G2, G2 beats G0

n = 9
A2 = [[0]*n for _ in range(n)]

for g in groups:
    a, b, c = g
    A2[a][b] = 1
    A2[b][c] = 1
    A2[c][a] = 1

# Cyclic between groups
for a in groups[0]:
    for b in groups[1]:
        A2[a][b] = 1  # G0 → G1
for b in groups[1]:
    for c in groups[2]:
        A2[b][c] = 1  # G1 → G2
for c in groups[2]:
    for a in groups[0]:
        A2[c][a] = 1  # G2 → G0

for i in range(n):
    for j in range(n):
        if i != j:
            assert A2[i][j] + A2[j][i] == 1, f"Error at ({i},{j})"

degs2 = [sum(A2[i]) for i in range(n)]
print(f"\nOut-degrees (cyclic inter-group): {degs2}")
print(f"Score sequence: {tuple(sorted(degs2))}")

print(f"\nCounting directed odd cycles...")
cycles2 = count_directed_odd_cycles(A2, n, max_length=9)
for k, count in sorted(cycles2.items()):
    print(f"  d_{k} = {count}")

alpha1_2 = sum(cycles2.values())
print(f"  α₁ = {alpha1_2}")

H2 = count_hp(A2, n)
print(f"  H = {H2}")
print(f"  This is the regular tournament (all degrees = 4)")

# ============================================================
# Exhaustive search at n=9: sample for α₁=3
# ============================================================

print(f"\n{'='*70}")
print("SAMPLE SEARCH for α₁=3 at n=9 (random tournaments)")
print(f"{'='*70}")

import random
random.seed(42)
n = 9
total_edges = n * (n - 1) // 2

found_alpha1_3 = []
sample_size = 10000

for trial in range(sample_size):
    bits = random.randint(0, 2**total_edges - 1)
    A_r = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A_r[i][j] = 1
            else:
                A_r[j][i] = 1
            idx += 1

    # Quick 3-cycle count
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A_r[i][j] and A_r[j][k] and A_r[k][i]:
            t3 += 1
        elif A_r[i][k] and A_r[k][j] and A_r[j][i]:
            t3 += 1

    if t3 <= 5:  # Only check full cycle counts if t3 is low
        cycles_r = count_directed_odd_cycles(A_r, n, max_length=5)  # skip 7,9 for speed
        alpha1_approx = sum(cycles_r.values())
        if alpha1_approx <= 5:
            # Now check 7-cycles too
            cycles_r = count_directed_odd_cycles(A_r, n, max_length=9)
            alpha1_full = sum(cycles_r.values())
            if alpha1_full == 3:
                found_alpha1_3.append((bits, t3, cycles_r))
                if len(found_alpha1_3) <= 5:
                    H_r = count_hp(A_r, n)
                    print(f"  Found α₁=3: t₃={t3}, d₅={cycles_r.get(5,0)}, d₇={cycles_r.get(7,0)}, d₉={cycles_r.get(9,0)}, H={H_r}")

    if trial % 2000 == 0:
        print(f"  ... {trial}/{sample_size} processed, {len(found_alpha1_3)} found")

print(f"\n  Total with α₁=3: {len(found_alpha1_3)}/{sample_size}")
if found_alpha1_3:
    # Check if all have disjoint cycles
    print(f"  All found cases:")
    for bits, t3, cyc in found_alpha1_3[:10]:
        print(f"    t₃={t3}, d₅={cyc.get(5,0)}, d₇={cyc.get(7,0)}, d₉={cyc.get(9,0)}")

# ============================================================
# The key question: at n=9, does α₁=3 → α₂≥2 or can α₂=0?
# ============================================================

print(f"\n{'='*70}")
print("EXPLICIT CHECK: Our n=9 construction")
print(f"{'='*70}")

# Rebuild our construction
n = 9
A = [[0]*n for _ in range(n)]
for g in groups:
    a, b, c = g
    A[a][b] = 1
    A[b][c] = 1
    A[c][a] = 1
for a in groups[0]:
    for b in groups[1]:
        A[a][b] = 1
    for c in groups[2]:
        A[a][c] = 1
for b in groups[1]:
    for c in groups[2]:
        A[b][c] = 1

# List all directed odd cycles
all_cycles = []
for k in range(3, n+1, 2):
    for verts in combinations(range(n), k):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            ok = True
            for idx in range(k):
                if A[order[idx]][order[(idx+1) % k]] != 1:
                    ok = False
                    break
            if ok:
                all_cycles.append(tuple(order))

print(f"  Total directed odd cycles: {len(all_cycles)}")
for c in all_cycles:
    vset = set(c)
    print(f"    {c} (vertices: {sorted(vset)}, length: {len(c)})")

# Build conflict graph
nc = len(all_cycles)
vsets = [frozenset(c) for c in all_cycles]
adj = [[False]*nc for _ in range(nc)]
for i in range(nc):
    for j in range(i+1, nc):
        if vsets[i] & vsets[j]:
            adj[i][j] = adj[j][i] = True

# Count edges in conflict graph = α₂
alpha2 = sum(1 for i in range(nc) for j in range(i+1, nc) if adj[i][j])
print(f"\n  Conflict graph Ω:")
print(f"    α₁ = {nc} vertices")
print(f"    α₂ = {alpha2} edges")
if alpha2 == 0:
    print(f"    → Ω has NO edges (all cycles vertex-disjoint)")
    print(f"    → I(Ω, 2) = (1+2)^{nc} = {3**nc}")
    print(f"    → H should be {3**nc}")
else:
    print(f"    → Ω has edges (some cycles share vertices)")

H = count_hp(A, n)
print(f"\n  H(T) = {H} (via Held-Karp)")
print(f"  OCF prediction = {3**nc if alpha2 == 0 else 'complex'}")
print(f"  MATCH: {H == 3**nc if alpha2 == 0 else 'N/A'}")

print(f"\n{'='*70}")
print("CONCLUSION")
print(f"{'='*70}")
if nc == 3 and alpha2 == 0:
    print(f"""
At n=9: α₁=3 with α₂=0 IS ACHIEVABLE.

Construction: 3 vertex-disjoint 3-cycles covering all 9 vertices,
with transitive inter-group structure.

This gives:
  H = (1+2)^3 = 27 (pure simplex^3 decomposition)
  Ω = 3 isolated vertices (no conflicts)

The α₁=3 → α₂≥2 gap CLOSES at n=9.

For n ≤ 8: α₁=3 forces α₂≥2 (verified computationally).
For n = 9: α₁=3 can have α₂=0 (by explicit construction).

The critical threshold is n=9 = 3×3 (three disjoint 3-cycles need exactly 9 vertices).
""")
elif nc == 3 and alpha2 > 0:
    print(f"  UNEXPECTED: α₁=3 but some cycles share vertices!?")
else:
    print(f"  UNEXPECTED: α₁ = {nc} ≠ 3. Construction error or additional cycles exist.")
