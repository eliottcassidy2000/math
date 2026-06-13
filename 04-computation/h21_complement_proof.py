#!/usr/bin/env python3
"""
h21_complement_proof.py — Prove remaining gaps for H=21 impossibility.

KEY DISCOVERY: At n=6, if {a,b,c} is a directed 3-cycle, is its
complement triple {d,e,f} always transitive? If so, alpha_2=0 always
at n=6, which gives a structural explanation for the T=10 gap.

For n≥7, disjoint 3-cycles exist on non-complementary triples.
We need different arguments there.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random

sys.stdout.reconfigure(line_buffering=True)

def get_directed_cycles(A, n):
    groups = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    groups[frozenset(verts)] += 1
    return groups

def compute_alpha(groups):
    vs_list = list(groups.items())
    n_vs = len(vs_list)
    alpha1 = sum(d for _, d in vs_list)
    alpha2 = 0
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if not (vs_list[i][0] & vs_list[j][0]):
                alpha2 += vs_list[i][1] * vs_list[j][1]
    return alpha1, alpha2

print("=" * 70)
print("H=21 COMPLEMENT PROOF")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: At n=6, test if cyclic triple → complement transitive
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: Complementary triple constraint at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

both_cyclic_count = 0
one_cyclic_count = 0
neither_count = 0

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    all_verts = set(range(6))
    for verts in combinations(range(6), 3):
        if verts[0] != 0:  # avoid double-counting complementary pairs
            continue
        vs = frozenset(verts)
        comp = all_verts - vs
        comp_list = sorted(comp)
        v = sorted(verts)

        # Check if vs is a 3-cycle (all edges go one way)
        score = A[v[0]][v[1]] + A[v[1]][v[2]] + A[v[2]][v[0]]
        vs_cyclic = (score in (0, 3))

        # Check if complement is a 3-cycle
        c = comp_list
        score_c = A[c[0]][c[1]] + A[c[1]][c[2]] + A[c[2]][c[0]]
        comp_cyclic = (score_c in (0, 3))

        if vs_cyclic and comp_cyclic:
            both_cyclic_count += 1
        elif vs_cyclic or comp_cyclic:
            one_cyclic_count += 1
        else:
            neither_count += 1

total_pairs = 2**ne * 10  # 10 complementary pairs per tournament
print(f"  Total complementary pair instances: {both_cyclic_count + one_cyclic_count + neither_count}")
print(f"  Both cyclic: {both_cyclic_count}")
print(f"  Exactly one cyclic: {one_cyclic_count}")
print(f"  Neither cyclic: {neither_count}")

if both_cyclic_count > 0:
    print(f"\n  COMPLEMENT LEMMA FAILS at n=6!")
    print(f"  But this means alpha_2 > 0 IS possible... checking which alpha_1")

    # Find examples
    examples = []
    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        all_verts = set(range(6))
        for verts in combinations(range(6), 3):
            vs = frozenset(verts)
            comp = all_verts - vs
            v = sorted(verts)
            c = sorted(comp)

            score = A[v[0]][v[1]] + A[v[1]][v[2]] + A[v[2]][v[0]]
            score_c = A[c[0]][c[1]] + A[c[1]][c[2]] + A[c[2]][c[0]]

            if score in (0,3) and score_c in (0,3):
                groups = get_directed_cycles(A, n)
                a1, a2 = compute_alpha(groups)
                if len(examples) < 5:
                    examples.append((bits, v, c, a1, a2))

    print(f"  Examples of both-cyclic complementary pairs:")
    for bits, v, c, a1, a2 in examples:
        print(f"    bits={bits}: {v} and {c} both cyclic, alpha_1={a1}, alpha_2={a2}")
else:
    print(f"\n  PROVED: No complementary triple pair is both cyclic at n=6!")
    print(f"  This means alpha_2=0 for all 3-cycle-only tournaments at n=6.")

# ═══════════════════════════════════════════════════════════════════
# Part 2: Full (alpha_1, alpha_2) landscape at n=6
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Full (alpha_1, alpha_2) landscape at n=6 ---")

a_landscape = Counter()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)
    a_landscape[(a1, a2)] += 1

print(f"  {'a1':>4} {'a2':>4} {'count':>6}  T=a1+2a2")
for (a1, a2) in sorted(a_landscape.keys()):
    T = a1 + 2*a2
    marker = " ← T=10!" if T == 10 else ""
    print(f"  {a1:4d} {a2:4d} {a_landscape[(a1,a2)]:6d}  T={T:3d}{marker}")

# Show which T values are achievable
T_vals = Counter()
for (a1, a2), count in a_landscape.items():
    T_vals[a1 + 2*a2] += count

print(f"\n  T-spectrum gaps: {[t for t in range(max(T_vals)+1) if t not in T_vals]}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: At n=7, sample the landscape near T=10
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: n=7 sampling near T=10 ---")

random.seed(42)
n = 7
edges7 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne7 = len(edges7)

near_t10 = Counter()  # (a1, a2) where T ∈ [8,12]
SAMPLES = 300000

for _ in range(SAMPLES):
    bits = random.randint(0, 2**ne7 - 1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges7):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)
    T = a1 + 2*a2
    if 6 <= T <= 14:
        near_t10[(a1, a2)] += 1

print(f"  n=7, {SAMPLES} samples, (alpha_1, alpha_2) with T ∈ [6,14]:")
for (a1, a2) in sorted(near_t10.keys()):
    T = a1 + 2*a2
    marker = " ← T=10!" if T == 10 else ""
    print(f"    ({a1:3d}, {a2:3d}): T={T:3d}, count={near_t10[(a1,a2)]:5d}{marker}")

t10_total = sum(v for (a1,a2), v in near_t10.items() if a1+2*a2==10)
print(f"\n  T=10 total: {t10_total}")
if t10_total == 0:
    print("  CONFIRMED: T=10 absent at n=7 (300k sample)")

# ═══════════════════════════════════════════════════════════════════
# Part 4: The achievable alpha_2 for each alpha_1 at n=7
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: alpha_2 spectrum by alpha_1 at n=7 (sample) ---")

full_spectrum_7 = defaultdict(Counter)
for _ in range(SAMPLES):
    bits = random.randint(0, 2**ne7 - 1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges7):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)
    full_spectrum_7[a1][a2] += 1

for a1 in sorted(full_spectrum_7.keys()):
    if a1 > 20:
        break
    a2_vals = full_spectrum_7[a1]
    a2_list = sorted(a2_vals.keys())
    total = sum(a2_vals.values())
    required_for_t10 = (10 - a1) / 2
    if a1 <= 10 and required_for_t10 == int(required_for_t10) and required_for_t10 >= 0:
        req_str = f"  need a2={int(required_for_t10)} for T=10"
        in_gap = int(required_for_t10) not in a2_vals
        req_str += " → IN GAP" if in_gap else " → ACHIEVABLE"
    else:
        req_str = ""
    print(f"  a1={a1:3d}: a2 ∈ {a2_list} (n={total}){req_str}")

print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
