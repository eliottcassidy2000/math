#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Investigate the α₁=3 gap — why can't a tournament have exactly 3 odd cycles?

Key finding from exploration: at n=3..6, the total number of directed odd cycles (α₁)
is never exactly 3. Values 0,1,2 occur; then it jumps to ≥4.

Hypothesis: α₁=3 is universally forbidden for all tournaments.
If true, this PROVES H=7 is impossible for all n.

Because H = 1 + 2α₁ + 4α₂ + ... = 7 requires:
- α₁=3, all higher=0 (only way since 4|remaining after 2α₁)
- But α₁=3 is forbidden → H=7 impossible. QED.

This script investigates:
1. Verify α₁=3 gap at n=7 (sampling)
2. Understand WHY: what constraint prevents exactly 3 odd cycles?
3. Check if similar gaps exist (α₁ values that are never achievable)
4. Extend to the H=21 forbidden value
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def count_3cycles(T):
    """Count 3-cycles using the score formula."""
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def count_5cycles(T):
    """Count directed 5-cycles by enumeration (exact)."""
    n = len(T)
    count = 0
    for verts in combinations(range(n), 5):
        # For 5 vertices, count directed Hamiltonian cycles
        # Fix first vertex, try all orderings of rest
        v0 = verts[0]
        for perm in permutations(verts[1:]):
            path = (v0,) + perm
            is_cycle = True
            for i in range(5):
                if T[path[i]][path[(i+1) % 5]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
                break  # one direction found; each undirected cycle has one valid direction
        # Actually each directed cycle is counted once per starting vertex
        # We fix v0 = smallest, so we count each cycle once
    return count

def count_all_odd_cycles(T, max_length=None):
    """Count all directed odd cycles up to given length."""
    n = len(T)
    if max_length is None:
        max_length = n

    counts = {}  # length -> count
    for length in range(3, max_length + 1, 2):
        if length > n:
            break
        count = 0
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                is_cycle = True
                for i in range(length):
                    if T[path[i]][path[(i+1) % length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    count += 1
                    break
        counts[length] = count
    return counts

def fast_alpha1(T):
    """Compute α₁ = total number of directed odd cycles."""
    cycles = count_all_odd_cycles(T)
    return sum(cycles.values()), cycles

# ─── Part 1: Verify α₁=3 gap at n=5,6 exhaustively ───
print("="*60)
print("PART 1: EXHAUSTIVE CHECK AT n=5")
print("="*60)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

alpha1_dist_5 = Counter()
alpha1_details_5 = defaultdict(list)
for T in all_tournaments(5):
    a1, cycles = fast_alpha1(T)
    alpha1_dist_5[a1] += 1
    if a1 == 3:
        alpha1_details_5[a1].append(cycles)
    elif a1 <= 4:
        key = (cycles.get(3,0), cycles.get(5,0))
        alpha1_details_5[a1].append(key)

print(f"α₁ distribution at n=5 (over {2**10} labeled tournaments):")
for a1 in sorted(alpha1_dist_5):
    print(f"  α₁={a1}: {alpha1_dist_5[a1]} tournaments")
print(f"\nα₁=3 count: {alpha1_dist_5[3]}")

# Show (c3, c5) breakdown for each α₁
print(f"\n(c3, c5) breakdown:")
for a1 in sorted(alpha1_details_5):
    if a1 <= 5:
        items = Counter(tuple(x) if isinstance(x, tuple) else x for x in alpha1_details_5[a1])
        print(f"  α₁={a1}: {dict(items)}")

# ─── Part 2: Exhaustive check at n=6 ───
print(f"\n{'='*60}")
print("PART 2: EXHAUSTIVE CHECK AT n=6")
print("="*60)

alpha1_dist_6 = Counter()
c3_to_alpha1_6 = defaultdict(Counter)
t0 = time.time()
for idx, T in enumerate(all_tournaments(6)):
    c3 = count_3cycles(T)
    # For efficiency, only compute full alpha1 for specific c3 values
    if c3 <= 4:  # low c3 — these are the ones where α₁=3 could happen
        a1, cycles = fast_alpha1(T)
        alpha1_dist_6[a1] += 1
        c3_to_alpha1_6[c3][a1] += 1
    else:
        # For high c3, α₁ >= c3 >= 5, can't be 3
        pass
    if idx % 10000 == 0 and idx > 0:
        elapsed = time.time() - t0
        print(f"  ...{idx}/{2**15} ({100*idx/2**15:.0f}%), elapsed={elapsed:.1f}s")

print(f"α₁ distribution at n=6 (c3 ≤ 4 only, over {sum(alpha1_dist_6.values())} tournaments):")
for a1 in sorted(alpha1_dist_6):
    print(f"  α₁={a1}: {alpha1_dist_6[a1]}")
print(f"\nα₁=3 among c3≤4: {alpha1_dist_6[3]}")

print(f"\nc3 -> α₁ map (c3 ≤ 4):")
for c3 in sorted(c3_to_alpha1_6):
    print(f"  c3={c3}: α₁ values = {dict(sorted(c3_to_alpha1_6[c3].items()))}")

# ─── Part 3: Sampling at n=7 ───
print(f"\n{'='*60}")
print("PART 3: SAMPLING AT n=7 (targeting α₁=3)")
print("="*60)

alpha1_dist_7 = Counter()
n_samples = 5000
n_with_c3_le_4 = 0
t0 = time.time()

for trial in range(n_samples):
    T = random_tournament(7)
    c3 = count_3cycles(T)
    if c3 <= 4:  # Only cases where α₁=3 is possible
        n_with_c3_le_4 += 1
        a1, cycles = fast_alpha1(T)
        alpha1_dist_7[a1] += 1
        if a1 == 3:
            print(f"  FOUND α₁=3! c3={c3}, cycles={cycles}")

elapsed = time.time() - t0
print(f"Sampled {n_samples} random tournaments at n=7 ({elapsed:.1f}s)")
print(f"  c3 ≤ 4: {n_with_c3_le_4} ({100*n_with_c3_le_4/n_samples:.1f}%)")
print(f"  α₁ distribution (c3≤4): {dict(sorted(alpha1_dist_7.items()))}")
print(f"  α₁=3 found: {alpha1_dist_7[3]}")

# ─── Part 4: WHY α₁=3 is forbidden ───
print(f"\n{'='*60}")
print("PART 4: STRUCTURAL ANALYSIS — WHY α₁=3 IS FORBIDDEN")
print("="*60)

# Key insight: α₁=3 requires exactly 3 odd cycles.
# The only way is c3=3, c5=c7=...=0  (since c3=2 gives α₁ ≤ 2+c5+...)
# or c3=2, c5=1, c7=c9=...=0
# or c3=1, c5=2, c7=...=0
# or c3=0, c5=3, c7=...=0
# etc.

# Let's check: at n=5, does c3=3 ALWAYS force c5≥1?
print("\nn=5: Testing c3=3 → c5≥1")
c3_3_with_c5_0 = 0
c3_3_total = 0
for T in all_tournaments(5):
    c3 = count_3cycles(T)
    if c3 == 3:
        c3_3_total += 1
        c5 = count_5cycles(T)
        if c5 == 0:
            c3_3_with_c5_0 += 1
print(f"  c3=3 tournaments: {c3_3_total}")
print(f"  c3=3 with c5=0: {c3_3_with_c5_0}")
print(f"  CONFIRMED: c3=3 → c5≥1 at n=5: {'YES' if c3_3_with_c5_0 == 0 else 'NO'}")

# At n=5: c3=2, can c5=1? (would give α₁=3)
print("\nn=5: Testing c3=2, c5=1 → α₁=3?")
c3_2_c5_1 = 0
for T in all_tournaments(5):
    c3 = count_3cycles(T)
    if c3 == 2:
        c5 = count_5cycles(T)
        if c5 == 1:
            c3_2_c5_1 += 1
print(f"  c3=2, c5=1: {c3_2_c5_1} tournaments")

# At n=5: c3=1, c5=2?
print("\nn=5: Testing c3=1, c5=2 → α₁=3?")
c3_1_c5_2 = 0
for T in all_tournaments(5):
    c3 = count_3cycles(T)
    if c3 == 1:
        c5 = count_5cycles(T)
        if c5 == 2:
            c3_1_c5_2 += 1
print(f"  c3=1, c5=2: {c3_1_c5_2} tournaments")

# At n=5: c3=0, c5=3?
print("\nn=5: Testing c3=0, c5=3 → α₁=3?")
c3_0_c5_3 = 0
for T in all_tournaments(5):
    c3 = count_3cycles(T)
    if c3 == 0:
        c5 = count_5cycles(T)
        if c5 == 3:
            c3_0_c5_3 += 1
print(f"  c3=0, c5=3: {c3_0_c5_3} tournaments")

# ─── Part 5: The (c3, c5, c7) → α₁ map at n=6 ───
print(f"\n{'='*60}")
print("PART 5: (c3, c5, c7) → α₁ MAP AT n=6 (sampling)")
print("="*60)

cycle_map_6 = defaultdict(Counter)
for trial in range(50000):
    T = random_tournament(6)
    c3 = count_3cycles(T)
    if c3 <= 4:
        a1, cycles = fast_alpha1(T)
        c5 = cycles.get(5, 0)
        cycle_map_6[(c3, c5)][a1] += 1

print("(c3, c5) -> α₁ distribution:")
for key in sorted(cycle_map_6):
    dist = dict(sorted(cycle_map_6[key].items()))
    print(f"  (c3, c5)={key}: α₁={dist}")

# ─── Part 6: Can c3=3 ever have c5=0 at n=6? ───
print(f"\n{'='*60}")
print("PART 6: c3=3, c5=0 AT n=6 (exhaustive)")
print("="*60)

c3_3_c5_0_n6 = 0
c3_3_total_n6 = 0
for T in all_tournaments(6):
    c3 = count_3cycles(T)
    if c3 == 3:
        c3_3_total_n6 += 1
        c5 = count_5cycles(T)
        if c5 == 0:
            c3_3_c5_0_n6 += 1
print(f"c3=3 tournaments at n=6: {c3_3_total_n6}")
print(f"c3=3, c5=0 at n=6: {c3_3_c5_0_n6}")

# ─── Part 7: What about H=21? ───
print(f"\n{'='*60}")
print("PART 7: WHY IS H=21 FORBIDDEN?")
print("="*60)
# H = 1 + 2α₁ + 4α₂ + 8α₃ + ... = 21
# → 2α₁ + 4α₂ + 8α₃ + ... = 20
# Possible: α₁=10, rest 0 (I=1+10x, H=21)
#   or α₁=8, α₂=1 (I=1+8x+x², H=1+16+4=21)
#   or α₁=6, α₂=2 (I=1+6x+2x², H=1+12+8=21)
#   or α₁=4, α₂=3 (I=1+4x+3x², H=1+8+12=21)
#   or α₁=2, α₂=4 (H=1+4+16=21)
#   or α₁=0, α₂=5 (H=1+0+20=21)
# But also: α₁=10, α₂=0, α₃=0 or α₁=6, α₂=1, α₃=1 (1+12+4+8=25 no)
# Let me enumerate systematically:
print("Decompositions of H=21 = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
print("Need: 2α₁ + 4α₂ + 8α₃ + 16α₄ + ... = 20")
print("With constraint: α_k ≤ C(α₁, k) (α_k counts independent k-sets)")
solutions = []
for a1 in range(21):
    rem = 20 - 2*a1
    if rem < 0:
        break
    if rem == 0:
        solutions.append((a1, 0, 0))
        continue
    for a2 in range(rem//4 + 1):
        rem2 = rem - 4*a2
        if rem2 < 0:
            break
        if rem2 == 0:
            if a2 <= a1*(a1-1)//2:  # α₂ ≤ C(α₁,2)
                solutions.append((a1, a2, 0))
            continue
        for a3 in range(rem2//8 + 1):
            rem3 = rem2 - 8*a3
            if rem3 == 0:
                if a2 <= a1*(a1-1)//2:
                    solutions.append((a1, a2, a3))

print(f"Possible (α₁, α₂, α₃) for H=21:")
for s in solutions:
    H_check = 1 + 2*s[0] + 4*s[1] + 8*s[2]
    print(f"  α₁={s[0]}, α₂={s[1]}, α₃={s[2]} → H={H_check}")

# Check at n=6: which α₁ values are achievable?
print(f"\nAchievable α₁ at n=6 (exhaustive, sampling high c3):")
alpha1_full_6 = Counter()
for T in all_tournaments(6):
    a1, _ = fast_alpha1(T)
    alpha1_full_6[a1] += 1
print(f"α₁ values: {sorted(alpha1_full_6.keys())}")
print(f"Full distribution: {dict(sorted(alpha1_full_6.items()))}")

# For H=21: need α₁=10 with α₂=0. Does α₁=10 ever have α₂=0?
# Or α₁=8, α₂=1. Etc.
# From the polynomial data at n=6: α₁=10 appears with I=1+10x+2x², giving H=29.
# So α₁=10 with α₂=0 would give I=1+10x, H=21.
# But from the data, I=1+10x never appears! The only α₁=10 entry has α₂=2.
# So the question is: can 10 odd cycles all be pairwise adjacent?

# ─── Part 8: Deeper structure at n=7 ───
print(f"\n{'='*60}")
print("PART 8: TARGETED n=7 — α₁ near 3")
print("="*60)

# Generate tournaments with small c3 at n=7 and count all odd cycles
alpha1_dist_7_targeted = Counter()
c3_alpha1_7 = defaultdict(Counter)
n_7_targeted = 50000
for trial in range(n_7_targeted):
    T = random_tournament(7)
    c3 = count_3cycles(T)
    if c3 <= 5:
        a1, cycles = fast_alpha1(T)
        alpha1_dist_7_targeted[a1] += 1
        c3_alpha1_7[c3][a1] += 1
        if a1 == 3:
            print(f"  *** FOUND α₁=3 at n=7! c3={c3}, cycles={cycles} ***")

print(f"\nTargeted n=7 sampling ({n_7_targeted} trials, c3≤5 subset):")
print(f"α₁ distribution: {dict(sorted(alpha1_dist_7_targeted.items()))}")
print(f"\nc3 -> α₁ at n=7:")
for c3 in sorted(c3_alpha1_7):
    print(f"  c3={c3}: {dict(sorted(c3_alpha1_7[c3].items()))}")

print("\n\nDONE.")
