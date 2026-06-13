#!/usr/bin/env python3
"""
Analyze the n=7 tournaments with alpha_1=3 (c3=3, c5=0, c7=0).

Questions:
1. What are their score sequences?
2. What is the structure of the 3 disjoint-ish triples?
3. Can we characterize them?
4. Verify: at n=5,6 alpha_1=3 truly never occurs (recap).

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def c3_from_scores(scores):
    n = len(scores)
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

def count_directed_cycles_on_subset(T, verts):
    k = len(verts)
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for vi in range(k):
            if not (mask & (1 << vi)):
                continue
            c = dp.get((mask, vi), 0)
            if c == 0:
                continue
            for ui in range(k):
                if mask & (1 << ui):
                    continue
                if T[verts[vi]][verts[ui]]:
                    key = (mask | (1 << ui), ui)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << k) - 1
    total = 0
    for vi in range(1, k):
        c = dp.get((full, vi), 0)
        if c > 0 and T[verts[vi]][verts[0]]:
            total += c
    return total

def count_ck(T, n, k):
    total = 0
    for combo in combinations(range(n), k):
        total += count_directed_cycles_on_subset(T, list(combo))
    return total

# ============================================================
# Analyze n=7 alpha_1=3 tournaments
# ============================================================
n = 7
m = 21
total_bits = 1 << m

print(f"n=7: Finding all tournaments with alpha_1=3...")

examples = []
score_dist = {}
triple_patterns = {}

for bits in range(total_bits):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 != 3:
        continue

    triples = get_cyclic_triples(T, n)
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)

    if span != 7:  # Only span=7 can have c5=0
        continue

    # Verify c5=0
    c5 = count_ck(T, n, 5)
    if c5 != 0:
        continue

    score_seq = tuple(sorted(scores))
    score_dist[score_seq] = score_dist.get(score_seq, 0) + 1

    # Analyze triple intersection pattern
    pw = []
    for i in range(3):
        for j in range(i+1, 3):
            pw.append(len(triples[i] & triples[j]))
    pattern = tuple(sorted(pw))
    triple_patterns[pattern] = triple_patterns.get(pattern, 0) + 1

    if len(examples) < 5:
        examples.append((bits, score_seq, [sorted(t) for t in triples], pattern))

print(f"\nTotal alpha_1=3 at n=7: {sum(score_dist.values())}")
print(f"\nScore sequence distribution:")
for s, cnt in sorted(score_dist.items()):
    print(f"  {s}: {cnt}")

print(f"\nTriple intersection patterns (pairwise sizes):")
for p, cnt in sorted(triple_patterns.items()):
    print(f"  {p}: {cnt}")

print(f"\nExamples:")
for bits, scores, trips, pw in examples:
    print(f"  bits={bits}, scores={scores}")
    print(f"    triples={trips}, pw_intersections={pw}")
    # Show the adjacency
    T = tournament_from_bits(n, bits)
    print(f"    Adjacency:")
    for i in range(n):
        print(f"      {i}: beats {[j for j in range(n) if T[i][j]]}")

# ============================================================
# Key question: do the 3 triples include a DISJOINT pair?
# ============================================================
print(f"\nDo any of the 3360 examples have a pair of disjoint triples?")
disjoint_count = 0
for bits in range(total_bits):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 != 3:
        continue
    triples = get_cyclic_triples(T, n)
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    if len(all_verts) != 7:
        continue

    has_disjoint = False
    for i in range(3):
        for j in range(i+1, 3):
            if len(triples[i] & triples[j]) == 0:
                has_disjoint = True
                break
        if has_disjoint:
            break
    if has_disjoint:
        disjoint_count += 1

print(f"  Tournaments with a disjoint triple pair: {disjoint_count}")

# ============================================================
# Verify at n=5,6: alpha_1=3 truly impossible
# ============================================================
for n in [5, 6]:
    m = n*(n-1)//2
    found = False
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i]) for i in range(n)]
        c3 = c3_from_scores(scores)
        if c3 > 3:
            continue
        # Check all odd cycles
        total = c3
        for k in range(5, n+1, 2):
            total += count_ck(T, n, k)
        if total == 3:
            found = True
            break
    print(f"\nn={n}: alpha_1=3 found? {found}")

print("\nDone.")
