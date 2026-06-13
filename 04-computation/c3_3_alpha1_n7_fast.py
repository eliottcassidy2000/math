#!/usr/bin/env python3
"""
Fast n=7 verification:
1. Among c3=3 tournaments: what are c5, c7? Does alpha_1=3 occur?
2. Among ALL tournaments: does alpha_1=3 ever occur?

Optimized: use Moon's formula to skip non-c3=3 tournaments fast.
For part 2: only compute full cycle counts for tournaments where c3<=3.

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations, permutations

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
    """Count directed Hamiltonian cycles on the given vertex subset."""
    k = len(verts)
    if k < 3:
        return 0
    # DP over Hamiltonian paths starting at verts[0]
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
    """Count directed k-cycles in T."""
    total = 0
    for combo in combinations(range(n), k):
        total += count_directed_cycles_on_subset(T, list(combo))
    return total

# ============================================================
# PART 1: c3=3 at n=7
# ============================================================
n = 7
m = n*(n-1)//2  # 21
total_bits = 1 << m

print(f"n=7: Scanning {total_bits} tournaments for c3=3...")

c3_3_count = 0
span_dist = {}
common_vertex_dist = {}
c5_dist = {}
c7_dist = {}
alpha1_values = {}

for bits in range(total_bits):
    if bits % 500000 == 0:
        print(f"  {bits}/{total_bits} ({100*bits//total_bits}%)", flush=True)

    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 != 3:
        continue

    c3_3_count += 1
    triples = get_cyclic_triples(T, n)

    # Span
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)
    span_dist[span] = span_dist.get(span, 0) + 1

    # Common vertex
    common = triples[0] & triples[1] & triples[2]
    cv = len(common)
    common_vertex_dist[cv] = common_vertex_dist.get(cv, 0) + 1

    # c5, c7
    c5 = count_ck(T, n, 5)
    c7 = count_ck(T, n, 7)

    c5_dist[c5] = c5_dist.get(c5, 0) + 1
    c7_dist[c7] = c7_dist.get(c7, 0) + 1

    a1 = 3 + c5 + c7
    alpha1_values[a1] = alpha1_values.get(a1, 0) + 1

print(f"\nc3=3 results at n=7:")
print(f"  Count: {c3_3_count}")
print(f"  Span dist: {span_dist}")
print(f"  Common vertex sizes: {common_vertex_dist}")
print(f"  c5 dist: {dict(sorted(c5_dist.items()))}")
print(f"  c7 dist: {dict(sorted(c7_dist.items()))}")
print(f"  alpha_1 values: {dict(sorted(alpha1_values.items()))}")
print(f"  alpha_1=3 occurs? {3 in alpha1_values}")

# ============================================================
# PART 2: Does alpha_1=3 occur for ANY tournament at n=7?
# ============================================================
print(f"\nPart 2: Checking all n=7 tournaments for alpha_1=3...")
print(f"  Strategy: alpha_1 = c3 + c5 + c7. Need c3+c5+c7=3.")
print(f"  If c3>=4, impossible. Only check c3 <= 3.")

alpha1_3_total = 0

for bits in range(total_bits):
    if bits % 500000 == 0:
        print(f"  {bits}/{total_bits} ({100*bits//total_bits}%)", flush=True)

    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)

    if c3 > 3:
        continue

    # Need c5 + c7 = 3 - c3
    target = 3 - c3
    if target < 0:
        continue

    c5 = count_ck(T, n, 5)
    if c5 > target:
        continue

    c7 = count_ck(T, n, 7)
    if c5 + c7 == target:
        alpha1_3_total += 1
        if alpha1_3_total <= 3:
            print(f"  FOUND: bits={bits}, c3={c3}, c5={c5}, c7={c7}")

print(f"\n  Total n=7 tournaments with alpha_1=3: {alpha1_3_total}")
print("Done.")
