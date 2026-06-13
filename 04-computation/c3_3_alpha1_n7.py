#!/usr/bin/env python3
"""
At n=7 with c3=3: verify that alpha_1 != 3 still holds,
even though triples can span > 5 vertices.

The question: does c3=3 always force some c_k >= 1 for k >= 5?
Or more precisely: is alpha_1 = c3 + c5 + c7 ever equal to 3?

Also: for the cases where triples span > 5 vertices, what is c5?
Specifically: when triples span 6 or 7 vertices, does c5 > 0?

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

def c3_from_scores(scores, n):
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def count_odd_cycles(T, n, max_k=None):
    """Count directed cycles of each odd length."""
    if max_k is None:
        max_k = n
    result = {}
    for k in range(3, max_k+1, 2):
        total = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
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
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    total += c
            result[k] = total
    return result

# Focus: n=7 exhaustive check
n = 7
m = n*(n-1)//2  # 21
total_bits = 1 << m  # 2097152

print(f"n=7: Exhaustive search ({total_bits} tournaments)")
print(f"Checking if alpha_1=3 ever occurs...")

c3_3_count = 0
alpha1_3_count = 0
c3_3_span_counts = {}
c3_3_odd_cycle_patterns = {}

# For c3=3 cases with large span, check c5
c3_3_with_c5_0 = 0
c3_3_with_c5_0_examples = []

for bits in range(total_bits):
    if bits % 500000 == 0:
        print(f"  Progress: {bits}/{total_bits} ({100*bits//total_bits}%)", flush=True)

    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores, n)

    if c3 != 3:
        continue

    c3_3_count += 1
    triples = get_cyclic_triples(T, n)
    assert len(triples) == 3

    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)
    c3_3_span_counts[span] = c3_3_span_counts.get(span, 0) + 1

    # Count c5 and c7
    cycles = count_odd_cycles(T, n)
    c5 = cycles.get(5, 0)
    c7 = cycles.get(7, 0)
    alpha1 = 3 + c5 + c7  # c3 = 3

    pattern = (c5, c7)
    c3_3_odd_cycle_patterns[pattern] = c3_3_odd_cycle_patterns.get(pattern, 0) + 1

    if alpha1 == 3:
        alpha1_3_count += 1
        print(f"  ALPHA1=3 FOUND! bits={bits}, c5={c5}, c7={c7}")

    if c5 == 0:
        c3_3_with_c5_0 += 1
        if len(c3_3_with_c5_0_examples) < 5:
            common = triples[0] & triples[1] & triples[2]
            c3_3_with_c5_0_examples.append((bits, span, len(common), [sorted(t) for t in triples], c7))

print(f"\nResults:")
print(f"  Total tournaments with c3=3: {c3_3_count}")
print(f"  alpha_1=3 found: {alpha1_3_count}")
print(f"  Span distribution: {c3_3_span_counts}")
print(f"  (c5, c7) patterns: {c3_3_odd_cycle_patterns}")
print(f"  c3=3 with c5=0: {c3_3_with_c5_0}")
if c3_3_with_c5_0_examples:
    print(f"  Examples with c5=0:")
    for bits, span, common_sz, trips, c7 in c3_3_with_c5_0_examples:
        print(f"    bits={bits}, span={span}, common={common_sz}, triples={trips}, c7={c7}")

# Also check: does alpha_1=3 occur for ANY c3 value?
print(f"\nChecking ALL tournaments for alpha_1=3...")
alpha1_3_any = 0
for bits in range(total_bits):
    if bits % 500000 == 0 and bits > 0:
        print(f"  Progress: {bits}/{total_bits}", flush=True)
    T = tournament_from_bits(n, bits)
    cycles = count_odd_cycles(T, n)
    a1 = sum(cycles.values())
    if a1 == 3:
        alpha1_3_any += 1
        if alpha1_3_any <= 3:
            print(f"  FOUND alpha_1=3: bits={bits}, cycles={cycles}")

print(f"  Total with alpha_1=3: {alpha1_3_any}")
print("Done.")
