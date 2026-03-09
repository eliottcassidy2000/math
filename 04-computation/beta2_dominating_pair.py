#!/usr/bin/env python3
"""
beta2_dominating_pair.py — Does every tournament have a dominating pair?

A dominating pair (v₁, v₂) has N⁺(v₁) ∪ N⁺(v₂) ∪ {v₁,v₂} = [n].
Equivalently: for every u ∉ {v₁,v₂}, either v₁→u or v₂→u.

If YES + full coverage implies cone success, then β₂ = 0 follows!

Author: opus-2026-03-08-S49
"""
import sys, time, random
from itertools import combinations
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


print("=" * 70)
print("DOMINATING PAIR EXISTENCE")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    has_dom_pair = 0
    no_dom_pair = 0
    min_dom_size = Counter()

    for bits in range(total):
        A = build_adj(n, bits)

        # Check all pairs
        found_pair = False
        for v1, v2 in combinations(range(n), 2):
            coverage = {v1, v2}
            for u in range(n):
                if A[v1][u] or A[v2][u]:
                    coverage.add(u)
            if len(coverage) == n:
                found_pair = True
                break

        if found_pair:
            has_dom_pair += 1
        else:
            no_dom_pair += 1
            # Check if single vertex dominates
            single = any(sum(A[v][u] for u in range(n) if u != v) == n-1 for v in range(n))
            if single:
                min_dom_size[1] += 1
            else:
                # Check triples
                found_triple = False
                for v1, v2, v3 in combinations(range(n), 3):
                    coverage = {v1, v2, v3}
                    for u in range(n):
                        if A[v1][u] or A[v2][u] or A[v3][u]:
                            coverage.add(u)
                    if len(coverage) == n:
                        found_triple = True
                        break
                if found_triple:
                    min_dom_size[3] += 1
                else:
                    min_dom_size['>3'] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.1f}s")
    print(f"  Has dominating pair: {has_dom_pair} ({100*has_dom_pair/total:.1f}%)")
    print(f"  No dominating pair: {no_dom_pair}")
    if no_dom_pair > 0:
        print(f"  Of those without pair: {dict(min_dom_size)}")

# Now test larger n by sampling
print(f"\n{'='*70}")
print("SAMPLED: DOMINATING PAIR AT n=7,8,9,10")
print("=" * 70)

for n in [7, 8, 9, 10, 15, 20]:
    samples = 10000 if n <= 10 else 1000
    t0 = time.time()
    has_pair = 0

    for _ in range(samples):
        A = build_random_adj(n)

        found = False
        for v1, v2 in combinations(range(n), 2):
            ok = True
            for u in range(n):
                if u == v1 or u == v2:
                    continue
                if not (A[v1][u] or A[v2][u]):
                    ok = False
                    break
            if ok:
                found = True
                break

        if found:
            has_pair += 1

    elapsed = time.time() - t0
    pct = 100*has_pair/samples
    print(f"  n={n}: {has_pair}/{samples} ({pct:.1f}%) have dominating pair ({elapsed:.1f}s)")

# Theoretical analysis
print(f"\n{'='*70}")
print("THEORETICAL: DOMINATING PAIR PROBABILITY")
print("=" * 70)
print("""
For a random tournament, the probability that vertex u is NOT dominated
by either v₁ or v₂ is P(u→v₁ AND u→v₂) = 1/4.

For a SPECIFIC pair (v₁,v₂), P(pair dominates) = (1 - 1/4)^{n-2} = (3/4)^{n-2}.

Number of pairs: C(n,2). So:
P(∃ dominating pair) ≥ 1 - C(n,2) * (1-(3/4)^{n-2})
Wait, that's not right. Let me use inclusion-exclusion or union bound.

P(NO dominating pair) ≤ Σ_{pairs} P(pair fails) - ... (by inclusion-exclusion)
Actually: P(NO dominating pair) ≤ Σ P(pair fails) - ... but this isn't tight.

Better: P(specific pair fails) = 1 - (3/4)^{n-2}.
P(ALL pairs fail) ≤ (1 - (3/4)^{n-2})^{C(n,2)} → 0 as n → ∞.

But this is only a probabilistic argument. We need a DETERMINISTIC guarantee.
""")

# Check: does every tournament on n ≥ 3 have a dominating pair?
# From the exhaustive data above, this should be clear.
# If not, check what fails.

print("CHECKING: Does the 3-cycle C₃ have a dominating pair?")
A_c3 = [[0,1,0],[0,0,1],[1,0,0]]  # 0→1→2→0
for v1, v2 in combinations(range(3), 2):
    coverage = {v1, v2}
    for u in range(3):
        if A_c3[v1][u] or A_c3[v2][u]:
            coverage.add(u)
    print(f"  ({v1},{v2}): coverage = {coverage}")

print("\nDone.")
