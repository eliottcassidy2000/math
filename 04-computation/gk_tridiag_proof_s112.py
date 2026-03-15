#!/usr/bin/env python3
"""
gk_tridiag_proof_s112.py — Why g_k is always cubic for k>=3
kind-pasteur-2026-03-15-S112

Key idea: Z_j = X_j - Y_j where X_j = 1[sigma(j+1)=sigma(j)+1],
Y_j = 1[sigma(j+1)=sigma(j)-1]. The covariance matrix is TRIDIAGONAL.

Hypothesis: the cubic universality follows from the Markov/tridiagonal
structure of the Z_j process. Specifically, for a 1D nearest-neighbor
process, connected correlation functions (cumulants) of order >= 4
that span > 3 consecutive sites must vanish.

Approach: directly compute E[prod_{j in S} Z_j] for all domino subsets S
by explicit enumeration of permutations, and analyze the structure.
"""

from itertools import permutations
from fractions import Fraction
from math import comb, factorial
from collections import defaultdict

def compute_expectations(n, max_k=None):
    """Compute E[prod_{j in S} Z_j] for all domino subsets S in [0..n-2].

    A domino subset is a union of adjacent pairs {j, j+1}.
    Return dict: frozenset(S) -> E[prod Z_j].
    """
    m = n - 1  # number of positions: indices 0 to n-2
    nfact = factorial(n)

    # Adjacent pairs: (j, j+1) for j = 0, ..., n-3
    pairs = [(j, j+1) for j in range(m - 1)]  # m-1 = n-2 pairs

    if max_k is None:
        max_k = len(pairs)

    # Compute Z values for all permutations
    z_vals = {}
    for perm in permutations(range(n)):
        z = []
        for j in range(m):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z.append(x - y)
        z_vals[perm] = z

    # Generate all subsets of pairs (up to max_k pairs)
    results = {}

    def recurse(start, chosen):
        # Compute expectation for current chosen set
        if chosen:
            positions = set()
            for p in chosen:
                positions.update(p)
            S = frozenset(positions)

            if S not in results:
                total = Fraction(0)
                for perm, z in z_vals.items():
                    prod = 1
                    for j in sorted(S):
                        prod *= z[j]
                    total += prod
                results[S] = total / nfact

        if len(chosen) >= max_k:
            return

        for i in range(start, len(pairs)):
            chosen.append(pairs[i])
            recurse(i + 1, chosen)
            chosen.pop()

    recurse(0, [])
    return results

def classify_domino_subset(S):
    """Classify a domino subset by its cluster structure.
    Returns tuple of cluster sizes."""
    positions = sorted(S)
    if not positions:
        return ()
    clusters = []
    current = [positions[0]]
    for p in positions[1:]:
        if p == current[-1] + 1:
            current.append(p)
        else:
            clusters.append(len(current))
            current = [p]
    clusters.append(len(current))
    return tuple(clusters)

def is_domino_subset(S):
    """Check if S is a valid union of adjacent pairs."""
    positions = sorted(S)
    # Each cluster must have even size (pairs of 2)
    # Wait no — a cluster of size 2s-1 represents s overlapping pairs
    # Actually, a "domino subset" is any subset that can be written as union of adjacent pairs
    # This means: for every position j in S, either j-1 or j+1 must also be in S
    for j in positions:
        if j-1 not in S and j+1 not in S:
            return False
    return True

print("="*70)
print("DIRECT COMPUTATION: E[prod Z_j] for domino subsets")
print("="*70)

for n in range(3, 8):
    print(f"\n--- n = {n} ---")
    max_k = (n-1) // 2
    exps = compute_expectations(n, max_k)

    # Group by |S| (=2k)
    by_size = defaultdict(list)
    for S, val in exps.items():
        if val != 0 and is_domino_subset(S):
            by_size[len(S)].append((S, val))

    for size in sorted(by_size.keys()):
        items = by_size[size]
        k = size // 2
        total = sum(v for _, v in items)
        print(f"  |S|={size} (k={k}): {len(items)} domino subsets, total = {total}")

        # Group by cluster type
        by_type = defaultdict(lambda: (0, Fraction(0)))
        for S, v in items:
            ctype = classify_domino_subset(S)
            count, tot = by_type[ctype]
            by_type[ctype] = (count + 1, tot + v)

        for ctype in sorted(by_type.keys()):
            count, tot = by_type[ctype]
            avg = tot / count if count > 0 else 0
            print(f"    cluster type {ctype}: {count} subsets, total={tot}, avg={avg}")

# Now the key test: for each n, compute the total |S|=2k contribution
# and verify it matches 2*g_k(n-2k)/(n)_{2k}
print("\n" + "="*70)
print("VERIFY g_k VALUES")
print("="*70)

def falling_factorial(n, k):
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

def g_true(k, m):
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k == 3: return Fraction(2*m**3 + m, 3)
    if k == 4: return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)
    return None

for n in range(3, 8):
    print(f"\nn={n}:")
    max_k = (n-1) // 2
    exps = compute_expectations(n, max_k)

    for k in range(1, max_k + 1):
        size = 2 * k
        total = sum(v for S, v in exps.items() if len(S) == size and is_domino_subset(S))
        m = n - 2*k
        predicted = 2 * g_true(k, m) / falling_factorial(n, 2*k) if g_true(k, m) is not None else None
        match = "OK" if predicted is not None and total == predicted else f"pred={predicted}"
        print(f"  k={k}: total={total}, predicted={match}")

# Key test: connected vs disconnected cluster contributions
print("\n" + "="*70)
print("CONNECTED CORRELATION ANALYSIS")
print("="*70)
print("For each cluster type, check if E[prod Z_j] factorizes.")

for n in range(5, 8):
    print(f"\nn={n}:")
    exps = compute_expectations(n, 3)

    # Find pairs of non-overlapping clusters
    # For each subset with 2+ clusters, check factorization
    for S, val in sorted(exps.items(), key=lambda x: len(x[0])):
        if not is_domino_subset(S) or val == 0:
            continue
        ctype = classify_domino_subset(S)
        if len(ctype) < 2:
            continue

        # Factor into cluster expectations
        positions = sorted(S)
        clusters = []
        current = [positions[0]]
        for p in positions[1:]:
            if p == current[-1] + 1:
                current.append(p)
            else:
                clusters.append(frozenset(current))
                current = [p]
        clusters.append(frozenset(current))

        # Product of individual cluster expectations
        cluster_exps = [exps.get(c, Fraction(0)) for c in clusters]
        product = Fraction(1)
        for ce in cluster_exps:
            product *= ce

        # Gap between clusters
        cluster_ranges = [(min(c), max(c)) for c in clusters]
        gaps = [cluster_ranges[i+1][0] - cluster_ranges[i][1] - 1 for i in range(len(clusters)-1)]

        ratio = val / product if product != 0 else "undef"
        if ratio != 1:
            print(f"  S={set(S)}, type={ctype}, gaps={gaps}: E={val}, prod={product}, ratio={ratio}")

print("\nDone!")
