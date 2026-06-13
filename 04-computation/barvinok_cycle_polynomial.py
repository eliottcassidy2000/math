#!/usr/bin/env python3
"""
barvinok_cycle_polynomial.py — Cycle cover polynomial for tournaments

Inspired by Barvinok (2014) "On Testing Hamiltonicity of Graphs".

Key insight: The permanent per(A) = Σ_σ Π a_{iσ(i)} decomposes by
number of cycles in σ:

    per(A) = Σ_{k=1}^{n} C_k(T)

where C_k(T) = Σ_{σ: c(σ)=k} Π a_{iσ(i)} counts cycle covers with exactly k cycles.

The CYCLE COVER POLYNOMIAL is:
    P_T(x) = Σ_{k=1}^{n} C_k(T) · x^k

So P_T(1) = per(A) and C_1(T) = ham(A) = number of Hamiltonian circuits.

BARVINOK'S THEOREM: For "softened" matrices (ε ≤ a_{ij} ≤ 1),
    (1/2r)(ε/n)^r · per(A) ≤ ham(A) ≤ per(A)
meaning the SINGLE-CYCLE contribution dominates the permanent.

QUESTION: For actual tournament adjacency matrices (ε=0), how does
the cycle cover polynomial behave? How does C_1/per vary across
isomorphism classes?

Also exploring:
- The "softened" permanent per(A_ε) as ε varies from 0 to 1
- Sinkhorn scaling of tournament matrices
- Connection to our H(T) = I(Ω(T), 2) formula
"""

import sys
import time
from itertools import permutations
from collections import defaultdict, Counter
from math import factorial, gcd
import numpy as np

# ================================================================
# CYCLE STRUCTURE OF PERMUTATIONS
# ================================================================

def cycle_count(sigma):
    """Count cycles in a permutation."""
    n = len(sigma)
    visited = [False] * n
    count = 0
    for i in range(n):
        if visited[i]:
            continue
        count += 1
        j = i
        while not visited[j]:
            visited[j] = True
            j = sigma[j]
    return count

def cycle_type(sigma):
    """Return sorted cycle lengths of a permutation."""
    n = len(sigma)
    visited = [False] * n
    lengths = []
    for i in range(n):
        if visited[i]:
            continue
        length = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = sigma[j]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


# ================================================================
# TOURNAMENT ENUMERATION
# ================================================================

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(arcs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj, bits

def canonicalize(adj, n, perms):
    """Canonical form under vertex relabeling."""
    best = None
    for p in perms:
        s = ''.join(str(adj[p[i]][p[j]]) for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def count_H_dp(adj, n):
    """Count Hamiltonian paths via Held-Karp DP."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)) or dp[S][v] == 0:
                continue
            for w in range(n):
                if not (S & (1 << w)) and adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full])


# ================================================================
# CYCLE COVER POLYNOMIAL
# ================================================================

def cycle_cover_polynomial(adj, n):
    """Compute the cycle cover polynomial P_T(x) = Σ C_k x^k.

    C_k = sum over permutations σ with exactly k cycles
           of Π_{i} A[i, σ(i)].

    Returns dict: k → C_k
    """
    coeffs = defaultdict(int)
    for sigma in permutations(range(n)):
        # Compute product Π A[i, σ(i)]
        product = 1
        for i in range(n):
            product *= adj[i][sigma[i]]
            if product == 0:
                break
        if product == 0:
            continue
        k = cycle_count(sigma)
        coeffs[k] += product
    return dict(coeffs)


def softened_permanent(adj, n, epsilon):
    """Compute permanent of softened matrix A_ε.

    A_ε[i][j] = A[i][j] if A[i][j]=1, else ε.
    (Diagonal stays 0 since i→i is not an arc.)
    """
    # Build softened matrix
    A_eps = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                A_eps[i][j] = epsilon  # or 0? Barvinok uses ε for ALL entries
            elif adj[i][j] == 1:
                A_eps[i][j] = 1.0
            else:
                A_eps[i][j] = epsilon

    # Compute permanent via Ryser's formula
    result = 0.0
    for mask in range(1, 1 << n):
        col_sums = [0.0] * n
        bits = 0
        m2 = mask
        while m2:
            j = (m2 & (-m2)).bit_length() - 1
            for i in range(n):
                col_sums[i] += A_eps[i][j]
            bits += 1
            m2 &= m2 - 1
        prod = 1.0
        for i in range(n):
            prod *= col_sums[i]
        result += ((-1)**(n - bits)) * prod
    return result * ((-1)**n)


def softened_ham(adj, n, epsilon):
    """Compute Hamiltonian permanent of softened matrix.

    ham(A_ε) = Σ_{σ∈H_n} Π A_ε[i,σ(i)]
    where H_n = single-cycle permutations.
    """
    result = 0.0
    for sigma in permutations(range(n)):
        if cycle_count(sigma) != 1:
            continue
        prod = 1.0
        for i in range(n):
            j = sigma[i]
            if i == j:
                prod *= epsilon
            elif adj[i][j]:
                prod *= 1.0
            else:
                prod *= epsilon
        result += prod
    return result


# ================================================================
# MAIN: COMPUTE FOR ALL TOURNAMENT CLASSES
# ================================================================

print("CYCLE COVER POLYNOMIAL FOR TOURNAMENT ISOMORPHISM CLASSES")
print("Inspired by Barvinok (2014)")
print("=" * 75)

for n in range(3, 7):
    t0 = time.time()
    m = n * (n - 1) // 2
    perms_list = list(permutations(range(n)))

    # Group tournaments by isomorphism class
    classes = {}
    for adj, bits in all_tournaments(n):
        canon = canonicalize(adj, n, perms_list)
        if canon not in classes:
            classes[canon] = adj

    class_data = []
    for canon in sorted(classes.keys()):
        adj = classes[canon]
        H = count_H_dp(adj, n)
        scores = score_seq(adj, n)
        poly = cycle_cover_polynomial(adj, n)
        per_A = sum(poly.values())
        ham_A = poly.get(1, 0)  # C_1 = Hamiltonian circuits

        class_data.append({
            'adj': adj, 'H': H, 'scores': scores,
            'poly': poly, 'per': per_A, 'ham': ham_A,
        })

    elapsed = time.time() - t0
    print(f"\nn={n}  ({len(class_data)} classes, {elapsed:.2f}s)")
    print(f"{'#':>3} {'H':>5} {'per(A)':>8} {'ham(A)':>8} {'ham/per':>8} {'scores':>16} {'C_k distribution'}")
    print("-" * 90)

    for ci, d in enumerate(sorted(class_data, key=lambda x: x['H'])):
        ratio = d['ham'] / d['per'] if d['per'] > 0 else 0
        poly_str = ' '.join(f'C{k}={v}' for k, v in sorted(d['poly'].items()))
        print(f"{ci:3d} {d['H']:5d} {d['per']:8d} {d['ham']:8d} {ratio:8.4f} {str(d['scores']):>16} {poly_str}")

    # Aggregate statistics
    total_per = sum(d['per'] for d in class_data)
    total_ham = sum(d['ham'] for d in class_data)
    print(f"\n  Total per(A) across classes: {total_per}")
    print(f"  Total ham(A) across classes: {total_ham}")
    if total_per > 0:
        print(f"  Global ham/per ratio: {total_ham/total_per:.4f}")

    # Check Barvinok's bound for softened matrices
    if n <= 5:
        print(f"\n  BARVINOK SOFTENING TEST (ε = n^{{-0.1}} = {n**(-0.1):.4f}):")
        eps = n ** (-0.1)
        r = int(4 * np.log(n) / eps**2) + 6
        print(f"  r = ⌊4 ln n / ε²⌋ + 6 = {r}")
        print(f"  Bound factor: (1/2r)(ε/n)^r = {1/(2*r) * (eps/n)**r:.2e}")

        for ci, d in enumerate(sorted(class_data, key=lambda x: x['H'])[:3]):
            per_eps = softened_permanent(d['adj'], n, eps)
            ham_eps = softened_ham(d['adj'], n, eps)
            bound = 1/(2*r) * (eps/n)**r * per_eps
            print(f"    Class #{ci}: per_ε={per_eps:.4f}, ham_ε={ham_eps:.4f}, "
                  f"bound={bound:.2e}, ham/per={ham_eps/per_eps:.4f}" if per_eps > 0 else
                  f"    Class #{ci}: per_ε={per_eps:.4f}")

    # C_1 / per(A) statistics
    ratios = [d['ham'] / d['per'] for d in class_data if d['per'] > 0]
    if ratios:
        print(f"\n  ham/per ratio: min={min(ratios):.4f}, max={max(ratios):.4f}, "
              f"mean={sum(ratios)/len(ratios):.4f}")

    # Which cycle count dominates the permanent?
    dominant = Counter()
    for d in class_data:
        if d['poly']:
            max_k = max(d['poly'], key=d['poly'].get)
            dominant[max_k] += 1
    print(f"  Dominant cycle count in per(A): {dict(sorted(dominant.items()))}")

print(f"""

INTERPRETATION:
  C_1(T) = ham(A) = number of Hamiltonian CIRCUITS in T
  H(T) = number of Hamiltonian PATHS in T
  per(A) = total cycle covers = C_1 + C_2 + ... + C_n

  Barvinok's theorem: for softened matrices, ham ≈ per up to O(exp(n^0.4)).
  For actual 0-1 tournament matrices: the ratio ham/per varies by class.

  KEY QUESTION: Is ham/per correlated with H, scores, or metagraph position?
  KEY QUESTION: Does the cycle cover polynomial P_T(x) classify tournaments
                more finely than H alone?
""")
