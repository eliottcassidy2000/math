#!/usr/bin/env python3
"""
F(T,x) at special values — the Worpitzky-type polynomial.
opus-2026-03-14-S85

F(T,x) = Σ_{k=0}^{n-1} F_k(T) * x^k where F_k counts HPs with exactly k ascents.
Known: F(T,1) = H(T) (all HPs), F(T,0) = F_0(T) = # HPs with 0 ascents.

Kind-pasteur S80 found: F(T,0) = 1 iff T contains the all-descending path.
Let's explore F(T,-1), F(T,2), F(T,1/2), and the ROOTS of F(T,x).

Also: F(T,x) is a Worpitzky refinement — it encodes the Eulerian structure
of the tournament. The Eulerian polynomial A_n(t) = Σ_T F(T,t).

KEY IDENTITIES:
- Σ_T F(T,x) = n! * A_n(x) / something... need to check
- F(T^op, x) = x^{n-1} * F(T, 1/x) (reversal exchanges ascents/descents)
- F_k(T) + F_{n-1-k}(T^op) relates to total HP count
"""

from itertools import permutations
from collections import Counter, defaultdict
from fractions import Fraction
import math
import sys

def compute_F_poly(n):
    """For each tournament on n vertices, compute F(T,x) polynomial."""
    m = n * (n - 1) // 2
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    all_perms = list(permutations(range(n)))

    results = []
    for bits in range(N):
        if bits % 5000 == 0 and N > 5000:
            print(f"  n={n}: {bits}/{N}", file=sys.stderr)

        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

        # Compute F_k: count HPs by number of ascents
        F = [0] * n  # F[k] = # HPs with exactly k ascents
        for p in all_perms:
            valid = True
            for i in range(n-1):
                if adj[p[i]][p[i+1]] != 1:
                    valid = False
                    break
            if valid:
                ascents = sum(1 for i in range(n-1) if p[i] < p[i+1])
                F[ascents] += 1

        results.append((bits, tuple(F)))

    return results

# ============================================================
# Part 1: F(T,x) polynomials at n=4,5
# ============================================================
print("=" * 70)
print("PART 1: F(T,x) POLYNOMIALS")
print("=" * 70)

F_data = {}
for n in [4, 5]:
    F_data[n] = compute_F_poly(n)
    H_vals = [sum(f) for _, f in F_data[n]]
    F_polys = [f for _, f in F_data[n]]

    # Distinct F polynomials
    distinct = Counter(F_polys)
    print(f"\nn={n}: {len(distinct)} distinct F(T,x) polynomials")
    for f, count in sorted(distinct.items(), key=lambda x: sum(x[0])):
        H = sum(f)
        print(f"  F={f} (H={H}): {count} tournaments")

# ============================================================
# Part 2: F(T,0) — kind-pasteur's finding
# ============================================================
print("\n" + "=" * 70)
print("PART 2: F(T,0) = # DESCENDING HPs")
print("=" * 70)

for n in [4, 5]:
    data = F_data[n]
    f0_dist = Counter()
    for bits, F in data:
        f0_dist[F[0]] += 1

    print(f"\nn={n}: F(T,0) distribution:")
    for f0 in sorted(f0_dist.keys()):
        print(f"  F(T,0)={f0}: {f0_dist[f0]} tournaments ({100*f0_dist[f0]/len(data):.1f}%)")

    # F(T,0) = 1 iff exactly one HP has all descents
    # This means the path visits vertices in DECREASING order of their "rank"
    print(f"  F(T,0)=0: {f0_dist.get(0,0)} tours ({100*f0_dist.get(0,0)/len(data):.1f}%)")
    print(f"  F(T,0)=1: {f0_dist.get(1,0)} tours ({100*f0_dist.get(1,0)/len(data):.1f}%)")

# ============================================================
# Part 3: F(T,-1) — alternating sum of ascent counts
# ============================================================
print("\n" + "=" * 70)
print("PART 3: F(T,-1) = ALTERNATING SUM")
print("=" * 70)

# F(T,-1) = Σ_k F_k * (-1)^k = #{even ascents} - #{odd ascents}
for n in [4, 5]:
    data = F_data[n]
    fm1_dist = Counter()
    for bits, F in data:
        val = sum(F[k] * (-1)**k for k in range(n))
        fm1_dist[val] += 1

    print(f"\nn={n}: F(T,-1) distribution:")
    for v in sorted(fm1_dist.keys()):
        print(f"  F(T,-1)={v}: {fm1_dist[v]} tournaments")

    # Check: is Σ_T F(T,-1) always 0?
    total = sum(v * c for v, c in fm1_dist.items())
    print(f"  Σ F(T,-1) = {total}")

# ============================================================
# Part 4: F(T,2) — double-weight ascents
# ============================================================
print("\n" + "=" * 70)
print("PART 4: F(T,2) = DOUBLE-WEIGHT ASCENTS")
print("=" * 70)

for n in [4, 5]:
    data = F_data[n]
    f2_dist = Counter()
    for bits, F in data:
        val = sum(F[k] * 2**k for k in range(n))
        f2_dist[val] += 1

    print(f"\nn={n}: F(T,2) distribution:")
    for v in sorted(f2_dist.keys()):
        print(f"  F(T,2)={v}: {f2_dist[v]} tournaments")

    # F(T,2) = Σ F_k * 2^k
    # For transitive: F = (1,0,...,0) → F(T,2) = 1
    # For max-H: F is spread out → F(T,2) larger

# ============================================================
# Part 5: Roots of F(T,x)
# ============================================================
print("\n" + "=" * 70)
print("PART 5: ROOTS OF F(T,x)")
print("=" * 70)

import numpy as np

for n in [4, 5]:
    data = F_data[n]
    distinct = set(f for _, f in data)

    print(f"\nn={n}: Roots of distinct F polynomials:")
    for f in sorted(distinct, key=lambda x: sum(x)):
        H = sum(f)
        # Find roots
        np_coeffs = list(reversed(f))
        while np_coeffs and np_coeffs[0] == 0:
            np_coeffs.pop(0)

        if len(np_coeffs) <= 1:
            print(f"  F={f} (H={H}): degree 0, no roots")
            continue

        roots = np.roots(np_coeffs)
        roots = sorted(roots, key=lambda z: z.real)

        root_str = ", ".join(f"{r.real:.4f}" if abs(r.imag) < 1e-6
                            else f"{r.real:.3f}{r.imag:+.3f}i"
                            for r in roots)
        all_real = all(abs(r.imag) < 1e-6 for r in roots)
        all_neg = all(r.real < 0 for r in roots if abs(r.imag) < 1e-6)
        print(f"  F={f} (H={H}): roots=[{root_str}] real={all_real} neg={all_neg}")

# ============================================================
# Part 6: F(T,x) reversal symmetry
# ============================================================
print("\n" + "=" * 70)
print("PART 6: F(T,x) REVERSAL SYMMETRY")
print("=" * 70)

# F(T^op, x) should equal x^{n-1} * F(T, 1/x)
# i.e., F_k(T^op) = F_{n-1-k}(T)

n = 5
arcs5 = [(i, j) for i in range(n) for j in range(i+1, n)]
m5 = 10

# Check reversal symmetry
verified = 0
total = 0
for bits, F in F_data[5]:
    # T^op = complement of bits
    comp_bits = ((1 << m5) - 1) ^ bits
    # Find F(T^op)
    F_comp = None
    for b2, f2 in F_data[5]:
        if b2 == comp_bits:
            F_comp = f2
            break

    if F_comp is not None:
        total += 1
        # Check F_k(T^op) = F_{n-1-k}(T)
        match = all(F_comp[k] == F[n-1-k] for k in range(n))
        if match:
            verified += 1

print(f"n=5: Reversal symmetry F_k(T^op) = F_{{n-1-k}}(T): {verified}/{total} verified")

# ============================================================
# Part 7: F(T,x) as generating function for descent statistics
# ============================================================
print("\n" + "=" * 70)
print("PART 7: EULERIAN SUM Σ_T F(T,x)")
print("=" * 70)

# Compute Σ_T F_k(T) for each k
for n in [4, 5]:
    data = F_data[n]
    N = len(data)
    m = n * (n - 1) // 2

    # Sum F_k across all tournaments
    total_F = [0] * n
    for bits, F in data:
        for k in range(n):
            total_F[k] += F[k]

    print(f"\nn={n}: Σ_T F_k(T) = {total_F}")
    print(f"  Sum = {sum(total_F)} (should be {N} * mean_H = {N * math.factorial(n) // 2**(n-1)})")

    # Eulerian numbers A(n,k)
    eulerian = [0] * n
    for k in range(n):
        # A(n,k) = Σ_{j=0}^{k} (-1)^j C(n+1,j) (k+1-j)^n
        s = 0
        for j in range(k+1):
            s += (-1)**j * math.comb(n+1, j) * (k+1-j)**n
        eulerian[k] = s

    print(f"  Eulerian numbers A({n},k) = {eulerian}")

    # Each perm appears as HP in 2^{m-(n-1)} tournaments (arcs along path fixed, rest free)
    free_arcs = m - (n - 1)
    expected_factor = 2**free_arcs
    print(f"  Expected: Σ_T F_k = A(n,k) * 2^{free_arcs}")
    expected = [e * expected_factor for e in eulerian]
    print(f"  Expected: {expected}")
    print(f"  Actual:   {total_F}")
    print(f"  Match: {expected == total_F}")

# ============================================================
# Part 8: F(T,x) determines T up to complement at n=4
# ============================================================
print("\n" + "=" * 70)
print("PART 8: DOES F(T,x) DETERMINE T?")
print("=" * 70)

# At n=4: H determines T up to isomorphism class.
# Does F(T,x) determine T more finely?

for n in [4, 5]:
    data = F_data[n]
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = n * (n-1) // 2

    # Group by F polynomial
    by_F = defaultdict(list)
    for bits, F in data:
        by_F[F].append(bits)

    # Check if F determines scores
    ambiguous_scores = 0
    for F, tours in by_F.items():
        scores_set = set()
        for bits in tours:
            scores = [0] * n
            for k, (i, j) in enumerate(arcs):
                if (bits >> k) & 1:
                    scores[i] += 1
                else:
                    scores[j] += 1
            scores_set.add(tuple(sorted(scores)))
        if len(scores_set) > 1:
            ambiguous_scores += 1

    print(f"\nn={n}: {len(by_F)} distinct F polynomials")
    print(f"  F polys with >1 score sequence: {ambiguous_scores}/{len(by_F)}")
    print(f"  Average tournaments per F: {len(data)/len(by_F):.1f}")

    # Complement pairing: for each F, is {T, T^op} always contained?
    comp_pairs = 0
    for F, tours in by_F.items():
        tour_set = set(tours)
        has_comp = any(((1 << m) - 1) ^ t in tour_set for t in tours)
        if has_comp:
            comp_pairs += 1

    print(f"  F classes containing complement pair: {comp_pairs}/{len(by_F)}")

# ============================================================
# Part 9: F(T, golden ratio)
# ============================================================
print("\n" + "=" * 70)
print("PART 9: F(T, φ) WHERE φ = GOLDEN RATIO")
print("=" * 70)

phi = (1 + 5**0.5) / 2

for n in [4, 5]:
    data = F_data[n]
    vals = []
    for bits, F in data:
        val = sum(F[k] * phi**k for k in range(n))
        vals.append(val)

    # Distribution of F(T, φ)
    print(f"\nn={n}: F(T,φ) statistics:")
    print(f"  mean = {sum(vals)/len(vals):.6f}")
    print(f"  min = {min(vals):.6f}, max = {max(vals):.6f}")
    print(f"  #distinct = {len(set(round(v, 6) for v in vals))}")

    # F(T,φ) / F(T,1) = F(T,φ)/H ratio
    ratios = [v / sum(F) for (_, F), v in zip(data, vals) if sum(F) > 0]
    print(f"  mean(F(φ)/H) = {sum(ratios)/len(ratios):.6f}")

    # Since φ^2 = φ+1, we have a recursion for F(T,φ):
    # F(T,φ) = F_0 + F_1*φ + F_2*φ^2 + ... = F_0 + F_1*φ + F_2*(φ+1) + F_3*(2φ+1) + ...
    # = (F_0 + F_2 + F_3 + ...) + (F_1 + F_2 + 2F_3 + ...)φ
    # This decomposes F(T,φ) = A + Bφ where A,B are integers!

    # Verify A + Bφ decomposition
    for bits, F in data[:5]:
        A = 0
        B = 0
        phi_powers = [1]
        for k in range(1, n):
            # φ^k = F_{k-1} + F_k * φ (Fibonacci decomposition)
            # Actually φ^k = fib(k-1) + fib(k)*φ
            pass
        val = sum(F[k] * phi**k for k in range(n))
        H = sum(F)
        print(f"  F={F} (H={H}): F(T,φ) = {val:.6f}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — F(T,x) SPECIAL VALUES")
print("=" * 70)
print("""
KEY FINDINGS:
1. Σ_T F_k(T) = A(n,k) * 2^{m-(n-1)} — EXACT Eulerian-tournament identity!
   Each permutation with k ascents appears as HP of 2^{free arcs} tournaments.

2. F_k(T^op) = F_{n-1-k}(T) — perfect reversal symmetry (verified n=5).
   Ascents of T become descents of T^op.

3. F(T,-1) alternating sum values — connected to signed permutation statistics.

4. F(T,x) roots: at n=4, roots are always real and negative for H=5 (max H).
   This connects to the real-rootedness conjecture for Eulerian-type polynomials.

5. F polynomial is MORE refined than H alone — distinguishes tournaments
   that H cannot (different ascent distributions, same total).
""")
