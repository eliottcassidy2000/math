#!/usr/bin/env python3
"""
H² - det(I+2A) ANALYSIS
opus-2026-03-13-S67k

OBSERVATION: H² - det(I+2A) appears always divisible by 8.
Also H² - det ≥ 0 always.

This script investigates:
1. Is (H² - det)/8 always a non-negative integer?
2. What does (H² - det)/8 count?
3. How does it relate to the channel decomposition?
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict

def adj_matrix(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def det_I2A(A, n):
    M = np.eye(n, dtype=float) + 2 * A.astype(float)
    return int(round(np.linalg.det(M)))

def fast_hash(A, n):
    scores = list(A.sum(axis=1))
    ns = []
    for i in range(n):
        o = sorted(scores[j] for j in range(n) if A[i][j])
        ins = sorted(scores[j] for j in range(n) if A[j][i])
        ns.append((scores[i], tuple(o), tuple(ins)))
    return tuple(sorted(ns))

def find_odd_cycles(A, n, max_len=None):
    if max_len is None:
        max_len = n if n % 2 == 1 else n - 1
    cycles = set()
    for length in range(3, max_len + 1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    mi = perm.index(min(perm))
                    cycles.add(tuple(perm[mi:] + perm[:mi]))
    return cycles

def conflict_graph_and_indpoly(cycles):
    cl = list(cycles)
    nc = len(cl)
    vs = [set(c) for c in cl]
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vs[i] & vs[j]:
                adj[i][j] = adj[j][i] = 1
    # Compute independence polynomial coefficients
    coeffs = [0] * (nc + 1)
    coeffs[0] = 1
    for size in range(1, nc + 1):
        cnt = 0
        for sub in combinations(range(nc), size):
            ok = True
            for a in range(len(sub)):
                for b in range(a+1, len(sub)):
                    if adj[sub[a]][sub[b]]:
                        ok = False
                        break
                if not ok: break
            if ok: cnt += 1
        coeffs[size] = cnt
    return coeffs

print("=" * 72)
print("H² - det(I+2A) ANALYSIS")
print("=" * 72)

print("\n--- Part 1: Is (H² - det)/8 always a non-negative integer? ---\n")

for n in range(3, 7):
    print(f"n={n}:")
    m = n*(n-1)//2
    hash_groups = defaultdict(list)
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    classes = []
    for h, A in hash_groups.items():
        H = count_hp(A, n)
        d = det_I2A(A, n)
        cycles = find_odd_cycles(A, n)
        coeffs = conflict_graph_and_indpoly(cycles) if cycles else [1]
        alpha1 = coeffs[1] if len(coeffs) > 1 else 0
        alpha2 = coeffs[2] if len(coeffs) > 2 else 0
        c3 = sum(1 for c in cycles if len(c) == 3)
        c5 = sum(1 for c in cycles if len(c) == 5)
        classes.append((H, d, alpha1, alpha2, c3, c5))

    classes.sort()

    all_div8 = True
    all_nonneg = True
    for H, d, a1, a2, c3, c5 in classes:
        diff = H*H - d
        div8 = diff // 8 if diff % 8 == 0 else None
        if diff % 8 != 0:
            all_div8 = False
        if diff < 0:
            all_nonneg = False

        div8_str = str(div8) if div8 is not None else "NOT INT"
        print(f"  H={H:3d}  det={d:6d}  H²={H*H:6d}  H²-det={diff:6d}  "
              f"(H²-det)/8={div8_str:>7s}  "
              f"α₁={a1:2d} α₂={a2:2d}  c3={c3:2d} c5={c5:2d}")

    print(f"  All divisible by 8: {all_div8}")
    print(f"  All non-negative: {all_nonneg}")
    print()

print("\n--- Part 2: What does (H²-det)/8 depend on? ---\n")

print("Testing: (H²-det)/8 vs α₂, α₁·(α₁-1)/2, c3·(c3-1)/2, etc.")
for n in range(3, 7):
    print(f"\nn={n}:")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        d = det_I2A(A, n)
        diff8 = (H*H - d) // 8
        cycles = find_odd_cycles(A, n)
        coeffs = conflict_graph_and_indpoly(cycles) if cycles else [1]
        a1 = coeffs[1] if len(coeffs) > 1 else 0
        a2 = coeffs[2] if len(coeffs) > 2 else 0

        # Test formulas
        test1 = a1 * (a1 - 1) // 2  # C(α₁, 2)
        test2 = a1 * a1  # α₁²
        test3 = a1 * (a1 + 1) // 2  # C(α₁+1, 2)
        test4 = H * (H - 1) // 2  # C(H, 2)

        # Since H = 1+2α₁+4α₂, H² = 1+4α₁+4α₁²+16α₂²+4α₂+16α₁α₂
        # H² - det = ?
        # If det = (1+2α₁)² - 8·something...
        # H² = (1+2α₁+4α₂)² = 1+4α₁+4α₁²+8α₂+16α₁α₂+16α₂²
        # If det ≈ (1+2α₁)² = 1+4α₁+4α₁²
        # Then H²-det ≈ 8α₂+16α₁α₂+16α₂² = 8α₂(1+2α₁+2α₂)

        formula_test = a2 * (1 + 2*a1 + 2*a2) if a2 > 0 else 0

        print(f"  H={H:3d}  (H²-det)/8={diff8:5d}  "
              f"α₁={a1:2d} α₂={a2:2d}  "
              f"C(α₁,2)={test1:4d}  α₁²={test2:4d}  "
              f"α₂(1+2α₁+2α₂)={formula_test:4d}  "
              f"match: {'✓' if diff8 == formula_test else '✗'}")

print("\n--- Part 3: Detailed formula search ---\n")

# Collect all data points
data = []
for n in range(3, 7):
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in hash_groups.items():
        H = count_hp(A, n)
        d = det_I2A(A, n)
        diff8 = (H*H - d) // 8
        cycles = find_odd_cycles(A, n)
        coeffs = conflict_graph_and_indpoly(cycles) if cycles else [1]
        a1 = coeffs[1] if len(coeffs) > 1 else 0
        a2 = coeffs[2] if len(coeffs) > 2 else 0
        a3 = coeffs[3] if len(coeffs) > 3 else 0
        data.append((n, H, d, diff8, a1, a2, a3))

# Try: (H²-det)/8 = f(α₁, α₂)
# H = 1 + 2α₁ + 4α₂ + 8α₃
# H² = 1 + 4α₁ + 4α₁² + 16α₂² + 64α₃² + 4α₂·8 + 4α₁·4α₂·2 + ...
# Let me just compute H² - (√det)² directly
# √det is some function of the tournament; det = (√det)²

# Better approach: compute H² - det algebraically
# H = 1 + 2α₁ + 4α₂ + 8α₃
# H² = 1 + 4α₁ + 4α₂·2 + (2α₁)² + 2·2α₁·4α₂ + (4α₂)² + ... cross terms
# Actually: H = 1 + 2S where S = α₁ + 2α₂ + 4α₃
# H² = 1 + 4S + 4S² = 1 + 4(α₁+2α₂+4α₃) + 4(α₁+2α₂+4α₃)²
# H² = 1 + 4α₁ + 8α₂ + 16α₃ + 4(α₁² + 4α₂² + 16α₃² + 4α₁α₂ + 8α₂α₃ + 16α₁α₃)

# So if det = 1 + 4α₁ + 4α₁², then:
# H²-det = 8α₂ + 16α₃ + 4(4α₂² + 4α₁α₂ + ...) = 8(α₂ + ...) + 4(...)
# This is getting messy. Let me just check numerically.

print("Checking: (H²-det)/8 vs various formulas:")
print(f"{'n':>2s} {'H':>4s} {'det':>6s} {'(H²-d)/8':>9s} {'α₁':>3s} {'α₂':>3s} {'α₃':>3s} "
      f"{'α₁α₂':>5s} {'C(α₁,2)':>7s} {'Δ':>5s}")

for n, H, d, diff8, a1, a2, a3 in data:
    a1a2 = a1 * a2
    ca1_2 = a1 * (a1-1) // 2

    # Try: diff8 = α₁·α₂ + something?
    # Or: diff8 = (H²-det)/8 = α₂·H/2 + ...
    # Let's just see what diff8 - a1*a2 looks like
    residual = diff8 - a1a2

    print(f"{n:2d} {H:4d} {d:6d} {diff8:9d} {a1:3d} {a2:3d} {a3:3d} "
          f"{a1a2:5d} {ca1_2:7d} {residual:5d}")

print("\n--- Part 4: Factor H²-det completely ---\n")

# Since H² - det ≥ 0 and H, √det are both odd:
# H - √det is even, H + √det is even
# H² - det = (H-√det)(H+√det)
# Both factors even → product divisible by 4
# Need to show divisible by 8

print("H²-det = (H-√det)·(H+√det):")
for n, H, d, diff8, a1, a2, a3 in data:
    sd = int(round(abs(d)**0.5))
    if sd*sd != abs(d):
        continue
    minus = H - sd
    plus = H + sd
    # Both even since H odd, sd odd → H±sd even
    print(f"  H={H:3d}  √det={sd:3d}  H-√det={minus:4d}  H+√det={plus:4d}  "
          f"product/4={minus*plus//4:5d}  mod4(H-√det)={minus%4}  mod4(H+√det)={plus%4}  "
          f"(H²-det)/8={diff8}")

print("""
ANALYSIS:
  H is always odd. √det is always odd.
  So H ± √det are both even.
  H² - det = (H-√det)(H+√det) = even × even.

  More precisely: H ≡ 1 (mod 2), √det ≡ 1 (mod 2)
  H - √det ≡ 0 (mod 2)
  H + √det ≡ 0 (mod 2)

  For divisibility by 8: need one of the factors ≡ 0 (mod 4).
  Since H ≡ √det (mod 2), they're both odd.
  H = 2k+1, √det = 2j+1, so H-√det = 2(k-j), H+√det = 2(k+j+1)
  (k-j)(k+j+1) = k²+k-j²-j = ... one of these is always even!
  Because (k-j) and (k+j+1) have different parity:
    k-j is even ⟺ k+j is even ⟺ k+j+1 is odd
    k-j is odd ⟺ k+j is odd ⟺ k+j+1 is even
  So EXACTLY ONE of (k-j), (k+j+1) is even.
  Therefore (k-j)(k+j+1) is even, and (H-√det)(H+√det) = 4(k-j)(k+j+1) is divisible by 8.

  PROVED: H² - det(I+2A) ≡ 0 (mod 8) for all tournaments!

  This is because H and √det are both odd, and the factorization
  (H-√det)(H+√det) = 2·even · 2·something = 4·even = 8k.
""")

print("\n--- Part 5: Summary ---\n")

print("PROVED IDENTITIES:")
print("  1. det(I+2A) is always a perfect square (HYP-788)")
print("  2. √det(I+2A) is always odd")
print("  3. H² - det(I+2A) ≡ 0 (mod 8)  [proved above, elementary]")
print("  4. H² ≥ det(I+2A)  [observed, not yet proved in general]")
print()
print("OPEN:")
print("  5. What does (H²-det)/8 count combinatorially?")
print("  6. Is (H²-det)/8 always a polynomial in α₁, α₂, α₃?")
print("  7. Does H² = det when T has no 3-cycles (only transitive)?")
print()

# Check Q7
print("Check Q7: does H=√det when c3=0?")
for n, H, d, diff8, a1, a2, a3 in data:
    cycles_info = (a1, a2, a3)
    sd = int(round(abs(d)**0.5))
    if a1 == 0:  # no cycles
        print(f"  n={n}, H={H}, √det={sd}, H=√det: {H==sd}")
    if diff8 == 0:
        print(f"  n={n}, H={H}, √det={sd}, H²=det: TRUE, α₁={a1} α₂={a2}")
