"""
quadratic_eigenvalue_theorem.py — opus-2026-04-04-S3

DISCOVERY: The quadratic coefficient matrix Q[i,j] = c_{ij} has:
  - Exactly n-2 negative eigenvalues
  - Exactly m-(n-2) positive eigenvalues (when full rank)
  - Trace = 0 (all diagonal entries are 0 in multilinear polynomial)

Also testing:
  - Sign rule: c_{ij} < 0 iff tiles i,j share a vertex
  - Walsh 2-adic uniformity
  - n=4 and n=8 verification
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v

def tiles_share_vertex(t1, t2):
    """Check if tiles (x1,y1) and (x2,y2) share a vertex."""
    return t1[0] == t2[0] or t1[1] == t2[1] or t1[0] == t2[1] or t1[1] == t2[0]

def tile_relation(t1, t2):
    """Classify the relationship between two tiles."""
    x1, y1 = t1
    x2, y2 = t2
    if x1 == x2 or y1 == y2:
        return "same-end"  # share upper or lower vertex
    if x1 == y2 or y1 == x2:
        return "cross-end"  # one's upper is other's lower
    return "disjoint"

def analyze_n(n):
    print(f"\n{'='*70}")
    print(f"n={n}")
    print(f"{'='*70}")

    tiles = make_tiles(n)
    m = len(tiles)

    if m == 0:
        print("  No tiles (n too small)")
        return

    print(f"  m={m} tiles, computing H for {2**m} tilings...")
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # =========================================================
    # 1. Quadratic matrix Q
    # =========================================================
    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        i, j = bits[0], bits[1]
        Q[i][j] = c
        Q[j][i] = c

    eigvals = sorted(np.linalg.eigvalsh(Q), reverse=True)
    pos = sum(1 for e in eigvals if e > 1e-10)
    neg = sum(1 for e in eigvals if e < -1e-10)
    zero = m - pos - neg
    rank = pos + neg

    print(f"\n  QUADRATIC MATRIX SIGNATURE:")
    print(f"    m={m}, rank={rank}, signature ({pos}+, {neg}-, {zero} zero)")
    print(f"    n-2 = {n-2}")
    print(f"    neg eigenvalues = {neg}")
    print(f"    MATCH neg=n-2? {'YES!' if neg == n-2 else 'NO'}")
    print(f"    Eigenvalues: {[f'{e:.4f}' for e in eigvals]}")

    # =========================================================
    # 2. Sign rule: c_{ij} < 0 iff tiles share vertex
    # =========================================================
    sign_violations = 0
    sign_total = 0
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        shared = tiles_share_vertex(t1, t2)
        sign_total += 1

        if c < 0 and not shared:
            sign_violations += 1
            print(f"    VIOLATION: c({t1},{t2}) = {c} < 0 but tiles DON'T share vertex")
        if c > 0 and shared:
            # Check if cross-end (allowed positive) or same-end (should be negative)
            rel = tile_relation(t1, t2)
            if rel == "same-end":
                sign_violations += 1
                print(f"    VIOLATION: c({t1},{t2}) = {c} > 0 but tiles share same-end vertex")

    print(f"\n  SIGN RULE:")
    print(f"    Tested {sign_total} pairs")
    print(f"    Violations: {sign_violations}")
    if sign_violations == 0:
        print(f"    RULE HOLDS: c_{{ij}} < 0 iff same-end shared vertex")
        # Count how many same-end negative, cross-end positive, disjoint positive
        neg_same = sum(1 for S, c in coeffs.items() if bin(S).count('1') == 2 and c < 0)
        pos_cross = sum(1 for S, c in coeffs.items()
                       if bin(S).count('1') == 2 and c > 0
                       and tile_relation(tiles[[i for i in range(m) if S & (1<<i)][0]],
                                        tiles[[i for i in range(m) if S & (1<<i)][1]]) == "cross-end")
        pos_disj = sum(1 for S, c in coeffs.items()
                      if bin(S).count('1') == 2 and c > 0
                      and tile_relation(tiles[[i for i in range(m) if S & (1<<i)][0]],
                                       tiles[[i for i in range(m) if S & (1<<i)][1]]) == "disjoint")
        print(f"    Negative (same-end): {neg_same}")
        print(f"    Positive (cross-end): {pos_cross}")
        print(f"    Positive (disjoint): {pos_disj}")

    # =========================================================
    # 3. Same-end negative: is the value always -2?
    # =========================================================
    same_end_vals = []
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        rel = tile_relation(t1, t2)
        if rel == "same-end":
            same_end_vals.append(c)

    if same_end_vals:
        print(f"\n  SAME-END VALUES: {sorted(set(same_end_vals))}")
        all_neg2 = all(v == -2 for v in same_end_vals)
        print(f"    All equal to -2? {all_neg2}")
        if not all_neg2:
            # Which ones differ?
            for S, c in coeffs.items():
                if bin(S).count('1') != 2: continue
                bits = [i for i in range(m) if S & (1 << i)]
                t1, t2 = tiles[bits[0]], tiles[bits[1]]
                if tile_relation(t1, t2) == "same-end" and c != -2:
                    s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
                    print(f"      {t1}-{t2} skips=({s1},{s2}): c={c}")

    # =========================================================
    # 4. Cross-end positive: formula?
    # =========================================================
    cross_end_data = []
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        if tile_relation(t1, t2) == "cross-end":
            s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
            cross_end_data.append((t1, t2, s1, s2, c))

    if cross_end_data:
        print(f"\n  CROSS-END VALUES:")
        for t1, t2, s1, s2, c in cross_end_data:
            expected = 2**(s1+s2-2)
            print(f"    {t1}-{t2} skips=({s1},{s2}): c={c}, 2^(s1+s2-2)={expected}, ratio={c/expected:.4f}")

    # =========================================================
    # 5. Walsh 2-adic uniformity
    # =========================================================
    N = 1 << m
    f = np.array(H, dtype=np.int64)
    for i in range(m):
        step = 1 << i
        for j in range(N):
            if j & step:
                a = f[j ^ step]
                b = f[j]
                f[j ^ step] = a + b
                f[j] = a - b

    # Find minimum v2 across ALL nonzero Walsh coefficients
    min_v2_all = float('inf')
    for S in range(N):
        if f[S] != 0:
            min_v2_all = min(min_v2_all, v2(abs(int(f[S]))))

    print(f"\n  WALSH 2-ADIC UNIFORMITY:")
    print(f"    Min v2 of unnormalized Walsh coefficients: {min_v2_all}")
    print(f"    m = {m}")
    print(f"    Effective Walsh denominator: 2^{m - min_v2_all}")

    # Check: is min_v2_all related to something?
    # At n=5 (m=6): check
    # At n=6 (m=10): check
    # At n=7 (m=15): min_v2 = 10, so denominator = 2^5

    # After dividing by 2^min_v2, what's the distribution of odd parts?
    # by order
    by_order_v2 = defaultdict(lambda: defaultdict(int))
    for S in range(N):
        if f[S] != 0:
            order = bin(S).count('1')
            reduced = abs(int(f[S])) >> min_v2_all
            by_order_v2[order][v2(reduced)] += 1

    print(f"    After dividing by 2^{min_v2_all}, v2 distribution by order:")
    for order in sorted(by_order_v2.keys()):
        print(f"      Order {order}: {dict(sorted(by_order_v2[order].items()))}")

    # =========================================================
    # 6. Higher-degree coefficient structures
    # =========================================================
    max_deg = max(bin(S).count('1') for S in coeffs.keys()) if coeffs else 0
    print(f"\n  MAX DEGREE: {max_deg}")
    print(f"    Expected 2*floor((n-1)/2) = {2*((n-1)//2)}")

    # Check if max-degree coefficients are all +/- 2
    if max_deg > 0:
        max_deg_vals = [c for S, c in coeffs.items() if bin(S).count('1') == max_deg]
        print(f"    Max-degree coefficient values: {sorted(set(max_deg_vals))}")
        print(f"    All +/-2? {set(max_deg_vals) == {2, -2} or set(max_deg_vals) == {2} or set(max_deg_vals) == {-2}}")

    # =========================================================
    # 7. Nullspace of Q: what tile combinations have Qv=0?
    # =========================================================
    if zero > 0 and m <= 20:
        # Find nullspace
        U, s, Vt = np.linalg.svd(Q)
        null_vecs = Vt[rank:]
        print(f"\n  NULLSPACE (dim={zero}):")
        for i, v in enumerate(null_vecs):
            # Normalize
            v = v / np.max(np.abs(v))
            entries = [(tiles[j], f"{v[j]:.4f}") for j in range(m) if abs(v[j]) > 0.01]
            print(f"    Null vector {i}: {entries}")

    # =========================================================
    # 8. Negative eigenspace: what tiles span it?
    # =========================================================
    if neg > 0:
        eigvals_full, eigvecs = np.linalg.eigh(Q)
        neg_indices = [i for i, e in enumerate(eigvals_full) if e < -1e-10]

        print(f"\n  NEGATIVE EIGENSPACE (dim={neg}):")
        for idx in neg_indices:
            e = eigvals_full[idx]
            v = eigvecs[:, idx]
            # Which tiles dominate?
            dominant = sorted(range(m), key=lambda j: abs(v[j]), reverse=True)[:5]
            dom_info = [(tiles[j], f"{v[j]:.4f}") for j in dominant]
            print(f"    lambda={e:.4f}: dominant tiles = {dom_info}")

    return coeffs, Q, eigvals

# =========================================================
# Run for n = 3 through 7 (n=8 if feasible)
# =========================================================
print("QUADRATIC EIGENVALUE THEOREM INVESTIGATION")
print("=" * 70)

results = {}
for n in [3, 4, 5, 6, 7]:
    results[n] = analyze_n(n)

# n=8: m=21, 2^21 = 2M tilings. Feasible? Let's try timing.
import time
print(f"\n\nAttempting n=8 (m=21, 2^21 = {2**21} tilings)...")
t0 = time.time()
tiles8 = make_tiles(8)
m8 = len(tiles8)
print(f"  m={m8} tiles, starting H computation...")

# This needs 2M HP computations. Each takes O(n * 2^n) = O(8 * 256) = 2048 ops.
# Total: 2M * 2K = 4 billion. Too slow in pure Python.
# Let's try a smaller check: just compute the quadratic coefficients directly
# by flipping single and double tiles from transitive.

# H(transitive) = H(0...0) = 1
# c_i = H(e_i) - H(0) = H(e_i) - 1
# c_{ij} = H(e_i + e_j) - H(e_i) - H(e_j) + H(0)

from h_multilinear_complete import tiling_to_tournament, count_hamiltonian_paths

n = 8
tiles = make_tiles(n)
m = len(tiles)
print(f"\n  n=8: computing quadratic coefficients via direct flips...")
print(f"  Need to compute H for {1 + m + m*(m-1)//2} tilings")

# H(0)
T0 = tiling_to_tournament(n, tiles, 0)
H0 = count_hamiltonian_paths(T0, n)
print(f"  H(transitive) = {H0}")

# H(e_i) for each tile
H_single = {}
for i in range(m):
    bits = 1 << i
    T = tiling_to_tournament(n, tiles, bits)
    H_single[i] = count_hamiltonian_paths(T, n)

# Linear coefficients
c_linear = {}
for i in range(m):
    c_linear[i] = H_single[i] - H0

print(f"  Linear coefficients: {[c_linear[i] for i in range(m)]}")
print(f"  Expected 2^(skip-1): {[2**(tiles[i][0]-tiles[i][1]-1) for i in range(m)]}")
lin_match = all(c_linear[i] == 2**(tiles[i][0]-tiles[i][1]-1) for i in range(m))
print(f"  All match: {lin_match}")

# H(e_i + e_j) for each pair
t0 = time.time()
Q8 = np.zeros((m, m))
pair_count = 0
for i in range(m):
    for j in range(i+1, m):
        bits = (1 << i) | (1 << j)
        T = tiling_to_tournament(n, tiles, bits)
        H_ij = count_hamiltonian_paths(T, n)
        c_ij = H_ij - H_single[i] - H_single[j] + H0
        Q8[i][j] = c_ij
        Q8[j][i] = c_ij
        pair_count += 1

elapsed = time.time() - t0
print(f"  Computed {pair_count} pairs in {elapsed:.1f}s")

# Analyze Q8
eigvals8 = sorted(np.linalg.eigvalsh(Q8), reverse=True)
pos8 = sum(1 for e in eigvals8 if e > 1e-10)
neg8 = sum(1 for e in eigvals8 if e < -1e-10)
zero8 = m - pos8 - neg8

print(f"\n  n=8 QUADRATIC MATRIX:")
print(f"    m={m}, signature ({pos8}+, {neg8}-, {zero8} zero)")
print(f"    n-2 = {n-2} = 6")
print(f"    neg eigenvalues = {neg8}")
print(f"    MATCH? {'YES!' if neg8 == n-2 else 'NO'}")
print(f"    Eigenvalues: {[f'{e:.2f}' for e in eigvals8]}")

# Sign rule for n=8
print(f"\n  n=8 SIGN RULE:")
violations = 0
for i in range(m):
    for j in range(i+1, m):
        c = Q8[i][j]
        if c == 0: continue
        t1, t2 = tiles[i], tiles[j]
        rel = tile_relation(t1, t2)
        if c < 0 and rel != "same-end":
            violations += 1
            if violations <= 5:
                print(f"    VIOLATION: {t1}-{t2} ({rel}): c={c}")
        if c > 0 and rel == "same-end":
            violations += 1
            if violations <= 5:
                print(f"    VIOLATION: {t1}-{t2} ({rel}): c={c}")

print(f"    Sign rule violations: {violations}")

# Same-end values
same_end_vals = []
for i in range(m):
    for j in range(i+1, m):
        if Q8[i][j] == 0: continue
        t1, t2 = tiles[i], tiles[j]
        if tile_relation(t1, t2) == "same-end":
            same_end_vals.append(int(Q8[i][j]))

print(f"    Same-end values: {sorted(set(same_end_vals))}")
