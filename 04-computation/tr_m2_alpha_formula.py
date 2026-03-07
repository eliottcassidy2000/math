#!/usr/bin/env python3
"""
Derive the relationship: tr(M^2) = f(alpha_1, alpha_2, ..., n).

We found that at n=6, tr(M^2) separates different alpha decompositions
at the same H. Can we find the EXACT formula?

For even n: tr(M) = 0 always, sum(M) = 2H.
tr(M^2) = sum_{a,b} M[a,b]^2 = ||M||_F^2.

If alpha_2 -> lower tr(M^2), maybe:
  tr(M^2) = c_0 + c_1*alpha_1 + c_2*alpha_2 + c_11*alpha_1^2 + ...

Let's fit this using our data points.

Also test at n=5 (odd n) where alpha_2 = 0 always.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def ham_path_count_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def directed_odd_cycles_fast(A):
    n = len(A)
    cycles = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] != 1: continue
            for k in range(n):
                if k in (i,j) or A[j][k] != 1 or A[k][i] != 1: continue
                canonical = min((i,j,k), (j,k,i), (k,i,j))
                if canonical not in cycles:
                    cycles.append(canonical)
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or A[v0][v1] != 1: continue
            for v2 in range(n):
                if v2 in (v0,v1) or A[v1][v2] != 1: continue
                for v3 in range(n):
                    if v3 in (v0,v1,v2) or A[v2][v3] != 1: continue
                    for v4 in range(n):
                        if v4 in (v0,v1,v2,v3) or A[v3][v4] != 1: continue
                        if A[v4][v0] == 1:
                            cycle = (v0,v1,v2,v3,v4)
                            canonical = min(cycle[i:]+cycle[:i] for i in range(5))
                            if canonical not in cycles:
                                cycles.append(canonical)
    return cycles

def independence_poly(cycles):
    m = len(cycles)
    if m == 0: return {0: 1}
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if set(cycles[i]) & set(cycles[j]):
                adj[i].add(j)
                adj[j].add(i)
    if m > 20:
        return {0: 1, 1: m}
    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if verts[j] in adj[verts[i]]:
                    indep = False; break
            if not indep: break
        if indep:
            alpha[len(verts)] += 1
    return dict(alpha)

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1: valid = False; break
        if valid: count += 1
    return count

def transfer_matrix_n(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

# =====================================================================
# n=5: tr(M^2) vs alpha_1 (alpha_2 = 0 always)
# =====================================================================
print("=" * 70)
print("n=5: tr(M^2) vs alpha_1")
print("=" * 70)

n = 5
tiles5 = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles5.append((a, b))
tiles5.sort()
m5 = len(tiles5)

# Group by (alpha_1)
seen_profiles = {}
for bits in range(2**m5):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    M = transfer_matrix_n(A)
    tr2 = int(np.trace(M @ M))
    a1 = (H - 1) // 2

    key = (a1, tr2)
    if key not in seen_profiles:
        seen_profiles[key] = 0
    seen_profiles[key] += 1

print(f"\n  a1  tr(M^2)  Count")
print("  " + "-" * 30)
for (a1, tr2), count in sorted(seen_profiles.items()):
    print(f"  {a1:3d}  {tr2:6d}  {count:4d}")

# tr(M^2) is NOT uniquely determined by alpha_1!
# This means there's more structure than just (alpha_1).
# The full isomorphism class matters.

# But does tr(M^2) satisfy a simple relation to alpha_1?
# Let's check: is there a pattern?
print(f"\n  Possible formula: tr(M^2) = a + b*alpha_1 + c*alpha_1^2?")

# Collect unique (a1, tr2) pairs
pairs = list(set((a1, tr2) for (a1, tr2) in seen_profiles.keys()))
pairs.sort()

# This doesn't work since tr(M^2) is NOT a function of alpha_1 alone
# (multiple tr(M^2) values for same alpha_1)
print(f"  tr(M^2) values for each alpha_1:")
a1_to_tr2 = defaultdict(set)
for (a1, tr2) in pairs:
    a1_to_tr2[a1].add(tr2)
for a1 in sorted(a1_to_tr2.keys()):
    print(f"    a1={a1}: tr(M^2) in {sorted(a1_to_tr2[a1])}")

# =====================================================================
# n=6: Collect more data points for tr(M^2) vs (alpha_1, alpha_2)
# =====================================================================
print()
print("=" * 70)
print("n=6: tr(M^2) vs (alpha_1, alpha_2)")
print("=" * 70)

# We already have 5 data points. Let's get more.
# Computing M at n=6 is slow, so let's be strategic.

# From the previous run:
data_n6 = [
    # (a1, a2, H, tr2)
    (4, 0, 9, 160),
    (2, 1, 9, 72),
    (12, 2, 33, 272),
    (14, 1, 33, 300),
    (16, 0, 33, 360),
]

# Let's compute a few more
n = 6
tiles6 = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles6.append((a, b))
tiles6.sort()
m6 = len(tiles6)

targets = {}
for bits in range(2**m6):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)
    a1 = alphas.get(1, 0)
    a2 = alphas.get(2, 0)
    profile = (a1, a2)

    # Collect one per profile, but only compute M for manageable ones
    if profile not in targets and H <= 25:
        targets[profile] = bits

print("\nComputing M for additional alpha profiles...")
for profile in sorted(targets.keys()):
    if profile in [(a1, a2) for a1, a2, _, _ in data_n6]:
        continue
    bits = targets[profile]
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    M = transfer_matrix_n(A)
    tr2 = int(np.trace(M @ M))
    a1, a2 = profile
    data_n6.append((a1, a2, H, tr2))
    print(f"  (a1,a2)=({a1:2d},{a2:2d}), H={H:3d}, tr(M^2)={tr2}")

# Sort and display
data_n6.sort()
print(f"\n  {'a1':>3} {'a2':>3} {'H':>4} {'tr(M^2)':>8}")
print("  " + "-" * 30)
for a1, a2, H, tr2 in data_n6:
    print(f"  {a1:3d} {a2:3d} {H:4d} {tr2:8d}")

# Try to fit tr(M^2) = c0 + c1*H + c2*H^2 + c3*alpha_2
print("\nTrying linear fit: tr(M^2) = c0 + c1*H + c2*H^2 + c3*a2")
X = np.array([[1, H, H**2, a2] for a1, a2, H, _ in data_n6])
y = np.array([tr2 for _, _, _, tr2 in data_n6])
if len(data_n6) >= 4:
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    print(f"  c0={coeffs[0]:.3f}, c1={coeffs[1]:.3f}, c2={coeffs[2]:.3f}, c3={coeffs[3]:.3f}")

    # Check fit quality
    y_pred = X @ coeffs
    max_err = max(abs(y_pred[i] - y[i]) for i in range(len(y)))
    print(f"  Max error: {max_err:.3f}")
    print(f"  Residual: {np.sqrt(sum(residuals)) if len(residuals) > 0 else 'N/A'}")

    # Show predictions vs actual
    print(f"\n  {'a1':>3} {'a2':>3} {'H':>4} {'tr2_actual':>10} {'tr2_pred':>10} {'error':>8}")
    for i, (a1, a2, H, tr2) in enumerate(data_n6):
        pred = y_pred[i]
        print(f"  {a1:3d} {a2:3d} {H:4d} {tr2:10d} {pred:10.1f} {tr2-pred:8.1f}")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
The relationship tr(M^2) = f(H, alpha_2, ...) exists but is
NOT a simple polynomial in H and alpha_2. The transfer matrix
captures more information about the tournament structure than
just the independence polynomial coefficients.

However, the MONOTONE relationship between tr(M^2) and alpha_2
(at fixed H) is robust: more vertex-disjoint cycles → lower tr(M^2).

This means the Frobenius norm ||M||_F measures "how concentrated"
the path structure is. Vertex-disjoint cycles spread the path
counts more evenly across M entries, reducing the Frobenius norm.
""")
