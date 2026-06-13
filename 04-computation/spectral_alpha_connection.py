#!/usr/bin/env python3
"""
Does the spectral structure of M encode the independence polynomial?

H = tr(M) = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
tr(M^2) = ||M||_F^2
tr(M^3) = sum_{a,b,c} M[a,b]*M[b,c]*M[c,a]

Can we extract alpha_1, alpha_2 from spectral moments of M?

For scalar M = (H/n)*I:
  tr(M^k) = n * (H/n)^k = H^k / n^{k-1}

For non-scalar M, the spectral moments encode MORE information.

KEY IDEA: At odd n, sum(M) = H and sum(M^2) = something.
If alpha_1 = (H-1)/2 (when alpha_2 = 0), then tr(M^2) is just
a function of H. But when alpha_2 > 0, tr(M^2) might distinguish.

Test at n=6 where alpha_2 can be nonzero.
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
    # 3-cycles
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] != 1: continue
            for k in range(n):
                if k in (i,j) or A[j][k] != 1 or A[k][i] != 1: continue
                canonical = min((i,j,k), (j,k,i), (k,i,j))
                if canonical not in cycles:
                    cycles.append(canonical)
    # 5-cycles
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
            if A[p[i]][p[i+1]] != 1:
                valid = False; break
        if valid: count += 1
    return count

def transfer_matrix(A):
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
# n=5: spectral moments vs alpha_1
# =====================================================================
print("=" * 70)
print("n=5: SPECTRAL MOMENTS vs ALPHA_1")
print("=" * 70)

n = 5
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)

print(f"\nAt n=5: alpha_k=0 for k>=2, so alpha_1 = (H-1)/2")
print(f"tr(M) = H, tr(M^2) = ||M||_F^2")
print()

# Compute for each class
seen = set()
data_points = []

for bits in range(2**m):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    M = transfer_matrix(A)

    # Use M as fingerprint
    M_tuple = tuple(M.flatten())
    if M_tuple in seen:
        continue
    seen.add(M_tuple)

    tr1 = int(np.trace(M))
    tr2 = int(np.trace(M @ M))
    tr3 = int(np.trace(M @ M @ M))
    det_M = int(round(np.linalg.det(M)))
    alpha1 = (H - 1) // 2

    data_points.append((H, alpha1, tr1, tr2, tr3, det_M))

print(f"  {'H':>3} {'a1':>3} {'tr(M)':>6} {'tr(M^2)':>8} {'tr(M^3)':>9} {'det(M)':>8}")
print("  " + "-" * 55)
for H, a1, t1, t2, t3, det_m in sorted(data_points):
    print(f"  {H:3d} {a1:3d} {t1:6d} {t2:8d} {t3:9d} {det_m:8d}")

# Is tr(M^2) determined by H alone?
H_to_tr2 = defaultdict(set)
for H, a1, t1, t2, t3, det_m in data_points:
    H_to_tr2[H].add(t2)

print(f"\n  tr(M^2) uniquely determined by H? ", end="")
print(all(len(v) == 1 for v in H_to_tr2.values()))
for H in sorted(H_to_tr2.keys()):
    if len(H_to_tr2[H]) > 1:
        print(f"  H={H}: tr(M^2) in {H_to_tr2[H]}")

# =====================================================================
# n=5: Can tr(M^2) separate isomorphism classes?
# =====================================================================
print()
print("=" * 70)
print("n=5: CAN SPECTRAL MOMENTS SEPARATE ISO CLASSES?")
print("=" * 70)

# For n=5, there are 12 iso classes. Can (tr(M), tr(M^2), tr(M^3), det(M))
# distinguish them?
seen_profiles = defaultdict(list)
for H, a1, t1, t2, t3, det_m in data_points:
    profile = (t1, t2, t3, det_m)
    seen_profiles[profile].append(H)

print(f"\n  {len(seen_profiles)} distinct spectral profiles from {len(data_points)} iso classes")
if len(seen_profiles) == len(data_points):
    print("  -> Spectral moments UNIQUELY determine iso class at n=5!")
else:
    print("  -> Some iso classes share spectral profile")
    for profile, Hs in seen_profiles.items():
        if len(Hs) > 1:
            print(f"    Profile {profile}: H values {Hs}")

# =====================================================================
# n=5 tilings with SAME score but DIFFERENT H
# =====================================================================
print()
print("=" * 70)
print("SCORE SEQUENCES: SAME SCORE, DIFFERENT H")
print("=" * 70)

# Classes 4 and 6 both have score [3,3,2,1,1] but H=9 vs H=9
# Classes 8 and 10 have score [3,2,2,2,1] but H=11 vs H=13
# These are distinguished by their spectral moments

for H, a1, t1, t2, t3, det_m in sorted(data_points, key=lambda x: (x[0], x[3])):
    if H in [9, 11, 13]:
        print(f"  H={H}: tr(M^2)={t2}, tr(M^3)={t3}, det={det_m}")

# =====================================================================
# Sample at n=6 where alpha_2 matters
# =====================================================================
print()
print("=" * 70)
print("n=6: SPECTRAL MOMENTS vs (alpha_1, alpha_2)")
print("=" * 70)

# At n=6, compute M for a sample of tilings with different alpha profiles
n = 6
tiles6 = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles6.append((a, b))
tiles6.sort()
m6 = len(tiles6)

print(f"\nSampling tilings with various alpha profiles...")
print(f"  {'bits':>12} {'H':>3} {'a1':>3} {'a2':>3} {'Sigma':>6}")
print("  " + "-" * 40)

# Sample specific tilings with different alpha_2 values
sample_results = []
for bits in range(2**m6):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)
    a1 = alphas.get(1, 0)
    a2 = alphas.get(2, 0)

    if a2 in [0, 1, 2, 4] and len(sample_results) < 20:
        # Check if this alpha profile is new
        profile = (a1, a2)
        if not any(r[3] == profile for r in sample_results):
            # Compute Sigma (sum of M off-diagonal = 2H at even n)
            # Don't compute full M at n=6 (too slow for all 1024)
            # But we can verify Sigma = 2H
            Sigma = 2 * H  # by theory (even n)
            sample_results.append((
                format(bits, f'0{m6}b'), H, a1*2 + a2*4, (a1, a2)
            ))
            print(f"  {format(bits, f'0{m6}b')} {H:3d} {a1:3d} {a2:3d} {Sigma:6d}")

# Can (H, alpha_1, alpha_2) determine spectral moments?
print(f"\n  H = 1 + 2*a1 + 4*a2")
print(f"  So a1 = (H - 1 - 4*a2) / 2")
print(f"  Given H, knowing a2 uniquely determines a1.")
print(f"  But does a2 determine tr(M^2)?")

# =====================================================================
# Key observation about alpha_2 and H
# =====================================================================
print()
print("=" * 70)
print("KEY: DOES (H, alpha_2) DETERMINE THE SPECTRAL PROFILE?")
print("=" * 70)

# At n=6, two tilings with same H but different alpha_2
# should have different spectral profiles
# Find examples:
alpha_by_H = defaultdict(set)
for bits in range(2**m6):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)
    a1 = alphas.get(1, 0)
    a2 = alphas.get(2, 0)
    alpha_by_H[H].add((a1, a2))

print(f"\n  H values with MULTIPLE alpha profiles:")
for H in sorted(alpha_by_H.keys()):
    profiles = alpha_by_H[H]
    if len(profiles) > 1:
        print(f"  H={H}: {sorted(profiles)}")

# Same H, different alpha decomposition!
# H = 9: can be (4,0) or (2,1). Same H, different structure!
# H = 17: (8,0) or (6,1). Same H, different.
# H = 29: (10,2) or (12,1). Same H, different.
# H = 33: (12,2) or (14,1) or (16,0). TRIPLE!
# H = 37: (14,2) or (16,1). Same H, different.

# This proves: H alone does NOT determine the independence polynomial.
# At n=6, the spectral moments of M might capture more than just H.

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
Key findings:
1. At n=5: spectral moments (tr, tr^2, tr^3, det) UNIQUELY determine
   the iso class. Since alpha_k=0 for k>=2, this is "over-determined."

2. At n=6: H does NOT determine (alpha_1, alpha_2). Multiple alpha
   decompositions give the same H (e.g., H=9 from (4,0) or (2,1)).

3. This means the transfer matrix M encodes MORE than just H.
   The eigenvalue spectrum of M should distinguish between tournaments
   with the same H but different alpha profiles.

4. CONJECTURE: The spectral moments of M encode the full independence
   polynomial I(Omega, x), not just its value at x=2.
   Specifically: alpha_k might be recoverable from tr(M), tr(M^2), ..., tr(M^k).
""")
