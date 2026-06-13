#!/usr/bin/env python3
"""
At n=6, H=33 has THREE alpha decompositions: (12,2), (14,1), (16,0).
Compute M for one representative of each to check if tr(M^2) separates them.

This tests the conjecture that spectral moments encode the independence
polynomial beyond just H = I(Omega, 2).

Also test H=9 with (4,0) and (2,1) decompositions.
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

def transfer_matrix_n6(A):
    """Transfer matrix for n=6 tournament."""
    n = 6
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

n = 6
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)

# Find one tiling for each alpha profile at H=33 and H=9
print("Finding representative tilings...")

targets = {
    9: [(4, 0), (2, 1)],
    33: [(12, 2), (14, 1), (16, 0)]
}

found = {}
for bits in range(2**m):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    if H not in targets:
        continue
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)
    a1 = alphas.get(1, 0)
    a2 = alphas.get(2, 0)
    profile = (a1, a2)
    key = (H, profile)
    if key not in found and profile in targets[H]:
        found[key] = bits
        print(f"  Found H={H}, (a1,a2)={profile}: tiling {format(bits, f'0{m}b')}")

    if all((H, p) in found for p in targets[H]) and all(
        all((h, p) in found for p in targets[h]) for h in targets
    ):
        break

# Now compute transfer matrix for each
print()
print("=" * 70)
print("TRANSFER MATRIX SPECTRAL COMPARISON")
print("=" * 70)

for H in sorted(targets.keys()):
    print(f"\n--- H = {H} ---")
    for profile in targets[H]:
        key = (H, profile)
        if key not in found:
            print(f"  (a1,a2)={profile}: NOT FOUND")
            continue

        bits = found[key]
        A = tiling_to_tournament(bits, n)

        print(f"\n  (a1,a2)={profile}, tiling={format(bits, f'0{m}b')}")
        print(f"  Computing transfer matrix...", flush=True)
        M = transfer_matrix_n6(A)

        tr1 = int(np.trace(M))
        tr2 = int(np.trace(M @ M))
        tr3 = int(np.trace(M @ M @ M))
        det_M = int(round(np.linalg.det(M)))
        sigma = int(M.sum() - np.trace(M))
        evals = sorted(np.linalg.eigvalsh(M).tolist(), reverse=True)

        print(f"  tr(M) = {tr1}, Sigma = {sigma}, sum(M) = {tr1 + sigma}")
        print(f"  tr(M^2) = {tr2}")
        print(f"  tr(M^3) = {tr3}")
        print(f"  det(M) = {det_M}")
        print(f"  eigenvalues = {[round(e, 3) for e in evals]}")
        print(f"  M symmetric? {np.allclose(M, M.T)}")

        # At even n: tr(M) = 0, sum(M) = 2H
        print(f"  Even n check: tr(M)={tr1} (should=0), sum(M)={tr1+sigma} (should={2*H})")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
If tr(M^2) differs between alpha profiles at the same H,
then the spectral structure of M encodes MORE than just H.
This would mean the transfer matrix "knows about" the
independence polynomial structure, not just its value at x=2.
""")
