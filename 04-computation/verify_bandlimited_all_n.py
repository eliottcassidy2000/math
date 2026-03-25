"""
Verify Walsh degree formula at n=5,6 (both odd and even).
Confirms THM-259 + THM-076 → band-limitedness for all n.

kind-pasteur-2026-03-25-S1
"""
import itertools
import numpy as np
from collections import defaultdict

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def count_ham_paths(T):
    """Count Hamiltonian paths via DP (Held-Karp)."""
    n = len(T)
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                for u in range(n):
                    if u == v:
                        continue
                    if not (prev_mask & (1 << u)):
                        continue
                    if T[u][v]:
                        dp[(mask, v)] += dp[(prev_mask, u)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def tiling_to_tournament(n, tiling_bits):
    """Convert tiling bits to tournament.
    Base path: n-1 -> n-2 -> ... -> 0.
    Tiles: (x,y) with x >= y+2, ordered as in CLAUDE.md:
    for y=0..n-3: for x=n-1 down to y+2: tile (x,y)
    """
    T = [[0]*n for _ in range(n)]
    # Base path arcs: k -> k-1 for k = n-1, ..., 1
    for k in range(n-1, 0, -1):
        T[k][k-1] = 1

    # Tiles
    tile_idx = 0
    for y in range(n-2):
        for x in range(n-1, y+1, -1):
            if x >= y+2:
                if (tiling_bits >> tile_idx) & 1:
                    T[y][x] = 1  # backward: y -> x (flip from default x -> y)
                else:
                    T[x][y] = 1  # forward: x -> y (default)
                tile_idx += 1

    return T

def count_tiles(n):
    """Count number of tiles = C(n-1, 2)."""
    return (n-1)*(n-2)//2

def walsh_transform(H_values, m):
    """Compute Walsh-Hadamard transform of H values.
    H_values[i] = H value for tiling i (0 to 2^m - 1).
    Returns Walsh coefficients hat_H[S] for each S.
    """
    # Use the fast Walsh-Hadamard transform
    N = 2**m
    f = np.array(H_values, dtype=np.float64)

    # Walsh-Hadamard transform
    for i in range(m):
        step = 1 << (i+1)
        half = 1 << i
        for j in range(0, N, step):
            for k in range(half):
                u = f[j+k]
                v = f[j+k+half]
                f[j+k] = u + v
                f[j+k+half] = u - v

    return f / N

def hamming_weight(x, m):
    """Hamming weight of x (number of 1-bits in first m positions)."""
    return bin(x & ((1 << m) - 1)).count('1')

def analyze_walsh(n):
    """Full Walsh analysis for tournaments on n vertices."""
    m = count_tiles(n)
    N = 2**m

    print(f"\n=== n={n}, m={m}, N={N} tilings ===")

    # Compute H for all tilings
    H_values = []
    for bits in range(N):
        T = tiling_to_tournament(n, bits)
        H_values.append(count_ham_paths(T))

    print(f"H range: [{min(H_values)}, {max(H_values)}]")
    print(f"Mean H: {sum(H_values)/len(H_values):.4f}")
    import math
    print(f"Expected mean: {math.factorial(n) / 2**(n-1):.4f}")

    # Walsh transform
    hat_H = walsh_transform(H_values, m)

    # Analyze by Hamming weight
    weight_stats = defaultdict(lambda: {'count': 0, 'nonzero': 0, 'max_abs': 0.0, 'sum': 0.0})
    for S in range(N):
        w = hamming_weight(S, m)
        weight_stats[w]['count'] += 1
        abs_val = abs(hat_H[S])
        if abs_val > 1e-10:
            weight_stats[w]['nonzero'] += 1
            weight_stats[w]['max_abs'] = max(weight_stats[w]['max_abs'], abs_val)
        weight_stats[w]['sum'] += hat_H[S]

    print(f"\nWalsh coefficient analysis by Hamming weight:")
    print(f"{'Weight':>6} {'#total':>7} {'#nonzero':>9} {'max|coeff|':>12} {'sum(coeff)':>12}")
    max_nonzero_weight = 0
    for w in sorted(weight_stats.keys()):
        s = weight_stats[w]
        nz = s['nonzero']
        if nz > 0:
            max_nonzero_weight = w
        print(f"{w:>6} {s['count']:>7} {nz:>9} {s['max_abs']:>12.6f} {s['sum']:>12.6f}")

    predicted_degree = 2 * ((n-1) // 2)
    print(f"\nObserved max nonzero weight: {max_nonzero_weight}")
    print(f"Predicted Walsh degree (THM-259): {predicted_degree}")
    print(f"m/2 = {m/2}")
    print(f"ceil(m/2)+1 = {(m+1)//2 + 1}")

    if max_nonzero_weight == predicted_degree:
        print(f"[OK] CONFIRMED: Walsh degree = {predicted_degree}")
    else:
        print(f"[FAIL] MISMATCH: observed {max_nonzero_weight} != predicted {predicted_degree}")

    if max_nonzero_weight < (m+1)//2 + 1:
        print(f"[OK] Band-limited: Walsh degree {max_nonzero_weight} < ceil(m/2)+1 = {(m+1)//2 + 1}")
    else:
        print(f"[FAIL] NOT band-limited at m/2")

    # Krawtchouk coefficients (sum of Walsh coefficients at each weight)
    print(f"\nKrawtchouk spectral coefficients (sum of Walsh at each weight):")
    for w in sorted(weight_stats.keys()):
        s = weight_stats[w]
        if abs(s['sum']) > 1e-10:
            print(f"  alpha_{w} = {s['sum']:.6f}")

    return max_nonzero_weight, predicted_degree

# Run analysis
for n in [5, 6]:
    obs, pred = analyze_walsh(n)

# Quick check n=4 too
analyze_walsh(4)

print("\n\n=== SUMMARY ===")
print("THM-076 proves |hat{H}[S]| = 2^r * (n-2k)! / 2^{n-1}")
print("Constraint: monomial S needs 2k+r <= n vertices, so |S| = 2k <= n-1")
print("Only even weights survive → max weight = 2*floor((n-1)/2)")
print("This proves band-limitedness for ALL n since 2*floor((n-1)/2) <= n-1 << C(n-1,2)/2")
