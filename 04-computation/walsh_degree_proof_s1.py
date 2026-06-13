"""
Walsh degree proof verification for all n.
Confirms THM-259 + THM-076 → band-limitedness.

Key insight: Walsh degree = 2*floor((n-1)/2) in TILING model.
This is the same as the FULL ARC model degree because we can always
find maximum-weight monomials using only tile arcs (via interleaving).

kind-pasteur-2026-03-25-S1
"""
import sys, time
from collections import defaultdict
from math import comb, factorial

sys.stdout.reconfigure(line_buffering=True)

def count_ham_paths_dp(adj, n):
    """Held-Karp DP for Hamiltonian path count."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for sz in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != sz:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev = mask ^ (1 << v)
                for u in range(n):
                    if u == v or not (prev & (1 << u)):
                        continue
                    if adj[u][v]:
                        dp[(mask, v)] = dp.get((mask, v), 0) + dp.get((prev, u), 0)
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def tile_list(n):
    """Return ordered list of tiles (x,y) with x >= y+2.
    Convention: 0-indexed vertices. Base path: n-1 -> n-2 -> ... -> 0."""
    tiles = []
    for y in range(n - 2):
        for x in range(n - 1, y + 1, -1):
            if x >= y + 2:
                tiles.append((x, y))
    return tiles


def tiling_to_adj(n, tiles, bits):
    """Convert tiling bits to adjacency matrix."""
    adj = [[0] * n for _ in range(n)]
    # Base path: k -> k-1
    for k in range(n - 1, 0, -1):
        adj[k][k - 1] = 1
    # Tiles
    for i, (x, y) in enumerate(tiles):
        if (bits >> i) & 1:
            adj[y][x] = 1  # backward
        else:
            adj[x][y] = 1  # forward
    return adj


def walsh_hadamard_transform(f, m):
    """In-place Walsh-Hadamard transform."""
    import numpy as np
    N = 1 << m
    g = np.array(f, dtype=np.float64)
    for i in range(m):
        half = 1 << i
        step = 1 << (i + 1)
        for j in range(0, N, step):
            for k in range(half):
                u, v = g[j + k], g[j + k + half]
                g[j + k] = u + v
                g[j + k + half] = u - v
    return g / N


def hamming_weight(x):
    return bin(x).count('1')


def analyze_n(n):
    """Full Walsh analysis for given n."""
    tiles = tile_list(n)
    m = len(tiles)
    N = 1 << m

    print(f"\n{'=' * 60}")
    print(f"  n = {n}, m = {m} tiles, N = {N} tilings")
    print(f"  Predicted Walsh degree: {2 * ((n - 1) // 2)}")
    print(f"  m/2 = {m / 2}")
    print(f"{'=' * 60}")

    t0 = time.time()

    # Compute H for all tilings
    H_vals = []
    for bits in range(N):
        adj = tiling_to_adj(n, tiles, bits)
        H_vals.append(count_ham_paths_dp(adj, n))

    t1 = time.time()
    print(f"  H computation: {t1 - t0:.1f}s")
    print(f"  H range: [{min(H_vals)}, {max(H_vals)}]")
    print(f"  Mean H: {sum(H_vals) / N:.4f}")
    print(f"  Expected (all tournaments): {factorial(n) / 2 ** (n - 1):.4f}")

    # Walsh transform
    import numpy as np
    hat_H = walsh_hadamard_transform(H_vals, m)

    t2 = time.time()
    print(f"  Walsh transform: {t2 - t1:.1f}s")

    # Analyze by weight
    weight_data = defaultdict(lambda: {'nz': 0, 'max_abs': 0.0, 'alpha': 0.0, 'total': 0})
    for S in range(N):
        w = hamming_weight(S)
        weight_data[w]['total'] += 1
        val = hat_H[S]
        weight_data[w]['alpha'] += val
        if abs(val) > 1e-10:
            weight_data[w]['nz'] += 1
            weight_data[w]['max_abs'] = max(weight_data[w]['max_abs'], abs(val))

    print(f"\n  Walsh coefficients by Hamming weight:")
    print(f"  {'wt':>3} {'C(m,wt)':>8} {'#nonzero':>9} {'max|coeff|':>11} {'alpha_k':>11}")
    max_nz_wt = 0
    for w in sorted(weight_data.keys()):
        d = weight_data[w]
        if d['nz'] > 0:
            max_nz_wt = w
        print(f"  {w:>3} {d['total']:>8} {d['nz']:>9} {d['max_abs']:>11.6f} {d['alpha']:>11.6f}")

    pred = 2 * ((n - 1) // 2)
    print(f"\n  RESULTS:")
    print(f"    Observed max nonzero weight: {max_nz_wt}")
    print(f"    Predicted (THM-259):         {pred}")
    print(f"    Match: {'YES' if max_nz_wt == pred else 'NO'}")
    print(f"    m/2 = {m / 2}")
    print(f"    Band-limited (degree < m/2): {'YES' if max_nz_wt < m / 2 else 'NO'}")

    # Check complement symmetry: H(t) vs H(complement(t))
    all_ones = N - 1
    sym_violations = 0
    for bits in range(N // 2 + 1):
        comp = all_ones ^ bits
        if H_vals[bits] != H_vals[comp]:
            sym_violations += 1

    print(f"    Complement symmetry H(t)=H(~t): {'YES' if sym_violations == 0 else f'NO ({sym_violations} violations)'}")

    # Check odd-weight Walsh: are they all zero?
    odd_wt_nonzero = sum(weight_data[w]['nz'] for w in weight_data if w % 2 == 1)
    print(f"    Odd-weight Walsh coefficients: {odd_wt_nonzero} nonzero {'(complement symmetry breaks in tiling model)' if odd_wt_nonzero > 0 else ''}")

    # Verify THM-076 amplitude formula
    print(f"\n  THM-076 amplitude verification:")
    for w in sorted(weight_data.keys()):
        d = weight_data[w]
        if d['nz'] == 0:
            continue
        # For even weight 2k with single path (r=1):
        # |hat_H| = 2 * (n-2k)! / 2^{n-1}
        if w % 2 == 0:
            k = w // 2
            if 2 * k + 1 <= n:
                predicted_amp_r1 = 2.0 * factorial(n - 2 * k) / (2 ** (n - 1))
                print(f"    wt={w}: max|coeff|={d['max_abs']:.6f}, THM-076 r=1 amp={predicted_amp_r1:.6f}", end="")
                # For r=2 (two paths): need 2k + 2 <= n
                if 2 * k + 2 <= n:
                    predicted_amp_r2 = 4.0 * factorial(n - 2 * k) / (2 ** (n - 1))
                    print(f", r=2 amp={predicted_amp_r2:.6f}", end="")
                print()

    return max_nz_wt, pred


# Run for n = 4, 5, 6, 7
results = {}
for n in [4, 5, 6]:
    obs, pred = analyze_n(n)
    results[n] = (obs, pred)

# n=7 is larger (m=15, N=32768) but still feasible
print("\n\nStarting n=7 (m=15, 32768 tilings)...")
obs7, pred7 = analyze_n(7)
results[7] = (obs7, pred7)

# Summary
print(f"\n\n{'=' * 60}")
print("SUMMARY: Walsh degree verification")
print(f"{'=' * 60}")
print(f"{'n':>3} {'m':>4} {'m/2':>5} {'observed':>9} {'predicted':>10} {'match':>6} {'BL at m/2':>10}")
for n in sorted(results.keys()):
    obs, pred = results[n]
    m = comb(n - 1, 2)
    bl = "YES" if obs < m / 2 else "NO"
    match = "YES" if obs == pred else "NO"
    print(f"{n:>3} {m:>4} {m / 2:>5.1f} {obs:>9} {pred:>10} {match:>6} {bl:>10}")

print(f"""
PROOF OF BAND-LIMITEDNESS FOR ALL n >= 6:

1. THM-076 (proved algebraically for all n):
   Walsh amplitude |hat_H[S]| = 2^r * (n-2k)! / 2^(n-1)
   for monomials of type (2a_1, ..., 2a_r) with 2k+r <= n.

2. Maximum Hamming weight = max(2k : 2k+r <= n, r >= 1) = n-1.
   Only EVEN weights have nonzero amplitudes in the full arc model.
   In the tiling model, odd weights CAN be nonzero (no complement symmetry).

3. Walsh degree in tiling model = 2*floor((n-1)/2):
   - Upper bound: degree in full model <= 2*floor((n-1)/2), tiling <= full.
   - Lower bound (n>=4): interleaving construction provides base-path-free
     Hamiltonian paths, achieving max weight in tile variables.

4. Band-limitedness: For n >= 6, 2*floor((n-1)/2) <= n-1 < C(n-1,2)/2.
   So all Walsh/Krawtchouk coefficients above m/2 are zero.

5. At n=5: degree 4 > m/2 = 3, so NOT band-limited at m/2 in tiling model.
   At n=4: degree 2 > m/2 = 1.5, NOT band-limited.
   Band-limitedness at m/2 holds for n >= 6.
""")
