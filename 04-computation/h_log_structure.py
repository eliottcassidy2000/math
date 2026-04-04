"""
h_log_structure.py — opus-2026-04-04-S3

The LOGARITHMIC structure of H(t).

Key discovery: log(H(t)) has FULL degree m (all tiles), while H has degree
only 2*floor((n-1)/2). The log is DENSER than H.

This means: H = exp(dense polynomial of degree m).
The sparsification from degree m to degree 2*floor((n-1)/2) happens through
the exponential — the OCF cycle packing structure IS the exponential.

The cluster expansion: log I(G, x) = Σ_{connected induced H ⊆ G} ...
For the independence polynomial: there's the Ursell function expansion.

SPECIFIC QUESTION: What are the log-linear coefficients?
log(H) has linear terms = log(c_0 + c_k) - log(c_0) = log(1 + c_k) ≈ c_k for small c_k.
But c_k can be large (up to 2^(n-2)), so log is significantly different from H.

DEEPER: log(H) = Σ_k λ_k t_k + higher order
where λ_k = log(H(e_k) / H(0)) = log((1 + c_k) / 1) = log(1 + 2^(skip-1))

For skip s: λ_k = log(1 + 2^(s-1))
  s=2: log(3) ≈ 1.099
  s=3: log(5) ≈ 1.609
  s=4: log(9) ≈ 2.197
  s=5: log(17) ≈ 2.833
  General: log(1 + 2^(s-1)) = (s-1)*log(2) + log(1 + 2^(1-s))
  ≈ (s-1)*log(2) for large s

So the log-linear coefficients are approximately (s-1)*log(2)!
This means: log H is approximately LINEAR in the skip values.

And the LOG-QUADRATIC coefficients measure the INTERACTION between tiles
in the log scale — the "mutual information" between tile flips.
"""
import numpy as np
from math import log, comb
from collections import defaultdict
from fractions import Fraction

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=int)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def compute_all_H(n):
    tiles = make_tiles(n)
    m = len(tiles)
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)
    return H, tiles, m

def mobius(vals, m):
    coeffs = {}
    for S in range(1 << m):
        val = 0.0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * vals[T]
            if T == 0: break
            T = (T - 1) & S
        if abs(val) > 1e-12:
            coeffs[S] = val
    return coeffs

def log_analysis(n):
    H, tiles, m = compute_all_H(n)
    logH = np.log(H.astype(float))

    print(f"\n{'='*60}")
    print(f"LOG STRUCTURE ANALYSIS n={n}, m={m}")
    print(f"{'='*60}")

    # Log-linear coefficients
    log_coeffs = mobius(logH, m)
    h_coeffs = mobius(H.astype(float), m)

    print(f"\nLINEAR COEFFICIENTS (H vs log H):")
    print(f"{'Tile':<10} {'skip':>5} {'c_k(H)':>10} {'λ_k(logH)':>12} {'log(1+c_k)':>12} {'match':>6}")
    for i in range(m):
        s = tiles[i][0] - tiles[i][1]
        c_h = int(h_coeffs.get(1 << i, 0))
        c_log = log_coeffs.get(1 << i, 0)
        expected = log(1 + c_h)  # log(H(e_i)/H(0))
        match = abs(c_log - expected) < 1e-8
        print(f"{str(tiles[i]):<10} {s:>5} {c_h:>10} {c_log:>12.6f} {expected:>12.6f} {'✓' if match else '✗':>6}")

    # Verify: λ_k = log(1 + 2^(s-1)) for all tiles
    all_match = True
    for i in range(m):
        s = tiles[i][0] - tiles[i][1]
        expected = log(1 + 2**(s-1))
        actual = log_coeffs.get(1 << i, 0)
        if abs(actual - expected) > 1e-8:
            all_match = False
    print(f"\nAll log-linear = log(1 + 2^(s-1)): {all_match}")

    # Log-quadratic: interaction strength
    print(f"\nQUADRATIC COEFFICIENTS (log H) — the 'interaction energy':")
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c_log = log_coeffs.get(S, 0)
            if abs(c_log) > 1e-10:
                t1, t2 = tiles[i], tiles[j]
                s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
                shared = set()
                if t1[0]==t2[0]: shared.add('upper')
                if t1[1]==t2[1]: shared.add('lower')
                if t1[0]==t2[1]: shared.add('cross')
                if t1[1]==t2[0]: shared.add('cross')
                disjoint = not ({t1[0],t1[1]} & {t2[0],t2[1]})
                rel = 'disjoint' if disjoint else str(shared) if shared else 'other'

                # Compare with H-quadratic
                c_h = h_coeffs.get(S, 0)

                # The log-quadratic has a nice interpretation:
                # λ_{ij} = log(H(ij)/H(0)) - log(H(i)/H(0)) - log(H(j)/H(0))
                # = log(H(ij)·H(0) / (H(i)·H(j)))
                # This is the LOG-ODDS RATIO (or "interaction energy" in Ising model)

                H00 = int(H[0])
                H10 = int(H[1 << i])
                H01 = int(H[1 << j])
                H11 = int(H[(1 << i) | (1 << j)])
                odds_ratio = (H11 * H00) / (H10 * H01) if H10 * H01 > 0 else 0
                log_or = log(odds_ratio) if odds_ratio > 0 else 0

                if n <= 5 or abs(c_log) > 0.5:
                    print(f"  {str(t1):<8}-{str(t2):<8} s=({s1},{s2}) {rel:<15}: "
                          f"λ₂={c_log:>8.4f}, c₂(H)={c_h:>6.0f}, "
                          f"odds={odds_ratio:.4f}, log(odds)={log_or:.4f}")

    # THE BIG PICTURE: log(H) as "energy landscape"
    print(f"\n{'='*60}")
    print(f"THE ENERGY LANDSCAPE")
    print(f"{'='*60}")

    # In the Ising model analogy:
    # log(H(t)) = Σ_i h_i t_i + Σ_{i<j} J_{ij} t_i t_j + higher order
    # h_i = "external field" (log-linear coefficient)
    # J_{ij} = "coupling" (log-quadratic coefficient)

    # The "mean field" approximation: E[H] ≈ exp(Σ h_i/2 + Σ J_{ij}/4)

    fields = [log_coeffs.get(1 << i, 0) for i in range(m)]
    couplings = []
    for i in range(m):
        for j in range(i+1, m):
            c = log_coeffs.get((1 << i) | (1 << j), 0)
            couplings.append(c)

    mean_field = np.exp(sum(f/2 for f in fields) + sum(c/4 for c in couplings))
    actual_mean = np.mean(H)
    print(f"  Mean field prediction: {mean_field:.2f}")
    print(f"  Actual mean H: {actual_mean:.2f}")
    print(f"  Ratio: {mean_field/actual_mean:.4f}")

    # The field strengths by skip
    skip_fields = defaultdict(list)
    for i in range(m):
        s = tiles[i][0] - tiles[i][1]
        skip_fields[s].append(fields[i])
    print(f"\n  Field strength by skip:")
    for s in sorted(skip_fields.keys()):
        vals = skip_fields[s]
        print(f"    skip {s}: {vals[0]:.4f} = log(1+2^{s-1}) = log({1+2**(s-1)})")

    # Coupling statistics
    if couplings:
        neg_c = [c for c in couplings if c < 0]
        pos_c = [c for c in couplings if c > 0]
        print(f"\n  Coupling statistics:")
        print(f"    Negative (repulsive): {len(neg_c)} couplings, "
              f"mean={np.mean(neg_c):.4f}" if neg_c else "    No negative couplings")
        print(f"    Positive (attractive): {len(pos_c)} couplings, "
              f"mean={np.mean(pos_c):.4f}" if pos_c else "    No positive couplings")
        print(f"    Zero: {len(couplings)-len(neg_c)-len(pos_c)}")

if __name__ == "__main__":
    log_analysis(4)
    log_analysis(5)
    log_analysis(6)
