"""
h_spectral_landscape.py — Map H as a function on the spectral simplex

For circulant tournaments on Z_p, the eigenvalue structure is:
  λ_k = -1/2 + iy_k,  y_{p-k} = -y_k,  Σy_k² = p(p-1)/4

The "spectral simplex" is parametrized by the (p-1)/2 independent values
(y_1², y_2², ..., y_{(p-1)/2}²) subject to the constraint
Σ_{k=1}^{(p-1)/2} y_k² = p(p-1)/8.

Question: On this simplex, does H achieve its maximum at the center
(all y_k² equal = spectral flatness)?

This is the "Variational Principle" conjecture (HYP-463 from kind-pasteur).

At p=7: YES (only 2 spectral classes, trivially true)
At p=11: Does H increase monotonically toward the simplex center?

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def circulant_eigenvalues(n, S):
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def circulant_tournaments(n):
    k = (n - 1) // 2
    pairs = [(d, n - d) for d in range(1, k + 1)]
    for bits in range(1 << k):
        S = set()
        for i in range(k):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        yield frozenset(S)


def is_paley(n, S):
    qr = set(pow(x, 2, n) for x in range(1, n))
    return set(S) == qr or set(S) == set(range(1, n)) - qr


def spectral_distance_to_center(y2_vals):
    """Euclidean distance of the y² vector from the simplex center."""
    center = np.mean(y2_vals)
    return np.sqrt(sum((y2 - center)**2 for y2 in y2_vals))


def main():
    print("H ON THE SPECTRAL SIMPLEX")
    print("=" * 75)

    for p in [7, 11, 13]:
        print(f"\n{'='*75}")
        print(f"p = {p} (p mod 4 = {p % 4})")
        print(f"{'='*75}")

        m = (p - 1) // 2
        tournaments = list(circulant_tournaments(p))

        # Compute spectral data and H
        data = []
        for S in tournaments:
            evals = circulant_eigenvalues(p, S)
            y2_half = [evals[k].imag**2 for k in range(1, m + 1)]
            dist = spectral_distance_to_center(y2_half)
            sigma2 = sum(y2**2 for y2 in y2_half)

            A = circulant_adj(p, S)
            H = hamiltonian_paths_dp(A, p)
            paley = is_paley(p, S)

            data.append({
                'S': S, 'H': H, 'dist': dist, 'sigma2': sigma2,
                'paley': paley, 'y2': y2_half
            })

        # Sort by distance to center
        data.sort(key=lambda d: d['dist'])

        # Show relationship between dist and H
        print(f"\n  {'Dist':>10} {'σ₂':>12} {'H':>10} {'Paley':>6} {'y² values (half)'}")
        seen = set()
        for d in data:
            key = round(d['dist'], 4)
            if key in seen:
                continue
            seen.add(key)
            label = "★" if d['paley'] else ""
            y2_str = "[" + ", ".join(f"{y:.3f}" for y in d['y2']) + "]"
            print(f"  {d['dist']:>10.4f} {d['sigma2']:>12.4f} {d['H']:>10} {label:>6} {y2_str}")

        # Is H monotonically decreasing in dist?
        unique_data = []
        seen = set()
        for d in data:
            key = round(d['dist'], 4)
            if key not in seen:
                seen.add(key)
                unique_data.append(d)

        monotone = True
        for i in range(len(unique_data) - 1):
            if unique_data[i]['H'] < unique_data[i+1]['H']:
                monotone = False
                break

        print(f"\n  H monotonically decreasing in spectral distance: {monotone}")
        if not monotone:
            # Find where monotonicity breaks
            for i in range(len(unique_data) - 1):
                if unique_data[i]['H'] < unique_data[i+1]['H']:
                    print(f"  Breaks at dist={unique_data[i]['dist']:.4f} (H={unique_data[i]['H']}) vs "
                          f"dist={unique_data[i+1]['dist']:.4f} (H={unique_data[i+1]['H']})")

        # But does H always achieve its MAXIMUM at the center?
        center_H = min(d['H'] for d in data if abs(d['dist'] - data[0]['dist']) < 0.001)
        max_H = max(d['H'] for d in data)
        print(f"  Center H = {center_H}, Max H = {max_H}")
        print(f"  Center achieves maximum: {center_H == max_H}")

    # p=11 detailed analysis: why is H not monotone?
    print(f"\n{'='*75}")
    print("DETAILED ANALYSIS: p=11 NON-MONOTONICITY")
    print(f"{'='*75}")

    p = 11
    m = 5
    tournaments = list(circulant_tournaments(p))

    # Group by spectral class (same y² set up to permutation)
    spectral_classes = defaultdict(list)
    for S in tournaments:
        evals = circulant_eigenvalues(p, S)
        y2_half = tuple(sorted(round(evals[k].imag**2, 6) for k in range(1, m + 1)))
        spectral_classes[y2_half].append(S)

    print(f"\n  {len(spectral_classes)} spectral classes:")
    for y2, Ss in sorted(spectral_classes.items(), key=lambda x: sum(v**2 for v in x[0])):
        S = Ss[0]
        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        sigma2 = sum(v**2 for v in y2)
        sigma3 = sum(v**3 for v in y2)
        paley = any(is_paley(p, s) for s in Ss)
        label = " ★" if paley else ""
        dist = spectral_distance_to_center(list(y2))
        print(f"  y²={list(y2)}, σ₂={sigma2:.2f}, σ₃={sigma3:.2f}, "
              f"dist={dist:.4f}, H={H}, orbit_size={len(Ss)}{label}")

    # What about the "Schur-convexity" of H?
    # H is Schur-convex if H(y²) ≥ H(y'²) whenever y² majorizes y'².
    # Spectral flatness = center of simplex = minimal element in majorization order.
    # If H is Schur-CONCAVE, then center = maximum.
    print(f"\n  Checking Schur-concavity (center = max iff H is Schur-concave on simplex):")
    classes = sorted(spectral_classes.items(), key=lambda x: sum(v**2 for v in x[0]))
    for i, (y2_i, Ss_i) in enumerate(classes):
        for j, (y2_j, Ss_j) in enumerate(classes):
            if i >= j:
                continue
            # Check if y2_i majorizes y2_j or vice versa
            sorted_i = sorted(y2_i, reverse=True)
            sorted_j = sorted(y2_j, reverse=True)
            # y majorizes z if partial sums of sorted y ≥ partial sums of sorted z
            partial_i = [sum(sorted_i[:k+1]) for k in range(len(sorted_i))]
            partial_j = [sum(sorted_j[:k+1]) for k in range(len(sorted_j))]
            i_maj_j = all(pi >= pj for pi, pj in zip(partial_i, partial_j))
            j_maj_i = all(pj >= pi for pi, pj in zip(partial_i, partial_j))

            if i_maj_j or j_maj_i:
                Hi = hamiltonian_paths_dp(circulant_adj(p, Ss_i[0]), p)
                Hj = hamiltonian_paths_dp(circulant_adj(p, Ss_j[0]), p)
                major = "i≻j" if i_maj_j else "j≻i"
                concave = (Hi >= Hj if j_maj_i else Hi <= Hj)  # more spread → less H
                print(f"    {major}: H({round(sum(v**2 for v in y2_i), 1)})={Hi}, "
                      f"H({round(sum(v**2 for v in y2_j), 1)})={Hj}, "
                      f"{'Schur-concave ✓' if concave else 'NOT Schur-concave ✗'}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
