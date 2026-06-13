#!/usr/bin/env python3
"""
paley_energy_minimizer.py -- kind-pasteur-2026-03-13-S60

KEY INSIGHT: Paley tournament MINIMIZES additive energy E(S) across all
orientations. Since disj3 = (p/4)*E + C(p) (THM-156), this means Paley
MINIMIZES disj3. But we know from S18h (BIBD analysis) that Paley
MAXIMIZES H by minimizing alpha_2 (=disj3)... wait, BIBD MINIMIZES alpha_2!

Let me check: does min E correspond to min or max disj3?
disj3 = (p/4)*E + C. Since p/4 > 0, min E => min disj3.
And BIBD (Paley) minimizes disj3 (from S18h analysis).
So: Paley minimizes BOTH E and disj3.

This means:
  min disj3 <=> min E(S) <=> max Fourier uniformity of S
  <=> QR set (Paley) has the most uniform Fourier spectrum

This is a NUMBER THEORY result: QR sets minimize additive energy among
all half-sets of Z_p* with the tournament orientation constraint.

Questions:
1. Does QR minimize E among ALL (p-1)/2-element subsets of Z_p*?
   (Known: QR has smallest E among sum-free-like sets? Check literature)
2. How does H relate to E? Is max H <=> min E?
3. Connection to Weil bound for character sums
4. Does interval MAXIMIZE E? (consistent with energy spectrum data)
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s)%p] = 1
    return A


def additive_energy(S, p):
    S_set = set(S)
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1
    return energy


def compute_H(A, p):
    """Compute H(T) = number of Hamiltonian paths for small p."""
    if p > 11:
        return None
    n = p
    # Held-Karp for Hamiltonian paths
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_c3_disj3(A, p):
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))
    c3 = len(c3_sets)
    disj3 = 0
    for i in range(c3):
        for j in range(i+1, c3):
            if not (c3_sets[i] & c3_sets[j]):
                disj3 += 1
    return c3, disj3


def analyze_H_vs_energy(p):
    """Full analysis: E, disj3, H for all orientations."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"H vs ENERGY ANALYSIS at p={p}, m={m}")
    print(f"{'='*70}")

    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        E = additive_energy(S, p)
        c3, disj3 = count_c3_disj3(A, p)
        H = compute_H(A, p)
        data.append((bits, S, E, c3, disj3, H))

    # Sort by energy
    data.sort(key=lambda x: x[2])

    # Show all unique (E, disj3, H) triples
    seen = {}
    for bits, S, E, c3, disj3, H in data:
        key = (E, disj3, H)
        if key not in seen:
            seen[key] = []
        seen[key].append(bits)

    print(f"\n  Unique (E, disj3, H) triples (sorted by E):")
    print(f"  {'E':>6} {'disj3':>6} {'H':>8} {'count':>6}")
    for key in sorted(seen.keys()):
        E, disj3, H = key
        print(f"  {E:6d} {disj3:6d} {H:8d} {len(seen[key]):6d}")

    # Check: min E => max H?
    E_vals = [d[2] for d in data]
    H_vals = [d[5] for d in data if d[5] is not None]
    if H_vals:
        r = np.corrcoef([d[2] for d in data if d[5] is not None],
                         [d[5] for d in data if d[5] is not None])[0, 1]
        print(f"\n  Corr(E, H) = {r:.6f}")

    # Identify Paley and Interval
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    for bits, S, E, c3, disj3, H in data:
        S_set = set(S)
        if S_set == QR:
            print(f"\n  PALEY: bits={bits:0{m}b}, S={S}, E={E}, disj3={disj3}, H={H}")
        if set(S) == set(range(1, m+1)):
            print(f"  INTERVAL: bits={'1'*m}, S={S}, E={E}, disj3={disj3}, H={H}")

    # Check: does E have the same ordering as -H (opposite)?
    # Build mapping: for each unique E value, what are the possible H values?
    E_to_H = defaultdict(set)
    for bits, S, E, c3, disj3, H in data:
        if H is not None:
            E_to_H[E].add(H)

    print(f"\n  E -> H mapping:")
    for E in sorted(E_to_H.keys()):
        H_set = sorted(E_to_H[E])
        print(f"    E={E}: H in {H_set}")

    # Is E a function of H? (i.e., does each H value correspond to unique E?)
    H_to_E = defaultdict(set)
    for bits, S, E, c3, disj3, H in data:
        if H is not None:
            H_to_E[H].add(E)

    E_determines_H = all(len(v) == 1 for v in E_to_H.values())
    H_determines_E = all(len(v) == 1 for v in H_to_E.values())
    print(f"\n  E determines H: {E_determines_H}")
    print(f"  H determines E: {H_determines_E}")

    return data


def energy_vs_fourier_uniformity(p):
    """Show that Paley minimizes E because QR has the most uniform Fourier spectrum.

    For S = QR (quadratic residues mod p), the Gauss sum gives:
      |S_hat(t)|^2 = p/4 for all t != 0  (when p = 3 mod 4)

    This is the MINIMUM of sum |S_hat(t)|^4 subject to sum |S_hat(t)|^2 = mp = p(p-1)/4,
    by the Cauchy-Schwarz / power mean inequality.
    """
    m = (p - 1) // 2
    print(f"\n{'='*70}")
    print(f"FOURIER UNIFORMITY ANALYSIS at p={p}")
    print(f"{'='*70}")

    omega = cmath.exp(2j * cmath.pi / p)

    # Compute |S_hat(t)|^2 for various orientations
    for name, get_S in [("Paley", lambda: sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)),
                         ("Interval", lambda: list(range(1, m + 1)))]:
        S = get_S()
        spectrum = []
        for t in range(p):
            val = sum(omega ** (s * t) for s in S)
            spectrum.append(abs(val)**2)

        print(f"\n  {name}: S = {S}")
        print(f"    |S_hat(0)|^2 = {spectrum[0]:.4f} (= m^2 = {m**2})")
        non_zero = spectrum[1:]
        mean_nz = sum(non_zero) / len(non_zero)
        var_nz = sum((x - mean_nz)**2 for x in non_zero) / len(non_zero)
        print(f"    |S_hat(t)|^2 for t!=0:")
        print(f"      min = {min(non_zero):.4f}")
        print(f"      max = {max(non_zero):.4f}")
        print(f"      mean = {mean_nz:.4f} (expected (p-1)/4 = {(p-1)/4:.4f})")
        print(f"      variance = {var_nz:.4f}")
        print(f"      std = {var_nz**0.5:.4f}")

        # For Paley at p=3 mod 4: |S_hat(t)|^2 = (p-1)/4 for all t!=0
        if name == "Paley":
            is_uniform = var_nz < 0.01
            print(f"    Uniform spectrum (Gauss sum): {is_uniform}")
            if is_uniform:
                print(f"    This is WHY Paley minimizes E: by C-S, sum x^2 >= (sum x)^2/n")
                print(f"    with equality iff all x equal.")

        E = additive_energy(S, p)
        L4 = sum(x**2 for x in spectrum)
        print(f"    E(S) = {E}")
        print(f"    sum |S_hat|^4 = {L4:.4f} (= p*E = {p*E})")


def weil_bound_connection(p):
    """Connect to the Weil bound for character sums.

    The Weil bound states: |sum_{x} chi(f(x))| <= (deg f - 1) * sqrt(p)
    for a non-trivial character chi and polynomial f.

    For the QR set S = {x : chi(x) = 1}, the Gauss sum is:
      S_hat(t) = sum_{s in S} omega^{st}
               = (1/2) * sum_{x=1}^{p-1} (1 + chi(x)) * omega^{xt}
               = (p-1)/2 * delta_{t,0} + (1/2) * sum_{x} chi(x) * omega^{xt}

    The inner sum is a Gauss sum: g(chi, t) = chi(t) * g(chi, 1) = chi(t) * g
    where g = sum_x chi(x) omega^x has |g| = sqrt(p).

    Therefore |S_hat(t)|^2 = |g|^2/4 = p/4 for t != 0.
    This gives E(QR) = (1/p)[m^4 + (p-1) * (p/4)^2] = (1/p)[m^4 + (p-1)*p^2/16].
    """
    m = (p - 1) // 2

    print(f"\n{'='*70}")
    print(f"WEIL BOUND CONNECTION at p={p}")
    print(f"{'='*70}")

    # Exact E(QR) from the uniform spectrum
    E_paley_exact = (m**4 + (p - 1) * (p / 4)**2) / p
    E_paley_computed = additive_energy(sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1), p)

    print(f"  E(QR) formula: (m^4 + (p-1)*p^2/16) / p")
    print(f"    = ({m**4} + {(p-1)*p**2//16}) / {p}")
    print(f"    = {E_paley_exact:.4f}")
    print(f"  E(QR) computed: {E_paley_computed}")
    print(f"  Match: {abs(E_paley_exact - E_paley_computed) < 0.01}")

    # Simplify: E(QR) = m^4/p + (p-1)p/16
    E_simple = m**4 / p + (p - 1) * p / 16
    print(f"\n  Simplified: E(QR) = m^4/p + (p-1)p/16")
    print(f"    = {E_simple:.4f}")
    print(f"    = {m**4}/p + {(p-1)*p//16}")

    # What about E(Interval)?
    S_int = list(range(1, m + 1))
    E_int = additive_energy(S_int, p)
    print(f"\n  E(Interval) = {E_int}")
    print(f"  E(Interval) - E(QR) = {E_int - E_paley_computed}")
    print(f"  Ratio E(Int)/E(QR) = {E_int/E_paley_computed:.6f}")

    # The gap is explained by the non-uniformity of interval Fourier spectrum
    omega = cmath.exp(2j * cmath.pi / p)
    spec_int = []
    for t in range(1, p):
        val = sum(omega ** (s * t) for s in S_int)
        spec_int.append(abs(val)**2)

    # Variance of |S_hat|^2 for interval
    mean_spec = sum(spec_int) / len(spec_int)
    var_spec = sum((x - mean_spec)**2 for x in spec_int) / len(spec_int)

    # E_excess = (1/p) * sum_t (|S_hat(t)|^2 - p/4)^2
    #          = (p-1) * var(|S_hat|^2) / p  (approximately)
    # Since E = (sum |S_hat|^4) / p and E_min = ((sum |S_hat|^2)^2/(p-1)) / p
    excess_from_var = (p - 1) * var_spec / p
    actual_excess = E_int - E_paley_computed

    print(f"\n  Interval spectrum variance = {var_spec:.4f}")
    print(f"  Excess energy from variance = {excess_from_var:.4f}")
    print(f"  Actual excess = {actual_excess}")
    print(f"  Match: {abs(excess_from_var - actual_excess) < 0.5}")


def H_energy_monotonicity(p):
    """Test whether H is a MONOTONE DECREASING function of E.

    If so: min E (Paley) => max H, max E (Interval) => min H among regular.
    This would prove Paley maximality for circulant tournaments!
    """
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"H-ENERGY MONOTONICITY at p={p}")
    print(f"{'='*70}")

    if p > 11:
        print(f"  [p too large for H computation]")
        return

    E_H_pairs = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        E = additive_energy(S, p)
        H = compute_H(A, p)
        E_H_pairs.append((E, H, bits))

    # Group by E
    E_groups = defaultdict(list)
    for E, H, bits in E_H_pairs:
        E_groups[E].append((H, bits))

    monotone = True
    prev_H_max = float('inf')
    for E in sorted(E_groups.keys()):
        H_vals = [h for h, b in E_groups[E]]
        H_min, H_max = min(H_vals), max(H_vals)
        is_const = (H_min == H_max)

        if H_min > prev_H_max:
            monotone = False

        prev_H_max = H_max

        print(f"  E={E:4d}: H in [{H_min}, {H_max}]{' (constant)' if is_const else ''}, n={len(H_vals)}")

    print(f"\n  Monotone decreasing: {monotone}")

    # Linear regression
    E_arr = np.array([e for e, h, b in E_H_pairs], dtype=float)
    H_arr = np.array([h for e, h, b in E_H_pairs], dtype=float)
    r = np.corrcoef(E_arr, H_arr)[0, 1]
    print(f"  Corr(E, H) = {r:.6f}")

    if np.std(E_arr) > 0:
        slope = np.cov(E_arr, H_arr)[0, 1] / np.var(E_arr)
        intercept = np.mean(H_arr) - slope * np.mean(E_arr)
        resid = H_arr - (slope * E_arr + intercept)
        print(f"  Linear: H = {slope:.4f}*E + {intercept:.2f}, max_resid = {np.max(np.abs(resid)):.4f}")


# ================================================================
# MAIN
# ================================================================

print("=" * 70)
print("PALEY AS ENERGY MINIMIZER: E(QR) <= E(S) for all orientations")
print("=" * 70)

# H vs Energy analysis (only feasible at small p)
for p in [7, 11]:
    analyze_H_vs_energy(p)

# Fourier uniformity
for p in [7, 11, 13]:
    energy_vs_fourier_uniformity(p)

# Weil bound connection
for p in [7, 11, 13, 19]:
    weil_bound_connection(p)

# H-Energy monotonicity
for p in [7, 11]:
    H_energy_monotonicity(p)

print("\nDONE.")
