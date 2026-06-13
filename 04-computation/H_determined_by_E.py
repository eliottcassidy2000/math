#!/usr/bin/env python3
"""
H_determined_by_E.py -- Is H(sigma) determined by additive energy E(S)?

DISCOVERY from ipr_energy_identity.py:
At p=7 and p=11, all orientations with the same E(S) have the SAME H.
This would mean H is a FUNCTION of E(S) alone.

If true, this reduces H-maximization to E(S)-maximization, which is a
well-studied problem in additive combinatorics. Interval maximizes E(S)
(contiguous sets have max additive energy by rearrangement inequality).

But WAIT: at p=7, max E corresponds to MIN H (anti-correlation)!
At p=11, the correlation is weak. At p=13, it's positive.
So E determines H but NOT monotonically.

This script verifies the E→H determination at all available primes
and investigates the exact functional form H = f(E).

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def additive_energy(S, p):
    S_set = set(S)
    e = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    e += 1
    return e


def representation_function(S, p):
    """r_S(n) = |{(a,b) in S x S : a + b = n mod p}|"""
    S_set = set(S)
    r = [0] * p
    for a in S:
        for b in S:
            r[(a + b) % p] += 1
    return r


def higher_energies(S, p):
    """Compute E_k(S) for k=2,3:
    E_2 = E(S) = sum r(n)^2 = standard additive energy
    E_3 = sum r(n)^3 (cubic energy)
    """
    r = representation_function(S, p)
    E2 = sum(x**2 for x in r)
    E3 = sum(x**3 for x in r)
    E4 = sum(x**4 for x in r)
    return E2, E3, E4


def main():
    print("=" * 70)
    print("IS H DETERMINED BY E(S)? (Additive Energy)")
    print("=" * 70)

    # PART 1: Full verification at p=7, 11, 13
    print("\n--- PART 1: FULL E -> H MAPPING ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}, m={m}, 2^m = {1 << m} orientations:")

        EH_data = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            E = additive_energy(S, p)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            EH_data.append((E, H, S))

        # Group by E
        E_groups = defaultdict(list)
        for E, H, S in EH_data:
            E_groups[E].append((H, S))

        # Check: is H constant within each E group?
        all_determined = True
        print(f"    {'E':>6} {'count':>6} {'H values':>20} {'determined':>12}")
        for E in sorted(E_groups):
            items = E_groups[E]
            H_vals = sorted(set(h for h, s in items))
            determined = len(H_vals) == 1
            if not determined:
                all_determined = False
            print(f"    {E:>6} {len(items):>6} {str(H_vals):>20} {'YES' if determined else 'NO ***'}")

        print(f"\n    H is {'FULLY' if all_determined else 'NOT'} determined by E(S) at p={p}")

    # PART 2: If H = f(E), what is the function?
    print("\n--- PART 2: THE FUNCTION H = f(E) ---")
    print("If H is determined by E, we can write H = a_0 + a_1*E + a_2*E^2 + ...")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        EH_map = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            EH_map[E] = H

        E_vals = sorted(EH_map.keys())
        H_vals = [EH_map[e] for e in E_vals]

        print(f"\n  p={p}: E -> H mapping:")
        for e, h in zip(E_vals, H_vals):
            print(f"    E={e:>6} -> H={h:>10}")

        # Check linearity
        if len(E_vals) >= 2:
            dE = [E_vals[i+1] - E_vals[i] for i in range(len(E_vals)-1)]
            dH = [H_vals[i+1] - H_vals[i] for i in range(len(H_vals)-1)]
            slopes = [dh/de if de != 0 else float('inf') for dh, de in zip(dH, dE)]
            print(f"    Slopes dH/dE: {[f'{s:.2f}' for s in slopes]}")

            if len(slopes) >= 2:
                is_linear = all(abs(s - slopes[0]) < 0.01 * abs(slopes[0]) for s in slopes)
                print(f"    Linear? {is_linear}")

    # PART 3: Is it the REPRESENTATION function that determines H?
    print("\n--- PART 3: REPRESENTATION FUNCTION r_S ---")
    print("r_S(n) = |{(a,b) in S x S : a+b = n mod p}|")
    print("E(S) = sum r_S(n)^2. But is H determined by the full r_S profile?")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}:")

        profiles = defaultdict(list)
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            r = representation_function(S, p)
            # The profile (sorted r values for n > 0) characterizes r_S up to relabeling
            profile = tuple(sorted(r[1:]))  # exclude r(0) = m (always)
            profiles[profile].append((H, S))

        print(f"    Distinct r-profiles: {len(profiles)}")
        for profile, items in sorted(profiles.items()):
            H_vals = sorted(set(h for h, s in items))
            count = len(items)
            E = sum(x**2 for x in profile) + m**2  # include r(0)=m
            print(f"    profile={profile}, count={count}, E={E}, H={H_vals}")

    # PART 4: Higher energies
    print("\n--- PART 4: HIGHER ENERGY INVARIANTS ---")
    print("E_2 = sum r(n)^2 (standard additive energy)")
    print("E_3 = sum r(n)^3 (cubic energy)")
    print("E_4 = sum r(n)^4 (quartic energy)")
    print("Question: do (E_2, E_3) determine H better?")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}:")

        data = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            E2, E3, E4 = higher_energies(S, p)
            data.append({'E2': E2, 'E3': E3, 'E4': E4, 'H': H, 'S': S})

        # Group by (E2, E3)
        E23_groups = defaultdict(list)
        for d in data:
            E23_groups[(d['E2'], d['E3'])].append(d['H'])

        print(f"    Distinct (E2,E3) pairs: {len(E23_groups)}")
        all_det = True
        for key in sorted(E23_groups):
            H_vals = sorted(set(E23_groups[key]))
            if len(H_vals) > 1:
                all_det = False
            print(f"    (E2={key[0]}, E3={key[1]}): H={H_vals}, count={len(E23_groups[key])}")

        print(f"    H determined by (E2, E3)? {all_det}")

    # PART 5: The Fourier spectrum determines everything
    print("\n--- PART 5: FOURIER SPECTRUM ANALYSIS ---")
    print("The sorted Fourier spectrum {|S_hat(k)|^2}_{k>0} should determine H")
    print("(up to automorphisms of Z_p)")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n  p={p}:")

        spec_groups = defaultdict(list)
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            # Fourier spectrum (k > 0)
            spectrum = []
            for k in range(1, (p+1)//2 + 1):
                val = sum(omega ** (k * s) for s in S)
                spectrum.append(round(abs(val)**2, 6))
            spec_key = tuple(sorted(spectrum))
            spec_groups[spec_key].append((H, S))

        print(f"    Distinct sorted Fourier spectra: {len(spec_groups)}")
        for spec, items in sorted(spec_groups.items()):
            H_vals = sorted(set(h for h, s in items))
            count = len(items)
            E_approx = round(sum(x**2 for x in spec) * 2 + m**4)  # approximate
            print(f"    spectrum={tuple(f'{x:.2f}' for x in spec[:4])}..., "
                  f"count={count}, H={H_vals}")

    # PART 6: The functional relationship H = f(E) at small primes
    print("\n--- PART 6: H vs E INTERPOLATION ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        EH_map = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            EH_map[E] = H

        E_sorted = sorted(EH_map.keys())
        H_sorted = [EH_map[e] for e in E_sorted]
        n_pts = len(E_sorted)

        print(f"\n  p={p}, m={m}:")
        print(f"    {n_pts} distinct (E, H) points")

        # Fit polynomial of degree n_pts - 1
        # Use Lagrange interpolation
        if n_pts <= 6:
            # Check if H is affine in E
            if n_pts >= 2:
                # Linear fit: H = a + b*E
                E_mean = sum(E_sorted) / n_pts
                H_mean = sum(H_sorted) / n_pts
                cov = sum((e - E_mean) * (h - H_mean) for e, h in zip(E_sorted, H_sorted))
                var_E = sum((e - E_mean)**2 for e in E_sorted)
                b = cov / var_E if var_E > 0 else 0
                a = H_mean - b * E_mean
                residuals = [h - (a + b*e) for e, h in zip(E_sorted, H_sorted)]
                max_res = max(abs(r) for r in residuals)

                print(f"    Linear fit: H = {a:.2f} + {b:.2f} * E")
                print(f"    Max residual: {max_res:.2f}")
                print(f"    Linear? {max_res < 0.01 * abs(H_sorted[0])}")

            # Check quadratic fit: H = a + b*E + c*E^2
            if n_pts >= 3:
                # Simple Vandermonde solve for 3 points
                # Actually let's just check residuals for first 3 points
                E0, E1, E2 = E_sorted[0], E_sorted[1], E_sorted[2]
                H0, H1, H2 = H_sorted[0], H_sorted[1], H_sorted[2]

                # Lagrange coefficients for quadratic
                L0 = lambda e: (e-E1)*(e-E2)/((E0-E1)*(E0-E2))
                L1 = lambda e: (e-E0)*(e-E2)/((E1-E0)*(E1-E2))
                L2 = lambda e: (e-E0)*(e-E1)/((E2-E0)*(E2-E1))

                quad = lambda e: H0*L0(e) + H1*L1(e) + H2*L2(e)
                quad_residuals = [h - quad(e) for e, h in zip(E_sorted, H_sorted)]
                max_quad_res = max(abs(r) for r in quad_residuals)

                print(f"    Quadratic fit residuals: {[f'{r:.1f}' for r in quad_residuals]}")
                print(f"    Max quadratic residual: {max_quad_res:.2f}")

    # PART 7: Connection to the E→H reversal via Walsh spectrum
    print("\n--- PART 7: WHY E anti-correlates with H at p=3 mod 4 ---")
    print("At p=7 (3 mod 4): max E = Interval, but max H = Paley")
    print("At p=13 (1 mod 4): max E = Interval = max H")
    print()
    print("Insight: Walsh degree-2 contribution f_2(sigma) depends on")
    print("the Legendre character chi at p=3 mod 4, creating sign reversal.")
    print("At p=1 mod 4, chi(-1)=+1 so no reversal; degree-4 takes over.")

    for p in [7, 11, 13]:
        m = (p - 1) // 2

        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        chi_neg1 = 1 if (p - 1) in QR else -1

        print(f"\n  p={p} ({p % 4} mod 4): chi(-1) = {chi_neg1}")
        print(f"    At p=3 mod 4: chi(-1)=-1, Paley tournament exists")
        print(f"    At p=1 mod 4: chi(-1)=+1, NO Paley tournament")
        print(f"    The sign reversal in Walsh degree-2 is controlled by chi(-1)")


if __name__ == '__main__':
    main()
