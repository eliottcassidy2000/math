#!/usr/bin/env python3
"""
walsh_overlap_bridge.py -- Connecting Walsh decomposition to overlap/co-occurrence structure

KEY QUESTION: How does the Walsh spectrum of H(sigma) relate to the overlap
structure of the conflict graph Omega(T)?

OBSERVATIONS FROM PRIOR ANALYSIS:
1. Paley co-occurrence variance = 0 (flat) vs Interval variance > 0
2. Paley Fourier L4/L2^2 ≈ 1/(p-1) (near-uniform) vs Interval ≈ 0.31
3. Additive energy: Paley < Interval at all tested primes
4. Paradox: Interval has MORE disjoint 3-3 pairs but FEWER H-paths at p=3 mod 4

HYPOTHESIS: The Walsh degree-2 coefficients h_hat[{i,j}] encode the
pairwise "chord interaction" between chords i and j. The degree-2 Walsh
energy is related to the co-occurrence variance of the connection set.

At p=3 mod 4 (Paley wins):
  - Walsh coefficients have UNIFORM magnitude within each degree (1 orbit)
  - Sign law holds => Paley achieves perfect alignment
  - Co-occurrence is flat => all chord pairs contribute equally

At p=1 mod 4 (Interval wins):
  - Walsh coefficients have MULTIPLE magnitudes (multiple orbits)
  - No sign law => no perfect alignment possible
  - Degree-4 energy dominates degree-2

This script:
1. Computes both Walsh decomposition and co-occurrence structure simultaneously
2. Tests if Walsh degree-2 energy = variance of co-occurrence
3. Tests if there's an exact formula relating them
4. Explores the interaction between chord pairs and higher-degree Walsh terms

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
from math import comb
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


def full_walsh(p):
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    H_dict = {}
    for bits in range(1 << m):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if sigma[i] == 1 else b)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_dict[sigma] = H

    n = 1 << m
    h_hat = {}
    for bits in range(n):
        S_idx = tuple(i for i in range(m) if bits & (1 << i))
        coeff = 0
        for sigma_bits in range(n):
            sigma = tuple((1 if sigma_bits & (1 << i) else -1) for i in range(m))
            H = H_dict[sigma]
            prod = 1
            for i in S_idx:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[S_idx] = coeff / n

    return h_hat, H_dict


def co_occurrence_structure(p, S):
    """Compute the 3-cycle co-occurrence structure on Z_p for connection set S."""
    A = build_adj(p, S)

    # Count 3-cycle vertex sets
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))

    # Co-occurrence by gap (using circulant symmetry: co_occur(u,v) depends on v-u mod p)
    gap_co = [0] * p
    for fs in c3_sets:
        verts = list(fs)
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                d = (verts[j] - verts[i]) % p
                gap_co[d] += 1
                d2 = (verts[i] - verts[j]) % p
                gap_co[d2] += 1

    # Normalize by p (circulant symmetry: each gap appears with multiplicity p)
    # Actually for circulant, co_occur(0,d) counts cycles containing both 0 and d.
    # By circulant symmetry, this equals co_occur(u, u+d) for any u.
    # Direct computation using vertex 0:
    gap_co_0 = [0] * p
    for fs in c3_sets:
        if 0 in fs:
            for v in fs:
                if v != 0:
                    gap_co_0[v] += 1

    return c3_sets, gap_co_0


def fourier_of_S(S, p):
    """Compute Fourier transform of indicator of S."""
    omega = cmath.exp(2j * cmath.pi / p)
    S_hat = []
    for k in range(p):
        val = sum(omega ** (k * s) for s in S)
        S_hat.append(val)
    return S_hat


def main():
    print("=" * 70)
    print("WALSH-OVERLAP BRIDGE ANALYSIS")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, p mod 4 = {p%4}")
        print(f"{'='*70}")

        # Compute Walsh decomposition
        h_hat, H_dict = full_walsh(p)

        # Walsh energy by degree
        energy_by_deg = defaultdict(float)
        for S_idx, c in h_hat.items():
            energy_by_deg[len(S_idx)] += c**2

        print(f"\n  Walsh energy by degree:")
        total_energy = sum(energy_by_deg.values())
        for d in sorted(energy_by_deg):
            pct = 100 * energy_by_deg[d] / total_energy
            print(f"    deg {d}: {energy_by_deg[d]:>14.2f} ({pct:.4f}%)")

        # Degree-2 Walsh matrix W2[i,j] = h_hat[{i,j}]
        print(f"\n  Degree-2 Walsh interaction matrix h_hat[{{i,j}}]:")
        W2 = [[0.0]*m for _ in range(m)]
        for i in range(m):
            for j in range(i+1, m):
                W2[i][j] = h_hat.get((i,j), 0)
                W2[j][i] = W2[i][j]

        for i in range(m):
            row = " ".join(f"{W2[i][j]:>+10.2f}" for j in range(m))
            print(f"    {row}")

        # Gap product classification
        print(f"\n  Gap product -> Walsh coefficient sign (degree 2):")
        for i in range(m):
            for j in range(i+1, m):
                c = h_hat.get((i,j), 0)
                if abs(c) < 0.01:
                    continue
                g_prod = ((i+1) * (j+1)) % p
                is_qr = g_prod in QR
                sign_c = "+" if c > 0 else "-"
                sign_chi = "QR" if is_qr else "NQR"
                match = (c > 0) == is_qr
                print(f"    {{i={i},j={j}}}: gaps=({i+1},{j+1}), "
                      f"prod={g_prod}, {sign_chi}, h={c:>+.2f} [{sign_c}] "
                      f"{'MATCH' if match else 'MISMATCH'}")

        # ====== CO-OCCURRENCE vs WALSH ======
        print(f"\n  --- Co-occurrence vs Walsh for each tournament ---")

        S_int = list(range(1, m + 1))
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            c3_sets, gap_co = co_occurrence_structure(p, S)

            # Gap co-occurrence values
            co_vals = [gap_co[d] for d in range(1, p)]
            co_mean = sum(co_vals) / len(co_vals)
            co_var = sum((c - co_mean)**2 for c in co_vals) / len(co_vals)

            # Fourier spectrum
            S_hat = fourier_of_S(S, p)
            fourier_sq = [abs(S_hat[k])**2 for k in range(p)]
            L4 = sum(f**2 for f in fourier_sq)
            L2 = sum(f for f in fourier_sq)

            print(f"\n    {name} (S={S}):")
            print(f"      c3 sets = {len(c3_sets)}")
            print(f"      Co-occ variance = {co_var:.4f}")
            print(f"      Fourier L4/L2^2 = {L4/L2**2:.6f}")

            # Compare: co-occurrence variance vs Walsh deg-2 energy
            walsh_deg2_energy = energy_by_deg.get(2, 0)
            print(f"      Walsh deg-2 energy = {walsh_deg2_energy:.4f}")

            # Check: is deg-2 energy related to variance of |S_hat|^2?
            fourier_sq_nonzero = [fourier_sq[k] for k in range(1, p)]
            fourier_mean = sum(fourier_sq_nonzero) / len(fourier_sq_nonzero)
            fourier_var = sum((f - fourier_mean)**2 for f in fourier_sq_nonzero) / len(fourier_sq_nonzero)
            print(f"      |S_hat|^2 variance = {fourier_var:.6f}")
            print(f"      |S_hat|^2 mean = {fourier_mean:.6f}")

        # ====== CHORD INTERACTION DECOMPOSITION ======
        print(f"\n  --- Chord interaction decomposition ---")
        # When we flip chord i: gap (i+1) <-> gap (p-i-1)
        # Effect on S_hat(k): delta = omega^{-k*(i+1)} - omega^{k*(i+1)} = -2i*sin(2pi*k*(i+1)/p)
        # This is a rank-1 update to the Fourier spectrum.
        # The Walsh degree-2 coefficient h_hat[{i,j}] measures the joint effect
        # of flipping chords i and j on H.

        omega = cmath.exp(2j * cmath.pi / p)
        for i in range(m):
            for j in range(i+1, m):
                c = h_hat.get((i,j), 0)
                if abs(c) < 0.01:
                    continue
                g_i = i + 1
                g_j = j + 1

                # The "spectral interaction" between chords i and j
                # When both flip: S_hat(k) changes by delta_i(k) + delta_j(k)
                # where delta_i(k) = omega^{-k*g_i} - omega^{k*g_i}
                # The quadratic term in the H change is proportional to
                # sum_k delta_i(k) * delta_j(k) * (something)

                spectral_int = 0
                for k in range(1, p):
                    di = omega**(-k*g_i) - omega**(k*g_i)  # = -2i*sin(2pi*k*g_i/p)
                    dj = omega**(-k*g_j) - omega**(k*g_j)
                    spectral_int += (di * dj.conjugate()).real
                spectral_int /= p

                print(f"    {{i={i},j={j}}}: h_hat = {c:>+10.2f}, "
                      f"spectral_int = {spectral_int:>+10.4f}")

        # ====== ADDITIVE ENERGY DECOMPOSITION ======
        print(f"\n  --- Additive energy and Walsh ---")
        # Additive energy E(S) = |{(a,b,c,d) in S^4 : a+b=c+d}| = sum_k |S_hat(k)|^4 / p
        # By Parseval: sum |S_hat(k)|^2 = m*p
        # The Fourier L4 norm controls additive structure.

        # For Walsh: h_hat[S] depends on S through the Fourier spectrum.
        # The key identity: for circulant tournaments,
        #   H(sigma) = f(S_hat(1), ..., S_hat(m))
        # where f is a symmetric function of the eigenspace contributions.

        # The Walsh coefficient at degree d measures the d-th order polynomial
        # interaction of the chord choices.

        # CONJECTURE: Walsh deg-2 energy = C * Var(|S_hat|^2) for some constant C?
        # Let's test numerically.

        # For each sigma, compute |S_hat(k)|^2 and H
        print(f"\n  --- |S_hat|^2 spectrum vs H correlation ---")
        specs = []
        H_vals_list = []
        for sigma, H in H_dict.items():
            S = []
            for i in range(m):
                if sigma[i] == 1:
                    S.append(i + 1)
                else:
                    S.append(p - i - 1)
            S = sorted(S)
            S_hat = fourier_of_S(S, p)
            spec = tuple(round(abs(S_hat[k])**2, 6) for k in range(1, m+1))
            specs.append(spec)
            H_vals_list.append(H)

        # Check: is H determined by the |S_hat|^2 spectrum?
        spec_H = defaultdict(set)
        for spec, H in zip(specs, H_vals_list):
            # Spectrum needs to be sorted (since H depends on multiset, not ordered tuple)
            sorted_spec = tuple(sorted(spec))
            spec_H[sorted_spec].add(H)

        n_determined = sum(1 for v in spec_H.values() if len(v) == 1)
        n_total = len(spec_H)
        print(f"    H determined by sorted |S_hat|^2 spectrum? "
              f"{n_determined}/{n_total} spectra have unique H")

        if any(len(v) > 1 for v in spec_H.values()):
            for spec, H_set in spec_H.items():
                if len(H_set) > 1:
                    print(f"    AMBIGUOUS: spectrum {spec[:3]}... -> H values {H_set}")

        # ====== DEGREE-4 WALSH AND HIGHER-ORDER INTERACTIONS ======
        print(f"\n  --- Degree-4 Walsh coefficients ---")
        deg4_items = [(S_idx, c) for S_idx, c in h_hat.items()
                      if len(S_idx) == 4 and abs(c) > 0.01]
        if deg4_items:
            # Group by magnitude
            mag_groups = defaultdict(list)
            for S_idx, c in deg4_items:
                mag = round(abs(c), 2)
                mag_groups[mag].append((S_idx, c))

            for mag in sorted(mag_groups, reverse=True):
                group = mag_groups[mag]
                n_pos = sum(1 for _, c in group if c > 0)
                n_neg = sum(1 for _, c in group if c < 0)
                print(f"    |h| = {mag}: {len(group)} coeffs ({n_pos}+, {n_neg}-)")

                # Check QR product structure
                for S_idx, c in group[:3]:  # show first 3
                    prod = 1
                    for i in S_idx:
                        prod = (prod * (i + 1)) % p
                    is_qr = prod in QR
                    print(f"      S={set(S_idx)}, prod gaps = {prod}, "
                          f"{'QR' if is_qr else 'NQR'}, h={c:>+.2f}")
        else:
            print(f"    No degree-4 coefficients (m={m} < 4)")

        # ====== THE BRIDGE: WALSH COEFFICIENTS AS DERIVATIVES ======
        print(f"\n  --- Walsh as partial derivatives of H ---")
        # h_hat[{i}] = (1/2) * (H(sigma_i=+1) - H(sigma_i=-1))_avg
        # h_hat[{i,j}] = (1/4) * avg_{sigma_{-ij}} [H(++)-H(+-)-H(-+)+H(--)]
        # This is the "interaction effect" of chords i and j.

        # Compute the interaction effect for each chord pair
        for i in range(m):
            for j in range(i+1, m):
                c = h_hat.get((i,j), 0)
                if abs(c) < 0.01:
                    continue

                # Average over all other chords
                plus_plus = 0
                plus_minus = 0
                minus_plus = 0
                minus_minus = 0
                count = 0

                for sigma, H in H_dict.items():
                    if sigma[i] == 1 and sigma[j] == 1:
                        plus_plus += H
                    elif sigma[i] == 1 and sigma[j] == -1:
                        plus_minus += H
                    elif sigma[i] == -1 and sigma[j] == 1:
                        minus_plus += H
                    else:
                        minus_minus += H

                n_each = (1 << m) // 4
                pp_avg = plus_plus / n_each
                pm_avg = plus_minus / n_each
                mp_avg = minus_plus / n_each
                mm_avg = minus_minus / n_each
                interaction = (pp_avg - pm_avg - mp_avg + mm_avg) / 4

                print(f"    chords ({i},{j}): "
                      f"avg H(++)={pp_avg:.1f}, H(+-)={pm_avg:.1f}, "
                      f"H(-+)={mp_avg:.1f}, H(--)={mm_avg:.1f}, "
                      f"interaction={interaction:.2f} [Walsh={c:.2f}]")

        # ====== OVERLAP WEIGHT vs WALSH ======
        print(f"\n  --- Overlap weight structure ---")
        # The overlap weight between two cycles C_i, C_j is |V(C_i) ∩ V(C_j)|.
        # For 3-cycles: overlap is 0, 1, or 2.
        # Overlap=0 means disjoint = independent in Omega.
        #
        # KEY: The Walsh degree-2 coefficient h_hat[{i,j}] measures how much
        # the joint chord choice (i,j) affects H. This is related to how
        # the cycle structure changes when we swap chords.
        #
        # When chord i flips: some cycles disappear, some appear.
        # When chords i and j BOTH flip: the interaction term captures
        # the "synergy" between these changes.
        #
        # HYPOTHESIS: |h_hat[{i,j}]| is large when chords i and j are
        # "structurally entangled" (many cycles use both gaps i+1 and j+1).

        # Compute: for each chord pair (i,j), count cycles using both gaps
        S_int = list(range(1, m + 1))
        A_int = build_adj(p, S_int)

        # Get all cycle vertex sets for Interval
        all_cycles = []
        max_k = p if p <= 11 else 7
        for k in range(3, max_k + 1, 2):
            for subset in combinations(range(p), k):
                verts = list(subset)
                n_cyc = count_ham_cycles_on_subset(A_int, verts)
                for _ in range(n_cyc):
                    all_cycles.append(frozenset(subset))

        if all_cycles:
            # For each cycle, find which gaps from S_int it uses
            gap_usage = defaultdict(set)  # gap -> set of cycle indices
            for idx, fs in enumerate(all_cycles):
                verts = sorted(fs)
                for vi in range(len(verts)):
                    for vj in range(len(verts)):
                        if vi == vj:
                            continue
                        d = (verts[vj] - verts[vi]) % p
                        if d in set(S_int):
                            gap_usage[d].add(idx)

            print(f"    Cycles using each gap (Interval):")
            for g in S_int:
                print(f"      gap {g}: {len(gap_usage[g])} cycles")

            # For each pair of gaps, count cycles using both
            print(f"\n    Cycle co-usage of gap pairs (Interval):")
            print(f"    {'':>6}", end="")
            for j in range(m):
                print(f" gap{j+1:>3}", end="")
            print()
            for i in range(m):
                print(f"    gap{i+1:>2}", end="")
                for j in range(m):
                    if i == j:
                        n_both = len(gap_usage[i+1])
                    else:
                        n_both = len(gap_usage[i+1] & gap_usage[j+1])
                    print(f" {n_both:>5}", end="")
                print()

            # Compare to Walsh deg-2
            print(f"\n    Gap co-usage vs Walsh h_hat[{{i,j}}]:")
            for i in range(m):
                for j in range(i+1, m):
                    n_both = len(gap_usage[i+1] & gap_usage[j+1])
                    c = h_hat.get((i,j), 0)
                    if abs(c) > 0.01:
                        print(f"      gaps ({i+1},{j+1}): {n_both} co-uses, "
                              f"h_hat = {c:>+.2f}")


def count_ham_cycles_on_subset(A, verts):
    """Count directed Ham cycles on a vertex subset."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        return fwd + bwd

    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


if __name__ == '__main__':
    main()
