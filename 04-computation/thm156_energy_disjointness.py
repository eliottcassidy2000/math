#!/usr/bin/env python3
"""
thm156_energy_disjointness.py -- kind-pasteur-2026-03-13-S60

THEOREM (THM-156): For circulant regular tournaments on Z_p,

    disj3(S) = (p/4) * E(S) + C(p)

where E(S) is the additive energy of the connection set S, and
C(p) = (K(p) + c3(p)*(c3(p)-1)) / 2 with:
  K(p) = -3p(p^2-1)(p^2-9)/320
  c3(p) = p(p^2-1)/24

Proof chain:
1. disj3 = -T5/2 + A(p)  where T5 = tr(A^5)/5 = # directed 5-cycles
   (from THM-155 applied to the Walsh level)
2. E*p = m^4 + m/8 + mp/4 + 2*S4  where S4 = sum_{t=1}^m D_t^4
   (from Parseval: E = (1/p)*sum|lambda_t|^4 and |lambda_t|^2 = 1/4 + D_t^2)
3. 5*T5 = m^5 - m/16 + 5mp/8 - 5*S4
   (from eigenvalue expansion of tr(A^5))
4. Eliminating S4 gives T5 = B(p) - p*E/2, hence disj3 = p*E/4 + const(p)

This gives an EXACT bridge between tournament cycle structure and additive combinatorics.
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
    """E(S) = |{(a,b,c,d) in S^4 : a+b=c+d mod p}|"""
    S_set = set(S)
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1
    return energy


def count_c3_disj3(A, p):
    """Count 3-cycle vertex sets and disjoint pairs."""
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


def directed_5cycles(A, p):
    """Compute T5 = tr(A^5)/5 = number of directed 5-cycles."""
    A_np = np.array(A, dtype=np.float64)
    A5 = np.linalg.matrix_power(A_np, 5)
    return int(round(np.trace(A5))) // 5


def const_C(p):
    """Compute the constant C(p) in disj3 = (p/4)*E + C(p).

    Derivation:
      disj3 = -T5/2 + A(p)       where A(p) = (K+c3*(c3-1))/2
      T5 = B(p) - p*E/2          where B(p) = m^5/5 + m^4/2 + m/20 + mp/4
      => disj3 = p*E/4 + A(p) - B(p)/2

    So C(p) = A(p) - B(p)/2
    = (K+c3*(c3-1))/2 - m^5/10 - m^4/4 - m/40 - mp/8
    """
    m = (p - 1) // 2
    c3 = p * (p**2 - 1) // 24
    K = -3 * p * (p**2 - 1) * (p**2 - 9) // 320
    A = (K + c3 * (c3 - 1)) / 2
    B_half = m**5 / 10 + m**4 / 4 + m / 40 + m * p / 8
    return A - B_half


def const_B(p):
    """B(p) = (2m^5 + 5m^4 + m/2 + 5mp/2) / 10 where m = (p-1)/2."""
    m = (p - 1) // 2
    # To avoid fractions, compute 10*B:
    # 10B = 2m^5 + 5m^4 + 5m + 5*5*m*p  ... no
    # 10B = 2m^5 + 5m^4 + 5m + 25mp ... no
    # Let me be precise: 10*B = 2m^5 + 5m^4 + 5m + 25mp
    # Actually: m/2 * 10 = 5m, and 5mp/2 * 10 = 25mp
    ten_B = 2*m**5 + 5*m**4 + 5*m + 25*m*p
    return ten_B / 10


def verify_thm156(p, verbose=True):
    """Verify THM-156 across all orientations."""
    m = (p - 1) // 2
    N = 1 << m
    C = const_C(p)
    c3_expected = p * (p**2 - 1) // 24

    if verbose:
        print(f"\n{'='*70}")
        print(f"THM-156 VERIFICATION at p={p}, m={m}")
        print(f"{'='*70}")
        print(f"  c3(p) = {c3_expected}")
        print(f"  C(p) = {C}")
        print(f"  Formula: disj3 = (p/4)*E + C(p)")

    max_err = 0
    n_checked = 0
    results = []

    limit = min(N, 128)  # cap for large p
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        c3, disj3 = count_c3_disj3(A, p)
        E = additive_energy(S, p)
        T5 = directed_5cycles(A, p)

        predicted = p * E / 4 + C
        err = abs(disj3 - predicted)
        max_err = max(max_err, err)
        n_checked += 1

        if err > 0.01 and verbose:
            print(f"  FAIL: bits={bits:0{m}b} disj3={disj3} predicted={predicted:.4f} err={err:.4f}")

        results.append((bits, E, T5, disj3, predicted))

    if verbose:
        print(f"\n  Checked {n_checked}/{N} orientations")
        print(f"  Max error = {max_err:.6f}")
        if max_err < 0.01:
            print(f"  *** THM-156 VERIFIED: disj3 = (p/4)*E(S) + {C} ***")
        else:
            print(f"  THM-156 FAILED at p={p}")

        # Show a few examples
        print(f"\n  Sample data:")
        print(f"    {'bits':>{m+2}}  {'E':>6}  {'T5':>6}  {'disj3':>6}  {'pred':>8}  {'err':>8}")
        for bits, E, T5, disj3, pred in results[:8]:
            print(f"    {bits:0{m}b}  {E:6d}  {T5:6d}  {disj3:6d}  {pred:8.2f}  {abs(disj3-pred):8.4f}")

    return max_err < 0.01, max_err


def closed_form_const(p):
    """Compute and verify the closed form for C(p)."""
    m = (p - 1) // 2
    c3 = p * (p**2 - 1) // 24
    K = -3 * p * (p**2 - 1) * (p**2 - 9) // 320

    # C(p) = (K + c3*(c3-1))/2
    C_formula = (K + c3 * (c3 - 1)) // 2

    # Alternatively, express as polynomial in p:
    # c3 = p(p^2-1)/24
    # c3*(c3-1) = p(p^2-1)/24 * (p(p^2-1)/24 - 1)
    #           = p(p^2-1)/24 * (p(p^2-1) - 24) / 24
    #           = p(p^2-1)(p^3 - p - 24) / 576

    # K = -3p(p^2-1)(p^2-9)/320

    # C = (K + c3(c3-1))/2
    # = (-3p(p^2-1)(p^2-9)/320 + p(p^2-1)(p^3-p-24)/576) / 2

    # Factor out p(p^2-1)/2:
    # C = p(p^2-1)/2 * [-3(p^2-9)/320 + (p^3-p-24)/576] / 2
    # Wait that's getting messy. Just compute numerically.

    return C_formula


def energy_range_analysis(p):
    """Analyze the range of additive energy across orientations."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n  Energy range at p={p}:")

    energies = []
    limit = min(N, 128)
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        E = additive_energy(S, p)
        energies.append(E)

    E_min, E_max = min(energies), max(energies)
    E_unique = sorted(set(energies))

    print(f"    E range: [{E_min}, {E_max}]")
    print(f"    E unique values: {E_unique}")
    print(f"    E step size: {[E_unique[i+1]-E_unique[i] for i in range(len(E_unique)-1)]}")

    # Paley energy (bits with S = QR set)
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_paley = sorted(QR & set(range(1, m+1)) | set(p-j for j in QR if j > m))
    # Actually, to find the Paley bits:
    paley_bits = 0
    for j in range(m):
        g = j + 1
        if pow(g, (p-1)//2, p) == 1:
            paley_bits |= (1 << j)
    S_paley_check = []
    for j in range(m):
        if paley_bits & (1 << j):
            S_paley_check.append(j + 1)
        else:
            S_paley_check.append(p - (j + 1))

    E_paley = additive_energy(S_paley_check, p)
    print(f"    Paley energy = {E_paley} (bits={paley_bits:0{m}b})")
    if N-1 < limit:
        print(f"    Interval energy = {energies[N-1]} (bits={'1'*m})")
    else:
        # Compute interval energy directly
        S_int = list(range(1, m + 1))
        E_int = additive_energy(S_int, p)
        print(f"    Interval energy = {E_int} (bits={'1'*m}, computed directly)")

    return E_unique


def thm155_bridge(p):
    """Show the connection chain: THM-155 <-> THM-156 <-> additive energy."""
    m = (p - 1) // 2
    c3 = p * (p**2 - 1) // 24
    K = -3 * p * (p**2 - 1) * (p**2 - 9) // 320

    print(f"\n{'='*70}")
    print(f"THM-155 <-> THM-156 BRIDGE at p={p}")
    print(f"{'='*70}")
    print(f"  THM-155: c5_dir + 2*ov1 + 2*ov2 = const  (K = {K})")
    print(f"  Rewrite: T5 + 2*(C(c3,2) - disj3) = K + c3*(c3-1)")
    print(f"           T5 - 2*disj3 = K + c3*(c3-1) - 2*C(c3,2)")
    print(f"           T5 - 2*disj3 = K - c3*(c3-1) + 2*c3*(c3-1) - 2*c3*(c3-1)/2")

    # Actually: overlap = ov1+ov2, and disj3 + overlap = C(c3,2) = c3(c3-1)/2
    # THM-155: T5 - 2*overlap = K_0  (for some constant K_0)
    # But wait, THM-155 as stated: c5 - 2*ov1 - 2*ov2 = K
    # where c5 is the number of 5-cycle VERTEX SETS? Or directed 5-cycles?

    # Let me just verify the chain numerically
    N = 1 << m
    print(f"\n  Numerical verification:")
    for bits in [0, 1, N-1]:
        if bits >= N:
            continue
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        c3_val, disj3 = count_c3_disj3(A, p)
        T5 = directed_5cycles(A, p)
        E = additive_energy(S, p)

        overlap = c3_val * (c3_val - 1) // 2 - disj3

        # THM-155: T5 - 2*overlap = ??
        lhs_155 = T5 - 2 * overlap
        # THM-156: disj3 = (p/4)*E + C
        C = const_C(p)
        lhs_156 = disj3 - p * E // 4  # should be C (approximately)

        print(f"    bits={bits:0{m}b}: T5={T5}, disj3={disj3}, E={E}")
        print(f"      T5 - 2*overlap = {lhs_155}")
        print(f"      disj3 - p*E/4 = {disj3 - p*E/4:.2f} (expected C={C})")


def explore_energy_spectrum_structure(p):
    """Deeper analysis: what determines the additive energy values?

    The additive energy E(S) = (1/p) * sum_t |S_hat(t)|^4
    = (1/p) * [m^4 + 2*sum_{t=1}^m (1/4 + D_t^2)^2]

    Since sum D_t^2 = mp/4 (constant), E depends on sum D_t^4.
    sum D_t^4 >= (sum D_t^2)^2 / m = (mp/4)^2 / m = m*p^2/16  (by Cauchy-Schwarz)
    with equality iff all D_t^2 are equal (i.e., D_t^2 = p/4 for all t).

    This gives E_min when D values are maximally spread (equal magnitudes).
    E_max when D values are maximally concentrated.
    """
    m = (p - 1) // 2

    print(f"\n{'='*70}")
    print(f"ENERGY SPECTRUM STRUCTURE at p={p}")
    print(f"{'='*70}")

    # Cauchy-Schwarz bound
    S2 = m * p / 4  # sum D_t^2
    S4_min = S2**2 / m  # by C-S
    S4_max_theory = S2 * max(abs(d) for d in [0])  # ... need actual bound

    # Compute actual D values for each orientation
    omega = cmath.exp(2j * cmath.pi / p)
    N = 1 << m

    all_D_profiles = {}
    limit = min(N, 64)
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        D_profile = []
        for t in range(1, m + 1):
            lam = sum(omega ** (s * t) for s in S)
            D_profile.append(round(lam.imag, 6))

        D2_profile = tuple(sorted([round(d**2, 4) for d in D_profile]))
        E = additive_energy(S, p)

        if D2_profile not in all_D_profiles:
            all_D_profiles[D2_profile] = (E, bits)

    print(f"\n  Unique D^2 profiles (sorted):")
    for prof in sorted(all_D_profiles.keys()):
        E, bits = all_D_profiles[prof]
        S4 = sum(d for d in prof)  # these are D^2 values, so sum = S2, but we want sum(D^4) = sum(d^2)
        S4_val = sum(d**2 for d in prof)
        print(f"    D^2={prof}  S4={S4_val:.2f}  E={E}")

    # What's special about the energy values?
    E_vals = sorted(set(v[0] for v in all_D_profiles.values()))
    print(f"\n  Energy values: {E_vals}")

    # Check: are they of the form m^2 + k*(p-1)/p or similar?
    if len(E_vals) > 1:
        diffs = [E_vals[i+1] - E_vals[i] for i in range(len(E_vals)-1)]
        print(f"  Energy diffs: {diffs}")


# =============================================================
# MAIN
# =============================================================

print("=" * 70)
print("THM-156: disj3 = (p/4) * E(S) + C(p)")
print("Exact bridge between 3-cycle disjointness and additive energy")
print("=" * 70)

# First: verify the const formula
print("\nConstant C(p) values:")
for p in [7, 11, 13, 17, 19, 23, 29, 31]:
    C = const_C(p)
    m = (p-1)//2
    c3 = p*(p**2-1)//24
    K = -3*p*(p**2-1)*(p**2-9)//320
    print(f"  p={p:2d}: m={m}, c3={c3:>5d}, K={K:>10d}, C={C:>12.2f}")

# Verify THM-156 at each prime
print()
all_pass = True
for p in [7, 11, 13, 17, 19]:
    passed, err = verify_thm156(p, verbose=True)
    if not passed:
        all_pass = False

if all_pass:
    print(f"\n{'*'*70}")
    print(f"*** THM-156 VERIFIED at all tested primes ***")
    print(f"*** disj3(S) = (p/4)*E(S) + C(p) EXACTLY ***")
    print(f"{'*'*70}")

# Energy range analysis
print("\n" + "=" * 70)
print("ENERGY RANGE ANALYSIS")
print("=" * 70)
for p in [7, 11, 13, 17]:
    energy_range_analysis(p)

# THM-155 bridge
for p in [7, 11]:
    thm155_bridge(p)

# Energy spectrum structure
for p in [7, 11, 13]:
    explore_energy_spectrum_structure(p)

print("\nDONE.")
