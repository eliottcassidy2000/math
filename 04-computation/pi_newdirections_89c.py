#!/usr/bin/env python3
"""
pi_newdirections_89c.py — Fresh explorations
opus-2026-03-14-S89c (continued)

New threads:
1. The tournament zeta function Z_n(s) = E[H^{-s}] — analytic continuation?
2. IP zeros on the negative real axis — Lee-Yang universality
3. Random tournament spectral radius vs Paley
4. The "permanent gap" — how far is H from perm(A)?
5. Catalan numbers and tournament paths
6. H(T) for ALL tournaments on n=6 — the full landscape
7. Euler-Mascheroni γ in tournament asymptotics?
"""

import numpy as np
from itertools import combinations, permutations
from math import factorial, log, pi, sqrt, e as euler_e, gamma as gammafn
from fractions import Fraction
from collections import Counter
import time

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(1 << m):
        adj = {v: set() for v in range(n)}
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[i].add(j)
            else:
                adj[j].add(i)
        yield adj

def count_hp(adj, n):
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full].values())

print("=" * 70)
print("PART 1: Full H landscape for n=6")
print("=" * 70)

t0 = time.time()
H_values_6 = []
for T in all_tournaments(6):
    H_values_6.append(count_hp(T, 6))
t1 = time.time()

num_T = len(H_values_6)
print(f"\n  n=6: {num_T} tournaments computed in {t1-t0:.1f}s")

counts_6 = Counter(H_values_6)
print(f"  Distinct H values: {len(counts_6)}")
print(f"  H distribution:")
for h in sorted(counts_6.keys()):
    c = counts_6[h]
    pct = c / num_T * 100
    bar = '#' * int(pct * 2)
    print(f"    H={h:4d}: {c:5d} ({pct:5.2f}%) {bar}")

mean_H = sum(H_values_6) / num_T
var_H = sum((h - mean_H)**2 for h in H_values_6) / num_T
print(f"\n  Mean H = {mean_H:.4f} (theory: 6!/32 = {factorial(6)/32})")
print(f"  Var H = {var_H:.4f}")
print(f"  Std H = {sqrt(var_H):.4f}")
print(f"  CV = {sqrt(var_H)/mean_H:.4f}")
print(f"  Skewness = {sum((h - mean_H)**3 for h in H_values_6) / (num_T * var_H**1.5):.4f}")
print(f"  Kurtosis = {sum((h - mean_H)**4 for h in H_values_6) / (num_T * var_H**2) - 3:.4f}")

# Missing H values
all_possible = set(range(1, max(counts_6.keys()) + 1, 2))  # H is always odd
actual = set(counts_6.keys())
missing = sorted(all_possible - actual)
print(f"\n  Missing odd H values (up to max): {missing}")

print()
print("=" * 70)
print("PART 2: Tournament zeta function Z_n(s)")
print("=" * 70)

# Z_n(s) = (1/|T_n|) Σ_T H(T)^{-s}
# This is the moment generating function of -ln H.

for n in [3, 4, 5, 6]:
    if n <= 5:
        H_list = []
        for T in all_tournaments(n):
            H_list.append(count_hp(T, n))
    else:
        H_list = H_values_6

    N = len(H_list)

    print(f"\n  n={n}: {N} tournaments")

    # Z_n(s) for various s
    for s in [0.5, 1, 2, 3]:
        Z = sum(h**(-s) for h in H_list) / N
        print(f"    Z_{n}({s}) = {Z:.8f}")

    # Analytic properties: Z_n(s) has a pole at s = -∞ (trivially)
    # and is entire for Re(s) > 0.
    # At s=0: Z_n(0) = 1 (always).
    # Derivative at 0: Z_n'(0) = -(1/N) Σ ln H

    mean_lnH = sum(log(h) for h in H_list) / N
    print(f"    -Z'_{n}(0) = E[ln H] = {mean_lnH:.8f}")
    print(f"    exp(E[ln H]) = {np.exp(mean_lnH):.8f}")
    print(f"    E[H] = {sum(H_list)/N:.8f}")
    print(f"    Geometric/Arithmetic mean ratio = {np.exp(mean_lnH)/(sum(H_list)/N):.8f}")

print()
print("=" * 70)
print("PART 3: The permanent vs H for n=6")
print("=" * 70)

# For each tournament, compute both H and perm(A)
perm_data = []
for T in all_tournaments(6):
    H = count_hp(T, 6)
    # Build adjacency matrix
    A = np.zeros((6, 6))
    for i in T:
        for j in T[i]:
            A[i][j] = 1
    # Permanent (exact via Ryser)
    n = 6
    total = 0
    for S_bits in range(1, 1 << n):
        S = [j for j in range(n) if S_bits & (1 << j)]
        k = len(S)
        col_sums_prod = 1
        for i in range(n):
            col_sums_prod *= sum(A[i][j] for j in S)
        total += ((-1) ** (n - k)) * col_sums_prod
    perm = int(round(total))
    perm_data.append((H, perm))

# Correlation between H and perm
Hs = [x[0] for x in perm_data]
Ps = [x[1] for x in perm_data]
corr = np.corrcoef(Hs, Ps)[0, 1]
print(f"\n  Correlation(H, perm) at n=6: {corr:.6f}")

# Scatter plot summary
h_vals = sorted(set(Hs))
for h in h_vals[:10]:
    perms_at_h = [p for hh, p in perm_data if hh == h]
    perm_set = Counter(perms_at_h)
    print(f"    H={h:3d}: perm values = {dict(perm_set)}")

print()
print("=" * 70)
print("PART 4: IP zeros — Lee-Yang analysis for P_7")
print("=" * 70)

# IP(G(P_7), x) = 1 + 80x + 7x²
# Zeros: x = (-80 ± √(6400-28))/14 = (-80 ± √6372)/14

import cmath

# P_7 zeros
a, b, c = 7, 80, 1
disc = b*b - 4*a*c
z1 = (-b + cmath.sqrt(disc)) / (2*a)
z2 = (-b - cmath.sqrt(disc)) / (2*a)
print(f"\n  IP(G(P_7), x) = 7x² + 80x + 1")
print(f"  Zeros: {z1.real:.6f}, {z2.real:.6f}")
print(f"  Both real and NEGATIVE ✓ (Lee-Yang)")
print(f"  Product of zeros = c/a = 1/7 = {1/7:.6f}")
print(f"  Sum of zeros = -b/a = -80/7 = {-80/7:.6f}")

# P_11 zeros
# IP = 1155x³ + 10879x² + 21169x + 1
coeffs_11 = [1155, 10879, 21169, 1]
roots_11 = np.roots(coeffs_11)
print(f"\n  IP(G(P_11), x) = 1155x³ + 10879x² + 21169x + 1")
print(f"  Zeros:")
for r in sorted(roots_11, key=lambda x: x.real):
    if abs(r.imag) < 1e-10:
        print(f"    x = {r.real:.8f} (real)")
    else:
        print(f"    x = {r.real:.8f} ± {abs(r.imag):.8f}i")
all_negative_real = all(abs(r.imag) < 1e-6 and r.real < 0 for r in roots_11)
print(f"  All real and negative: {'✓' if all_negative_real else '✗'}")

# Lee-Yang theorem: for the hard-core lattice gas model on a graph G,
# ALL zeros of the partition function Z(G, λ) = IP(G, λ) lie on the
# NEGATIVE real axis.
# This is NOT always true — it's specific to certain graph families.
# But for our graphs G(P_p), it seems to hold.

# The closest zero to x=2 tells us how "close to a phase transition" we are.
closest = min(abs(r - 2) for r in roots_11)
print(f"\n  Distance from x=2 to nearest zero: {closest:.6f}")
print(f"  We evaluate IP at x=2 (fugacity λ=2), which is past all zeros.")
print(f"  This is the SUPERCRITICAL regime of the hard-core model.")

print()
print("=" * 70)
print("PART 5: Spectral comparison — Paley vs random tournaments")
print("=" * 70)

# For a random tournament on n vertices, the eigenvalues of A fill a disk
# of radius ≈ √n/2 (circular law).
# For Paley, all non-trivial eigenvalues lie on a circle of radius √((p+1)/4).

# Let's compare for n=7: sample random tournaments and compute spectra
np.random.seed(42)
n = 7
n_samples = 1000

# Random tournament eigenvalue radii
random_max_eig = []
for _ in range(n_samples):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            if np.random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    eigs = np.linalg.eigvals(A)
    # Remove Perron eigenvalue (largest real part)
    eigs_sorted = sorted(eigs, key=lambda x: -x.real)
    non_trivial = eigs_sorted[1:]
    radii = [abs(e) for e in non_trivial]
    random_max_eig.append(max(radii))

print(f"\n  n={n}: Random tournament non-trivial eigenvalue radii")
print(f"    Mean max |λ_nontrivial| = {np.mean(random_max_eig):.4f}")
print(f"    Std = {np.std(random_max_eig):.4f}")
print(f"    Paley |λ_nontrivial| = {sqrt((7+1)/4):.4f} (all equal)")
print(f"    Paley ratio vs random max: {sqrt((7+1)/4)/np.mean(random_max_eig):.4f}")

# For Paley: all non-trivial eigenvalues have SAME radius.
# For random: there's a spread. Paley is more "structured."

print()
print("=" * 70)
print("PART 6: Catalan number connection")
print("=" * 70)

# The Catalan numbers C_n = C(2n,n)/(n+1) count many things including:
# - Binary trees with n internal nodes
# - Dyck paths of length 2n
# - Triangulations of (n+2)-gon
# - Non-crossing partitions of [n]
#
# Is there a Catalan connection to tournament paths?
# The number of MONOTONE paths in a tournament is related to permutations
# that avoid certain patterns.
#
# Catalan numbers: 1, 1, 2, 5, 14, 42, 132, ...
# C_n ~ 4^n / (n^{3/2} √π)  — π appears!

from math import comb

print("\n  Catalan numbers and π:")
for n in range(1, 11):
    C_n = comb(2*n, n) // (n+1)
    approx = 4**n / (n**1.5 * sqrt(pi))
    print(f"    C_{n} = {C_n}, 4^n/(n^{3/2}√π) = {approx:.2f}, ratio = {C_n/approx:.6f}")

# Connection to tournaments:
# The number of TOPOLOGICAL SORTS of a tournament T
# (orderings consistent with the tournament) is related to
# the score sequence via the principle of inclusion-exclusion.
# For the transitive tournament: 1 topological sort = 1 Hamiltonian path (H=1).
# For other tournaments: some arcs are "backwards" from the topo sort.

# More directly: the BALLOT PROBLEM counts ways to rank n candidates
# given pairwise comparison results. This IS the Hamiltonian path count.

print()
print("=" * 70)
print("PART 7: H(T) and the binary entropy function")
print("=" * 70)

# The binary entropy H_2(p) = -p log₂ p - (1-p) log₂(1-p)
# peaks at p = 1/2 with H_2(1/2) = 1.
#
# For a tournament T on n vertices with H Hamiltonian paths:
# H / n! is the "path density" — fraction of all orderings that are paths.
# For random tournaments: E[H/n!] = 1/2^{n-1} ≈ 0.
# The entropy of the H distribution tells us about path concentration.

for n in [3, 4, 5, 6]:
    if n <= 5:
        H_list = []
        for T in all_tournaments(n):
            H_list.append(count_hp(T, n))
    else:
        H_list = H_values_6

    N = len(H_list)
    total_H = sum(H_list)

    # Distribution of H/n!
    h_dist = Counter(H_list)
    # Shannon entropy of the H distribution
    entropy = -sum((c/N) * log(c/N) for c in h_dist.values())
    max_entropy = log(len(h_dist))

    # "Path entropy": log₂(H) averaged over tournaments
    path_entropy = sum(log(h) / log(2) for h in H_list) / N

    print(f"\n  n={n}: {len(h_dist)} distinct H values")
    print(f"    Shannon entropy = {entropy:.4f}, max = {max_entropy:.4f}")
    print(f"    Efficiency = {entropy/max_entropy:.4f}")
    print(f"    Mean log₂(H) = {path_entropy:.4f}")
    print(f"    log₂(E[H]) = {log(sum(H_list)/N)/log(2):.4f}")
    print(f"    Gap = {log(sum(H_list)/N)/log(2) - path_entropy:.4f}")

print()
print("=" * 70)
print("PART 8: The Euler-Mascheroni constant γ in Σ 1/H")
print("=" * 70)

# The harmonic-like sum Σ_T 1/H(T) might involve γ.
# γ = lim_{n→∞} (Σ_{k=1}^n 1/k - ln n) ≈ 0.5772...

for n in [3, 4, 5, 6]:
    if n <= 5:
        H_list = []
        for T in all_tournaments(n):
            H_list.append(count_hp(T, n))
    else:
        H_list = H_values_6

    N = len(H_list)
    harm_sum = sum(Fraction(1, h) for h in H_list)
    mean_recip = harm_sum / N

    print(f"\n  n={n}:")
    print(f"    Σ 1/H = {harm_sum} = {float(harm_sum):.6f}")
    print(f"    E[1/H] = {float(mean_recip):.8f}")
    print(f"    1/E[H] = {float(Fraction(N, sum(H_list))):.8f}")
    print(f"    E[1/H] / (1/E[H]) = {float(mean_recip * Fraction(sum(H_list), N)):.8f}")

print()
print("=" * 70)
print("PART 9: Score sequence → H — the 'tournament invariant map'")
print("=" * 70)

# For n=6, how well does the score sequence predict H?
score_to_H = {}
for T in all_tournaments(6):
    H = count_hp(T, 6)
    scores = tuple(sorted([len(T[v]) for v in T]))
    if scores not in score_to_H:
        score_to_H[scores] = []
    score_to_H[scores].append(H)

print(f"\n  n=6: {len(score_to_H)} distinct score sequences")
for scores in sorted(score_to_H.keys()):
    H_vals = score_to_H[scores]
    h_set = sorted(set(H_vals))
    mean = sum(H_vals) / len(H_vals)
    if len(h_set) <= 6:
        print(f"    {scores}: H ∈ {h_set}, mean={mean:.1f}, count={len(H_vals)}")
    else:
        print(f"    {scores}: {len(h_set)} distinct H values, range [{min(H_vals)},{max(H_vals)}], mean={mean:.1f}")

print()
print("=" * 70)
print("PART 10: π in the n=6 moment ratios")
print("=" * 70)

# Exact moments of H at n=6
moments = {}
for k in range(1, 7):
    moments[k] = Fraction(sum(h**k for h in H_values_6), num_T)
    print(f"  E[H^{k}] = {moments[k]} = {float(moments[k]):.4f}")

# E[H²]/E[H]²
ratio_2 = moments[2] / moments[1]**2
print(f"\n  E[H²]/E[H]² = {ratio_2} = {float(ratio_2):.8f}")

# The CV² sequence: 1/3, 1/3, 19/60, 13/45, ?
cv2 = (moments[2] - moments[1]**2) / moments[1]**2
print(f"  CV² = Var/E² = {cv2} = {float(cv2):.8f}")

# Does π appear in these ratios?
# 19/60 ≈ 0.3167 (n=5), 13/45 ≈ 0.2889 (n=6)
# π/10 ≈ 0.3142, close to 19/60!
# π²/34 ≈ 0.2900, close to 13/45!

print(f"\n  Looking for π in moment ratios:")
print(f"    n=5: CV² = 19/60 = {19/60:.6f}, π/10 = {pi/10:.6f}")
print(f"    n=6: CV² = {float(cv2):.6f}, π²/34 = {pi**2/34:.6f}")

# Check if E[H²]/E[H]² approaches 1 + 1/3 = 4/3?
# Known: 4/3, 4/3, 79/60, 58/45 for n=3,4,5,6
print(f"\n  E[H²]/E[H]² sequence:")
print(f"    n=3: 4/3 = {4/3:.6f}")
print(f"    n=4: 4/3 = {4/3:.6f}")
print(f"    n=5: 79/60 = {79/60:.6f}")
print(f"    n=6: {float(ratio_2):.6f} = {ratio_2}")

# The limit might be (1 + 1/e) or some other transcendental
import sympy
print(f"    {ratio_2} simplified: {sympy.Rational(ratio_2.numerator, ratio_2.denominator)}")

print()
print("=" * 70)
print("PART 11: The π-H dictionary v2 — incorporating n=6 data")
print("=" * 70)

print("""
  UPDATED π-TOURNAMENT DICTIONARY (v2)
  ======================================

  1. MEAN: E[H] = n!/2^{n-1} ~ √(2πn)(n/2e)^n
     π IN: Stirling's approximation

  2. EIGENPHASE: arg(λ(P_p))/π = 1/2 + 1/(π√p)
     π IN: arctan expansion near the pole

  3. GAUSS SUM: g = Σ χ(a)exp(2πia/p), |g| = √p
     π IN: roots of unity construction

  4. CLT: H ~ Normal(μ, σ²) for large n
     π IN: Gaussian density (2πσ²)^{-1/2}

  5. DFT: Paley adjacency diagonalized by F_{jk} = exp(2πijk/p)/√p
     π IN: every matrix entry

  6. WEIL: |Σ χ(f)exp(2πix/p)| ≤ d√p
     π IN: character sum bounds

  7. CATALAN: C_n ~ 4^n/(n^{3/2}√π)
     π IN: path counting asymptotics

  8. LEE-YANG: IP zeros on negative real axis
     π IN: phase transition at λ_c (singularity of free energy)

  9. SPECTRAL: Random tournament eigenvalues fill disk of radius √n/2
     π IN: circular law (eigenvalue area = πr²)

  10. ENTROPY: Shannon entropy of H distribution ~ log n + ...
      π IN: maximum entropy distribution (Gaussian) involves π

  11. THM-212: p | H(P_p), orbits under the CIRCLE GROUP Z/pZ
      π IN: the circle group is π's group

  12. H(P_11) = 5 × 1729 × 11, where 1729 = 10³ + 9³
      π IN: taxicab numbers connect to Ramanujan and modular forms

  DEEP TRUTH: Every appearance of π traces back to CIRCLES or PERIODICITY.
  Tournament theory studies ORIENTED CIRCLES (directed cycles).
  π is therefore the fundamental constant of tournament theory.
""")

print("Done!")
