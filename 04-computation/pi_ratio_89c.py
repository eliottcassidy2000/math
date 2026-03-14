#!/usr/bin/env python3
"""
pi_ratio_89c.py — Is H(P_p)/E[H] → e?
opus-2026-03-14-S89c

Data: 2.0, 2.4, 2.44, 2.53, 2.56 for p=3,7,11,19,23
ln: 0.693, 0.875, 0.892, 0.927, 0.939

If limit is e (=2.718...), ln(limit) = 1.
The ln values approach 1 from below. What's the rate?

Also: compare with doubly-regular tournament heuristic.
A doubly regular tournament on p vertices has exactly (p-3)/4
common out-neighbors for each pair. This extreme regularity
might give H in terms of a "doubly regular permanent".
"""

from math import factorial, log, pi, sqrt, e as euler_e, lgamma
from fractions import Fraction

print("=" * 70)
print("PART 1: Is ln(H/E[H]) → 1?")
print("=" * 70)

data = [
    (3, 3),
    (7, 189),
    (11, 95095),
    (19, 1172695746915),
    (23, 15760206976379349),
]

print(f"\n  {'p':>3} | {'ln(H/E[H])':>12} | {'1 - ln(H/E[H])':>14} | {'gap × √p':>10}")
print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}")

for p, H in data:
    mean = factorial(p) / 2**(p-1)
    ratio = H / mean
    lr = log(ratio)
    gap = 1 - lr
    print(f"  {p:3d} | {lr:12.8f} | {gap:14.8f} | {gap * sqrt(p):10.6f}")

# If gap ~ C/√p, then ln(H/E[H]) = 1 - C/√p + ...
# and H/E[H] = e × exp(-C/√p) → e.
# Let's fit C:
print(f"\n  Fitting gap ≈ C / √p:")
for i, (p, H) in enumerate(data):
    mean = factorial(p) / 2**(p-1)
    lr = log(H / mean)
    gap = 1 - lr
    C = gap * sqrt(p)
    print(f"    p={p}: C = {C:.6f}")

# C values: 0.531, 0.329, 0.359, 0.318, 0.294
# Not perfectly constant, still varying. Maybe gap ~ C/p^α?
print(f"\n  Fitting gap ≈ C / p^α (linear regression in log-log):")
import numpy as np
ps = [p for p, _ in data]
gaps = []
for p, H in data:
    mean = factorial(p) / 2**(p-1)
    lr = log(H / mean)
    gaps.append(1 - lr)

log_ps = np.log(ps)
log_gaps = np.log(gaps)
# Fit: log_gap = log_C - α × log_p
alpha, log_C = np.polyfit(log_ps, log_gaps, 1)
C = np.exp(log_C)
print(f"    α = {-alpha:.4f}, C = {C:.4f}")
print(f"    gap ≈ {C:.3f} / p^{-alpha:.3f}")

# Predicted at p=31 (next Paley prime ≡ 3 mod 4):
gap_31 = C * 31**alpha
print(f"\n    Predicted gap at p=31: {gap_31:.6f}")
print(f"    Predicted ln(H/E[H]) at p=31: {1 - gap_31:.6f}")
print(f"    Predicted H/E[H] at p=31: {np.exp(1 - gap_31):.6f}")

# Is the limit actually e?
print(f"\n  If limit is e = {euler_e:.6f}:")
print(f"  Current values: {[f'{H/factorial(p)*2**(p-1):.4f}' for p, H in data]}")
print(f"  Need to reach: {euler_e:.4f}")
print(f"  Gap from e: {[f'{euler_e - H/factorial(p)*2**(p-1):.4f}' for p, H in data]}")

print()
print("=" * 70)
print("PART 2: Alternative limit candidates")
print("=" * 70)

# What if the limit is not e but something else?
# Candidates: e, 8/3 ≈ 2.667, e^{π/(2e)} ≈ ?, sqrt(2π/e) ≈ ?
import math

candidates = [
    ("e", euler_e),
    ("8/3", 8/3),
    ("π - 1/2", pi - 0.5),
    ("√(2π)", sqrt(2*pi)),
    ("e^{1-1/e}", euler_e**(1-1/euler_e)),
    ("3 - 1/e", 3 - 1/euler_e),
    ("5/2", 2.5),
    ("1 + √2", 1 + sqrt(2)),
    ("2 + 1/π", 2 + 1/pi),
    ("e × (1 - 1/(2e))", euler_e * (1 - 1/(2*euler_e))),
]

last_ratio = data[-1][1] / (factorial(data[-1][0]) / 2**(data[-1][0]-1))
print(f"\n  Last computed ratio (p=23): {last_ratio:.6f}")
print(f"\n  {'Candidate':>25} | {'Value':>10} | {'Distance':>10}")
print(f"  {'-'*25}-+-{'-'*10}-+-{'-'*10}")
for name, val in sorted(candidates, key=lambda x: abs(x[1] - last_ratio)):
    dist = val - last_ratio
    print(f"  {name:>25} | {val:10.6f} | {dist:+10.6f}")

print()
print("=" * 70)
print("PART 3: Doubly regular structure and Hamiltonian paths")
print("=" * 70)

# For a doubly regular tournament on p vertices:
# - Each vertex has outdegree (p-1)/2
# - Each pair of vertices has exactly (p-3)/4 common out-neighbors
#
# This means: starting from any vertex v, we can go to (p-1)/2 neighbors.
# From each neighbor u, we can go to (p-3)/4 common out-neighbors of u
# (minus vertices already visited).
#
# For a HEURISTIC estimate:
# Step 1: p choices for first vertex
# Step 2: (p-1)/2 choices for second vertex
# Step 3: from the perspective of the current vertex u, about
#   (p-1)/2 - 1 remaining out-neighbors, minus already visited
# But the doubly regular property means the correlation between
# successive steps is controlled.

# The "mean-field" estimate:
# H ≈ p × (p-1)/2 × ((p-3)/2)/(p-2) × ((p-5)/2)/(p-3) × ...
# This is the "greedy" estimate assuming at each step we have
# outdegree × (remaining-1)/(total-1) choices.

# Actually for random d-regular tournament (d = (p-1)/2):
# Mean H = p! × (d/p-1))^{p-1} × correction
# Hmm, this isn't right either.

# For a doubly regular tournament:
# The permanent of the adjacency matrix relates to something...

# Let's compute a different quantity: the "expected next step" approach
for p in [3, 7, 11, 19, 23]:
    d = (p-1)//2  # out-degree
    lambda_common = (p-3)//4 if p >= 7 else 0  # common out-neighbors

    # Greedy estimate: at step k (having visited k vertices),
    # expected out-neighbors among unvisited ≈ d × (p-k)/(p-1)
    # So H_greedy ≈ p × Π_{k=1}^{p-1} d×(p-k)/(p-1)
    # = p × (d/(p-1))^{p-1} × Π_{k=1}^{p-1} (p-k)
    # = p × (d/(p-1))^{p-1} × (p-1)!
    # = p! × (d/(p-1))^{p-1}

    # But d = (p-1)/2, so d/(p-1) = 1/2
    # H_greedy = p! × (1/2)^{p-1} = p!/2^{p-1} = E[H]! Just the mean!

    # So the greedy estimate IS the mean. The doubly regular property
    # must give an EXCESS over the mean.

    # The excess might come from the positive correlation between steps.
    # In a doubly regular tournament, knowing that v→u tells you
    # exactly how many of u's out-neighbors are also out-neighbors of v.

    # The number of common out-neighbors is λ = (p-3)/4.
    # The expected number for random: d²/(p-1) = (p-1)²/(4(p-1)) = (p-1)/4
    # So λ = (p-3)/4 < (p-1)/4. Doubly regular has FEWER common
    # out-neighbors than random! (By 1/2.)

    # This means the directions "spread out" more uniformly,
    # which should INCREASE the number of Hamiltonian paths.

    H = dict(data).get(p, 0)
    mean = factorial(p) / 2**(p-1)

    print(f"\n  p={p}: d=(p-1)/2={d}, λ=(p-3)/4={(p-3)/4:.1f}")
    print(f"    Random λ_expected = (p-1)/4 = {(p-1)/4:.1f}")
    print(f"    Deficit = {(p-1)/4 - (p-3)/4:.1f}")
    print(f"    H/E[H] = {H/mean:.6f}" if H > 0 else "    H not computed")

print()
print("=" * 70)
print("PART 4: The skew-adjacency matrix and Pfaffian")
print("=" * 70)

# The skew-adjacency matrix S of a tournament T:
# S_{ij} = 1 if i→j, -1 if j→i, 0 if i=j
# S is skew-symmetric: S^T = -S
#
# For even n: det(S) = Pf(S)² where Pf is the Pfaffian
# For odd n: det(S) = 0 (skew-symmetric matrices have zero det when n is odd)
#
# All our Paley primes p are odd, so det(S(P_p)) = 0 trivially.
# But det(S) for the n-1 minor might be nonzero.

# For the Paley tournament, S = A - A^T = 2A - J + I
# where J is all-ones and A is adjacency.
# Eigenvalues of S: for the trivial eigenvalue direction (all-ones):
# S × 1 = (A - A^T) × 1 = d×1 - (p-1-d)×1 = (2d-p+1)×1 = 0
# So 0 is always an eigenvalue with eigenvector 1.
# Other eigenvalues: λ_k(S) = λ_k(A) - λ_k(A)^* = 2i×Im(λ_k(A))
# For Paley: λ_k(A) = (-1 ± i√p)/2
# So λ_k(S) = ±i√p

# Therefore S has eigenvalues: 0 (once), i√p ((p-1)/2 times), -i√p ((p-1)/2 times)
# |det(S)| = 0 (because of the zero eigenvalue)

# But the (p-1)×(p-1) principal minor:
# det(S_{00}) = product of nonzero eigenvalues / ???
# Actually the cofactor matrix might be more interesting.

import numpy as np

for p in [3, 7, 11, 19]:
    # Build skew-adjacency
    S = np.zeros((p, p))
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    for i in range(p):
        for j in range(p):
            if i == j:
                continue
            if (j - i) % p in qr:
                S[i][j] = 1
            else:
                S[i][j] = -1

    # Eigenvalues
    eigs = np.linalg.eigvals(S)
    # Should be 0, ±i√p
    print(f"\n  P_{p}: S eigenvalues (sorted by Im):")
    eigs_sorted = sorted(eigs, key=lambda x: x.imag)
    for e in eigs_sorted[:3]:
        print(f"    {e.real:+.4f} + {e.imag:+.4f}i")
    print(f"    ... (then symmetric)")

    # Minor: delete row 0, col 0
    S_minor = S[1:, 1:]
    det_minor = np.linalg.det(S_minor)
    print(f"    det(S[1:,1:]) = {det_minor:.1f}")

    # det of (p-1)×(p-1) skew matrix: if p-1 is even, det = Pf²
    if (p-1) % 2 == 0:
        # Pfaffian
        pf = round(sqrt(abs(det_minor)))
        print(f"    Pfaffian ≈ ±{pf}")
        print(f"    p^{(p-1)//2} = {p**((p-1)//2)}")
        print(f"    |det|/p^{(p-1)//2} = {abs(det_minor)/p**((p-1)//2):.6f}")

print()
print("=" * 70)
print("PART 5: H(P_p) mod various quadratic residues")
print("=" * 70)

# For each p, check H mod q for all QR q of p
for p, H in data:
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    qr_sorted = sorted(qr)
    nqr = sorted(set(range(1, p)) - qr)

    print(f"\n  p={p}: QR = {qr_sorted[:10]}{'...' if len(qr_sorted) > 10 else ''}")
    print(f"    H mod QR: {[(q, H % q) for q in qr_sorted[:8]]}")
    print(f"    H mod NQR: {[(q, H % q) for q in nqr[:8]]}")

    # How many QR divide H?
    qr_divides = [q for q in qr_sorted if H % q == 0]
    nqr_divides = [q for q in nqr if H % q == 0]
    print(f"    QR that divide H: {qr_divides}")
    print(f"    NQR that divide H: {nqr_divides}")

print()
print("=" * 70)
print("PART 6: The path-count generating function over primes")
print("=" * 70)

# Dirichlet series: L_H(s) = Σ_p H(P_p) / p^s
# where the sum is over Paley primes p ≡ 3 mod 4

print("\n  Partial sums of Σ H(P_p) / p^s:")
for s in [1, 2, 3]:
    total = sum(H / p**s for p, H in data)
    print(f"    s={s}: Σ = {total:.6e}")

# The series diverges for s=1 (since H ~ p!/2^p and p^s is polynomial).
# For s > p, it converges trivially.
# More interesting: Σ ln(H(P_p)) / p^s
print(f"\n  Partial sums of Σ ln(H(P_p)) / p^s:")
for s in [1, 2]:
    total = sum(log(H) / p**s for p, H in data)
    print(f"    s={s}: Σ = {total:.6f}")

# Even more interesting: Σ (H/E[H]) / p^s
print(f"\n  Partial sums of Σ (H(P_p)/E[H]) / p^s:")
for s in [1, 2]:
    total = sum((H / (factorial(p) / 2**(p-1))) / p**s for p, H in data)
    print(f"    s={s}: Σ = {total:.6f}")

print()
print("=" * 70)
print("PART 7: Connection to random matrix theory")
print("=" * 70)

# The Paley tournament adjacency matrix has spectrum
# {(p-1)/2, (-1±i√p)/2 each (p-1)/2 times}
# The non-trivial eigenvalues lie on a circle of radius √((p+1)/4)
# This is like a GUE matrix where eigenvalues lie on a curve.
#
# In random matrix theory, the circular law says eigenvalues of
# random non-Hermitian matrices fill a disk.
# But Paley eigenvalues all lie on a SINGLE CIRCLE — perfectly structured.
#
# The phase distribution is uniform (each root of unity gives one eigenvalue).
# This is the most ORDERED possible spectrum — antithesis of random.

print("  Paley adjacency spectrum:")
print("    - Trivial: λ₀ = (p-1)/2 (the Perron eigenvalue)")
print("    - Non-trivial: on circle |λ| = √((p+1)/4)")
print("    - Phases: 2πk/p for k = 1,...,p-1, mapped through χ_p")
print()
print("  This is maximally STRUCTURED — all eigenvalues on one circle.")
print("  Random tournament eigenvalues fill a DISK of radius ~√n/2.")
print("  Paley concentrates all spectral weight on the boundary circle.")
print()
print("  The spectral radius √((p+1)/4) ≈ √p/2 matches the random")
print("  tournament spectral radius of ~√n/2 — same scale, different geometry.")

# For random tournaments: eigenvalues fill disk of radius ≈ √n/2
# by the circular law.
# Paley: eigenvalues on circle of radius √(p+1)/4 ≈ √p/2.
# Same RADIUS, but Paley has all eigenvalues on the boundary!

# This concentration on the boundary means maximum "spectral energy"
# at each eigenvalue, which relates to path counts.

print()
print("=" * 70)
print("DONE — The π exploration continues...")
print("=" * 70)
