#!/usr/bin/env python3
"""
analytic_crossover.py — Analytic prediction of the Paley→Interval crossover

Combining THM-141/142 (kind-pasteur S59) with the Walsh difference formula (S65):

H(Int) - H(Pal) = 2·Σ_{NQR products} Ĥ[S]
                 = sum of degree-2 + degree-4 + ... contributions

At the LOWEST ORDER (OCF with only 3-cycles):
  H ≈ 1 + 2α₁ + 4α₂ + 8α₃ + ...

The key quantities:
  α₁(Pal) - α₁(Int) = number of extra odd cycles in Paley
  α₂(Int) - α₂(Pal) ≥ Δ_3,3 = p(p-1)(p+1)(p-3)/192  [THM-142, for 3-3 pairs only]

The crossover occurs when:
  2·Δα₁ < 4·Δα₂ + 8·Δα₃ + ...
  i.e., 2·(α₁(P)-α₁(I)) < 4·(α₂(I)-α₂(P)) + higher order

Can we predict the critical p from these?

Author: opus-2026-03-12-S67
"""

import numpy as np
from itertools import product as iprod

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def tournament_adjacency(sigma, p):
    m = (p-1)//2
    n = p
    A = np.zeros((n, n), dtype=int)
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

print("=" * 70)
print("ANALYTIC CROSSOVER PREDICTION")
print("=" * 70)

# THM-142: exact 3-3 disjointness excess
print("\n--- THM-142: 3-3 Disjointness Excess ---")
for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    if p % 4 == 3:
        m = (p-1)//2
        delta_33 = p * (p-1) * (p+1) * (p-3) // 192
        # c3 for regular tournament on p vertices
        c3 = p * (p*p - 1) // 24
        # α₁ (total odd cycles) is harder — but c3 is the 3-cycle count
        print(f"  p={p:>2}: Δ(3,3 disjoint) = {delta_33:>10}, "
              f"c3 = {c3:>6}, "
              f"4·Δ/c3 = {4*delta_33/c3:.2f}")

# The ratio 4·Δ_{3,3}/c₃ measures how much the 3-cycle disjointness
# advantage contributes to H relative to the total 3-cycle count.
# But we also need Δα₁.

print("\n--- Known Δα₁ data ---")
# From THM-138:
data = {
    7: {"alpha1_P": 80, "alpha1_I": 59, "alpha2_P": 7, "alpha2_I": 14, "H_P": 189, "H_I": 175},
    11: {"alpha1_P": 21169, "alpha1_I": 18397, "alpha2_P": 10879, "alpha2_I": 11110,
         "H_P": 95095, "H_I": 93027},
}
# At p=19: only alpha1 known
data[19] = {"alpha1_P": 130965270477, "alpha1_I": 126443605257,
            "H_P": 1172695746915, "H_I": 1184212824763}

for p in sorted(data.keys()):
    d = data[p]
    da1 = d["alpha1_P"] - d["alpha1_I"]
    dH = d.get("H_P", 0) - d.get("H_I", 0)
    print(f"  p={p}: Δα₁ = {da1:>15}, 2·Δα₁ = {2*da1:>15}, H(P)-H(I) = {dH:>15}")
    if "alpha2_P" in d:
        da2 = d["alpha2_I"] - d["alpha2_P"]
        print(f"         Δα₂ = {da2:>15}, 4·Δα₂ = {4*da2:>15}")
        # Residual from higher order
        residual = dH - 2*da1 + 4*da2
        print(f"         residual (8·Δα₃+...) = {residual}")

# The crossover formula:
# H(P) - H(I) = 2·Δα₁ - 4·Δα₂ - 8·Δα₃ - ...
# For Paley to win: 2·Δα₁ > 4·Δα₂ + 8·Δα₃ + ...

print("\n--- Growth rates ---")
# Δα₁ grows as some power of p. Let's fit.
ps = [7, 11, 19]
da1s = [80-59, 21169-18397, 130965270477-126443605257]
for i, p in enumerate(ps):
    print(f"  p={p}: Δα₁ = {da1s[i]}, log Δα₁ / log p = {np.log(da1s[i])/np.log(p):.2f}")

# The THM-142 disjointness excess Δ_{3,3} = p(p-1)(p+1)(p-3)/192 ~ p⁴/192
# If Δα₂ ≈ Δ_{3,3} (3-3 pairs dominate), then:
#   4·Δα₂ ≈ 4p⁴/192 = p⁴/48

# For the crossover: 2·Δα₁ = p⁴/48
# Need to know the growth rate of Δα₁.

# From the trace formula:
# α₁(P) = Σ_k c_k(P) and α₁(I) = Σ_k c_k(I)
# Δα₁ = Σ_k [c_k(P) - c_k(I)]
# Each c_k is bounded by the trace formula: c_k = (1/p)Σ λ_j^k
# For Paley: λ_j = ±√((p+1)/4) → |λ_j|^k ~ ((p+1)/4)^{k/2}
# For Interval: dominant eigenvalue |μ₁| ≈ m·2/π ≈ p/π

# At length k=3: c₃ is the SAME for both (both regular)
# At length k=5: c₅(P) > c₅(I) (Paley has ~10% more)
# The bulk of Δα₁ comes from MIDDLE length cycles

# Rough scaling: Δα₁ ~ Σ_k [c_k(P) - c_k(I)]
# c_k(P) - c_k(I) ~ (1/p)[(p-1)(√p/2)^k - 2(p/π)^k cos(kθ) + ...]
# For moderate k: Paley term ~ (p-1)(p/4)^{k/2}/p grows exponentially
# Interval term: |μ₁|^k/p ~ (p/π)^k/p also grows, but |μ₁|^k >> (√p)^k for large k

# Dominant contribution to Δα₁ is from short-to-medium cycles where Paley leads.

print("\n--- Trace formula analysis of Δα₁ ---")
for p in [7, 11, 19]:
    m = (p-1)//2
    # Paley eigenvalues: all have |λ| = √((p+1)/4) (from Weil bound)
    pal_mag = np.sqrt((p+1)/4)
    # Interval eigenvalues from DFT
    omega = np.exp(2j * np.pi / p)
    int_eigs = []
    for k in range(1, p):
        lam = sum(omega**(j*k) for j in range(1, m+1))
        int_eigs.append(abs(lam))
    int_top = max(int_eigs)

    print(f"\n  p={p}: |λ_Pal| = {pal_mag:.3f}, |μ₁_Int| = {int_top:.3f}, ratio = {int_top/pal_mag:.3f}")

    # Trace comparison for each odd k
    for k in [3, 5, 7, 9]:
        if k > p:
            break
        # Paley: tr(A^k)/p = Σ λ_j^k / p
        # For Paley: all |λ| = pal_mag, but phases matter
        # Actually for Paley circulant: λ_j = Σ_{a∈QR} ω^{ja}
        # |λ_j| = √p for p≡3 mod 4 (Gauss sum magnitude)
        qr = sorted(a for a in range(1, p) if legendre(a, p) == 1)
        if p % 4 == 3:
            pal_eigs = []
            for j in range(1, p):
                lam = sum(omega**(j*a) for a in qr)
                pal_eigs.append(lam)

            trace_pal = sum(e**k for e in pal_eigs).real / p
            trace_int = sum(complex(sum(omega**(j*s) for s in range(1, m+1)))**k
                           for j in range(1, p)).real / p
            diff = trace_pal - trace_int
            print(f"    k={k}: tr_P/p={trace_pal:.1f}, tr_I/p={trace_int:.1f}, Δ={diff:.1f}")

# =============================================================
# THE CROSSOVER ESTIMATE
# =============================================================
print("\n" + "=" * 70)
print("CROSSOVER ESTIMATE")
print("=" * 70)

# From the data:
# p=7:  H(P)-H(I) = +14  (Paley wins)
# p=11: H(P)-H(I) = +2068  (Paley wins, barely)
# p=19: H(P)-H(I) = -11.5B  (Interval wins)

# The ratio H(P)/H(I):
for p, hp, hi in [(7, 189, 175), (11, 95095, 93027), (19, 1172695746915, 1184212824763)]:
    ratio = hp / hi
    print(f"  p={p}: H(P)/H(I) = {ratio:.6f}, log ratio = {np.log(ratio):.6f}")

# Linear interpolation in log(H_P/H_I) vs 1/p
# p=7: log ratio = +0.077
# p=11: log ratio = +0.022
# p=19: log ratio = -0.0098
# The crossover (log ratio = 0) occurs at approximately...

from scipy.optimize import brentq

ps_data = [7, 11, 19]
log_ratios = [np.log(189/175), np.log(95095/93027), np.log(1172695746915/1184212824763)]

# Fit quadratic in 1/p
x = [1/p for p in ps_data]
y = log_ratios
coeffs = np.polyfit(x, y, 2)
poly = np.poly1d(coeffs)

# Find zero
try:
    p_c_inv = brentq(poly, 0.01, 0.2)
    p_c = 1 / p_c_inv
    print(f"\n  Quadratic fit: log(H_P/H_I) = {coeffs[0]:.4f}/p² + {coeffs[1]:.4f}/p + {coeffs[2]:.6f}")
    print(f"  Predicted crossover: p_c ≈ {p_c:.1f}")
except:
    print("  Could not find crossover via quadratic fit")

# Linear fit
coeffs_lin = np.polyfit(x[:2] + [x[2]], y[:2] + [y[2]], 1)
print(f"  Linear fit: log(H_P/H_I) = {coeffs_lin[0]:.4f}/p + {coeffs_lin[1]:.6f}")
if coeffs_lin[0] != 0:
    p_c_lin = -coeffs_lin[1] / coeffs_lin[0]
    if 0 < p_c_lin < 1:
        print(f"  Linear predicted crossover: 1/p_c = {p_c_lin:.4f}, p_c ≈ {1/p_c_lin:.1f}")

# =============================================================
# The scaling argument
# =============================================================
print("\n" + "=" * 70)
print("SCALING ARGUMENT")
print("=" * 70)
print("""
From the data and THM-142:

  Δα₁ = α₁(P) - α₁(I) grows as ~ p^{a} for some a
  Δα₂ ≥ p(p-1)(p+1)(p-3)/192 ~ p⁴/192  [3-3 pairs only, THM-142]

  H(P) - H(I) = 2·Δα₁ - 4·Δα₂ - 8·Δα₃ - ...

The crossover is when the POSITIVE 2·Δα₁ term is overcome by the
NEGATIVE -4·Δα₂ and higher terms.

From the trace formula at p=19:
  Δα₁ ≈ 4.5 × 10⁹, grows roughly as (p/e)^{p/2}
  But so does α₁ itself, and the RATIO Δα₁/α₁ → 0 as p → ∞

In contrast, Δα₂/α₂ stays bounded away from 0 (the disjointness
gap is STRUCTURAL, not just a fluctuation).

The competition between exponentially growing terms with different
bases is what determines the crossover.

For p ≡ 3 mod 4 primes: crossover between p=11 and p=19.
The exact critical prime is BETWEEN 11 and 19 — likely p=13 or p=17
(but these are ≡ 1 mod 4, where Paley doesn't exist!).

This suggests: for p ≡ 3 mod 4, the crossover occurs at the first
prime p > 11, which is p = 19. The transition is "sharp" — there
is no intermediate prime to test.

For p ≡ 1 mod 4: Interval always wins (no Paley eigenvector advantage).
""")

# =============================================================
# NEW: Can we prove Interval wins for ALL p ≥ 19 (p ≡ 3 mod 4)?
# =============================================================
print("=" * 70)
print("TOWARDS A PROOF: INTERVAL WINS FOR p ≥ 19")
print("=" * 70)
print("""
STRATEGY:

1. At degree 2 (J matrix):
   Paley eigenvalue λ_P = +7 (p=7), +561 (p=11), -544M (p=19)
   λ_P changes sign → degree-2 ALSO favors Interval at p=19!

   If λ_P < 0 for all p ≥ 19 (p ≡ 3 mod 4), then EVERY Walsh degree
   favors Interval, and the proof is immediate.

2. Can we compute λ_P analytically?
   λ_P = (1/m) σ_P^T J σ_P

   From the Walsh decomposition:
   J[i,j] = Ĥ[{i,j}] = interaction between chords i and j

   σ_P^T J σ_P = Σ_{i,j} J[i,j] χ(i) χ(j)
   = Σ_{QR pairs} J[i,j] + Σ_{NQR pairs} J[i,j] - 2·Σ_{mixed pairs} J[i,j]
   = (Σ_all J[i,j]) - 4·Σ_{mixed} J[i,j]

   At p=13: Σ_{mixed} = 0 exactly (mixed pairs cancel)
   At p=7,11: Σ_{mixed} < 0 (all negative)
   At p=19: Σ_{mixed} > 0 (from HYP-512 data, λ_P < 0)

   So the question is: when does Σ_{mixed} J[i,j] change sign?

3. The J matrix entries J[i,j] involve the FULL Hamiltonian path counting
   function, which is hard to compute analytically. But the LINEARIZED
   version (from cycle counts only) is computable!

   From S66: the linearized Ĥ is COMPLETELY WRONG (nonlinear corrections
   dominate by 100:1 ratio). So cycle counts alone don't determine J.

4. ALTERNATIVE APPROACH: Use the SPECTRAL concentration bound.
   IPR(Interval) → 1/3 as p → ∞
   IPR(Paley) → 0 as p → ∞

   If IPR determines H ranking for large p, and we can prove this
   rigorously, then Interval wins for all sufficiently large p.

   Combined with exhaustive verification at p=19 (done), this gives
   Interval wins for all p ≥ 19 with p ≡ 3 mod 4.
""")

print("DONE.")
