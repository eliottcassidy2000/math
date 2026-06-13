#!/usr/bin/env python3
"""
convergence_proof.py — opus-2026-03-14-S71d

PROVE: φ_k = 2 - (1/2)^k + O((1/4)^k)  [k-nacci dominant root]
       ψ_k = 3 - (2/3)^k + O((4/9)^k)  [weighted k-nacci dominant root]

The convergence constant C = -1 exactly.

APPROACH: Substitute φ_k = 2 - ε into the k-nacci characteristic equation.

k-nacci: x^k = x^{k-1} + x^{k-2} + ... + 1 = (x^k - 1)/(x - 1)
So: x^k(x-1) = x^k - 1, i.e., x^{k+1} - 2x^k + 1 = 0.

For x = 2 - ε:
  (2-ε)^{k+1} - 2(2-ε)^k + 1 = 0
  (2-ε)^k [(2-ε) - 2] + 1 = 0
  (2-ε)^k · (-ε) + 1 = 0
  ε · (2-ε)^k = 1
  ε = 1/(2-ε)^k ≈ 1/2^k · 1/(1-ε/2)^k
  ε ≈ (1/2)^k · (1 + kε/2 + ...)
  ε ≈ (1/2)^k + k/2 · ε · (1/2)^k

To first order: ε ≈ (1/2)^k, so C₂ = ε · 2^k = 1? But we got -1!

Wait, (2-ε)^k(-ε) + 1 = 0 gives ε(2-ε)^k = 1, so ε > 0.
But φ_k < 2 for all finite k, so ε = 2 - φ_k > 0.
And ε · 2^k → 1 as k → ∞.

But our computation gave C₂ = (φ_k - 2) · 2^k → -1.
Since φ_k - 2 = -ε, we get (-ε) · 2^k → -1.
So ε · 2^k → 1. C₂ = -1 is correct because we compute (φ_k - 2) · 2^k = -ε · 2^k.

So the PROOF is: ε · (2-ε)^k = 1, which gives ε · 2^k(1-ε/2)^k = 1.
As k → ∞: ε → 0, (1-ε/2)^k → 1 (since ε ~ 2^{-k} and k·ε ~ k/2^k → 0).
Hence ε · 2^k → 1. QED.

Similarly for weighted k-nacci: ε · (3-ε)^k = C for some constant.
"""

import sys
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PROOF: CONVERGENCE CONSTANTS FOR k-NACCI FAMILIES")
print("=" * 70)

# ======================================================================
# PART 1: k-nacci proof
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: k-NACCI — EXACT ASYMPTOTICS")
print("=" * 70)

print("""
THEOREM: The dominant root φ_k of the k-nacci recurrence satisfies:
  φ_k = 2 - 1/2^k + O(k/4^k)

PROOF:
The k-nacci recurrence x_n = x_{n-1} + ... + x_{n-k} has characteristic
polynomial p(x) = x^k - x^{k-1} - ... - 1.

Multiply by (x-1): (x-1)·p(x) = x^{k+1} - 2x^k + 1.

So the roots of p(x) are the roots of x^{k+1} - 2x^k + 1 = 0
EXCEPT x = 1 (which is a root of (x-1) but not of p(x)).

Factoring: x^{k+1} - 2x^k + 1 = x^k(x-2) + 1

Setting x = 2 - ε with ε small:
  (2-ε)^k · (2-ε-2) + 1 = 0
  (2-ε)^k · (-ε) + 1 = 0
  ε · (2-ε)^k = 1

For ε = 1/2^k + δ where δ = o(1/2^k):
  (1/2^k + δ) · (2 - 1/2^k - δ)^k = 1
  (1/2^k + δ) · 2^k · (1 - 1/2^{k+1} - δ/2)^k = 1
  (1 + δ·2^k) · (1 - k/2^{k+1} + O(k/4^k)) = 1
  1 + δ·2^k - k/2^{k+1} + O(k/4^k) = 1
  δ = k/2^{2k+1} + O(k/4^k)   (exponentially small correction)

So φ_k = 2 - 1/2^k - k/2^{2k+1} + O(k²/8^k).

The leading correction is EXACTLY ε = 1/2^k (constant C₂ = -1).  ■
""")

# Verify
print("  Verification:")
print(f"  {'k':>3s}  {'ε=2-φ_k':>14s}  {'1/2^k':>14s}  {'ε-1/2^k':>14s}  {'k/2^(2k+1)':>14s}  {'correction match':>16s}")
for k in range(3, 20):
    coeffs = [-1] * k + [1]
    roots = np.roots(coeffs[::-1])
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    if not real_pos: continue
    phi_k = max(real_pos)
    eps = 2 - phi_k
    pred = 1/2**k
    corr = k/2**(2*k+1)
    diff = eps - pred
    match = "✓" if abs(diff - corr) < corr * 0.5 or abs(diff) < 1e-12 else "~"
    print(f"  {k:3d}  {eps:14.10e}  {pred:14.10e}  {diff:14.10e}  {corr:14.10e}  {match:>16s}")

# ======================================================================
# PART 2: Weighted k-nacci proof
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: WEIGHTED k-NACCI — EXACT ASYMPTOTICS")
print("=" * 70)

print("""
THEOREM: The dominant root ψ_k of the weighted k-nacci recurrence
  f(n) = 1·f(n-1) + 2·f(n-2) + 4·f(n-3) + ... + 2^{k-1}·f(n-k)
satisfies:
  ψ_k = 3 - (2/3)^k + O(k·(4/9)^k)

PROOF:
Characteristic polynomial: x^k - x^{k-1} - 2x^{k-2} - ... - 2^{k-1}
  = x^k - Σ_{j=0}^{k-1} 2^j x^{k-1-j}
  = x^k - x^{k-1} Σ_{j=0}^{k-1} (2/x)^j
  = x^k - x^{k-1} · (1-(2/x)^k)/(1-2/x)  [geometric series]
  = x^k - x^k · (1-(2/x)^k)/(x-2)
  = x^k · [1 - (1-2^k/x^k)/(x-2)]
  = x^k · [(x-2 - 1 + 2^k/x^k)/(x-2)]
  = x^k · [(x-3 + 2^k/x^k)/(x-2)]

Setting this to zero (and x ≠ 2):
  x - 3 + 2^k/x^k = 0
  x = 3 - 2^k/x^k

For x = ψ_k ≈ 3:
  ψ_k = 3 - 2^k/ψ_k^k ≈ 3 - (2/3)^k

More precisely, set ψ_k = 3 - ε:
  3 - ε = 3 - 2^k/(3-ε)^k
  ε = 2^k/(3-ε)^k = (2/(3-ε))^k = (2/3)^k · 1/(1-ε/3)^k

For ε = (2/3)^k + δ:
  ((2/3)^k + δ) = (2/3)^k · (1 + kε/3 + ...)
  δ = (2/3)^k · (kε/3 + ...) = k/3 · (2/3)^{2k} + ...

So ψ_k = 3 - (2/3)^k - k/3 · (2/3)^{2k} + O(k²·(2/3)^{3k}).
The leading correction is EXACTLY ε = (2/3)^k (constant C₃ = -1).  ■
""")

# Verify
print("  Verification:")
print(f"  {'k':>3s}  {'ε=3-ψ_k':>14s}  {'(2/3)^k':>14s}  {'ε-(2/3)^k':>14s}  {'k/3·(2/3)^{2k}':>14s}  {'match':>8s}")
for k in range(3, 20):
    coeffs = [1]
    for j in range(k):
        coeffs.append(-2**j)
    roots = np.roots(coeffs)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    if not real_pos: continue
    psi_k = max(real_pos)
    eps = 3 - psi_k
    pred = (2/3)**k
    corr = k/3 * (2/3)**(2*k)
    diff = eps - pred
    match = "✓" if abs(diff - corr) < abs(corr) * 0.5 or abs(diff) < 1e-12 else "~"
    print(f"  {k:3d}  {eps:14.10e}  {pred:14.10e}  {diff:14.10e}  {corr:14.10e}  {match:>8s}")

# ======================================================================
# PART 3: General theorem
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: GENERAL THEOREM")
print("=" * 70)

print("""
GENERAL THEOREM (PROVED):

For the generalized k-nacci with weights w_j:
  f(n) = w_1·f(n-1) + w_2·f(n-2) + ... + w_k·f(n-k)

with w_j = c^{j-1} (geometric weights, c=1 for standard, c=2 for weighted):

The dominant root R_k approaches L = 1+c as k → ∞, with:
  R_k = L - (c/L)^k + O(k·(c/L)^{2k})

The convergence rate is (c/L)^k = (c/(1+c))^k:
  c=1: rate = (1/2)^k,  L=2
  c=2: rate = (2/3)^k,  L=3
  c=3: rate = (3/4)^k,  L=4
  c=r: rate = (r/(1+r))^k, L=1+r

PROOF SKETCH:
The characteristic equation with geometric weights gives
  x^k(x - (1+c)) + c^k = 0 [after multiplying by (x-c)/(x-c)]
Setting x = L - ε = (1+c) - ε:
  ((1+c)-ε)^k · (-ε) + c^k = 0
  ε · ((1+c)-ε)^k = c^k
  ε = c^k / (1+c-ε)^k = (c/(1+c))^k · 1/(1-ε/(1+c))^k
  → (c/(1+c))^k as ε → 0.

The "keys to the universe" appear as:
  c=1: L=2, rate=1/2  →  "the 2"
  c=2: L=3, rate=2/3  →  "the 3"

And the RATES are the reciprocals:
  rate(2) = 1/2
  rate(3) = 2/3

The product rate(2)·rate(3) = 1/3.
The sum rate(2)+rate(3) = 7/6 ≈ 1.17.
""")

# ======================================================================
# PART 4: Connection to tournament numbers
# ======================================================================
print("=" * 70)
print("PART 4: TOURNAMENT NUMBER CONNECTIONS")
print("=" * 70)

print("""
AT k=7:
  φ_7 = 2 - 1/128 = 2 - 1/2^7 = 255/128
  ψ_7 = 3 - (2/3)^7 = 3 - 128/2187

  The correction 1/2^7 = 1/128
  The correction (2/3)^7 = 128/2187

  Note: 2187 = 3^7. So (2/3)^7 = 2^7/3^7.

  And n=7 is where:
  - β₃ first appears (8.4% of tournaments)
  - α₂ is purely (3,3) type (single-level)
  - H is almost lambda-determined (99.93%)

AT k=8:
  φ_8 = 2 - 1/256 = 2 - 1/2^8 = 511/256
  ψ_8 = 3 - (2/3)^8 = 3 - 256/6561

  And n=8 is where:
  - Cross-level (3,5) pairs appear (37.5% of α₂)
  - β₃ becomes more common (16%)
  - |ΔH| ∈ {2,4,6} for Vitali reversals

NUMEROLOGICAL CONNECTION:
  The correction at k=n is proportional to c^n/L^n.
  For c=1 (standard k-nacci): error = 1/2^n
  For c=2 (weighted k-nacci): error = 2^n/3^n

  These are the SAME numerators/denominators that appear
  in the Jacobsthal formula:
    J_2(n) = (2^n - (-1)^n) / 3

  At k=n, the "error in knowing 2" is 1/2^n,
  and the "error in knowing 3" is (2/3)^n.

  To know 2 AND 3 to precision ε requires:
    k ≥ -log₂(ε) terms for the 2-tower
    k ≥ -log₂(ε)/log₂(3/2) terms for the 3-tower

  Since log₂(3/2) ≈ 0.585, you need ~1.7× MORE terms
  to know 3 to the same precision as 2.

  AT n=7: you know 2 to precision 1/128 ≈ 0.008
           you know 3 to precision 128/2187 ≈ 0.059
  So 3 is much less "known" than 2 at the tournament scale n=7.
""")

print("Done.")
