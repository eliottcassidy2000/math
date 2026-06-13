#!/usr/bin/env python3
"""
gk_naive_formula_89c.py — Test naive formula g_k(m) = Σ C(k-1,r-1)·C(m,r)·2^{r-1}
opus-2026-03-15-S89c

Motivated by:
  h_s(n) = 2/(n)_{2s} for all cluster sizes s
  Tilings with r clusters: count = C(k-1,r-1)·C(m,r)
  If weight per tiling is 2^r/(n)_{2k}, then g_k = Σ C(k-1,r-1)·C(m,r)·2^{r-1}
"""

from math import comb
from fractions import Fraction

def g_naive(k, m):
    """Naive formula: g_k(m) = Σ_{r=1}^{k} C(k-1,r-1)·C(m,r)·2^{r-1}"""
    return sum(comb(k-1, r-1) * comb(m, r) * 2**(r-1) for r in range(1, k+1))

def g_true(k, m):
    """True g_k from verified polynomials."""
    if k == 1: return m
    if k == 2: return m * m
    if k == 3: return (2*m**3 + m) // 3
    if k == 4: return (10*m**3 - 33*m**2 + 50*m - 24) // 3
    if k == 5: return (388*m**3 - 2040*m**2 + 3431*m - 1776) // 3
    if k == 6: return (69660*m**3 - 380445*m**2 + 653748*m - 342960) // 3
    if k == 7: return (19826270*m**3 - 109486152*m**2 + 189674605*m - 100014720) // 3
    return None

print("="*70)
print("NAIVE vs TRUE g_k VALUES")
print("="*70)

for k in range(1, 8):
    print(f"\nk={k}:")
    for m in range(1, 8):
        naive = g_naive(k, m)
        true = g_true(k, m)
        diff = naive - true if true is not None else "?"
        ok = "✓" if diff == 0 else f"diff={diff}"
        print(f"  m={m}: naive={naive:>10}, true={true:>10}  {ok}")

# The corrections: true - naive
print("\n" + "="*70)
print("CORRECTION: E_k(m) = g_k(m) - naive_k(m)")
print("="*70)

for k in range(4, 8):
    print(f"\nk={k}: E_{k}(m) =")
    corrections = []
    for m in range(1, 10):
        true = g_true(k, m)
        naive = g_naive(k, m)
        if true is not None:
            corr = true - naive
            corrections.append((m, corr))
            print(f"  E_{k}({m}) = {corr}")

    # Check if correction is a polynomial in m
    # For k=4: correction should be -(m-1)(m-2)(m-3)(m-4)/3
    if k == 4:
        print(f"\n  Check: -(m-1)(m-2)(m-3)(m-4)/3 =")
        for m, corr in corrections:
            pred = -((m-1)*(m-2)*(m-3)*(m-4)) // 3
            match = "✓" if pred == corr else f"✗ (pred={pred})"
            print(f"    m={m}: {pred}  {match}")

# Factor out the degree-4 corrections
print("\n" + "="*70)
print("DEGREE-4 CORRECTION ANALYSIS")
print("="*70)

# For k=4: correction = -(m-1)(m-2)(m-3)(m-4)/3 = -4!C(m-1,4)/3 = -8C(m-1,4)
# Check: C(m-1,4) for m=1..7: 0,0,0,0,1,5,15
# -8C(m-1,4): 0,0,0,0,-8,-40,-120
# Difference (naive-true) for k=4: 0,0,0,0,8,40,120 → correction is -8C(m-1,4) ✓

for k in range(4, 8):
    print(f"\nk={k}: correction pattern analysis")
    corrections = []
    for m in range(1, 10):
        true = g_true(k, m)
        naive = g_naive(k, m)
        if true is not None:
            corrections.append(true - naive)

    # Try to identify the correction polynomial
    # Correction starts being nonzero at m = ?
    first_nonzero = next((i+1 for i, c in enumerate(corrections) if c != 0), None)
    if first_nonzero:
        print(f"  First nonzero at m={first_nonzero}")

    print(f"  corrections: {corrections}")

    # For k=4: corrections involve C(m-1,4) (nonzero for m>=5)
    # For k=5: likely involves C(m-1,4) and C(m-1,5) terms
    # Check if corrections[:4] are all 0 (they should be for r≤3 formula working)
    zeros_at_start = sum(1 for c in corrections if c == 0)
    print(f"  Leading zeros: {zeros_at_start}")

# Now: the "binomial" basis expansion
# g_k(m) = Σ_{r≥1} a_{k,r} · C(m,r)
# where a_{k,r} = C(k-1,r-1)·2^{r-1} for r ≤ 3 (from the naive formula)
# but a_{k,r} differs for r ≥ 4

print("\n" + "="*70)
print("BINOMIAL BASIS EXPANSION: g_k(m) = Σ a_{k,r}·C(m,r)")
print("="*70)

for k in range(1, 8):
    true_vals = [g_true(k, m) for m in range(1, max(k+3, 8))]
    if any(v is None for v in true_vals):
        continue

    # Convert to binomial basis using forward differences
    # g_k(m) = Σ_{r=0}^{deg} Δ^r g_k(0) · C(m,r)
    # where Δ^r f(0) = Σ_{j=0}^{r} (-1)^{r-j} C(r,j) f(j)

    # Compute g_k(0), g_k(1), ..., g_k(k+2) for the forward differences
    vals = []
    for m in range(0, k + 4):
        v = g_true(k, m)
        if v is None:
            # Use the polynomial to compute
            break
        vals.append(v)

    # Forward difference table
    deltas = [vals[0]]
    curr = list(vals)
    for r in range(1, len(curr)):
        next_diff = [curr[i+1] - curr[i] for i in range(len(curr)-1)]
        if next_diff:
            deltas.append(next_diff[0])
        curr = next_diff

    print(f"\nk={k}: g_{k}(m) = ", end="")
    terms = []
    for r, d in enumerate(deltas):
        if d != 0:
            naive_coeff = comb(k-1, r-1) * 2**(r-1) if r >= 1 else 0
            if r == 0:
                terms.append(f"{d}·C(m,0)")
            else:
                match = "✓" if d == naive_coeff else f"(naive={naive_coeff})"
                terms.append(f"{d}·C(m,{r}) {match}")
    print(" + ".join(terms))

print("\nDone!")
