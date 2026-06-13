#!/usr/bin/env python3
"""
P_11 Betti number prediction from Ω dims.

Known Ω dims (k=1 eigenspace): [1, 5, 20, 70, 205, 460, 700, 690, ?, ?, ?]
Constraint: alt sum = 1

Per-eigenspace chain complex:
  ∂_d: Ω_d → Ω_{d-1}

The boundary maps have ranks r_d = rank(∂_d: Ω_d → Ω_{d-1}).
Then: β_d = (Ω_d - r_d) - r_{d+1} = ker_d - im_{d+1}

Constraints:
  r_d ≤ min(Ω_d, Ω_{d-1})
  β_d ≥ 0
  Σ (-1)^d β_d = Σ (-1)^d Ω_d = 1 (alternating sum)
  Σ (-1)^d Ω_d = Σ (-1)^d (r_d + β_d + r_{d+1})

For P_7: Ω dims = [1, 3, 6, 9, 9, 6, 3], boundary ranks r = [0, 1, 3, 6, 5, 3, 3]
Then: ker_d = Ω_d - r_d, β_d = ker_d - r_{d+1}
  β_0 = 1-1 = 0... wait, P_7 per-eigenspace (k≠0) has β = (0,0,0,0,1,0,0)

Actually for the per-eigenspace:
  Ω_0 = 1, r_0 = 0 (∂_0 = 0), r_1 = rank(∂_1: Ω_1 → Ω_0)
  ker_0 = Ω_0 - 0 = 1
  β_0 = ker_0 - r_1

For k≠0, the eigenspace k=1 at P_7:
  The ∂ map preserves eigenspace, and β_0^(k) = 0 for k≠0.

Let me just compute what's constrainted by the Ω dims alone.

If β is concentrated at dimension d (all other β = 0), then:
  Ω_d = r_d + β_d + r_{d+1}
  For all i ≠ d: Ω_i = r_i + r_{i+1} (since β_i = 0)

This gives a system of equations. Let me solve for each possible d.
"""
import numpy as np

# Known Ω dims
# omega = [1, 5, 20, 70, 205, 460, 700, 690, Ω8, Ω9, Ω10]
# Alt sum = 1, so Ω8 - Ω9 + Ω10 = 300

# For P_7 reference:
omega_7 = [1, 3, 6, 9, 9, 6, 3]
print("P_7 per-eigenspace (k=1):")
print(f"  Ω = {omega_7}")

# If β concentrated at d=4: β_4=1, all others 0
# Then: r_0=0, r_1=1, r_2=3-1=2... let me work through
# Actually: r_0=0, and Ω_0=r_0+r_1 (since β_0=0) → r_1 = Ω_0 = 1
# Ω_1 = r_1 + r_2 → r_2 = Ω_1 - r_1 = 3-1 = 2
# Ω_2 = r_2 + r_3 → r_3 = Ω_2 - r_2 = 6-2 = 4
# Ω_3 = r_3 + r_4 → r_4 = Ω_3 - r_3 = 9-4 = 5
# Ω_4 = r_4 + β_4 + r_5 → r_5 = Ω_4 - r_4 - β_4 = 9-5-1 = 3
# Ω_5 = r_5 + r_6 → r_6 = Ω_5 - r_5 = 6-3 = 3
# Ω_6 = r_6 + 0 (no r_7, since dim 7 has Ω_7=0 for p=7)
# r_6 = 3 ≤ Ω_6 = 3 ✓

print(f"\n  β concentrated at d=4:")
r = [0]
for d in range(1, 7):
    if d <= 4:
        r_next = omega_7[d-1] - r[d-1]
    elif d == 5:
        r_next = omega_7[4] - r[4] - 1  # β_4 = 1
    else:
        r_next = omega_7[d-1] - r[d-1]
    r.append(r_next)

# Cleaner: for β at position d_beta with value 1
def solve_ranks(omega, d_beta, beta_val=1):
    """Given Ω dims and β_d = beta_val (all others 0),
    compute boundary ranks r_1, ..., r_n."""
    n = len(omega)
    # r_0 = 0 by convention
    r = [0] * (n + 1)
    for d in range(1, n + 1):
        if d - 1 < n:
            if d - 1 == d_beta:
                r[d] = omega[d-1] - r[d-1] - beta_val
            else:
                r[d] = omega[d-1] - r[d-1]

    # Verify all ranks non-negative
    valid = all(r[d] >= 0 for d in range(n + 1))
    # Verify last rank is consistent (r_n should match Ω_{n-1})
    return r, valid

print("\nP_7: Testing all possible concentration dimensions:")
for d_beta in range(7):
    r, valid = solve_ranks(omega_7, d_beta, 1)
    if valid:
        print(f"  d={d_beta}: r = {r[1:]}, valid={valid}")

# ===== P_11 =====
print(f"\n\n{'='*70}")
print("P_11: PREDICTING β CONCENTRATION")
print("="*70)

omega_11_partial = [1, 5, 20, 70, 205, 460, 700, 690]
# Need Ω_8, Ω_9, Ω_10 with Ω_8 - Ω_9 + Ω_10 = 300

# For each possible d_beta, compute what Ω_8, Ω_9, Ω_10 would need to be
print("\nFor each concentration dimension d, assuming β_d = 1:")

for d_beta in range(11):
    # Try to solve: starting from r_0 = 0, compute ranks
    # For dims 0-7 we know Ω, for 8-10 we need to determine

    # First compute r_1 through r_8 from known Ω dims
    r = [0]
    valid = True
    for d in range(1, 9):
        if d - 1 < len(omega_11_partial):
            if d - 1 == d_beta:
                r_next = omega_11_partial[d-1] - r[d-1] - 1
            else:
                r_next = omega_11_partial[d-1] - r[d-1]
            if r_next < 0:
                valid = False
                break
            r.append(r_next)

    if not valid:
        print(f"  d={d_beta}: INVALID (negative rank before dim 8)")
        continue

    # Now for dims 8-10, we need to determine Ω_8, Ω_9, Ω_10
    # and also r_9, r_10, r_11
    # Constraints: Ω_8 - Ω_9 + Ω_10 = 300, all non-negative
    # Ω_d = r_d + (1 if d==d_beta else 0) + r_{d+1}

    # For d_beta ≤ 7: all dims 8-10 have β=0
    if d_beta <= 7:
        # r_8 is computed above, Ω_8 = r_8 + r_9
        # r_9 = Ω_8 - r_8, r_10 = Ω_9 - r_9
        # We need Ω_8, Ω_9, Ω_10 satisfying the alt sum
        # AND r_9 ≥ 0, r_10 ≥ 0, r_11 ≥ 0
        # r_11 = Ω_10 - r_10 ≥ 0
        # Also: dim 11 has Ω_11 = 0, so r_11 = Ω_10 - r_10 and Ω_10 = r_10 + r_11
        # So Ω_10 = r_10 + r_11 ≥ 0 (trivially)
        # And we need Ω_10 ≥ 0, Ω_9 ≥ 0, Ω_8 ≥ 0

        # Express: Ω_8 = r_8 + r_9, Ω_9 = r_9 + r_10, Ω_10 = r_10 + r_11
        # So Ω_8 - Ω_9 + Ω_10 = r_8 + r_9 - r_9 - r_10 + r_10 + r_11 = r_8 + r_11 = 300
        # Therefore r_11 = 300 - r_8

        r_8_val = r[8]
        r_11_val = 300 - r_8_val

        if r_11_val < 0:
            print(f"  d={d_beta}: INVALID (r_11 = {r_11_val} < 0)")
            continue

        print(f"  d={d_beta}: r[1:9]={r[1:]}, r_8={r_8_val}, r_11={r_11_val}")
        print(f"    Ω_8 = r_8 + r_9 (free), Ω_9 = r_9 + r_10 (free), Ω_10 = r_10 + r_11")
        print(f"    Constraint: r_11 = {r_11_val}, all r ≥ 0")

    elif d_beta == 8:
        r_8_val = r[8]
        # Ω_8 = r_8 + 1 + r_9 → r_9 = Ω_8 - r_8 - 1
        # Ω_9 = r_9 + r_10, Ω_10 = r_10 + r_11
        # Ω_8 - Ω_9 + Ω_10 = (r_8+1+r_9) - (r_9+r_10) + (r_10+r_11)
        # = r_8 + 1 + r_11 = 300 → r_11 = 299 - r_8
        r_11_val = 299 - r_8_val
        print(f"  d={d_beta}: r[1:9]={r[1:]}, r_8={r_8_val}, r_11={r_11_val}")

    elif d_beta == 9:
        r_8_val = r[8]
        # Ω_8 = r_8 + r_9, Ω_9 = r_9 + 1 + r_10, Ω_10 = r_10 + r_11
        # Ω_8 - Ω_9 + Ω_10 = r_8 - 1 + r_11 = 300 → r_11 = 301 - r_8
        r_11_val = 301 - r_8_val
        print(f"  d={d_beta}: r[1:9]={r[1:]}, r_8={r_8_val}, r_11={r_11_val}")

    elif d_beta == 10:
        r_8_val = r[8]
        # Ω_10 = r_10 + 1 + r_11
        # Ω_8 - Ω_9 + Ω_10 = r_8 + 1 + r_11 = 300 → r_11 = 299 - r_8
        r_11_val = 299 - r_8_val
        print(f"  d={d_beta}: r[1:9]={r[1:]}, r_8={r_8_val}, r_11={r_11_val}")

# ===== Summary =====
print(f"\n\nSUMMARY:")
print(f"The boundary ranks r_1,...,r_8 are determined by Ω_0,...,Ω_7 and d_beta.")
print(f"For d_beta = p-3 = 8:")
r = [0]
for d in range(1, 9):
    if d - 1 == 8:
        r_next = omega_11_partial[d-1] - r[d-1] - 1
    else:
        r_next = omega_11_partial[d-1] - r[d-1]
    r.append(r_next)
print(f"  r = {r}")
print(f"  Need r_11 = 299 - r_8 = 299 - {r[8]} = {299 - r[8]}")

print("\nDone.")
