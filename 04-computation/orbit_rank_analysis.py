#!/usr/bin/env python3
"""
Analysis of orbit complex rank patterns for Paley tournaments.

Known orbit data:
  P_7 (m=3): Omega_orb = [1,1,2,3,3,2,1], R_orb = [0,0,1,1,2,1,1,0], β_orb = [1,0,0,0,0,0,0]
  P_11 (m=5): Omega_orb = [1,1,4,14,41,92,140,138,90,36,6], R_orb = [0,0,1,3,11,30,61,78,60,30,6,0]
               β_orb = [1,0,0,0,0,1,1,0,0,0,0]

Questions:
1. Is there a pattern in R_orb values?
2. Can we predict R_orb for P_19 (m=9)?
3. Why does β_m^orb = (m-3)/2?

Key observation: acyclicity holds for d < m and d > m+1.
So R_d^orb is determined by the recursion R_d = Omega_{d-1}^orb - R_{d-1} for d=1..m
and R_d = Omega_{d-1}^orb - R_{d-1} for d=m+2..2m (from the top).

This means R_orb is determined by Omega_orb, except for the "gap" at d=m,m+1.

opus-2026-03-13-S71b
"""

# Known data
data = {
    7: {
        'm': 3,
        'Omega_orb': [1, 1, 2, 3, 3, 2, 1],
        'R_orb': [0, 0, 1, 1, 2, 1, 1, 0],
        'beta_orb': [1, 0, 0, 0, 0, 0, 0],
    },
    11: {
        'm': 5,
        'Omega_orb': [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6],
        'R_orb': [0, 0, 1, 3, 11, 30, 61, 78, 60, 30, 6, 0],
        'beta_orb': [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    }
}

for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    R = d['R_orb']
    beta = d['beta_orb']

    print(f"=== P_{p} (m={m}) ===")
    print(f"Omega_orb: {O}")
    print(f"R_orb:     {R}")
    print(f"β_orb:     {beta}")

    # Bottom recursion: R_1 = 0, R_d = Omega_{d-1} - R_{d-1} for d ≤ m
    print(f"\nBottom recursion (d=0..m):")
    R_bottom = [0]  # R_0
    for dd in range(1, m + 1):
        r_next = O[dd - 1] - R_bottom[-1]
        R_bottom.append(r_next)
        match = "✓" if r_next == R[dd] else "✗"
        print(f"  R_{dd} = Omega_{dd-1} - R_{dd-1} = {O[dd-1]} - {R_bottom[-2]} = {r_next} {match}")

    # Top recursion: R_{2m+1} = 0, R_d = Omega_d - R_{d+1} for d ≥ m+2
    print(f"\nTop recursion (d=2m..m+2):")
    R_top = {2*m + 1: 0}
    for dd in range(2*m, m + 1, -1):
        r_prev = O[dd] - R_top[dd + 1]
        R_top[dd] = r_prev
        match = "✓" if r_prev == R[dd] else "✗"
        print(f"  R_{dd} = Omega_{dd} - R_{dd+1} = {O[dd]} - {R_top[dd+1]} = {r_prev} {match}")

    # Gap at d=m and d=m+1
    print(f"\nGap at d={m},{m+1}:")
    R_m = R_bottom[m]
    R_m2 = R_top.get(m + 2, R[m+2])
    budget = O[m] - R_m
    budget2 = O[m + 1] - R_m2
    print(f"  R_m = {R_m}")
    print(f"  R_{{m+2}} = {R_m2}")
    print(f"  Budget at d=m: Omega_m - R_m = {budget}")
    print(f"  Budget at d=m+1: Omega_{{m+1}} - R_{{m+2}} = {budget2}")
    print(f"  R_{{m+1}} = Budget - β_m = {budget} - {beta[m]} = {budget - beta[m]}")
    print(f"  R_{{m+1}} check = Budget2 - β_{{m+1}} = {budget2} - {beta[m+1]} = {budget2 - beta[m+1]}")
    print(f"  Match: R_{{m+1}} = {R[m+1]}")
    print()

    # Alternating sum formula for R_m
    print(f"R_m as alternating sum of Omega_orb:")
    alt_sum = 0
    terms = []
    for dd in range(m):
        sign = (-1) ** (m - 1 - dd)
        alt_sum += sign * O[dd]
        terms.append(f"({sign:+d})*{O[dd]}")
    print(f"  R_m = {' + '.join(terms)} = {alt_sum}")
    print(f"  Actual R_m = {R_m}")

    # R_{m+2} as alternating sum from top
    print(f"\nR_{{m+2}} as alternating sum:")
    alt_sum_top = 0
    terms_top = []
    for dd in range(2*m, m+1, -1):
        sign = (-1) ** (dd - m - 2)
        alt_sum_top += sign * O[dd]
        terms_top.append(f"({sign:+d})*{O[dd]}")
    print(f"  R_{{m+2}} = {' + '.join(terms_top)} = {alt_sum_top}")
    print(f"  Actual R_{{m+2}} = {R_m2}")
    print()

    # Symmetry analysis: R_d^orb vs R_{2m+1-d}^orb
    print(f"R_d symmetry (R_d vs R_{{2m+1-d}}):")
    for dd in range(2*m + 2):
        r1 = R[dd]
        r2 = R[2*m + 1 - dd] if 0 <= 2*m + 1 - dd < len(R) else "N/A"
        print(f"  R_{dd}={r1}, R_{2*m+1-dd}={r2}", end="")
        if isinstance(r2, int):
            print(f"  ratio={r1/r2:.4f}" if r2 != 0 else "  (r2=0)", end="")
        print()

    print()

# Now: relationship between Omega_orb and Omega (full)
print("=== OMEGA ORBIT vs FULL ===\n")
for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    print(f"P_{p} (m={m}): Omega_orb = {O}")

    # Check if Omega_orb[d] has a pattern
    # Omega_d / m for d >= 1
    # For P_7: Omega = [1,3,6,9,9,6,3]
    # For P_11: Omega = [1,5,20,70,205,460,700,690,450,180,30]
    if p == 7:
        Omega_full = [1,3,6,9,9,6,3]
    else:
        Omega_full = [1,5,20,70,205,460,700,690,450,180,30]

    print(f"  Omega_full = {Omega_full}")
    print(f"  Omega_full/m: {[Omega_full[dd]//m for dd in range(len(Omega_full))]}")
    print(f"  Check: {[Omega_full[dd] % m for dd in range(1, len(Omega_full))]}")
    print()

# COMBINATORIAL ANALYSIS: what IS β_m^{orb}?
print("=== β_m^{orb} = (m-3)/2 ANALYSIS ===\n")
print("Values:")
for m in range(3, 20, 2):
    p = 2*m + 1
    beta = (m - 3) // 2
    print(f"  m={m}, p={p}: β_m^orb = {beta}")

print(f"\nThis is:")
print(f"  β_m^orb = (m-3)/2 = (p-7)/4")
print(f"  β_m = m*(m-3)/2")
print(f"  β_{{m+1}} = m*(m-3)/2 + (p-1) = m*(m-3)/2 + 2m = m*(m+1)/2 = C(m+1,2)")
print()

# KEY INSIGHT: The orbit complex has 2m+1 = p dimensions.
# It's acyclic except at d=m,m+1.
# The Euler characteristic is 1.
# So β_m = β_{m+1} = β_m^{orb}.
# And β_m^{orb} is the "topological complexity" of the orbit complex.

# The orbit complex is a chain complex with dimensions:
# 1, 1, O_2/m, O_3/m, ..., O_{2m}/m
# For P_11: 1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6
# Sum = 563

# It's acyclic up to d=m-1 (bottom), and acyclic from d=m+2 (top).
# The gap at d=m gives β_m = β_{m+1} = 1.

# Why 1? The bottom alternating sum R_m = Σ (-1)^{m-1-d} Omega_d^{orb}
# And the top alternating sum R_{m+2} = Σ (-1)^{d-m-2} Omega_d^{orb} for d=m+2..2m.

# β_m = Omega_m^{orb} - R_m - R_{m+1}
# For P_11: 92 - 30 - 61 = 1. The CLOSENESS of R_m + R_{m+1} to Omega_m gives β=1.

# R_m (bottom) = 1 - 1 + 4 - 14 + 41 = 31. Wait, that should be 30...
# For P_11: alternating from bottom: R_1=0, R_2=Ω_1-R_1=1, R_3=Ω_2-R_2=4-1=3,
# R_4=Ω_3-R_3=14-3=11, R_5=Ω_4-R_4=41-11=30. ✓

# R_m = Ω_{m-1} - Ω_{m-2} + Ω_{m-3} - ... + (-1)^{m-1} Ω_0
#      = Σ_{d=0}^{m-1} (-1)^{m-1-d} Ω_d^{orb}

# For m=5: R_5 = Ω_4 - Ω_3 + Ω_2 - Ω_1 + Ω_0 = 41 - 14 + 4 - 1 + 1 = 31.
# But actual R_5 = 30! Discrepancy = 1. Where does the -1 come from?

print("=== ALTERNATING SUM DISCREPANCY ===\n")
for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    R = d['R_orb']

    # R_m via recursion
    R_m_recursion = R[m]

    # R_m via alternating sum formula
    alt_sum = sum((-1)**(m-1-dd) * O[dd] for dd in range(m))
    print(f"P_{p}: R_m={R_m_recursion}, alt_sum={alt_sum}, diff={alt_sum - R_m_recursion}")

    # The discrepancy is due to Ω_0^{orb} = 1 but the recursion starts with R_0=0.
    # R_1 = Ω_0 - R_0 = 1 - 0 = 1? But R_1 = 0!
    # Wait: R_1 = rank(∂_1^{orb}). The boundary ∂_1 maps the single edge orbit to
    # the single vertex: ∂(s_1) = () - () = 0 (face_0 and face_1 are both ()).
    # So R_1 = 0, meaning the edge is a cycle, H_1^{orb} = 0 only if ∂_2 is surjective.
    #
    # So the recursion is NOT R_1 = Ω_0 - R_0 (which would give 1).
    # R_1 = 0 because ∂_1 is zero. Then β_0 = Ω_0 - R_0 - R_1 = 1 - 0 - 0 = 1.
    # But β_1 = Ω_1 - R_1 - R_2 = 1 - 0 - R_2. Need R_2 = 1 for β_1 = 0.
    #
    # The acyclicity at d=1 says: β_1 = 0, so R_2 = Ω_1 - R_1 = 1 - 0 = 1. ✓
    # At d=2: β_2 = 0, so R_3 = Ω_2 - R_2 = Ω_2^{orb} - 1.
    # ...
    # At d=m-1: β_{m-1} = 0, so R_m = Ω_{m-1} - R_{m-1}.

    # Now: the alternating sum R_m = Σ (-1)^{m-1-d} Ω_d^{orb} for d=1..m-1
    # (starting from R_1 = 0, not from R_0 = 0)
    # R_2 = Ω_1 - R_1 = Ω_1 = 1
    # R_3 = Ω_2 - R_2 = Ω_2 - Ω_1 + 0 = Ω_2 - 1
    # ...
    # R_m = Σ_{d=1}^{m-1} (-1)^{m-1-d} Ω_d^{orb}

    alt_sum_correct = sum((-1)**(m-1-dd) * O[dd] for dd in range(1, m))
    print(f"  Correct alt_sum (d=1..m-1): {alt_sum_correct}")
    print(f"  Match: {'✓' if alt_sum_correct == R_m_recursion else '✗'}")

    # R_{m+2} from top
    R_m2 = R[m+2] if m+2 < len(R) else 0
    alt_top = sum((-1)**(dd-m-2) * O[dd] for dd in range(m+2, 2*m+1))
    print(f"  R_{{m+2}}: recursion={R_m2}, alt_sum(d=m+2..2m)={alt_top}")

    # Budget
    budget_m = O[m] - R_m_recursion
    budget_m1 = O[m+1] - R_m2
    R_m1 = R[m+1]
    print(f"  Budget_m = Ω_m - R_m = {budget_m}")
    print(f"  Budget_{{m+1}} = Ω_{{m+1}} - R_{{m+2}} = {budget_m1}")
    print(f"  R_{{m+1}} = {R_m1}")
    print(f"  β_m = budget_m - R_{{m+1}} = {budget_m - R_m1}")
    print(f"  β_{{m+1}} = budget_{{m+1}} - R_{{m+1}} = {budget_m1 - R_m1}")
    print()

print("\n=== BUDGET FORMULA ===\n")
print("For the orbit complex:")
print("  R_m = Σ_{d=1}^{m-1} (-1)^{m-1-d} Ω_d^{orb}")
print("  R_{m+2} = Σ_{d=m+2}^{2m} (-1)^{d-m-2} Ω_d^{orb}")
print("  Budget_m = Ω_m^{orb} - R_m")
print("  Budget_{m+1} = Ω_{m+1}^{orb} - R_{m+2}")
print()
print("  β_m = β_{m+1} since Euler char = 1 and all other β = 0.")
print("  So: Budget_m + Budget_{m+1} = R_{m+1} + R_{m+1} + β_m + β_{m+1}")
print("      = 2R_{m+1} + 2β_m")
print("  And: R_{m+1} = Budget_m - β_m = Budget_{m+1} - β_{m+1}")
print("  So: β_m = (Budget_m - Budget_{m+1}) / 2 + ... hmm no.")
print()

# Actually: chi_orb = 1 means β_0 + Σ_d (-1)^d β_d = 1
# With only β_0, β_m, β_{m+1} nonzero:
# 1 + (-1)^m β_m + (-1)^{m+1} β_{m+1} = 1
# So (-1)^m (β_m - β_{m+1}) = 0, hence β_m = β_{m+1}.
# And Budget_m = R_{m+1} + β_m, Budget_{m+1} = R_{m+1} + β_{m+1} = R_{m+1} + β_m.
# So Budget_m = Budget_{m+1}!

print("CHECK: Budget_m = Budget_{m+1}?")
for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    R = d['R_orb']
    R_m = R[m]
    R_m2 = R[m+2]
    budget_m = O[m] - R_m
    budget_m1 = O[m+1] - R_m2
    print(f"  P_{p}: Budget_m={budget_m}, Budget_{{m+1}}={budget_m1} "
          f"{'✓ EQUAL' if budget_m == budget_m1 else '✗ NOT EQUAL'}")

print()
print("If Budget_m = Budget_{m+1} = B, then:")
print("  β_m = B - R_{m+1}")
print("  β_{m+1} = B - R_{m+1}")
print("  So we need to determine R_{m+1} (or equivalently β_m).")
print()
print("R_{m+1} is NOT determined by acyclicity alone.")
print("It depends on the actual rank of ∂_{m+1}^{orb} restricted to Omega_{m+1}^{orb}.")
print()
print("The question becomes: what is the CORANK of ∂_{m+1}^{orb} at degree m+1?")
print("corank = dim(ker) = Omega_{m+1}^{orb} - R_{m+1}")
print()
for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    R = d['R_orb']
    ker_m1 = O[m+1] - R[m+1]
    print(f"  P_{p}: ker(∂_{{m+1}}^orb) = {O[m+1]} - {R[m+1]} = {ker_m1}")
    print(f"    = R_{{m+2}} + β_{{m+1}} = {R[m+2]} + {d['beta_orb'][m+1]} = {R[m+2] + d['beta_orb'][m+1]}")
