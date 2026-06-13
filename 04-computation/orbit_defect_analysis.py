#!/usr/bin/env python3
"""
Deep analysis of orbit complex defect.

The orbit complex has ∂_1 = 0 (both faces are the same vertex).
This means R_1^{orb} = 0 instead of the "generic" value 1.

In the FULL k=0 complex:
  ∂_1(s) = () - () = 0 for ALL s ∈ QR.
  So R_1^{(0)} = 0 and β_0^{(0)} = 1 (connected), β_1^{(0)} = m (all edges are cycles).
  Wait, β_1^{(0)} should be 0 for Paley tournaments!

Actually: R_1^{(0)} = 0 (since ∂_1 is identically zero for k=0).
  Then β_1 = Ω_1 - R_1 - R_2 = m - 0 - R_2.
  For β_1 = 0: R_2 = m.
  This means ∂_2 maps ONTO the entire Ω_1 space (rank m).

In the orbit complex:
  R_1^{orb} = 0, β_1^{orb} = Ω_1^{orb} - 0 - R_2^{orb} = 1 - R_2^{orb}.
  For β_1^{orb} = 0: R_2^{orb} = 1. ✓ (Confirmed for both P_7 and P_11)

  But now R_2^{orb} = 1 means the acyclicity recursion gives:
  R_3^{orb} = Ω_2^{orb} - R_2^{orb} = Ω_2^{orb} - 1

  Compare with the "generic" case where R_1 = 1:
  R_2_generic = Ω_1 - R_1 = 1 - 1 = 0
  R_3_generic = Ω_2 - R_2 = Ω_2
  R_4_generic = Ω_3 - R_3 = Ω_3 - Ω_2
  ...

  But the actual orbit complex has R_1 = 0:
  R_2 = Ω_1 - R_1 = 1 - 0 = 1
  R_3 = Ω_2 - R_2 = Ω_2 - 1
  R_4 = Ω_3 - R_3 = Ω_3 - Ω_2 + 1
  ...

  The difference: R_d^{actual} - R_d^{generic} = (-1)^{d+1} for d=1..m.
  This is EXACTLY the rank shift from the full complex!

  At d=1: R_1^{actual} - R_1^{generic} = 0 - 1 = -1 = (-1)^2
  At d=2: R_2^{actual} - R_2^{generic} = 1 - 0 = +1 = (-1)^3
  At d=3: R_3^{actual} - R_3^{generic} = (Ω_2-1) - Ω_2 = -1 = (-1)^4
  ...
  At d=m: R_m^{actual} - R_m^{generic} = (-1)^{m+1}

  So the orbit rank shift is IDENTICAL to the eigenspace rank shift!

  And at d=m+1:
  R_{m+1}^{actual} = Budget - β_m^{orb}
  R_{m+1}^{generic} = Budget_generic (from the other direction)
  The deficit at d=m+1 determines β_m^{orb}.

  KEY: What is R_m^{generic}?
  R_m^{generic} = Σ_{d=0}^{m-1} (-1)^{m-1-d} Ω_d^{orb}
  R_m^{actual} = Σ_{d=1}^{m-1} (-1)^{m-1-d} Ω_d^{orb}
  Difference = (-1)^{m-1} · Ω_0 = (-1)^{m-1} · 1

  Since m is odd: (-1)^{m-1} = 1, so R_m^{generic} - R_m^{actual} = 1.
  (And R_m^{actual} = R_m^{generic} - 1.)

  So Budget_m = Ω_m - R_m = Ω_m - (R_m^{generic} - 1) = (Ω_m - R_m^{generic}) + 1

  Similarly for the top:
  R_{m+2}^{generic} = Σ_{d=m+2}^{2m} (-1)^{d-m-2} Ω_d = R_{m+2}^{actual}
  (no defect from the top since the top recursion starts from R_{2m+1} = 0 normally)

  Budget_{m+1} = Ω_{m+1} - R_{m+2}

  And Budget_m = Budget_{m+1} (proved above).

  So: Ω_m - R_m^{generic} + 1 = Ω_{m+1} - R_{m+2}
  i.e., Ω_m - R_m^{generic} = Ω_{m+1} - R_{m+2} - 1

  This means the "generic" complex (with ∂_1 ≠ 0) would have:
  Budget_m^{generic} = Ω_m - R_m^{generic} = Budget_m - 1

  And β_m^{generic} = Budget_m^{generic} - R_{m+1}^{generic}
  But R_{m+1}^{generic} might also differ.

opus-2026-03-13-S71b
"""

# Let me compute all these quantities

data = {
    7: {
        'm': 3,
        'Omega_orb': [1, 1, 2, 3, 3, 2, 1],
        'R_orb': [0, 0, 1, 1, 2, 1, 1, 0],
    },
    11: {
        'm': 5,
        'Omega_orb': [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6],
        'R_orb': [0, 0, 1, 3, 11, 30, 61, 78, 60, 30, 6, 0],
    }
}

for p, d in data.items():
    m = d['m']
    O = d['Omega_orb']
    R = d['R_orb']

    print(f"=== P_{p} (m={m}) ===\n")

    # "Generic" recursion (as if R_1 = Ω_0 = 1)
    R_gen = [0, 1]  # R_0, R_1 generic
    for dd in range(2, m + 1):
        R_gen.append(O[dd-1] - R_gen[-1])

    print(f"R_actual: {R[:m+1]}")
    print(f"R_generic: {R_gen}")
    print(f"Difference: {[R[dd] - R_gen[dd] for dd in range(m+1)]}")
    print(f"Expected: {[(-1)**(dd+1) if dd > 0 else 0 for dd in range(m+1)]}")
    print()

    # Budget analysis
    R_m_actual = R[m]
    R_m_generic = R_gen[m]
    R_m2 = R[m+2]

    print(f"R_m: actual={R_m_actual}, generic={R_m_generic}, diff={R_m_actual - R_m_generic}")
    print(f"Budget: actual={O[m] - R_m_actual}, generic={O[m] - R_m_generic}")
    print(f"R_{{m+2}} = {R_m2}")

    # Euler characteristic of the "generic" orbit complex
    # If R_d^{gen} followed normal acyclicity (R_1=1), then:
    # β_d^{gen} = 0 for all d: completely acyclic.
    # Because the chain complex with ∂_1 rank 1 is acyclic.
    # R_{m+1}^{gen} = Budget^{gen} = O[m] - R_m^{gen}
    R_m1_gen = O[m] - R_m_generic
    print(f"R_{{m+1}}^gen = Budget^gen = {R_m1_gen}")

    # Check from top
    R_m1_from_top = O[m+1] - R_m2
    print(f"R_{{m+1}} from top = O_{{m+1}} - R_{{m+2}} = {R_m1_from_top}")
    print(f"Actual R_{{m+1}} = {R[m+1]}")
    print()

    # The defect:
    # R_{m+1}^{actual} = R_{m+1}^{gen} + Δ where Δ = ?
    # β_m = Budget - R_{m+1} = (O[m] - R_m) - R_{m+1}
    #     = (O[m] - R_m^{gen} + 1) - R_{m+1}
    # If R_{m+1} = R_{m+1}^{gen} + something...
    #
    # Actually: Budget^{actual} = O[m] - R_m^{actual} = O[m] - R_m^{gen} + 1
    #         = Budget^{gen} + 1
    # And β_m = Budget^{actual} - R_{m+1}
    # If the complex were fully acyclic (generic): β_m^{gen} = 0, R_{m+1}^{gen} = Budget^{gen}
    # But Budget^{actual} = Budget^{gen} + 1
    # If R_{m+1} = R_{m+1}^{gen} (rank doesn't change from shifting ∂_1):
    #   β_m = Budget^{gen} + 1 - Budget^{gen} = 1
    # This gives β_m^{orb} = 1!
    # But for P_7: β_m^{orb} = 0, so this doesn't work directly.

    print(f"If R_{{m+1}} = R_{{m+1}}^gen: β_m would be {O[m] - R_m_actual - R_m1_gen}")
    print(f"Actual β_m^orb = {O[m] - R_m_actual - R[m+1]}")
    print(f"Deficit: R_{{m+1}} - R_{{m+1}}^gen = {R[m+1] - R_m1_gen}")
    print()

    # CHECK: is the deficit at d=m+1 equal to (-1)^{m+2}?
    deficit_m1 = R[m+1] - R_m1_gen
    print(f"Deficit at d=m+1: {deficit_m1}")
    print(f"(-1)^{{m+2}} = {(-1)**(m+2)}")
    print()

print("=== SUMMARY ===\n")
print("The ∂_1 = 0 in the orbit complex creates a rank-1 shift that propagates:")
print("  R_d^{actual} = R_d^{generic} + (-1)^{d+1} for d=1..m")
print()
print("At d=m+1, the shift also propagates but the amount differs:")
print("  P_7:  deficit = 1 (so β_3^orb = 0)")
print("  P_11: deficit = 0 (so β_5^orb = 1)")
print()
print("For P_7 (m=3): R_{m+1}^{gen} = O_3 - R_3^{gen} = 3 - 2 = 1")
print("  R_{m+1} = 2, so the rank INCREASES by 1 compared to generic")
print("  This means ∂_4 has rank 2, using up all of Budget=2, leaving β_3=0")
print()
print("For P_11 (m=5): R_{m+1}^{gen} = O_5 - R_5^{gen} = 92 - 31 = 61")
print("  R_{m+1} = 61, SAME as generic. So the rank shift stops at d=m.")
print("  β_5 = Budget - R_{m+1} = 62 - 61 = 1")
print()
print("CONCLUSION: For P_7, the shift continues to d=m+1 (alternating +1).")
print("For P_11, the shift STOPS at d=m.")
print("The question: when does the shift stop?")
print()
print("For m=3: shift continues through d=4. Total shift accumulation = 0.")
print("  (shifts: -1, +1, -1, +1 → cancel)")
print("For m=5: shift stops at d=m. Total shift at d=m+1 is 0.")
print("  Budget_actual = Budget_gen + 1, but R_{m+1} = R_{m+1}^gen, so β_m = 1.")
