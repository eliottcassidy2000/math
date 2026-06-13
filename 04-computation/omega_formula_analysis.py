#!/usr/bin/env python3
"""
Analyze k=0 eigenspace Omega formulas for Paley tournaments.

DISCOVERED: Ω_d = m(m-1)(2m-3)/2 for d=3.

Proof via Gauss sums:
  A_3 = m³ - N_3(0) where N_3(0) = #{s_i ∈ QR : s1+s2+s3 ≡ 0}
  N_3(0) = (p²-1)/8 = m(m+1)/2 (via Gauss sum evaluation)
  junk_rank_3 = 2m(m-1)
  Ω_3 = A_3 - junk_rank_3 = m³ - m(m+1)/2 - 2m(m-1) = m(m-1)(2m-3)/2

Goal: find formula for Ω_4 and higher.

opus-2026-03-13-S71c
"""
from math import comb
import sys
sys.stdout.reconfigure(line_buffering=True)

# Known data
data = {
    7:  [1, 3, 6, 9, 9, 6, 3],
    11: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30],
    19: [1, 9, 72, 540, 3753],
}

print("=" * 70)
print("OMEGA FORMULA VERIFICATION")
print("=" * 70)

for p, omega in sorted(data.items()):
    m = (p - 1) // 2
    print(f"\nP_{p} (m={m}, p={p}):")
    print(f"  Omega = {omega}")

    # Verify formulas
    formulas = []
    for d, od in enumerate(omega):
        f = {}
        f['d'] = d
        f['omega'] = od

        if d == 0:
            f['formula'] = 1
            f['name'] = '1'
        elif d == 1:
            f['formula'] = m
            f['name'] = 'm'
        elif d == 2:
            f['formula'] = m * (m - 1)
            f['name'] = 'm(m-1)'
        elif d == 3:
            f['formula'] = m * (m - 1) * (2 * m - 3) // 2
            f['name'] = 'm(m-1)(2m-3)/2'
        else:
            f['formula'] = None
            f['name'] = '?'

        match = (f['formula'] == od) if f['formula'] is not None else None
        formulas.append(f)
        status = '✓' if match else ('✗' if match is not None else '?')
        print(f"  d={d}: Ω={od:>8}, formula={f['name']:>20} = {str(f['formula']):>8} {status}")

# Prove N_3(0) via Gauss sums
print("\n" + "=" * 70)
print("N_3(0) = #{(s1,s2,s3) ∈ QR³ : s1+s2+s3 ≡ 0 mod p}")
print("=" * 70)

print("""
PROOF (via Gauss sums):

Let χ = Legendre symbol mod p, G = Gauss sum = Σ_{a≠0} χ(a) ζ^a.
For p ≡ 3 mod 4: G² = -p.

Σ_{s∈QR} ζ^{ts} = (χ(t)·G - 1) / 2

N_3(0) = (1/p) Σ_t (Σ_{s∈QR} ζ^{ts})³
       = m³/p + (1/8p) Σ_{t≠0} (χ(t)G - 1)³

Expanding:
  (χ(t)G - 1)³ = χ(t)³G³ - 3χ(t)²G² + 3χ(t)G - 1

Summing over t ≠ 0:
  Σ χ(t)³ = Σ χ(t) = 0     (since χ³ = χ)
  Σ χ(t)² = Σ 1 = p-1       (since χ² = 1)
  Σ χ(t) = 0

So Σ_{t≠0} (χ(t)G-1)³ = -3G²(p-1) - (p-1) = (p-1)(3p-1)

N_3(0) = m³/p + (p-1)(3p-1)/(8p)
       = [(p-1)³ + (p-1)(3p-1)] / (8p)
       = (p-1)(p²+p) / (8p)
       = (p²-1)/8
       = m(m+1)/2  ✓  (since m = (p-1)/2)
""")

# Verify computationally
print("Computational verification of N_3(0):")
for p in [7, 11, 19, 23, 31, 43]:
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    if p % 4 != 3:
        continue
    m = len(qr)

    count = 0
    for s1 in qr:
        for s2 in qr:
            for s3 in qr:
                if (s1 + s2 + s3) % p == 0:
                    count += 1

    expected = m * (m + 1) // 2
    print(f"  P_{p} (m={m}): N_3(0) = {count}, m(m+1)/2 = {expected}, match = {count == expected}")

# Now analyze Ω_4
print("\n" + "=" * 70)
print("OMEGA_4 ANALYSIS")
print("=" * 70)

# A_4 = m⁴ minus sum-constraint failures (PIE)
# Need to subtract sequences where partial sums collide or hit 0
# Constraints: s1+s2+s3 ≠ 0, s1+s2+s3+s4 ≠ 0, s2+s3+s4 ≠ 0
# (All other distinctness conditions are automatic from -1 ∉ QR)

print("""
A_4 = #{(s1,s2,s3,s4) ∈ QR⁴ : all partial sums distinct}

Automatic conditions (from -1 ∉ QR for p ≡ 3 mod 4):
  s1 ≠ 0, s1+s2 ≠ 0, s2 ≠ 0, s3 ≠ 0, s4 ≠ 0, s3+s4 ≠ 0

Non-automatic conditions:
  (i)   s1+s2+s3 ≠ 0
  (ii)  s1+s2+s3+s4 ≠ 0
  (iii) s2+s3+s4 ≠ 0

By PIE: A_4 = m⁴ - |S_i ∪ S_ii ∪ S_iii| where S_x = violations of condition x.
""")

for p in [7, 11, 19, 23, 31, 43]:
    if p % 4 != 3:
        continue
    qr = sorted(set((a*a)%p for a in range(1,p)))
    m = len(qr)
    qr_set = set(qr)

    # Count A_4 directly
    count_A4 = 0
    count_S1 = 0  # s1+s2+s3 ≡ 0
    count_S2 = 0  # s1+s2+s3+s4 ≡ 0
    count_S3 = 0  # s2+s3+s4 ≡ 0
    count_S12 = 0  # both (i) and (ii)
    count_S13 = 0  # both (i) and (iii)
    count_S23 = 0  # both (ii) and (iii)
    count_S123 = 0  # all three

    for s1 in qr:
        for s2 in qr:
            t2 = (s1 + s2) % p
            for s3 in qr:
                t3 = (t2 + s3) % p
                b1 = (t3 == 0)
                for s4 in qr:
                    t4 = (t3 + s4) % p
                    b2 = (t4 == 0)
                    t_234 = (s2 + s3 + s4) % p
                    b3 = (t_234 == 0)

                    # Also check s1+s2+s3+s4 ≠ s1 (i.e., s2+s3+s4 ≠ 0) => covered by b3
                    # And s1+s2+s3+s4 ≠ s1+s2 (i.e., s3+s4 ≠ 0) => automatic
                    # And s1+s2+s3 ≠ s1 (i.e., s2+s3 ≠ 0) => automatic

                    has_bad = b1 or b2 or b3
                    if not has_bad:
                        count_A4 += 1
                    if b1: count_S1 += 1
                    if b2: count_S2 += 1
                    if b3: count_S3 += 1
                    if b1 and b2: count_S12 += 1
                    if b1 and b3: count_S13 += 1
                    if b2 and b3: count_S23 += 1
                    if b1 and b2 and b3: count_S123 += 1

    by_pie = m**4 - (count_S1 + count_S2 + count_S3) + (count_S12 + count_S13 + count_S23) - count_S123
    n30 = m * (m + 1) // 2

    print(f"\nP_{p} (m={m}):")
    print(f"  m⁴ = {m**4}")
    print(f"  |S_i| = {count_S1} = N_3(0)·m = {n30}·{m} = {n30*m} {'✓' if count_S1 == n30*m else '✗'}")
    print(f"  |S_ii| = {count_S2}")
    print(f"  |S_iii| = {count_S3} = m·N_3(0) = {m*n30} {'✓' if count_S3 == m*n30 else '✗'}")
    print(f"  |S_i∩S_ii| = {count_S12}")
    print(f"  |S_i∩S_iii| = {count_S13}")
    print(f"  |S_ii∩S_iii| = {count_S23}")
    print(f"  |S_i∩S_ii∩S_iii| = {count_S123}")
    print(f"  A_4 = {count_A4}, by PIE = {by_pie} {'✓' if count_A4 == by_pie else '✗'}")
    print(f"  A_4/m = {count_A4/m:.1f}")

    # Now what's the formula for S_ii = N_4(0)?
    # N_4(0) = #{(s1,...,s4) ∈ QR⁴ : s1+...+s4 ≡ 0}
    # By Gauss sums: similar to N_3(0) but with 4th power
    n40_expected = (m**4 + (p-1)*(3*p-1)*(m-1) + 0) // p  # guess
    # Actually let's compute via Gauss sums
    # N_4(0) = m⁴/p + (1/16p) Σ_{t≠0} (χ(t)G-1)⁴
    # (χ(t)G-1)⁴ = χ⁴G⁴ - 4χ³G³ + 6χ²G² - 4χG + 1
    #             = G⁴ - 4χG³ + 6G² - 4χG + 1
    # Σ_t: G⁴(p-1) + 0 + 6G²(p-1) + 0 + (p-1)
    # = (p-1)[G⁴ + 6G² + 1]
    # G² = -p, G⁴ = p²
    # = (p-1)[p² - 6p + 1]
    n40_gauss = (m**4 * 1 + (p-1)*(p*p - 6*p + 1) // 16) // 1  # need to be careful
    # Actually: N_4(0) = m⁴/p + (p-1)(p²-6p+1)/(16p)
    n40_exact = m**4 // 1  # placeholder
    # Let me just compute it numerically
    print(f"  N_4(0) = {count_S2}")
    # Check formula: (p-1)(p³-9p+16) / (16p)?
    # Numerator check:
    # N_4(0) = [m⁴·16 + (p-1)(p²-6p+1)] / (16p)
    # = [(p-1)⁴/16·16 + (p-1)(p²-6p+1)] / (16p)... this is getting messy
    # Let me just verify the Gauss sum formula directly
    numerator = m**4 * 16 + (p-1)*(p*p - 6*p + 1)
    if numerator % (16) == 0:
        n40_formula = numerator // 16
        # But wait, we divided by p already in the formula
        # N_4(0) = [16·m⁴ + (p-1)(p²-6p+1)] / (16p) doesn't work if not divisible
        # Actually N_4(0) = m⁴/p + (p-1)(p²-6p+1)/(16p)
        # = [16m⁴ + (p-1)(p²-6p+1)] / (16p)
        num = 16 * m**4 + (p-1)*(p*p - 6*p + 1)
        if num % (16*p) == 0:
            n40_f = num // (16*p)
            print(f"  N_4(0) formula = {n40_f} {'✓' if n40_f == count_S2 else '✗'}")
        else:
            print(f"  N_4(0) formula: {num}/{16*p} = {num/(16*p):.4f} (not integer?)")

    if p >= 31:
        break  # skip large primes to save time
