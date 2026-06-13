#!/usr/bin/env python3
"""
Verify β_9 = 3 for P_19 using the stacking trick at degree m = 9.

Instead of computing Ω_9 (which requires 6M orbit reps),
use the stacking trick at degree m = 9 to directly compute R_9^{stacked}.

R_9 = rank([C_9; B_9]) - rank(C_9)

But this requires boundary B_9 (mapping d=9 to d=8) PLUS constraint C_9.

Alternative: use the FACT that R_d (from bottom) = R_d (actual) for d ≤ m
when all β = 0 for d < m. We already have R_9 from bottom recursion = 268340.

Then β_9 = Ω_9 - R_9 - R_10.
We need both Ω_9 and R_10.

R_10 = rank([C_10; B_10]) - rank(C_10) — this is the rank of ∂_10 on Ω_10.
Or equivalently: R_10 from TOP recursion: R_10 = Ω_10 - R_11 (if β_10 = 0...
but β_10 might be nonzero!)

Actually wait: β_{m+1} = β_m (proved). So β_10 = β_9.
And for d > m+1: β_d = 0 (acyclicity above).

So the top recursion gives:
R_{2m} = Ω_{2m}
R_{2m-1} = Ω_{2m-1} - R_{2m}
...
R_{m+2} = alternating sum of Ω_{m+2},...,Ω_{2m}

And: R_{m+1} = Ω_{m+1} - R_{m+2} - β_{m+1}

This requires Ω for the TOP HALF of the sequence.

For P_19 (m=9, degrees 0..18):
Bottom: Ω_0,...,Ω_8 → R_1,...,R_9
Top: Ω_{11},...,Ω_{18} → R_{18},...,R_{11}
Gap: degrees 9, 10 (= m, m+1)

So we need EITHER:
1. Ω_9 and the full top half, OR
2. R_10 directly from stacking trick at d=10

Option 2: R_10 from stacking trick requires orbit reps at d=10 (even MORE than d=9).

Let me try option 3: compute the top half Ω using DUALITY.

For the Paley tournament, there might be a DUALITY relating Ω_d and Ω_{2m-d}.

P_7 (m=3): Ω_orb = [1, 1, 2, 3, 3, 2, 1]
Ω_d vs Ω_{2m-d}: (1,1), (1,2), (2,3), (3,3), (3,2), (2,1), (1,1)
NOT palindromic.

P_11 (m=5): Ω_orb = [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6]
Ω_d vs Ω_{2m-d}:
  d=0↔10: (1,6)
  d=1↔9: (1,36)
  d=2↔8: (4,90)
  d=3↔7: (14,138)
  d=4↔6: (41,140)
  d=5↔5: (92,92)
Not palindromic either.

Actually, what IS palindromic?
For P_7: |A_d|/m = [1, 1, 3, 7, 13, 15, 9]
For P_11: |A_d|/m = [1, 1, 5, 22, 86, 286, 794, 1747, 2879, 3149, 1729]
These are NOT palindromic either.

Let me check the FULL A_d values.
For P_11: A_orb = [1, 1, 5, 22, 86, 286, 794, 1747, 2879, 3149, 1729]
Hmm, the last value 1729 = 12^3 + 1 = Hardy-Ramanujan number. Coincidence.

Check: Σ (-1)^d |A_d|^{orb} = 1 - 1 + 5 - 22 + 86 - 286 + 794 - 1747 + 2879 - 3149 + 1729
= (1+5+86+794+2879+1729) - (1+22+286+1747+3149)
= 5494 - 5205 = 289 = 17^2.
Hmm, 289 = 17^2. And p = 11. Not obvious connection.

Let me try another approach: use the RATIO Ω_d/Ω_{2m-d}.

P_11:
d=1/d=9: Ω_1/Ω_9 = 1/36
d=2/d=8: Ω_2/Ω_8 = 4/90 = 2/45
d=3/d=7: Ω_3/Ω_7 = 14/138 = 7/69
d=4/d=6: Ω_4/Ω_6 = 41/140

No clean ratios.

OK, let me try a COMPLETELY different approach.

For the FULL complex (not orbit), we can use the JH-J method from earlier sessions
to compute ranks in eigenspaces. But that requires building the boundary matrices.

Actually, we can just compute β_9 from the FORMULA, assuming THM-130 is correct:
β_9 = m(m-3)/2 = 9·6/2 = 27.
And β_9^{orb} = (m-3)/2 = 3.

Let me verify CONSISTENCY: if β_9^{orb} = 3, what do we need?

From bottom: R_9 = 268340
Budget_m = Ω_9 - R_9 = Ω_9 - 268340
β_9 = Budget - R_10

From top: need Ω_{10},...,Ω_{18}. If I can compute these, I can verify.

But computing Ω_d for large d is the bottleneck.

ALTERNATIVE: Can I compute Ω_{2m-d} from Ω_d using some symmetry?

The operation T ↦ T^op (complement) maps QR ↦ NR. For the diff-seq complex,
the complement exchanges QR and NR. But diff-seqs use only QR elements.
The complement duality β(T) = β(T^op) is a theorem (proved earlier).

For P_p: T^op has connection set NR instead of QR. But P_p is self-complementary
(up to relabeling). So β(P_p) = β(P_p^{op}).

But this doesn't give a duality between Ω_d and Ω_{2m-d}.

What about REVERSING diff-seqs? If σ = (s_1,...,s_d), then σ^{rev} = (s_d,...,s_1).
The reversed sequence has partial sums P'_k = Σ_{i=d-k+1}^d s_i = P_d - P_{d-k}.
σ ∈ A_d iff all P_k distinct and nonzero.
σ^{rev} ∈ A_d iff all P'_k distinct and nonzero.
P'_k = P_d - P_{d-k}. Since P_0=0, P'_d = P_d - P_0 = P_d.
P'_k = P_d - P_{d-k}. Distinct iff P_{d-k} distinct iff original P distinct ✓.
P'_k ≠ 0 iff P_d ≠ P_{d-k}, which is true since P_i are all distinct.
BUT we also need P'_k ≠ P_d (wait, we need all P'_k distinct AND nonzero).
P'_0 = P_d - P_d = 0 ← this is the "starting vertex" and must not appear again.
P'_k = P_d - P_{d-k} ≠ 0 iff P_d ≠ P_{d-k}, which is true for k < d. ✓

So σ ∈ A_d ⟺ σ^{rev} ∈ A_d. The reversal map is a bijection on A_d.
|A_d| is invariant under reversal.

But does reversal induce a map on Ω_d? The constraint matrix involves junk faces.
If face_i(σ) is junk, is face_{d-i}(σ^{rev}) also junk?
face_i(σ) = (s_1,...,s_{i-1}+s_i,...,s_d), checking if s_{i-1}+s_i ∈ QR.
face_{d-i}(σ^{rev}) involves the reversed sequence at position d-i.
This might not give the same merged element.

So reversal preserves |A_d| but may not preserve Ω_d. Not a useful symmetry.

CONCLUSION: Without computing Ω_9 (or the top half), we cannot directly verify
β_9 = 3 from the orbit computation alone. The P_19 computation is too large
for current methods.

However, the formula β_m = m(m-3)/2 is STRONGLY supported by:
1. P_7 (m=3): β = 0 ✓
2. P_11 (m=5): β = 5 ✓
3. β_{m+1} = C(m+1,2) and χ = p ✓ for both
4. The eigenspace mechanism (rank shift + budget equality) ✓
5. The combinatorial interpretation (diagonals of m-gon) ✓
"""

print("P_19 β_9 verification status:")
print()

# Known orbit values
omega = [1, 1, 8, 60, 417, 2648, 15140, 76474, 331958]
print(f"Ω_orb = {omega} (d=0..8)")
print()

# Bottom recursion
R = [0, 0]
for d in range(2, len(omega)+1):
    if d-1 < len(omega):
        R.append(omega[d-1] - R[-1])
print(f"R (bottom): {R}")
print()

# Verification status
print(f"R_9 = Ω_8 - R_8 = {omega[8]} - {R[8]} = {R[9]}")
print()
print(f"To compute β_9 = Ω_9 - R_9 - R_10:")
print(f"  R_9 = {R[9]} (known)")
print(f"  Ω_9 = ? (need 6M orbit reps, ~20min+ computation)")
print(f"  R_10 = ? (need top-half Ω values)")
print()
print(f"Predicted (THM-130):")
print(f"  β_9^orb = (m-3)/2 = 3")
print(f"  β_9 = m(m-3)/2 = 27")
print(f"  β_10 = m(m+1)/2 = 45")
print(f"  χ = 1 - 27 + 45 = 19 = p ✓")
print()

# What we CAN verify: the alternating sum from bottom gives R_9 consistently
print("Bottom recursion consistency:")
alt_sum = sum((-1)**(8-d) * omega[d] for d in range(1, 9))
print(f"  R_9 (alt sum) = {alt_sum}")
print(f"  R_9 (recursion) = {R[9]}")
print(f"  Match: {'✓' if alt_sum == R[9] else '✗'}")
print()

# Budget prediction
print("If β_9 = 3:")
print(f"  Budget = R_9 + β_9 + R_10")
print(f"  Budget = Ω_9 - R_9 = Ω_9 - {R[9]}")
print(f"  R_10 = Budget - β_9 = Ω_9 - {R[9]} - 3 = Ω_9 - {R[9]+3}")
print()

# Can we estimate Ω_9 from the growth pattern?
print("Ω growth pattern:")
for d in range(1, len(omega)):
    ratio = omega[d] / omega[d-1] if omega[d-1] > 0 else 0
    print(f"  Ω_{d}/Ω_{d-1} = {omega[d]}/{omega[d-1]} = {ratio:.4f}")

# Extrapolate Ω_9 from the ratio trend
# Ratios: 1.0, 8.0, 7.5, 6.95, 6.35, 5.72, 5.05, 4.34
ratios = [omega[d]/omega[d-1] for d in range(1, len(omega))]
print(f"\nRatios: {[f'{r:.3f}' for r in ratios]}")

# The ratios are roughly decreasing by ~0.6 per step
# Next ratio ≈ 4.34 - 0.7 ≈ 3.6
# Ω_9 ≈ 331958 * 3.6 ≈ 1,195,050
# But the ratio decrease is not linear; let me try quadratic fit

print("\nRatio differences:")
ratio_diffs = [ratios[i+1] - ratios[i] for i in range(len(ratios)-1)]
print(f"  {[f'{d:.3f}' for d in ratio_diffs]}")

# Crude estimate: next diff ≈ same as last
# next_ratio ≈ 4.34 + (-0.71) ≈ 3.63
est_omega_9 = int(omega[8] * 3.63)
print(f"\nCrude estimate: Ω_9 ≈ {omega[8]} × 3.63 ≈ {est_omega_9}")

# From the partial streaming computation: at 3.8M/6M reps,
# rank was 3.245M, so null space = 555K. If this scales to 6M:
# Ω_9 ≈ 6M - rank_final. Rank was growing at about (3.245/3.8) ≈ 85.4% of reps.
# If this rate continues: final rank ≈ 0.854 * 6M = 5.1M
# Ω_9 ≈ 6M - 5.1M = 0.9M
# Actually the rate was decreasing (rank/reps): 86.4% at 3.5M, 85.4% at 3.8M
# Let's extrapolate: at 6M, maybe ~82%: rank ≈ 4.9M, Ω_9 ≈ 1.1M

print(f"\nFrom streaming data: rank/reps was ~85% at 3.8M/6M reps")
print(f"Extrapolating: Ω_9 ≈ 1,000,000 to 1,200,000")

# Compare with P_11 pattern:
# P_11: Ω_9/Ω_8 = 36/90 = 0.40
# P_19 at d=9: Ω_9/Ω_8 = ?/331958
# For P_11, the ratio at d=m=5 was Ω_5/Ω_4 = 92/41 = 2.24, then drops fast
# P_19 at d=m=9 should have similar behavior to P_11 at d=m=5
# Actually, d=9=m for P_19, so we're at the "peak" ratio, not far past it
# P_11 at d=5 (=m): ratio 2.24. P_19 at d=8 (=m-1): ratio 4.34
# P_11 at d=6 (=m+1): ratio 1.52
# So P_19 at d=9 (=m): ratio ≈ 2-4? Hard to extrapolate.

print(f"\nCOMPARISON: P_11 ratio at d=m = {92/41:.3f}, at d=m+1 = {140/92:.3f}")
print(f"P_19 ratio at d=m-1 = {omega[8]/omega[7]:.3f}")
print(f"P_19 ratio at d=m should be lower (approaching peak), estimated 2.5-4")
print(f"Ω_9 estimate: {omega[8]}×3 = {omega[8]*3} to {omega[8]}×4 = {omega[8]*4}")
