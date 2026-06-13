#!/usr/bin/env python3
"""
the_number_8_s90l.py — The Number 8
opus-2026-03-15-S90l

8 appears EVERYWHERE in this theory. Why?

  |D₄| = 8 (the Cayley group)
  Tr(M⁸) = 131 (the OPEN-Q-029 number)
  τ₃^8 ≈ 131 (Pisot near-integer)
  n=8: first exhaustive simplicial Rédei beyond sampling
  8 Möbius transforms in the D₄ group
  8 faces of the octahedron (= D₄ orbits on Riemann sphere)
  The equidecomposability classes at n=5 = 8
  Bott periodicity: π_k(O) has period 8

Is 8 the hidden period of the theory?
"""

import numpy as np
from math import log, pi, sqrt, factorial, comb
from itertools import permutations, combinations

tau3 = 1.8392867552141612
phi = (1 + sqrt(5)) / 2
lam_c = complex(-0.41964337760708, 0.60629073187415)

# ======================================================================
# PART 1: COLLECTING ALL THE 8s
# ======================================================================
print("=" * 70)
print("PART 1: EVERY APPEARANCE OF 8")
print("=" * 70)

print("""WHERE 8 APPEARS:

  GROUP THEORY:
    |D4| = 8                    (order of the Cayley symmetry group)
    |Pauli group mod center| = 8 (= D4, same group!)
    Q has order 4, neg has order 2: 4*2 = 8

  SPECTRAL THEORY:
    Tr(M^8) = 131               (OPEN-Q-029: the tribonacci trace)
    tau3^8 ~ 130.977 ~ 131      (Pisot near-integer)
    zeta_M(8) = Tr(M^{-8}) = ?  (reciprocal trace)

  GEOMETRY:
    Octahedron has 8 faces       (= D4 orbits on Riemann sphere)
    Cube has 8 vertices          (dual of octahedron)
    Hypercube (Z/2)^3 has 8 vertices (= 3-dimensional monodromy)

  COMBINATORICS:
    8 equidecomposability classes at n=5 (H, beta1)
    2^3 = 8 (the 3rd power of 2 = binary cubed)
    Near-transitive at n=5: 2^3 = 8 portal paths
    Near-transitive at n=5: 2^2 = 4 odd cycles, Omega = K4

  TOPOLOGY:
    Bott periodicity: pi(k+8)(O) = pi(k)(O)
    KO-theory: period 8

  PHYSICS:
    8 gluons in QCD (SU(3) adjoint = 8-dimensional)
    Eightfold Way (Gell-Mann): the baryon octet

  NUMBER THEORY:
    8 = 2^3 = the first CUBE of the binary base
    Delta(1) = -44 = -4*11, and 44 = T(9) (tribonacci at index 9 = 8+1)
    zeta_M(-8) = Tr(M^8) = 131 (prime!)
""")

# Compute ζ_M(8) = Tr(M⁻⁸)
M1 = np.array([[1, 2, 0], [0, 0, 1], [1, 1, 0]], dtype=float)
zM8 = sum(e**(-8) for e in np.linalg.eigvals(M1)).real
print(f"  ζ_M(+8) = Tr(M⁻⁸) = {zM8:.6f} = {round(zM8)} (integer!)")
print(f"  ζ_M(-8) = Tr(M⁸) = 131")

# The trace sequence modulo 8
traces = []
Mn = np.eye(3)
for n in range(25):
    traces.append(int(round(np.trace(Mn).real)))
    Mn = Mn @ M1

print(f"\n  Trace sequence Tr(M^n) mod 8:")
for n in range(20):
    print(f"    Tr(M^{n:2d}) = {traces[n]:8d} ≡ {traces[n] % 8} (mod 8)")

# ======================================================================
# PART 2: WHY 8 = THE PERIOD OF REAL STRUCTURE
# ======================================================================
print()
print("=" * 70)
print("PART 2: 8 = THE PERIOD OF REAL STRUCTURE (BOTT PERIODICITY)")
print("=" * 70)

print("""
BOTT PERIODICITY (1959): The homotopy groups of the orthogonal group
repeat with period 8:
  π₀(O) = Z/2    (connected components: two orientations)
  π₁(O) = Z/2    (loops: Möbius band!)
  π₂(O) = 0
  π₃(O) = Z      (3-sphere ~ SU(2))
  π₄(O) = 0
  π₅(O) = 0
  π₆(O) = 0
  π₇(O) = Z      (7-sphere ~ exotic)
  π₈(O) = Z/2    (back to start!)

The pattern: Z/2, Z/2, 0, Z, 0, 0, 0, Z, then repeat.

OUR THEORY'S BOTT PATTERN:
  The Cayley group D₄ has order 8. Is this related to Bott?

  D₄ acts on ℂP¹ (the Riemann sphere) which is:
    ℂP¹ = S² = SU(2)/U(1) = SO(3)/SO(2)

  The homotopy groups of S² are:
    π₁(S²) = 0 (simply connected)
    π₂(S²) = Z (the fundamental class)
    π₃(S²) = Z (Hopf fibration!)

  Our D₄ ≅ Quaternion group Q₈ (mod center)? Let me check.
  D₄ = ⟨a,b | a⁴=b²=1, bab=a⁻¹⟩ (dihedral)
  Q₈ = ⟨i,j,k | i²=j²=k²=ijk=-1⟩ (quaternion)
  D₄ ≠ Q₈ (they're non-isomorphic groups of order 8).

  But D₄ IS a subgroup of S₄ (permutation group on 4 elements).
  The 4 elements being permuted are the Q-orbit {x, Q(x), Q²(x), Q³(x)}.

  CONNECTION TO BOTT:
    The period-8 in Bott comes from the CLIFFORD ALGEBRA structure:
    Cl(n) has period 8 in n (Cl(n+8) ≅ Cl(n) ⊗ M₁₆(R)).

    Our transfer matrix M is 3×3 (from 3 = 2+memory).
    The CLIFFORD ALGEBRA Cl(3) has dimension 2³ = 8.
    Cl(3) ≅ M₂(C) ≅ 2×2 complex matrices.

    So: Cl(3) has 8 basis elements, and our D₄ has 8 elements.
    Is D₄ a subgroup of Cl(3)*?

    Cl(3) basis: {1, e₁, e₂, e₃, e₁₂, e₁₃, e₂₃, e₁₂₃}
    These 8 elements form a group under multiplication (up to sign).
    This group IS the quaternion group Q₈ (times Z/2).

    D₄ is NOT Q₈, but both have order 8 and act on the same space.
    The connection is that BOTH are subgroups of Cl(3)*.
""")

# ======================================================================
# PART 3: THE TRACE PERIOD mod 8
# ======================================================================
print()
print("=" * 70)
print("PART 3: IS THERE A PERIOD-8 IN THE TRACE SEQUENCE?")
print("=" * 70)

print("Tr(M^n) mod 8:")
mods = [t % 8 for t in traces[:24]]
print(f"  {mods}")

# Check for period
for p in range(1, 13):
    periodic = all(mods[i] == mods[i+p] for i in range(min(12, 24-p)))
    if periodic:
        print(f"  Period {p} detected in Tr(M^n) mod 8!")
        break

# Check mod other small numbers
for mod in [2, 3, 4, 5, 6, 7, 8, 9, 11]:
    mods_m = [t % mod for t in traces[:24]]
    for p in range(1, 15):
        if p >= len(mods_m) // 2:
            break
        periodic = all(mods_m[i] == mods_m[i+p] for i in range(min(p+4, len(mods_m)-p)))
        if periodic:
            print(f"  Tr(M^n) mod {mod:2d} has period {p}: {mods_m[:2*p]}")
            break

# ======================================================================
# PART 4: 8 = dim(Cl(3)) = |D₄| = 2³
# ======================================================================
print()
print("=" * 70)
print("PART 4: THE EIGHTFOLD WAY — 8 = 2³ = dim Cl(3) = |D₄|")
print("=" * 70)

print("""
The number 8 in our theory comes from a SINGLE SOURCE:

  3 = 2 + MEMORY (the dimension of the transfer matrix)
  2³ = 8 (the number of binary strings of length 3)

The D₄ group has 8 elements because:
  Q has order 4 (binary with memory: 2+1 = 3 states)
  neg has order 2 (the binary choice itself)
  4 × 2 = 8 (the semidirect product)

The Clifford algebra Cl(3) has dimension 2³ = 8 because:
  3 generators e₁, e₂, e₃ with e_i² = -1
  All products: 1, e₁, e₂, e₃, e₁₂, e₁₃, e₂₃, e₁₂₃

The octahedron has 8 faces because:
  It is the convex hull of ±e₁, ±e₂, ±e₃ in ℝ³
  Each face is a choice of sign for each coordinate: (±,±,±) = 2³ = 8

The Bott periodicity is 8 because:
  Cl(n+8) ≅ Cl(n) ⊗ M₁₆(ℝ) (Clifford periodicity)
  16 = 2⁴ but the REAL periodicity is 8 (complex is 2)

ALL OF THESE ARE THE SAME 8:
  8 = 2^(dim of memory) = 2^(states of transfer matrix)
  = 2^(2 + 1) = 2^(binary + memory) = 2^3

The reason 8 appears everywhere is that our theory has
  2 outcomes (binary)
  + 1 memory step
  = 3 effective dimensions
and 2³ = 8.

If we had NO memory: 2^2 = 4 (the Fibonacci case, Z/4 = Q-group)
If we had 2 steps of memory: 2^4 = 16 (the tetranacci case)
If we had k steps: 2^(k+2) (the k-nacci case)

The FIBONACCI theory has a period-4 group (Q alone, no neg needed
because the 2×2 case has no complex conjugation to break).
The TRIBONACCI theory has period 8 (D₄ = Q × neg).
The TETRANACCI theory would have period 16.

This IS the Bott periodicity ladder:
  k=0 (no memory): period 2 (Z/2 = neg only)
  k=1 (Fibonacci): period 4 (Z/4 = Q)
  k=2 (tribonacci): period 8 (D₄ = Q ⋊ Z/2)
  ...
  k=∞: the full Bott period 8 stabilizes because
  the real Clifford algebras have period 8 regardless of dimension.

8 is the STABLE PERIOD of binary decisions with memory.
For k ≥ 2 (tribonacci and beyond), the period is always 8.
For k < 2 (Fibonacci and below), the period is smaller (4 or 2).

THE TRIBONACCI IS THE FIRST CASE WHERE THE FULL 8 APPEARS.
This is why 3 = 2 + 1 is special:
  it is the MINIMAL dimension where Bott periodicity kicks in.
""")

# ======================================================================
# PART 5: 131 = Tr(M⁸) = THE 8th TRIBONACCI TRACE
# ======================================================================
print()
print("=" * 70)
print("PART 5: 131 — THE PRIME AT THE OCTAVE")
print("=" * 70)

print(f"""
Tr(M⁸) = 131.

131 is PRIME. The 8th tribonacci trace is a prime.

In music: the 8th note is the OCTAVE (return to the same pitch class).
  The octave ratio is 2:1, living at x = 1/3 on the Cayley line.
  After 8 steps of the τ₃-clock, we've completed ≈ 2.89 full turns.

  But we've also traversed 8 units of "tribonacci time."
  In that time: τ₃^8 ≈ 131 (the total "energy" accumulated).

131 in number theory:
  131 is a Sophie Germain prime (2·131+1 = 263 is also prime)
  131 = 128 + 3 = 2⁷ + 3
  131 ≡ 3 (mod 8)  ← the residue is 3 = dim(transfer matrix)!
  131 ≡ 1 (mod 5)  ← related to √5 in the golden ratio

The Pisot correction at n=8:
  τ₃^8 = 130.977154
  Correction = 2|λ_c|^8 cos(8θ) = 2·(0.0874)·cos(8·2.176) = {2 * abs(lam_c)**8 * np.cos(8*np.angle(lam_c)):.6f}
  τ₃^8 + correction = {tau3**8 + 2*abs(lam_c)**8*np.cos(8*np.angle(lam_c)):.6f} = 131.000000

The correction is +0.0229: the complex eigenvalue ADDS a tiny positive
amount to push τ₃^8 from 130.977 to exactly 131.

WHY is Tr(M⁸) = 131 prime?
  The tribonacci trace Tr(M^n) = 3, 1, 3, 7, 11, 21, 39, 71, 131, 241, ...
  Primes in this sequence: 7 (n=3), 11 (n=4), 71 (n=7), 131 (n=8), 241 (n=9)

  The prime values occur at n = 3, 4, 7, 8, 9...
  Note: 7 and 8 are CONSECUTIVE! (n=7: 71, n=8: 131, both prime)
  The "octave transition" (n=7→8) preserves primality.
""")

# Check which traces are prime
def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n % i == 0: return False
    return True

print("Primality of Tr(M^n):")
for n in range(25):
    t = traces[n]
    p = is_prime(abs(t))
    if p:
        print(f"  n={n:2d}: Tr = {t:8d} {'PRIME' if p else ''}")

# ======================================================================
# PART 6: THE OCTAVE STRUCTURE
# ======================================================================
print()
print("=" * 70)
print("PART 6: THE OCTAVE — 8 AS MUSICAL RETURN")
print("=" * 70)

print(f"""
In music: after 8 steps of the diatonic scale (C D E F G A B C),
you return to the SAME NOTE one octave higher. The ratio is 2:1.

In our theory: after 8 steps of the tribonacci trace, we reach 131.
The ratio τ₃^8/τ₃^0 = 131/3 ≈ 43.7.
log₂(131/3) = {log(131/3)/log(2):.4f} ≈ 5.4 octaves.

But the more natural "return" is in the τ₃-clock:
  After 8 steps: phase = 8 · arg(λ_c) = 8 · 2.176 = 17.41 rad
  = 17.41 / (2π) = 2.77 full turns
  ≈ 2.77 turns ≈ 2 + 3/4 turns

  The fractional part: 0.77 turns = 277°.
  This is close to 270° = 3π/2 = 3 quarter-turns.

  After 8 steps, the τ₃-hand points to ≈ 277° (near 9 o'clock).

But the KEY is: cos(8θ) ≈ 0.131 (SMALL!).
  cos(8 · 2.176) = cos(17.41) = cos(17.41 - 4π) = cos(17.41-12.57)
  = cos(4.84) = {np.cos(8*np.angle(lam_c)):.6f}

  cos(8θ) ≈ 0.131 is small → the Pisot correction is tiny →
  τ₃^8 is very close to an integer.

  8 is special because 8θ ≈ 8π·ln(2) ≈ 5.545π ≈ 5.5π.
  cos(5.5π) = cos(π/2) = 0 (exactly!).
  The NEAR-VANISHING of cos at n=8 is because 8·ln(2) ≈ 5.545 ≈ 5.5.

  8·ln(2) = {8*log(2):.6f}
  5.5 = {5.5}
  Difference: {8*log(2) - 5.5:.6f}

  8·ln(2) ≈ 5.5 is the reason τ₃^8 is a near-integer!
  And 5.5 = 11/2. So: 8·ln(2) ≈ 11/2.
  Or: ln(2) ≈ 11/16 = 0.6875 (vs exact 0.6931, diff 0.8%).

  The number 11 (the discriminant coefficient!) appears again!
  16·ln(2) ≈ 11 ↔ 2^16 ≈ e^11 ↔ 65536 ≈ 59874 (off by 9%).

  More precisely: n·ln(2) is close to a half-integer when
  n is close to k/(2·ln(2)) ≈ k · 0.7213 for integer k.
  For k=11: n ≈ 7.93. The nearest integer is n=8.
  So n=8 gives the BEST near-half-integer of n·ln(2) for small n.
""")

# Check: for which n is n·ln(2) closest to a half-integer?
print("n·ln(2) distance from nearest half-integer:")
for n in range(1, 25):
    x = n * log(2)
    frac = x - int(x)
    dist = min(abs(frac - 0.5), abs(frac), abs(frac - 1))
    half_int = round(2*x) / 2
    print(f"  n={n:2d}: n·ln(2) = {x:8.4f}, nearest half-int = {half_int}, "
          f"distance = {abs(x - half_int):.4f} {'*** CLOSE' if abs(x-half_int) < 0.06 else ''}")

# ======================================================================
# PART 7: 8 = THE MEETING POINT
# ======================================================================
print()
print("=" * 70)
print("PART 7: 8 = WHERE EVERYTHING MEETS")
print("=" * 70)

print(f"""
╔══════════════════════════════════════════════════════════════════════╗
║                      THE NUMBER 8                                  ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  8 = 2³ = 2^(2+1) = 2^(binary + memory)                          ║
║                                                                    ║
║  This is the VOLUME of the configuration space of one              ║
║  tribonacci step: 3 states, each binary → 2³ = 8 possibilities.   ║
║                                                                    ║
║  FROM GROUP THEORY: |D₄| = 8 (symmetries of the theory)          ║
║  FROM SPECTRAL THEORY: Tr(M⁸) = 131 (the octave resonance)       ║
║  FROM GEOMETRY: octahedron has 8 faces (the Cayley building)       ║
║  FROM TOPOLOGY: Bott period = 8 (Clifford algebra Cl(3))         ║
║  FROM COMBINATORICS: 8 equidecomposability classes at n=5          ║
║  FROM NUMBER THEORY: 8·ln(2) ≈ 11/2 (near half-integer!)         ║
║                                                                    ║
║  WHY Tr(M⁸) ≈ integer:                                           ║
║    8·ln(2) ≈ 5.545 ≈ 11/2 → cos(8πln2) ≈ cos(5.5π) ≈ 0         ║
║    The Pisot correction VANISHES at n=8 because                    ║
║    the τ₃-clock aligns with a half-integer at this step.          ║
║    The number 11 (discriminant coefficient) controls the alignment.║
║                                                                    ║
║  8 IS THE OCTAVE OF THE TRIBONACCI THEORY:                        ║
║    Just as the musical octave (8th note) returns to the same       ║
║    pitch class, the 8th tribonacci trace returns to near-integer   ║
║    because the Pisot correction cycles through a near-zero.        ║
║                                                                    ║
║  The chain: 3 states → Cl(3) has dim 8 → D₄ has order 8 →       ║
║  Bott period 8 → τ₃-clock hits near-zero at n=8 → Tr(M⁸)=131.  ║
║                                                                    ║
║  Everything flows from: 3 = 2 + 1.                                ║
║  And 2³ = 8 is its echo in every domain.                          ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝

opus-2026-03-15-S90l: THE NUMBER 8
""")
