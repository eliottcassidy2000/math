#!/usr/bin/env python3
"""
moonshine_tournament_s123.py — opus-2026-03-21-S123

MONSTROUS MOONSHINE AND TOURNAMENT THEORY

The "+1" pattern:
  196884 = 196883 + 1  (moonshine: j-coefficient = Monster rep dim + 1)
  1729   = 1728   + 1  (tournaments: taxicab = j(i) + Rédei quantum)

Both are cases where a MODULAR QUANTITY exceeds an ALGEBRAIC QUANTITY by 1.
In moonshine: the 1 is the trivial representation.
In tournaments: the 1 is the Rédei quantum g(2) = 1.

THIS SESSION: explore the moonshine-tournament parallel creatively.
Not claiming tournaments are related to the Monster — but the STRUCTURAL
PATTERN of "modular exceeds algebraic by 1" may have a common origin.

Author: opus-2026-03-21-S123
"""

import sys
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  MONSTROUS MOONSHINE AND TOURNAMENT THEORY")
print("  opus-2026-03-21-S123")
print("=" * 72)

# ========================================================================
# PART 1: The "+1" pattern
# ========================================================================
print("\n" + "=" * 72)
print('  PART 1: THE "+1" PATTERN')
print("=" * 72)
print(f"""
  MONSTROUS MOONSHINE:
    j(q) = q^(-1) + 744 + 196884q + 21493760q² + ...

    196884 = 196883 + 1
    where 196883 = dim of smallest faithful Monster rep
    and 1 = dim of trivial rep

    The j-function coefficient DECOMPOSES as a sum of Monster
    representation dimensions, with the trivial rep adding "+1."

  TOURNAMENT THEORY:
    1729 = 1728 + 1
    where 1728 = j(i) = 12³ = modular j-invariant at complement
    and 1 = g(2) = Rédei quantum

    The taxicab number DECOMPOSES as a modular invariant plus
    the tournament quantum, with the Rédei parity adding "+1."

  THE PARALLEL:
    Both involve the j-function (j-coefficient in moonshine, j(i) in tournaments).
    Both have a "+1" from a TRIVIAL/UNIVERSAL contribution.
    Both connect algebraic structure to modular forms.

  But the SCALE is vastly different:
    Monster: |M| ≈ 8×10^53, dim(V_1) = 196884
    Tournaments: |S_n| = n!, dim = C(n,2)

  The connection is at the STRUCTURAL level, not the numerical level.
""")

# ========================================================================
# PART 2: j-function coefficients and tournament numbers
# ========================================================================
print("=" * 72)
print("  PART 2: j-FUNCTION COEFFICIENTS")
print("=" * 72)

# j(q) = q^{-1} + 744 + 196884q + 21493760q^2 + 864299970q^3 + ...
j_coeffs = [1, 744, 196884, 21493760, 864299970, 20245856256]

print(f"  j(q) = 1/q + 744 + 196884q + 21493760q² + ...")
print(f"\n  j-coefficients: {j_coeffs}")

# Check divisibility by tournament numbers
tournament_nums = {
    3: "cycle atom",
    7: "forbidden atom",
    21: "C(7,2) = arcs at n=7",
    42: "Hurwitz = 2×3×7",
    12: "T(5) = iso classes at n=5",
    1729: "taxicab",
}

print(f"\n  Divisibility of j-coefficients by tournament numbers:")
print(f"  {'coeff':>12s} {'value':>12s} | {'mod 3':>6s} {'mod 7':>6s} {'mod 21':>6s} {'mod 42':>6s} {'mod 12':>6s}")
for i, c in enumerate(j_coeffs):
    name = f"c_{i-1}" if i > 0 else "c_{-1}"
    print(f"  {name:>12s} {c:>12d} | {c%3:>6d} {c%7:>6d} {c%21:>6d} {c%42:>6d} {c%12:>6d}")

print(f"""
  OBSERVATION: 744 = 8 × 93 = 8 × 3 × 31. Divisible by 3 but not 7.
  196884 = 196884. 196884 mod 7 = {196884 % 7}. Divisible by 7? {196884 % 7 == 0}.
  196884 mod 42 = {196884 % 42}. Divisible by 42? {196884 % 42 == 0}.

  196884 = 4 × 49221 = 4 × 3 × 16407 = 12 × 16407.
  So 196884 = 12 × 16407. Divisible by 12 = T(5)!
  And 196884 / 12 = 16407 = 3 × 5469 = 3 × 3 × 1823.
  So 196884 = 36 × 5469.

  196884 mod 1729 = {196884 % 1729}. So 196884 = 113 × 1729 + {196884 - 113*1729}.
  Actually: 113 × 1729 = {113*1729}. 114 × 1729 = {114*1729}. So 196884/1729 ≈ {196884/1729:.2f}.
  Not an integer multiple.
""")

# ========================================================================
# PART 3: The "+1" as a universal pattern in representation theory
# ========================================================================
print("=" * 72)
print('  PART 3: THE "+1" AS UNIVERSAL PATTERN')
print("=" * 72)
print(f"""
  The "+1" pattern appears throughout mathematics:

  1. MOONSHINE: 196884 = 196883 + 1
     (j-coefficient = Monster dimension + trivial rep)

  2. TOURNAMENTS: 1729 = 1728 + 1
     (taxicab = j(i) + Rédei quantum)

  3. RÉDEI: H(T) = 1 + 2K for some integer K ≥ 0
     (HP count = 1 + even part — the "+1" is the base path)

  4. OCF: H = I(Ω, 2) = 1 + 2α₁ + 4α₂ + ...
     (independence polynomial starts with 1 = empty set contribution)

  5. EULER: e^(iπ) + 1 = 0
     (the "+1" bridges exponential and additive worlds)

  6. FERMAT: a^p ≡ a (mod p) → a^(p-1) ≡ 1 (mod p)
     (the "+1" from dividing out the trivial orbit)

  THE COMMON THREAD: In each case, the "+1" represents the
  TRIVIAL CONTRIBUTION that must be added to the NON-TRIVIAL part.

  In moonshine: 196883 is the smallest non-trivial Monster rep.
  The "+1" is the trivial rep ρ₁. Together: the j-coefficient
  decomposes into group actions.

  In tournaments: 1728 is the j-invariant (the "modular background").
  The "+1" is the Rédei quantum (the "tournament departure").
  Together: 1729 is the tournament's modular address.

  The STRUCTURAL PATTERN: when a modular object (j-function, cusp form)
  meets an algebraic object (representation, tournament), the relationship
  is always "modular = algebraic + trivial" — the "+1" quantifies the gap
  between continuous (modular) and discrete (algebraic) descriptions.
""")

# ========================================================================
# PART 4: The moonshine module and the tournament module
# ========================================================================
print("=" * 72)
print("  PART 4: THE MOONSHINE MODULE V♮ AND THE TOURNAMENT MODULE")
print("=" * 72)
print(f"""
  MOONSHINE MODULE V♮:
  - A vertex operator algebra (VOA)
  - Aut(V♮) = Monster group M
  - Graded: V♮ = ⊕_n V_n, with dim(V_n) = j-coefficient c_n
  - V₋₁ = trivial (dim 1), V₁ = ρ₁ ⊕ ρ₁₉₆₈₈₃ (dim 196884)

  THE TOURNAMENT MODULE (speculative):
  What would a "tournament VOA" look like?

  Consider: for each n, the space of tournament invariants.
  - T_n = {{all tournament invariants at order n}} is a vector space
  - dim(T_n) depends on the number of iso classes
  - The "grading" is by tournament order: T = ⊕_n T_n

  The GENERATING FUNCTION of dim(T_n):
    Sum dim(T_n) q^n = Sum (iso classes at n) q^n
    = 1 + 1 + 1 + 2 + 4 + 12 + 56 + 456 + ... (A000568)

  Is there a modular form lurking here?

  A000568: 1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, ...
  The growth rate is roughly 2^(C(n,2))/n!.

  Compare with j-coefficients: 1, 744, 196884, 21493760, ...
  These grow MUCH faster (exponentially in n).

  The tournament iso class sequence is NOT a modular form.
  But it might be a COMPONENT of a modular form — just as
  the Monster rep dimensions are components of j-coefficients.

  THE MOONSHINE DECOMPOSITION:
    c_n = sum_i a(n,i) * dim(rho_i)
  where rho_i are Monster irreps and a(n,i) are non-negative integers.

  A TOURNAMENT DECOMPOSITION (hypothetical):
    dim(T_n) = sum_k b(n,k) * f(k)
  where f(k) is some modular quantity and b(n,k) are structural constants.

  This is highly speculative. But the PATTERN — algebraic structure
  decomposing into modular coefficients — is the same.
""")

# ========================================================================
# PART 5: The 744 constant and tournament theory
# ========================================================================
print("=" * 72)
print("  PART 5: THE CONSTANT 744")
print("=" * 72)

# j(q) = 1/q + 744 + ... The constant 744 in the j-function
# is the unique value that makes j have no constant term in the
# normalized Hauptmodul.

print(f"""
  The constant term of j is 744 = 8 × 93 = 8 × 3 × 31.
  744 = 2³ × 3 × 31.

  In tournament theory:
  744 / 42 = {744/42:.4f} (not clean)
  744 mod 7 = {744 % 7}
  744 mod 21 = {744 % 21}
  744 = 4 × 186 = 4 × 6 × 31 = 24 × 31.
  24 = 4! = the order of S_4 (permutations of the n=4 tournament).

  So 744 = |S_4| × 31. And 31 is the 11th prime... hmm.

  More interestingly: 744 = 12 × 62 = T(5) × 62.
  And 62 = 2 × 31. Not obviously meaningful.

  THE HONEST ASSESSMENT: 744 doesn't have an obvious tournament
  interpretation. The moonshine numbers are controlled by the
  Monster group, not by tournament theory. The connection between
  the two is at the level of the j-FUNCTION, not the j-COEFFICIENTS.
""")

# ========================================================================
# PART 6: The vertex operator algebra parallel
# ========================================================================
print("=" * 72)
print("  PART 6: VOA STRUCTURE AND TOURNAMENT ALGEBRA")
print("=" * 72)
print(f"""
  A VERTEX OPERATOR ALGEBRA has:
  1. A graded vector space V = ⊕_n V_n
  2. A state-field correspondence Y(v, z) for each state v
  3. A vacuum vector |0⟩
  4. A Virasoro element ω

  The TOURNAMENT ALGEBRA (from S111) has:
  1. A graded space: tournaments at each order n
  2. A "nesting operator" N(T, bit): adds src+snk shell
  3. A "vacuum": the empty tournament (H = 1)
  4. A "Virasoro": the curvature g(x) = x³-x²-x-1

  THE PARALLEL:
  VOA vacuum ↔ transitive tournament (H = 1 = "nothing")
  VOA grading ↔ tournament order n
  VOA fields ↔ tournament polynomial evaluations at x
  VOA Virasoro ↔ curvature classifier g(x)

  The NESTING from S110 (H = c₀ - H_inner + c₂×sbs) is like
  the VOA STATE-FIELD CORRESPONDENCE: each inner tournament
  becomes a "field" that modifies the outer tournament's H.

  The PERIOD-2 PROPERTY (S111: two nestings restore the H sign)
  is like the VOA LOCALITY: fields at the same point commute
  (up to a sign from the Z₂ grading).

  This is NOT a proof that tournaments form a VOA.
  It IS a structural parallel that suggests tournaments
  have a VERTEX ALGEBRAIC structure at some level.
""")

# ========================================================================
# PART 7: The graded dimensions
# ========================================================================
print("=" * 72)
print("  PART 7: GRADED DIMENSIONS — TOURNAMENTS VS j")
print("=" * 72)

# Tournament iso class counts (A000568)
T_n = [1, 1, 1, 2, 4, 12, 56, 456, 6880]

# j-function coefficients (n=-1, 0, 1, 2, ...)
j_n = [1, 744, 196884, 21493760, 864299970]

# Monster rep dimensions (first few)
monster_reps = [1, 196883, 21296876, 842609326]

print(f"  n  | T(n) tournaments | j-coeff c_n | Monster dim d_n")
for i in range(min(len(T_n), len(j_n)+1)):
    t = T_n[i] if i < len(T_n) else "?"
    j = j_n[i-1] if 0 < i <= len(j_n) else ("1/q" if i==0 else "?")
    m = monster_reps[i-1] if 0 < i <= len(monster_reps) else "?"
    print(f"  {i:2d} | {str(t):>16s} | {str(j):>11s} | {str(m):>15s}")

print(f"""
  The three sequences have COMPLETELY DIFFERENT growth rates:
  T(n): roughly 2^(n²/2) / n!
  c_n: roughly exp(4π√n)
  d_n: controlled by Monster rep theory

  They are not the same sequence. But they MAY share structural
  properties (like the "+1" pattern) that come from their shared
  relationship to modular forms.
""")

# ========================================================================
# PART 8: The modular discriminant and the tournament discriminant
# ========================================================================
print("=" * 72)
print("  PART 8: DISCRIMINANTS — MODULAR AND TOURNAMENT")
print("=" * 72)
print(f"""
  MODULAR DISCRIMINANT: Δ(τ) = η(τ)²⁴
    Weight 12 cusp form. Vanishes at the cusp (q → 0).
    Delta = q * product (1-q^n)^24 for n=1,2,3,...

  TOURNAMENT DISCRIMINANT (kind-pasteur S18e): Δ_T = (H-1)/2
    Non-negative integer. Vanishes for transitive tournament (H=1).
    Δ_T = α₁ + 2α₂ + 4α₃ + ... (OCF expansion in base 2)

  THE PARALLEL:
  Modular Δ = cusp form (vanishes at cusp, nonzero elsewhere)
  Tournament Δ_T = cusp tournament analog (vanishes at transitive, positive elsewhere)

  BOTH measure "departure from the degenerate limit."
  The modular cusp = the transitive tournament.
  The modular Δ coefficient c(n) counts something structural.
  The tournament Δ_T counts cycle packings (in base 2).

  THE 24: η²⁴ means the modular discriminant involves 24th powers.
  24 = |S_4| = the permutation group at n=4 (last trivial tournament order).
  Also: 24 = 2³ × 3 = the order of the (2,3) Coxeter group.
  And: Klein quartic has 24 heptagons.

  MOONSHINE + TOURNAMENT SYNTHESIS:
  In moonshine: the Virasoro central charge is c = 24 (for V♮).
  In tournaments: the "central charge" is... the number of arcs?
  C(n,2) at n=7 is 21 = 3 × 7. Not 24.
  But C(5,2) = 10, and 24 = C(5,2) + C(5,2) + 4 = 24... reaching.

  THE HONEST CONCLUSION:
  The moonshine-tournament parallel is at the STRUCTURAL LEVEL:
  - Both have "+1" patterns (trivial representation / Rédei quantum)
  - Both are controlled by modular forms (j-function / curvature g(x))
  - Both have "cusp" limits (trivial representation / transitive tournament)
  - Both involve the number 12 (Monster weight / T(5) iso class count)

  But tournaments are NOT Monster moonshine. They are {'{2,3,7}'} moonshine —
  a much smaller group (PSL(2,7) = 168 instead of |M| ≈ 8×10^53).
  The tournament's modular structure is governed by the (2,3,7) triangle
  group, not the Monster.

  Perhaps: TOURNAMENT MOONSHINE is a "baby moonshine" for PSL(2,7)
  analogous to how Mathieu moonshine is a smaller version for M₂₄.
""")

print("\nDone. opus-2026-03-21-S123")
