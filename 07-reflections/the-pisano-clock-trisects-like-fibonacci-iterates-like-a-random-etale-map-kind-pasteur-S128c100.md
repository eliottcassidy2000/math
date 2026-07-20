# The Pisano clock: the counterexample trisects like Fibonacci and iterates like a random étale map

**kind-pasteur-2026-07-20-S128c100** (HYP-8145) · owner: *"see how this could
possibly relate to the fact that the last digit of the fibbonacci numbers is
periodic every 60 values."*

## 1. The verdict up front

The relation is real, and it is exactly one level deep. The Fibonacci last-digit
period 60 and the Jacobian counterexample share their TRIPLING LAW — the same
polynomial identity governs both — but NOT their iteration arithmetic: the
counterexample's own "Pisano period" exists (new instrument, computed:
**π_F(10) = 24**) and obeys no Pisano-style law. The analogy lives in the fiber
direction and dies in the orbit direction, and both halves are informative.

## 2. The mechanism bridge (proved)

Binet + αβ = −1 give, verified symbolically:

    L_{3n} = L_n³ − 3(−1)ⁿ L_n          (Lucas index-tripling)
    even n:  L_{3n}/2 = T₃(L_n/2)        (EXACTLY the Chebyshev trisection)
    F_{3n} = 5F_n³ + 3(−1)ⁿ F_n          (the U-side companion)
    T₃ ∘ T₃ = T₉                          (the semigroup step)

THM-1335 says the counterexample's fibers solve 4T³ − 3T = m — the trisection
normal form. So **the identity that trisects the modulus angle of the Jacobian
counterexample is the identity that trisects Fibonacci/Lucas indices.** The
iteration ladder F → F∘F → … (degrees 3, 9, 27, klein-S327's tower) is the
Chebyshev index semigroup T₃ → T₉ → T₂₇, which on the Fibonacci side is
n → 3n → 9n. Depressed cubics everywhere: F's fiber cubics have no square term
(the trace law), and the tripling laws have no square term — same normal form,
same reason (trace zero across conjugates).

## 3. The 60-decode in the fleet's groups (exact facts, honestly graded)

π(10) = lcm(π(2), π(5)) = lcm(3, 20), and each factor lands in a group the
fleet already owns:

- **mod 2:** the Fibonacci shift M = [[1,1],[1,0]] has order 3 in
  GL₂(F₂) ≅ S₃ — a 3-CYCLE in the very group that is the counterexample's
  monodromy (forced by klein's Smith rule at degree 3). Mechanically: x²−x−1
  mod 2 = Φ₃, the golden ratio becomes a primitive cube root of unity.
- **mod 5** (the ramified prime of ℚ(√5)): char poly = (x−3)², order
  20 = ord(3)·5 = 4·5, and M⁵ = 3I — **projectively an order-5 pentagon
  rotation** in PGL₂(F₅), whose simple core PSL₂(F₅) ≅ A₅ has order 60
  (klein's "three sixties," T1534).
- So 60 = lcm(a 3-rotation living in S₃, a dressed 5-rotation living over A₅),
  and A₅ = ⟨3-rotation, 5-rotation⟩ IS the icosahedral group. **The Fibonacci
  last-digit clock is assembled from exactly the two rotation orders whose
  closure is the Smith-allowed quintic monodromy.**

Grading: the group facts are exact; the "sameness" of the S₃'s (Fibonacci's
GL₂(F₂) vs F's monodromy) is a shared-normal-form observation, not a shared
cause — the PROVED shared cause is §2's tripling identity. The shaped program
conjecture (wild, arbiter = the realization program): the degree-5 Keller rung
with A₅ monodromy — the icosahedral analog of the trisection counterexample —
would carry its modulus arithmetic on ℚ(√5), whose mod-10 Frobenius clock IS
Pisano-60. Fibonacci's 60 is the shadow, at the quintic rung, of what the
χ(−L) law is at the cubic rung.

## 4. π_F: the Pisano period of the counterexample (new instrument, computed)

Definition: π_F(n) = lcm of the cycle lengths of F on its eventual core in
(ℤ/n)³ (the core is the stabilized image; F restricts to a bijection there).
Fibonacci's π(n) is the linear, det = −1, dimension-2 case of the same
definition. Data (full table in the .out):

    π_F(2) = 1 (core 2 — the det ≡ 0 mod 2 degeneracy), π_F(5) = 24,
    **π_F(10) = 24** = lcm(π_F(2), π_F(5)) (CRT verified); π_F(7) = 24,
    π_F(11) = 900, π_F(13) = 792, π_F(17) = 46800, π_F(19) = 817236,
    π_F(23) = 60060, π_F(29) = 212520, π_F(31) = 260.

Findings, honestly reported:
- **No Pisano-style law.** Pisano has π(p) | p−1 (split) or | 2(p+1) (inert) —
  a torus/global-field law. π_F's cycle lengths include primes dividing no
  p^k ± 1 for small k (p = 19 has a 47-cycle; 47 ∤ p±1, p²±1, p²+p+1). Core
  sizes scale like ~√(p³) with coefficient ~2 — random-map birthday statistics,
  not toral structure.
- **The contrast IS the finding:** Fibonacci's clock is lawful because M is
  linear and unimodular (det −1 — the "trivial Keller map" of the plane); F's
  clock is lawless because F is genuinely nonlinear with det −2 (bad prime 2:
  π_F(2) = 1, the mod-2 collapse). The Chebyshev inheritance lives entirely in
  the FIBER/modulus direction (§2), not the orbit direction: **F trisects like
  Fibonacci and iterates like a random étale map.**
- Curios, flagged not claimed: π_F(10) = 24 (the η²⁴/Leech constant of the
  repo's pentagonal hub); π_F(5) = π_F(7) = 24 = 4!.

## 5. Prior threads connected

klein T1534 ("the three sixties": ord₁₀₀₁(2) = π(10) = |A₅| = 60), T1532/T1533
(Fibonacci skip-row sums, Zeckendorf exchange walk), the pentagonal/η²⁴ hub,
and the Smith d=5 row (D₅/F₂₀/A₅/S₅ allowed — dihedral-iff-odd includes D₅,
and |F₂₀| = 20 = π(5), same order different group, numerology-graded). The
realization program now has an arithmetic face: each realized monodromy rung
should carry a characteristic modulus field — ℚ(√−L(target)) fiberwise at the
S₃ rung, conjecturally ℚ(√5)-flavored at the A₅ rung — and the rung's mod-p
fiber statistics are that field's Frobenius clock. Fibonacci's 60 is the
oldest such clock in mathematics.
