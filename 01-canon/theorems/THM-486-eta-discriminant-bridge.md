# THM-486: the η-discriminant bridge — the code discriminant is the pentagonal product, and the extremal-enumerator correction is secular (−42m), not two-rate

**Status:** PROVED (the −42m law, part B) + VERIFIED EXACT (parts A,C: independent
agent's η/θ computation, all asserts passed; my extremal-enumerator build,
n = 24m up to m = 154, exact integers). Scripts
`04-computation/extremal_enumerator_bridge_kps3_0611.py`,
`verify_eta_dictionary_agentr.py` (+ .out).
**Source:** kind-pasteur-2026-06-11-S3 (HYP-2419/2420). Bridges THM-485 (the
pentagonal Lyapunov theme) to THM-481/484 (the Gleason/Paley code line) and the
MOS extremal-enumerator theory. **CORRECTS HYP-2420's mechanism guess.**

## A. The code discriminant IS the modular discriminant (VERIFIED)

Under Construction A (W ↦ Θ via x = θ₃(q²), y = θ₂(q²)):
* W_{ê₈} = x⁸ + 14x⁴y⁴ + y⁸ ↦ E₄(q²) (verified to q¹⁰⁰).
* **The Gleason second-axis discriminant P₂₄ = x⁴y⁴(x⁴−y⁴)⁴ ↦ 16·η(q²)²⁴ = 16·Δ**
  — the CODE discriminant is exactly the modular discriminant, the 24th power of
  the pentagonal product Δ = q·Π(1−qⁿ)²⁴ (proof chain θ₂θ₃θ₄ = 2η³, θ₃⁴ = θ₂⁴+θ₄⁴,
  θ₄ = η²/η(q²) — each step series-verified; the (4,4,16) exponent pattern is
  exactly what lands on level 1 / SL₂(ℤ), vs (8,8,2) → level 2).
* W_{g₂₄} = W_{ê₈}³ − 42·P₂₄ (exact in ℤ[x,y]).

So the same pentagonal product that drives the partition function (THM-485,
1/Π(1−qⁿ)) is the second generator of the Type II weight-enumerator ring, in its
24th power. The Gleason ring is ℂ[E₄-image, Δ-image] = ℂ[W_{ê₈}, P₂₄], and P₂₄ is
η²⁴ in disguise.

## B. The extremal correction is SECULAR: c₁(m) = −42m (PROVED)

The extremal Type II enumerator at length n = 24m is
  W_n = Σ_{j=0}^{m} c_j · W_{ê₈}^{3(m−j)} · P₂₄^{j},  c₀ = 1,
the c_j fixed by the extremal conditions A₄ = A₈ = … = A_{4m} = 0 (minimum weight
4m+4). The basis is lower-triangular (P₂₄^j starts at y^{4j}), so the c_j solve
forward.

**Theorem.** The leading discriminant-correction coefficient is exactly
  **c₁(m) = −42m**  for every m ≥ 1.
Proof. The A₄ = 0 condition reads c₀·[A₄ of W_{ê₈}^{3m}] + c₁·[A₄ of P₂₄·W_{ê₈}^{3(m−1)}] = 0.
W_{ê₈}^{3m} = (1 + 14φ + φ²)^{3m} with φ = y⁴-monomial; its linear-in-φ
coefficient is 3m·14 = 42m. P₂₄ = x⁴y⁴(x⁴−y⁴)⁴ has lowest term x²⁰y⁴ (coefficient
+1), so [A₄ of P₂₄·(higher)] = +1. Hence 42m·1 + c₁ = 0, c₁ = −42m. ∎
(Verified m = 1..11: −42, −84, −126, …, −462 exactly; full c_j vectors at
n = 24,48,72,96,120 in the .out.)

**This is the MOS secular factor, made exact and elementary.** The leading
correction grows LINEARLY in m (rate 0 on the log scale, ratios c₁(m+1)/c₁(m) =
(m+1)/m → 1) — it is a polynomial/secular cocycle, not an exponential. The full
deeper coefficients c_j carry the saddle.

## C. CORRECTION to HYP-2420: the MOS crossover is same-rate secular, NOT two competing rates

HYP-2420 hypothesized the Mallows–Odlyzko–Sloane negativity (extremal Type II
enumerators acquire a negative coefficient at large length) as a *two competing
exponential rates* crossover. **This is WRONG** (literature audit, MOS 1975 §5 +
Jenkins–Rouse 2011 + Zhang 1999): the two boundary coefficients α_{μ+1}, α_{μ+2}
of the extremal form are Bürmann–Lagrange extractions of the SAME kernel
H(q)^r = Π(1−q^r)^{−24·s} — **same saddle, same exponential base** c₁ ≈ 69.1
(= F(y₀), the free-energy minimum of F(y) = e^{2πy}H(e^{−2πy})) — differing only
by a polynomially-growing SECULAR prefactor (24μ + 744, the j-function constant).
The sign flips when the linear secular factor 24μ crosses the constant amplitude
ratio (~1.64×10⁵), giving μ ≈ 6800 for lattices (exact dim 163264, Jenkins–Rouse)
and, for codes, the EXACT first Type II negativity **n = 3696** (Zhang; residue-0
branch; residue-16 survives to 3928). My build reproduces n = 3696 exactly.
The −42m of part B is the code-side avatar of this secular factor.

So the honest mechanism: ONE pentagonal-driven exponential rate (the η^{−24}
saddle), with positivity decided sub-exponentially by a secular cocycle. This is
the deterministic cousin of THM-485's sign-rigidity, where ONE rate (the
pentagonal product's boundary singularity) is modulated by the SIGNS to give
either subexponential (Euler) or exponential (everything else) growth.

## D. The Lyapunov family η^{−b} (honest framing)

The partition function (b = 1: 1/Π(1−qⁿ), subexponential e^{π√(2n/3)}) and the
MOS Bürmann kernel (b = 24: η^{−24}, exponential saddle c₁ ≈ 69.1) are the b = 1
and b = 24 members of one family of "η^{−b} growth constants." THM-485's all-plus
constant γ₊ = 0.548 (root of Σ(x^{g_k}+x^{ḡ_k}) = 1) is the b = 1, sign-flipped
analog; the random-sign γ_pent = 0.206 its stochastic version. The shared object
across the whole dispatch is **the pentagonal product Π(1−qⁿ)**: its reciprocal
powers η^{−b} furnish the growth engines, and the SIGNS decorating it (Euler's
alternation for partitions; the c_j secular alternation for extremal codes) decide
positivity. This is the precise sense in which "random-sign Lyapunov × pentagonal
× self-dual codes" is one circle of ideas, not three.

## Honesty

- c₁(m) = −42m is exact and proved; the deeper c_j are computed, not given closed
  form (the saddle lives there). The "η^{−24} saddle c₁ ≈ 69.1" is MOS's, quoted
  from the audited literature, not recomputed here to high precision.
- Part A's identities are the agent's exact q-series computation (asserts passed,
  brute-force lattice cross-checks); P₂₄ ↦ 16Δ and W_{g₂₄} = W_{ê₈}³ − 42P₂₄ are
  the load-bearing facts I use.
- The η^{−b} "family" (D) is a framing that is precise at b = 1, 24; intermediate
  b is not developed.

**Cross-refs:** THM-485 (the Lyapunov/rigidity side), THM-481/484 (Gleason/Paley),
THM-487 (the 72 frame, which uses part C's n = 3696 ≫ 72), HYP-2419 (PROVED here),
HYP-2420 (CORRECTED here), MOS 1975 / Zhang 1999 / Jenkins–Rouse 2011.
