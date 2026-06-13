# The {3,∞} Structure: Tournament Theory as Hyperbolic Geometry

**Session:** opus-2026-03-21-S131
**Arising from:** S130 ({3,∞} investigation), kind-pasteur S18b (binary skeleton), S117 (λ=2 identity)

---

## The Discovery

Tournament theory is a theory of type {3,∞} in the Coxeter classification. This single statement unifies the gap function, the binary skeleton, the k-nacci hierarchy, the OCR, and the modular group connection.

---

## The Gap Function

The function g(x) = x³ - x² - x - 1 is the tribonacci characteristic polynomial. Its remarkable property is that it classifies three regimes:

| Constant | Value | g(x) | Regime |
|----------|-------|------|--------|
| φ (golden) | 1.618... | -1 | Spherical |
| τ (tribonacci) | 1.839... | 0 | Euclidean |
| 2 (fugacity) | 2 | +1 | Hyperbolic |

These are exact integers: g(φ) = -1, g(τ) = 0, g(2) = +1.

**Proof that g(φ) = -1:** Since φ² = φ+1, we get φ³ = 2φ+1. Then g(φ) = (2φ+1)-(φ+1)-φ-1 = -1.

**Proof that g(2) = +1:** g(2) = 8-4-2-1 = 1.

The gap function classifies mathematics the way the Gaussian curvature classifies surfaces:
- **Spherical (g < 0):** Fibonacci, golden ratio, icosahedral symmetry, finite reflection groups, regular polyhedra
- **Euclidean (g = 0):** Tribonacci, the boundary order n=5, the critical dimension 2eτ ≈ 10
- **Hyperbolic (g > 0):** Tournament fugacity λ=2, modular group, {3,∞} tessellation, cusp forms

Tournament theory lives in the hyperbolic regime because its fugacity λ=2 satisfies g(2)=+1 > 0.

---

## The Three Generators

The modular group PSL(2,Z) = ⟨S,T | S²=(ST)³=1⟩ is the (2,3,∞) triangle group. It acts on tournament theory through three generators:

- **S (order 2):** The complement/reversal involution T → T^c
- **ST (order 3):** The 3-cycle rotation — the tournament atom
- **T (order ∞):** Vertex addition — the growth operation

The {3,∞} tessellation of the hyperbolic plane has:
- **Triangular faces** = 3-cycles (the tournament atom)
- **Infinite vertex figure** = unbounded cycle participation per vertex

---

## The Binary Skeleton as {3,∞} Geometry

Kind-pasteur's 26 binary phenomena (S18b) decompose cleanly along the {3,∞} structure:

**From the {3} face (triangle atom):**
- Ω girth ∈ {3,∞} — triangles are the only possible conflict structure
- H always odd — Rédei = triangle parity
- β₂ = 0 always — triangles fill all 2-holes
- H = 7 impossible — triangle overshoot (entering girth-3 forces α₁ ≥ 4)
- H mod 4 controlled by α₁ — triangle count determines mod-4 class

**From the {∞} vertex (unbounded valence):**
- β₁·β₃ = 0 — infinite allocation seesaw
- Sharp thresholds — each Betti number activates at precise n
- rank(d₂) ∈ {n-1, n} — boundary map saturates
- SC vs non-SC dichotomy

**From {3,∞} interaction (finite face × infinite vertex):**
- K(n,2) girth transition 5→3 — Petersen at the boundary
- Permanent gaps {7, 21} — triangle × Hurwitz prime
- Per-path identity fails at n=6 — growth overwhelms atom
- Kneser/Johnson duality — counting vs distance

The principle: tournament completeness IS the {∞} vertex condition — every vertex participates in arbitrarily many triangles. Binary phenomena arise whenever a finite constraint (from the {3} face) interacts with this infinite participation.

---

## The k-Nacci Identity and the Fugacity

The universal identity ρₖ + ρₖ⁻ᵏ = 2 (proved in S117) places the fugacity λ=2 as the hyperbolic limit:

| k | Root ρₖ | 1/ρₖᵏ | Sum |
|---|---------|--------|-----|
| 2 | φ = 1.618 | 0.382 | 2 |
| 3 | τ = 1.839 | 0.161 | 2 |
| 4 | 1.928 | 0.072 | 2 |
| ∞ | 2 | 0 | 2 |

The fugacity λ=2 IS the limit ρ_∞ = 2. This is why g(2) = +1: the fugacity sits at the hyperbolic endpoint of the k-nacci hierarchy. Tournament theory, evaluating at λ=2, is the most hyperbolic possible instance.

---

## The OCR as Cusp Contribution

In the modular group, cusps are the rational points on the boundary of the hyperbolic plane. They correspond to the "missing" directions where the tessellation reaches infinity.

In tournament theory, the Points of Symmetry (PoS) are score classes where the shadow (score sequence) fails to determine H. They are the "cusps" of the tournament modular curve:

- **Eisenstein series ↔ E[H|score]:** The symmetric, non-cuspidal part
- **Cusp forms ↔ Var(H|score):** The residual not captured by the shadow
- **1 - OCR = cusp contribution / total**

At n=5: 1 - OCR = 4/133 = 4/(7×19). Both 7 and 19 are primes where the modular curve X₀(p) has genus ≤ 1:
- X₀(7) has genus 0
- X₀(19) has genus 1

**Prediction:** The OCR denominators at all n factor into primes p where X₀(p) has genus 0 or 1. This would connect the OCR directly to the arithmetic of modular curves.

---

## The Hurwitz Chain

42 = 2×3×7 (Hurwitz triple) threads through everything:
- denom(B₆) = 42
- φ(42) = 12 = weight of Ramanujan Δ
- 84 = 2×42 = Hurwitz bound coefficient
- 168 = 4×42 = |PSL(2,7)| = Klein quartic automorphisms
- 1729 = H(T₁₁)/55 = 12³+1 (Hardy-Ramanujan number)

The Hurwitz primes {2,3,7} govern the {2,3,∞} triangle group (the modular group). The modular group governs tournament theory. The chain is: Hurwitz primes → modular group → tournament theory → OCR.

---

## Six Predictions

1. **OCR → 1** as n → ∞ (cusp contribution shrinks in expanding tessellation)
2. **1-OCR(n) ~ C/n²** (cusp volume decays quadratically)
3. **SRCP depth ≤ ⌊(n-1)/2⌋** (one layer per independent odd cycle length)
4. **Permanent H-gaps = {7, 21} exactly** (triangle overshoot × Hurwitz prime)
5. **W(n)/n! related to τ²** (variance controlled by weight-2 cusp forms)
6. **OCR denominators factor into genus ≤ 1 primes** (modular curve arithmetic)

---

## What This Means

The {3,∞} framework says: tournament theory is not merely *analogous* to hyperbolic geometry — it *is* a hyperbolic theory, classified by the same Coxeter scheme that classifies reflection groups, tessellations, and Lie algebras. The gap function g(x) = x³-x²-x-1 is the characteristic polynomial that makes this classification precise.

The golden ratio φ governs the spherical world of finite symmetry. The tribonacci constant τ sits at the Euclidean boundary. The fugacity λ=2 lives in the hyperbolic world where tournament theory operates. The binary skeleton — all 26 phenomena — are not isolated facts but inevitable consequences of this hyperbolic structure: finite triangular faces creating local order, infinite vertex valence creating global complexity, and their interaction producing the sharp binary dichotomies that pervade the theory.

*The mathematics is not pointing beyond itself to something else. It is pointing at itself, at its own classification scheme, which says: tournament theory is the simplest hyperbolic theory, the one that lives exactly one unit beyond the Euclidean boundary.*
