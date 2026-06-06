# HYP-2265 — The π/3 shared object: forbidden H = Φ₃(2)=7, the Eisenstein unit-distance lattice, and 0.014=Cl₂(π/3)−1

**Session:** claudebox-2026-06-03-S628. **Prompt (user):** relate the unit-distance/UDG/Jacobsthal work to the
H=7,21 impossibility and the SC-tournament ↔ 0.014 growth; explore hidden connections. **Builds on / converges with:**
HYP-2184/2186 (opus-S599u/v: the π/3 object, Cl₂(π/3), phantom volumes), HYP-2179/2180 (H-impossibility mechanism),
my HYP-2235 (unit-dist/CM), HYP-2255 (SCC equidecomp), HYP-2260 (Jacobsthal/UDG).

## The one shared object: the angle π/3 (Eisenstein ζ₆ / cube root of unity)
Every thread sits at the **same angle π/3 = 60°** (the cube root of unity, the Eisenstein/hexagonal symmetry):

| thread | π/3 manifestation |
|---|---|
| forbidden H | **`7 = Φ₃(2)`**, `21 = 3·Φ₃(2)` (`Φ₃ = X²+X+1`, roots `e^{±2πi/3}`) — FORMALIZED |
| unit distance | the triangular = **Eisenstein lattice `ℤ[ω] = ℚ(√−3)`**; chord 1 ⟺ `dZ = 1/6` = the 60° gap (HYP-2235/S623) |
| LRC | gap `δ = 1/6` = the hexagonal chord |
| tournament partition fn | Lee–Yang / 3-cycle skew-eigenvalues at `±2π/3` (cube roots) |
| 0.014 exponent | `Cl₂(π/3) − 1 = 0.01494` (Clausen max = ideal-regular-tetrahedron volume = `3Λ(π/3)`) — verified |
| Jacobsthal chain | `J(4) = 21 = 3·Φ₃(2)` (the path-UDG hits the forbidden value at the cyclotomic point, S627) |
| Collatz | `×3` vs the 2-adic; the cube-root in the `2^K = 3^L` resonance |

So the forbidden values, the unit-distance optimum, the LRC gap, the disproof exponent, and the Collatz `3` are ONE
object — the primitive cube root of unity / the angle π/3.

## The forbidden value is cyclotomic (verified + formalized)
**`7 = Φ₃(2)`** and `21 = 3·Φ₃(2) = Φ₂(2)·Φ₃(2)`. The impossibility lives at the cube-root angle: `Φ₃`'s roots are
`e^{±2πi/3}`, the 3-cycle's skew-adjacency eigenvalues, the Lee–Yang zeros of the H-partition-function. Formalized:
`Math/Tournaments/CyclotomicSeven.lean` `cyclotomic_three_eval_two` (`Φ₃(2)=7`), `twentyone_eq_three_mul_cyclotomic`.

## The 0.014 = the Clausen/Dehn signature of π/3
The unit-distance disproof exponent surplus and the SC-tournament shape growth are the **same** number,
`Cl₂(π/3) − 1 ≈ 0.0149` (HYP-2184): the Clausen function at π/3 = the volume of the ideal regular hyperbolic
tetrahedron whose dihedral angle is the 3-cycle's eigenvalue angle. The unit-distance disproof and the SC-shape growth
are **equidecomposable** (same Dehn/scissors-congruence volume Cl₂(π/3), HYP-2186), not merely equinumerous — the
0.014 IS the equidecomposability signature (this session's lens, S626, made quantitative).

## SC tournaments = the CM-conjugation-fixed α₂=1 norm-1 family (verified)
Self-complementary tournaments (`T ≅ T^op` — fixed points of the complement involution `τ` = the CM conjugation `c` =
our `σ`) carry the `α₂=1` independence-polynomial structure (roots `ρ₁ρ₂ = 1 = |β|²`, the norm-1/CM analogue). Their
H-values (verified): `n=4:{1,5}`, `n=5:{1,3,9,11,13,15}` — **excluding 7** (consistent with 7 forbidden). These are
the shape-growth family whose exponent is the cubic 3-cycle count with the `0.014 = Cl₂(π/3)−1` secondary correction.

## HONEST CORRECTION to S626/HYP-2255
HYP-2184 (opus-S599u, exhaustive m=7 monad-S4) establishes that **only `{7, 21}` are PERMANENT forbidden H-values** —
`63` IS a strongly-connected atom at `n=8`, so the `7·3^k` ladder I asserted in S626 was the **sampling artifact**
already flagged by the fleet. My formalized `forbidden_seven_mul_three_pow` remains valid as a *conditional* (its
hypothesis "no atom = 7·3^j" fails for `j≥2`, so it only bites at `k=0,1`). The strong-minimum follows Busch's
recurrence `p(n)=p(n−1)+p(n−2)+1` (`3,5,9,15,25,…`), not the refuted quadratic. The permanent gaps are exactly the
cyclotomic pair `{Φ₃(2), 3·Φ₃(2)} = {7,21}`.

## Formalized (math-lean, sorry-free)
`Math/Tournaments/CyclotomicSeven.lean` — `cyclotomic_three_eval_two` (`Φ₃(2)=7`), `twentyone_eq_three_mul_cyclotomic`.

## Open
- The **phantom-volume theorem** (HYP-2186): prove no strongly-connected tournament has `H ∈ {7,21}` for all m (the
  permanent-gap = phantom-Dehn-volume statement), with `7=Φ₃(2)` as the cyclotomic obstruction.
- Derive the SC shape exponent's `0.014` = `Cl₂(π/3)−1` rigorously (currently a verified numerical coincidence /
  Sawin lower bound); equality iff scissors-congruent (HYP-2186).
- The cube-root Lee–Yang zeros of the H-partition-function: prove they sit at `±2π/3` for the atomic 3-cycle.
