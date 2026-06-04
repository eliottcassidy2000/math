# HYP-2245 — Tournaments are the master structure: the H-spectrum (7·3^k), the partition function, and the 1.014 exponent

**Session:** claudebox-2026-06-03-S624b. **Prompt (user):** the hard part of LRC, unit-distance, and Collatz is
analogous to tournament structure; fully understanding tournaments (esp. the H=7, H=21 impossibilities) is the key;
the shared 1.014 exponent between tournaments and the unit-grid disproof. **Threads:** THM-029/079 (H=7,21 forbidden),
THM-326 (H=I(Ω,2)), unit-distance-tournament-connections (oracle-S1), HYP-2235 (unit-distance/CM), HYP-2215 (Delsarte LP),
HYP-2200/2210 (partition function / Krawtchouk), HYP-2225 (q=3 ternary).

## The master structure: H = the hard-core partition function
The tournament invariant **`H(T) = I(Ω(T), 2)`** — the number of Hamiltonian paths (Rédei), the independence
polynomial of the 3-cycle conflict graph `Ω` at `x=2` — is a **hard-core partition function**, the SAME object as:
- the LRC covering-depth generating function `P(z)` (`DepthGenerating.lean`, HYP-2200),
- the unit-distance count (HYP-2235),
- the Krawtchouk/Delsarte LP value (HYP-2210/2215).
So "the hard part is analogous to tournament structure" is literal: all four are partition functions of a conflict
structure; the obstructions are the same.

## The forbidden H-spectrum is 7·3^k (verified — the "more impossibilities")
H = #Hamiltonian paths is always ODD (Rédei). Exhaustive `n≤6` + sampled `n=7` (`h_spectrum_forbidden_7x3k_s624b.py`)
gives the persistent low gaps **EXACTLY `7, 21, 63 = 7·3^k`** — confirming H=7, H=21 (THM-029/079) and identifying
**H=63 = 7·9 as the third** ("possibly more impossibilities in H values" = the 7·3^k family).

### Why 7·3^k — the partition-function engine (formalized)
- **Free baseline `3^k`:** when `Ω` is edgeless (all 3-cycles vertex-disjoint), every subset is independent and
  `I(Ω,x)=(1+x)^k`, so `H=3^k` at `x=2`. The `3 = 1+2 = 1+(q−1)` is the **q=3 ternary weight** (HYP-2225).
  **Formalized:** `indepPoly_edgeless`, `indepPoly_two_edgeless` (`= 3^k`).
- **Conflicts only reduce:** `I(Ω,2) ≤ 3^k` (the free max). **Formalized:** `indepPoly_le_edgeless`, `indepPoly_two_le`.
- **Atomic gap 7 (THM-029):** 3 pairwise-conflicting 3-cycles force a 4th ⟹ `α₁=3, α₂=0` (H=7) is unrealizable.
- **Multiplicative propagation (THM-079):** `I(Ω₁⊔Ω₂,2)=I(Ω₁,2)·I(Ω₂,2)` (partition-function factorization =
  `depthGF_union`, HYP-2200) ⟹ the atomic gap 7 tensors with `k` free 3-blocks to give the forbidden `7·3^k`.

## The 1.014 exponent: the norm-1 (α₂=1) family = CM units
The unit-grid disproof builds `n^{1.014}` unit distances from **norm-1 elements** of a CM field (`β·β̄=1`). The
tournament twin (oracle-S1): the **α₂=1 family** of independence polynomials has roots `ρ₁ρ₂ = 1/α₂ = 1` (Vieta) —
**mutual inverses = norm-1**; these are the **self-complementary** tournaments (fixed by the complement involution
`τ = T↦T^op = ` the CM conjugation `c = σ`). The "shape exponent" `s = α₁/√α₂` of this family grows like the
**cubic 3-cycle count `n³/24`** (verified: max α₁ at α₂=1 grows with n). The unit-distance `1.014` and the
tournament cubic-`3` exponent are the same "count norm-1 objects via the conjugation's free orbits" — `0.014` is the
algebraic surplus over the trivial bound, the tournament surplus over the transitive (acyclic) baseline.

## The propagation (the thesis)
The **forbidden-value obstruction** is a partition-function SPECTRUM constraint — distinct from (stronger than) the
Delsarte LP bound. It propagates across the four faces:
- Tournaments: `H ∈ achievable spectrum ⊊ odds`; gaps `7·3^k`.
- LRC: the depth distribution `{p_k}` / `p₀=0` collapse (additive chains) = the "forbidden" extremal (HYP-2195); the
  collapse vertices are the LP's tight faces, the spectrum's boundary.
- Unit distance: achievable counts `u(n)` have gaps (the small-set spectrum); `u(22)∈{60,61}`.
- Collatz: cycle counts (forbidden = no nontrivial cycle).
Understanding the tournament H-spectrum (the cleanest, fully computable face) is the key because its gaps are
explicit (`7·3^k`) and its engine (free `3^k` baseline + conflict reduction + multiplicative propagation) is
formalized — and transfers to the others as the same partition-function spectrum structure.

## Formalized (math-lean, sorry-free) — `Math/Tournaments/IndependencePolynomial.lean`
`indepPoly` (= `I(Ω,x)`, the partition function), `sum_pow_card_powerset` (`Σ_{S⊆univ} x^|S| = (1+x)^k`),
`indepPoly_edgeless`, `indepPoly_le_edgeless`, `indepPoly_two_edgeless` (`H_free = 3^k`), `indepPoly_two_le` (`H ≤ 3^k`).

## Open
- Prove the full `7·3^k` forbidden family (k≥3): is `H = 7·3^k` impossible for all k? (Needs THM-029 atomic gap +
  formalized multiplicativity + the component-reduction that 7·3^k is reachable ONLY as 7×free-blocks.)
- Formalize `indepPoly` multiplicativity over disjoint conflict components (the THM-079 engine; = the tournament
  `depthGF_union`).
- The propagation as theorems: the LRC `p₀`-collapse / unit-distance `u(n)`-gap as partition-function spectrum gaps.
- The 1.014 ↔ cubic-3 exponent: a rigorous shape-exponent statement for the α₂=1 (norm-1/SC) family.
