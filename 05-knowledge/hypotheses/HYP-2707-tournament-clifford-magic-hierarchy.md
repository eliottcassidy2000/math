---
id: HYP-2707
title: The tournament tractability hierarchy IS the Clifford-magic (Gottesman-Knill) hierarchy — c3 is a GF(2) quadratic form (Clifford), H is magic
status: OPEN; the degree facts VERIFIED (c3 deg 2, c5 deg 4); the Gauss-sum/stabilizer-rank tooling proposed
source: kind-pasteur-2026-06-20-S22
depends_on:
  - THM-554   # score partition function Z_n
  - THM-555   # score->OCF wall: c3 is the LAST score-determined OCF datum
related:
  - THM-027   # H-maximizer = regular/Paley
  - THM-064   # imaginary evaluation p0(2)=(-2)^m W(i/2)
---

# HYP-2707 - Tournaments are quantum circuits: the cut/cycle split is the Clifford/magic split

## The connection (new external lens: Gottesman-Knill / stabilizer formalism)

A tournament on n vertices is a point of `GF(2)^{C(n,2)}` (one bit per arc). The tiling model splits
this into the **cut space** (the n-1 base-path bits = a spanning tree) and the **cycle space** (the
`F=C(n-1,2)` tile bits) — the standard GF(2) Cut ⊕ Cycle decomposition of the edge space of K_n.

**Gottesman-Knill / stabilizer dictionary.** Over GF(2): *quadratic forms* (and their Gauss sums) are
exactly the Clifford-simulable, polynomial-time level; *higher-degree* phase functions are "magic"
(non-stabilizer), requiring the T (π/8) gate and exponential (stabilizer-rank) resources. CLAIM:

> **The tournament invariant hierarchy graded by polynomial degree in the tile bits IS the
> Clifford -> magic hierarchy. The score->OCF wall (THM-555) is the stabilizer->magic boundary.**

## The degree facts (VERIFIED, n=5,6,7)

- **scores** `s_v` = **degree 1** (affine) in the bits — each tile adds +1 to one endpoint (THM-554).
- **`c3 = C(n,3) - sum_v C(s_v,2)`** = **degree 2** = a GF/integer **QUADRATIC FORM** in the bits
  (max nonzero mixed finite-difference order = 2, verified n=5,6,7). Its distribution = the score
  partition function `Z_n` evaluated/Gauss-summed — the **Clifford-tractable** part. `c3 = alpha_1` is
  the leading OCF term.
- **`c5`** = **degree 4** (verified n=6); **`alpha_2`, higher `c_k`, and H** = higher degree = **magic**.
  `c_{2k+1}` has degree ~`2k`; H mixes all degrees.

So the EXACT reason "c3 is poly-time / score-determined but H is not" (THM-555, THM-442 "H not
cell-affine") is the EXACT reason "Clifford is classically simulable but T-augmented circuits are not":
**degree-2 (quadratic/Gauss-sum/stabilizer) is the tractable ceiling; degree-3+ is magic.** The
"magic degree" of a tournament invariant := its polynomial degree in the tile bits.

## PROVED CORE (kps-S22, exact n=4..8) — the c3-parity IS a Clifford Gauss sum

The c3-parity `q(b) = c3(b) mod 2` is a GF(2) quadratic (Boolean degree-2) function of the F tile bits.
Its associated bilinear form `B(x,y)=q(x+y)+q(x)+q(y)` has GF(2) **rank `r = 2*floor((n-1)/2)`** (= n-1 for
odd n, n-2 for even n): r = 2,4,4,6,6 for n=4..8. The standard Boolean Gauss sum then gives
`E[(-1)^c3] = (1/2^F) sum_b (-1)^{q(b)} = 2^{-r/2} = 2^{-floor((n-1)/2)}` — which is EXACTLY THM-555's
c3-parity formula. So THM-555's `E[(-1)^c3]` is the **Gottesman-Knill rank formula** for the Clifford
amplitude of the c3-quadratic-form, and the c3-parity is computable in poly(F) by GF(2) rank (no
2^F enumeration). The "symplectic rank" of the c3 magic-free core is n-1 (odd) / n-2 (even).

## Why this is more than analogy (precise hooks)

- `Z_n = (prod x_v) prod_{tiles}(x_a+x_b)` is literally a **partition function / Gauss sum**; setting
  `x_v = i^{(.)}` (cf THM-064's imaginary evaluation `p0(2)=(-2)^m W(i/2)`) turns the c3-quadratic-form
  sum into a **Clifford amplitude** (a GF(2) Gauss sum), which Gottesman-Knill computes in poly time.
- The H-maximizer = regular/Paley (THM-027) is the **most "magic"** tournament: maximal cycle-space
  (degree-high) content at fixed quadratic (score) content — the analogue of a maximal-magic state.

## Proposed tooling (the payoff for the repo's H-computation goal)

1. **Gauss-sum algorithm for the c3/score layer** (Clifford): compute the full c3-distribution as a
   GF(2) quadratic-form Gauss sum in poly(F) — sharper than the `Z_n` state-DP (which is exp in the
   score box). TEST: does the standard quadratic-form-rank Gauss-sum formula reproduce the c3-distribution?
2. **Stabilizer-rank decomposition of H:** write the magic (degree>=3) content of H as a sum of `chi`
   quadratic (Clifford) pieces; then `H` costs poly * 2^{chi}. The "stabilizer rank" `chi(T)` = the
   tournament's magic content; conjecture `chi` is small for near-transitive, maximal for Paley.
3. **Magic monotone:** is there a tournament invariant (a "mana"/stabilizer-Renyi analogue) that lower-
   bounds H-hardness and is itself score/c3-computable?

## Tests / next
- VERIFY the c3 Gauss-sum formula (quadratic-form rank => exact c3-distribution, poly time) vs `Z_n`.
- Compute the "magic degree" spectrum of `c_5, c_7, alpha_2` and the H-polynomial's degree profile.
- Connect to THM-064 (imaginary evaluation) as the literal Clifford amplitude.
- The road-coloring / Gibbs / Feynman-propagator / cat-map siblings: see the S22 connection-mining
  workflow and reflection `tournaments-are-quantum-circuits-clifford-cut-magic-cycle`.
