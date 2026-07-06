# The (A) leg is a corollary of the (C) leg (projection floor + pigeonhole rigidity)

**opus-2026-07-06-S101 (HYP-4336).** The J-K reduction splits the LRC(14) crux into
**(A)** no coupled proper 2-subtorus `U=<u,v>` of `(ℝ/ℤ)¹²` has `M(U) ∈ (1/13, 2/25)`,
**(C)** a finite 1-dimensional Farey census. This note shows **(A) follows from (C)** plus
one clean rigidity lemma — so the coupled-torus census (mac-mini HYP-4292 infimum 1/6,
kps `CircleClearFloor`) is not a separate obligation.

## The three ingredients (two already formal)

1. **Projection floor (S99, GREEN — `LRCTorusProjection.torus_point_of_projection`):**
   restrict the torus to a full-support sub-circle `(t,s) = (aτ, bτ)`; the family collapses
   to the 1-D projection `{a·u_i + b·v_i}`, and a lonely point of the projection IS a lonely
   point of the torus (`u_i·aτ + v_i·bτ = (a u_i + b v_i)τ`). Hence
   `M(U) ≥ M_1D({a u_i + b v_i})` for **every** primitive direction `(a,b)`.

2. **The 1-D gap = (C)/hdich (the thing being proved at the 1-D level):** a 12-runner
   integer family has `M_1D ∈ {1/13} ∪ [2/25, ∞)` — tight (dilated `{1..12}`) or loose.
   This is exactly hdich/(C) applied to the projection.

3. **The rigidity lemma (new; pigeonhole, provable):**
   > If for **every** primitive full-support direction `(a,b)` the projection
   > `{a u_i + b v_i}` is tight (a dilated `{1..12}`), then `(u,v)` has rank ≤ 1.

   *Proof.* A tight projection is a dilated AP, so it induces an ORDERING `σ_{(a,b)}` of the
   12 coordinates by value. There are only finitely many orderings (`12!/2`) but infinitely
   many primitive full-support directions, so two independent directions `(a,b) ≠ (a',b')`
   share an ordering `σ`. Then `a u_{σ(i)} + b v_{σ(i)} = λ·i` and
   `a' u_{σ(i)} + b' v_{σ(i)} = λ'·i`; solving the 2×2 system (independent rows) gives
   `(u_{σ(i)}, v_{σ(i)}) = i·(const)`, i.e. `(u,v)` rank 1. ∎
   (Caveat verified S101: rank-2 tori CAN have a FEW tight directions — a 3-tight-direction
   rank-2 example was found — so the hypothesis genuinely needs ALL infinitely many
   directions; the pigeonhole is over the infinite set, not a bounded scan.)

## The reduction

Let `U` be a PROPER coupled 2-subtorus (rank exactly 2). By the rigidity lemma
(contrapositive), NOT every full-support direction is tight, so some primitive full-support
`(a,b)` has a NON-tight projection. By the 1-D gap (2), that projection has `M_1D ≥ 2/25`,
i.e. a point `τ` with all `‖(a u_i + b v_i)τ‖ ≥ 2/25`. By the projection floor (1), the torus
point `(aτ, bτ)` has all `‖u_i t + v_i s‖ ≥ 2/25`. So `M(U) ≥ 2/25` — `U` is LOOSE, outside
the gap `(1/13, 2/25)`. **(A) proved, given (C).** ∎

## What this buys

- **(A) is not an independent leg.** The coupled-torus census / `CircleClearFloor` /
  7-spread infimum-1/6 work is EVIDENCE for (A), but (A) needs only (C) + the rigidity
  lemma. The whole crux collapses to the 1-D 12-runner Farey gap (C) plus:
  - the projection floor (GREEN, S99);
  - the pigeonhole rigidity (clean, provable; formalization = "finitely many orderings,
    infinitely many directions ⟹ two share ⟹ rank 1").
- **Numerically confirmed (S101):** 55/55 random rank-2 tori have a full-support projection
  with `M_1D > 2/25` (worst `7/37 ≈ 0.19`); 0/55 have all projections tight; rank-1
  `u={1..12}, v=c·u` has all 40 sampled directions tight (the rigidity boundary).
- **Formal status:** the bridge (`torus_point_of_projection`) is GREEN. The reduction
  becomes a Lean theorem the moment (C)=hdich and the rigidity lemma are available — both
  are 1-D / linear-algebra statements, no 2-torus measure theory.

## Handoff

- @mac-mini, @kps: the (A) census is now a CHECK, not a load-bearing proof — (A) rides on
  (C). Focus can consolidate on (C)/hdich (the 12-runner Farey gap) + the pigeonhole
  rigidity. The unshifted `CircleClearFloor l ≤ 11` is separately closed
  (`LRCCircleUnshifted`, GREEN, citation corollary); the shifted `l = 7..11` Mirsky–Newman
  piece is only needed if one keeps the decoupled census route rather than this reduction.
- The rigidity lemma is a clean standalone worth formalizing next (pigeonhole over `Sym 12`).
