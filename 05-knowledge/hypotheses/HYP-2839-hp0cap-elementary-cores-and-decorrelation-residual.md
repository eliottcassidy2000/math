---
id: HYP-2839
title: hp0cap (the wide cover bound p0(E)<=cap) -- elementary Lean cores formalized; the binding bound = the analytic decorrelation residual
status: CORES FORMALIZED (sorry-free) + p0 WIRED CONCRETE INTO SKELETON (LRC14.p0 := slowμ(coverSet), full skeleton builds sorry-free); binding k=8..12 = analytic residual (resonance bound p0<=p0_decorr), isolated
source: kind-pasteur-2026-06-22-S31
related:
  - OPEN-Q-108
  - HYP-2832   # the witness/p0 unification (floor reduces to hp0cap)
  - THM-534    # decorrelated closed form / moment dual
  - THM-546    # far-element decorrelation (comb) bound
---

# HYP-2839 -- hp0cap: elementary cores + the decorrelation residual

## What hp0cap is
`hp0cap = (p0(E) <= cap_k for every wide k-set E)` where `p0(E) = meas(coverSet E) =
meas{x : the phases {frac(e*x):e in E} hit all 6 inner sectors [j/7,(j+1)/7)}`. It is
the sector route's deep node AND (via HYP-2832, MISTAKE-084) the single analytic input
the witness floor needs (`p0 < cap => witnessG2 > 0`, positivity is all hpartA needs).

## Attacked it "like hpartA's reach core": ELEMENTARY CORES, sorry-free (LRCCoverBound.lean)
- `coverSet_mono` / `slowμ_coverSet_mono`: p0 is MONOTONE in the speed list (more speeds
  => easier to cover). So the wide bound's worst case is the full speed set.
- `six_le_card_of_coverSet_mem`: **six disjoint inner sectors need six distinct speeds**
  (an injective `Fin 6 -> E.toFinset` witness map; one phase lies in one half-open sector).
- `slowμ_coverSet_eq_zero_of_card_lt_six`: **|E|<6 => coverSet=empty => p0=0 <= cap.** The
  cover analogue of mac-mini's hnu1 max-gap pigeonhole. Closes hp0cap for k<6.
- `slowμ_coverSet_lt_cap_of_decorrelation` / `slowμ_coverSet_lt_cap`: the RESIDUAL ISOLATION
  -- hp0cap reduces (split by cluster size: k<6 trivial, k>=6 via decorrelation) to a SINGLE
  analytic hypothesis, the resonance bound `p0 <= p0_decorr`. (The hp0cap analogue of
  `witness_pos_from_wide_bound`.)
All axioms propext/Classical/Quot only; lake build EXIT=0.

## Codex S86g strict-cover handoff to the concrete witness floor
`LRCWitnessFloorConcrete.witness_pos_from_strict_cover_bound` now consumes the
strict output of `slowμ_coverSet_lt_cap` directly:

> `p0(E) < cap_k` and `cap_k <= meas(G_P)` imply
> `0 < meas((coverSet E)^c ∩ safeSet P)`.

This removes a small formal mismatch: the p0 positivity route does not need the
dual cap floor to be strict and does not need an explicit `delta` at the
concrete carrier layer.  The quantitative margin theorem remains useful for
finite-ruler error budgets, but the positive-PartA route can now use the strict
hp0cap residual plus non-strict `hmeasGP` exactly as stated.

## Codex S86g dense-complement proxy handoff
The same strict-cover output now feeds the formal max-gap proxy.  For anchored
`0 ∈ E`, `LRCDenseCovers.coverSet_compl_subset_denseSet_compl` proves
`(coverSet E)^c ⊆ (denseSet E)^c`, and
`LRCWitnessFloorConcrete.dense_compl_witness_pos_from_strict_cover_bound`
combines this with measure monotonicity:

> `p0(E) < cap_k`, `cap_k <= meas(G_P)`, and `0 ∈ E` imply
> `0 < meas((denseSet E)^c ∩ safeSet P)`.

This is deliberately weaker than claiming the full `goodSet` readout.  It
records the proved lower-carrier path from hp0cap to a `Dense17`-complement
event, leaving the sorted cyclic-gap equivalence with `goodSet`/`witnessG2` as
the next formal node.

The quantitative version is also formalized:
`dense_compl_witness_margin_from_wide_bound` carries the p0 margin itself,
`delta <= meas((denseSet E)^c ∩ safeSet P)`, from
`p0(E) <= cap_k - delta`, `cap_k <= meas(G_P)`, and `0 ∈ E`.

The phase-gap version narrows the remaining `goodSet` readout further:
`exists_phase_arc_empty_of_not_dense` proves the finite cyclic-gap step,
`(denseSet E)^c ⊆ phaseGapSet E`, and
`phaseGap_witness_margin_from_wide_bound` carries the same p0 margin to
`meas(phaseGapSet E ∩ safeSet P)`.  The open quotient is now the
speed-difference identity into `GoodSet.goodSet`, not the max-gap/empty-arc
finite geometry.

That quotient is now closed too: `LRCGoodSet.phaseGapSet_subset_goodSet` and
`denseSet_compl_subset_goodSet` prove the concrete `GOOD` readout, and
`LRCWitnessFloorConcrete.goodSet_witness_margin_from_wide_bound` carries the p0
margin to `meas(goodSet E ∩ safeSet P)`.  Thus hp0cap plus `hmeasGP` feeds the
actual slow-time `GOOD∩G_P` witness carrier; Part A still has to transfer that
slow-time carrier to the finite-ruler witness.

The shape-level readout is now bridged as well.  `LRCEventMeasureBridge` proves
`shape_goodSet_witness_margin_from_wide_bound` and
`shape_goodSet_witness_pos_from_strict_cover_bound`, so once a shape supplies
lists `Eof/Pof` and the readout equality
`witnessG2 = meas(goodSet(Eof) ∩ safeSet(Pof))`, the same p0/cap hypotheses give
`delta <= witnessG2` or `0 < witnessG2` directly.  This leaves the residual
hp0cap proof itself, the cap floor, and finite-ruler Part A; it is no longer a
slow-time event quotient problem.

The finite Part-A handoff now consumes that shape-level readout directly.
`LRCWitnessPartA.finite_witness_pos_from_goodSet_margin_shapes` and its uniform
arc-bound variant combine the readout equality, anchored `Eof`, p0 margin,
cap floor, and the `rhoK/arcCount/Vmax` error budget to produce positive
finite witness density.  The companion
`lrc14_from_finite_partA_goodSet_margin_shapes` packages this into the
conditional LRC14 assembly.  This does not prove hp0cap, hmeasGP, or the
finite-ruler approximation; it removes the extra abstract `delta<=witnessG2`
glue between those nodes and Part A.

After S31b made p0 a concrete named atom, the same handoff also has a
`DenseCovers.p0` surface.  `finite_witness_pos_from_goodSet_p0_margin_shapes`,
its uniform arc-bound variant, and
`lrc14_from_finite_partA_goodSet_p0_margin_shapes` accept
`DenseCovers.p0 (Eof s) <= cap s - delta s` and unfold p0 internally to the
cover-set measure.  This keeps hp0cap stated in the live skeleton vocabulary
while still feeding the actual measurable cover event into the goodSet/Part-A
route.

After the mac-mini S28/HYP-2840 pull, this handoff should be read as neutral
between hp0cap proof routes.  The older route proves the p0 margin through the
decorrelation residual; HYP-2840 proposes the cleaner THM-534 `L_y` route
(`p0 <= L_y` plus scalar extremality).  Either route supplies the same concrete
`DenseCovers.p0` margin expected by the Part-A wrappers.

## The residual (the genuinely analytic part, NOT elementary)
The binding cases k=8..12 route through the decorrelated closed form (THM-534, kps-S24):
> `p0(E) <= p0_decorr(E) = sum_t P_t^{(r)} p_t(B)`  [RESONANCE BOUND -- the residual]
> `p0_decorr(E) <= Q(k-1)`                          [finite combinatorial check, VERIFIED]
> `Q(k-1) < cap_k`                                  [exact rational margins 0.13-0.25, VERIFIED]
The middle + last are finite/rational (machine-checkable). The FIRST -- `p0 <= p0_decorr`,
i.e. the far phases decorrelate -- is the analytic core: THM-546's comb bound |Delta_w| <=
kappa V(E')/(pi^2 w) (single-far) + the joint Mordell-Tornheim/ET-Koksma tail (multi-far).
This is the genuine residual, exactly analogous to hpartA's Vmax-embedding + equidistribution.

## Why no elementary core bounds the binding cases
Marginal uniformity (`meas{frac(w x) in I} = |I|`, the decorrelation building block) is true
(mathlib `AddCircle.measurePreserving_zsmul`) but gives only a union bound (k/7 > 1, useless);
the binding bound REQUIRES the JOINT decorrelation (independence of far phases), which IS the
analytic content. So like hpartA, the elementary cores + residual isolation are the formalizable
progress; the binding bound stays a cited classical analytic input.

## Status / next
hp0cap cores DONE (sorry-free). The remaining residual = the resonance bound `p0 <= p0_decorr`
(THM-546 + Tornheim tail). To fully formalize: (i) marginal uniformity (w-fold cover or AddCircle
bridge), (ii) the decorrelated IE identity P_t^{(r)}, (iii) the comb/Tornheim tail -- a multi-session
analytic-infrastructure effort. -> OPEN-Q-108, HYP-2832, THM-534, THM-546.
