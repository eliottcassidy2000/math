---
id: HYP-3761
title: THREE results toward tight-set completeness for n=14 (the {AP,GW} census), the multi-patch moduli, and the off-class double-heal. (1) NO OFF-CLASS DOUBLE-HEAL FLOOR-ACHIEVER: skip v, add 2v reaches the LRC floor 1/n ONLY at n=2 mod6 (always v=n-2), verified n=5..30; off-class the best double-heal is M=2/(2n-1)>1/n -- so klein-S49's n=2 mod6 criterion is SHARP even allowing any skip position, not just v=n-2. (2) MULTI-PATCH MODULI: every tight set is a residue-multi-patch of the AP (residues mod n = {1..n-1} complete, or one-gap {1..n-1}\{k}, by the HYP-+2914 necessary conditions); the j-patch moduli generate the near-AP family -- j=0 AP, j=1 GW (UNIQUE single-swap, verified adds<=100), j=2 cross-types. (3) COMPLETENESS REDUCTION: the AP is the UNIQUE RESIDUE-COMPLETE tight set (no lifted residue-complete set is tight, verified <=2 lifts m<=3) -- so every NON-AP tight set has the one-gap structure; in the one-gap case the duplicating lift is small (drop 12 => uniquely GW at lift m=1). Reduces {AP,GW} completeness to bounding the lifts in the one-gap multi-patch case (the residual scale-separation gap).
status: (1) VERIFIED n=5..30 (clean negative). (2) structural + verified: j=1 uniquely GW (adds<=100); j=2 EMPTY for n=14 (0 tight of 36270 sets, adds<=44). (3) sub-results VERIFIED (AP unique residue-complete <=2 lifts; one-gap lift bounded to m=1). So {AP,GW} is COMPLETE for n=14 within all tested bounds. Does NOT fully PROVE completeness -- reduces it to bounding lifts/patch-count in the one-gap case (open, the scale-separation residual of HYP-+2909). Builds on klein-S48/S49, opus-S1, HYP-+2913/+2914.
source: mac-mini-2026-06-30-S63
depends_on:
  - HYP-+2914  # the tight-set necessary conditions (residues = complete or one-gap) -- the reduction rests on this
  - HYP-3753   # klein-S49 patch-tuning unification (the double-heal / n=2 mod6 criterion this sharpens)
related:
  - HYP-3750   # klein-S48 tight classification (GW-type + cross-type)
  - HYP-+2913  # three-gap/Steinhaus tight-locus + census (the completeness route)
  - HYP-+2909  # tightness forces apex-7 binding pair + denom 14; the scale-separation residual
  - HYP-3749   # difference-closed = AP (the residue-complete analog)
  - HYP-3760   # my S62 confirmation/extension of the classification
results:
  - 04-computation/tight_multipatch_and_completeness_macmini_20260630.py
  - 05-knowledge/results/tight_multipatch_completeness_macmini_20260630.out
  - 05-knowledge/results/tight_double_heal_offclass_macmini_20260630.out
  - 05-knowledge/results/tight_residue_complete_liftbound_macmini_20260630.out
---

# HYP-3761 -- tight completeness, the multi-patch moduli, and the off-class double-heal

Three sub-questions on the tight-set census `{AP, GW}` for `n=14` (M(S)=1/n; the difference-closed tight sets
are the dilated APs, HYP-3749; GW = the skip-(n-2)/double-patch, klein-S49 HYP-3753).

## (3) The off-class double-heal -- a clean negative
**Question:** does the double-heal (skip `v`, add `2v`) ever reach the floor `1/n` when `n ≢ 2 (mod 6)`?
**Answer: NO** (verified `n = 5..30`). `{1..n-1}\{v} ∪ {2v}` has `M = 1/n` **only** at `n = 8,14,20,26`
(`n ≡ 2 mod 6`), always at `v = n-2`. Off-class, the best double-heal (over all `v`) plateaus at
`M = 2/(2n-1) > 1/n` and never touches the floor. So klein-S49's `n ≡ 2 (mod 6)` criterion is **sharp not only
for `v=n-2` but for every skip position** -- no other doubling rescues the floor off the residue class. (The
double heals the resonance hole `1/v` down to `2/(2n-1)`, but only the Jacobsthal-gated `v=n-2, n≡2 mod 6`
takes it all the way to `1/n`.)

## (2) The multi-patch moduli
The tight-set necessary conditions (HYP-+2914) force the residues mod `n` to be **complete `{1..n-1}`** or
**one-gap `{1..n-1}\{k}`**. Hence **every tight set is a residue-multi-patch of the AP**: pick residues
(complete or one-gap) and lift some to `r + n·m`. The `j`-patch moduli (skip `j` residues, add `j` lifted
elements) generate the near-AP family:
- **`j=0`**: the AP `{1..n-1}` (residue-complete, no lift).
- **`j=1`**: for `n=14`, the **unique** single-swap tight set is GW `{1..11,13,24}` (drop 12, duplicate residue
  10 via `24=10+14`) -- verified with additions up to **100**, far beyond `2n`.
- **`j=2`**: cross-types exist at `n=8` (`{1,4,5,6,7,11,13}` = drop residue 2, duplicate residue 5) but
  **NOT at `n=14`**: the `j=2` enumeration (adds `<= 44`, **36270 sets tested**) finds **0 tight**. So `n=14`
  -- unlike `n=8` -- has no cross-type 2-patch tight set within the bound.

So the moduli space is `(dropped residues) x (duplicated residues) x (lift vector)`, and the tight locus is the
(finite, per bound) set of lift vectors on which `M` collapses to `1/n`.

## (1) Completeness reduction
Two verified sub-results narrow `{AP, GW}` completeness:
- **The AP is the UNIQUE residue-complete tight set.** No residue-complete lift (all 13 residues present, one
  speed each, some lifted to `r+14m`) is tight except the all-zero-lift AP -- verified over all `<= 2` lifted
  residues with `m <= 3`. So **every non-AP tight set has the one-gap residue structure** (a genuinely dropped
  residue), matching GW.
- **The one-gap lift is bounded.** For the one-gap structure `drop 12, duplicate residue d via d+14m`, the only
  tight instance is `d=10, m=1` = GW; no other `d`, and no larger `m`, is tight. The duplicating lift is small.

**Net:** within all tested bounds, `{AP, GW}` is COMPLETE for `n=14`: residue-complete `=>` uniquely AP
(`<=2` lifts, `m<=3`); `j=1` single-swap `=>` uniquely GW (adds `<=100`); `j=2` `=>` **none** (36270 sets,
adds `<=44`). So `{AP, GW}` completeness reduces to **bounding the lifts / patch-count in the one-gap case** --
showing no one-gap tight set uses a large lift or `>=2` patches. The empirical evidence is strong (no escapee
up to adds 44-100); the residual is the large-lift / high-multiplicity one-gap tail, exactly the
scale-separation gap flagged in HYP-+2909.

## Honest scope
(1) is a complete verified negative (`n<=30`). (2) is the structural moduli description + the `j=1` uniqueness
(adds `<=100`) and the `j=2` enumeration (bounded). (3) reduces completeness to the one-gap lift bound via two
verified sub-results (AP unique residue-complete; one-gap lift small), but does NOT close the lift bound -- that
remains the open scale-separation residual. No claim of a full completeness proof; a genuine reduction + the
sharp off-class answer.
