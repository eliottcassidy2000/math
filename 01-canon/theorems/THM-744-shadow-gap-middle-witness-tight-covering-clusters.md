---
id: THM-744
title: The shadow-gap middle-witness — for a speed set C with an even element, if max(C) < 6·(smallest even in C) then t=1/2+δ is lonely for {1}∪C for every δ∈(1/(14 e), 3/(7 max C)) (e = smallest even), an explicit non-empty good interval in the MIDDLE. This is the first UNCONDITIONAL middle-reaching witness for TIGHT covering clusters (ratio < 6) — exactly the residual that the isolated-far disc_v certificate and the simultaneous multi-peel (THM-735) do not reach; it attacks the AP-stability / shadow-gap rigidity (S296) at k=2, p=1/2
status: PROVED (klein-2026-07-13-S297; elementary, 6-line proof). Verified NG-dense on covering families ({1,90..101}, {1,30..40,42}, spread tight clusters) — the interval is all-good, zero violations. Handles every covering {1}∪C whose cluster C has max/min < 6 (the tight-cluster residual). Residual after it: covering C with ratio ∈ [6,13] and no isolated far element (a thinner slice; opus-S271 true-disc per-family, or a sharper multi-resonance witness).
source: klein-2026-07-13-S297
depends_on:
  - THM-523   # the non-covering witness t=1/q; THM-744 is the covering-side companion (a MIDDLE resonance witness)
related:
  - HYP-6590  # klein-S296 (the shadow-gap mechanism this makes quantitative)
  - HYP-6600  # klein-S297 (this theorem)
  - THM-735   # kps simultaneous multi-peel (bounded body + ≤6 far) — THM-744 handles the DISJOINT tight-cluster case
  - THM-731   # klein disc_v certificate (isolated far element) — THM-744 handles the NON-isolated tight case
  - THM-724   # covering-min M≥14/183>1/14, so for covering L>0 ⟺ the needed M>1/14
---

# THM-744 — the shadow-gap middle-witness for tight covering clusters

The quantitative form of the S296 shadow-gap mechanism, at the resonance `k=2, p=1/2`. It is the first
*unconditional* proof that a class of covering `{1}∪C` reaches the middle `[1/14,13/14]` (hence `L>0`),
attacking the AP-stability rigidity directly.

## Statement

Let `C⊂ℤ_{≥2}` be finite with at least one even element; put `e = min{c∈C : c even}` and `m = max C`.
**If `m < 6e`**, then for every `δ ∈ (1/(14e),\ 3/(7m))` (a non-empty interval) the time `t = 1/2 + δ`
satisfies `‖c t‖ ≥ 1/14` for every `c ∈ {1}∪C`. Hence `M({1}∪C) ≥ 1/14`, and the good set contains the
interval, so `L({1}∪C) ≥ 3/(7m) − 1/(14e) > 0`.

## Proof (elementary)

Fix `δ ∈ (1/(14e), 3/(7m))`, `t = 1/2 + δ`. For `cδ < 1/2` (which holds since `cδ < m·3/(7m) = 3/7`):

- **Even `c` (`c ≥ e`):** `c/2 ∈ ℤ`, so `‖ct‖ = ‖c/2 + cδ‖ = ‖cδ‖ = cδ` (as `cδ<1/2`). Since
  `δ > 1/(14e) ≥ 1/(14c)`, `cδ > 1/14`; and `cδ ≤ 3/7 < 1/2`. So `‖ct‖ = cδ ∈ (1/14, 3/7)`. ✓
- **Odd `c` (`c ≤ m`):** `c/2 ∈ ℤ + 1/2`, so `‖ct‖ = ‖1/2 + cδ‖ = 1/2 − cδ` (as `cδ<1/2`). Since
  `δ < 3/(7m) ≤ 3/(7c)`, `cδ < 3/7`, so `‖ct‖ = 1/2 − cδ > 1/2 − 3/7 = 1/14`. ✓
- **Speed `1`:** `t ∈ (1/2, 1/2+3/7) = (1/2, 13/14)`, so `‖t‖ = 1−t ∈ (1/14, 1/2) > 1/14`. ✓

The interval is non-empty iff `1/(14e) < 3/(7m)`, i.e. `m < 6e`. ∎

## Why it matters — the tight-cluster residual

`e` is the smallest even, so `m < 6e` is essentially "`max/min < 6`": **tight covering clusters**. These
are exactly the covering `{1}∪C` that the two class-level tools miss:
- the **isolated-far disc_v certificate** (THM-731/732, kps-THM-736) needs `c_max ≫` the rest — it fails
  for a packed cluster;
- the **simultaneous multi-peel** (THM-735) needs a body of `≥7` speeds and `≤6` far elements — for
  `{1}∪(12\text{-speed cluster})` the only small element is the outlier `1`, whose peel has all-energy
  disc, and there are `12>6` "far" speeds.

THM-744 handles them with an explicit rational-neighbourhood witness at `t=1/2`. E.g. `{1,90,…,101}`
(`e=90, m=101, m/e=1.12<6`): lonely on `t ∈ (1/2+1/1260,\ 1/2+3/707) = (0.5008, 0.5042)` — the S289
"decisive counterexample" is dispatched by a two-line inequality. Census: every sampled tight covering
`{1}∪C` (`min≥15`, moderate spread) satisfies `m<6e` and is certified.

## The AP-stability connection, and the honest residual

THM-744 is the covering-side companion of THM-523: THM-523 gives the middle witness `t=1/q` for
*non-covering* sets; THM-744 gives the middle witness `t=1/2+δ` for *covering tight clusters*, using the
`k=2` shadow gap. The factor `6 = (14−2)/2` is optimal for this method (other resonances `p=a/k` give
`(14−k)/k ≤ 6`; large-`k` resonances have vanishing margin). So the surviving residual is:
**covering `{1}∪C` with `C` of ratio `∈[6,13]` and no isolated far element** — a thinner slice than the
whole non-isolated stratum, handled per-family by opus-S271's true disc; its uniform class-level closure
(or a multi-resonance sharpening beyond factor 6) is the remaining rigidity.

## Addendum (klein-S298) — the PARITY-SPLIT refinement pushes the factor past 6 (up to 13 on the even branch)

The `m<6e` bound is not tight: at `t=1/2` the *odd* and *even* speeds constrain the gap differently, and
only the odds carry the factor `6`. Precisely, with `e=`min-even, `E=`max-even, `o=`max-odd (`o=0` if none):

> **Refined THM-744.** If `o < 6e` **and** `E < 13e`, then `t=1/2+δ` is lonely for `{1}∪C` for every
> `δ ∈ (1/(14e),\ \min(3/(7o),\ 13/(14E)))` (non-empty). Hence `M({1}∪C)≥1/14`, `L>0`.

*Proof.* Same as the main proof, with the even bound sharpened: an even `c` has `‖ct‖=cδ`, which stays
`≥1/14` for `δ>1/(14e)` **and** stays `≤13/14` (i.e. good) until `δ=13/(14c)`, tightest `δ<13/(14E)`; the
odd `c` gives `‖ct‖=1/2−cδ>1/14` for `δ<3/(7o)`. Non-empty iff `1/(14e)<\min(3/(7o),13/(14E))`. ∎

So the odds cap the factor at `6`, but the evens only at `13`. Consequences:
- **All-even `C`** (`o=0`, the odd constraint vacuous): factor **13** — `E<13e`, matching the THM-405
  bounded-ratio bound (an all-even cluster at `1/2` is a `C/2` first-window near `0`).
- **Odd-tight clusters** (large evens but `o<6e`): covered up to `E<13e`. Verified: 112 sampled covering
  clusters with `max/e` up to `8.7` that fail `m<6e` are certified by the refined bound (interval all-good,
  1315/1315).

**Multi-resonance (numerical, klein-S298).** Taking the widest gap over resonances `p=a/k` (`k=2,3,…`),
each with its own residue-split refinement, certifies far more of the ratio-`[6,13]` band than any single
`k`: `k=2,3` alone reach 96% at ratio 6–8, 82% at 8–11, 85% at 13–20 (vs 100% true `L>0`). The clean
single-factor statement is the `k=2` parity split above (`6` odd / `13` even); a uniform closed-form past
that is a *per-cluster* choice of resonance, not one constant — the residual `[6,13]`-ratio clusters with
an unfavorable parity/residue profile still need opus-S271's per-family true disc or the general
equidistribution. The optimal single-point factor is `6` (odds at `t=1/2` sit at the maximal distance
`1/2`); the parity split is what lifts the *even* half to `13`.

*(klein-S298, HYP-6610; verified `04-computation/lrc14_multiresonance_klein_S298.py`.)*

*Files: `04-computation/lrc14_shadowgap_thm_klein_S297.py`, `lrc14_multiresonance_klein_S298.py` (+.out).
HYP-6600, HYP-6610. Attacks the S296 shadow-gap rigidity (HYP-6590); companion to THM-523; fills the
tight-cluster gap between THM-731 (isolated far) and THM-735 (multi-peel).*
