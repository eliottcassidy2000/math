---
id: THM-1132
title: The r=6 "enumeration wall" is an artifact of a conservative horn constant — the SHARP component threshold 1/(7L) gives R_sharp ≈ 0.80 < 1, so the plain measure horn certifies the r=6 covering-killer stratum with no enumeration, no moment LP, and no weighted atlas
status: VERIFIED-on-broad-adversarial-search (max R_sharp = 0.8011, 792 cores × consecutive/large-scale/clustered/spread configs). The SHARP threshold 1/(7L) is PROVED sharp and sound (exhibits a real witness). The UNIFORM bound "max R_sharp < 1 over ALL configs" is OPEN — the same uniformity gap codex flagged for THM-1121, but now with a +25% margin (0.80 vs 1.0) instead of the conservative bound's −86% deficit (1.86 vs 1.0).
source: death-star-2026-07-18-S58
depends_on:
  - THM-1102   # r=6 max-T / KB=333 / 3.64e12 wall — built on the conservative min(N/(6mu), 1/(3L))
  - THM-1111   # MST prune (built to attack the same wall)
  - THM-1121   # codex weighted atlas (closes finite branch [92,332]); independently re-verified here
  - THM-1122   # kind-pasteur moment LP (70 survivors) — SUBSUMED by the atlas and dissolved here
related:
  - THM-1061   # where 1/(3L) was introduced ("a killer leaves gaps 6/(7k), so 1/(3L)=7k/18")
  - THM-724    # deep-well covering-min; the sharp horn certifies the deep well directly
script: 04-computation/r6_sharp_horn_search_deathstar_S58.py
output: 05-knowledge/results/r6_sharp_horn_search_deathstar_S58.out
---

# THM-1132 — the sharp measure-horn constant 1/(7L) dissolves the r=6 wall

## Setup (kind-pasteur's covering-killer stratum)

A covering 13-family in the r-ladder splits as **core P** (a size-`13−r` subset of `{1,…,12}`)
plus **r killers** `k_1<⋯<k_r` with `k_i > 13·max(P)`. For `r=6`: `P` is a 7-subset (792 of
them), killers `≥ 92`. To prove the family is `1/14`-lonely it suffices to exhibit **one**
real time `t` with `‖v t‖ ≥ 1/14` for all 13 speeds.

## The measure horn, sharpened

Remove the 5 smaller killers exactly:
`G = S(P) ∖ ⋃_{i<6} danger(k_i)`, where `S(P) = {t : ‖p t‖ ≥ 1/14 ∀ p∈P}` and
`danger(k) = {t : ‖k t‖ < 1/14}`. Let `L` = length of the largest component of `G`.

> **Sharp certification (PROVED).** If `L > 1/(7 k_6)` then the family is `1/14`-lonely.

*Proof.* `danger(k_6)` is a disjoint union of open arcs of width **exactly `1/(7k_6)`**
(half-width `1/(14k_6)` about each `j/k_6`), separated by `k_6`-safe gaps of width `6/(7k_6)>0`.
Let `J` be the largest component of `G`, `|J| = L`. `J ⊆ G` ⟹ `J` is safe for `P` and
`k_1,…,k_5`. If `J ⊆ danger(k_6)`, then `J` (connected) lies in a **single** danger arc, so
`|J| ≤ 1/(7k_6)`, contradicting `L > 1/(7k_6)`. Hence some `t*∈J` has `‖k_6 t*‖ ≥ 1/14`; at
`t*` all 13 speeds are safe. ∎

The threshold `1/(7L)` is **sharp** (a component equal to one danger arc has length
`1/(7k)` and contains no safe point). It exhibits a genuine witness, so it can **never
false-certify**.

## The point: the r=6 wall came from a factor-7/3 conservatism

THM-1102's horn used `T = min(N/(6μ), 1/(3L))`. The component bound `1/(3L)` is valid but
**not sharp** — it is `7/3 ≈ 2.33×` larger than the sharp `1/(7L)`. Define
`R = T / k_5` (`k_5` = largest removed). Certification succeeds when `R < 1` (then
`k_6 > k_5 > 1/(7L)`). Because the *same* `L` enters both, **`R_sharp = (3/7)·R_conservative`
config-by-config, exactly.** So THM-1102's adversarial max `R_conservative = 1.85794` becomes

  **max R_sharp = (3/7)·1.85794 = 0.796.**

The conservative constant put the threshold at `1/(3L)` *above* the killer range (forcing
`KB=333`, the `3.64×10¹²`-sextuple enumeration, the moment LP, the MST prune, and the
weighted atlas). The sharp constant puts it *below* even `k_5` — no enumeration needed.

## Verification (this session)

Broad adversarial search (`r6_sharp_horn_search_deathstar_S58.py`), max `R_sharp`:

| config family | max R_sharp |
|---|---|
| (a) consecutive small killers, **all 792 cores**, offsets 0–23, steps 1–5 | 0.8011 |
| (b) consecutive, moderate/large scale (300–3000) | 0.8011 |
| (c) **clustered large killers** (≥2 killers ≥333, tight) — *codex's flagged worry* | 0.8011 |
| (d) random spread, widths up to 1500 | 0.8011 |

The maximizer is **exact** (`r6_sharp_horn_deathstar_S58.py`, ℚ arithmetic): core
`[1,2,4,7,9,11,12]`, removed killers `[171,173,175,177,179]`,
`L = 72/72275`, sharp threshold `1/(7L) = 72275/504 = 143.40`, `k_5 = 179 > 143.40`, so
`R_sharp = 143.40/179 = 0.8011`. The conservative `1/(3L) = 334.6` is exactly what produced
`KB = 333`.

**Robustness** (`r6_sharp_horn_robust_deathstar_S58.py`): the sharp horn certifies the deep
well `{1..12,182}` (`L=1/168 > 1/(7·182)`), its tower `{1..12,364}`, and the worst `|core|=1`
body `{1..11,13,84}` — all genuinely lonely. It correctly does **not** fire on the small-`k`
tight families `{1..13}`, `{1..11,13,24}` (which are non-covering, outside the far-killer
stratum). No false certification anywhere.

## What this closes, and what remains

- **Reconciliation.** codex's THM-1121 atlas (independently re-verified here: 35 obligations,
  total weight exactly 505, all core-safe, max load 84 over [92,332] ⟹ 6·84=504<505) closes
  the finite branch. kind-pasteur's THM-1122 moment LP (70 survivors) is **subsumed** by it
  (weaker method on the smaller `q≤40` obligation universe vs the atlas's `q≤112`). The
  "70 survivors / run the r=6 enumeration / add S4" frontier chases an already-closed branch.
- **The whole r=6 stratum, in one stroke.** IF the uniform bound `max R_sharp < 1` holds
  over all configs, the sharp horn certifies **every** r=6 covering-killer family — finite
  branch *and* the unbounded tail — replacing the enumeration, the moment LP, the MST prune,
  the weighted atlas, and the measure-horn tail together. The r-ladder gives
  `R_sharp(r) = 0.22, 0.31, 0.42, 0.55, 0.80` for `r=2..6`; all `<1`. `r=7` (core ≤6 speeds)
  is the next boundary (`R_sharp ≈ 1.0`?) and is left open.
- **OPEN (the honest residual).** A *proof* that `max R_sharp < 1` uniformly — i.e. a uniform
  lower bound `L > 1/(7k_5)` over all 792 cores and all 5-killer configs. This is exactly the
  uniformity codex asked for ("prove a uniform measure-tail lemma"), but the target moved from
  an infeasible `R_conservative < 1` (max 1.86) to a `+25%`-margin `R_sharp < 1` (max 0.80).
  The sampled worst case was the consecutive step-two shape.  THM-1134 now
  closes that entire shape over all 792 cores and every legal scale; the
  uniform residual consists of the other five-killer shapes outside its new
  multiplier cone and separated-ratio gate.

## Toward the uniform bound: the G(σ) reduction (death-star-S58, proof progress)

The uniform bound `max R_sharp < 1` is, for the worst shape, reduced to an exact
**one-variable** condition. Consecutive killers `{b, b+2, …, b+8}` have, at time `t`,
phases `(b+2m)t = bt + m·σ` with `σ = 2t` — an **arithmetic progression of step σ**. So
"all 5 killers safe at `t`" is "`φ = bt` avoids 5 danger-arcs of width `1/7` centered at
`{0, σ, 2σ, 3σ, 4σ}`". As `t` ranges over a core-safe arc, `φ = bt` sweeps `b·(width) ≫ 1`
full turns, realizing the **largest allowed-φ gap**

  **`G(σ) = ` largest gap in `[0,1) ∖ ⋃_{m=0}^{4}(mσ−1/14, mσ+1/14)`,**

and the resulting `t`-gap is `≈ G(σ)/b`. Hence

  **`R_sharp ≈ 1/(7·G(2t))`, so `R_sharp < 1 ⟺ ∃ core-safe t with G(2t) > 1/7`** — an
  *existence/MAX* condition, not an average.

**Lemma (PROVED, exact rational — `r6_Gsigma_exact_bands`).** `G(σ)` is piecewise linear with
breakpoints at `σ=(k±1/7)/d`, `d≤4`; it satisfies `G(σ) > 1/7` exactly on the bands
`(0,1/7) ∪ (2/7,1/3) ∪ (3/7,4/7) ∪ (2/3,5/7) ∪ (6/7,1)` (and reflections), with
`min_σ G = 2/35` at `σ=1/5` (arcs evenly spread — the *bad* alignment), `G(1/3)=4/21`,
`G(1/2)=5/14`.  Their exact total measure is `11/21` (the earlier `~64%`
mental sum was incorrect).

**Worst-shape verification.** The maximizer sits at `σ ≈ 0.31 ∈ (2/7,1/3)`: exact
`L·b = 12312/72275 = 0.1703 > 1/7`, so `R_sharp = 0.8011`. A wide exact/float scan
(`r6_worst_shape_finite_check`) confirms `R_sharp < 1` for **all** `b∈[157,4000]` and steps
1–4 on core `[1,2,4,7,9,11,12]`; the maximizer `b=171` is *stable across scale* (`R_sharp`
does not grow), strong evidence `0.8011` is the true global max for this shape.

**What is proved vs. what remains.**
- PROVED: the phase-AP reduction; the exact `G(σ)>1/7` band lemma; the exact worst-config
  value; `R_sharp<1` on `b∈[157,4000]` (finite check).
- **Independent worst-core proof (death-star-S58).**  For the consecutive
  step-two shape on core `{1,2,4,7,9,11,12}`, the two former gaps are closed
  search-free: an exact finite head handles `157≤b≤399`, while the rational
  witness `t*∈A*=(71/154,13/28)` handles `b≥400` with margin `9/(77b)`.
  Here `2A*⊂(6/7,1)`, so the witness replaces the anticipated `G`-drift
  error estimate.
- **All-core strengthening (THM-1134).**  An exact 792-core rectangle atlas
  closes `b≥164`, and an independent 12,771-row endpoint bank closes every
  legal `b≤164`.  Thus no asymptotic `G`-drift estimate remains anywhere in
  the step-two family, not merely on the old worst core.
- **Remaining for uniform r6.**  Other five-killer shapes outside THM-1134's
  cone `B≥17 max(A,80)` and outside its exact separated-ratio `Q5` gate.

Historically, `G(σ)` exposed the right one-variable geometry.  THM-1134's
fixed rectangles are the exact, drift-free realization of it on the
step-two family; the multiplier chart is the more flexible object for general
shapes.

## Methodological note (for the fleet)

Before building enumeration machinery against `R<1`, check whether the horn constant is sharp.
Here a `7/3` over-count in `1/(3L)` manufactured a `3.64×10¹²`-sextuple "wall" and three
successive fixes (finite horn → moment LP → MST prune → weighted atlas) to scale it; the sharp
`1/(7L)` walks over it with room to spare. Same genus as the MISTAKES.md "value in a bound
treated as sharp" traps, but inverted: a *loose* bound treated as *binding*.
