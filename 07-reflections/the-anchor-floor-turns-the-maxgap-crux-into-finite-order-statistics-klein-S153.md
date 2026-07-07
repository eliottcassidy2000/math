---
source: klein-2026-07-07-S153 (HYP-4748)
status: sharpening + constructive certificate (not a proof). Extends kps-S57's density-floor
  reduction (inf_E E[maxgap] > 1/7): a FINITE ANCHOR SET turns the global max-gap order statistic
  into a computable lower bound; the AP is the minimizer with a COMFORTABLE 48% margin; and the
  origin saturates the AP's max gap (E[gap@0]=E[maxgap] for the AP, exactly).
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - max-gap
  - paley-zygmund
  - three-gap
  - anchor-floor
---

> **⛔ CORRECTION BANNER (klein-2026-07-07-S154; death-star-S1 HYP-4777, opus-S133; MISTAKE-118).**
> Two claims below are REFUTED: (1) "the AP is the minimizer of E[maxgap]" — FALSE, exact: the primitive
> saturated family `{2,4,…,24,13}` has `E[maxgap] = 145091/720720 ≈ 0.20131 < 93/440 ≈ 0.21136` (AP), and
> jump-move descent reaches `≈ 0.202005` at `{1,3,5,6,7,8,9,10,11,13,15,20,29}`; this session's descent
> used only ±1/±2 nudges and could not leave the AP basin. The margin over 1/7 is ≈ **41%**, not 48%, and
> the reduced open target is the DIRECT `inf_E E[maxgap] > 1/7`, not AP-minimality. (2) the "bonus lead"
> `inf_E E[gap@0] ≈ 0.156 > 1/7` — FALSE: direct minimization of `E[gap@0]` gives 0.137 (kps-S58) /
> 0.134 (opus-S133 `{6..58}`), both < 1/7. SURVIVING: the anchor floor (pointwise, correct), the AP
> origin-saturation identity `E[gap@0](AP) = E[maxgap](AP)`, the PZ-ceiling and mod-14-seam guardrails.
> NOTE ALSO (S154): μ_{1/7} IS still AP-minimized (death-star re-check) — the two functionals genuinely
> have different minimizers.

# The anchor floor: turning the max-gap crux into finite order statistics

Owner: work the real frontier, think Paley–Zygmund; pull in the fleet and extend. The fleet has
converged the Route-1 density floor onto one clean analytic target (**kps-S57, HYP-4747**, via a
reverse-Markov / Paley–Zygmund step): `μ_{1/7} = P_x(maxgap{frac(e·x)} > 1/7) ≥ (7/6)(E[maxgap]−1/7)`,
so the floor `μ_{1/7} > 0` reduces to

> **`inf_E E[maxgap] > 1/7`**  over cluster shapes `E` (`k = 8..13`),

with kps's partial `E[origingap] = 1/7` giving only the **non-strict** `E[maxgap] ≥ 1/7`. kps left the
strict margin ("max beats typical") open. This session gives it a **constructive, finite handle** and
confirms the margin is comfortable.

## The anchor floor (constructive lower bound)

For any finite anchor set `A ⊂ ℝ/ℤ`, the max gap dominates the gap at any single anchor, pointwise:

    maxgap(x) ≥ max_{a∈A} gap_a(x)   ⟹   E[maxgap] ≥ E[ max_{a∈A} gap_a ].

As `A` densifies, the RHS ↑ `E[maxgap]`. So `E[maxgap]` is the supremum of a **hierarchy of finite
order statistics**, each a computable integral — the global max gap becomes `max` over finitely many
**local** gaps. Numerically (exact-in-the-limit grid, `lrc14_maxgap_anchor_floor_klein_S153`), a tiny
anchor set already recovers the whole floor:

| k | `E[maxgap]` (AP) | `{0}` | `{0,½,¼,¾}` | `{j/8}` |
|---|---|---|---|---|
| 13 (AP) | 0.2114 | 0.2114 | 0.2114 | 0.2114 |
| 13 (spread) | 0.235–0.240 | 0.156–0.161 | 0.229–0.235 | =E[maxgap] |

`|A| = 8` matches `E[maxgap]` to grid precision for every family. This converts kps-S57's open
`inf_E E[maxgap] > 1/7` into the finite-anchor target `inf_E E[max_{a∈A} gap_a] > 1/7` — a `max`
over a fixed handful of local-gap integrals, each amenable to a moment computation.

## The AP is the minimizer, with a comfortable margin

Adversarial descent (`lrc14_maxgap_minimizer_klein_S153`, 60 random starts + local descent, `k=13`):
**nothing found below the AP**; the descent *converges to* `{1..13}` from random starts.

> **`inf_E E[maxgap] = E[maxgap](AP) = 0.2114 > 1/7 = 0.1429`, margin `+0.0685` (48% above threshold).**

So the density-floor crux is a **comfortable-margin extremal problem**, not a razor edge — the
route is robust. (Cross-confirms opus-S131's "AP uniquely minimizes `μ_{1/7}`" for the `E[maxgap]`
functional, extending the exhaustive `k≤10` evidence to a `k=13` descent.)

## The origin saturates the AP's max gap (clean structural fact)

For the AP, `E[gap@0] = E[maxgap]` **exactly** (0.21136 = 0.21136). Since `gap@0 ≤ maxgap` pointwise,
equality of means forces `gap@0(x) = maxgap(x)` for a.e. `x`:

> **For the AP orbit `{frac(i·x): i=1..13}`, the origin ALWAYS lies in the maximal gap.**

(This is the three-distance fact that the longest gaps of `{α,2α,…,Nα}` abut `0`.) Two consequences:
- It explains *why the AP's tight witness centers on the origin* (`t = 1/14`) and why the AP is the
  extremal shape — its single largest gap is "spent" straddling `0`.
- For the AP the single anchor `{0}` already achieves the full floor; the multi-anchor set is only
  needed to catch the *spread* families, whose max gap wanders off the origin.

## Guardrail: `E[gap@a]` is NOT structure-independent

Tempting (by analogy with the structure-free pairwise overlap `3/196`, opus-S131) to guess
`E[gap@a] = 2/(k+1) = 1/7`. **False.** Measured, `E[gap@a]` ranges `0.14–0.21`, depending on both the
anchor and the family (AP: `E[gap@0]=0.211` but `E[gap@generic]=0.162`). kps-S57's "origin-gap `1/7`"
holds only under its specific reduced normalization, not for the raw family orbit. Do not assume a
uniform origin-gap identity.

## Honest ledger

- **New / this session:** the **anchor-floor** lower bound `E[maxgap] ≥ E[max_{a∈A} gap_a]` (finite
  `A` recovers the floor; converts kps-S57's target to a finite order statistic); the **AP is the
  `E[maxgap]`-minimizer** with a **48% margin** (descent, `k=13`); the **origin saturates the AP's
  max gap** (`gap@0=maxgap` a.e.); the **`2/(k+1)` guardrail**.
- **Does NOT claim:** a proof of `inf_E E[maxgap] > 1/7`. That still needs either the AP-minimizer
  extremal (opus's `μ_{1/7}` program, exhaustive only to `k≤10`) or a proof that `inf_E E[max_{a∈A}
  gap_a] > 1/7` for a fixed finite `A` (the cleaner sub-target this session opens). Bonus lead:
  `inf_E E[gap@0] ≈ 0.156 > 1/7` in every sample — a *single*-gap sub-target worth a moment bound.
- **Extends / credits:** kps-S57 (the reverse-Markov reduction to `E[maxgap]`, HYP-4747), opus-S131
  (`μ_{1/7}` AP-minimizer, exact `477/1078`), opus-S132 (three-gap governing frame), mac-mini-S39
  (single-scale no-escape). Stays on the **sup** side (survives MISTAKE-117).

- **Files:** `04-computation/lrc14_maxgap_anchor_floor_klein_S153.py`,
  `lrc14_maxgap_minimizer_klein_S153.py` (+ `.out`); the earlier PZ-ceiling / relation-lattice probes
  `lrc14_paley_zygmund_good_set_klein_S153.py`, `lrc14_relation_lattice_seam_klein_S153.py`.
- **Next:** prove `inf_E E[max_{a∈{0,½}} gap_a} > 1/7` (a 2-anchor moment bound) or `inf_E E[gap@0] >
  1/7` — either closes kps-S57's target without the full AP-extremal.
