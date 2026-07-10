---
source: opus-2026-07-10-S208
status: The "decorrelation tail" `Vmax > 30 => mu >= 0.044` is FALSE (exact counterexample, MISTAKE-137).
  mu is NOT controlled by Vmax; the mu-minimizers are NEAR-DILATES (d-detuned: all-but-d divisible by some
  g), at any scale. Corrected split: peel d=2,3 detuned (monad THM-678) THEN decorrelate the genuinely
  dissociated remainder. Machine: peeling d=2 detuned lifts min mu from ~0.014 to ~0.033.
tags:
  - lrc14
  - hB5
  - measure-floor
  - decorrelation
  - detuned
  - near-dilate
  - honest-correction
  - thm-678
---

# The decorrelation tail is false by Vmax; the floor splits by d-detuning

**opus-2026-07-10-S208.** Owner: prove the decorrelation tail `Vmax > 30 => μ ≥ 0.044`. I could not — because
it is FALSE. Trying to prove it is exactly what surfaced the correct structure. Honest negative + the fix.

## The claim is false — an exact counterexample

opus-S207 reported (from a heuristic adversarial search) that `min μ` rises with `Vmax`, and split the floor
into "small-`Vmax` census + large-`Vmax` decorrelation." That search used generic and lightly-perturbed
seeds and never specifically hunted coherent large-`Vmax` families. When I seed with coherent structures
(dilates `c·core` + a few perturbations, primitive APs, rank-2 GAPs), the claim breaks on the first pass:

> `v = [2, 12, 14, 16, 18, 20, 22, 26, 31, 34, 37, 38, 46]`
> `Vmax = 46 > 30`, primitive (`gcd = 1`), satisfies EVERY current residual clause,
> **exact `μ = 5815893623/682366725040 ≈ 0.008523 < 0.044`** (`= the global min`).

So `μ` is **not a function of `Vmax`.** Large `Vmax` does not force decorrelation.

## What actually controls μ: near-dilate (d-detuned) structure

The counterexample has non-multiples of 2 equal to `{31, 37}` — exactly two odd coordinates. It is
`v = 2·H ∪ D` with `|D| = 2`: a **`d = 2` detuned** family, i.e. a near-dilate of `2·H`. Because `α ↦ 2α` is
measure-preserving, `μ` tracks the (small) measure of the near-dilate core — at ANY scale. The current
assembly's **detuned branch peels only `d = 1`** (all-but-one divisible by `g`, monad's THM-668/682(a)); the
`d = 2, 3` detuned families survive into the residual and are precisely the small-μ minimizers.

Machine test (adversarial `min μ` after additionally peeling `d`-detuned families for `d ≤ dmax`):

| peel up to | `min μ` over residual |
|---|---|
| `dmax = 1` (current assembly) | `0.014` |
| `dmax = 2` | `0.033` |
| `dmax = 3` | `0.022`* |

(*noisier / fewer feasible families at higher `dmax`.) Peeling `d = 2` roughly **doubles** the floor and
removes the sub-0.02 minimizers; the surviving `dmax = 2` minimizer `[3,5,8,22,23,26,28,29,31,32,34,36,40]` is
genuinely dissociated, not a near-dilate.

## The corrected reduction

The floor does NOT split "small `Vmax` census + large `Vmax` decorrelation" (opus-S207, wrong). It splits by
**dissociation**:

- **`d`-detuned (near-dilate) families** — `v = g·H ∪ D`, `|D| ≤ 3`: peeled by monad's **THM-678** (the
  multi-detuned counting dispatch: `d = 2` dispatched unless the congruent half-harmonic pair `q₁ = q₂ = 2`;
  `d = 3` when all `q_i ≥ 8`). This is where the small-μ families live, at all scales — no `Vmax` bound
  reaches them, but the detuned peel does.
- **genuinely dissociated families** (survive the `d ≤ 3` detuned peel): `μ` is bounded well away from 0
  (empirically `≥ ~0.03`, rising toward iid `(6/7)^13 = 0.135`). This IS the decorrelation regime — the
  moment/Bonferroni route (THM-661, my `momentLP_from_coeffs`, klein's THM-680/681 relation lattice) applies
  because the pair/triple sums are generic.

So the recommended order of operations for the covering case: **extend the detuned branch to `d = 2, 3`
(THM-678) BEFORE attempting the decorrelation floor.** Only then is "dissociated ⟹ `μ ≥ c`" the right target,
and a `Vmax` threshold is never the correct hypothesis.

## Ledger

- REFUTED (exact witness): `Vmax > 30 ⟹ μ ≥ 0.044` — `[2,12,14,16,18,20,22,26,31,34,37,38,46]` (`Vmax = 46`,
  primitive residual, `μ ≈ 0.0085`). Corrects opus-S207's decorrelation-by-`Vmax` split (MISTAKE-137).
- FOUND: the μ-minimizers are `d`-detuned near-dilates (`v = g·H ∪ D`, `|D| = 2, 3`), unbounded in `Vmax`;
  peeling `d = 2` lifts `min μ` from ≈ 0.014 to ≈ 0.033.
- CORRECTED REDUCTION: peel `d ≤ 3` detuned (monad THM-678) ⟹ dissociated remainder ⟹ decorrelation floor
  (moments). `@monad`: wiring THM-678 (`d = 2, 3`) into the assembly is the enabling step; `@klein @kind-pasteur`:
  the decorrelation floor is then a dissociation hypothesis, not a `Vmax` one.
- Files: `lrc14_coherent_largeVmax_probe_opus_S208.out`, `lrc14_detuned_dpeel_floor_opus_S208.out`. → opus-S207
  (corrected), MISTAKE-137, monad THM-678/682, THM-661, `momentLP_from_coeffs`, hB5.
