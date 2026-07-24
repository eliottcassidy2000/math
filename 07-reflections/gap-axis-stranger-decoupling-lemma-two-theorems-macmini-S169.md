---
source: mac-mini-2026-07-23-S169 (Opus 4.8)
status: TWO THEOREMS, PROVED (fully exact rational arithmetic, no floating point in the decision path),
  both via one new reusable lemma. (1) SGC'(13) holds on the single-perturbation family: no
  ({1..13}\{j}) u {w} has gap in (1/14, 3/41); extremal 3/41 at {1..11,13,36}. (2) TIGHT CLASSIFICATION on
  that family: the only primitive tight sets are the AP {1..13} and the Goddyn-Wong sporadic
  {1..11,13,24} -- a rigorous exhaustive instance of OPEN-Q-108's {AP,GW} conjecture. Both follow from an
  explicit GAP-AXIS STRANGER-DECOUPLING LEMMA that makes the stranger search FINITE with a computable bound.
tags: [lrc, lrc14, spectral-gap, sgc, stranger-decoupling, tight-locus, open-q-108, theorem, exact]
related: [THM-518, THM-522, THM-523, THM-541, kps-S133, opus-S4, mac-mini-S169-sgc-refutation]
---

# A gap-axis stranger-decoupling lemma, and two theorems it proves

**mac-mini-2026-07-23-S169.** Continuation of my SGC(13) refutation (band = Ostrowski rungs `k/(13k+1)`,
corrected buffer `1/574`). Here I make the corrected conjecture a **theorem on the single-perturbation
family**, and get a **tight-locus classification** on the same family, both from one lemma.

Notation: `gap(S) = max_τ min_{v∈S} ‖vτ‖`; `f_C(τ) = min_{v∈C} ‖vτ‖`; tight ⟺ `gap = 1/14`.

## The lemma (new, explicit, and the engine for both results)

> **LEMMA (gap-axis stranger decoupling).** Let `C` be a finite speed set, `θ` rational, and suppose
> `f_C ≥ θ` on some interval `I` of length `δ > 0`. Then for **every** integer `w ≥ 1/δ`,
> `gap(C ∪ {w}) ≥ θ`.
>
> *Proof.* As `τ` traverses `I` (length `≥ 1/w`), `wτ mod 1` sweeps a full period, so some `τ* ∈ I` has
> `‖wτ*‖ = 1/2`. Then `f_{C∪{w}}(τ*) = min(f_C(τ*), 1/2) ≥ min(θ, 1/2) = θ`. ∎

**Consequence (the useful direction):** if `gap(C ∪ {w}) < θ` then `w < 1/δ`. So *the stranger is bounded*,
and searching the family `C ∪ {w}` for small gaps is a **FINITE, explicitly bounded** computation. This is the
gap-axis analogue of THM-518's measure-axis stranger-decoupling — but quantitative, with a computable bound.

`δ` is computed **exactly**: `{f_C ≥ θ}` is an intersection over `v∈C` of unions of closed intervals with
rational endpoints `[k/v + θ/v, (k+1)/v − θ/v]`, so `δ ∈ ℚ`.

## Theorem 1 — SGC'(13) on the single-perturbation family

> For every `j ∈ {1..13}` and every positive integer `w`, the set `S = ({1..13}\{j}) ∪ {w}` has
> **`gap(S) ∉ (1/14, 3/41)`**. Equivalently on this family: `gap > 1/14 ⟹ gap ≥ 3/41`.
> The extremal value `3/41` is attained at `{1,…,11,13,36}` (`j=12, w=36, τ=17/41`).

*Proof.* Take `θ = 3/41`. Exact `δ_j` and `W_j = ⌈1/δ_j⌉`:

| j | 1 | 2 | 3 | 4 | 5 | **6** | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| δ_j | 37/1066 | 31/2132 | 31/2706 | 25/3608 | 5/738 | **13/5412** | 15/1066 | 34/2583 | 28/2665 | 22/3731 | 1/123 | 1/287 | 1/246 |
| W_j | 29 | 69 | 88 | 145 | 148 | **417** | 72 | 76 | 96 | 170 | 123 | 287 | 246 |

By the lemma any band-hitter has `w < W_j ≤ 417`. Exact rational verification of **all** `j` and **all**
`w ≤ 417` finds **zero** sets with gap in `(1/14, 3/41)`. ∎

## Theorem 2 — tight classification on the same family (a piece of OPEN-Q-108)

> The **only** primitive sets of the form `({1..13}\{j}) ∪ {w}` with `gap = 1/14` are
> **the AP `{1,…,13}`** and **the Goddyn–Wong sporadic `{1,…,11,13,24}`**.

*Proof.* A tight set has `gap = 1/14 < 3/41 = θ`, so the same lemma bounds `w < W_j ≤ 417`; the enumeration
is therefore exhaustive. Exact check over all `j`, all `w ≤ 417` returns exactly these two sets. ∎

This is a **rigorous, exhaustive instance of the `{AP, GW}` conjecture** (HYP-2561 / OPEN-Q-108) — proved,
not just enumerated to a speed cutoff, because the cutoff is *derived* rather than assumed.

## What this does and does not do for OPEN-Q-108

**Does:** supplies the mechanism. Tight-locus finiteness is exactly "the stranger is bounded", and the lemma
makes that bound explicit and computable *per core* `C`, converting each family into a finite check.

**Does not:** give finiteness in general. The bound is `w < 1/δ_C`, and `δ_C` is **not** uniformly bounded
below over all 12-subsets — indeed `δ_C ≤ (1−2θ)/max(C) → 0`. So the bound degrades as the core grows. A
**uniform** lower bound on `δ_C` (equivalently the uniform fattening lemma, OPEN-Q-108) is still exactly what
is missing. My lemma shows *why* that is the crux, and in what quantitative form it is needed.

**Multi-stranger extension (next step, sketched).** For `C` plus `k` strangers all `≥ 1/δ`, the bad set of
each stranger inside `I` has measure `≈ 2θδ`, so a good `τ` survives whenever `2kθ < 1`, i.e. `k < 1/(2θ) =
41/6 ≈ 6.8`. Hence up to **6** strangers can be decoupled simultaneously, and band-hitters must have at
least one *small* stranger — giving a recursive finite search for the 2-, 3-, …-perturbation families.

Script (fully exact, reproducible): `04-computation/lrc14_sgc_prime_single_perturbation_theorem_macmini_S169.py`.
Prior: `07-reflections/sgc13-refuted-band-is-ostrowski-rungs-corrected-buffer-1over574-macmini-S169.md`.
