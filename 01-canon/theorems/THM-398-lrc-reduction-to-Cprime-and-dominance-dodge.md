---
id: THM-398
title: LRC reduces to C′ (multiple-of-n ⟹ loose); the dominance-dodge lemma
status: PROVED (the reduction and the dodge lemma); C′ itself remains CONJECTURE
source: opus-2026-06-03-S572
depends_on:
  - THM-369   # divisibility / 1-clock witness
  - LRC(n-1)  # proven in the literature for n-1 <= 13
related:
  - HYP-2102  # the reduction, discovered S571
  - HYP-2097  # the 64-class worry container
  - HYP-2095  # the lift lemma (paired/anchored split)
---

# THM-398 — LRC reduces to C′, and the dominance-dodge lemma

## 0. Setup and definitions

Throughout `n ≥ 2`; a **speed set** `S = {v_1,…,v_{n-1}}` is `n-1` distinct positive
integers (the moving runners; the observer is runner `0`). For `t ∈ ℝ/ℤ` write
`‖x‖ = dist(x, ℤ)`. Define

```
M(S) = max_{t} min_{v∈S} ‖v t‖,
```

the maximal gap (max-collar). **LRC(n)** is the statement `M(S) ≥ 1/n` for every such
`S`. It is a theorem for `n ≤ 13` (literature) and open at `n = 14`.

The **safe set at level `1/n`** is `G(S) = { t : ‖v t‖ > 1/n for all v ∈ S }`, an open
subset of the circle. Note `μ(G(S)) > 0 ⟺ M(S) > 1/n` (one direction: a point of `G`
has all `‖v t‖ > 1/n` so `M > 1/n`; conversely the optimum `t*` with
`min ‖v t*‖ = M > 1/n` has an open neighbourhood still in `G` by continuity). We call
`S` **loose** if `M(S) > 1/n` and **tight** if `M(S) = 1/n`.

> **C′ (conjecture).** If `n | v` for some `v ∈ S`, then `M(S) > 1/n` (`S` is loose).

C′ was observed in S564 as a *property of tight sets* ("tight ⟹ no multiple of `n`").
Here we prove its converse direction is the entire conjecture.

## 1. Lemma A (the 1-clock witness; restatement of THM-369)

**If no element of `S` is a multiple of `n`, then `M(S) ≥ 1/n`.**

*Proof.* Evaluate at `t = 1/n`. For each `v ∈ S`, write `r = v mod n`. Since `n ∤ v`,
`r ∈ {1,…,n-1}`, so `‖v/n‖ = ‖r/n‖ = min(r, n-r)/n ≥ 1/n`. Thus
`min_{v∈S} ‖v/n‖ ≥ 1/n`, whence `M(S) ≥ 1/n`. ∎

## 2. Theorem (the reduction): **C′ ⟹ LRC(n)**

*Proof.* Let `S` be any speed set.
- If no `v ∈ S` is a multiple of `n`: `M(S) ≥ 1/n` by Lemma A.
- If some `v ∈ S` is a multiple of `n`: `M(S) > 1/n ≥ 1/n` by C′.

Either way `M(S) ≥ 1/n`; this is LRC(n). ∎

**Remark (why this is a genuine reduction).** The hypothesis C′ constrains only the
*proper subclass* of speed sets that contain a multiple of `n` — precisely the sets
on which the elementary `1/n`-clock witness fails. All of LRC's difficulty is thereby
isolated into a single structural statement about the distinguished runner that
vanishes on the whole `n`-clock. (Verified: 0 tight-with-multiple over exhaustive
small boxes and large samples, `n = 4..14`, every multiplier size — HYP-2102.)

**Remark (the Vitali handoff — S551o).** This split *is* the "Vitali wall" of S551o:
LRC = (positive-measure bulk, settled by *measure*) ∪ (measure-zero core, settled by
*construction*), with the Vitali set marking the handoff. THM-398 *locates the
handoff*: it is the equation **`n | v`**. Configs with `n ∤ v` are handled by Lemma A
— the `t=1/n` construction, which is *measure-blind* and so reaches the measure-zero
core (the worry-set has no multiple of `n`, S564). Configs with `n | v` are pushed to
the *measure* side (C′). On that side the danger of `v=nw` is a **periodic,
bounded-eccentricity arc family** (a genuine Vitali cover), so the **Vitali covering
lemma / Lebesgue density theorem** is the natural tool — see §3–§4 and HYP-2104.

## 3. Lemma B (the dominance-dodge): a partial proof of C′ — **and more**

For a speed set `S` and `v ∈ S` write `V'(v) = max(S \ {v})`.

> **Lemma B.** Assume LRC(n−1). If some `v ∈ S` satisfies `v > (n-1)·V'(v)`, then
> `M(S) > 1/n` (`S` is loose). *(No divisibility hypothesis: this holds for any
> dominant runner, multiple of `n` or not.)*

*Proof.* Put `S' = S \ {v}`, with `|S'| = n-2` distinct positive integers, and
`V' = V'(v) = max S'`. By LRC(n−1) applied to `S'` (proven for `n-1 ≤ 13`),
`M(S') ≥ 1/(n-1)`, so there is `t₀` with `‖u t₀‖ ≥ 1/(n-1)` for every `u ∈ S'`.

*The S′-window.* For `|t - t₀| < δ` with
```
δ := 1 / ( n (n-1) V' ),
```
each `u ∈ S'` satisfies `‖u t‖ ≥ ‖u t₀‖ - u·|t-t₀| > 1/(n-1) - V'·δ
   = 1/(n-1) - 1/(n(n-1)) = 1/n`.
So on the open interval `I = (t₀ - δ, t₀ + δ)` every runner of `S'` is `> 1/n`.

*Dodging v.* The danger set of `v` is `D_v = { t : ‖v t‖ < 1/n }`, a union of open
arcs of radius `ρ := 1/(n v)` centred at the points `k/v`, separated by safe gaps.
An interval can lie inside `D_v` only if it lies inside a *single* arc (a gap is
`v`-safe, so an interval spanning a gap meets the `v`-safe set). Hence if
`|I| = 2δ > 2ρ`, i.e. **`δ > ρ`**, then `I ⊄ D_v`, so `I` contains a point with
`‖v t‖ > 1/n` as well — and an open neighbourhood of it, since `I ∩ (circle ∖ D̄_v)`
is open and nonempty.

*The threshold.* `δ > ρ ⟺ 1/(n(n-1)V') > 1/(nv) ⟺ v > (n-1)V'`, exactly the
hypothesis. On the resulting open set every runner of `S` exceeds `1/n`, so
`μ(G(S)) > 0` and `M(S) > 1/n`. ∎

**Corollary B1 (large multiples).** If `n | v` and `v > (n-1)·V'(v)` then `S` is
loose. In particular every multiple-of-`n` config whose multiple dominates the rest
is handled — and for `n = 14` this draws on the *proven* `LRC(13)`.

**Corollary B2 (one dominant runner ⟹ loose).** Any speed set with a single runner
exceeding `(n-1)×` all others is loose, regardless of divisibility. *(Verified:
1500/1500 loose at each of `n = 6,8,10,12,14`, dominant runners chosen with arbitrary
residues — `lrc_dodge_formalization_s572.py`.)* This is the **clean general fact**
that the divisibility-flavoured C′ partial is a special case of: the dodge is about
**dominance**, not about multiples of `n`.

## 4. The sharpened criterion and the residual dichotomy

Inspecting the proof, Lemma B only used that *some* component interval of `G(S')`
(the safe set of the `n-2` runners at level `1/n`) is longer than one danger arc of
`v`. This is weaker than `v > (n-1)V'`:

> **Criterion B′.** If `G(S\{v})` has a connected component of length `> 2ρ = 2/(nv)`,
> then `S` is loose. *(Proof identical: that long interval cannot sit inside one
> `v`-arc, so it meets `G(S)`.)*

**Vitali-covering iff (HYP-2104).** Because `G(S\{v})` is a finite union of intervals
and `D_v`'s gaps are open and nonempty, *an interval of `G(S\{v})` lies in `D_v` iff
it fits inside a single arc*. So Criterion B′ is the Vitali-covering direction, and a
multiple-of-`n` config is tight (measure-0) *only if* every component of `G(S\{v})` is
both short **and** arc-aligned. Quantified: B′ already proves looseness for **72%
(n=6) → 96.8% (n=14)** of multiple-of-`n` configs, and the all-short residual is
**never tight** (0 across n=6..14) — `lrc_vitali_covering_residual_s573.py`.

This yields the honest **dichotomy** for a multiple-of-`n` config `S` (`v = nw`):

1. **Long-interval case** — `G(S\{v})` has a component longer than `2/(n·nw) = 2/(n²w)`:
   `S` is loose by Criterion B′. *(Always holds when `v` dominates; Lemma B.)*
2. **All-short case** (the **residual**) — every component of `G(S\{v})` is `≤ 2/(n²w)`:
   the dodge of a single interval is not guaranteed. Here `S` is still loose
   (verified), but for the equidistribution reason
   ```
   μ(G(S)) = μ(G(S')) − μ(G(S') ∩ D_v),   μ(G(S') ∩ D_v) ≈ (2/n)·μ(G(S')) < μ(G(S')),
   ```
   i.e. the single arithmetic progression of thin arcs `{k/(nw)}` (period `1/(nw)`,
   total danger `2/n`) cannot *align to cover* the fixed union of intervals `G(S')`.
   Proving this non-covering is the remaining open core of C′ — a discrepancy /
   three-distance statement about one AP against a fixed open set.

**The residual is sharpest at `v = n` (`w = 1`):** the arcs sit on the `n`-clock
points `k/n` (radius `1/n²`), and one must find a point in some clock gap
`(k/n + 1/n², (k+1)/n − 1/n²)` where all `n-2` other runners exceed `1/n`. The gap
*midpoint* `(2k+1)/(2n)` (where `v=n` is maximally safe, at `1/2`) suffices only about
half the time (773/1499 at `n=6`), so the witness genuinely ranges over the gap
interior — confirming the residual is a true interval-search, not a fixed sub-clock.

## 5. Status ledger

| Statement | status |
|---|---|
| Lemma A (no-multiple ⟹ `M ≥ 1/n`) | **PROVED** (THM-369) |
| Reduction **C′ ⟹ LRC(n)** | **PROVED** (§2) |
| Lemma B / Cor B1 (dominant multiple dodge, uses LRC(n−1)) | **PROVED** (§3) |
| Cor B2 (one dominant runner ⟹ loose) | **PROVED** (§3), verified 1500×5 |
| Criterion B′ (long component ⟹ loose) | **PROVED** (§4) |
| C′ in the all-short / small-multiple case | **OPEN** (the equidistribution residual) |

So LRC(14) now sits on a single open assertion: *the thin evenly-spaced danger arcs
of a multiple of 14 cannot cover the safe set of the other twelve runners* (the
all-short case), with the long-interval/dominant case fully proved from the literature's
LRC(13).

**Artifacts:** `04-computation/lrc_dodge_formalization_s572.py` (+`.out`),
`lrc_lift_Cprime_residual_s571.py`, `lrc_lift_lemma_measure_bound_s571.py`.
See reflection `07-reflections/lrc-formalizing-the-Cprime-reduction-and-dominance-dodge-s572.md`.
