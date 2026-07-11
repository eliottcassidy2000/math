---
id: THM-707
title: "B5 = liveCount − deep-coverage penalty (the EXACT form of THM-671's B5 ≤ liveCount), via the partial alternating binomial identity T5(n) = −C(n−1,5); and the CLEAN-RULER reduction it yields — the single Lean obligation hB5 reduces to 'every residual covering family has a ruler q with a live multiplier and no multiplier covering ≥ 6 runners' (maxBand ≤ 5 ⟹ penalty 0 ⟹ B5 = liveCount > 0). The binding near-AP {1..12,26} is discharged by the clean pair-sum ruler q = 27 (maxBand 4). This fleshes out the two-routes-one-ladder architecture (HYP-5995): the clean-ruler condition is itself a seven-sector danger-arc occupancy statement, closing the fine-scale certificate onto the coarse-scale ladder."
status: PROVED (the identity T5(n) = Σ_{d=0}^5 (−1)^d C(n,d) = −C(n−1,5) is the elementary partial-alternating-binomial-sum identity, verified n = 0..9; the reformulation B5 = liveCount − Σ_{bandCount≥6} C(bandCount−1,5) and the clean-ruler corollary are one line each from it; Lean-formalizable). The DOWNSTREAM reduction "every residual covering family HAS a clean ruler" is the (now much more transparent) restatement of hB5 — still open, but reduced from a depth-5 signed Bonferroni to a shallow-coverage + liveness condition.
source: kind-pasteur-2026-07-11-S127 (cont.28) — fleshing out the HYP-5995 architecture
depends_on:
  - THM-671   # klein-S210: B5 ≤ liveCount (the ≤ direction; this is its exact form)
  - THM-701   # the peel (the unbounded half of the architecture)
related:
  - HYP-5995  # two-routes-one-ladder: B5 = the seven-sector moment ladder at scale 1/14
  - THM-703   # mac-mini: the coarse-scale factorial-moment majorant (same ladder, scale 1/7)
  - LRC14CompletionAudit  # hB5 = the single remaining Lean obligation
external: the partial alternating binomial sum Σ_{d=0}^k (−1)^d C(n,d) = (−1)^k C(n−1,k).
---

# THM-707 — B5 = liveCount − penalty, and the clean-ruler architecture

## Setup (the single Lean obligation)

`LRC14CompletionAudit`: LRC(14) is fully formalized, kernel-pure, foundational-axioms-only, modulo
`LRCUpTo13` (cited) and one obligation `hB5`: every residual covering family `v` has some `q` with
`B5(v,q) > 0`, where — writing `bandCount(v,q,p) = #{i : (v_i·p mod q) ∉ [q/14, 13q/14]}` (runners *not*
lonely at time `p/q`) —
> `B5(v,q) = Σ_{d=0}^{5} (−1)^d S_d`,  `S_d = Σ_{p=1}^{q−1} C(bandCount(v,q,p), d)`.
`THM-671` gives `B5 ≤ liveCount` and `B5 > 0 ⟹` a live multiplier ⟹ lonely.

## The exact identity

**Lemma (partial alternating binomial sum).** For `n ≥ 0`, `T5(n) := Σ_{d=0}^{5} (−1)^d C(n,d) = −C(n−1,5)`
(with `C(−1,5) = −1`). Hence
> `T5(0) = 1`, `T5(n) = 0` for `1 ≤ n ≤ 5`, and `T5(n) = −C(n−1,5) < 0` for `n ≥ 6`.
(Verified `n = 0..9`: `1, 0,0,0,0,0, −1,−6,−21,−56`.)

**Theorem (B5 exact).** Summing `T5` over the multipliers,
> **`B5(v,q) = liveCount(q) − PENALTY(v,q)`,  `PENALTY(v,q) := Σ_{p : bandCount(p) ≥ 6} C(bandCount(p)−1, 5)`.**

This is the *equality* form of THM-671's inequality `B5 ≤ liveCount`: the deficit is exactly the
`C(·−1,5)`-weighted mass of **deeply covered** multipliers (those where `≥ 6` of the 13 runners
simultaneously fall in the `1/7` danger arc `[−q/14, q/14]`). Multipliers with `1 ≤ bandCount ≤ 5`
contribute **zero** — the depth-5 Bonferroni is exact through 5-fold overlaps.

## The clean-ruler reduction (the architecture)

**Corollary (clean ruler).** If a ruler `q` has `liveCount(q) ≥ 1` and `max_p bandCount(v,q,p) ≤ 5`, then
`PENALTY = 0` and
> **`B5(v,q) = liveCount(q) > 0`.**

So the opaque depth-5 signed obligation collapses to a transparent, positive, Lean-friendly condition:

> **`hB5` ⟸ every residual covering family has a *clean ruler* `q` — one with a live multiplier and no
> multiplier covering `≥ 6` runners.**

Both clauses are **seven-sector occupancy statements** at the fine ruler `q`, because
`bandCount(v,q,p) = #{i : v_i p mod q ∈ the 1/7 danger arc}` is exactly the seven-sector danger-arc
occupancy:
- `liveCount ≥ 1` ⟺ some multiplier leaves the danger arc **empty** (a lonely time);
- `maxBand ≤ 5` ⟺ the danger arc **never holds ≥ 6** of the 13 runners (shallow coverage).

This is why HYP-5995's "two routes, one ladder" holds with equality here: the fine-scale certificate `B5`
*is* the coarse-scale seven-sector occupancy, and its positivity is governed by the same danger-arc counting.

## The binding case, discharged

The adversarial floor of `max_q B5` over residual covering families sits at the **near-AP `{1..12, 26}`**,
and it has a clean ruler: `q = 27 = 1 + 26` (a **pair-sum ruler**) with `maxBand = 4 ≤ 5` and `liveCount = 2`,
so `B5 = 2 > 0`. General bounded residuals have clean rulers with more room (`live = 6, maxBand ≤ 5` at
`q = 14` for many). (The *full* AP `{1..13}` has **no** clean ruler — best is `q = 14`, `live = 6`,
`penalty = 1`, `B5 = 5` — but it is **not residual**: it is discharged by the exact `M(AP) = 1/14` branch of
the grand assembly. So the clean-ruler reduction is asked only of the genuinely residual families, exactly
where it holds.)

## The assembled architecture (fleshed out)

> **`hB5` ⟸ [THM-701: peel far elements down to bounded cores] + [clean ruler on the bounded residual].**
>
> **Clean ruler on the bounded residual ⟸ `liveCount(q) ≥ 1` ∧ `maxBand(q) ≤ 5` for some `q`** — two
> seven-sector danger-arc occupancy conditions, binding at the near-AP, met by a pair-sum ruler.

The remaining open content is a single, transparent seven-sector statement: *every bounded residual covering
family admits a ruler at which the `1/7` danger arc is empty at one multiplier and never holds ≥ 6 runners.*
No depth-5 signed cancellation remains to reason about — THM-707 removed it.

## Why it matters

- **Removes the Bonferroni opacity.** `B5 > 0` looked like a delicate signed alternating-sum positivity; it
  is exactly `liveCount > (deep-coverage penalty)`, and on clean rulers it is simply `liveCount > 0`.
- **Sharpens THM-671 to an equality**, and the penalty term names precisely what obstructs positivity:
  multipliers of coverage-depth ≥ 6.
- **Lean-tractable.** The identity is an elementary induction; `measureFloor`-style real analysis is not
  needed for this node — the clean-ruler corollary is a decidable per-`q` check plus an existence claim.
- **Closes the fine scale onto the coarse ladder**, making HYP-5995's unification an equality, not an analogy.

## Files

`04-computation/lrc14_B5_architecture_kps_S127.py` (+`.out`): the identity `T5 = −C(n−1,5)`, the reformulation
`B5 = liveCount − penalty` (verified `== ` direct), and the clean ruler for the binding near-AP. Companion:
`lrc14_B5_moment_ladder_kps_S127.py`, `lrc14_B5_adversarial_floor_kps_S127.py`.
