---
id: THM-668
title: The detuned-harmonic dispatch (d = 1) — every 13-family of the form g·H ∪ {δ} (g ≥ 2, |H| = 12, g ∤ δ) has M(v) ≥ 1/13 > 1/14, unconditionally (no covering, primitivity, or ratio hypotheses); closes the entire monad-S2 composition residual and the detuned-harmonic slice of the ratio > 13 realization residual
status: PROVED + FORMALIZED (Lean: TournamentH7/LRCDetunedDispatch.lean, lonely14_of_detuned, KERNEL-PURE [propext, Classical.choice, Quot.sound]; the formalization uses the triangle shortcut — Bezout branch shift c with c·δ ≡ d·⌊(g/d)/2⌋ (mod g) in the quarter window, then max of two branches ≥ 1/8 — instead of the full coset walk; conclusion min(1/13, 1/8) = 1/13, same as the prose). Machine-verified: explicit witnesses constructed on the monad-S2 residual instances and on ratio > 13 detuned constructions; exact rational clearances ≥ 1/13 confirmed (companion .out).
source: monad-explorer-2026-07-09-S3 (HYP-5727) — executing the S2 handoff ("the gcd-subgroup dispatch").
depends_on: []   # LRC(13): named citation per CLAUDE.md policy (Sungkawichai–Trakulthongchai; LRC13Citation.lean)
related:
  - LRCClusterGcd (kps-S18)   # the mechanism's in-repo precedent: 1/d-periodic margins + tooth pigeonhole, used contrapositively there (no-margin ⟹ gcd bound); here the positive direction at threshold 1/14
  - spread13_lonely / LRCHlargeRoute / LRCPureClusterCorner (kps-S28, death-star-S1)  # subsume ratio ≤ 13; this theorem eats the detuned-harmonic slice of the ratio > 13 residual they leave
  - THM-666 (mac-mini-S65, pair-sum ruler)  # cross-validation: the S2 residual witness 48/161 = 96/322 sits on the pair-sum ruler 322 = 154 + 168
  - HYP-5717 / THM-667  # the φ-interval composition whose named residual this dispatches
---

# THM-668 — the detuned-harmonic dispatch (d = 1)

## Statement

Let `v` be a family of 13 nonzero integer speeds of the form

> `v = g·H ∪ {δ}`,  where `g ≥ 2`, `H` is a set of 12 nonzero integers, and `g ∤ δ`.

Then

> **`M(v) = sup_τ min_i ‖v_i τ‖ ≥ min(1/13, 1/4) = 1/13 > 1/14`.**

In particular every such family is 14-lonely. No covering, primitivity, scale, or ratio
hypothesis is needed.

## Proof

**Step 1 (harmonic witness — the LRC(13) citation).** By LRC(13), applied to the 12
nonzero speeds `H`, there is a real `u*` with `‖m·u*‖ ≥ 1/13` for every `m ∈ H`.

**Step 2 (branch pigeonhole).** Consider the `g` branch times `τ_j = (u* + j)/g`,
`j = 0, …, g−1`. Every harmonic runner sees only `u*`:
`‖(g m)·τ_j‖ = ‖m·u* + m j‖ = ‖m·u*‖ ≥ 1/13` for all `j`. The detuned phases
`δ·τ_j = (δ u* + δ j)/g mod 1` run, as `j` varies, over the coset
`δu*/g + (gcd(δ,g)/g)·Z mod 1` — equally spaced points at spacing
`s = gcd(δ, g)/g ≤ 1/2` (the inequality because `g ∤ δ` gives `g/gcd(δ,g) ≥ 2`).
A coset of spacing `s` has a point within `s/2` of the antipode `1/2`, so some branch
`j*` has `‖δ·τ_{j*}‖ ≥ 1/2 − s/2 ≥ 1/4`.

**Step 3.** At `τ* = (u* + j*)/g`: the twelve harmonic clearances are `≥ 1/13` and the
detuned clearance is `≥ 1/4`, so `min_i ‖v_i τ*‖ ≥ 1/13`. ∎

## Scope and context

1. **It closes the monad-S2 residual entirely.** Every composition-invisible instance
   found by the S2 adversarial batteries has the form `14·(12 multipliers) ∪ {one
   detuned element}` (`d = 1` by construction: one perturbed harmonic) — e.g.
   `{14,…,70, 83, 98,…,182}` = `14·({1..13}∖{6}) ∪ {83}`: here the theorem's witness
   uses the LRC(13) time of `{1,…,5,7,…,13}` and one of 14 branches. (The measured
   exact `M = 2/23 = 0.0870` is consistent with the guaranteed `≥ 1/13 = 0.0769`.)
2. **It bites inside the ratio > 13 residual.** `spread13_lonely` (ratio ≤ 13 ⟹ lonely)
   subsumes detuned harmonics built on `14·{1..13}` (ratio exactly 13); the dispatch has
   no ratio hypothesis, so it also covers e.g. `14·{1,…,11, 20} ∪ {δ}` (ratio 20) —
   the detuned-harmonic slice of the "speeds-ratio > 13 only" realization residual
   (death-star-S1).
3. **Mechanism precedent (honest):** the 1/g-periodicity of the harmonic margin plus a
   pigeonhole over its `g` copies is exactly the engine of kps-S18's cluster-gcd ladder
   (`LRCClusterGcd`), used there contrapositively (a no-margin family bounds the gcd of
   any 11-subfamily). This theorem is the positive direction at the LRC(14) threshold,
   where it is three lines.
4. **The open generalization (d ≥ 2).** For `v = g·H ∪ D`, `|D| = d ≥ 2`, the branch
   vector `(δ₁ j, …, δ_d j) mod g` must hit the product safe box `∏[1/14, 13/14]`
   shifted by `(δ_i u*/g)`; the branch subgroup can be degenerate (e.g. all
   `δ_i` congruent mod `g` move diagonally), and then the rescue is one extra
   constraint on `u*` (a difference-speed condition), i.e. an LRC-type system with
   `12 + (d−1)` mixed-threshold constraints. Finite and decidable per `(g, D mod g)`;
   left open. Note covering forces structure on `D` (the detuned elements must carry
   every divisor `q` the harmonics miss), which shrinks the genuinely-open zone.
5. **Lean shape:** `LRC13Citation` + `nearInt` algebra + a `Finset.exists` pigeonhole —
   one page; the witness is explicit given the citation's `u*`.

## Verification & files

`04-computation/lrc14_detuned_dispatch_monad_S3.py` (+ `.out`): explicit witness
construction (exact rational `u*` from the exact 12-runner `M`-computation, branch
selection, clearance table) on the S2 residual instances and ratio > 13 constructions;
the union-coverage hunt ([φ-interval composition] ∪ [spread13] ∪ [this dispatch]) over
primitive covering adversaries, with the `d ≥ 2` frontier mapped.
