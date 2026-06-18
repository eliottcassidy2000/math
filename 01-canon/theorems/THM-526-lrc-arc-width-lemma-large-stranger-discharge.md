---
id: THM-526
title: The arc-width lemma for LRC — a non-tautological discharge of the large-parked-runner direction (if the core's widest level-1/14 safe arc exceeds the parked runner's danger-tooth width 1/(98m), then M(core ∪ {14m}) ≥ 1/14, with an explicit witness), closing LRC(14) unconditionally for the family of (size-12 cores ⊆ {1..13} covering 2..13) ∪ {14m} over all m
status: PROVED (the lemma: elementary interval-covering, proof below; tighter threshold m>1/(98·W) verified vs dense check). The sub-family closure is PROVED (lemma for m≥4 + exhaustive small-m check m≤3, 0 counterexamples). Co-derived with the S4 creative-reframing workflow (Angle B / "homogenization"), which obtained the same lemma with a weaker threshold m>1/(14·W).
source: mac-mini-2026-06-17-S4 (creative reframings of the last LRC(14) inequality)
depends_on:
  - THM-523   # covering-set reduction (the only hard case is covering 13-sets)
  - THM-524   # binding-pair reduction (M attained at a pair crossing)
  - THM-525   # easy-dominates-hard (the parked-runner / core framing)
related:
  - HYP-2578  # the reframings ledger: sawtooth-exclusion, off-grid fence, even-clearing, dip-law collapse, the clearing-crossing reformulation
  - HYP-2577  # the peeling-chain decomposition (this lemma supplies its missing rigorous large-m piece)
  - THM-518   # measure-side stranger-decoupling (this is the GAP-side, with an explicit threshold)
external: Lonely Runner Conjecture; proven for ≤7 runners (Barajas–Serra), open for 13.
---

# THM-526 — The arc-width lemma: discharging the large parked runner

**Context.** LRC(14) reduces (THM-523/524/525) to: every covering 13-set `S` has
`M(S) = max_τ min_{v∈S} ‖vτ‖ ≥ 1/14`. Such `S` contains a "parked" runner `w = 14m`
(`w ≡ 0 mod 14`, pinned at the observer at every grid time). Writing `S = A ∪ {w}`,
the genuine non-tautological difficulty was: bound how much the parked runner can pull
`M` below the core gap `M(A)`. This lemma settles the **large-`m`** half outright.

## The lemma (PROVED)

Let `A` be any finite set of positive integers with `M(A) > 1/14`, and let
`G_A = {τ ∈ [0,1) : ‖aτ‖ ≥ 1/14  ∀ a ∈ A}` be its **level-1/14 safe set** — a finite
union of closed arcs of positive total measure (positive because `M(A) > 1/14`). Let
`W(A)` be the width of its **widest arc**.

> **Lemma.** If `W(A) > 1/(98m)`, then `M(A ∪ {14m}) ≥ 1/14`.
> Hence for every `m ≥ M₀(A) := ⌊1/(98·W(A))⌋ + 1`, looseness is automatic.

**Proof.** The parked runner's danger set at level 1/14 is
`D = {τ : ‖14m·τ‖ < 1/14} = ⋃_k ( k/(14m) − 1/(196m),  k/(14m) + 1/(196m) )`:
`14m` open "teeth", each of full width `t := 1/(98m)`, with centers spaced
`1/(14m) = 7t` apart — so consecutive teeth are separated by **safe gaps of width `6t`**.
Let `I` be a widest arc of `G_A`, of width `W(A) > t`. Since each tooth has width `t < W(A)`,
`I` cannot lie inside any single tooth; as the teeth are isolated (gaps between them),
`I` must contain a point `τ₀ ∉ D`, i.e. `‖14m·τ₀‖ ≥ 1/14`. (Quantitatively, the safe
measure inside `I` is `≥ W(A) − (⌊W(A)/(7t)⌋+1)·t > 0` whenever `W(A) > t`.) Because
`τ₀ ∈ I ⊆ G_A`, also `‖aτ₀‖ ≥ 1/14` for all `a ∈ A`. Therefore
`min_{v ∈ A∪{14m}} ‖vτ₀‖ ≥ 1/14`, so `M(A ∪ {14m}) ≥ 1/14`. ∎

**Why it is non-tautological.** `W(A)` is computed from the **core `A` alone**, with no
reference to `S` or to `M(S)`; the conclusion `M(S) ≥ 1/14` is *derived*, with an explicit
witness `τ₀` (any safe point of the widest arc; e.g. the gap-center
`(2k+1)/(2·14m)`, where the parked runner is at distance ½). This is the gap-side analogue
of the measure-side stranger-decoupling (THM-518), but with an **explicit, computable
threshold** `M₀(A)` rather than an asymptotic `m → ∞`.

## What it closes (PROVED sub-family, unconditional)

For every size-12 core `A ⊆ {1,…,13}` that covers `2..13`, the worst widest-arc width is
`W(A) = 5/1848`, attained at the drop-6 core `A = {1,2,3,4,5,7,8,9,10,11,12,13}`, giving
threshold `M₀ = ⌊1848/490⌋+1 = 4`. So:

> **Corollary.** For every size-12 core `A ⊆ {1,…,13}` covering `2..13` and every `m ≥ 1`,
> `M(A ∪ {14m}) ≥ 1/14`. *Proof:* `m ≥ 4` by the lemma; `m ∈ {1,2,3}` by exhaustive exact
> check (0 counterexamples; in fact `M(A∪{14m}) = M(A) ≥ 2/23` there — the small parked
> runner does not bind). ∎

This is a genuine infinite sub-family of covering 13-sets for which **LRC(14) is now
proved unconditionally** — the canonical hard cores `{1..13}\{j} ∪ {14m}` among them.

## Honest scope and the remaining hard regime

`W(A)` is **not uniform** over all cores: it shrinks toward 0 as core speeds grow
(e.g. `A = {1..11,13,1400}` has `W(A) ≈ 6·10⁻⁴`), so `M₀(A)` is unbounded across the full
covering family. The lemma therefore does **not** prove LRC(14): the surviving hard regime
is exactly **large-speed cores in small-`m` resonance**, where the parked runner *is* the
binding partner and `W(A)` is too small for the threshold to bite. Two further pieces would
finish the program: (i) bound the essential core speeds (compactness — the gap-side of
THM-522; the S4 finite-reduction probe confirmed the inf `= 7/89` stabilizes once
speeds ≤ 84, but found the "min at smallest representative" map is not monotone, so the
bound is empirical); (ii) the small-`m` window per bounded core (finite). The sharpest
distilled target (S4 synthesis) is **clearing-crossing existence**: for the binding pair
`D = flank + w`, `M(S) = j/D` with `j = (flank·num mod D)`, and the non-tautological content
is that the crossing index `j` is forced `≥ D/14` — i.e. some pair-crossing clears all 13
runners at level ≥ 1/14.

(Runner-up exact identity, S4: `M({1,…,13, 14m}) = m/(14m+1) < 1/14` — a **14**-element
set, the cleanest proof that the 13-element cardinality bound is essential.)
