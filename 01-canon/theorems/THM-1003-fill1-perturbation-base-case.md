# THM-1003 — The fill-1 perturbation lemma (base case): the far-element regime closes by a direct witness (boxeph-2026-07-17-S82)

**Status:** PROVED (elementary, exact) **and FORMALIZED in Lean, kernel-pure**
(`TournamentH7/LRCFill1Perturbation.lean`, `LonelyRunner.fill1_perturbation`, axioms
`[propext, Classical.choice, Quot.sound]`, no `sorry`/`native_decide`; corpus-registered).
Verified numerically on the deep-well and residue extremals and the single-killer family
{1..12,w} (`lrc_fill1_perturbation_boxeph_S82.py` + `.out`). Discharges
the "under-filled circle" clause of the resonance-fill crux
([[the-resonance-fill-profile-one-lens-for-every-lrc-face-boxeph-S81]], [[THM-998-farey-circle-deep-law]])
in its base case: **fill = 1 with a stranded runner dominant over the body.**

## Setup

`‖x‖ =` distance from `x` to `ℤ`. Threshold `1/N` (LRC(N)); `N=14` is the live case. A family
`V ⊂ ℤ_{>0}` is **`1/N`-lonely** if `∃ t, ∀v∈V: ‖v t‖ ≥ 1/N`. Fill `f_b(V) = #{v∈V : b∣v}`.

**Positioning (Lemma 0).** For `b ≥ 2`, `gcd(a,b)=1`, `t_0 = a/b`: if `b∣v` then `‖v t_0‖ = 0`;
if `b∤v` then `v a mod b ∈ {1,…,b-1}`, so `‖v t_0‖ ≥ 1/b`.

## Statement

> **Theorem (fill-1 base case).** Let `V ⊂ ℤ_{>0}`, `b ∈ {2,…,N-1}`, `gcd(a,b)=1`. Suppose
> circle `b` has **fill one**: a unique `v* ∈ V` with `b ∣ v*`, and every other `v ∈ V` has
> `b ∤ v`. Let `B := max{v ∈ V : b ∤ v}` (the body maximum). If
> `` `b·B ≤ (N-b)·v*` ``  (equivalently `B/v* ≤ (N-b)/b`),
> then `t = a/b + 1/(N v*)` witnesses `1/N`-loneliness: `min_{v∈V} ‖v t‖ ≥ 1/N`, and equals
> `1/N` at `v*` (and at any body runner attaining `B` with equality in the hypothesis).

## Proof

Minimal kick `s := 1/(N v*)`, `t = a/b + s`.

- **Stranded runner.** `b∣v* ⟹ v* a/b ∈ ℤ`, so `‖v* t‖ = ‖v* a/b + v* s‖ = ‖(ℤ) + 1/N‖ = 1/N`
  (using `1/N ≤ 1/2`). ✓

- **Body runner** `v` (`b∤v`). The hypothesis gives `v ≤ B ≤ (N-b)v*/b < N v*/2`, so
  `v s = v/(N v*) < 1/2`, whence `‖v s‖ = v s`. By the reverse triangle inequality on `ℝ/ℤ`
  (`‖x+y‖ ≥ ‖x‖ − ‖y‖`) and Lemma 0,
  `‖v t‖ ≥ ‖v a/b‖ − ‖v s‖ ≥ 1/b − v/(N v*).`
  This is `≥ 1/N` iff `1/b − 1/N ≥ v/(N v*)` iff `(N-b)v* ≥ b v`, which holds for every body
  runner because `v ≤ B` and `b B ≤ (N-b)v*`. ✓ ∎

## Reach and tightness (exact, verified)

- **Closes the single-killer far-element families.** `{1,…,12, w}` with `182 ∣ w` is covering;
  at `b = 13` (`f_13 = 1`, `v* = w`, `B = 12`) the condition is `156 ≤ w`, automatic since
  `w ≥ 182`. **Both covering-min extremals** land: deep well `{1..12,182}` at `b=13`
  (`t = 197/2548`, min `= 1/14`); residue `{1..11,13,84}` at `b=12`
  (`v*=84`, `B=13`, `156 ≤ 168`, `t = 33/392`, min `= 1/14`).
- **Tight.** For `{1..12,w}` at `b=13` the naive witness fails exactly below the boundary:
  `w=143` gives min `71/1001 < 1/14`; `w=156` (`=b·B`) gives min `= 1/14` (the body runner `12`
  also lands exactly on threshold). The inequality `b·B ≤ (N-b)v*` is sharp.
- **Unconditional sub-case.** If `b ≤ N/2` and `v* = max V`, then `B ≤ v*` and `(N-b)/b ≥ 1`,
  so the hypothesis is automatic — every fill-1 circle at `b ≤ 7` whose stranded runner is the
  maximum is dodgeable with no size check.

## Significance

- **The far-element regime of the crux is elementary.** klein-S288/kind-pasteur-S128 close the
  far-element regime analytically (disc_v; threshold `v* ≳ 83` for the deep well). THM-1003 closes
  it by a **direct rational witness** — cruder constant (`v* ≥ b·B/(N-b) = 156`) but no measure
  theory, and **Lean-trivial** (one reverse-triangle inequality; a candidate companion to
  `LRCLiveCountLonely`/`LRCC8Consecutive`). Complementary, not competing.
- **It is exactly the base of the deflation recursion.** Perturbing around circle `b` rescales
  the stranded runners `{v : b∣v}` by `1/b` into a sub-family `{v/b}` that must itself be lonely
  (Mode-B descent). Fill 1 = a **single** stranded speed = the trivial base case (one speed is
  always lonely). Fill `≥ 2` recurses into a smaller LRC problem.
- **The residual is pinned to the finish map.** What THM-1003 does NOT reach: (i) the
  **bounded/compact core** (`v*` too small, e.g. `{1..12,14}`), and (ii) **multi-killer** families
  (`{1..10,13,22,84}` — no single fill-1 circle with a dominant stranded runner). These are exactly
  the surviving strata of Route B (finish-map §"bounded-`Vmax` compact core" + `k≤8` multi-killer),
  where THM-724/726 rigidity and the uniform exposure bound live.

Related: [[THM-998-farey-circle-deep-law]], [[THM-724]], [[THM-733]] (kind-pasteur far-element),
[[THM-731]]/[[THM-732]] (disc_v analytic far-element), HYP-7315, HYP-7325.
