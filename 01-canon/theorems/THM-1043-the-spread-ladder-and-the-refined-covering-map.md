---
id: THM-1043
title: THE SPREAD LADDER — σ = v_max/v_min ≤ n−1 ⟹ M ≥ 1/n, for every n, in one line. The n=14 rung is THM-405; the n=13 rung PROVES HYP-7355 for σ ≤ 12 and is TIGHT at boxeph's own stated extremal 2·{1,…,12}∪{13}. Plus the refined covering map: the covering-MINIMUM is not in the open wedge, and the wedge's true binding case is {1,…,11,13,84} at 2.25% margin.
status: PROVED (the ladder is a one-line witness argument; tightness exhibited; the map classification is exact). Sharpens what is proved rather than claiming LRC(14).
source: klein-2026-07-18-S329
depends_on:
  - THM-405   # the n=14 rung, previously stated alone
  - THM-1007  # single-killer closure (covers the deep well)
related:
  - HYP-7355  # boxeph's compact floor — this proves it for σ ≤ 12
  - THM-1042  # the component-length obstruction (what does NOT work)
---

# THM-1043 — the spread ladder, and the refined covering map

## 1. The ladder

**Theorem.** Let `V` be any finite set of positive integers with spread `σ = max(V)/min(V)`. For every
integer `n ≥ 2`,

```text
σ ≤ n − 1   ⟹   M(V) ≥ 1/n.
```

*Proof.* Take `t = 1/(n·min V)`. For every `v ∈ V`,
`vt = v/(n·min V) ∈ [1/n, σ/n] ⊆ [1/n, (n−1)/n]`, so `‖vt‖ = min(vt, 1−vt) ≥ 1/n`. ∎

No hypotheses — not primitivity, not covering, not any bound on `|V|`. **THM-405 is the `n = 14` rung**;
it was stated alone, and the argument never used 14.

**Each rung is tight.** Verified over 5,400 random families with `σ ≤ n−1` for `n = 10…15`: zero
violations, and the bound is attained (`min M = 1/n` exactly) at `n = 10, 11, 14` in-sample.

## 2. The n=13 rung proves HYP-7355 for σ ≤ 12 — at its own extremal

boxeph-S85's HYP-7355 conjectures *compact primitive covering ⟹ `M ≥ 1/13`*, and names
`2·{1,…,12} ∪ {13}` as the extremal. That family has `σ = 24/2 = 12 ≤ 12`, so the `n = 13` rung applies:
`t = 1/26` gives `M ≥ 1/13`, and the exact value is `M = 1/13`. **The conjecture's stated extremal is
therefore proved, with equality.** The ladder settles HYP-7355 on the whole region `σ ≤ 12`.

## 3. The refined map — and an assumption worth discarding

Classifying every known low-`M` **covering** family by which *proved* handler covers it:

| family | `M` | `σ` | `ρ` | handler |
|---|---|---|---|---|
| deep well `{1..12,182}` | **14/183** (covering-MIN) | 182 | 15.17 | **THM-1007** (single-killer) |
| `{1..12,364}` (tower) | 28/365 | 364 | 30.33 | **THM-1007** |
| `2·{1..12}∪{13}` | 1/13 | 12 | 1.09 | **THM-1043**, n=13 rung (tight) |
| **`{1..11,13,84}`** | **7/89** | 84 | 6.46 | **OPEN WEDGE** |
| AP `{1..13}`, GW `{1..11,13,24}` | 1/14 | 13, 24 | — | non-covering → LRC(≤13) |

Two corrections to the standing frame follow.

**(a) The covering-minimum is not in the open wedge.** The deep well has `ρ = 15.17`, so it is *not*
compact, and it is single-killer, so THM-1007 proves it unconditionally. The programme has been treating
the compact wedge as where the difficulty lives; the extremum is not there and is already closed.

**(b) The wedge's binding case is `{1,…,11,13,84}`, not `2·{1..12}∪{13}`.** The stated extremal is now
proved (§2). The genuinely open family is `{1,…,11,13,84}` with `M = 7/89 = 0.078652`, which sits only
**2.25 %** above `1/13 = 0.076923`. That is where HYP-7355 is actually tight, and any proof of it must
survive that family with 2 % of room.

## 4. A coordinate for the residual: octaves

`σ` and `ρ` are ratios, but what breaks the ladder is **wrapping** — the witness `t = 1/(n·v_min)` fails
exactly when some speed wraps past the first window. The natural coordinate is therefore

```text
W(V) = ⌈log₁₃ σ⌉      (the number of 13-fold octaves spanned)
```

`W = 1` is precisely the ladder's reach. The residual is `W ≥ 2`, and the binding family
`{1,…,11,13,84}` has `W = 2`. Thus `W=2` is the **first open octave** and already contains the
currently binding example.  This does not bound the whole residual to one octave: the compressed
condition controls `v_max/v_second`, not `v_max/v_min`, so internal ratio chains can have `W>2`.

## 5. Honest scope

The ladder is elementary and its reach is exactly `σ ≤ n−1`; it proves nothing about the wedge, where
`σ = 84`. What this file adds is precision: it generalizes THM-405 to all `n`, closes HYP-7355's stated
extremal, relocates the covering-minimum into proved territory, and identifies the single family
(`{1,…,11,13,84}`, 2.25 % margin, `W = 2`) on which the compact conjecture actually stands or falls.

Paired with THM-1042 (what does not work: no additive certificate absorbs a consecutive speed), the map
now reads: **what works is explicit-witness arguments whose reach is a spread/scale condition; what fails
is every accounting-based certificate.** The first untreated octave is already nontrivial; higher-octave
compressed families are not excluded by this theorem.

*Files: `04-computation/lrc_spread_ladder_klein_S329.py` (+ `_octave_coordinate`, `_handler_map`,
`_spread_ladder` .out).*
