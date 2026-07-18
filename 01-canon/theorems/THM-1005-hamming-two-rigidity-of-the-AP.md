---
id: THM-1005
title: Hamming-2 rigidity of the AP — replacing ANY TWO elements of {1,…,12} by any two new speeds still forces M ≥ 2/25. The inverse/rigidity theorem for CRUX (C), proved at Hamming radius ≤ 2.
status: PROVED (complete, three regimes with overlapping ranges: exact finite check for both new speeds ≤ 68; joint interval-survival tail for both > 67.24; mixed-case tail for the large one > 52.53). Hamming radius ≤ 2 only — radius ≥ 3 open.
source: klein-2026-07-17-S313e
depends_on:
  - THM-1004  # Hamming-1 rigidity — the radius-1 case and the tail lemma reused here
  - THM-1002  # pair-sum denominator bound (exact evaluator) + gap parametrization
  - THM-999   # death-star Lemma A (bounded denominator)
external: LRC(11)/LRC(12) SETTLED — give M(B) ≥ 1/11 (10 speeds) and ≥ 1/12 (11 speeds), which force L > 0 below.
related:
  - HYP-6820  # n=12 uniformity audit / AP-uniqueness (the sporadic branch)
  - HYP-7310  # AP-uniqueness census
---

# THM-1005 — Hamming-2 rigidity of the AP

Owner directive: *prove the inverse/rigidity theorem `M near 1/13 ⟹ A near AP`*. THM-1004 proved it at
Hamming radius 1. **This proves radius 2**: you may move *two* elements of the AP anywhere and you still
cannot approach `1/13` — you land at `2/25` or above.

`AP = {1,…,12}`, `M(A) = max_t min_{v∈A} ‖vt‖`, `δ = 2/25`.

## Statement

Let `j ≠ k ∈ AP` and let `w₁ < w₂` be positive integers not in `AP`. Put
`A = (AP ∖ {j,k}) ∪ {w₁,w₂}`. Then

```text
M(A) ≥ 2/25.
```

Hence **no family within Hamming distance 2 of the AP lies in the gap `(1/13, 2/25)`**, and combined with
THM-1004 the AP is an isolated tight point of the whole radius-≤2 neighbourhood.

## Extremals (two distinct ones)

`2/25` is attained, and by *inequivalent* configurations:

```text
{1,…,11, 24}                              (radius 1, THM-1004)
{1,2,3,5,7,8,9,10,11,12, 17, 19}          (radius 2: remove 4,6 — add 17,19)
```

so the bound is sharp at radius 2 as well, and the extremal set is **not** unique — a fact any eventual
full inverse theorem must accommodate.

## Proof

Write `B = AP ∖ {j,k}` (10 speeds) and `G_S(δ) = {t : ‖vt‖ ≥ δ ∀ v ∈ S}`.

**The tail lemma (from THM-1004 Step 2).** If `I ⊆ G_S(δ)` is an interval of length `L` and `w` is a new
speed, the bad set `{‖wt‖ < δ}` meets `I` in measure `≤ 2δL + 2δ/w`. Hence for `r` new speeds
`w_1,…,w_r`, a point of `I` survives all of them as soon as

```text
2rδL + 2δ·Σ_i (1/w_i)  <  L.
```

**Regime 1 — both new speeds large.** With `r = 2` and `δ = 2/25` the condition reads
`Σ 1/w_i < 17L/4`, which holds whenever both `w_i > 8/(17L)`. Since `B` has 10 speeds, LRC(11) gives
`M(B) ≥ 1/11 > δ`, so `L_{j,k} := ` (largest interval of `G_B(δ)`) `> 0`, computed exactly in ℚ over all
`C(12,2)=66` pairs:

```text
W_joint = max_{j,k} 8/(17 L_{j,k}) = 8000/119 = 67.23…
```

So if `w₁, w₂ > 68` the claim holds unconditionally.

**Regime 2 — mixed (one small, one large).** Fix `w₁ ≤ 68` and set `B' = B ∪ {w₁}` (11 speeds). LRC(12)
gives `M(B') ≥ 1/12 > δ`, so `G_{B'}(δ)` has an interval of length `L' > 0`. Applying the tail lemma with
`r = 1` (as in THM-1004) the remaining speed is absorbed once `w₂ > 4/(21 L')`. Over all
`66 × 56 = 3,696` pairs `(j,k)` and `w₁ ∈ [13,68]` (no degenerate `L' = 0` case occurs):

```text
W₂ = max 4/(21 L') = 11400/217 = 52.53…   (worst: j,k = 4,6 and w₁ = 38)
```

So any `w₂ > 68 (> W₂)` is absorbed.

**Regime 3 — both small.** For `13 ≤ w₁ < w₂ ≤ 68` and all 66 pairs `(j,k)` — **101,640 families** — the
predicate `M(A) ≥ 2/25` was verified by exact rational arithmetic (early-exit scan over the pair-sum
denominators of THM-1002 §1). **Zero violations**; the minimum is exactly `2/25`.

Since `w ∉ AP` forces `w ≥ 13`, and `68 > max(W_joint, W₂)`, the three regimes overlap and exhaust all
`(w₁,w₂)`. ∎

## Scope — what is still open

- **Radius ≥ 3 is not covered.** The tail lemma degrades as `r` grows: it needs `2rδL < L`, i.e.
  `r < 1/(2δ) = 6.25`, so this argument can in principle reach radius ≤ 6, but each radius needs its own
  (rapidly growing) finite check — `C(12,r)` pairs times a `w`-box. Radius 3 is the next rung.
- **It says nothing about families far from the AP**, which is where CRUX (C) actually lives: by
  THM-1002 §4b the gap is `{val/(13val−s) : 1 ≤ s < val/2}` with `val` unbounded, so a full inverse
  theorem needs control of *all* configurations, not a neighbourhood of the AP.
- Together THM-1004 + THM-1005 establish the **local** form of the stability bound
  (`non-AP ⟹ M ≥ 2/25` within Hamming distance 2); the **global** form remains CRUX (C) / HYP-6820.

*Files: `04-computation/lrc14_rigidity_hamming2_klein_S313e.py` (+ `_hamming2`, `_hamming2_tail`,
`_hamming2_ext`, `_hamming2_W2` .out).*
