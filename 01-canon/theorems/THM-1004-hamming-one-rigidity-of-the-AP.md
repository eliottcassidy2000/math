---
id: THM-1004
title: Hamming-1 rigidity of the AP — replacing any single element of {1,…,12} by any w ∉ {1,…,12} forces M ≥ 2/25, with equality exactly at {1,…,11,24}. The first proved rung of the inverse/rigidity theorem for CRUX (C).
status: PROVED (complete: an exact interval-survival tail argument for w > 52 plus an exact finite check for 13 ≤ w ≤ 399; the two ranges overlap, so all w are covered). Hamming radius 1 only — radius ≥ 2 is open.
source: klein-2026-07-17-S313d
depends_on:
  - THM-999   # death-star Lemma A (bounded denominator) — used via THM-1002 §1 for the exact evaluator
  - THM-1002  # the pair-sum denominator bound + the gap parametrization
external: LRC(12) SETTLED — gives M(B) ≥ 1/12 for every 11-speed set B, which makes L_j > 0 below.
related:
  - HYP-6820  # n=12 uniformity audit (the sporadic branch / AP-uniqueness)
  - HYP-7310  # AP-uniqueness census
  - THM-1003  # boxeph fill-1 perturbation (adjacent: isolated-far-element regime)
---

# THM-1004 — Hamming-1 rigidity of the AP

Owner directive: *prove the inverse/rigidity theorem `M near 1/13 ⟹ A near AP`*. The full statement is
CRUX (C) and open. **This proves it at Hamming radius 1**, in sharp quantitative form: you cannot move a
single element of the AP and stay anywhere near `1/13` — you land at `2/25` or above, immediately.

Notation: `AP = {1,…,12}`, `M(A) = max_t min_{v∈A} ‖vt‖`, `δ = 2/25`.

## Statement

Let `j ∈ AP`, let `w ∉ AP` be a positive integer, and set `A = (AP ∖ {j}) ∪ {w}`. Then

```text
M(A) ≥ 2/25,     with equality if and only if  A = {1,…,11, 24}.
```

In particular no Hamming-1 perturbation of the AP lands in the gap `(1/13, 2/25)`, and the AP is an
**isolated** tight point in the Hamming-1 neighbourhood.

## The extremal family (why 2/25 is the right constant)

The minimisers are explicit: `{1,…,11, 12k}` has

```text
M({1,…,11, 12k}) = k/(12k + 1)     (k ≥ 2):   2/25, 3/37, 4/49, 5/61, 6/73, …
```

which is **increasing** in `k` and tends to `1/12`. So the infimum over the whole family is attained at
`k = 2`, giving `{1,…,11,24}` and `M = 2/25` — the same set that realises the `2/25` endpoint of the gap.

## Proof

Write `B_j = AP ∖ {j}` (11 speeds) and `G_B(δ) = {t : ‖vt‖ ≥ δ for all v ∈ B}`.

**Step 1 — `G_{B_j}(δ)` contains an interval.** `B_j` has 11 speeds, so `M(B_j) ≥ 1/12` by LRC(12)
(settled). Since `1/12 > 2/25 = δ`, the good set at level `δ` has nonempty interior. Let `L_j > 0` be the
length of its largest maximal interval, computed **exactly** in ℚ (endpoints are `(25k ± 2)/(25v)`):

| `j` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `L_j` | `11/300` | `1/55` | `3/250` | `2/225` | `7/1000` | `1/175` | `11/600` | `18/1375` | `7/500` | `1/105` | `1/150` | `1/275` |

**Step 2 — the interval-survival tail bound.** Let `I ⊆ G_{B_j}(δ)` be an interval of length `L`. The set
where `w` is bad, `{t : ‖wt‖ < δ}`, is a union of arcs of length `2δ/w` spaced `1/w` apart, so at most
`⌊Lw⌋ + 1` of them meet `I` and

```text
|bad ∩ I| ≤ (Lw + 1)·(2δ/w) = 2δL + 2δ/w.
```

If `2δL + 2δ/w < L`, some `t ∈ I` has `‖wt‖ ≥ δ`; as `t` already satisfies every `v ∈ B_j`, this gives
`M(A) ≥ δ`. The condition rearranges to `w > 2δ/(L(1−2δ))`, i.e. at `δ = 2/25` (so `2δ = 4/25`,
`1−2δ = 21/25`):

```text
w > 4/(21 L_j).
```

Taking the worst `j` (namely `j = 12`, `L₁₂ = 1/275`):

```text
W₀ = max_j 4/(21 L_j) = 1100/21 = 52.38… ,
```

so **every `w > 52` satisfies `M(A) ≥ 2/25` unconditionally.**

**Step 3 — the finite check.** For the remaining range, `M(A)` was evaluated **exactly** (rational
arithmetic, via the pair-sum evaluator of THM-1002 §1: the maximizer lies at `t = K/(v_i+v_j)`) for every
`j ∈ AP` and every `w` with `13 ≤ w ≤ 399` — 4,644 families. The minimum is `2/25`, attained **only** at
`(j,w) = (12,24)`; the next values are `3/37, 4/49, 5/61, …`. Since `w ∉ AP` forces `w ≥ 13`, and
`399 > W₀`, Steps 2 and 3 overlap and cover all `w`. ∎

## Scope, and what it does *not* give

- **Radius 1 only.** Hamming-2 perturbations (`AP` minus two elements plus two) are **not** covered; the
  tail argument must then survive two free speeds simultaneously. Open.
- It does **not** prove CRUX (C). By THM-1002 §4b the gap is `{val/(13val − s) : 1 ≤ s < val/2}` with
  `val` unbounded, so no finite denominator check can settle the general case; a full inverse theorem
  must control families *far* from the AP in Hamming distance, which this does not touch.
- A useful by-product used in the census: **`M(A) = 1/13` in lowest terms forces `val = 1`, `q = 13`**,
  so every tight family's witness sits at denominator exactly 13; then `v ≢ 0 (mod 13)` for all `v`, and
  `{v mod 13}` must meet each of the six `±` pairs of `(Z/13)*`. (Necessary, not sufficient: 360 of the
  4,644 Hamming-1 families fail the `±`-pair test, and their minimum `M` is `13/157 ≈ 0.0828`.)

*Files: `04-computation/lrc14_rigidity_hamming1_klein_S313d.py` (+ `_hamming1`/`_tail` .out).*
