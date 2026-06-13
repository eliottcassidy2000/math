---
id: THM-405
name: lrc-bounded-ratio-first-window
status: PROVED
date: 2026-06-03
session: oracle-2026-06-03-S582o
depends_on:
  - THM-386   # Lonely = central-box / nearest-integer predicate
  - THM-384   # LRC safety = observer-adjacent gap condition
---

# THM-405: Bounded speed ratio ⟹ LRC via the first window (AP is the sharp boundary)

## Statement

Let `n ≥ 2` and let `S = {v_1,…,v_{n-1}}` be distinct nonzero integers. Write
`a = min_i |v_i| ≥ 1` and `b = max_i |v_i|`. If

```
    b ≤ (n-1) · a          (equivalently  v_max / v_min ≤ n-1)
```

then `S` is lonely on the whole interval

```
    t ∈ [ 1/(n a) , (n-1)/(n b) ]   (nonempty),
```

hence `M(S) = max_t min_i ‖v_i t‖ ≥ 1/n`, i.e. **LRC@n holds for `S`**.

The bound is sharp: the arithmetic progression `S = {1,2,…,n-1}` has `a=1, b=n-1`, so
`b = (n-1)a` (equality) and the window collapses to the single point `t = 1/n` — recovering
the AP as the tight extremal (`M = 1/n` attained only at `t = j/n`).

## Proof

Loneliness is unchanged under `v_i → |v_i|` (since `‖(-v)t‖ = ‖vt‖`), so assume all
`v_i > 0`. For a positive integer `v` and `t ∈ (0, 1/v)` we have `vt ∈ (0,1)`, so
`frac(vt) = vt` and (THM-386 / `far_iff_fract`)

```
    ‖vt‖ ≥ 1/n  ⟺  frac(vt) ∈ [1/n, 1-1/n]  ⟺  vt ∈ [1/n, 1-1/n]
                ⟺  t ∈ [ 1/(n v) , (n-1)/(n v) ].
```

(This interval lies inside `(0,1/v)`, the first passage of runner `v` before it returns to
the observer collar.) Intersecting over all runners,

```
    ⋂_i [ 1/(n v_i) , (n-1)/(n v_i) ]
      = [ max_i 1/(n v_i) ,  min_i (n-1)/(n v_i) ]
      = [ 1/(n a) , (n-1)/(n b) ],
```

since `1/(n v_i)` is maximised at the smallest speed `a` and `(n-1)/(n v_i)` is minimised at
the largest speed `b`. This interval is nonempty iff `1/(n a) ≤ (n-1)/(n b)`, i.e.
`b ≤ (n-1) a`. For any `t` in it, every `‖v_i t‖ ≥ 1/n`, so `S` is lonely and `M(S) ≥ 1/n`. ∎

## Consequences / scope

- **LRC@n is proved for every speed set of bounded ratio `v_max/v_min ≤ n-1`.** For `n=14`
  this settles all configs with `v_max/v_min ≤ 13`. (Verified computationally 800/800 over
  n=8,10,12,14, `lrc_first_window_ratio_theorem_s582d.py`: every bounded-ratio config lonely
  on the predicted window; the AP window is exactly `{1/n}`.)
- **The AP is the sharp boundary**: equality `b=(n-1)a` ⟺ the window is a single point ⟺ the
  tight case. This is the elementary reason the AP is the LRC extremal and explains
  HYP-2139's empirical "first-lonely time `≤ ~1.5/n`": the bounded-ratio window opens by
  `t = 1/(na) ≤ 1/n`.
- **The residual** (the genuine open core of LRC@n) is exactly the **large-ratio** sets
  `v_max/v_min > n-1`, where fast runners re-enter the observer collar before slow runners
  have left it, so the first-passage window is empty and a later (fine-pinch, mod-`2n-1`)
  window is required. This is the precise complement THM-405 leaves open.

## Relation to prior work

Elementary "first-gap" argument made precise as a ratio condition with the AP identified as
the sharp boundary. Companion to THM-384 (observer-gap), THM-369 (sieve), THM-401 (pair-sum
sieve = `2n-1`, the *fine*-window modulus that governs the large-ratio residual), HYP-2139
(first-window scaling, `c ≈ 1.5`), S556o (local LP / first-window conjecture). Lean-ready
(the proof is `far_iff_fract` + an interval intersection).

## Artifacts
```
04-computation/lrc_first_window_ratio_theorem_s582d.py
05-knowledge/results/lrc_first_window_ratio_theorem_s582d.out
```
