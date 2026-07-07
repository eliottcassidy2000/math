---
source: kind-pasteur-2026-07-07-S57
status: creative reduction — Route 1's density floor reduces, by elementary reverse-Markov, to
  a mean-max-gap bound inf_E E[maxgap] > 1/7 (comfortable margin), simplifying opus's PZ route
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - paley-zygmund
  - max-gap
  - three-distance
---

# The density floor reduces to a mean max-gap (elementary reverse-Markov)

Owner: *work the crux creatively, think Paley–Zygmund.* opus-S131's PZ reduction is
`μ_{1/7}(E) = P_x(maxgap{frac(e_i x)} > 1/7) ≥ E[U]`, `U = Σ_gaps (gap−1/7)_+`, with the open
lemma `inf_E E[U] > 0` — obstructed by the triple+ inclusion–exclusion overlaps. Working it, I
found a cleaner route to the *same* conclusion that collapses to a **single order statistic**.

## Two elementary reductions to the mean max-gap

- **Pointwise:** `U(x) = Σ (gap − 1/7)_+ ≥ (maxgap − 1/7)_+ ≥ maxgap − 1/7`, so
  **`E[U] ≥ E[maxgap] − 1/7`** (a lower bound on opus's `E[U]` by one term).
- **Reverse-Markov (bypasses `E[U]` entirely):** `maxgap ∈ [0,1]`, so for `a = 1/7`,
  `E[maxgap] ≤ a + P(maxgap > a)·(1−a)`, i.e.

  > **`μ_{1/7}(E) = P(maxgap > 1/7) ≥ (E[maxgap] − 1/7)/(1 − 1/7) = (7/6)(E[maxgap] − 1/7).`**

Either way, the whole density-floor positivity reduces to

> **`inf_E E[maxgap] > 1/7`** — a **mean max-gap** statement, one order statistic.

## Why this is the better target: a comfortable margin

Verified (grid + adversarial descent, 13-speed families):

| family | `μ_{1/7}` | `E[maxgap]` | `(7/6)(E[mg]−1/7)` |
|---|---|---|---|
| AP `{1..13}` | 0.442 | 0.211 | 0.080 |
| `2·AP` (tight, `M=1/14`) | 0.444 | 0.211 | 0.080 |
| prim-sat `2·{1..12}∪{13}` | 0.501 | 0.201 | 0.068 |
| adversarial min-`E[maxgap]` | 0.674 | **0.203** | 0.070 |
| random | 0.985 | 0.240 | 0.114 |

**`inf_E E[maxgap] ≈ 0.203`, margin `+0.06` above `1/7 = 0.143`** — robust to adversarial
descent, and *comfortable* (not razor-thin, unlike the raw loneliness `M = 1/14`). Note even
the **tight** family `2·AP` (`M=1/14`) has `E[maxgap] = 0.211`; the density floor is not tight
where the raw problem is, which is exactly why it is a workable route. Compared with `E[U]`
(full inclusion–exclusion, triples the obstruction), `E[maxgap]` is a single statistic whose
positivity margin is large and stable.

## Rigorous partials toward `inf E[maxgap] > 1/7`

Two clean lower bounds, both a little short — the shortfall is the actual content:

- **Length-biased:** `maxgap ≥ Σ g_i²` (max ≥ weighted mean) `≥ (Σg_i)²/k = 1/k = 1/13`
  (Cauchy–Schwarz). Empirically `E[Σg²] ≈ 0.14–0.16` — approaches but does not clear `1/7`.
- **Origin gap:** `maxgap ≥ gap_0 = min_i frac(e_i x) + (1 − max_i frac(e_i x))`, so
  `E[maxgap] ≥ E[min] + 1 − E[max] = 1 − E[range]`. For `k` uniform points `E[range] =
  (k−1)/(k+1) = 6/7`, giving `E[gap_0] = 2/(k+1) = 1/7` — **exactly the threshold**. Empirically
  `E[gap_0] ≈ 0.15–0.21 ≥ 1/7`, so the phases are *at least as spread as uniform* on average.

So the remaining content is precisely **"the max gap beats the length-biased / origin gap"**
by a fixed margin, uniformly in `E` — a **three-distance / discrepancy** statement about the
orbit `{frac(e_i x)}`, not an inclusion–exclusion of triples. This is a classical-flavoured
target (the max gap of `{n α}` is governed by continued fractions; here averaged over the
"multiplier" `x` for a fixed integer set `E`), and it carries a comfortable margin.

## Status

- **Rigorous & elementary:** `μ_{1/7} ≥ (7/6)(E[maxgap] − 1/7)` and `E[U] ≥ E[maxgap] − 1/7`.
- **Reduced target:** `inf_E E[maxgap] > 1/7` (mean max-gap; empirically `≈ 0.20`, margin `0.06`).
- **Partials:** `E[maxgap] ≥ 1/13` (length-biased); `E[maxgap] ≥ E[gap_0] ≥ 1/7`-borderline
  (origin gap / `E[range] ≤ 6/7`).
- **Open:** the uniform strict `inf E[maxgap] > 1/7` — a discrepancy/three-distance statement,
  cleaner than `E[U]`'s triples, and with margin.
- **Does NOT prove LRC(14).**

## Ledger

- **Files:** `lrc_maxgap_reduction_kps_S57.py` (+ out). HYP-4747.
- **Pointers:** opus-S131 (PZ reduction `μ_{1/7} ≥ E[U]`, the E[U] first-moment route);
  kps-S54 (`μ_{1/7}` census, `2·AP` tight); three-distance theorem / discrepancy of `{nα}`.
