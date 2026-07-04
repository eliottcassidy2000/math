# The multi-tightener shift obstruction: an extremity lemma, a compactness bound, and exactly why it does not close

*opus-2026-07-03-S61. The owner asked me to prove confinement for the multi-tightener case — THM-612's
open gap (`primitive tight ⟹ q*=14`, proved by mac-mini for one tightener, open for ≥2). I did not close
it. I extended the shift obstruction as far as it honestly goes: an extremity lemma and a compactness
bound for f=2, and a clean account of the exact place the argument runs out. Careful, given MISTAKE-097.*

> **Convergence note (read first).** Same-prompt, mac-mini reached this first and further: their S32
> anti-correlation IS the extremity lemma below, their S33 **Lemma D** (now in THM-612) gives the same
> "tighteners bounded by the even part" compactness AND the switch-point divisibility `w_i|w_j` that
> reduces `f=2` to a finite per-`U` check — which I did not get. So what follows is an *independent
> convergence*, not a new result; the canonical, stronger line is THM-612 Lemma D. My one useful artifact
> is the independent search confirmation; my compactness form `w_i ≤ 12 u_max` is a cleaner-but-weaker
> version of mac-mini's.

## What confinement needs, and what the shift argument sees

A primitive tight family with `q*=28` (`m=2`) splits as `S = 2U ∪ F`, `F` odd tighteners. On the U-loose
region `R = {g_U(2t) > 1/14}` — which is `(+1/2)`-invariant and positive-measure (`U` loose by LRC≤13) —
the even runners are strictly safe, so tightness forces the tighteners to **cover `R`**. The one clean
handle is that an odd `w` reflects: `‖w(t+1/2)‖ = 1/2 − ‖wt‖`, so `D_w ∩ (D_w+1/2) = ∅` — a tightener
covers at most one of each pair `{t, t+1/2}`. That reproves `f ≥ 2` (Lemma C) and gives `f ≥ 7·meas(R)`.

## The extremity lemma and the compactness bound (f=2, proved)

For `f=2` the reflection sharpens to a pointwise dichotomy. Both `t` and `t+1/2` lie in `R` and must be
covered; the reflection then forces, at every `t∈R`, **exactly one tightener `≤ 1/14` and the other
`≥ 3/7`** (extremity). Because `‖w_1 t‖` is continuous and skips the gap `(1/14, 3/7)` on `R`, every
component of `R` is *single-tightener* (entirely `w_1`-danger or entirely `w_2`-danger), the two types
swapped by `+1/2`. A single-type component sits in one danger arc, of length `1/(7w_i)`.

Applied to the component around `U`'s global-max point (length `≥ (M(U)−1/14)/u_max` by the Lipschitz
slope bound, THM-613's idea), this yields the **compactness bound**
`w_1, w_2 ≤ 12·u_max`: a `q*=28`, `f=2`, primitive tight family, if it existed, would have its tighteners
bounded by `12×` the even part's scale. The family is forced *compact* — a real reduction of the search
space, and (with the mod-28 grid) a step toward a finite check.

## Exactly why it does not close (the honest part)

The shift argument uses only one consequence of tightness: **`F` covers `R`.** That is necessary but far
from sufficient. If `U`'s loose region is a single arc (and its `+1/2` copy), covering `R` by two
tighteners is *structurally easy* — pick `w_1` dangerous on the arc, `w_2` on its shift. Nothing in the
covering-of-`R` picture forbids it. The obstruction that actually makes `q*=28` primitive tight families
non-existent lives in the conditions the shift argument **cannot see**:

- `M(S) ≤ 1/14` **everywhere off `R`** (the global loose condition, not just on `R`);
- the max is **attained at denominator exactly 28** (the `q*` condition);
- **primitivity** beyond "`E` even".

These are exactly the LRC-hard, integer/global constraints — the same wall as everywhere in this problem.
The single search over 938 structured even-block+odd-tightener families returns 0, matching mac-mini: the
families don't exist, but the reason is global, not the local covering the shift controls. And `m ≥ 3` is
untouched — the clean `1/2 − ‖wt‖` reflection is special to `m=2`.

## Status

- **Proved (elementary, verified):** `f ≥ 2`; `f ≥ 7 meas(R)`; the Extremity Lemma (f=2); single-tightener
  components; the compactness bound `w_i ≤ 12 u_max` (f=2). THM-614.
- **Confirmed:** independent exact-`M`/`q*` search, 938 families, 0 primitive tight with `q*>14`.
- **NOT proved (= the gap):** full confinement `primitive tight ⟹ q*=14`. The residual is the
  global-tightness + attained-denominator + primitivity conditions, and all `m ≥ 3`. Confinement stays a
  CONJECTURE (THM-612), and with it the tight-locus rigidity and the measure floor.

I flag the non-closure plainly. The contribution is a proper, clearly-scoped subset of the gap — the
shift obstruction pushed to its limit — not the proof requested.

Related: THM-612 (the tower + Lemma C, the gap this extends), THM-610 (`14|q*`), THM-613 (the slope idea
reused for the compactness bound), HYP-4062/kps (tight-locus rigidity — what confinement would give),
HYP-2913 (`g(14)≤3`, the other gap), MISTAKE-097 (why I flag the non-closure). Script:
`lrc14_confinement_setup_opus_S61.py`. THM-614, HYP-4066.
