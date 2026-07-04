# The covering-min needs razor-sharpness only at the deep well — every other covering family has 35/16287 of slack

*klein-2026-07-04-S129 (HYP-4090). Owner: investigate further and keep improving the proof. This
session ran into heavy fleet convergence (opus-S71 independently confirmed multi-swap ≥ 14/183;
mac-mini-S40 already had the 2-point equioscillation; mac-mini-S44 extended THM-618 to the whole
single-killer stratum). So most of what I re-derived is confirmation. The genuinely useful new
piece is a **reframing of the residual**, plus an honestly-documented dead-end.*

## The useful new piece: the residual only needs a NON-sharp bound

The covering-min crux is `M(S) ≥ 14/183` for every 13-element covering system `S`. Exhaustively over
the 509 minimal-tightener covering families (klein-S128 enumeration), exact `M`:

> **The deep well `{1,…,12,182}` is the UNIQUE family attaining `14/183`. Every other covering
> family is `≥ 7/89`.** (Second value `7/89`, gap `7/89 − 14/183 = 35/16287 ≈ 0.00215`.)

And the deep well is a **single-killer** family (`{1,…,12}` + the single killer `182 = 13·14`),
which is **already PROVED** by THM-618 (mac-mini) — the runner-1 ↔ killer equioscillation
`M({1,…,12,X}) = X/(13(X+1))`, minimized at `X=182`. mac-mini-S44 extended this to the whole
single-killer stratum (`0/8410` below `14/183`).

**Consequence for the open residual (m=2 folding).** The razor-sharp value `14/183` is attained at
exactly one family, and that family is already proved. So *every remaining argument* — the `m=2`
folding (opus), any multi-swap or split-killer bound — only has to clear `14/183` for families that
are **actually `≥ 7/89`**, i.e. with `35/16287` of slack in every inequality. **The residual is a
non-sharp bound.** This matters because the folding/pigeonhole arguments lose constant factors and
cannot hit a razor edge, but they *can* clear a target with a definite margin. Concretely: opus's
`m=2` folding needs `min_t max(g_E, Ψ) ≥ 14/183`; since the only `m=2` family that could be tight
(the deep well) is peeled off to THM-618, the folding may prove the slacker `≥ 7/89` (or just
`≥ 14/183` with room) for the rest.

Where `7/89` itself lives: `{1,…,11,13,84}` (drop-12 residue-liar), **also already proved** (kps
`LRCResidueLiar.lean`). So the two lowest rungs of the spectrum are both proved; the residual is
"nothing sneaks below the proved rungs," a gap statement with margin.

## Confirmations (credit to the fleet)

- **2-point equioscillation is universal.** Over the 509 families, the covering-min is a 2-point
  equioscillation between a small runner and a killer/tightener in `53/60` of the lowest cases (the
  rest are 3-point). This confirms mac-mini-S40's Chebyshev-equioscillation observation across the
  whole stratum. New detail: the **binding small runner varies** — runner-1 for the deep well
  (THM-618), runner-5 for `{1,…,11,13,84}` (kps residue-liar at `t*=37/89`), runner-2/3/5/7/11 for
  others — so THM-618 and the residue-liar are the *same* equioscillation mechanism with different
  binding runners. Confirms opus-S70's "parametric residue-formula ladder" route (not Delsarte).
- **Multi-swap ≥ 14/183, deep well unique** — opus-S71 established this independently the same
  window; my S128+S129 enumerations agree. The killer-offset (THM-618) is single-killer-specific;
  the multi-swap global argmax is a higher offset (opus-S71 / mac-mini non-convexity).

## A clean universal upper bound (minor)

For every covering family, `M(S) ≤ 1/(base_len(S)+1)` where `base_len` is the longest contiguous
`{1,…,r} ⊆ S` (proof: the sub-family `{1,…,r}` has `M = 1/(r+1)`, and adding runners only lowers
`M`). Verified `509/509`. The deep well saturates the extreme: longest possible base `r=12`
(covering forbids `{1,…,13}`, since then `q=14` is uncovered), so `M ≤ 1/13`, and THM-618's offset
`1/2379` lands it exactly at `14/183`. Equivalently, at the optimal witness `a/Q` the ratio
`Q/r ≤ 183/14 = 13.07`, maximized (uniquely) by the deep well.

## Dead-end (documented, per the "never waste ideas" mandate)

I conjectured the covering-min is **monotone in base-length** `r` (which would have explained
S128's swap-depth monotonicity via THM-618's mechanism). **REFUTED.** The per-`r` minima are *not*
monotone: `r=9 → 7/89 (0.0787)` but `r=10 → 1/12 (0.0833)` (up), `r=11 → 7/89` again, `r=12 →
14/183`. Base-length is **not** the order parameter; a family with a longer contiguous base can have
a *larger* covering-min than one with a shorter base (the killer offset depends on the killer/base
arithmetic, not just base length). S128's swap-depth monotonicity is a property of the
minimal-tightener enumeration per depth, not a base-length law. Do not pursue base-length as the
covering-min order parameter.

## Honest scope

No new theorem; no Lean. This session is mostly confirmation (heavy convergence with opus-S71,
mac-mini-S40/S44). The transferable takeaways: (1) the **residual needs only a non-sharp bound**
(`35/16287` slack; the sharp point is the proved deep well) — a genuine simplification for opus's
`m=2` folding; (2) the binding small runner varies, unifying THM-618 and the residue-liar as one
equioscillation; (3) base-length monotonicity is a dead-end. The open core is unchanged: the `m=2`
folding (opus), now known to have definite slack.

## Links

- Scripts: `04-computation/lrc14_equioscillation_structure_klein_S129.py` (+ `.out`),
  `lrc14_baselength_killeroffset_klein_S129.py` (+ `.out`). Exact. HYP-4090.
- Builds on / credits: THM-618 (mac-mini, single-killer + killer-offset), mac-mini-S40 (2-point
  equioscillation), mac-mini-S42/S44 (m≥3 shift-pigeonhole + single-killer stratum), opus-S70
  (Delsarte dead → parametric route), opus-S71 (multi-swap ≥ 14/183), kps `LRCResidueLiar`/
  `LRCOneSwapLadders` (proved low rungs), klein-S128 (isolated spectrum, gap 35/16287), klein-S126
  (even-part gap, the `m=2` folding input). Open core: the `m=2` folding (opus), with slack.
