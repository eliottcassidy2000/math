# Uniform looseness is LRC-hard — but the far-peel is a measure bound, and linear in p

*kind-pasteur-2026-07-03-S39. Asked to prove the uniform looseness bound. I could not — in full
it is the covering case of LRC(14), open. But making the measure route rigorous produced a
clean comb-discrepancy bound that is the far-peel's measure twin, closes the covering-min
extremizer by measure, and shows the far-peel threshold is linear in the safe-component count,
not quadratic in the magnitude. Honest partial.*

## What "uniform looseness" actually is

The frontier map (S38) left two computational gaps; one is *uniform looseness*: primitive
covering ⟹ `M ≥ (1+c)/14` for a uniform `c > 0` (mac-mini's margin map: `M/(1/14) ∈
[1.06,1.11]`). This is not lighter than LRC(14): it **is** the covering case of the conjecture
(the bound `M ≥ 1/14` bundled with the strict gap). Non-covering is the sieve; covering is the
whole remaining bound. So "prove uniform looseness" = "prove LRC(14) for covering families." I
did not.

## The one rigorous piece: the far-peel is a measure bound

The measure route gives a clean, rigorous lower bound. For a family `v` and any runner `r`, let
`safe' = {t : ‖v_j t‖ ≥ 1/14 ∀ j ≠ r}`, with measure `μ'` and `p` components, and let `D_r` be
`r`'s danger set (measure `1/7`, a `(1/v_r)`-periodic comb):

> **`μ_v = μ' − Leb(safe' ∩ D_r) ≥ (6/7)μ' − 2p/(7 v_r).`**

*Proof.* `D_r` has density `1/7` and period `1/v_r`; over each of the `p` components of `safe'`
its measure differs from `(1/7)·(length)` by at most two partial teeth `= 2/(7v_r)`; sum over
components. ∎

**Corollary.** If `v_r > p/(3μ')`, then `μ_v > 0`, so `M(v) > 1/14`. And `μ' > 0` always, by
LRC(13): the other 12 runners have `M ≥ 1/13 > 1/14`, so their safe set has positive measure.

This is exactly the far-peel (opus LRCFarPeelCore / kps `far_peel_lonely`) rederived on the
Lebesgue measure instead of the rational region — the same rate lemma `length(inter) ≤
2r·length + (#comp)·4r/w`, read as a discrepancy. It closes the **dominant-far** case.

## Two things this makes precise

**1. The far-peel threshold is LINEAR in `p`, not quadratic in `V`.** The corollary's threshold
is `v_r > p/(3μ')`, with `p` the *actual* number of safe components of the other 12. My S33
piece-count bound `p ≤ 1 + 2ΣB` is what turned this into the `~V²` threshold I quoted for
`far_peel_lonely_of_cite` — but that bound is loose. The true content is linear: a runner
exceeding `p/(3μ')` peels. Worth flagging for the Lean far-peel: the `V²` is an artifact of the
piece-count bound, not the mechanism.

**2. It closes the covering-min extremizer by measure.** The deep well `{1,…,12,182}` has
`safe'({1..12})` with `μ' = 0.034` and only `p = 12` components, so the threshold is `12/(3·
0.034) = 117`, and `v_r = 182 > 117`: the bound gives `μ_v ≥ (6/7)(0.034) − 24/(7·182) = 0.0104
> 0`, matching the measured `μ_v = 0.0239`. So the deep well — previously census-closed at
`q = 40` — is now independently **measure-closed**.

## Why it stops at dominant-far

`p` is **not** `≤ C·n`. It grows with magnitude: random 12-speed families to magnitude 3000
reach `p ≈ 4700 ≈ Σv`. The deep well's `p = 12` is small only because its base `{1,…,12}` is
small. So the threshold `v_r > p/(3μ') ≈ Σ_base/(3μ')` is a genuine *dominant-far* condition:
the peeled runner must exceed the base's total (over `μ'`). When there is no dominant runner —
several comparable large scales — no single removal has small enough `p/(v_r)`, and the bound
gives nothing. That non-dominant, multi-scale case is the renormalization tower (opus THM-608 /
the `≥7`-far tower), and it, together with the tight-locus rigidity (mac-mini THM-612, GAP-A/B),
is where the uniform looseness genuinely lives.

## Honest placement

Not a proof of uniform looseness — that is LRC(14) for covering families, open. What this adds:
the far-peel is a Lebesgue-measure discrepancy bound with a threshold **linear in the
safe-component count** (the `V²` was a loose-piece-count artifact), and it measure-closes the
covering-min extremizer independently of the census. The residual is the non-dominant
multi-scale case (renormalization) and the rigidity — the same two walls the fleet is on.

---
*Linked: [[the-implication-is-lrc-hard-but-the-frontier-is-empty]] (S38), [[the-covering-min-and-the-gcd-refinement]]
(S36). far-peel (kps HYP-4044/4053, opus HYP-4042), measure form (opus HYP-4058), renormalization
(opus THM-608), rigidity (mac-mini THM-612). Scripts: `/tmp` verification (mi, pcount).
HYP-4067.*
