# Making the gap-growth rigorous — and finding it is false: the confinement gap is uniform, not a magnitude bound

*opus-2026-07-03-S63. The owner asked me to make the S62 "gap grows with u_max" rigorous — which would
bound the even part and close the m=2, f=2 confinement. Attempting the proof, I first stress-tested the
premise on the commensurate families random sampling misses, and it collapsed: the gap does not grow. This
is a correction of my own S62 (MISTAKE-101), caught before anyone built on it. Honest, per the discipline.*

## The premise was false

S62 claimed the f=2 tightening gap `gap(U) = min_{odd w_1,w_2}(M(2U ∪ {w_1,w_2}) − 1/14)` rises with
`u_max`, over *random* 11-runner even parts. Random subsets are incommensurate, and those do have growing
gap. But the near-extremal (small-gap) families are **commensurate dilated APs**, which random sampling
never hits. Exact arithmetic on them:

> `U = c·{1..11}`:  `M(2c·{1..11} ∪ {w_1,w_2}) = 1/12` exactly,  so  `gap = 1/12 − 1/14 = 1/84 ≈ 0.0119`,
> **constant** for `c = 1,2,3,4,5` (`u_max = 22c → ∞`), every family **primitive**.

The gap is flat at `1/84` at every scale. The reason: `M(2c·{1..11}) = M({1..11}) = 1/12` (scale-invariant
AP), and two odd tighteners **cannot reduce it** — they fail to cover the point achieving `1/12`, so
`M(S) = 1/12` regardless of `u_max`. The tighteners are useless on these families; the "gap" is just how
loose the AP even part already is.

So **`u_max` is not bounded by gap-growth** — dilated APs are near-extremal at arbitrarily large magnitude.
S62's "the finite per-U check terminates" is wrong and retracted (MISTAKE-101). It is the same failure mode
as MISTAKE-095/096/098: a random/generic sampler under-represents the maximally-commensurate family and
shows a false trend. I should have tested the dilated AP before claiming a magnitude bound.

## The correct picture

The min-gap extremizers are **dilated APs at all scales** — magnitude-unbounded. So the confinement fact,
if it holds, is a **uniform positive gap**

> `inf_U gap(U) > 0`  (min `≈ 1/84`, extremized by `c·{1..11}` for all `c`),

not a bound on `u_max`. "Bound `v_max(U)`" is the **wrong target**: there is nothing to bound — the
extremizers run off to infinity with the gap pinned at `1/84`. mac-mini's finite-per-`U` check (Lemma D)
is finite *for each* `U`, but the set of near-extremal `U` is infinite (all dilations), so it does **not**
become globally finite by bounding magnitude. Confinement must be proved as a uniform gap, scale-invariantly
— which is the same shape as everything else in this problem (measure floor, covering-min): a
*scale-invariant positive infimum*, extremized by an AP, whose proof is the tight-locus rigidity.

This actually *simplifies the map*, honestly: it says the `m=2, f=2` residual is **not** a separable
"bound the magnitude then check finitely" problem (as S62 and the Lemma-D framing hoped), but a genuine
uniform-gap / rigidity statement — consistent with mac-mini's own S34 reframing ("confinement is a
covering-min piece, not an independent gap"). The dilated APs are the reason: they are the AP again, dilated,
exactly the imprimitive/AP locus the whole rigidity is about.

## Status

- **Refuted (exact):** S62/HYP-4068 "gap grows with `u_max` ⟹ `u_max` bounded" — `gap(c·{1..11}) = 1/84`
  constant, primitive, `u_max → ∞`. MISTAKE-101 logged.
- **Corrected picture:** confinement (m=2, f=2) is a uniform positive gap `inf_U gap(U) > 0` (≈ 1/84),
  extremized by unbounded dilated APs; "bound `u_max`" is not the path.
- **Stands:** the per-family reduction `w_i ≤ 12 u_max` (S61/THM-614) — a correct bound on the tighteners
  in terms of `u_max`; it simply does not bound `u_max`.
- **Open (unchanged):** confinement itself, now correctly framed as a scale-invariant uniform gap = the
  tight-locus AP rigidity (HYP-4062, THM-612).

I flag plainly: this session's deliverable is a *negative* result — the requested gap-growth does not exist
— caught before it propagated. The value is redirecting the fleet off "bound `u_max`" and onto the uniform
gap, and the exact dilated-AP extremizer that shows why.

Related: MISTAKE-101 (the artifact), HYP-4068 (refuted), THM-612 Lemma D + mac-mini S34 (the finite-per-U
framing this corrects, and the covering-min reframing it agrees with), HYP-4062/kps (the AP rigidity the
uniform gap reduces to), THM-614/HYP-4066 (the `w_i ≤ 12 u_max` reduction that stands). Script:
`lrc14_gap_growth_refuted_dilated_AP_opus_S63.py`. No HYP reservation (a correction).
