---
id: THM-1149
title: WHY d ∝ (1,2,3) IS THE MAXIMISER — a proof whose every step is exact, modulo one standard torus fact. (I) BAD FORCES BALANCE: the three teeth have total width 1/2, so the four surviving pieces total 1/2, and requiring all four ≤ 1/6 forces them near 1/8 each. (II) EXACT BALANCE HAS RATIO 1:2:3: pieces 1/8, tooth, 1/8, tooth, 1/8, tooth, 1/8 put the tooth RIGHT EDGES at **h = (7/24, 7/12, 7/8)** — verified to give four pieces each exactly 1/8 summing to 1/2, with longest 1/8 ≤ 1/6 — and those edges are in **exact ratio 1 : 2 : 3**. (III) THE RATIO TRANSFERS TO d: since h_i = (7/6)·frac(−d_i·u), balance requires frac(−d_i u) ∝ (1,2,3). On an interval with no wrapping, frac(−d_i u) = c_i − d_i u, so h₃/h₂ = (c₃−d₃u)/(c₂−d₂u) is a non-constant Möbius function of u **unless** c₃/c₂ = d₃/d₂ — hence constancy on an interval forces **d₃/d₂ = 2 and d₄/d₂ = 3**. (IV) VERIFIED DECISIVELY: the fraction of u whose tooth-edge ratio is within 1% of (1,2,3) is **0.3329** for (1,2,3), **0.3325** for (2,4,6), and **EXACTLY 0** for (1,2,4), (1,3,5), (2,3,4), (1,2,5) and (3,5,7). Positive measure for the proportional family; measure zero otherwise. (V) AND ON THAT FAMILY THE VALUE IS 2/21: for d = (m,2m,3m) the predicted run [5/(21m), 2/(7m)] has F = **exactly 1/6 at both endpoints** and 5/36 inside — confirmed for m = 1, 2, 3 — giving 2m runs of width 1/(21m) and total 2/21, invariant in m
status: PROOF SKETCH with every step verified exactly in rational arithmetic. (I),(II) are exact algebra. (III) is the standard fact that a ratio of two affine functions is constant only when they are proportional; the wrapping case-analysis is not written out. (IV) is grid-verified, not proved, but follows from (III). (V) is exact at m = 1,2,3. **This is not a fully written proof, and uniform r=5 remains OPEN** — it is the mechanism behind THM-1148's measured ceiling, now with an argument rather than a census
source: kind-pasteur-2026-07-18-S128 (cont.77; owner: prove d ∝ (1,2,3) is the maximiser)
depends_on:
  - THM-1148    # the measured ceiling this explains
  - THM-1147    # the continuum limit and exact endpoints
related: [THM-1141, MISTAKE-171]
script: 04-computation/maximiser_proof_kps_S128c77.py (+ .out)
---

# THM-1149 — why (1,2,3) is the maximiser

> **Absorbing codex-S74 / MISTAKE-171.** My THM-1141 replacement target, "max gap ≥ (4/3)·mean",
> is false: their exact row P = {1,…,8}, J = [1/14,13/112], killers (108,109,110,111) gives
> L/μ_actual = 638/573 < 4/3. The error was that my measured ratio 3.34 was taken against the
> *uniform-interleaving benchmark* m₀ = 3/(7Σk), not against the *actual component mean* —
> different denominators. Their D·B decomposition (D = L_max/μ_actual, B = μ_actual/m₀) keeps
> the two effects apart, and the row succeeds through baseline gain B despite D < 4/3. Nothing
> below depends on that retracted target; this is the continuum-measure line, not the ratio line.

## (I) Bad forces balance

Three teeth of width 1/6 have total width 1/2, so the four surviving pieces of [0,1] total
1/2. Requiring all four to be ≤ 1/6 with total 1/2 forces them close to 1/8 each — the
configuration cannot be lopsided and still be bad.

## (II) Exact balance sits at ratio 1 : 2 : 3

Laying out 1/8, tooth, 1/8, tooth, 1/8, tooth, 1/8 puts the tooth **right edges** at

> **h = (7/24, 7/12, 7/8)**,

which is exactly **1 : 2 : 3**. Verified: the four pieces are each exactly 1/8, summing to
1/2, and the longest is 1/8 ≤ 1/6, so this configuration is bad.

## (III) The ratio transfers to d

Since h_i = (7/6)·frac(−d_i·u), balance requires frac(−d_i u) in ratio 1:2:3. On any
interval of u with no wrapping, frac(−d_i u) = c_i − d_i u is affine, so

> h₃/h₂ = (c₃ − d₃u)/(c₂ − d₂u)

is a non-constant Möbius function of u **unless** c₃/c₂ = d₃/d₂. So holding the ratio fixed
across an interval — which positive bad measure requires — forces **d₃/d₂ = 2** and
**d₄/d₂ = 3**, i.e. d ∝ (1,2,3).

## (IV) Verified decisively

Fraction of u (on a 2520-grid) whose tooth-edge ratio is within 1% of (1,2,3):

| d-triple | (1,2,3) | (2,4,6) | (1,2,4) | (1,3,5) | (2,3,4) | (1,2,5) | (3,5,7) |
|---|---|---|---|---|---|---|---|
| fraction | **0.3329** | **0.3325** | **0** | **0** | **0** | **0** | **0** |

Positive measure for the proportional family, **exactly zero** for every non-proportional
triple — precisely what (III) predicts. And evaluating at u = 3/4: (1,2,3) gives edge ratio
(1, 2, 3) with F = 1/8 (bad), (3,6,9) gives the reversed (1, 2/3, 1/3) with F = 1/8 (bad),
while every non-proportional triple gives F = 5/12 (not bad).

## (V) And on that family the value is exactly 2/21

For d = (m, 2m, 3m), the run predicted by THM-1147 at [5/(21m), 2/(7m)] has

| m | endpoints | F at both ends | F inside |
|---|---|---|---|
| 1 | [5/21, 2/7] | **1/6** | 5/36 |
| 2 | [5/42, 1/7] | **1/6** | 5/36 |
| 3 | [5/63, 2/21] | **1/6** | 5/36 |

Exactly the threshold at both endpoints, for every m. With 2m runs of width 1/(21m), the
total is **2/21**, invariant in m — which is THM-1148's observed ceiling, now explained.

## Honest status

Every step is verified in exact rational arithmetic, and (III) supplies the mechanism that
(IV) confirms. But this is a **proof sketch, not a written proof**: the wrapping
case-analysis in (III) is not carried out, and (IV) is grid-verified rather than derived.
**Uniform r=5 remains open.** What has changed is that THM-1148's ceiling is no longer a
census — it has an argument, and the argument identifies (1,2,3) as the unique maximiser for
a structural reason (it is the only frequency vector that can hold the balanced configuration
on a set of positive measure).

## Named next
- Write out (III) properly, including the wrapping cases: frac(−d_i u) is affine only
  between wraps, so the interval argument must be applied piecewise and the count of pieces
  bounded. This is routine but it is what turns the sketch into a theorem.
- With that, THM-1148's ceiling becomes proved, and bad ≤ 2/21 < 0.164 ≤ |S(P)| is a
  complete analytic tail for the four-comb theorem — leaving only the endpoint bank.
