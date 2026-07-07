---
id: THM-641
title: The WEIGHTED ANCHORED PAIR-MASS LAW ("weighted cherry") — exact four-Bernoulli-term closed form for meas{frac(q1 x) ∈ (a, a+α], frac(q2 x) ∈ (b, b+β]} at arbitrary widths α,β and anchors a,b; THM-638 (both signs) is the (θ,θ) specialization; apex-invisibility (s|cq ⟹ zero correction) holds at ALL anchors and weights; cross-anchor corrections can be NEGATIVE — the mechanism behind the {0,1/2} 2-anchor tail
status: PROVED (elementary: CRT offset sweep + Bernoulli periodization + multiplication theorem; half page below). VERIFIED exact: 156/156 random (q≤24 coprime × widths {1/14,1/7,3/14,2/7} × anchors {0,1/2,1/3}) closed-vs-brute; 489/489 reproduces THM-638 same-sign at θ=1/7; 277/277 mixed-sign closed-vs-brute AND brute-vs-THM-638-mixed. Convention pinned by sign-forcing probes (an earlier φ-sign ambiguity was caught by exactly this verification loop).
source: mac-mini-2026-07-07-S45 (HYP-4947; owner items "2-anchor joint tail" + "closed-form weighted cherry theorem")
depends_on:
  - THM-638   # klein-S156: the equal-width anchored-at-0 law (both signs) — now a corollary
related:
  - HYP-4791/HYP-4831  # klein: Hunter floor, cherry/moment-LP, the pairwise barrier
  - THM-637   # opus: Farey roof; apex-7 invisibility (generalized here)
  - HYP-4947  # this session
external: Bernoulli polynomial multiplication theorem (Raabe); CRT/Bezout.
---

# THM-641 — The weighted anchored pair-mass law

## Statement

Let `q1, q2 ≥ 1` be coprime integers, `N = q1 q2`, widths `α, β ∈ (0,1]` and anchors
`a, b ∈ [0,1)` (all rational for exactness, real in general). Let
`B̄₂(t) = frac(t)² − frac(t) + 1/6` (the periodized second Bernoulli polynomial) and
`φ = q1·b − q2·a`. Then

> `M(q1,q2; α,β; a,b) := meas{x ∈ [0,1) : frac(q1 x) ∈ (a, a+α], frac(q2 x) ∈ (b, b+β]}`
> `= α·β + [ −B̄₂(φ + q1 β) + B̄₂(φ) + B̄₂(φ + q1 β − q2 α) − B̄₂(φ − q2 α) ] / (2 q1 q2).`

For `gcd(q1,q2) = g > 1`, substitute `y = g x mod 1` (measure-preserving):
`M(q1,q2;…) = M(q1/g, q2/g;…)` with the same windows.

## Proof

**1. Pullback and offset sweep.** `frac(q1 x) ∈ (a, a+α]` iff `x` lies in one of the `q1`
intervals `I_{m1} = ((m1+a)/q1, (m1+a+α)/q1]`, each of length `A = α/q1`; similarly `q2`
intervals `J_{m2}` of length `B = β/q2`. So `M = Σ_{m1,m2} len(I_{m1} ∩ J_{m2})`. The
relative offset of the pair `(m1,m2)` is
`δ = (m2+b)/q2 − (m1+a)/q1 = ((q1 m2 − q2 m1) + φ)/N`,
and since `gcd(q1,q2)=1`, the integer `q1 m2 − q2 m1` sweeps every residue mod `N`
exactly once (Bezout/CRT). Hence, with `ov(δ) = len((0,A] ∩ (δ, δ+B])` interpreted on the
circle,
> `M = Σ_{j=0}^{N−1} ov((j + φ)/N)`.

**2. The overlap trapezoid and its Bernoulli form.** `ov` is the 1-periodic trapezoid:
slope jumps `+1, −1, −1, +1` at the corners `c ∈ {−B, 0, A−B, A}` (this corner/jump
multiset is the same whether `A ≥ B` or `A < B`), with mean `∫₀¹ ov = A·B`. On the
circle, `(B̄₂/2)'' = 1 − δ_ℤ` (as distributions: `B̄₂' = 2t−1` jumps by `−2` at integers).
A 1-periodic piecewise-linear function is determined by its slope jumps and mean, so
> `ov(δ) = A·B − Σ_c (Δs_c / 2) · B̄₂(δ − c)`
(both sides periodic, equal means since `B̄₂` has zero mean, equal second derivatives
since `Σ_c Δs_c = 0` kills the constant).

**3. Collapse by the multiplication theorem.** Raabe's identity
`Σ_{j mod N} B̄₂((t+j)/N) = B̄₂(t)/N` gives
`Σ_j B̄₂((j+φ)/N − c) = B̄₂(φ − N c)/N`. Summing step 2 over the offsets of step 1:
`M = N·A·B − (1/2N) Σ_c Δs_c B̄₂(φ − N c)`. Now `N·A·B = αβ`, and the corner values
scale to `N·(−B) = −q1β`, `0`, `N(A−B) = q2α − q1β`, `N·A = q2α`; with jumps
`(+1,−1,−1,+1)` and `B̄₂` even this is exactly the displayed formula. ∎

## Corollaries

**(C1) THM-638, both signs.** At `(α,β,a,b) = (θ,θ,0,0)` with `θ = c/s`: `φ = 0` and the
formula evaluates to `θ² + min(r1,r2)(s − max(r1,r2))/(s² q1 q2)`, `r_i = c q_i mod s`
(0 ↦ s) — klein's same-sign law (489/489 exact). The mixed-sign law is the anchoring
`a = −θ` (window `(−θ, 0]`): reproduces `θ² − min(r1 r2, (s−r1)(s−r2))/(s² q1 q2)`
(277/277 exact). THM-638 is the two-point specialization of a four-parameter law.

**(C2) Apex-invisibility at ALL anchors and weights.** If `q2 α ∈ ℤ` (e.g. `α = c/s` and
`s | c q2`), then `B̄₂(φ + q1β − q2α) = B̄₂(φ + q1β)` and `B̄₂(φ − q2α) = B̄₂(φ)`, so the
correction is IDENTICALLY ZERO — for every width `β` and both anchors. This generalizes
the THM-637/THM-638 apex-7 vanishing rows (7 | q1q2 at θ = 1/7) from the equal-width
anchored-at-0 case to the full four-parameter family: apex speeds are invisible to pair
statistics no matter how the windows are weighted or anchored.

**(C3) Cross-anchor anticorrelation — the {0,1/2} mechanism.** At `θ = 1/7` the exact
`(a,b) ∈ {0,1/2}²` table (correction × 2N shown; script section iv) has NEGATIVE entries
exactly at the even/odd interplays, e.g. `(q1,q2) = (1,2)`: `(0,0) ↦ +10/49` but
`(0,1/2) ↦ −4/49`; `(3,4)`: `(0,1/2) ↦ −17/49`. Since the 2-anchor tail lemma consumes
`P(A_0 ∪ A_{1/2}) = p_0 + p_{1/2} − P(A_0 ∩ A_{1/2})` and the joint term's pair layer is
built from these cross-anchor masses, **anticorrelated pairs make the union LARGER than
independent** — the quantitative reason the `{0, 1/2}` anchor pair is the robust choice
(boxeph's 0.187) and the exact input layer for proving `inf_E P(max(gap@0, gap@1/2) >
1/7) ≥ T_k` (the fleet's load-bearing 2-anchor lemma). The law supplies every pair mass
of that Bonferroni/Hunter computation in closed form, for any window-weight optimization.

**(C4) Full pair-marginal data.** Sweeping `(α,β)` gives the COMPLETE joint law of any
two phases `(frac(q1 x), frac(q2 x))` in closed form — the maximal "pairwise
information". NB (honest): a first attempt to feed width-swept data into a
pairwise-marginal LP at klein-S159's barrier shape produced a non-comparable relaxation
(my 3-bin event algebra ≠ klein's 128-atom LP; my LP min = 0, meaning my encoding loses
their structural constraints). Whether full pair marginals move the 0.1233 barrier is
OPEN and now a well-posed handoff: rerun klein's exact LP with THM-641's width-swept
constraints added.

## Verification & files
`04-computation/lrc14_weighted_cherry_law_macmini_S45.py` (+ `.out`): convention pinned
by sign-forcing probes; broad exact verification as in the status line. The φ-sign
ambiguity (B̄₂ evenness makes symmetric probes insensitive) is itself a recorded lesson:
reconciliation probe sets must include asymmetric anchors.
