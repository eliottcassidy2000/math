---
id: HYP-2635
status: PARTIAL — p_6 maximality VERIFIED (closed form); HYP-2607 crux consolidated, tools mapped, leads identified; LRC(14) NOT proved
source: kind-pasteur-2026-06-19-S11
related:
  - HYP-2607   # the crux: consec maximizes L_y = E[g(N)]
  - THM-534    # the moment-LP dual (per-E bound)
  - THM-531    # AP-orbit invariance
  - HYP-2606   # relation-lattice covolume (consec minimizes it)
  - HYP-2610   # stranger-contraction (peeling)
  - HYP-2604   # consec maximizes meas(S7)
---

# HYP-2635: the HYP-2607 frontier consolidated — the p_6 lemma, the cleaner target p_0, and what is eliminated

A kind-pasteur-S11 session + a 4-angle adversarial workflow (coupling / wide-bound / literature /
three-distance, all verify-confirmed) sharpened the LRC(14) endgame. Summary for the next agent.

## 1. A clean PROVED-style sub-lemma: consec maximizes the all-missed term `p_6`

`p_6(E) := meas{x : frac(e x) ∈ [0,1/7) ∀ e∈E}` (all six sectors {1..6} missed; the `N=6` term,
weight 1 in `L_y` for k=8). VERIFIED exact (k=6,7,8; gcd(E)=1):
```
  p_6(consec_k) = 1/(7(k-1))   (k=6,7,8 -> 1/35, 1/42, 1/49),
  and consec UNIQUELY maximizes p_6 (0 sets tie or beat).
```
Reason (proof sketch): the all-missed event near x=0 is the interval `[0, 1/(7·max E))`, length
`1/(7·max E)`; consec MINIMIZES `max E = k-1` among k distinct nonneg integers with 0, giving the
largest near-0 interval `1/(7(k-1))`; the wrapping contributions never lift any E above this.
Dilation-invariant (`p_6(cE)=p_6(E)`), so equality holds for all dilations of consec.
**This is a self-contained extremality where consec's "smallest max element" is decisive** — a model
for the kind of argument the full crux needs. (Script `04-computation/lrc14_p6_maximality_kps.py`.)

## 2. The cleaner target: consec maximizes `p_0 = meas(S7)` ITSELF (not just `L_y`)

VERIFIED exact (angle A + its adversarial verify, k=8 and k=9, bounded spread): **consec STRICTLY
maximizes `p_0 = meas(S7)`** (0 beat, 0 non-dilation ties). Exact: `p_0(consec_8)=481/1470`,
`p_0(consec_9)=2447/5880`; and `p_0(consec_k) < cap_k` directly. So the cleanest sufficient lemma is
HYP-2604 ("consec maximizes meas(S7)"), STRONGER and cleaner than HYP-2607's `L_y`; it bypasses the
moment dual entirely. The mechanism (verified): competitors LEAK `p_0` mass into `p_1` (where `g=0`).
Full exact consec distributions (angle D, verified):
```
  consec_8: p = [481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49]
  consec_9: p = [2447/5880, 653/2940, 27/196, 23/245, 13/196, 9/196, 1/56]
```

## 3. The correlation-boost frame (corrected independent limit)

`L_y` measures POSITIVE CORRELATION among orbit points. The fully-independent limit, with the
CORRECT exponent `k-1` (e=0 is pinned in sector 0, only the other `k-1` generators are free):
```
  L_y^inf = Σ_r y_r C(6,r)(1-r/7)^{k-1}:  k=8: 40573/823543≈0.0493;  k=9≈0.1513;  k=10≈0.2188.
```
All `<< cap` (margins ≥ 0.33). consec (the AP, MAX correlation via the densest relation lattice,
HYP-2606) sits at ≈0.358 (k=8). So the danger is ENTIRELY near-AP; wide/dissociated → independent
limit (huge margin). The crux is exactly "the AP maximizes the (signed) correlation `L_y`."

## 4. What is ELIMINATED (do not retry — all verify-confirmed)

- **Coupling / convex-order / FOSD / stochastic-dominance on `N_E`**: ALL fail. Reasons: `g` is not
  convex (k=8: `g=(1,0,0,1/10,0,0,1)`, `Δ²g` sign `+,+,−,+,+`); the moment weights `y_r` ALTERNATE in
  sign (no separable/monotone-order argument); and `L_y` is DILATION-INVARIANT (no stochastic-dominance
  argument can distinguish a set from its dilate). This closes the "AP-orbit-majorization" hope.
- **Dissociation/spread DICHOTOMY (the HYP-2610 peel as stated)**: FAILS. Many primitive wide-spread
  (>20) k=8 sets have NO dissociated stranger — every element sits in a height-2 relation (3873/20000
  sampled; e.g. `[0,3,5,16,28,30,33,35]`, `[0,4,12,15,20,21,25,31]`). So "peel a dissociated stranger
  until bounded spread" does not cover all sets.
- **Literature sieve**: Sungkawichai–Trakulthongchai (arXiv:2604.23906) proves only k≤12 (≤13 runners);
  LRC(14)=k=13 is NOT covered (genuinely open, not preempted). Their method is a DIVISIBILITY-SIEVE
  barrier with no measure content — does not transfer to bounding `L_y`/`meas(S7)`.
- (Earlier: PD/Gram/Fredholm kernel; Selberg minorant; absolute |corr|≤C·W; local monotone descent.)

## 5. The LIVE leads (what MIGHT work)

- **The additive-energy reframe of the wide gap.** "No dissociated stranger ⟹ every element in a
  height-2 relation ⟹ HIGH ADDITIVE ENERGY ⟹ (Freiman) E ⊆ a low-dimensional generalized AP (GAP)."
  GAPs are multi-dim dilated APs; AP-orbit-invariance (THM-531) + dimension should bound `L_y` on them.
  This would replace the failed dichotomy with: {dissociated stranger → peel} ∪ {high energy → GAP,
  bounded by AP-invariance}. NEW concrete target — needs a quantitative Freiman + a GAP `L_y` bound.
- **Tao's second-moment / cluster method** (the only literature tool that structurally matches): bound
  the variance of the number of hit sectors `H(x)=7−N(x)` and conclude `P(H=7)=p_0` is maximized by the
  AP. Worth importing in detail.

## 6. The synthesis verdict — the precise remaining gap (the "third pocket")

The workflow synthesis (verify-confirmed) pins it:
- **Angle B is closest** — the only angle with a CLOSURE ARCHITECTURE: Part 1 (independent limit
  `L_y^inf` exact, ≪ cap) + Part 2 (bounded-spread finite check, consec unique max < cap). What is
  missing is one explicit step bridging them.
- **The gap is the THIRD POCKET:** primitive, wide-spread, **no dissociated stranger** (rich relation
  lattice Λ(E) yet large diameter — partial-AP / Sidon-complement sets, e.g. `[0,3,5,16,28,30,33,35]`).
  These are covered only EMPIRICALLY (they cap at ≈0.18/0.28/0.38 for k=8/9/10, margin ≥0.20). To close:
  bound `Σ_{m≠0,m∈Λ(E)} |ĝ(m)|` (the signed correction `L_y(E)−L_y^inf`) for sets whose Λ(E) is rich.
  The dissociation-peel (HYP-2610) cannot reach them; the absolute envelope DIVERGES (MISTAKE-078); so
  the bound MUST be SIGNED.
- **The binding constraint is the k=9 cap-margin = 0.001384** (consec is only that far below cap_9; the
  nearest competitor `{0,…,7,9}` is 0.005587 below consec). The TIGHTNESS lives entirely at consec
  (bounded spread, in the DONE finite check); the third pocket has margin ≥0.20. So the rigorous gap is
  a LOOSE bound on a generous-margin region — the obstruction is purely the SIGNED cancellation over a
  rich lattice, not tightness.

## 7. Route convergence — this meets codex's reciprocal-tail route

The third pocket (signed correction over a rich relation lattice Λ(E)) is EXACTLY codex's open target
(HYP-2632/2633/2634): "the integer relation-lattice lift samples the finite lanes fairly, then signed
Abel/summation-by-parts inside additive-frequency shells before taking absolute values." My empty-sector
`p_0`/`L_y` route and codex's reciprocal-hyperplane route are two coordinates on the SAME final object:
the signed sum `Σ_{m∈Λ(E)} ĝ(m)` over the relation lattice of a wide, high-additive-energy set. Both
need: a residue-lift/equidistribution statement making the signed shell-sums controllable. **This is the
single shared crux of both active LRC(14) threads.**

## Honest status
LRC(14) NOT proved. The crux is the single inequality "consec maximizes `meas(S7)`" (cleaner than
`L_y`); all majorization/coupling/dichotomy tools are dead; the live routes are the additive-energy/GAP
structure of the wide non-AP sets and Tao's second-moment. Scripts: `lrc14_empty_sector_distribution_kps.py`,
`lrc14_p6_maximality_kps.py`, `lrc14_p0_maximality_test_kps.py` (+ workflow verify scripts).
