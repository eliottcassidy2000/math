---
id: HYP-3779
title: THREE MORE REGULARIZATION AVATARS of the LRC14 margin (deepening S67 HYP-3774) -- topological, exact, and p-adic -- with multiple proof-directions. The seed (S67): margin(n)=n/Phi6-1/n=-12 s(n,Phi6)/n^2, s(n,Phi6)=-T/(12T+6)->-1/12=zeta(-1), T=1+..+(n-1). Beyond the Ramanujan/archimedean avatar (S67), THREE more, all verified: (AVATAR 1, SPECTRAL/TOPOLOGICAL) the Dedekind sum has the cotangent form s(h,k)=(1/4k) sum_j cot(pi j/k)cot(pi h j/k), which is the APS eta-invariant / Hirzebruch signature defect of the LENS SPACE L(k,h); so the LRC14 margin s(14,183)=-91/1098 IS the eta-invariant of the 3-manifold L(183,14), and Dedekind-Rademacher reciprocity = the APS cobordism/gluing formula. (AVATAR 2, EULER-MACLAURIN) the margin is EXACTLY the Bernoulli-B2 (2nd-order) remainder of the speed-sum T: margin=-12 s/n^2 and n^2*margin -> -12 zeta(-1)=1, with B2=1/6 the E-M coefficient. (AVATAR 3, p-ADIC) Kubota-Leopoldt zeta_p(-1)=-(1-p)/12; at the apex prime 7, zeta_7(-1)=1/2, which DISAGREES with the archimedean zeta(-1)=-1/12 -- the un-regularizable residual (f14 at the 7-cusp) is exactly where the 7-adic and archimedean regularizations SPLIT. SYNTHESIS: the LRC14 margin has THREE regularization avatars (Ramanujan/archimedean -1/12, spectral/eta-invariant of L(183,14), p-adic zeta_p(-1)); they AGREE on the bulk but SPLIT at the apex prime 7 = the residual = the genus-1 cusp form. PROOF-DIRECTIONS: (1) bound the margin via lens-space eta-invariant positivity + cobordism (reciprocity); (3) characterize the residual as the 7-adic/archimedean regularization discrepancy (the "wild at 7" part).
status: AVATARS EXACT/verified. Avatar 1 (cotangent=sawtooth=exact Dedekind sum; = lens-space eta-invariant, classical Hirzebruch-Zagier/APS) verified n=4,7,14. Avatar 2 (margin = -12 s/n^2, n^2*margin->1) exact. Avatar 3 (zeta_p(-1)=-(1-p)/12, zeta_7(-1)=1/2) exact (Kubota-Leopoldt). The three-avatar synthesis + the two proof-directions are reframings/directions, NOT closed proofs. Deepens S67 HYP-3774 (does not repeat it: those were the Mobius/24/hexagonal/Faulhaber/bulk-residual results).
source: mac-mini-2026-06-30-S70
related:
  - HYP-3774   # S67: the zeta-regularization carrier (Mobius interpolation, 24=psi(14), hexagonal anomaly) -- this DEEPENS it
  - HYP-3768   # S64: margin = order-6 Dedekind sum; klein-S56 proved the closed form
  - HYP-3586   # X0(14) apex cusp d=7, the f14 cusp-form obstruction (the residual)
  - HYP-2808   # Dedekind-Rademacher reciprocity (the far-coherence tool) = APS cobordism here
  - HYP-3771   # the crystallographic (2,3,p) spine; the apex-7 hyperbolic
results:
  - 04-computation/regularization_avatars_eta_padic_macmini_20260630.py
  - 05-knowledge/results/regularization_avatars_eta_padic_macmini_20260630.out
---

# HYP-3779 -- three more regularization avatars of the margin (topological, exact, p-adic)

The owner's regularization seed, revisited to go **deeper** than S67 (HYP-3774). Seed: `margin(n) = n/Phi_6 -
1/n = -12 s(n,Phi_6)/n^2`, with `s(n,Phi_6) = -T/(12T+6) -> -1/12 = zeta(-1)` and `T = 1+...+(n-1)`, `Phi_6 =
2T+1`. S67 gave the Ramanujan/archimedean avatar. Here are **three more**, each verified.

## Avatar 1 -- SPECTRAL / TOPOLOGICAL: the margin is a lens-space eta-invariant
The Dedekind sum has the **cotangent (spectral) form**
`s(h,k) = (1/4k) sum_{j=1}^{k-1} cot(pi j/k) cot(pi h j/k)` (verified `= sawtooth = exact` for `n=4,7,14`).
This cotangent sum is the **APS eta-invariant / Hirzebruch signature defect of the lens space `L(k,h)`**
(a `zeta`-regularized trace of the signature operator; Hirzebruch-Zagier, Atiyah-Patodi-Singer). Hence

> the LRC14 margin `s(14,183) = -91/1098` **is the eta-invariant of the 3-manifold `L(183,14)`**,

and **Dedekind-Rademacher reciprocity** (the repo's far-coherence tool, HYP-2808) **is the APS cobordism/gluing
formula** (the lens spaces co-bound). So the margin is not just an arithmetic regularization -- it is a
**topological/spectral** one. *Proof-direction:* bound the margin via eta-invariant positivity and the
cobordism relation; the cobordism gluing is the natural home for the CRT-linkage the residual needs.

## Avatar 2 -- EULER-MACLAURIN: the margin is the B2-remainder of the speed-sum
`margin = -12 s(n,Phi_6)/n^2` (exact), and `n^2 * margin -> -12 zeta(-1) = 1` (verified `0.977, 0.995, ... ->
1`). The margin is exactly the **Bernoulli-`B_2` (second-order) Euler-Maclaurin remainder** of the finite
speed-sum `T`: `B_2 = 1/6`, `zeta(-1) = -B_2/2 = -1/12`, and the `-1/12` is the E-M coefficient. So "finite
carries the actual sum, asymptotic carries the regularization" (owner) is literally the Euler-Maclaurin
expansion of `sum k`.

## Avatar 3 -- p-ADIC: the apex prime 7 regularizes DIFFERENTLY (the residual)
Kubota-Leopoldt: `zeta_p(-1) = -(1-p)/12`. At the **apex prime 7**: `zeta_7(-1) = -(1-7)/12 = 1/2`, which
**disagrees** with the archimedean `zeta(-1) = -1/12`. The Euler factor `(1-p)` is nontrivial exactly at the
covering prime's own place. So:

> the un-regularizable residual (the genus-1 cusp form `f_14` at the 7-cusp) is exactly **where the 7-adic and
> archimedean regularizations split**.

*Proof-direction:* characterize the residual as the "wild at 7" part -- the bulk is the away-from-7
(archimedean `-1/12`) regularization; the hard core is the discrepancy `zeta_7(-1) - zeta(-1) = 1/2 + 1/12 =
7/12` localized at the apex.

## Synthesis -- three avatars, one split
The LRC14 margin has **three regularization avatars**: the Ramanujan/archimedean `zeta(-1) = -1/12` (S67), the
**spectral eta-invariant** of `L(183,14)` (Avatar 1), and the **p-adic** `zeta_p(-1)` (Avatar 3). They agree on
the **bulk** (all encode the `-1/12`-type second-order behavior away from the apex) and **split at the apex
prime 7** -- which is precisely the residual / the cusp form `f_14`. This upgrades the S67 "regularizable bulk
+ un-regularizable residual" dichotomy to a **concrete triple**: the residual is the place the three
regularizations fail to coincide, and (Avatar 3) that place is the 7-adic one.

## Honest scope
Avatars 1-3 are exact/classical identities (Dedekind = lens-space eta-invariant; Euler-Maclaurin;
Kubota-Leopoldt `zeta_p(-1)`), verified. The three-avatar synthesis and the two proof-directions (eta-invariant
+ cobordism bound; residual = 7-adic/archimedean split) are reframings and directions, **not** closed proofs.
This DEEPENS S67 (HYP-3774) -- which had the Mobius interpolation, `24 = psi(14)`, the hexagonal anomaly, and
the Faulhaber/bulk-residual split -- with three genuinely different regularization lenses (topological, exact,
p-adic), not a repeat.
