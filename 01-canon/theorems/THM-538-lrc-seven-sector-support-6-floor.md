---
id: THM-538
title: "The support-6 floor for the LRC(14) seven-sector signed Fourier sum — in meas(S7(E)) = M7(k) + Σ_{0≠n∈Λ°(E)} K(n), the kernel K(n) VANISHES unless the relation n has ≥6 coordinates that are nonzero and ≢0 mod 7; the shortest (support ≤5) relations contribute EXACTLY zero, explaining the F3 'absolute 5× lossy' obstruction (the signed sum annihilates all support-≤5 mass)"
status: RESOLVED/CORRECTED (CASE-thm538-support6-floor-zero-padding CONCEDED by author kps-2026-06-19; reconfirmed independently: K(1,1,−1,0,0,0,0)=+0.00066≠0). The support-6 floor is **PROVED for the active-coordinate sum Q(n)** but **FALSE for the zero-padded kernel K(n)** that appears in the measure (the zero coords contribute |T|-dependent factors (1−|T|/7)^z that break the cancellation). The correct kernel structure is the coset factorization K(n)=D7(n mod 7)/∏n_j (HYP-2646), with the correction conditionally convergent and ruled by support-6 relation DENSITY R6 (HYP-2645), NOT a floor. Short relations (support 2–5) DO contribute (the AP's are support-3-dominated). The stranger-contraction / far-element plateau (HYP-2610/2644) is the live wide-spread route and is UNAFFECTED.
source: kind-pasteur-2026-06-19-S9 (workflow angle wide-spread-signed-weyl)
depends_on:
  - HYP-2606  # the signed offset-relation-lattice Fourier identity meas(S7)=M7(k)+Σ K(n)
  - THM-534   # the seven-sector moment / S7 object
  - THM-503   # the seven-vanishing ĉ_T(7m)=0 (apex prime 7)
related:
  - HYP-2610  # the stranger-contraction / multiplicative decoupling (the wide-spread reduction this feeds)
  - HYP-2608  # the wide-spread-bound route
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case, 13 speeds); inclusion–exclusion; Poisson summation on the relation lattice.
---

# THM-538 — The support-6 floor of the seven-sector signed sum

## Setup (HYP-2606)
For an integer co-offset set `E` (0∈E, |E|=k), the seven-sector cover measure has the exact signed
Fourier expansion over the offset relation lattice `Λ°(E) = {n∈ℤ^{k−1}: Σ_j n_j e_j = 0}` (rank k−2):
> `meas(S7(E)) = M7(k) + Σ_{0≠n∈Λ°(E)} K(n)`,  `K(n) = Σ_{T⊆{1,…,6}} (−1)^{|T|} ∏_j ĉ_T(n_j)`,
where `ĉ_T` = Fourier coefficient of `1_{[0,1)∖∪_{j∈T}[j/7,(j+1)/7)}`, `ĉ_T(0)=1−|T|/7`, and (THM-503)
`ĉ_T(7m)=0`. `M7(k)=K(0)=Σ_{t=0}^{6}(−1)^t C(6,t)(1−t/7)^{k−1}` is the iid (large-spread) limit.

## The theorem — CORRECTED (kps-2026-06-19, resolving CASE-thm538): floor holds for Q, NOT for K

> ⚠️ **The originally-stated support-6 floor "K(n)=0 for support<6" is FALSE for the zero-padded kernel
> K(n) that appears in the measure expansion** (CASE-thm538-support6-floor-zero-padding, CONCEDED by the
> author kps-S9 after independent reconfirmation: K(1,1,−1,0,0,0,0)=+0.00066≠0 at k=8, a support-3
> relation; the AP's 1+2−3=0 relation genuinely contributes). The bug: the proof's `C(U)=0` step dropped
> the **zero-coordinate factors `ĉ_T(0)=1−|T|/7`, which depend on |T|** and weight the alternating T-sum
> by `(1−|T|/7)^z` (z=#zeros), breaking the `(1−1)^{6−|U|}=0` cancellation.

> **CORRECT statement (PROVED).** The floor holds for the **ACTIVE-COORDINATE sum** `Q(n) :=
> Σ_T(−1)^{|T|}∏_{j∈U}ŝ_T(n_j)` (the |U|-dimensional kernel with NO zero padding): `Q(n)=0` for
> support |U|<6 (the proof below is valid for Q), nonzero at |U|=6. The exhaustive "support≤5 ⟹
> |·|=5e-17" computed Q, not K. **THM-538 conflated Q (active) and K (zero-padded).**

**Proof (for Q).** `Q(n)=Σ_T(−1)^{|T|}∏_{j∈U}ŝ_T(n_j)`; for a fixed sector-assignment the T-sum factors
through `C(U)=Σ_{T⊇U}(−1)^{|T|}=(−1)^{|U|}(1−1)^{6−|U|}=0` for |U|<6. ∎ (For Q the zero coords are simply
absent, so no `(1−|T|/7)^z` weight — which is exactly why the argument works for Q but NOT for K.)

## The CORRECT kernel structure (HYP-2646, supersedes the floor framing)
The right object is the **coset/reciprocal factorization** (kps-S12, verified 1.6e-19):
`K(n) = D7(n mod 7)/∏_j n_j`, `D7` a FINITE mod-7 character coefficient, NONZERO on all 46656 cosets.
So `corr(E)=Σ_{c∈(F_7^*)^6} D7(c)·S_c(E)`, `S_c=Σ_{n≡c}1/∏n_j` — the cancellation is the alternating sign
of `Re D7(c)` (|Re D7|≤0.1431) + the CONDITIONAL convergence of the reciprocal sums, NOT support-6
sparsity. Short relations (support 2–5) DO contribute. The wide-spread ruler is the support-6 relation
**DENSITY R6** (HYP-2645: R6=982/924/546/414/156 tracks corr 0.303/0.184/0.213/0.0096/0.0005), and the
convergent representation is the finite x-cell integral (HYP-2645) / far-element plateau (HYP-2644), NOT
the box-truncated lattice sum (which is only conditionally convergent — MISTAKE-078).

## Why it matters

- **It explains the F3 "absolute 5× lossy" obstruction (HYP-2606).** The absolute bound `Σ|K(n)|`
  counts the support-2..5 mass that the SIGNED sum annihilates exactly (52% of the absolute mass at the
  AP, k=8). So the signed structure is not a luxury — the cancellation is total below support 6, and the
  surviving sum is intrinsically a ≥6-body object.
- **It removes all short relations**, so the wide-spread correction is controlled by the **6-fold
  additive energy** `R6(E,L)=#{support-≥6 relations, |n|∞≤L}` (corr is monotone in R6; dissociated/Sidon
  sets have R6 small and corr≈0). This de-risks the Minkowski/successive-minima count.
- **Companion envelope (LEMMA B, PROVED):** `|ŝ_T(n)| = |sin(πn/7)|/(π|n|)` (independent of T, =0 iff
  7|n); sharp per-coordinate `|ĉ_T(n)| ≤ 0.6973/|n|` (2.7× better than the triangle bound 1.862).

## Scope — what it does and does NOT close

THM-538 does **not** close the cap bound by itself: the absolute support-≥6 residual is still ~1.2 ≫
`δ_8 = meas(S7(consec_8)) − M7(8) = 0.303` at the AP. **BUT this is not the operative gap:** the AP (and
all bounded-spread shapes) are handled by the **exact finite check** (consec is the verified argmax,
0 violations), not by this Fourier bound. The Fourier sum is needed only for **WIDE** shapes, where the
stranger-contraction (HYP-2610 / kps-S9: peel the largest offset, multiply by `1−r/7`; fast convergence
`N≥61`) makes the correction collapse to the bounded core's value (wide ceiling ≈0.21 ≪ cap 0.38, margin
0.16). So THM-538 is the structural backbone of the **wide-spread** half (with the contraction), while the
**tight** margin (0.023 at k=8) lives entirely in the finite check. The remaining analytic content is the
explicit single-offset Weyl error `|Ŵ(N)| ≤ C(E')/N` (1-dimensional) + the bounded-spread finite max.
→ HYP-2606, HYP-2610, HYP-2608, THM-534, OPEN-Q-108.
