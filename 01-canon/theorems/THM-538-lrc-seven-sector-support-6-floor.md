---
id: THM-538
title: "The support-6 floor for the LRC(14) seven-sector signed Fourier sum — in meas(S7(E)) = M7(k) + Σ_{0≠n∈Λ°(E)} K(n), the kernel K(n) VANISHES unless the relation n has ≥6 coordinates that are nonzero and ≢0 mod 7; the shortest (support ≤5) relations contribute EXACTLY zero, explaining the F3 'absolute 5× lossy' obstruction (the signed sum annihilates all support-≤5 mass)"
status: PROVED (elementary algebra: C(U)=Σ_{T⊇U}(−1)^{|T|}=(−1)^{|U|}(1−1)^{6−|U|}=0 for |U|<6; the support-2 kernel matrix is identically zero). VERIFIED exhaustively (support ≤5 max |K|=5e-17 incl 7-multiples; nonzero at support 6). Does NOT by itself close the cap bound (absolute support-≥6 residual ~1.2 > δ_8=0.303 at the AP — but the AP is bounded-spread, handled by the exact finite check, NOT this bound).
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

## The theorem (PROVED)

> **Support-6 floor.** `K(n) = 0` unless `n` has **at least 6** coordinates that are both nonzero and
> not multiples of 7. Equivalently: only relations supported on ≥6 distinct offsets (with non-7
> coordinates) contribute; **every support-≤5 relation — in particular the shortest (λ₁) relations —
> contributes EXACTLY zero.**

**Proof.** Write `ĉ_T(n_j) = δ_{n_j,0}·(1−|T|/7) + [n_j≠0]·ŝ_T(n_j)` and expand `K(n)` over which
coordinates are nonzero. Let `U = supp(n)` (the nonzero coordinates; `n_j≡0 mod 7` give `ŝ_T(n_j)=0` by
the 7-vanishing, so effectively drop them too). For a fixed assignment of sectors, the `T`-sum factors
through `C(U) := Σ_{T⊇U}(−1)^{|T|}` over the sectors NOT pinned by the nonzero coordinates. Carrying out
the alternating sum,
`C(U) = Σ_{T⊇U}(−1)^{|T|} = (−1)^{|U|}\sum_{S⊆\{1..6\}∖U}(−1)^{|S|} = (−1)^{|U|}(1−1)^{6−|U|} = 0`
unless `|U| = 6`. Hence a surviving term needs all six sectors pinned, i.e. `≥6` genuine (nonzero,
non-7) coordinates. In particular the support-2 kernel `W[j,l]=Σ_{T⊇\{j,l\}}(−1)^{|T|}` is identically
zero. ∎ (Verified exhaustively: support ≤5 gives max `|K(n)|=5·10⁻¹⁷`, including zeros and 7-multiples;
first nonzero at support 6.)

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
