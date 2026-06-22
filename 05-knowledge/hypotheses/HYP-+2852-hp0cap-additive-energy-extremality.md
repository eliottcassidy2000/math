---
id: HYP-+2852
title: hp0cap extremality "consec maximizes L_y" = [APs maximize additive energy A(E) -- classical/rigorous] + [L_y monotone in A(E) -- the bridge]; connects the LRC(14) sector-cover extremality to classical additive combinatorics
status: CONCEPTUAL FRAME ESTABLISHED + VERIFIED (k=8,9,10). Piece (a) rigorous (classical). Piece (b) the bridge = the remaining analytic gap, now cleanly isolated.
source: mac-mini-2026-06-22-S29
related:
  - THM-534   # p0<=L_y; residual = consec maximizes L_y (this DECOMPOSES that residual)
  - HYP-2840  # the L_y route + sharp arc-complexity
  - HYP-2637  # higher Freiman dimension LOWERS p0 (= low additive energy => low L_y, the same frame)
  - THM-531   # AP-orbit invariance (APs tie at the max -- consistent: scale-invariance)
---

# HYP-+2852 -- hp0cap: the additive-energy extremality frame

## The decomposition
hp0cap (via THM-534) = "consec maximizes L_y(E)" where L_y = Σ_r y_r S_r is the moment-LP
dual. This residual DECOMPOSES into two pieces via the **additive energy**
`A(E) = #{(a,b,c,d) ∈ E⁴ : a+b=c+d} = Σ_s r_E(s)²`:

**(a) A(E) ≤ A(consec_k) [APs MAXIMIZE additive energy among k-sets] -- CLASSICAL, rigorous.**
A(consec_k) = additive energy of [k] = Σ_s r(s)² (triangular r(s)=1,2,..,k,..,2,1).
Maximality is the standard rearrangement fact (AP peaks the representation function;
Sidon sets minimize it). VERIFIED exhaustive (max(E)<=k+4): consec is THE argmax,
A(consec_k) = 344, 489, 670 for k=8,9,10.

**(b) L_y(E) is MONOTONE in A(E) [the bridge] -- empirical, the remaining gap.**
VERIFIED: over 900 random configs (k=8,9,10, spreads up to 5k), ZERO have L_y > L_y(consec)
and ZERO have A > A(consec). L_y tracks A near-perfectly: (A,L_y) sorted is co-monotone
(consec 489/0.493 → Sidon 153/0.148 at k=9). consec is the GLOBAL max of L_y, on the
upper envelope of the (A, L_y) scatter, which is increasing in A.

## Why this is the right frame
- **Unifies with HYP-2637** (Freiman dimension penalty): high dimension ⟺ low additive
  energy ⟺ low L_y. AP = dimension 1 = max energy = max L_y. The "dimension penalty" IS
  the additive-energy monotonicity.
- **Fourier mechanism:** L_y(E) = L_y_decorr(k) + Σ_{relations h: Σ h_e e=0} (sector-Fourier
  weights). The DOMINANT relations are additive quadruples (a+b=c+d), contributing ∝ A(E).
  So L_y = decorr + w·A + (higher-order relations); the leading term is the additive energy.
- **Reduces the crux to a CLASSICAL question**: combined with (a) [rigorous], the bridge (b)
  "L_y ↑ A" is the single remaining analytic statement for the bounded extremality.

## Status / next
The bounded-spread "consec maximizes L_y" (my half of the kps division) reduces to the
bridge (b). To close: (i) the exact functional L_y(A) (is L_y affine in A, or which
weighted energy is exact?); (ii) bound the higher-order relation terms. → THM-534, HYP-2840,
HYP-2637. Far-element half is kps (resonance/dilation, S32).

Files: lrc14_hp0cap_additive_energy_macmini_S29.py (to save).
