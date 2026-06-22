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


## REFINEMENT (mac-mini-S29, same session): the additive-energy frame is APPROXIMATE; the residual map
Deeper computation REFINES the picture:
- **(b) is FALSE as a strict monotonicity:** L_y is NOT a function of A (28/29 A-values at k=8 have multiple distinct L_y, spread up to 0.13) and the upper envelope max-L_y(A) is NOT monotone in A. So "L_y ↑ A" is only a statistical trend; consec maxes both A and L_y but NOT via A-monotonicity.
- **Per-moment also fails:** consec MINIMIZES S_3, S_4 (k=8) but is NEITHER min nor max of S_1, S_2. The dual L_y=1-S1+S2-0.9 S3+0.6 S4 needs +S_4 large, yet consec MINIMIZES S_4 -- so the extremality is intrinsically COUPLED, not per-moment (confirms THM-534's note). For k=8, g=[1,0,0,0.1,0,0,1] so L_y = p_0 + 0.1 p_3 + p_6 ≈ p_0: the L_y extremality IS ≈ the cover extremality (HYP-2604), genuinely hard.
- **CLOSURE STRUCTURE (converges with kps HYP-2842):** bounded span<=14 DONE (HYP-2830 exhaustive); WIDE span>=15 is the residual, needs E-ADAPTIVE equidistribution (kps HYP-2842 refuted fixed-center; the correct route = three-distance = my L_y/decorrelation). 
- **MIDDLE-residual quantified (k=9):** max L_y over PRIMITIVE span>=15 = 0.4588 (margin 0.035 to cap), BINDING at the **1far config {0..7,21}** (consec core + one far element). The fully-wide configs are far safer (L_y 0.16-0.24). So the binding wide case = 1far.
- **The gap:** for 1far, L_y_decorr(consec core) ~ 0.44 (margin 0.05) + peel deviation Δ_w; actual Δ_w ~ 0.019 but the bound (6/49)V/w ~ 0.16 is **5-8x too loose**. Closing needs a sharper equidistribution bound on Δ_w (the arc-phase cancellation THM-546 left "not in closed form"). This is the shared crux with kps's resonance/marginal-uniformity work.


## BREAKTHROUGH (mac-mini-S29): square-root cancellation in the peel deviation => feasible wide closure
The peel deviation Delta_w = L_y(C u {w}) - L_y_decorr(C) (THM-546 Abel form:
Delta_w = -(1/w) Sum_j Sum_{arcs(a,b) in B_j} [F_j(wa)-F_j(wb)], |F_j|<=3/49) has
**SQUARE-ROOT CANCELLATION**, not the worst-case V-bound:
- F_j is the CENTERED sawtooth (mean 0). Evaluated at the scattered points {wa mod 1}
  over the ~V arcs, the sum has mean 0 + fluctuation ~ sqrt(V)*std(F_j), NOT V*||F_j||_inf.
- MEASURED |Delta_w * w|: consec_8 (V=28) -> 0.34 (= sqrt(28)*0.06); dilated (V=84) -> 0.78;
  2block -> 0.61. ALL ~ C*sqrt(V) << (6/49)V = 3.4-10. The 5-8x THM-546 looseness EXPLAINED
  as the difference between L2 (sqrt) and L1 (linear) discrepancy.
- For BOUNDED cores (span<=14), V<=~56 so sqrt(V)<=7.5 and |Delta_w*w| <= ~0.5 -- a SMALL
  CONSTANT. The peel cutoff w* = |Delta_w*w|/margin <= 0.5/0.047 ~ 11-17 (FEASIBLE).

## WIDE RESIDUAL decomposes by #far elements -- ALL CLOSE (k=9, verified float search)
- 0far (span<=14): DONE (HYP-2830 exhaustive).
- **1far** (bounded core + 1 far w): max L_y = 0.4588, MARGIN 0.035, binding = consec_8 + far@21.
- **2far**: max L_y = 0.3826, MARGIN 0.112 (more decorrelation = safer).
- 3+far: safer still (L_y 0.16-0.24). So the BINDING wide case is 1far, margin 0.035.

## CLOSURE PATH (the clean structure)
1far closure = [L_y_decorr(bounded core) <= cap-margin: a BOUNDED-CORE finite check] +
[|Delta_w*w| <= C*sqrt(V) < margin*w for w > cutoff~17: the sqrt-cancellation peel] +
[finite check w in [15,17]]. The remaining ANALYTIC gap: rigorous sqrt-cancellation
|Delta_w*w| <= C*sqrt(V) for the specific far w (generic w OK via Erdos-Turan/L2; RESONANT w
= kps's HYP-2842 domain). This CONNECTS my L_y peel route to kps's resonance/equidistribution
work: the closure = [sqrt-cancellation, generic] + [resonant-center handling, kps]. -> THM-546, HYP-2842.
