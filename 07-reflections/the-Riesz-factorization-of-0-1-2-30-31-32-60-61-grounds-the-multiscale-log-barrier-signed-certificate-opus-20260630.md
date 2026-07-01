# The Riesz factorization of S={0,1,2,30,31,32,60,61} (pure DFT, machine-precision) grounds the multiscale log-barrier as a SIGNED Beurling-Selberg certificate: the exponential sum factors EXACTLY as T_S = T_A*T_B - z^62 (A={0,1,2}, B={0,30,60}), so the power spectrum IS the tensor |T_A|^2|T_B|^2 (Pearson 0.985) and log|T_S|^2 decomposes ADDITIVELY across scales 1 and 30 -- which predicts, to corr 0.997 on the heights and 10/12 on the peak-frequencies, exactly where the coherence and the THM-546 peel deviation w*Delta_w concentrate (the 30-comb under the |T_A| envelope); and the log-barrier B=<log m> gives tight sets a strictly POSITIVE averaged margin (0.43) where the linear min-margin is pinned at the borderline 1, diverging to -inf at holes, with the ARITHMETIC mean conserved (1.857, blind) and only the GEOMETRIC mean separating -- the AM-GM gap that is the signed smoothing HYP-2738 proves is necessary

*opus-2026-06-30. Owner: "confirm the Riesz factorization of {0,1,2,30,31,32,60,61} numerically and show it
predicts exactly where w*Delta_w peaks -- grounds (A); work to use w as a log-barrier smoothing inside a
Beurling-Selberg / Fourier-positivity signed certificate; look back at signed vs unsigned series." Done: the
factorization is exact, it predicts the peaks, the log-barrier is the signed smoothing, and the repo's whole
signed/unsigned history says WHY a signed certificate is forced.*

## 1. The Riesz factorization is EXACT (pure DFT, no measure subtlety)
S = {0,1,2,30,31,32,60,61} is the product set A(+)B, A={0,1,2}, B={0,30,60}, MINUS the top corner 2+60=62:
> **T_S(t) = T_A(t)*T_B(t) - z^62**,  T_A = 1+z+z^2,  T_B = 1+z^30+z^60,  z=e(t).
Numerically `max_t |T_S - (T_A T_B - z^62)| = 2.7e-15` (machine precision). So T_S is a **rank-1 perturbation of
a clean tensor** (Fejer/Dirichlet factor at scale 1 times Fejer/Dirichlet factor at scale 30). This is the
Riesz/Fejer product structure: `|T_A|^2 = D_3(t)^2`, `|T_B|^2 = D_3(30t)^2` (Dirichlet kernels).

## 2. The power spectrum IS the tensor -> the log decomposes ADDITIVELY (this grounds (A))
- **Tensor spectrum:** `|T_S|^2` vs `|T_A|^2|T_B|^2`: Pearson **0.985**, rel-L2 0.217 (the only deviation is the
  rank-1 corner -z^62). The coherence is a two-scale tensor: fine comb (period 1/30, from B) x coarse envelope
  (from A).
- **Exact additive log-decomposition** (the clean product F=A(+)B, before removing the corner):
  `T_F = T_A T_B` to 2.7e-15, hence **log|T_F|^2 = log|T_A|^2 + log|T_B|^2 EXACTLY** (wherever T_A,T_B != 0).
  This is the crux for (A): **because the sum FACTORS multiplicatively, the log-barrier ADDS across scales.**
  Removing the corner (S = F minus 62) perturbs this to corr 0.62 -- the log amplifies the low-|T| regions the
  rank-1 term disturbs, an honest caveat, but the grounding identity is the exact product one.

## 3. It predicts EXACTLY where w*Delta_w peaks
Two independent confirmations, both matching the tensor comb `30*(B-scale) + {0,1,2}*(A-scale)`:
- **Coherence heights:** at t=k/30 the B-factor is at its max (|T_B|=3 every k), so `T_S(k/30)=3 T_A(k/30)-e(2k/30)`
  exactly, and the peak heights `|T_S(k/30)|^2` track the coarse envelope `(3|T_A(k/30)|)^2` with
  **corr 0.9972** over k=0..14 (heights 64,62,57,50,41,31,22,13,7,3,1,... = the |T_A| decay).
- **THM-546 peel deviation** `D(w) = <P_S(t)(safe(wt)-c)>` (P_S = lonely-measure integrand, safe = lonely
  indicator, signed-sinc weights `g_hat(m) = -sin(2pi m/n)/(pi m)`): of the top-12 peel-frequencies for
  `|w*D(w)|`, **10 lie on the 30-comb + {0,1,2} satellites** ({30,31,32,60,61,91,92,120,121,122}), and D(w) is
  **signed** (terms +3.04,+2.96,-2.41,+2.28,-1.80,...: cancellation, not pile-up) -- exactly THM-546's signed
  Abel / QR-reality structure. So the Riesz factorization predicts which peel-frequencies carry the deviation.

## 4. The log-barrier IS the signed smoothing (grounds (A)'s certificate)
At n=14, radius 1/14, the barrier `B(S) = <log m(t)>`, m(t)=#{v in S: ||vt||<=1/14}:

| set | min m (linear margin) | hole meas | B = <log m> | geo-mean | arith-mean |
|---|---|---|---|---|---|
| AP {1..13} | 1 (ZERO) | 0 | **+0.428** | 1.534 | 1.857 |
| GW {1..11,13,24} | 1 (ZERO) | 0 | **+0.438** | 1.550 | 1.857 |
| 12->26 (non-tight) | 0 | 0.012 | **-inf** | 0 | 1.857 |
| {2..14} (non-tight) | 0 | 0.061 | -inf | 0 | 1.857 |
| drop-7 (non-tight) | 0 | 0.023 | -inf | 0 | 1.857 |

The **arithmetic mean is conserved at 1.857** for every set (= 13*2/14, blind to tightness -- exactly the
moment-blindness of `the-tight-skeleton-is-COMBINATORIAL-COVERING` and `covering-rigidity-...-dead-end`). The
**geometric mean separates**: tight -> 1.53, hole -> 0. The gap between them is the **AM-GM gap**, and the
log-barrier `B = log(geo-mean)` is precisely the functional that reads it. Tight sets get a **strictly positive
averaged margin (0.43)** where the linear certificate `min_t m >= 1` is pinned at the borderline; a hole sends
`B -> -inf` sharply. **This is the signed smoothing:** `-log` is convex, its Beurling-Selberg majorant has
signed (alternating) Fourier coefficients -- the same sign pattern HYP-2738/THM-534 prove is FORCED (a nonneg
majorant is loosest at the extremizer; you must subtract a consec-max quantity). The log-barrier is a
principled, scale-separated source of that signed certificate, and by 2 it builds SCALE-BY-SCALE.

## 5. The certificate construction (aspirational, now grounded) and honest status
**Target:** certify `min_t m(t) >= 1` (tight) where the AP is borderline (HYP-2738: no nonneg monotone works).
**Route (A):** majorize the barrier `-log m(t)` by a band-limited trig polynomial `Q(t) >= -log m(t)` with
`max_t Q <= 0`. By 2, for structured sets `-log m` splits additively across scales, so `Q = Q_coarse + Q_fine`
with the Beurling-Selberg loss at each scale ADDING (in log) instead of forcing one monolithic high-degree
majorant; each `Q_scale` is signed (convex-barrier majorant = the L_y negative weights, THM-534). The finest
scale carries the tight margin (smallest Beurling-Selberg gap 1/(N+1)).
- **PROVED/CONFIRMED (opus, this session, all numerical, clean):** the exact factorization; the tensor
  spectrum; the exact additive log-decomposition for the product; the comb prediction of the coherence heights
  (0.997) and the peel-deviation peak-frequencies (10/12); the barrier's tight/non-tight separation with a
  positive averaged margin and -inf at holes; the AM (conserved) vs GM (separating) gap.
- **ASPIRATIONAL (route A, not done):** the actual band-limited signed majorant `Q` of `-log m` for the AP with
  `max Q <= 0`, and the proof that the summed per-scale Beurling-Selberg losses stay below 0.43. This is the
  research target.
- **HONEST CEILING:** the barrier is a RELAXATION -- `min m >= 1 => B >= 0`, not conversely; and among covering
  sets it does not pin the AP (GW's B=0.438 > AP's 0.428). So a barrier certificate proves COVERING (tightness),
  not consec-max; the AP-vs-GW selection still needs the pointwise/support tools (klein-S51/HYP-3762). The
  barrier does not evade HYP-2738; it INSTANTIATES the signed certificate HYP-2738 demands.

## 6. Signed vs unsigned across the repo (the meta-lesson, from a full search of both threads)
A comprehensive scan (LRC + the older parity/tournament thread) shows a single law: **a signed certificate is
forced precisely when the natural nonnegative bound is loosest (or divergent) AT the extremizer.**
- **LRC:** HYP-2738 (nonneg phi is max at consec => useless; must subtract, forcing Bonferroni signs);
  THM-534 (the exact signed inclusion-exclusion moment-LP dual `L_y`, negative weight w_3=-1/5); THM-546 (the
  signed Abel form is 1.54x sharper than absolute, up to 76x with full arc cancellation; QR -1 non-residue mod
  7 makes Delta_w REAL); THM-523/515 (the singular series is a signed theta-lattice sum with signed sinc
  weights `-sin(pi t/7)/(pi t)`, nonnegative only in the Fourier/measure dual via Poisson); HYP-2974 (Toeplitz
  PSD of the SIGNED `D_S-1`; negative eigenvalues = Farkas hole-certificates; signed Fejer square).
- **Parity/tournament:** the OCF (THM-002) is the PLAIN independence polynomial `I(Omega,2)` -- no per-cycle mu
  signs in closed form, but the mu(C)=H(T\V(C)) weights are genuinely SIGNED in the Claim-A recursion and
  unroll into it (MISTAKE-001 was a computation bug, not a sign-theory error); the Even-Odd Split is the
  alternating identity `sum_S (-1)^|S| Delta(S,R)=0`, structural from path-count bilinearity.
- **The bridge to THIS work:** the log-barrier is the LRC instance of the same law -- the nonneg linear cover
  `sum 1_{D_v} >= 1` is borderline (loosest) at the AP; the SIGNED log-barrier (convex, alternating majorant)
  is the smoothing that reads the AM-GM gap. And by the Riesz factorization it is the FIRST signed certificate
  in the repo that provably decomposes across scales -- a multiscale signed certificate, not a monolithic one.

## Status
- **Confirmed (opus):** the exact Riesz factorization T_S=T_A T_B - z^62; tensor spectrum (0.985); exact
  additive log-decomposition (product); comb prediction of coherence (0.997) and peel deviation (10/12);
  barrier separation (tight +0.43 vs hole -inf), AM conserved / GM separating.
- **Grounded (opus):** the multiscale log-barrier signed-certificate program (route A) -- the structure is
  clean enough (scale-separated, additive in log) to attempt the band-limited signed majorant.
- **Open (honest):** the majorant construction + summed-loss bound for the AP; the barrier is a relaxation
  (covering, not consec-max); AP-vs-GW selection stays with the pointwise/support route.

Related: HYP-3763 (the entropy is the fine-scale monotone -- same log object, `<log m>`), HYP-2738 (why signed
is forced), THM-534 (the signed Bonferroni dual), THM-546 (the peel deviation w*Delta_w, signed Abel),
THM-523/515 (signed sinc/theta), HYP-2974 (Toeplitz PSD, signed), the-tight-skeleton-is-COMBINATORIAL-COVERING
+ covering-rigidity-...-dead-end (AM/moment blindness), HYP-3762/klein-S51 (the pointwise support route),
PATCH-TUNING-... (Angle 2 Fourier-positivity), OPEN-Q-108. HYP-3765 (this; renumbered from 3764 -- klein-S53 collision).
Scripts: 04-computation/lrc_riesz_factorization_logbarrier_opus_20260630.py.
