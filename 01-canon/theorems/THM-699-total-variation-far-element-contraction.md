---
id: THM-699
title: The total-variation far-element contraction — the far-offset correction to every seven-sector functional collapses by integration by parts with the PROVEN explicit rate C(E')/b, C(E') = (Σ_{e∈E'} e)·Σ_t C(7,t)t²(7−t)/7; the far-element plateau (HYP-2610/2644) is a theorem; iterated over the scale hierarchy it gives the explicit multi-scale decorrelation, closing the scale-separated branch of the wide-spread rows WITHOUT the divergent absolute lattice count
status: PROVED (one-line Riemann–Stieltjes integration by parts per sector atom; constants pure arithmetic; verified exact-rational with ZERO fitting: 40 random atoms worst |exact|/bound = 0.031, full measS7 on 4 families × 4 decades of b all inside bound, 3-scale iterate inside bound). HONEST SCOPE: the constant is universal but crude (≈4032·Σe' vs klein's measured ≤0.97) — it proves the PLATEAU AND RATE; the non-scale-separated wide middle (all ratios bounded, primitive) remains the open extremal core (it is exactly where the census-not-proof gap lives).
source: mac-mini-2026-07-09-S65 (cont.28-math, 2026-07-11)
depends_on:
  - THM-538   # the corrected kernel structure (short relations DO contribute — this bypasses the lattice sum entirely)
  - THM-534   # the seven-sector object
related:
  - HYP-2610/2644  # the stranger-contraction / far-element plateau — now PROVED with explicit rate
  - HYP-2653/2671  # the C(k) discrepancy program — this is its closed form
  - THM-687 (klein) # the measured C/V two-scale transfer — this is its proven-constant twin
  - opus-S216 reflection # the conditional-convergence wall — bypassed: TV needs no lattice count
---

# THM-699 — the total-variation far-element contraction

## Statement

Let E' be co-offsets, b a far offset, T ⊆ {0..6} a sector set, and
A_T(E) = ∫₀¹ ∏_{e∈E} 1_{frac(ex) ∉ U_T} dx the seven-sector atom. Then

>  |A_T(E' ∪ {b}) − (1 − |T|/7)·A_T(E')|  ≤  V_T(E')·‖H_T‖_∞ / b,

with V_T(E') ≤ 2|T|·Σ_{e∈E'} e (the total variation of the E'-product: each offset
contributes ≤ 2|T| jumps per unit period) and ‖H_T‖_∞ ≤ |T|(7−|T|)/14 (the sup of the
mean-zero antiderivative of the b-factor). Summing the inclusion–exclusion over the 127
sector atoms:

>  |meas(S7(E' ∪ {b})) − plateau(E')| ≤ C(E')/b,   C(E') = (Σ_{e∈E'} e)·Σ_{t=1}^{7} C(7,t)·t²(7−t)/7.

*Proof.* F = the E'-product is a step function with jump set {m/e + s/(7e)}; h(x) =
1_{frac(x)∉U_T} is 1-periodic with mean 1−|T|/7; H = ∫(h − h̄) is 1-periodic, piecewise
linear with slopes |T|/7 and −(1−|T|/7), so ‖H‖_∞ ≤ |T|(7−|T|)/14. Then
∫F(x)(h(bx) − h̄)dx = −(1/b)·Σ_i J_i·H(b x_i) by parts, and |Σ J_i| ≤ V(F). ∎

## Why this matters (the route map)

opus-S216 proved the t≥3 character route dead and localized the last analytic gap to the
signed relation-lattice sum, whose ABSOLUTE version diverges (conditional convergence — no
Minkowski point count can close it per se). This theorem BYPASSES the lattice: the far
element's entire signed cross-sum is captured by one integration by parts, signs and all —
the boundary term IS the conditionally-convergent tail, summed. Iterating over a scale
hierarchy e_{σ(1)} ≪ … gives the explicit multi-scale decorrelation (klein's THM-687/688
C/V transfers with PROVEN constants), closing the scale-separated branch of the k=8,9,10
wide-spread rows. The remaining extremal core is precisely the bounded-ratio primitive
middle — the same locus every route isolates (near-dilates/near-APs), now with a proved
moat around it.

## Files
04-computation/lrc14_tv_contraction_macmini_S65cont28.py (+ .out): exact-rational verification,
zero fitting — atoms (40 random, worst ratio 0.031), full measS7 (4 families × b ∈ {100..3000}),
3-scale iterate.
