---
id: THM-750
renumber_note: "was THM-748 (claimed opus-S283); klein-S303 independently renumbered their shadow-gap THM-744 -> 748 (crossing my claim) and took 749; ceding to their completed canon+Lean rename, mine -> 750"
title: THE CLOSED BUDGET of the tail lane -- W(L-Area) = PHI(W) + QPOT(W) + KAP(W) EXACTLY (Fraction-equal vs the engine at all battery points, both shapes): the wedge sums, the exact quadratic pot, and closed-form per-event kappa's (death/birth/swap/cut, incl. static-arc-0 partners); full-period scans: signed envelope 0.0720 / 0.3486, all-W termwise bound 3.01/W; every constant a rational or a completed scan -- U1 discharged, no analytic unknown remains in the lane
status: PROVED (the identity, referee-gated: EXACT Fraction equality at W in {600,601,700,977 / 930,1001}) + SCANNED-EXACT (full periods Q = 8820 / 97020; extremes exact-confirmed) + the all-W extension by termwise decay (N_max = 3.0164 / 3.0093)
source: opus-2026-07-14-S283 (owner prompt: run the U1 scan and close the lane's budget)
depends_on:
  - THM-745/746/747 (the identity chain, phase sums, periodicity); THM-742/743 (assembly); S273 sweeps (below-threshold coverage)
related:
  - MISTAKE-142 (the unsound-order lesson -- here every order is an exact identity)
  - klein-S302 (HYP-6650 fleet triangulation -- this lane's contribution is now a closed exact object)
---

# THM-750 -- the closed budget

## The exact identity (referee-gated)

For W >= Wz (every segment has a full-inside crossing; Wz = 924 / 588):

>  W (L - Area)  =  PHI(W) + QPOT(W) + KAP(W)   EXACTLY,

- PHI: per-segment wedge sums (j/W) Sum psi(h_m) over full-inside crossings (arithmetic series,
  wrap-free by the strong no-wrap lemma);
- QPOT: the exact quadratic pot Sum orient h_m j^2/(W(W+j)) (per-crossing f - g =
  (j/W)psi(h) + h j^2/(W(W+j)) with NO remainder);
- KAP: closed-form per-event corrections at the segment-end events, phi = {u_event W}:
  death: W delta (phi-x)_+/((W+j_up)(W+j_lo)) - delta phi^2/(2W);  birth: mirror;
  swap: o[sigma_valid - x - (j_A phi^2 - j_B(1-phi)^2)/(2W)];  cut: per-run clipped forms;
  static-arc-0 partners for the origin-band events; phi = 0 events contribute 0 (boundary
  degeneracy -- covered by the full-inside sums).

REFEREE: exact Fraction equality W(L_engine - Area) = PHI + QPOT + KAP at all six battery
points (shape 2: W = 600, 601, 700, 977; shape 1: W = 930, 1001).  Unmatched ends: 0.

## The scans (U1 discharged)

| shape | Q | signed envelope max|E| (exact-confirmed) | at W | all-W termwise N_max | at W |
|---|---|---|---|---|---|
| 1 | 97020 | **0.071977** (W=1052) / +0.070180 (W=1485) | -- | **3.0093** | 1094 |
| 2 | 8820 | **0.348627** (W=729) / -0.347155 (W=606) | -- | **3.0164** | 1054 |

>  **|L - Area| <= 3.02/W for ALL W >= Wz** (termwise decay at fixed phase class);
>  **|L - Area| <= 0.072 / 0.349 / W on the scanned window** (the sharp envelope).

Margins at Wz itself: 23x / 9x below Area (and improving as 1/W).  Combined with the S273
sweeps (every W <= 1948 / 2676 verified exactly) and THM-743 (sound for W >= 339/513), the two
shape families are certified at EVERY W with the tail's exact envelope known.

## The closure statement

Every quantity in the lane's error budget is now (i) an exact rational (the identity's
coefficients, Xi_sv, the event data), or (ii) a completed finite scan (S(W), T(W), the total
envelope E(W), N(W)).  U1 of THM-747 is DISCHARGED: no analytic unknown remains in the tail
lane.  What remains outside the lane: the compact core (kps) and the general multi-speed
equidistribution (klein-S300) -- as mapped in THM-747.

## Debugging ledger (methodological, for the Lean pass)

The referee-driven build caught three real errors: (1) events at phi = 0 (vertex exactly on a
strip boundary) must contribute 0 -- the post-event segment's full-inside range already covers
the strip (first symptom: pred = 2x true on aligned-pair swap strips at divisible W); (2) cut-run
reconstruction must pad with the STATIC arc-0 boundaries (runs bounded by the origin band);
(3) Python int/int division leaked a float 0.0 into the Fraction chain (max(x,0) vs max(x,F(0)))
-- 2^-55 noise masquerading as a formula error.  The per-strip bisection against exact
f(m), g(m) localized (1)-(2) in one pass.

## Files

04-computation/lrc14_closed_budget_thm748_opus_S283.py (+.out, incl. scan appendices).
