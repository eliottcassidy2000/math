---
id: HYP-3844
title: The K=0 final-window lemma -- d=1 emptiness (PROVED, formalized), diam<=28 => K=0 (PROVED, formalized), exposure => band-critical residue system at Q>=29 (PROVED, exact finite test); universality for covering sets REFUTED (real exposed window kinks exist at large diameter)
status: MIXED -- (a,b,c) PROVED (b,c arithmetic shells FORMALIZED in Lean, build clean); universal K=0 for covering sets REFUTED by explicit examples; the per-set band test is cheap and exact
source: klein-2026-07-01-S89
script: 04-computation/k0_final_window_lemma_klein.py (+ .out); Lean: 04-computation/lean/TournamentH7/TournamentH7/LRCFinalWindowBand.lean (builds, sorry-free)
related:
  - HYP-3841 (tangent-ladder; this lemma is its last-rung input)
  - THM-592 (kink taxonomy), MISTAKE-093 (the r*=2/29 correction found en route)
  - mac-mini HYP-3850(4) (K=0 genericity empirics -- consistent)
  - opus HYP-3902 / doc sec 7.2 (three-rational certificate; this supplies its K-test)
---

# HYP-3844: the K=0 final-window lemma

## Proved

(a) [THM-592, cited] Concave (overtaking/engulfment) kink radii of Lambda_S are d/(w-v),
convex (gap-death) radii are d/(v+w), v < w in S, d >= 1.

(b) **d=1 emptiness + diameter corollary** (FORMALIZED, LRCFinalWindowBand.lean):
r0 = d/(w-v) in (1/15, 1/14) iff 14d < w-v < 15d; for d=1 the band (14,15) is EMPTY;
hence any concave window kink has d >= 2 and w-v >= 29, and **diam(S) <= 28 => K = 0 on
(1/15, 1/14)**. Covers the entire 11-core census (range <= 14), AP/GW, all opus tower
patterns -- the regime the census-exhaustiveness doc needs. Bands: {29}, {43,44},
{57,58,59}, ... (band d has d-1 members; also formalized).

(c) **Exposure => band-critical residue system** (proved this session): if the crossing
of same-side endpoints of a/v, b/w (d = bv - aw) is EXPOSED at r0 = d/(w-v), then the
coincidence point is x0 = (b-a)/(w-v) EXACTLY (key identity: (b-a)/(w-v) - a/v =
d/(v(w-v)) = r0/v), with ||v x0|| = ||w x0|| = r0 and m_S(x0) = r0. In lowest terms
x0 = B/Q: m_S(B/Q) = D/Q with Q/D in (14,15), Q >= 29, all residues of S*B mod Q at
distance >= D, v = w mod Q both attaining. **An exact finite per-set test** (enumerate
band pairs, check m_S at the crossing rationals) -- the K-entry generator for opus's
three-rational certificate (m, m', K).

## Refuted (the honest boundary)

Universal K=0 for covering sets is FALSE: the adversarial hunt (planted band differences)
found REAL exposed concave window kinks, e.g. S = [9,12,16,20,26,33,42,44,54,70,77,86,107]
kinks at r0 = 5/74 ((33,107), d=5, w-v=74 in band), slopes -3.769 -> -3.811. So the
ladder's last rung is defect-free by THEOREM only in the bounded-diameter regime; beyond
it the band test is the (cheap) certificate. Deep well: all 14 window candidates
((d=12, w-v=170..179)) are SUBMERGED (m(x0) < r0), K=0 confirmed structurally.

## Also en route: MISTAKE-093

The (d)(i) reconciliation probe found GW's kink at 2/29 is REAL: my S88 r* = 1/15 was
wrong; identity window is [2/29, 1/14]; mac-mini S94(3) confirmed in full. Filed
MISTAKE-093; HYP-3843 corrected in place.
