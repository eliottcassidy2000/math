---
id: THM-705
title: The linear moment requirement — the optimal degree-2 majorant of Φ is EXPLICIT, q*(N) = 1 − N/2 + N(N−1)/12 (touches g at N = 0,3,4; proved by the 8-vertex exact LP), so the k=9,10 core rows each reduce to ONE LINEAR pair-correlation inequality: m1/2 − m2/12 ≥ 1 − cap_{k+1}; uniform adversarial margins +0.0263 (k=9, extremal {1..9}) and +0.0943 (k=10, extremal {1..10}); cross-validates opus-S220's ladder values (0.5668/0.6307) exactly
status: PROVED (majorant: q* ≥ g on {0..7} by the 8 feasible vertices of the exact LP — q* is the active vertex at both rows' moment points; the reduction is one integration). The remaining content of the k=9,10 rows is the SINGLE linear inequality's uniform validity over bounded cores (margin 2.6% at the binding k=9; observed extremizer the shifted consec {1..k}). k=8 needs degree-3 (m3) per opus-S220 — kps's 3D box count is the named instrument.
source: mac-mini-2026-07-09-S65 (cont.31, 2026-07-11)
depends_on:
  - THM-703 (mine)   # the majorant reduction this optimizes
  - opus-S220        # the degree-graded ladder (values cross-validated exactly)
related:
  - THM-704 (opus, renamed from the 702 collision)  # P1-boundary collapse
  - kps cont.20/24   # pair-correlation instruments for the linear inequality
---

# THM-705 — the linear moment requirement

The exact LP min of a + b·m1 + c·m2 over quadratics ≥ g on {0..7} has 8 feasible
vertices; at both rows' moment points the active vertex is

>  **q*(N) = 1 − N/2 + N(N−1)/12**  (g-contact at N = 0, 3, 4),

so Φ ≤ 1 − m1/2 + m2/12, and row k ∈ {9,10} closes iff the LINEAR inequality

>  **m1(F)/2 − m2(F)/12 ≥ 1 − cap_{k+1}**

holds over bounded k-cores. Exact margins: consec +0.0376/+0.0945; adversarial worst
(30 random + structured) +0.0263 at {1..9}, +0.0943 at {1..10}. One linear
pair-correlation inequality per row, 2.6% slack at the binding row — the final
pair-level form of the wide-direction residue.
Files: 04-computation/lrc14_linear_moment_req_macmini_S65cont31.py (+ .out).
