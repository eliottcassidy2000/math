---
id: THM-712
title: The k=8 cubic base form — the optimal degree-3 majorant of g on {0..7} is EXPLICIT (touches N = {0,1,3,7}; q₃ = 1 − (2/3)N + (47/252)N(N−1) − (5/252)N(N−1)(N−2), exact 70-vertex LP), giving the k=8 requirement (2/3)m1 − (47/252)m2 + (5/252)m3 ≥ 1 − cap₉ with m3 entering POSITIVELY (triples need only a LOWER bound); adversarial hunt (50 families, every enemy class): worst bound 0.4483 at SHIFTED-CONSEC {1..8} vs cap₉ = 0.4943 — unconditional margin +0.0459, TWICE the k=9 margin
status: MAJORANT PROVED (6 feasible vertices of the exact deg-3 LP; the active vertex at consec touches {0,1,3,7}); consec_8 optimal bound 40561/92610 = 0.4380, margin +0.0563; the unconditional inf with minimizer {1..8} is the sharp conjecture (0 violations over 50 adversarial families; dilate-invariant; far elements handled by THM-710's eigen-transfer with m3 → (4/7)m3). TOGETHER WITH THM-711: both wide-direction base checks are now single compact moment inequalities with conjectured-exact shifted-consec minimizers and margins +0.026/+0.046 — and the m3 sign means kps's 3D box count is needed only as a LOWER-bound instrument, the easier direction.
source: mac-mini-2026-07-09-S65 (cont.36, 2026-07-11)
depends_on:
  - THM-703/705/710/711 (the ladder machinery)
related:
  - opus-S220 (deg-graded ladder), kps 3D box (now lower-bound-only), klein-S252 (pigeonhole DE-CITED in Lean — the citation assembly is down to two)
---
# THM-712 — the cubic base form
Exact LP over quadruples: 6 feasible vertices; active at consec: q₃(N) = 1 − (2/3)N +
(47/252)N(N−1) − (5/252)N(N−1)(N−2) ≥ g on {0..7} (contact {0,1,3,7}). Requirement:
(2/3)m1 − (47/252)m2 + (5/252)m3 ≥ 1 − cap₉. Hunt: 50 families (consec/shifted, doubling,
near-AP, mod-7, random, far-mix): worst 0.4483 at {1..8}, margin +0.0459 unconditional.
Files: 04-computation/lrc14_cubic_base_macmini_S65cont36.py (+ .out).
