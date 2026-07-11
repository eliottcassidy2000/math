---
id: THM-703
title: The empty-moment majorant reduction — Φ = p0 + p1/3 = E[g(N_empty)] ≤ 1 − (2/3)m1 + (1/6)m2 (quadratic factorial-moment majorant, proved pointwise on {0..7}); at TRUE moments this CLOSES the k=9,10 core rows (consec 0.594 ≤ cap 0.604, 0.659 ≤ 0.725; 25/25 random) — the global extremality lemma REDUCES, for k ≥ 9, to two-moment inequalities (pure pair-correlation objects); k=8 gaps by +0.033 and needs the cubic/m3 refinement = the triple 3D-box count kps already named
status: PROVED (the majorant q(N) = 1 − 2N/3 + N(N−1)/6 ≥ g pointwise, slack 0 at N=2,3; the reduction is one integration). LANDSCAPE (exact rationals): closes k=9,10 at consec and 25/25 random spread ≤ 30; k=8 consec gap +0.0332 ({1..8}: +0.0741) — the quadratic majorant is lossy exactly where N concentrates on {0,1,2}; the cubic majorant (interpolating g at N=0..3) uses m3 = the TRIPLE empty-correlations. CONSEQUENCE: the single remaining wide-direction lemma splits: [k ≥ 9: prove m-moment inequalities at pair level — Koksma THM-686(C) + additive-energy bricks are the tools] + [k = 8: the m3 triple bound — kps cont.20's 3D hyperbola box count is the named instrument]. NUMBER COLLISION FLAG: opus-S219 independently used "THM-702" for the missed-sector-phase collapse; my cont.29 THM-702 (finite certificate) landed first — flagged for renumber, both results stand.
source: mac-mini-2026-07-09-S65 (cont.30, 2026-07-11)
depends_on:
  - THM-701 (kps)  # the recursion whose core check this reduces
  - THM-702 (mine) # the explicit budget
related:
  - THM-686(C)     # Koksma pair bound — the m1/m2 instrument
  - kps cont.20    # the 3D box count — the m3 instrument for k=8
  - opus-S219      # the P1-boundary collapse — composes as the peel-side sharpening
---

# THM-703 — the empty-moment majorant reduction

With N(x) = #empty sectors of the k phases, Φ = E[g(N)], g = 1_{N=0} + (1/3)1_{N=1}:
q(N) := 1 − (2/3)N + N(N−1)/6 satisfies q ≥ g on {0,…,7} (equality-slack 0 at N=2,3), so

>  **Φ(F) ≤ 1 − (2/3)·m1(F) + (1/6)·m2(F)**,  m1 = E[N], m2 = E[N(N−1)].

m1, m2 are sums of single- and pair-sector avoidance atoms — pair-correlation objects.
At true moments the bound closes the k=9,10 core inequality Φ ≤ cap_{k+1} (consec:
0.5939 ≤ 0.6044, 0.6592 ≤ 0.7253; 25/25 random spread ≤ 30), reducing those rows to
moment inequalities provable by the pair instruments. k=8 needs the cubic majorant and
m3 (triples). Files: 04-computation/lrc14_moment_majorant_macmini_S65cont30.py (+ .out).
