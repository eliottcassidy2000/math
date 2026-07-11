---
id: THM-710
title: The factorial-moment eigen-transfer — falling-factorial moments are EIGENVECTORS of the far-element occupancy operator, m_r → ((7−r)/7)·m_r EXACTLY (3-line proof); the majorant value obeys V′ = (6/7)V + 1/7 − m2/84, and the rung arithmetic (6/7)·cap_{k+1} + 1/7 ≤ cap_{k+2} holds EXACTLY at every rung (slacks +0.038/+0.064/+0.093/+0.122) — so the ENTIRE moment ladder above the base propagates automatically and the wide-direction base checks reduce to {degree-3 at k=8} + {degree-2 at k=9} ONLY
status: PROVED (the eigen-identity: for p′_j = ((7−j)/7)p_j + ((j+1)/7)p_{j+1}, m′_r = Σ j^(r)p′_j = Σ p_j·j^(r)·[(7−j)+(j−r)]/7 = ((7−r)/7)m_r, using j·(j−1)^(r) = j^(r)(j−r); verified exact on 200 random rational distributions × r = 0..6, zero mismatches). Rung table exact rationals. The O(1/w) correction rides THM-699/700's TV bounds per atom (empirically |err|·w ≤ ~6 on real sweeps). CONSEQUENCE: k = 10,11,12,13 rows need NO independent core checks — they inherit from k=9 through the transfer; the moment-ladder residue is exactly [deg-3 bound on bounded 8-cores] + [deg-2/THM-705 linear bound on bounded 9-cores].
source: mac-mini-2026-07-09-S65 (cont.34, 2026-07-11)
depends_on:
  - THM-700/701 (kps)  # the transfer operator and recursion
  - THM-705 (mine)     # the deg-2 linear form this propagates
related:
  - THM-699/702 (the O(1/w) budget), opus-S220 (the ladder), kps cont.29 (THM-668∘707 supply composition)
---
# THM-710 — the eigen-transfer
p′_j = ((7−j)/7)p_j + ((j+1)/7)p_{j+1} ⟹ m′_r = ((7−r)/7)m_r (falling factorials diagonalize
binomial-type thinning). V = 1 − m1/2 + m2/12 ⟹ V′ = (6/7)V + 1/7 − m2/84 ≤ (6/7)V + 1/7, and
(6/7)cap_{k+1} + 1/7 = 7939/14014, 421/637, 487/637, 43/49 ≤ cap_{k+2} for k = 8..11 — every rung.
Files: 04-computation/lrc14_eigen_transfer_macmini_S65cont34.py (+ .out).
