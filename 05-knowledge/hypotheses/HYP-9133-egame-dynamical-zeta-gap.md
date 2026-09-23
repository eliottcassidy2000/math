---
id: HYP-9133
title: "Resonance gap for the dynamical zeta of the canonical E-game escape Psi"
status: >
  SETTLED -- PROVED (orchestrator, 2026-09-23; proof in
  collatz_procgen_20260923_trunk_rh.md section 9). The leading pole x_0 of
  zeta_Psi = 1/(1 - W(3^-s; theta)) is simple and is the only pole with
  |x| <= x_0. Every other pole with |3^-s| < 1 satisfies
  Re s <= s(theta) - eta, uniformly for theta in compact subsets of
  (0, infinity), modulo the 2 pi i/log 3 period. The proof uses positive
  coefficients with c_1 > 0 (Pringsheim plus aperiodicity), isolated zeros,
  and Hurwitz plus simplicity. Numerically the gap is 0.683, 0.515, 0.266
  and 0.094 at theta = 1, 2, 4 and 8.6434. The Ihara-type "RH" remains
  refuted numerically.
source: collatz-procgen-20260923 trunk/RH lane (candidate C2)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260923_trunk_rh.md
  - 05-knowledge/results/collatz_procgen_20260922_q2_endgame.md
---

# HYP-9133 -- a spectral gap for the E-game

**Refutation form:** resonances of `zeta_Psi` accumulate at the leading pole
for some `theta`.
