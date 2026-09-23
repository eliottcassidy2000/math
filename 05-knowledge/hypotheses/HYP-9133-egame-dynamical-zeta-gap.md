---
id: HYP-9133
title: "Resonance gap for the dynamical zeta of the canonical E-game escape Psi"
status: >
  OPEN; floating-point evidence. Every non-leading pole of zeta_Psi(.; theta)
  with |3^-s| < 1 satisfies Re s <= s(theta) - eta, with eta > 0 uniform for
  theta in compact subsets of (0, infinity). Here zeta_Psi is the E-game's
  own dynamical zeta, in which the 1/2-escape price enters as a factor
  1 + 2^theta. Numerics: the gap is 0.64150 - 0.23547 = 0.406 at theta*,
  and 1 - 0.90647 = 0.094 at theta_L. The trunk/RH lane found this to be the
  only precise form of the "hostile point 1/2 vs critical line" analogy. The
  naive Ihara-type Riemann hypothesis for zeta_Psi is REFUTED numerically.
source: collatz-procgen-20260923 trunk/RH lane (candidate C2)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260923_trunk_rh.md
  - 05-knowledge/results/collatz_procgen_20260922_q2_endgame.md
---

# HYP-9133 -- a spectral gap for the E-game

**Refutation form:** resonances of `zeta_Psi` accumulate at the leading pole
for some `theta`.
