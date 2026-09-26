---
id: HYP-9142
title: "Conjecture R (Robin inequality): the survivors of the reflected barrier at slope level m are at most a constant times the undecided words with peak below 3^(m+1), uniformly in L and m; this implies the peak-discounted private pairing price pi_L <= O(L) rho^peak_L"
status: >
  PROVED UP TO A POLYNOMIAL FACTOR (2026-09-26), independently by two
  sessions; OPEN with a constant K_0.
  (a) codex crossroads223 (crossroads223_20260926_robin.md and _bridge.md;
  THM-4488 reserved): N_m(L) <= 4e(m+3) L A_(m+5)(L).
  (b) collatz-procgen robin lane (procgen_robin_20260926_robin_inequality.md;
  orchestrator-audited, including an independent interval-arithmetic
  re-check): N_m(L) <= (m+2)^17 A_(m+2)(L). This is Theorem 1 (a sine
  supersolution at width m+2 with the exact eigenvalue kappa(pi/(m+2)))
  together with Theorem 2 (a sine subsolution plus bridge lemmas).
  Both give the private pairing price pi_L <= poly(L) rho^peak_L and the
  sharp second-order term ln rho^peak_L = -(1-H)L ln 2 - kappa_3 L^(1/3)
  + O(log L), elementary, without Mogul'skii.
  Earlier status: OPEN; FINITE-EXACT support. For 1 <= L <= 180 and all m,
  max N_m(L)/A_(m+1)(L) = 1.003255 (at (L,m) = (50,7)); with K = 2 the
  maximum is 1.0000889. The ratio is <= 1.0005 at L = 256, 512, 1024
  (m <= 60). So K_0 = 1 is false but K_0 = 1.0033 holds on the range.
  PROVED consequence (Theorem C of the pairpeak note): Conjecture R implies
  pi_L <= (K_0 3^K (27L+18) + 4L 2^(-L)) rho^peak_L, i.e. HYP-9140 for the
  private price. Modulo Mogul'skii (THM-4480 Theorem 4) it also gives the
  sharp constant ln pi_L = -(1-H)L ln 2 - kappa_3 L^(1/3)(1+o(1)).
source: collatz-procgen-20260922 session, pairpeak lane (2026-09-26), section 5 of procgen_pairpeak_20260926_pairing_peak_price.md
depends_on:
  - 05-knowledge/hypotheses/HYP-9140-pairing-price-is-peak-discounted.md
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md
---

# HYP-9142 -- the Robin inequality for the barrier walk

**Definitions.**
* `A_M(L)` counts the undecided words (every prefix slope `3^(e_j)/2^j > 1`) whose peak slope is `< 3^M`.
* `N_m(L)` counts the survivors of the reflected barrier `Pi^(m)`. This construction flips every odd point of its own orbit at slope level `>= m - 1`, so the top boundary is partially absorbing: in the tilted measure a zone visit has weight `2(1-c) = 0.738`.

**Conjecture.** `sup_(L >= 1, m >= 1) N_m(L)/A_(m+1)(L) < infinity`.

**Reading.** A reflected (Robin-type) top boundary costs at most a constant compared with a hard wall one unit higher. This is a discrete boundary-value comparison for a random walk in a strip, of the kind whose continuum analogue compares Robin and Dirichlet principal eigenvalues. It isolates the analytic core of the private pairing price. The consistent price `delta_L` (HYP-9140) additionally needs interference control.
