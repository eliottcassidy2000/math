---
id: THM-3332
title: "AMM 12592: hypothesis S reduces to a one-row certificate (S-cone-fc); S verified through R = 32768; C* < (5+3g)/(3+2g) modulo certificate entry; and the exact-floor integer point exists at R = 256"
status: >
  PROVED + VERIFIED-EXACT / CROSS-AUDIT IN FLIGHT (final-referee lane running;
  E1's own lemma referee ran 600-trial exact certificates per lemma on an
  independent clamp implementation; the prior wave's fourth-implementation
  referee passed 108/108 checks on everything it reached).
  (S1, PROVED) Exact post-feed magnitude law: all-negative junk evolves by
  a'_t = max(0, (K_delta a)_t - 2C(d'-1,t)), a'_0 = max(0, a_0 - 2) — cell 0
  is an autonomous exact 2/row clock; subsumes negativity invariance and
  upgrades cap-ratio damping to an exact map. (S2, PROVED) Comparison
  principle: the map is cellwise monotone, so any capturing majorant proves
  capture. (S3, PROVED) The drain deadline is unconditionally non-binding.
  (S4, PROVED) Cell-1 spill is in-box iff a_0 <= d'-1 — explaining the
  |j_0| = d +- O(10) feed-end edge law (feed hands over at the
  drain-reinjection threshold).
  (S-cone-fc, PROVED — the main result) If ONE post-feed row satisfies
  (F1) negativity + support, (F2) a_0 <= d-1, (F3)
  2a_{t-1} + a_{t-2} <= 2C(d-1,t) for t in [2, m+2], then all cells are
  non-increasing, the support freezes, death is impossible, per-cell
  extinction clocks hold, and capture occurs by row R-2: S(R) FOLLOWS FROM A
  ONE-ROW CERTIFICATE. A-priori corollary (exact-rationally certified for
  all dyadic 2^9..2^40, both eps): F3 alone bounds all magnitudes
  (a_t <= C(d-1,t+1)) and the extinction clocks fit with ~35% slack, so the
  hypothesis is F1 ^ F2 ^ F3 alone for R >= 512.
  (VERIFIED-EXACT) fcscan 18/18 (R = 128..32768 x eps in {1/32, 1/16}): a
  certifying row exists at feed-end + 0..8 rows EVERYWHERE, margins ~38%;
  R = 32768 closed at BOTH eps (captures 19865/20185, each EXACTLY
  i_pf + a_0/2 - 1, the second predicted before the run finished; ballot
  debt (R-2)/2 exact; zero re-ignitions) — so S(1/32) and S(1/16) hold for
  ALL dyadic 128 <= R <= 32768.
  (CONSTANT, PROVED) 1 + gamma* + eps* = (5 + 3g)/(3 + 2g), certified
  rational bracket (1.6191617801, 1.6191618342): the E-lin limit constant
  is a single fraction in g = gamma*. Consequence chain: certificate ENTRY
  for all large R (the one remaining gap, lane running) => S(eps) for
  eps > eps* => C* <= (5+3g)/(3+2g) < 1.6191619 UNCONDITIONALLY.
  (RECORD CORRECTION + NEW POINT, VERIFIED twice) The (256, D0 = 0)
  exact-floor epoch is CLOSED (desc1 witness, sha 5950bd42, referee
  ALL-PASS twice; "D0*(256) = 1" was a plain-rule threshold, not a
  feasibility threshold). Its transportation form is exact: all 58,837
  cells f = (C - delta)/2 integer in [0, C] — the transportation polytope at
  the EXACT floor profile contains an INTEGER point at R = 256. The
  exact-floor existence frontier is now (512, D0 = 0), OPEN both directions
  (84/84 rule deaths and an action-alphabet schedule beam all die — search
  negatives, not evidence; the exact LP at 512 is out of current budget).
  (BARRIER LEMMAS, PROVED) L1: (t2 - 1) + F0 = d; L2: dominant front load
  = C(d, F0) = exactly half the box width (the front is born at box scale);
  L3: depth-rate law (A_{t-1}/A_t)/(C(d,t-1)/C(d,t)) = (R-1-t)/t ->
  g/(1-g) in (1.4874828, 1.4874887); L4: influence cone — junk at cell c
  after k rows depends only on cells [c-2k, c], so the feed cannot touch
  the front for ~0.2R rows. Window-locality alone is VACUOUS as a barrier
  class (a width-2 cap-reading rule can replay any precomputed global
  schedule), so the barrier is stated as two labeled conjectures
  (E3-BC(param), E3-BC(mech)) — not canon.
source: boxeph-2026-08-04-multifront-S
depends_on:
  - THM-3331
  - THM-3330
  - THM-3002
related:
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
notes:
  - 05-knowledge/results/amm12592-S-invariant-cone-proof-boxeph.md
  - 05-knowledge/results/amm12592-localrule-barrier-and-gap-boxeph.md
output: 05-knowledge/results/amm12592_S_cone_constants_boxeph.out
---

# THM-3332 -- S-cone-fc, the 32768 closures, and the 256 floor point

Statements, proofs, the 74 certificate files, and the barrier note are under
the frontmatter paths; every claim labeled, every computation exact. The
program state after this theorem:

1. **Upper bound:** `C* < (5+3g)/(3+2g) = 1.61916...` holds modulo ONE
   remaining lemma — certificate ENTRY (the feed phase delivers an
   F1 ^ F2 ^ F3 row for every large dyadic R), currently verified at all
   18 scale-eps pairs with the certifying row within 8 rows of feed-end.
   Entry-proof lane running.
2. **Attainment:** T(n) <= n + 1 + floor(gamma* n) + 14 for all n <= 2047;
   exact-floor closures now reach R = 256.
3. **Golden:** the exact-floor existence frontier is (512, D0 = 0); the
   general-class floor (MISTAKE-361 repair) remains the other half of
   C* = log_5(5 phi^2).
