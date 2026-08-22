---
id: THM-3331
title: "AMM 12592: the E-lin conditional theorem (C* < 1.62924 under one hypothesis S), the closed-form feed-survival constant eps*, the LIFT theorem, and the super-block death law"
status: >
  PROVED + VERIFIED-EXACT / AWAITING CROSS-SESSION HOSTILE AUDIT. Referee:
  12/12 witnesses ALL-PASS on a FOURTH independent implementation
  (streaming O(M)-memory verifier via t = (1-x)/x); all core theorems below
  CONFIRMED with two corrections of record (G3 seeded-chain valid rows are
  46/95 ~ 0.09R, not 25/42 — a majorant-truncation artifact, theorem sound;
  an immaterial i_feed off-by-one convention, re-instantiated with no
  downstream impact).
  (LIFT, PROVED) Degree-raising delta'_t = delta_t + delta_{t-1} maps an
  admissible block at d to an admissible block at d+1 representing the same
  polynomial: epoch feasibility is MONOTONE in D0, and any closing
  construction may replace rule A.
  (B, PROVED) Feed-phase survival: i_feed = floor((R(1-g) - D0)/(1+g))
  exactly (never an integer: phi^{2q} = 5^p is impossible), and for
  D0 >= eps* R in the window R/2 < d_0 < 2R/3, rule A cannot die at any row
  <= i_feed, with the NEW closed-form constant
  eps* = 2(1 - g - g^2)/(3 + 2g) = 0.0211741... (g = gamma*). Constant
  ordering: eps_phi = 1/phi - g = 0.02005 < eps* = 0.02117 < eps_inf
  (measured ~ 0.025).
  (C, PROVED) Post-feed one-step lemmas: negativity invariance; exact
  cap-ratio damping; two-cell front freeze; ballot-debt identity
  j_0 = e_i(1) = (2-R) + 2#minus2rows with 2/row drain and a NEVER-BINDING
  deadline (i_feed < (R-2)/2).
  (S, HYPOTHESIS — the single open item) For dyadic R >= 32768 at
  D0 = ceil(eps R), the autonomous post-feed flow captures by row R-2.
  VERIFIED-EXACT for ALL dyadic 128 <= R <= 16384 at eps = 1/32 and 1/16
  (16/16 closures, debt (R-2)/2 exact), plus probes to eps ~ 0.4. Feed-end
  state (16/16): all-negative contiguous junk block on cells [0, m],
  m ~ 0.93 sqrt(R); cells >= 1 have R-INDEPENDENT amplitudes
  (max log2 ~ 10-11 over seven doublings); |j_0| = d +- <= 11 (a second,
  unexplained edge law at the T7 drain edge).
  (CONDITIONAL BOUND, arithmetic referee-confirmed)
  under S(1/16): C* < 435287/262144 = 1.66049...;
  under S(1/32): C* < 427095/262144 = 1.62924...;
  and S for every eps > eps* would give C* <= 1 + gamma* + eps* = 1.61916...
  (G1, PROVED) Alternation-transport calculus: for j_t = (-1)^{d-t} a_t,
  (K_delta * j)_c = (-1)^{d-c} (D^{1+delta} a)_c exactly (D = backward
  difference); L1(K*j) = 2 x (one-signed mass of D^{1+delta} a) — the exact
  survival budget; the naive envelope form |j| <= M => |K*j| <= |D^2 M| is
  FALSE (counterexample recorded).
  (SUPER-BLOCK LAW, VERIFIED-EXACT 12/12) Death <=> formation of a
  front-region same-sign junk run whose min |value| exceeds the entire
  remaining cap tail; forms by row ~ 0.09R (long before the visible march),
  persists with zero recoveries; closures have zero super rows. Bulk-rule
  sweep: best local variant (desc1) improves D0* to 0/0/4/14/36/87/<=188
  (R = 128..8192) — real gains, but the LINEAR signature persists in the
  whole clamp-local class; attainment improves to
  T(n) <= n + 1 + floor(gamma* n) + 14 for all n <= 2047 (R = 1024 closes
  at D0 = 14, referee-verified).
source: boxeph-2026-08-03-multifront-estimateE
depends_on:
  - THM-3330
  - THM-3329
  - THM-3002
  - THM-3026
related:
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
notes:
  - 05-knowledge/results/amm12592-Elin-theorem-boxeph.md
  - 05-knowledge/results/amm12592-golden-transient-bound-boxeph.md
  - 05-knowledge/results/amm12592_bulkrule_design_sweep_boxeph.md
referee: 04-computation/amm12592_estimateE_referee_boxeph.py
output: 05-knowledge/results/amm12592_estimateE_referee_boxeph.out
---

# THM-3331 -- E-lin conditional theorem, eps*, LIFT, and the super-block law

Full statements and proofs are in the three companion notes (frontmatter),
each claim labeled PROVED / VERIFIED-exact / CONJECTURED, all computations
exact, and the whole wave re-audited by a fourth-implementation hostile
referee (64 PASS / 42 CONFIRMED / 0 FAIL).

The remaining program for the golden bound:

1. **Prove S(eps)** (post-feed collapse; in progress — invariant-cone and
   super-block-impossibility routes): S for eps > eps* yields
   `C* <= 1 + gamma* + eps* = 1.61916...` UNCONDITIONALLY.
2. **The local-rule barrier**: every clamp-local rule tried has linear
   slack; exact-floor witnesses exist at R <= 128 — the gap between rule
   dynamics and epoch feasibility is the true golden frontier.
3. **Lower bound**: the deadline-bounded routing window (MISTAKE-361
   repair) remains the other half of C* = log_5(5 phi^2).
