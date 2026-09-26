---
id: HYP-9140
title: "The pairing family also catches at the peak: delta_L = rho^peak_L L^(O(1)), so the pairing price of L-step provability is rho_L exp(-Theta(L^(1/3))), answering THM-4478's P1 negatively for the pairing family too"
status: >
  OPEN HYPOTHESIS for the consistent price delta_L; EMPIRICAL support for L = 8..64.
  UPDATE 2026-09-26 (later): the PRIVATE version, pi_L <= poly(L) rho^peak_L,
  is PROVED independently by crossroads223 (bridge note, THM-4488 reserved)
  and by our robin lane (Corollary 4). Our proof rests on the pairpeak
  note's Theorems B-C. The consistent price delta_L (interference between
  sources) remains OPEN.
  UPDATE 2026-09-26 (pairpeak lane, procgen_pairpeak_20260926_pairing_peak_price.md):
  the rationale below is CORRECTED. A single pairing flip does not give a
  fresh continuation. By the coupling Lemma 1 the flipped orbit (v-1)/2 and
  the unflipped T(v) re-merge at time r+3 whenever T(v)'s word begins
  1^r 00, which happens for half of all odd v (orchestrator spot-check to
  v = 4*10^5). A single highest flip rescues only 0.33, 0.25, 0.20 of
  undecided sources at L = 16, 32, 48.
  PROVED for the private price pi_L (each source its own flips, no
  consistency): pi_L = rho_L exp(-Theta(L^(1/3))), so P1 is negative
  there. The private version of this hypothesis follows from Conjecture R
  (HYP-9142).
  The real obstacle for delta_L is consistency: undecided orbits merge late
  into earlier frozen certificates. A consistent analysable construction
  gives valid sections (n <= 10^6 for L = 8..64; n <= 10^7 for L = 16, 32)
  at 2.0-5.3 x 2 rho^peak, a slowly growing factor; the peak lane's greedy
  is about 2x cheaper. Whether interference can be avoided at polynomial
  cost is exactly this hypothesis.
  Earlier status line: THM-4480 proves
  rho^peak_L / M*_L <= delta_L <= 2 rho_L, with the upper bound from
  THM-4475's construction. A partner-isolated greedy peak catch (flip the
  highest admissible odd path point v = 2 mod 3) gives valid members of P_L
  on n <= 10^6. Its flip density is 1.19-1.51 x 2 rho^peak for L = 8..32,
  while its ratio to 2 rho_L falls from 0.21 to 0.047.
source: collatz-procgen-20260922 session, peak lane (2026-09-26), section 7 of procgen_peak_20260926_peak_discounted_price.md
depends_on:
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md
---

# HYP-9140 -- the pairing price is peak-discounted

**Statement.** `delta_L <= poly(L) rho^peak_L`, where `delta_L` is the least flip density of a member of the pairing family (THM-4470/4475) in which every `n >= 3` falls below itself within `L` steps.

**Why it is plausible** (original rationale; CORRECTED 2026-09-26: a single flip re-merges with probability 1/2, so the freshness argument below applies only after the coupling has broken, and a construction needs repeated or barrier-type flips; see the pairpeak note).
* **Freshness.** A pairing flip at an odd point `v` replaces the step factor `3/2` by `1/2`. By Lemma F of the peak note (post-flip freshness), the next parities after a flip are uniform and independent of the prefix. So after a flip at the peak, the continuation is a fresh Collatz word with drift `ln(3/4)/2` per step.
* **Budget.** It has to fall by about `ln(w*/3) = O(L^(1/3))` for the words that carry `rho^peak`. It does so with probability tending to 1 if the flip leaves `>> L^(1/3)` steps.

**What a proof must control.**
* Flips that are not at stopping times.
* Partner interactions: the flip pushes `2i` up to `3i`. The `v = 2 (mod 3)` isolation device of THM-4475 keeps partners out of other orbits.
* Frozen-pair blocking, i.e. the greedy never getting stuck. The data show two stranding failures of weaker versions, at `n = 890154` and `735217`; the isolated version has none up to `10^6`.

**Consequence.** Together with THM-4480 it would make `delta_L = rho_L exp(-Theta(L^(1/3)))`. The literal P1 of THM-4478 (for the pairing family) would then be negative, and the three settings would read:
* arbitrary edits and pairings: peak-discounted;
* the periodic cube (THM-4479): polynomially sharp.
