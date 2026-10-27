---
id: HYP-9422
title: "After THM-4814: (1) every regular tournament on 2m+1 >= 7 vertices is exactly m black native moves from a singleton plus a strong core; (2) the nuclear-norm maximizers at n = 1, 2 (mod 4) have two distinct nonzero singular values and are one-vertex deletions / bordered forms of the odd-order maximizers; (3) the exception count of P_q against the join/C3 tree grammar is q^2/4 - o(q^2), no bank beating the transitive parse asymptotically; (4) order-11 cycle-index types would deliver the order-12 ledger (A280772(12), c(12)) from A000568(12) = 154,108,311,168"
status: >
  OPEN (klein-2026-10-27). (1) lower bound PROVED (THM-4814 (D)); witnesses FINITE-EXACT at n = 7, 9 (every regular class)
  and for P_11 plus twelve random regular 11-tournaments; the all-order statement needs a construction (alternate the
  singleton as start and end, keep the core strong, keep every path black); n = 5 is excluded by the pure-blue atom T5.
  (2) FINITE-EXACT at n = 5, 6, 9 (ratios to the Cauchy-Schwarz bound 0.9789, 0.9636, 0.9842; at n = 9 all 216 maximizers
  share the singular values sqrt7 (x6), sqrt15 (x2)); the proposed description is untested at n = 10.
  (3) heuristic: the spectral bound gives q^2/8 for every fixed bank and the score bound q^2/8 for linear parses, while
  f(P_q) = q^2/4 - O(q^{3/2}) for the transitive parse; no construction with fewer than (1/4 - o(1)) q^2 exceptions is known.
  (4) a task: 903,753,248 classes at order 11 by gentourng in parts (hours on eight cores), then species11.py at n = 12.
session: klein-2026-10-27 (codex/session-overnight-20261027)
source: 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md (§2–§4)
scripts:
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/blackcore.py (+ receipts/blackcore7-*, blackcore9-*, blackcore11paley-*)
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/nucmax.py, nucmax9sv.py, parse_bound.py
related:
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 05-knowledge/hypotheses/HYP-9420-strong-untouchables-1-7-21-complement-graph-components-paley-depth-quadratic-bound.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
  - 05-knowledge/hypotheses/HYP-9413-minimum-intransitivity-chain-family-fibonacci-primes-and-the-intransitivity-count-asymptotics.md
---

# HYP-9422 — Black distance of regular tournaments, nuclear-norm maximizers off the Hadamard orders, and the next ledger

1. **Exact black distance.** For every regular tournament on `2m+1 >= 7` vertices there is a word of
   `m` black native moves ending at a tournament with a vertex of score `0` or `2m` whose deletion is
   strong, and no shorter word (black or blue) reaches such a tournament. The lower bound is THM-4814
   (D); the geodesics alternate the singleton as start and end of the paths (the "two distinct path
   starts"); every regular class at orders 7 and 9, `P_11` and twelve random regular 11-tournaments
   have such words (found within 60 candidates per level). Test next: all 1,223 regular classes at
   order 11; whether the first move can always be chosen with a prescribed singleton `a` (true in every
   witness so far: the search fixed `a = 0`).

2. **Maximizers off the Hadamard orders.** For `n ≡ 1, 2 (mod 4)` the maximum of `ν` is attained by
   tournaments with exactly two distinct nonzero singular values; at `n = 9` they are `√7` (×6) and
   `√15` (×2) (216 classes). Conjecturally these are vertex deletions or borderings of the order-`n∓1`
   extremals; `n = 10` (9,733,056 classes, numpy at 10 × 10, about an hour) decides the next case.

3. **No bank beats the transitive parse.** The minimum number of exceptions of `P_q` against any
   bounded-arity substitution parse (the join/`C3` grammar of session73, or any fixed bank) is
   `q²/4 − o(q²)`, the same as the feedback number of the transitive parse. The spectral bound of
   THM-4814 (C) gives `q²/8`, the score bound gives `q²/8` on linear parses only; the gap between
   `1/8` and `1/4` is the open quantity.

4. **The order-12 ledger.** The order-11 cycle-index types (a parts run of `gentourng 11` through
   `tourn_types`) would give, with `A000568(12) = 154,108,311,168` and `A002785(12) = 1,492,288`, the
   values `A280772(12)`, `SC_prime(12)`, `c(12)`, `l(12)` by the same identities; the prediction from the
   growth of non-rigid primes (895 at 10, about 12,200 at 11 by the Burnside defect) is that fewer than
   two hundred thousand order-11 primes have a nontrivial automorphism.
