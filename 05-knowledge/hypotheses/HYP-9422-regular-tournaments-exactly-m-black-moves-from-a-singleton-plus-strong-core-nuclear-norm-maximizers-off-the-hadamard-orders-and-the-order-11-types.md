---
id: HYP-9422
title: "After THM-4814: (1) every regular tournament on 2m+1 >= 7 vertices is exactly m black native moves from a singleton plus a strong core; (2) the nuclear-norm maximizers at n = 1, 2 (mod 4) have two distinct nonzero singular values and are one-vertex deletions / bordered forms of the odd-order maximizers; (3) the exception count of P_q against the join/C3 tree grammar is q^2/4 - o(q^2), no bank beating the transitive parse asymptotically; (4) order-11 cycle-index types would deliver the order-12 ledger (A280772(12), c(12)) from A000568(12) = 154,108,311,168"
status: >
  OPEN clauses (2), (4), and prescribed-singleton strengthening of (1); Session75 reconciliation 2026-10-27.
  (1) PROVED for every regular tournament on 2m+1 >= 7 vertices: unprescribed singleton plus one strong core has
  exact black native distance m, by the independently audited named-start construction. The n = 5 pure-blue
  exception remains. This does not prescribe which singleton is detached.
  (2) FINITE-EXACT at n = 5, 6, 9 (ratios to the Cauchy-Schwarz bound 0.9789, 0.9636, 0.9842; at n = 9 all 216 maximizers
  share the singular values sqrt7 (x6), sqrt15 (x2)); the proposed description is untested at n = 10.
  (3) PROVED for fixed-b disjoint bounded-arity substitution parses and individual named arc disagreements;
  the matching leading upper constant assumes transitive baselines are available. Variable banks, physical
  port identifications, block-exception counts and native-move distances are not identified with that scope.
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

1. **Exact black distance — PROVED, 2026-10-27.** For every regular tournament on
   `2m+1 >= 7` vertices, the black native distance to a singleton plus one strong
   core is exactly `m`. The [Session75 named-start proof](../results/session75_overnight_black_root_spectral_20261027/session75-black-core/BLACK_CORE_GEODESICS.md)
   and [independent adversarial audit](../results/session75_overnight_black_root_spectral_20261027/session75-black-audit/ADVERSARIAL_REPORT.md)
   supply the all-order construction; THM-4814(D)'s lower bound and finite
   witnesses remain valid. The proof retains one reversal-gauge endpoint and
   two distinct *other named starts*, not merely alternating the eventual
   singleton between native start/end positions. This guards both blackness
   and a terminal with only one singleton. The pure-blue order-five exception
   remains excluded.

   **Still OPEN:** choosing the detached singleton in advance. The construction
   selects it through its first black path and does not prove the prescribed-
   singleton strengthening. The [original filing](../results/session75_overnight_black_root_spectral_20261027/session75-incoming/HYP-9422.md)
   is preserved as dated history, including its complete n7/n9 and sampled n11
   witness evidence. The balanced even-order extension is separate, not part
   of this original hypothesis.

2. **Maximizers off the Hadamard orders.** For `n ≡ 1, 2 (mod 4)` the maximum of `ν` is attained by
   tournaments with exactly two distinct nonzero singular values; at `n = 9` they are `√7` (×6) and
   `√15` (×2) (216 classes). Conjecturally these are vertex deletions or borderings of the order-`n∓1`
   extremals; `n = 10` (9,733,056 classes, numpy at 10 × 10, about an hour) decides the next case.

3. **No bank beats the transitive parse — PROVED in the fixed-bank scope,
   2026-10-27.** For each fixed effective arity bound `b>=2`, the
   [audited all-shape spectral theorem](../results/session75_overnight_black_root_spectral_20261027/session75-spectral/RESULT.md)
   gives, for finite-field Paley order `q`,

       e_min >= max(0,ceil((q(q-1)-C_b q^(3/2) log2(q))/4)).

   Here `C_b` is the explicit entropy constant in that theorem. Counts are
   individual disagreeing unordered named pairs, children are disjoint modules,
   all internal tournaments are retained, and arbitrary parse shapes/labels are
   allowed. If the grammar includes transitive baselines (as join/C3 does),
   averaging total orders gives `e_min <= floor(q(q-1)/4)`. Hence
   `e_min=q^2/4-o(q^2)` for fixed `b`. This closes the earlier gap in the leading
   constant, not a claim of equality with the finite transitive optimum.

   The [independent audit](../results/session75_overnight_black_root_spectral_20261027/session75-spectral-audit/REPORT.md)
   and [incoming reconciliation](../results/session75_overnight_black_root_spectral_20261027/session75-incoming/INTEGRATION_NOTE.md)
   retain the scope: no arbitrary growing bank, undecomposed large module,
   physical identification of shared ports, block-exception currency, or native-
   move distance is inferred. A bank without transitive baselines needs its
   own upper bound/admissible orders. THM-4814(C)'s weaker earlier lower bound
   and the [original filing](../results/session75_overnight_black_root_spectral_20261027/session75-incoming/HYP-9422.md)
   remain historical provenance, not counterexamples to this resolution.

4. **The order-12 ledger.** The order-11 cycle-index types (a parts run of `gentourng 11` through
   `tourn_types`) would give, with `A000568(12) = 154,108,311,168` and `A002785(12) = 1,492,288`, the
   values `A280772(12)`, `SC_prime(12)`, `c(12)`, `l(12)` by the same identities; the prediction from the
   growth of non-rigid primes (895 at 10, about 12,200 at 11 by the Burnside defect) is that fewer than
   two hundred thousand order-11 primes have a nontrivial automorphism.
