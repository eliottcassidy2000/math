---
id: HYP-9435
title: "Open boundaries of the marked-carrier ladder: (1) full necessity of the Euler-characteristic criterion (RESOLVED, now THM-4835 part 4); (2) the Maxwell ladder j >= 2: marking jn - C(j+1,2) edges reduces a tournament by j vertices, but the naive equivariant count is neither necessary nor sufficient (pre-registered prediction FAILED; the failures sit at subgroups of order divisible by 3); (3) 'every non-linear carrier has an involution witness' (REFUTED at n = 10, as predicted); (4) whether counts can coincide where A-sets differ (OPEN; never for n <= 8)."
status: >
  MIXED. Independently audited (audit_20261027), with repairs applied (MISTAKE-689).
  - (1) RESOLVED, PROVED by the S8 referee (audit N1) and adopted into THM-4835 part 4: linear trijectivity holds iff
    chi(G^K) = 1 for all K, iff chi(G^g) = 1 for all elements g.
  - (2) OPEN. The pre-registered naive criterion FAILED:
    - j = 2: 0/10, 0/144, and 6/3817 at n = 4, 5, 6;
    - j = 3: 2/50 at n = 5, and 12/2740 at n = 6.
    The referee recomputed the j = 2, 3 tallies independently. Rigidity of the rung j = 1 is now proved (THM-4835
    part 4).
  - (3) REFUTED by the referee's n = 10 carrier, confirming the session's prediction "FALSE in general". Involution
    witnesses exist for every non-linear carrier with n <= 8. The minimal counterexample therefore has n = 9 or 10;
    n = 9 was not enumerated.
  - (4) OPEN, with FINITE-EXACT evidence: no coincidence for n <= 8 (referee census, 7031 carriers).
session: kind-pasteur-2026-10-27-S8
source: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md
scripts:
  - 04-computation/tournament_carriers/maxwell_ladder.py (+ .out, .json): (2)
  - 04-computation/tournament_carriers/tiers_and_sporadics.py (+ .out, .json): (3), (4)
  - 04-computation/tournament_carriers/audit_20261027/ (necessity_lemmas.py, hyp3_check.py, maxwell.py, p2_anatomy_check.py, census.py, analyse_census.py)
related: [THM-4835, THM-4836, session80_tournament_symmetry_20261027, session81_counterexample_families_20261028, Connelly-Fowler-Guest-Schulze-Whiteley (symmetry-extended Maxwell counts)]
---

# HYP-9435 — Open boundaries of the marked-carrier ladder

## (1) Full necessity — RESOLVED

The referee's lemmas L1–L4 are reproduced in THM-4835 part 4:
- L1, twins: isomorphic components must have `χ = 0`.
- L2, sign: `δ_K >= 0` for every K.
- L3, solvable steps: `φ_K ∈ {-1, 0, 1}` is constant along prime-index steps.
- L4, cyclic to all: fixed sets contain the 2-core.

**Consequences.**
- The criterion is an exact characterisation.
- The cyclic version of the criterion suffices.
- Linear carriers have an A-fixed vertex and no twinned components of nonzero Euler characteristic.

## (2) The Maxwell ladder (OPEN; the naive form is refuted)

- **Setup.** Marking an oriented subgraph with `jn - C(j+1,2)` edges leaves exactly `C(n-j,2)` free pairs. These are
  the edge counts of minimally rigid frameworks in dimension j.
- **Proved.**
  - Sufficiency: A fixes at least j vertices and `|E^K| = j|V^K| - C(j+1,2)` for all K.
  - Necessity at 2-subgroups.
- **Where the naive count fails (audit m2, corrected wording).**
  - Every failing subgroup has order divisible by 3: orders 3, 6, 9 or 18.
  - Each failure has an alternative fixed count: `x ∈ {0, 1}` at an empty fixed set (e.g. a fixed-point-free `Z_3`
    for j = 2, or the whole group `A ≅ Z_6` for j = 3, n = 5), or `x = 3` at `|V^K| = 3`, `|E^K| = 0`.
  - The 3 "count without linearity" cases have no failing subgroup. There A fixes only `2 < j` vertices.
- **Question.** Characterise linear j-jectivity for `j >= 2`. The rung `j = 1` is rigid: THM-4835 part 4.
- **Analogy, typed.** The fixed-subgraph counts parallel the symmetry-extended Maxwell rule of
  Fowler–Guest–Connelly–Schulze–Whiteley (permutation-character version). No implication is claimed.

## (3) Involution witnesses — REFUTED

- **The claim.** For `n <= 8`, every non-linear carrier has an involution defect (1841/1841).
- **The counterexample (referee).** `D = K_4(C_3 + source) + K_1 + K_2 + vec P_3`, with n = 10:
  - `A = Z_3`, with no involutions;
  - `δ(Z_3) = 3`;
  - not linear and not weak (`audit_20261027/hyp3_check.out`).
- **The session's own prediction.** The chiral theta construction (n = 19) also works: `A = Z_3`, `δ = 6`.
- **Open.** Whether the minimal size is 9 or 10.

## (4) Counts versus A-sets (OPEN)

- **Question.** Is there an oriented carrier whose marked count equals `|Tour(X)/A|` for some A-set X of size
  `n-1`, while the A-sets differ?
- **Evidence.** None for `n <= 8` (7031 carriers).
- **Context.** Session80's 81-vertex cycle-index collision is the analogue one level up.
