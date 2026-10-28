---
id: HYP-9435
title: "Open boundaries of the marked-carrier ladder: (1) full necessity of the Euler-characteristic criterion at every subgroup; (2) the Maxwell ladder j >= 2 (marking jn - C(j+1,2) edges reduces a tournament by j vertices), where the naive equivariant count is neither necessary nor sufficient; (3) involution witnesses for non-linearity, predicted to fail at larger n through chiral odd-order symmetry; (4) whether counts can ever coincide where A-sets differ."
status: >
  OPEN. FINITE-EXACT evidence; one pre-registered prediction FAILED and is recorded as such. Not yet audited.
  - (1) Part 4 of THM-4835 proves necessity for 2-subgroups and for subgroups of odd prime order. The census agrees at
    every subgroup for all 1250 oriented carriers with n <= 7. Open for subgroups of odd composite order and for
    mixed subgroups whose fixed sub-carrier is not reached by a normal series of prime-order steps. The restriction to
    the full group matters: K4 + 3K1 is linear for Z3 although chi(G^{Z3}) = -2 (THM-4835 (8)).
  - (2) Pre-registered P2 (PREREGISTRATION.md, pushed in efdf0b0a5f before the run) predicted "linear j-jectivity iff
    |E^K| = j|V^K| - C(j+1,2) for all K" for j = 2 (n <= 6) and j = 3 (n = 5, 6). It FAILED:
    - j = 2: 0 disagreements at n = 4, 5; 6 of 3817 at n = 6;
    - j = 3: 2 of 50 at n = 5; 12 of 2740 at n = 6 (9 linear without the count, 3 count without linearity).
    Every disagreement is at an order-3 subgroup with an alternative fixed count x != |V^K| - j, or (3 cases) at a
    group fixing only 2 < j vertices.
  - (3) Every non-linear carrier with n <= 7 has an involution defect (215/215).
  - (4) For n <= 7 no carrier has its marked count equal to some |Tour(X)/A| while its A-set differs; the 5 non-weak
    carriers fail at the count level too.
session: kind-pasteur-2026-10-27-S8
source: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md
scripts:
  - 04-computation/tournament_carriers/maxwell_ladder.py (+ .out, .json): (2)
  - 04-computation/tournament_carriers/tiers_and_sporadics.py (+ .out, .json): (3), (4)
  - 04-computation/tournament_carriers/carrier_census.py: (1)
related: [THM-4835, THM-4836, session80_tournament_symmetry_20261027, session81_counterexample_families_20261028, Connelly-Fowler-Guest-Schulze-Whiteley (symmetry-extended Maxwell counts)]
---

# HYP-9435 — Open boundaries of the marked-carrier ladder

## (1) Full necessity (conjecture)

**Conjecture.** For every oriented excess -1 carrier D with full automorphism group A: if D is linearly trijective,
then `χ(G^K) = 1` for **every** subgroup `K <= A`.

**Proof route.**
- Write `K` as an extension, with a normal subgroup N of prime order inside a subnormal series.
- The sub-carrier `D^N` is linearly trijective for `K/N` by restricting the `A x Z2` isomorphism to N-fixed points.
- The prime-order argument of THM-4835 (4)(ii) needs transpositions of isomorphic components *of `D^N`*. These need
  not be restrictions of automorphisms of D; this is the gap.

## (2) The Maxwell ladder (open problem with a refuted naive form)

- **Setup.** Marking an oriented subgraph with `jn - C(j+1,2)` edges leaves exactly `C(n-j,2)` free pairs. These are
  the edge counts of minimally rigid frameworks in dimension j: trees (j = 1), Laman graphs (j = 2), and so on.
- **What is proved.**
  - Sufficiency holds when A fixes at least j vertices and `|E^K| = j|V^K| - C(j+1,2)` for all K.
  - Necessity holds for 2-subgroups (index-2 chain argument).
- **What fails.** The naive equivalence was refuted by the pre-registered test: odd-order subgroups admit alternative
  fixed counts. For example, a fixed-point-free `Z_3` can act on X with one fixed point when `j = 2`.
- **Question.** Characterise linear j-jectivity for `j >= 2`.
- **Conjectured rigidity of j = 1.** j = 1 is the unique rung where the equivariant count is exactly the criterion.
  This is proved at prime-order subgroups and holds in the census.
- **Analogy, typed.** The fixed-subgraph counts parallel the symmetry-extended Maxwell rule of
  Fowler–Guest–Connelly–Schulze–Whiteley (scalar version: permutation characters only). No implication is claimed.

## (3) Involution witnesses

- **Observed.** For `n <= 7`, every non-linear carrier has an involution defect.
- **Prediction: FALSE in general.** A candidate counterexample: a component containing a theta graph whose three
  paths carry a directed triangle through their midpoints, together with pairwise non-isomorphic, rigidly oriented
  tree components making `χ = 1`.
  - Chirality removes every involution. The `Z_3` rotation then has `δ = 6`, which no X can match, so the carrier is
    non-linear (and non-weak) with no involution defect.
- **Status.** The construction is to be built and tested (n ≈ 20).

## (4) Counts versus A-sets

- **Question.** Is there an oriented carrier whose marked count equals `|Tour(X)/A|` for some A-set X of size `n-1`,
  although the A-sets `{T ⊇ D}` and `Tour(X)` are not isomorphic?
- **Observed.** None for `n <= 7`.
- **Context.** Session80's 81-vertex cycle-index collision is the analogous phenomenon one level up: equal counts,
  different groups.
