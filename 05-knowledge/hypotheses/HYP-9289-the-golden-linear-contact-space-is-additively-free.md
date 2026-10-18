---
id: HYP-9289
title: "The golden linear contact space is additively free: for every finite source A' + uC' + vC' with A', C' in K = Q(zeta15, sqrt(-11)), the golden contact d = 1 + zeta5 and clock_duality_oct15's eight literal ports, every six-colour port pattern that extends to the unit graph of A' extends to the whole source by an additive colouring c_A(a) + c_1(b) + c_2(c) mod 6; hence six-colour forcing needs edges outside the linear space K + uK + vK, and the natural product coordinate v^2 supplies none"
status: >
  OPEN (CONJECTURE), with a proved reduction and FINITE-EXACT evidence.
  PROVED (elementary): the additive extension lemma, and its reformulation as a cyclic-shift correspondence colouring of
  the golden two-distance graph on A'.
  FINITE-EXACT evidence:
  - all 570 patterns of clock_duality's 5400-vertex receipt extend additively with c_1 = c_2, checked on all 30,856 edges;
  - six larger sources (up to 10,597 first-coordinate points or 2,297,178 mixed edges) have no additively blocked
    pattern;
  - the v^2 census: every unit difference with a nonzero v^2 part is a pure v^2-axis unit.
session: opus-2026-10-17-S41
source: 05-knowledge/results/golden_additive_boundary_20261017.md
scripts:
  - 04-computation/experiments/golden_additive_colourings_20261017.py, golden_vsquare_unit_census_20261017.py, golden_forcing_search_20261017.py (+ .out, .json)
related: [clock_duality_oct15 joint_contact_field_20261017f.md (G5, section 5), shared_golden_contact_20261017g.md, THM-4661, THM-4662, HYP-9247, THM-4614]
---

# HYP-9289 — the golden linear contact space is additively free

## Statement

Let `K = Q(zeta15, sqrt(-11))`, `d = 1 + zeta5`, and let `u`, `v` be clock_duality's golden contact units, so that
`|d + u + v| = 1`. By clock_duality's classification (G5 and section 5 of joint_contact_field_20261017f.md), the unit
differences of `K + uK + vK` are the three axes and the contact family `q(d, 1, 1)`, `|q| = 1`.

**Conjecture.** For all finite `A', C' ⊂ K`, with `A'` containing the eight literal ports and `0 ∈ C'`, every
six-colour port pattern that extends to the unit graph of `A'` extends to `A' + uC' + vC'` additively.

**Weaker form.** The same holds with `c_1 = c_2`, i.e. one colouring of the shared coefficients.

**Consequence.** Every six-colour forcing construction in the golden contact field needs edges outside `K + uK + vK`.
Before building a graph, any new coordinate should be census-tested for cross edges, as was done for `v^2`.

## Evidence

All pre-screens below run in seconds to minutes.

| source | `\|A'\|` | `\|C'\|` | mixed edges | additively blocked patterns |
|---|---|---|---|---|
| clock_duality's receipt (base) | 54 | 10 | 56 | 0 (also verified on the full 5400-vertex graph) |
| `A' = P + M + dM` | 364 | 10 | 2,142 | 0 |
| `A' = A + d·U(C)` | 1,099 | 10 | 5,302 | 0 |
| `A' = P + M + dM`, `C' = C + U(C)` | 364 | 174 | 2,297,178 | 0 |
| `A' = P + M + dM`, `C' = M + M + pentagon` | 364 | 29 | 72,358 | 0 |
| two contact steps | 10,597 | 10 | 87,846 | see the .out file |

**The naive reduction fails.** The golden pentagon's diagonals are `d`-steps, so the two-distance graph contains a
`K_5` on the ports and blocks 465 of the 570 patterns. The coefficient colourings evade this by shifting the
constraints.

## Tests

* **Rigid coefficients.** Coefficient sets whose unit graph forbids homomorphism-like colourings (many distinct shifts
  per direction). These are the most likely counterexamples.
* **Multiple golden contacts** (`d = 1 + zeta5` together with `1 + zeta5^2`). This requires first a classification
  census for the joint linear space.
* **A proof.** Choose the coefficient colourings so that the shifted two-distance colouring of `A'` reduces to the unit
  graph of `A'`.
