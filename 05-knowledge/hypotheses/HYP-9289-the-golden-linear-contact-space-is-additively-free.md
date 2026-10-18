---
id: HYP-9289
title: "The golden linear contact space is additively free: for every finite source A' + uC' + vC' with A', C' in K = Q(zeta15, sqrt(-11)), the golden contact d = 1 + zeta5 and clock_duality_oct15's eight literal ports, every six-colour port pattern that extends to the unit graph of A' extends to the whole source by an additive colouring c_A(a) + c_1(b) + c_2(c) mod 6; hence six-colour forcing needs edges outside the linear space K + uK + vK, and the natural product coordinate v^2 supplies none"
status: >
  OPEN (CONJECTURE), with a proved reduction and FINITE-EXACT evidence.
  PROVED (elementary): the additive extension lemma (clock_duality's coordinate-sum method, graded_tail_carrier_20261017b.md),
  and its reformulation as a Z/6 colouring of the golden two-distance graph on A' with forbidden difference sets (a
  correspondence colouring when the sets are singletons).
  FINITE-EXACT evidence:
  - all 570 patterns of clock_duality's 5400-vertex receipt extend additively with c_1 = c_2, checked on all 30,856 edges;
  - four larger linear sources (up to 1,099 first-coordinate points or 2,297,178 mixed edges; contact2 with 10,597 points was built but its pre-screen did not finish) have no additively blocked
    pattern;
  - the v^2 census: every unit difference with a nonzero v^2 part is a pure v^2-axis unit.
session: opus-2026-10-17-S41
source: 05-knowledge/results/golden_additive_boundary_20261017.md
scripts:
  - 04-computation/experiments/golden_additive_colourings_20261017.py, golden_vsquare_unit_census_20261017.py, golden_forcing_search_20261017.py (+ .out, .json)
related: [clock_duality_oct15 marked_contacts_root_current_20261017e.md (G1, section 5), joint_contact_field_20261017f.md (section 5), graded_tail_carrier_20261017b.md (coordinate-sum additive CSPs, prior method), shared_golden_contact_20261017g.md, golden_square_contexts_20261017h.md (concurrent: G union (G + v^2) Cartesian, additive formula (3), degree-192 escape outside the linear space), THM-4661, THM-4662, HYP-9247, THM-4614]
---

# HYP-9289 — the golden linear contact space is additively free

## Statement

Let `K = Q(zeta15, sqrt(-11))`, `d = 1 + zeta5`, and let `u`, `v` be clock_duality's golden contact units, so that
`|d + u + v| = 1`. By clock_duality's classification (lemma G1 of marked_contacts_root_current_20261017e.md section 5, hypotheses certified in joint_contact_field_20261017f.md section 5), the unit
differences of `K + uK + vK` are the three axes and the contact family `q(d, 1, 1)`, `|q| = 1`.

**Conjecture.** For all finite `A', C' ⊂ K`, with `A'` containing the eight literal ports and `0 ∈ C'`, every
six-colour port pattern that extends to the unit graph of `A'` extends to `A' + uC' + vC'` additively.

**Weaker form.** The same holds with `c_1 = c_2`, i.e. one colouring of the shared coefficients.

**Consequence.** Every six-colour forcing construction in the golden contact field needs edges outside `K + uK + vK`.
For example, clock_duality's degree-192 direction `e` of golden_square_contexts_20261017h.md lies outside it. THM-4668
(PROVED) shows `e` is Cartesian on `K + uK + vK + eK`, so adjoining `e` over `K` preserves additive freedom whenever `k >= χ(L)` and the ports lie in one `e`-layer (it adds a
Cartesian factor), and no `e`-edge changes the old coordinates.
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
| two contact steps | 10,597 | 10 | 87,846 | not finished |

**The naive reduction fails.** The golden pentagon's diagonals are `d`-steps, so the two-distance graph contains a
`K_5` on the ports and blocks 465 of the 570 patterns. The coefficient colourings evade this by shifting the
constraints.

**Beyond the linear space.** For the coefficient sets D(A), D(C), D(C) ∪ D(L), the product coordinates `ue` and `ve` add no cross edges (FINITE-EXACT census), and `e` is Cartesian
(THM-4668). The coordinates `uv` and `e^2` add one contact family each:
- `uv`: `q (0, 1, 1, conj d)`, the reflected contact;
- `e^2`: `q (d, 1, 0, conj d)`, which shares the `(A, u)`-step with the `v`-contact.

The coupled source `A + uC + vC + e^2 C` is still additively free (golden_forcing_search4_20261017.py).

**S42 (opus-2026-10-17-S42; golden_forcing_coupled_20261017.py, golden_reflection_forcing_20261017.py, CaDiCaL 1.9.5).**
All three known families were put into one five-coordinate source `A' + uC' + vC' + (uv)C' + e^2 C'`:
- the `v`-contact, the reflected contact and the `e^2` contact;
- the `v`- and `e^2`-contacts share their `(A, u)`-step;
- the reflected and `v`-contacts share their `(u, v)`-step.

Every case below blocks **0 of 570** patterns, with shared and with independent coefficient colourings:

| `A'` | `C'` | mixed patterns (`v` / reflected / `e^2`) |
|---|---|---|
| base | base | 112 / 16 / 84 |
| base | MM (Eisenstein-rigid, 29 points) | 2,072 / 296 / 224 |
| dM | MM | 144,716 / 296 / 3,488 |
| contact1 | MM | 329,340 / 296 / 9,664 |

**Subgraph caveat.** There is no census of differences with both `uv`- and `e^2`-parts, so these graphs are subgraphs
of the unit-distance graphs on the same addresses. A block would still need certification on the full graph.

**Reflection-forcing experiment.** If every proper colouring of the `uv`-layer `W` showed all five nonzero
`conj(d)q`-differences, the reflected contact would force odd reflection (zero shifts), and the `v`-contact would
reduce to plain 6-colouring of the two-distance graph, which blocks 465 patterns.
- Enlarging `W` to `C' + conj(d)({0} ∪ U)` (211 points) and then once more (2,113 points, 460-896 pairs per
  direction) still blocks 0 of 570.
- In the solver's colourings of `W`, every direction keeps at least one difference value missing (set sizes 4-5).

This adds a precise design target to the **Tests** below: a *difference-forcing* layer.

## Tests

* **Rigid coefficients.** Coefficient sets whose unit graph forbids homomorphism-like colourings (many distinct shifts
  per direction). These are the most likely counterexamples.
* **Multiple golden contacts** (`d = 1 + zeta5` together with `1 + zeta5^2`). This requires first a classification
  census for the joint linear space.
* **A proof.** Choose the coefficient colourings so that the shifted two-distance colouring of `A'` reduces to the unit
  graph of `A'`.
