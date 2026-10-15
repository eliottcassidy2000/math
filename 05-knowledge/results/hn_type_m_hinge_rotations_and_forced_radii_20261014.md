# Hadwiger–Nelson type-M graphs: hinge rotations, forced radii and mono-set rigidity

Session `mac-mini-2026-10-14-wilson`, 2026-10-14. Exploratory and not independently audited; nothing here is canon.

Typing:
* Graph constructions and unit distances are **FINITE-EXACT** (exact arithmetic in a real multiquadratic field).
* Colourability verdicts are **SOLVER** verdicts: CaDiCaL 1.9.5 through pysat. No DRAT certificates were produced.

Scripts and outputs: `04-computation/experiments/hn_records_20261014/` (see its `README.md`).

## Why this target

Note §3 of `repair_algebra_minus_one_and_wilson_dictionary_20261014.md` ranks a constructive Hadwiger–Nelson record as the most realistic notable result among the owner's four problems. The current records:

* **Smallest 5-chromatic unit-distance graph:** 509 vertices, 2,442 edges (Parts, arXiv:2010.12665).
* **Smallest Moser-spindle-free one:** 1,441 vertices (Heule 2021; corroborated by Haugland, arXiv:2608.04542; the venue "Geombinatorics 31(2)" is not confirmed).

Parts builds his graphs as **type M**, `G = L ∪ ρS`. Here `L` and `S` are Minkowski sums of rotated hexagonal wheels, and `ρ = exp(i arccos(7/8))` rotates `S` about a shared vertex. He writes that "working constructions of type M with other rotations are not known yet". This note tests that question.

## Setting

* **Pieces** (Parts' notation, reproduced exactly):
  * `η = exp((i/2) arccos(5/6)) = (√33 + i√3)/6`, and `H` is the hexagonal wheel.
  * `H^m = ∪_(|α|≤m) η^α H`.
  * `V25 = H¹ ∪ (i/√3)H`, `V31 = H¹ ∪ i√3·(η^(−1)H ∪ ηH)`, `V37 = H¹ ∪ (i/√3)H¹`, `V49 = H² ∪ (i/√3)H¹`.
  * `L727 = V49 ⊕ V37` and `S361 = V31 ⊕ V25`. Sizes 727 and 361 match Parts. S361 needs his second definition of V31.
* **Hinge rotations.**
  * For a squared radius `t`, `w_t = exp(iθ_t)` with `cos θ_t = 1 − 1/(2t)` and `sin θ_t = √(4t − 1)/(2t)`, so `|p − w_t p| = 1` whenever `|p|² = t`.
  * `ρ = w_4`. The rotation field is `Q(√−(4t − 1))`, matching the corpus's rotation-field law (`hn_rotation_field_law_s688.py`).
  * In the field `Q(√−3, √−11)` of the pieces, `w_7/3 = (11 + 5i√3)/14`, `w_3 = η²` and `w_7 = (13 + 3i√3)/14` stay inside; `ρ` adds `√−15`.
* **Forced radius.** `t` is *forced* in a piece `B ∋ 0` if every proper 4-colouring of `B` gives some vertex with `|p|² = t` the colour of `0`.
* **Mono-set.** The *mono-set* of a colouring is the set of such vertices.

## Findings

1. **Reproduction.** `L727 ∪ ρS361` has 1,087 vertices and 5,892 edges and is not 4-colourable (SOLVER, 13 s). Parts' 509-vertex graph splits as an L-subgraph of 374 and an S-subgraph of 136. His search also adds base orbits from a reserve of about 400 candidates, so that graph need not lie inside this one (not checked).
2. **Forced radii** (`forced_radii_scan.out`, 17 pieces; SOLVER).
   * `t = 4` (radius 2, 18-point orbit, rotation `ρ`) is never forced in a single piece.
     * Among the two-fold sums it is forced in all except `V25 ⊕ V25` and `V31 ⊕ V31`.
     * It is also forced in the three-fold sums `H1 ⊕ H1 ⊕ H1`, `H1 ⊕ H1 ⊕ V25` and `H2 ⊕ H1 ⊕ H1`.
   * `t = 4/9` (6 points, rotation field `Q(√−7)`) and `t = 20/9` (12 points, `Q(√−71)`) are forced only in `V49 ⊕ V37` and `V49 ⊕ V49`.
   * No other rational radius is forced. That includes `t = 7/3`, whose rotation stays in the Moser field.
3. **Mono-set rigidity at radius 2** (`mono_sets_*.out`; SOLVER, by blocking-clause enumeration).
   * **L727** admits exactly 19 realizable mono-sets on its 18-point radius-2 orbit: six of size 6, twelve of size 12, and the full orbit. Each minimal one is three antipodal pairs, for example angles `{60°, 136.78°, 163.22°} ∪ (that + 180°)`.
     * So every 4-colouring of L727 contains three collinear monochromatic triples `p, 0, −p`, with `|p| = 2`.
   * **S361** admits 19,916 mono-sets. The smallest are single antipodal pairs at the six non-hexagonal angle classes.
   * No single rotation joins every pair of L727's minimal mono-sets by a radius-2 unit chord (`θ_4 = 28.955°`). So the type-M contradiction is not carried by the reference edges `p ~ ρp` alone; the other cross edges are essential.
4. **Rotation scan** (`rotation_scan_1.out`, `rotation_scan_2.out`). 250 type-M graphs `A ∪ w_t B`, over 25 base pairs and 10 radii (SOLVER, all decided within the conflict budget).
   * Only `t = 4` (Parts' `ρ`) produces 5-chromatic graphs. The radii `t = 4/9, 11/9, 4/3, 5/3, 20/9, 7/3, 3, 5, 7` give 4-colourable graphs for all 25 pairs.
   * This includes `t = 4/9` with both pieces rigid and forcing it, e.g. `(V31 ⊕ V49) ∪ w_4/9 (V49 ⊕ V37)`.
   * At `t = 4`, the graph is 5-chromatic exactly when one piece is `V49 ⊕ V37` or `V49 ⊕ V49` (rigid) and the other forces radius 2. Two pieces that merely force radius 2, such as `(V37 ⊕ V37) ∪ ρ(V31 ⊕ V25)`, give 4-colourable graphs.
5. **Minimisation baseline** (`minimize.py`; SOLVER).
   * Core-guided random deletion on `L727 ∪ ρS361` gives near-critical subgraphs of 599–617 vertices, about 450 from L727 and 150 from `ρS361`. That is far above Parts' 509, whose split is 374 + 136.
   * After a fix (the symmetry-breaking triangle's vertices are retested under another fixed triangle), four seeds give genuinely vertex-critical subgraphs of 614, 611, 617 and 618 vertices, each with 448–462 from L727 and 156–166 from `ρS361` (`minimize_L727_S361_rho.log`, `crit_*.json` with exact coordinates).
   * Parts' 509 splits as 374 + 136, so greedy deletion leaves about 75 surplus vertices on the L side and 25 on the S side. The union of the four critical sets has 784 vertices.
6. **Exact smallest-MUS search does not scale.** pysat's OptUx (implicit hitting sets) found the smallest non-3-colourable subgraph of V25 (9 vertices) instantly. On `H1 ⊕ H1` it ran more than 10 minutes without finishing. Parts' 509 came from a structured hitting-set search with symmetry and orbit compression, not from a generic smallest-MUS solver.

7. **The spindle-free graphs are type M with Parts' hinge** (`spindle_free_*.log`; SOLVER).
   * In the rebuild, `ψ = φ_0^(−1) ρ`. Its 120 new edges all have angular gap `θ_4 = arccos(7/8) = 28.955°`. Only 72 of them join radius-2 points `z` and `ψw` (for example `2φ_0^(a−1) ~ 2ψφ_0^a`). The other 48 join `|z|² = 2/3` to `|ψw|² = 8/3`, 24 in each direction; the distance works because `2/3 + 8/3 − 2·(4/3)·(7/8) = 1`. In `M_3 ∪ ψP_1` the 48 cross edges split as 24 at radius 2 and 24 between `8/3` and `2/3` (audit B, `hn_cross_edges.out`). So case 10 is a type-M graph with Parts' ρ-hinge, in the field `Q(i, √2, √3, √5)`, but the hinge is not carried by radius 2 alone: the hinge principle must account for the `2/3`–`8/3` edges too. (The pre-audit text said all 120 join radius-2 points, MISTAKE-610.)
   * **Forcing.**
     * `M_1 ⊕ M_1` (2,593 vertices) does not force radius 2.
     * `P_1 = M_2 ⊕ M_1^clock` (15,865 vertices; `M_1^clock = {0} ∪ {φ_0^a}`) forces it on its 24-point orbit `2φ_0^a`.
     * `M_3` (32,257 vertices) forces it on its 72-point orbit.
   * **4-colourable reductions.** `(M_1 ⊕ M_1) ∪ ψ(M_1 ⊕ M_1)` (5,185), `M_3 ∪ ψ(M_1 ⊕ M_1)` (34,849), `P_1 ∪ ψP_1` (31,729) and `(M_2^c ⊕ M_1) ∪ ψ(M_2^c ⊕ M_1)` (11,569) are all 4-colourable.
   * **Independent confirmation (audit B, kissat 4.0.4 on its own exact rebuild in `Q(ζ₂₄)(√5)`):** `M_3 ∪ ψP_1` (48,121 vertices, 385,296 edges) is not 4-colourable (8,786 s); `P_1` forces radius 2 on its 24-point orbit (427 s); `M_3 ∪ ψ(M_1 ⊕ M_1)` and `P_1 ∪ ψP_1` are 4-colourable, with explicit colourings checked edge by edge (`audit_B_20261014/hn_*_kissat404.out`).
   * As with L727, a non-forcing piece kills 5-chromaticity, and two forcing pieces are not enough. The rigidity of `M_3` was not measured. Only the mixed pair `M_3 ∪ ψP_1` was run (`spindle_free_mixed.log`). `P_1 ∪ ψM_3` is congruent to it under complex conjugation, since `M_3` and `P_1` are both conjugation-invariant (checked exactly by audit B), so it needs no separate test.

8. **Later experiments (SOLVER; added after the audit-B scope, not audited).**
   * **The hybrid fails.** Parts' rigid `L727` hinged by `ρ` to the clock-family forcing piece `P_1` (16,591 vertices, field `Q(√2, √3, √5, √11)`) is **4-colourable**. So rigid plus forcing is not sufficient across families.
     * The ρ-edges only join radius-2 points with equal angles, or angles differing by `2θ_4 = 57.91°`.
     * `L727` and `P_1` share only the six hexagonal directions, while `L727`'s minimal mono-sets use the η-directions `60k ± 16.78°`.
     * Even within `L727 ∪ ρS361`, colour-0 sets alone do not force the contradiction: S361's minimal mono-sets (antipodal η-pairs) miss `L727`'s along both angle offsets. The obstruction lives in the full four-colour interface.
   * **Small clock pieces do not force.** `clock ⊕ clock`, `clip_1(clock ⊕ clock) ⊕ clock` and `(clock ⊕ clock ⊕ clock) ∩ {r ≤ 2}` (289–1,273 vertices) do not force radius 2. `P_1` (15,865) does.
   * **Exact smallest MUS is out of reach.** OptUx on the 784-vertex union of the four critical sets did not finish in 45 minutes.
   * **Symmetric minimisation loses here.**
     * `L727 ∪ ρS361` is C6-invariant: one centre plus 181 orbits of size 6.
     * Orbit-level deletion gave an orbit-critical subgraph of 116 orbits (691 vertices). Vertex-level refinement brought it to a vertex-critical **683**, worse than asymmetric greedy deletion (611). This is Parts' "overcompression".
   * **Clipped forcing pieces (SOLVER, CaDiCaL 1.9.5; `spindle_free_batch.log`).** Replace `P_1` by `P' = {p ∈ P_1 : |p|² <= r} ∪ {radius-2 orbit}`. `P'` still forces radius 2, and `M_3 ∪ ψP'` stays not 4-colourable:
     * `r = 3.5`: `|P'| = 15,289`, graph 47,545 vertices / 374,880 edges (4,568 s);
     * `r = 3.0`: `|P'| = 14,329`, graph 46,585 vertices / 364,272 edges (1,776 s);
     * `r = 2.5`: `|P'| = 13,153`, graph **45,409 vertices / 351,816 edges** (3,884 s).
     * Replacing `M_2` by `clip(M_1 ⊕ M_1, 0.8)` loses the radius-2 orbit altogether, so it cannot force radius 2; that branch is dead.
     * **Clipping the rigid side** (`spindle_free_m3side.py`, `spindle_free_m3side.log`; SOLVER, mac-mini-2026-10-14-tourn). `M_3' = {p ∈ M_3 : |p|² <= r} ∪ {radius-2 orbit}` against `ψP'` with `P'` the `r = 3.0` piece (14,329). At `r = 3.5` the graph has **45,241 vertices / 340,584 edges and is not 4-colourable** (CaDiCaL 1.9.5, 3,922 s), the smallest so far in this family. The `r = 3.0` and `r = 2.7` runs were killed by a machine reboot and restarted; see the log. These are steps of a few percent, far from Heule's 1,441.
   * **The spindle-free graph is C24-invariant.** `M_3 ∪ ψP_1` is one centre plus 2,005 orbits of size 24, so symmetric subgraphs have `1 + 24m` vertices. `1441 = 1 + 24·60` is NUMEROLOGY unless Heule's graph is itself C24-symmetric, which was not checked.
     * Orbit-level minimisation (`orbit_core.py`) runs slowly: CaDiCaL with 2,006 assumption literals had not finished its first core after more than an hour, against 25 minutes for the plain solve.
     * The next step is phase-1 cores without assumptions, e.g. by fixing a random half of the orbits off and testing.

## Reading (HEURISTIC)

* **Hinge principle.** In this family a type-M graph is 5-chromatic when one piece is rigid at the hinge radius and the other forces it. Rigid here means a tiny mono-set family: 19 patterns, against about 20,000 for the flexible piece.
* **Forcing is not enough.** At `t = 4/9`, both pieces are rigid in the forcing sense and the graph is still 4-colourable. The hinge angle must also line up with the cross edges that the rigid piece's orbits create.
* **Why ρ.** Its angle `θ_4 = 28.955°` meets the hexagonal and `η`-orbits (angles `60k` and `60k ± 16.78°`) of both pieces in auxiliary unit chords. That is Parts' auxiliary orbit `(0, 0, 2, 6)` at distances `√11/2 ± √3/6`.
* **The geometric meaning of rigidity.** Every 4-colouring of L727 has three collinear monochromatic triples `(p, 0, −p)` with `|p| = 2`. These are the "mono-pairs of length 4" behind de Grey's and Parts' type-A arguments.

## Plan for a record attempt (weeks, not this session)

1. **Reproduce Parts' search.**
   * Hyperedges: minimal vertex sets whose removal makes `L727 ∪ ρS361` 4-colourable, searched by orbit under the graph's symmetry, with reserve orbits added as in his expansion step.
   * Then minimal hitting sets in increasing size.
   * Validate by recovering a 509-vertex graph. Parts estimates 100 laptop hours from scratch.
2. **Search beyond Parts.**
   * His M6B and M6C subtypes, which he says were less explored.
   * Pieces with other rigid orbits, for example sums containing V49 twice.
   * Second hinges: `L ∪ ρS ∪ ρ'S'`, or hinges at a non-central shared vertex.
3. **Spindle-free variant.**
   * Rebuild the Voronov–Neopryatnaya–Dergachev graphs in `Q(i, √2, √3, √5)`, which have no `√11` and so no Moser spindle. Heule reduced these from 64,513 to 1,441 vertices.
   * Then apply the same minimiser. This record is older and less contested.
   * **Recipe** (arXiv:2106.11824, series 2, case 10, read from its HTML; field primes `2, 3, 5`):
     * `φ_0 = e^(iπ/12) = (√6+√2)/4 + i(√6−√2)/4` and `φ_1 = √6/3 + i√3/3`;
     * `M_1 = {0} ∪ {φ_0^a φ_1^b : 0 ≤ a ≤ 23, b ∈ {−1, 0, 1}}`, which has 73 points;
     * `M_2 = clip(M_1 ⊕ M_1, 1)`, keeping `|z| ≤ 1`, then `M_3 = M_2 ⊕ M_1`, which has 32,257 points;
     * `W = M_3 ∪ ψM_3` with `ψ = φ_0^(−1)(7/8 + i√15/8)`, giving 64,513 vertices.
   * The paper reports Glucose 4 needing about 1,400 s to show such a graph is not 4-colourable. Its Proposition 6 proves the graphs are spindle-free.
   * **Done this session** (`spindle_free_build.py`).
     * The rebuild is exact: `|M_2| = 865`, `|M_3| = 32,257`, `|W| = 64,513`.
     * `|E| = 542,472`, which is the paper's 542,352 common edges plus the 120 new edges of case 10.
     * CaDiCaL 1.9.5 confirms it is not 4-colourable (1,974 s; SOLVER).
     * The next step is minimisation towards Heule's 1,441 (`minimize_npz.py`, `spindle_free_pieces.py`).
4. **Certificates.** Any candidate record needs a DRAT proof of non-4-colourability (CaDiCaL can emit one, checked by `drat-trim`), exact coordinates and an independent re-check of all unit distances.
