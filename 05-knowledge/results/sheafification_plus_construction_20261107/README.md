# Sheafification, the plus construction applied twice, and tournament structure

opus-2026-11-07-S85. Status: PROVED (THM-4892), FINITE-EXACT checks, ANALOGY where typed below. No new H-value; Chan-Pak
OPEN.

The owner's directive: merge the prior work on sheafification, and say how Grothendieck's plus construction, which must
be applied twice, relates to tournament structure.

## 0. Outcome

1. **Prior work** (section 1). There is no earlier precise use of sheafification or of the plus construction on
   tournaments.
   - The March 2026 "sheafification" texts are metaphors, and their "HP presheaf" has no restriction maps.
   - MISTAKE-529 retracted the related "tournament topos" claims.
   - The repository's genuine descent theorems on tournaments are THM-796, THM-801 and THM-818 (Čech gluing on the
     fixed-path tiling cube, with Z/2 torsors). The H-spectrum lane has proved gluing inequalities and gluing laws
     (Sessions 195 and 202) and torsor counting (Session 210), but never in sheaf language.
2. **Three tournament sites** on which the plus construction can be computed exactly (THM-4892):
   - **The pair site.** Tournaments are the sheafification of linear orders. One plus step suffices, and cyclic triangles
     are where gluing fails.
   - **The strong-component site.** Its sheaves are exactly the presheaves that are multiplicative over strong
     components. Hamiltonian paths are its global sections, via the sheaf of 2-coloured odd-cycle packings (OCF).
   - **The permutohedron's face poset.** THM-4891's relaxations are the sheafifications for coarser topologies.
3. **Where "twice" lives.** On these sites the second plus step only ever acts on data at the empty tournament.
   - Everything natural (packings, linear orders, tournaments, the projective carrier) has a single empty structure, so
     it sheafifies in one step.
   - The canonical two-step case is global data. A single global bit, such as the converse T -> T^op, becomes one bit
     per strong component after exactly two steps.
   - H is invariant under the resulting group of componentwise converses.
4. **The +2 operations.**
   - Padding adds two points to the strong-component site, and the gate collapses everything into one point.
   - mac-mini's THM-4900 shows the gate is literally the square of the continued-fraction map y -> 1 + 1/y on the
     marked-resolvent coordinate.
   - "Plus applied twice" against "+2" is therefore a precise coincidence of forms, typed ANALOGY. No functor
     identifies them.

## 1. Prior work, merged

This section is from a repository-wide search (git grep over all tracked files except agent inboxes and `.scratch`).

| thread | where | what it is | status |
|---|---|---|---|
| S71 "descent chain" (2026-03-14/15) | `04-computation/beyond_ouroboros_S71u.py` (sheafification section), `symbolic_abstract_S71q.py`, `tau_projective_duality_S71s.py`, `beyond_equation_S71x.py`, `category_deep2_S71o.py`, `new_object.py`; HYP-1479, HYP-1511 | "Tournament theory IS the sheafification of HP counting"; an "HP presheaf" F(S) = {HPs of T[S]} with H = global sections | metaphor. As defined the HP presheaf has no restriction maps: an HP of T restricted to S is a linear order, usually not an HP of T[S]. category_deep2's "sheaf condition verified" test is `if True:`. MISTAKE-529 (2026-08-29) retracted the tournament-topos synthesis. |
| Tiling-cube descent (codex S9-S15, July 2026) | THM-796, THM-801, THM-813, THM-818; HYP-6825 | X_n (tournaments containing the path 1 -> ... -> n) = (X_(n-1) x over X_(n-2) of X_(n-1)) x F_2(apex): two endpoint-deleted faces glued over their core, plus the one tile (n, 1) that neither face sees. The complement-line tower maps two-to-one onto the fibre product of lower lines (a C_2 endpoint-phase torsor). THM-801 gives three-face Čech descent with Z/2 holonomy: at n = 5, 32 of 64 compatible triples do not glue. | PROVED. This is genuine sheaf-type gluing, about tiling addresses rather than H. |
| Metagraph ultrafilters (2026-06) | HYP-2245 | principal ultrafilters on the tiling cube do not descend to the isomorphism quotient | OPEN |
| LRC sheaf lane (2026-06/07) | HYP-2101, 2127, 2902, 3006, 3234, 3301, 3311, 3493, 3766, ...; THM-2542, 2591, 2657, 2751 (Čech holonomy, coboundary, cocycle, Mayer-Vietoris) | sheaves on Lonely Runner charts | proved or synthesis; no tournament H content |
| H-spectrum gluing, without sheaf vocabulary | session202 `overlap-gluing` (H(T) >= sum over pi of h_X(pi) h_Y(pi)); session195 `pathcover-algebra` (exact cut gluing law; no fixed-dimensional bilinear counting algebra); session165 `fi_interfaces` (path-cover interfaces in the presheaf topos [FI, Set]); session210 (permutohedron face-sink products, the projective carrier and its relative-colour torsor, the gate) | gluing of path data across cuts and overlaps | PROVED (audited) |
| The +2 / x2 grid | THM-4900 (mac-mini, 2026-11-06); session211 `recursion-bridge`; session194 `two-sided-padding` | padding and gate act as n -> n + 2 and doubling as n -> 2n, with Baumslag-Solitar d.a = a^2.d; the gate on s(T) is the square of the continued-fraction map | PROVED |
| Numerology | HYP-1193 (Quillen K-groups "generate" 7, 21), HYP-1231 ("d^2 = 0") | | stale or heuristic. The only Quillen-plus mention (`operads_derived_23.py:207`) has no tournament link. |

External context, for the pair-site reading only:
- Abramsky-Brandenburger, *The sheaf-theoretic structure of non-locality and contextuality* (New J. Phys. 13, 2011).
  Local sections that do not glue are contextuality.
- Sargsyan, *Localizing preference aggregation conflicts* (arXiv:2512.02416). Discrete sheaves on pairwise preferences;
  Condorcet cycles are the gluing failures.

Neither work treats Hamiltonian paths or the H-spectrum.

## 2. The plus construction

For a presheaf F, F^+(S) collects the families of local sections on a covering that agree on overlaps, identified when
they agree on a finer covering.
- **Step one** can only make F separated: sections that agree locally become equal. It also adds glued families, but
  these may still fail to glue uniquely one level further.
- **Step two** then produces a sheaf.

The textbook case needing both steps is the constant presheaf on a space:
- step one fixes the value on the empty set;
- step two turns "one value" into "one value per connected component".

## 3. Three tournament sites (THM-4892)

**Pair site (A).** Covers must contain every pair.
- A tournament is exactly a compatible family of orders on pairs, so Tour = L^+ = a(L) for linear orders L.
- The transitive tournaments are the sections that glue to a global order.
- Gluing first fails on cyclic triangles. This is Condorcet's paradox in sheaf form.
- Here one step always suffices, for every presheaf, because the covering pieces (pairs) have no proper covers.

**Strong-component site (B), (C), (D).** For a fixed tournament T, a cover of S must contain every strong component of
T[S].
- The points are the strong subsets. Order joins play the role of disjoint unions.
- Sheaves are exactly the presheaves that are multiplicative over strong components.
- **2-coloured odd-cycle packings.** They form a sheaf with exactly H(T) global sections; this is OCF, THM-002.
  - Multiplicativity H(A => B) = H(A) H(B) is the sheaf condition.
  - Rédei's theorem is the unique fixed point of the global colour swap.
- **The spectrum.** {H(T)} is the multiplicative closure of the point spectrum (strong tournaments). So 7 and 21 are
  statements about points (THM-1370 h-spectrum-omits-7-21-all-n, THM-4094).
- **Session 210's projective carrier** (packings modulo the global swap) is not separated.
  - Its separation defect is exactly the relative-colour torsors.
  - One plus step replaces the join law q_A + q_B + 2 q_A q_B by q_A + q_B + q_A q_B.

**Permutohedron face site (E).**
- Coherent face orientations, prod Tour(blocks), form a sheaf for the topology of all 2-faces.
- Sheafifying for coarser topologies gives exactly THM-4891's relaxations: edges plus squares, edges plus hexagons, and
  edges alone.
- So the gaps 7 and 21 belong to the finest topology. Each coarser sheafification has sections with 7 and 21 sinks.

## 4. Where "twice" lives (typed)

| "applied twice" | where | type |
|---|---|---|
| a(F) = F^(++); on the pair and strong-component sites the second step acts only on data at the empty tournament | THM-4892 (A1), (B4) | PROVED |
| one global bit (the converse) becomes one bit per strong component after exactly two steps; H is invariant under the resulting (Z/2)^c | THM-4892 (B5) | PROVED |
| Session 210's projective carrier needs one step, because its value at the empty tournament is a point; its defect is the torsor | THM-4892 (D) | PROVED |
| the gate is the square of y -> 1 + 1/y on y = (s - 1)/2 | THM-4900 (mac-mini) | PROVED there; the parallel with (+)(+) is ANALOGY |
| padding adds two points to the strong-component site, and the gate collapses everything into one; HC(Gamma T) = H(T) | THM-4892 (F); Session 210 | PROVED facts; "plus twice = +2" is ANALOGY, with no functor between them |
| THM-796: each level glues two faces over their core, plus one torsor bit | THM-796, THM-801 | PROVED there. Read as iterating a plus construction over a non-transitive coverage, one level per step, this needs n - 2 steps, not 2. |
| OCF's "two colours per odd cycle" | THM-002; THM-4892 (C) | PROVED. The 2 is the colour set of a sheaf, not |Omega_Set|; MISTAKE-529 forbids that bridge. |

**Reading.**
- In tournament structure, the separate-then-glue shape of sheafification appears as two levels:
  - the strong components (points, where H lives, and where the gaps 7 and 21 are decided);
  - their transitive order (gluing by order join, where H multiplies).
- The plus construction needs a second pass only for data attached to the whole tournament rather than to its points.

## 5. Consequences for the spectrum problem

- The sheaf view localises the problem: Chan-Pak is a statement about the point spectrum {H(T) : T strong} under
  multiplication. This is THM-4094's strong-lane reduction, now as the sheaf condition for Pk_2.
- The gate turns every tournament into a point: Gamma(T) is strong and HC(Gamma T) = H(T).
  - THM-4892 (F) re-derives Session 210's gate formula from OCF. The gate's new odd cycles are exactly the odd-sized
    path pieces: H(Gamma T) = H(T) + 2 sum over |S| odd of H(T[S]) H(T[V - S]).
- Nothing here produces new H-values.
- The sheafification view does not see why 7 and 21 are missing. It only says where to look: at points.

## 6. Checks and reproduction

Run in this folder, standard library only (Python 3):
- `python sheaf_check.py` gives:
  - (1) the pair site for n <= 5: L^+ has 1, 2, 8, 64, 1024 sections = 2^C(n,2), and the image of L has n! elements;
  - (2) pullback stability of J_ss: 0 violations in 300 random cases;
  - (4) the constant presheaf's two steps: 1,072 subsets checked, 688 of them needing both steps;
  - (3) and (5): OCF and multiplicativity on 1,559 tournaments, |Q| = (H + 1)/2, and the torsor fibres.
- `python gate_packing_check.py`: the gate identity of (F) on 1,399 tournaments (all with n <= 5, plus random n = 6, 7).
- `python q_torsor_check.py`: (D) on 80 order joins of 2-3 strong pieces. It checks the join laws, the 2^(m-1) fibres
  and non-separation.

## 7. Next obligations

- A site on which Hamiltonian paths themselves, not packings, form a presheaf with honest restriction maps. The
  candidate is the position site of intervals of [n].
  - There the presheaf of injective T-paths is separated, and its sheafification is the T-walks.
  - Inclusion-exclusion over vertex subsets (Karp, Bax) recovers H from walk counts. This is a Möbius inversion between
    a separated presheaf and its sheafification; formalise it as such.
- Session 202's overlap-gluing inequality and Session 195's no-bilinear-algebra theorem, restated as statements about
  the failure of a path presheaf to be a sheaf.
- Whether mac-mini's continued-fraction "half-gate" y -> 1 + 1/y is realised by a one-vertex operation on marked
  tournaments. If it is, the analogy of section 4 would have a concrete half-step.

## 8. Audit

Pending at this revision.
