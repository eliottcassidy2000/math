# Sheafification, the plus construction applied twice, and tournament structure

opus-2026-11-07-S85.
- Status: PROVED (THM-4892, audited PASS WITH CORRECTIONS), FINITE-EXACT checks, ANALOGY where typed below.
- No new H-value. Chan-Pak is OPEN.

The owner's directive: merge the prior work on sheafification, and say how Grothendieck's plus construction, which must
be applied twice, relates to tournament structure.

## 0. Outcome

1. **Prior work** (section 1). There is no earlier precise use of sheafification or of the plus construction on
   tournaments.
   - The March 2026 "sheafification" texts are metaphors, and their "HP presheaf" has no restriction maps. HYP-1479 and
     HYP-1511, still marked CONFIRMED in the historical index, are retyped as metaphors (MISTAKE-731).
   - MISTAKE-529 retracted the related "tournament topos" claims.
   - The repository's genuine descent theorems on tournaments are THM-796, THM-801 and THM-818: Čech gluing on the
     fixed-path tiling cube, with Z/2 torsors.
   - The H-spectrum lane has proved gluing inequalities and gluing laws (Sessions 195 and 202) and torsor counting
     (Session 210), but never in sheaf language.
   - THM-512 (the Arrow-Condorcet bridge) is the in-repository precursor of the pair-site reading.
2. **Three tournament sites** on which the plus construction can be computed exactly (THM-4892):
   - **The pair site.** Tournaments are the sheafification of linear orders. One plus step suffices, and cyclic
     triangles are where gluing fails (Condorcet).
   - **The strong-component site.** By the comparison lemma, its sheaves are the presheaves on the poset of strong
     subsets, which is to say the presheaves that are multiplicative over strong components. H is the number of global
     sections of the sheaf of 2-coloured odd-cycle packings (OCF). No bijection with Hamiltonian paths is claimed.
   - **The permutohedron's face poset.** THM-4891's relaxations are the sheafifications for finer topologies: fewer
     2-faces are required, so there are more covering sieves.
3. **Where "twice" lives.**
   - On the strong-component site the second plus step can act only through the empty set, the one object covered by
     the empty sieve. On the face site the vertices play this role; on the pair site the second step never acts.
   - Everything natural (packings, linear orders, tournaments, the swap quotient of packings) has a single empty
     structure, so it sheafifies in one step.
   - The canonical two-step case is global data. A single global bit, such as the converse T -> T^op, becomes one bit
     per strong component after exactly two steps.
   - The resulting group of componentwise converses acts on the sheafified tournaments and preserves H (classical).
4. **The +2 operations.**
   - Padding adds two points to the strong-component site. The gate makes the whole tournament a point.
   - mac-mini's THM-4900 shows that the gate acts on the marked-resolvent invariant y = (s(T) - 1)/2, where
     s(T) = F_T(1), as the square of the continued-fraction map y -> 1 + 1/y.
   - "Plus applied twice" and the "+2" operations share only the word and the number 2. No mechanism or functor
     relates them (verbal ANALOGY; cf. MISTAKE-529: equal numerals do not identify mechanisms).

## 1. Prior work, merged

This section is from a repository-wide search: git grep over all tracked files except agent inboxes and `.scratch`. The
audit extended the search (section 8).

| thread | where | what it is | status |
|---|---|---|---|
| S71 "descent chain" (2026-03-14/15) | `04-computation/beyond_ouroboros_S71u.py` (sheafification section), `symbolic_abstract_S71q.py`, `tau_projective_duality_S71s.py`, `beyond_equation_S71x.py`, `category_deep2_S71o.py`, `new_object.py`; HYP-1479, HYP-1511 | "Tournament theory IS the sheafification of HP counting"; an "HP presheaf" F(S) = {HPs of T[S]} with H = global sections | metaphor. As defined the HP presheaf has no restriction maps: an HP of T restricted to S is a linear order, usually not an HP of T[S]. category_deep2's "sheaf condition verified" test is `if True:`. MISTAKE-529 (2026-08-29) retracted the tournament-topos synthesis. HYP-1479 and HYP-1511 are retyped in the historical index (MISTAKE-731). |
| S71 topos texts and the arc-flip sheaf (2026-03 and 2026-06) | `euler_identity_S71w.py` (PART 7), `functors_S71y.py` (PART 5), `tournament_category_23.py` (item 5) | the presheaf topos Set^{TOURN^op}; "F(U) = {HPs of tournaments in U} defines a sheaf" on the arc-flip poset, whose Čech cohomology is said to obstruct the global HP count | metaphor. The two topos texts use different categories TOURN. euler_identity_S71w's TOURN includes embeddings and quotients, and that text calls the internal logic intuitionistic. functors_S71y calls TOURN a groupoid ("all dualities are invertible") and the topos Boolean. Neither defines a coverage or computes a sheafification. The arc-flip "sheaf" is asserted without restriction maps or a coverage. |
| Arrow-Condorcet bridge | THM-512 | identifies the tiling cube with the social-choice cube, transitive tournaments with rational preferences, and the cyclic triangle with the Condorcet paradox | PROVED. This is the in-repository precursor of THM-4892 (A). |
| Tiling-cube descent (codex S9-S15, July 2026) | THM-796, THM-801, THM-813, THM-818; HYP-6825 | X_n (tournaments containing the path 1 -> ... -> n) = (X_(n-1) x over X_(n-2) of X_(n-1)) x F_2(apex): two endpoint-deleted faces glued over their core, together with a free bit for the one tile (n, 1) that neither face sees. The complement-line tower maps two-to-one onto the fibre product of lower lines (a C_2 endpoint-phase torsor). THM-801 gives three-face Čech descent with Z/2 holonomy: at n = 5, 32 of 64 compatible triples do not glue. | PROVED. This is genuine sheaf-type gluing, about tiling addresses rather than H. |
| Metagraph ultrafilters (2026-06) | HYP-2245 | principal ultrafilters on the tiling cube do not descend to the isomorphism quotient | OPEN |
| LRC sheaf lane (2026-06/07) | HYP-2101, 2127, 2902, 3006, 3234, 3301, 3311, 3493, 3766, ...; THM-2542, 2591, 2657, 2751 (Čech holonomy, coboundary, cocycle, Mayer-Vietoris) | sheaves on Lonely Runner charts | proved or synthesis; no tournament H content |
| H-spectrum gluing, without sheaf vocabulary | session202 `overlap-gluing` (H(T) >= sum over pi of h_X(pi) h_Y(pi)); session195 `pathcover-algebra` (exact cut gluing law; no fixed-dimensional bilinear counting algebra); session165 `fi_interfaces` (path-cover interfaces in the presheaf topos [FI, Set]); session210 (permutohedron face-sink products, the projective carrier and its relative-colour torsor, the gate) | gluing of path data across cuts and overlaps | PROVED (audited) |
| The +2 / x2 grid | THM-4900 (mac-mini, 2026-11-06); session211 `recursion-bridge`; session194 `two-sided-padding` | padding and gate act as n -> n + 2 and doubling as n -> 2n, with Baumslag-Solitar d.a = a^2.d; on y = (s(T) - 1)/2 the gate is the square of the continued-fraction map y -> 1 + 1/y | PROVED |
| Numerology | HYP-1193 (Quillen K-groups "generate" 7, 21), HYP-1231 ("d^2 = 0") | | stale or heuristic. The only Quillen-plus mention (`operads_derived_23.py:207`) has no tournament link. |

External context, for the pair-site reading only:
- Abramsky-Brandenburger, *The sheaf-theoretic structure of non-locality and contextuality* (New J. Phys. 13 (2011)
  113036, arXiv:1102.0264). Local sections that do not glue are contextuality.
- Sargsyan, *Localizing preference aggregation conflicts: a graph-theoretic approach using sheaves* (arXiv:2512.02416).
  - Cellular sheaves on a graph of voters. The stalks are total orders on each voter's visible alternatives, with edge
    stalks on the overlaps.
  - Condorcet-type configurations give empty merged stalks ("local impossibilities"), that is, no global section.
  - The paper uses no Grothendieck topology, sheafification, tournaments or Hamiltonian paths.
- The audit's literature search found no prior statement of (A) on the pair site or of (C). Its nearest finds were
  Abramsky et al., *Contextuality, cohomology and paradox* (arXiv:1502.03097), and arXiv:2601.07283 ("Condorcet's
  paradox as non-orientability"). The content of (A) is folklore; only the J_2 formulation is new.

## 2. The plus construction

For a presheaf F, F^+(S) collects the families of local sections on a covering that agree on overlaps. Two families are
identified when they agree on a finer covering.
- **Step one** always yields a separated presheaf: sections that agree locally become equal. It also adds the families
  that agree on overlaps. A compatible family of F^+-sections still need not come from a single F-family, because
  gluings may fail to exist one level further. Uniqueness already holds, since F^+ is separated.
- **Step two** then produces a sheaf.

The textbook case needing both steps is the constant presheaf on a space:
- step one fixes the value on the empty set;
- step two turns "one value" into "one value per connected component".

## 3. Three tournament sites (THM-4892)

**Pair site (A).** Covers must contain every pair.
- A tournament is exactly a compatible family of orders on pairs, so Tour = L^+ = a(L) for linear orders L.
- The transitive tournaments are the sections that glue to a global order.
- Gluing first fails on cyclic triangles. This is Condorcet's paradox in sheaf form (THM-512).
- Here one step always suffices, for every presheaf. The covering pieces (pairs) and their overlaps (singletons and the
  empty set) have only the trivial cover.
- On the strong-component site, by contrast, the overlap of two components is the empty set, which is covered by the
  empty sieve. That is exactly where a second step can act.

**Strong-component site (B), (C), (D).** For a fixed tournament T, a cover of S must contain every strong component of
T[S].
- The irreducible objects are the strong subsets. They are the points of the sheaf topos, which by the comparison lemma
  is the presheaf topos on the poset of strong subsets.
- Order joins play the role of disjoint unions.
- Sheaves are exactly the presheaves that are multiplicative over strong components.
- **2-coloured odd-cycle packings.** They form a sheaf with exactly H(T) global sections (OCF, THM-002).
  - The multiplicativity H(A => B) = H(A) H(B) is classical (THM-1862). Via OCF it is the cardinality of the sheaf
    condition.
  - Rédei's theorem is the unique fixed point of the global colour swap. This is Grinberg-Stanley's parity consequence
    of OCF.
- **The spectrum.** {H(T)} is the multiplicative closure of the point spectrum (strong tournaments). This is classical
  (THM-1862; THM-1370 (1)), so 7 and 21 are statements about points.
- **The swap quotient of packings.** Its global sections index the cells of Session 210's projective carrier, and it is
  not separated.
  - Its separation defect is exactly the relative-colour torsors: (Z/2)^(m-1) over a packing meeting m components.
  - For strong A and B, one plus step replaces the join law q_A + q_B + 2 q_A q_B by q_A + q_B + q_A q_B.
  - In general |Q^+(A => B)| = |Q^+(A)| |Q^+(B)|. For example A = C3 => C3 and B = C3 give 7, not 9.

**Permutohedron face site (E).**
- Coherent face orientations, prod Tour(blocks), form a sheaf for the topology that requires all 2-faces. That is the
  coarsest of the four topologies considered.
- Sheafifying for finer topologies gives exactly THM-4891's relaxations. The topologies require edges plus squares,
  edges plus hexagons, and edges alone.
- So the gaps 7 and 21 hold for the sheaf O itself, for the coarsest topology. Each finer sheafification has sections
  with 7 and 21 sinks.

## 4. Where "twice" lives (typed)

| "applied twice" | where | type |
|---|---|---|
| a(F) = F^(++). The second step can act only through overlaps of generators that are covered by smaller sieves: the empty set on the strong-component site, the vertices on the face site, nothing on the pair site | THM-4892 (Background, (A1), (B4)) | PROVED |
| one global bit (the converse) becomes one bit per strong component after exactly two steps; the resulting (Z/2)^c acts on the sheafified tournaments and preserves H | THM-4892 (B5) | PROVED (the H-invariance is classical) |
| the swap quotient of packings needs one step, because its value at the empty tournament is a point; its defect is the torsor | THM-4892 (D) | PROVED |
| the gate acts on y = (s(T) - 1)/2 as the square of y -> 1 + 1/y | THM-4900 (mac-mini) | PROVED there; the parallel with (+)(+) is a verbal ANALOGY |
| padding adds two points to the strong-component site; the gate makes the whole tournament a point; HC(Gamma T) = H(T) | THM-4892 (F); Session 210 | PROVED facts; "plus twice = +2" is a verbal ANALOGY, with no functor between them |
| THM-796: each level glues two faces over their core, times a free apex bit; the C_2 torsor arises at the line level | THM-796, THM-801 | PROVED there. As iterated gluing it has one level per added vertex. No site is defined, so it is not a plus construction (ANALOGY at most). |
| OCF's "two colours per odd cycle" | THM-002; THM-4892 (C) | PROVED. The 2 is the colour set of a sheaf, not \|Omega_Set\|; MISTAKE-529 forbids that bridge. |

**Reading.** In tournament structure, the separate-then-glue shape of sheafification appears as two levels:
- the strong components, which are the points: H lives there, and the gaps 7 and 21 are decided there;
- their transitive order, which glues by order join: H multiplies there.

The plus construction needs a second pass only for data attached to the whole tournament rather than to its points.

## 5. Consequences for the spectrum problem

- **The sheaf view localises the problem.** Chan-Pak is a statement about the point spectrum {H(T) : T strong} under
  multiplication.
  - This multiplicative reduction is classical (THM-1862; THM-1370 (1)). The sheaf condition for Pk_2 is its
    packing-level form.
  - THM-4094's two-lane reduction uses more (the finite carry atoms 49, 63, 343) and is not re-derived here.
- **The gate.** It makes every tournament the top point of a larger site: Gamma(T) is strong and HC(Gamma T) = H(T).
  - THM-4892 (F) re-derives Session 210's gate formula from OCF. The gate's new odd cycles are exactly the odd-sized
    path pieces (Session 210's Theorem 2): H(Gamma T) = H(T) + 2 sum over |S| odd of H(T[S]) H(T[V - S]).
- Nothing here produces new H-values.
- The sheafification view does not see why 7 and 21 are missing. It only says where to look: at points.

## 6. Checks and reproduction

Run in this folder, standard library only (Python 3):
- `python sheaf_check.py` gives:
  - (1) the pair site for n <= 5: L^+ has 1, 2, 8, 64, 1024 sections = 2^C(n,2), and the image of L has n! elements;
  - (2) pullback stability of J_ss: 0 violations in 300 random cases;
  - (4) the constant presheaf's two steps: 1,072 subsets checked, 688 of them needing both steps;
  - (3) and (5): OCF and multiplicativity on 1,559 tournaments, and |Q| = (H + 1)/2.
- `python gate_packing_check.py`: the gate identity of (F) on 1,399 tournaments (all with n <= 5, plus random n = 6, 7).
- `python q_torsor_check.py`: the join laws of (D) on 80 order joins of 2-3 strong pieces, and non-separation.

**Limits of these scripts** (the audit's MINOR-14). They check cardinalities only:
- check (1) uses the count formula 2^C(n,2), not an enumeration of matching families;
- the "torsor fibre" tests are identities, true for every input;
- the sheafified size of Q is defined by the product formula, so the sheafified join law is checked against itself;
- transitivity of J_ss and Lemma 0 are not tested.

The structural checks are the audit's (`audit/`, standard library plus numpy for the face site):
- `aud_sites_colimit.py`:
  - the plus construction as the colimit over all covering sieves;
  - the J_2 and J_ss axioms;
  - Lemma 0, (A), (B) and the constant presheaf (n <= 6);
  - the comparison formula on 8,410 (object, presheaf) pairs, 402 of the presheaves non-separated.
- `aud_ocf_q_gate.py`:
  - the explicit Pk_2 sheaf bijection (34,311 pairs);
  - the fibres of Q and the (Z/2)^M action (55,613 fibres);
  - the (D)4 counterexample;
  - the gate formula, Session 210 (5) and Lemma 3;
  - padding and (B5).
- `aud_perm_site.py`:
  - the J_X axioms and the finer/coarser order;
  - O separated and sheaf;
  - the global sections at n = 3, 4, 5;
  - the n = 4 sink spectra and THM-4891's witness tables.
- `aud_gate_points.py`: the points of the sites of T, Sigma(T) and Gamma(T), which grow by 2 and by 2^n + 1.

## 7. Next obligations

- **A path site.** Find a site on which Hamiltonian paths themselves, not packings, form a presheaf with honest
  restriction maps. Session 210 lists an explicit bijection between paths and packings as open.
  - The candidate is the position site of intervals of [n]. There the presheaf of injective T-paths is separated, and
    its sheafification is the T-walks.
  - Inclusion-exclusion over vertex subsets (Karp, Bax) recovers H from walk counts. This is a Möbius inversion between
    a separated presheaf and its sheafification; formalise it as such.
- **Gluing failures as sheaf failures.** Restate Session 202's overlap-gluing inequality and Session 195's
  no-bilinear-algebra theorem as statements about the failure of a path presheaf to be a sheaf.
- **A half-gate.** Is mac-mini's continued-fraction "half-gate" y -> 1 + 1/y realised by a one-vertex operation on
  marked tournaments? If so, the analogy of section 4 would have a concrete half-step.
- **Intermediate face topologies.** THM-4891 leaves open the relaxations between tournaments and single face types,
  such as all hexagons plus some squares. In the language of (E), these are the topologies J between J_{squares,
  hexagons} and J_{hexagons}. Which of their sheafifications exclude 7 and 21?

## 8. Audit

Independent audit (background agent, 2026-11-07): **PASS WITH CORRECTIONS**. The report is `audit/AUDIT.md`.
- **The sheaf theory is correct.** All three topologies are Grothendieck topologies. Lemma 0 holds as a colimit over
  every covering sieve. Every computation was re-done in independent code.
- **MAJOR-1.** (E) had finer and coarser inverted: the relaxations are sheafifications for finer topologies. Fixed in
  THM-4892, this note and both INDEX rows.
- **MAJOR-2.** The sheafified join count in (D)4 is false for non-strong factors (counterexample C3 => C3 with C3: 7,
  not 9). It is restricted to strong factors, and the general product law is added.
- **MAJOR-3.** "Hamiltonian paths as global sections" overclaimed. It now reads "H as the number of global sections",
  and no bijection is claimed.
- **MAJOR-4.** Several parts typed new were known: Rédei via the swap (Grinberg-Stanley), the spectrum closure
  (THM-1862, THM-1370 (1)), the gate cycles (Session 210 Theorem 2), the two-component torsor (Session 210), and the
  Condorcet reading (THM-512). THM-4892's status block now uses the audit's typing.
- **Minors and upgrades.** All 16 minors and 4 upgrades are applied. The main ones:
  - the comparison lemma (points = strong subsets);
  - the exact criterion for the second step;
  - (B5) as an action;
  - the gate site gains 2^n + 1 points;
  - credits and the README wording.
- **Ledger.** The errors and the historical retype of HYP-1479 and HYP-1511 are recorded as MISTAKE-731.
