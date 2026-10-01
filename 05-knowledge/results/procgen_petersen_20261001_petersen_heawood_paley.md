# The Petersen and Heawood families in Paley coordinates; a tournament form of the linear Conway–Gordon–Sachs theorem; why the Paley parity bridge cannot prove Collatz; and an all-odd tournament on 14 vertices

**Lane:** procgen "petersen" lane, session `collatz-procgen-20260922`, 2026-10-01. Resumed twice (a reboot, then a network
interruption); nothing from the first instance survived.
**Owner's prompt (verbatim part):** "regarding 'the residues mod 7 one way and the non-residues the other way, the length-3
uniqueness, the Fano/Singer facts, and the precise symmetry split. Of Paley's 21 symmetries, the 3 that fix {1, 2, 4} are
the cycle's own dynamics, and the 7 translations have no Collatz counterpart.' consider the family of 7 petersen graphs. try
proving collatz using the paley parity bridge."
**Scope (orchestrator):** the Petersen family taken literally, the Heawood family, oriented Delta-Y and THM-4524's arc
parity, and a Collatz attempt from an angle different from S15's twentieth note (whose snark analysis and Propositions 3–4
are not repeated).
**Code:** `04-computation/experiments/procgen_petersen_20261001_run.py` (runner; stdout only; ends `ALL CHECKS PASSED`),
helpers `procgen_petersen_20261001_lib.py`, C engines `procgen_petersen_20261001_hp.c` (exact Hamiltonian-path and arc
counts), `procgen_petersen_20261001_par.c` (the same mod 2, bitset) and `procgen_petersen_20261001_anti.c` (mod 2 for
anti-circulant tournaments, one array, N ≤ 26).
**Output:** `05-knowledge/results/procgen_petersen_20261001.out` (runner), and
`05-knowledge/results/procgen_petersen_20261001_census26.out` (the N = 26 census, script
`procgen_petersen_20261001_census26.py`, 18 minutes).
**Audit:** independent audit OWED. No git state was changed by this lane; no HYP/THM files were created.

Labels: PROVED, FINITE-EXACT (exhaustive computation), VERIFIED (random or sampled checks), CITED, EXPLAINED COINCIDENCE,
ANALOGY, NUMEROLOGY, REFUTED, OPEN.

## Status

| # | Claim | Label |
|---|---|---|
| 1 | The Delta-Y/Y-Delta closure of K6 has exactly 7 graphs, all with 15 edges: K6, P7, K3,3,1, P8, K4,4−e, P9, P10 = Petersen. The Delta-Y Hasse diagram is a tree with 6 arrows; all members except K3,3,1 descend from K6 by Delta-Y alone. Heawood family: 20 graphs, 14 Delta-Y descendants of K7. K3,3,1,1 family: 58 and 26 | FINITE-EXACT; the counts agree with the literature (CITED) |
| 2 | Heawood graph = Delta-Y of K7 along the 7 translates of {1,2,4} = split (bipartite double) of the Paley tournament P7. Delta-Y along a set S of Fano lines gives 10 distinct members, one per GL(3,2)-orbit of S | PROVED + FINITE-EXACT |
| 3 | In P7 − 0 (underlying graph K6) the 4 translates of {1,2,4} avoiding 0 (t ∈ {0} ∪ QR7) are directed 3-cycles forming a Pasch configuration. Delta-Y on k of them gives K6, P7, P8, P9, P10 (k = 0..4). The 3 translates through 0 (t ∈ NQR7) become the perfect matching {u, 3u} of the Petersen graph. K4,4−e = Delta-Y on the real code {1,2,4} and the 2-adic code {3,5,6} of the trivial cycle. K3,3,1 = Y-Delta at the meeting point of two Pasch lines | PROVED |
| 4 | Petersen = Heawood minus one point-vertex, with the three resulting degree-2 vertices suppressed. Every Petersen-family member except K3,3,1 is such a vertex-deletion shadow of a Heawood-family graph; K3,3,1 is a shadow of the K3,3,1,1 family | PROVED (first part) + FINITE-EXACT |
| 5 | No Petersen-family member has a symmetry of order 7 (Aut orders 720, 36, 72, 72, 8, 12, 120). In the Heawood family only K7 and the Heawood graph do. A bijection "7 translations ↔ 7 family members" is induced by nothing | FINITE-EXACT; the bijection is NUMEROLOGY |
| 6 | Deleting the vertex 0 of P7 kills exactly the 7 translations and keeps exactly the Frobenius: Aut(P7 − 0) = Stab(0) = Stab({1,2,4}) = ⟨x ↦ 2x⟩, the trivial cycle's dynamics. P7 − 0 is the Paley tournament on the six codes of the trivial cycle (three real, three 2-adic; the deleted 0 is the code of the fixed point 0), and it is THM-4524's first all-odd tournament | PROVED; the Collatz reading is an EXPLAINED COINCIDENCE |
| 7 | Of the 28560 Fano colourings of the Petersen graph in Paley coordinates, 24 are Frobenius-equivariant. The 6 that use only 4 lines use exactly the pencil through the deleted point 0 plus the trivial-cycle line {1,2,4}; this is forced. The other 18 use all 7 lines | PROVED + FINITE-EXACT |
| 8 | Arc-HP parity over all 2^15 orientations of each member. Rédei's "H odd for every orientation" holds only for K6 (in general exactly for complete graphs). All-odd orientations exist only for K6 (the 240 copies of P7 − v). No orientation of a Delta-Y move preserves either property | FINITE-EXACT + PROVED; a preserved HP parity is REFUTED |
| 9 | The parity that IS preserved across the family is the Conway–Gordon–Sachs sum: for every embedding, the linking numbers of the disjoint cycle pairs add up to an odd number | CITED (Sachs; Conway–Gordon) + VERIFIED on 180 random PL embeddings of each member |
| 10 | Theorem L (linear K6). Two triangles are linked iff exactly one vertex of one triangle is "uniform-opposite" in the Gale-dual tournament T*; equivalently iff a ±1 height walk crosses level 3. Hence there are 1 or 3 linked pairs, and 3 iff T* ≅ C3[TT2,TT2,TT2], the circular one of the two 6-vertex H-maximizers (H = 45). P7 − v, the other maximizer (all arcs odd), is not circular and is the Gale tournament of no linear K6. Theorem L′: for n + 3 points in R^n (n odd) the linked pairs of complementary simplices have |lk| = 1, and their number is odd and at most r (r odd) or r − 1 (r even), with r = (n+3)/2 | PROVED (the 1-or-3 count and linear CGS are known: Hughes; Huh–Jeon; Bogdanov–Matushkin; we did not find the tournament form, the walk, the general counts or the C3[TT2,TT2,TT2] identification in print) + FINITE-EXACT (n = 3: all 18720 totally cyclic sign patterns, exact arithmetic) + VERIFIED (n = 5, 7, 9, random) |
| 11 | Every Cayley digraph of an odd-order abelian group, with any connection set, has every arc on an even number of directed HPs. So every "code digraph" of a cycle, every dynamics-invariant tournament on an odd cycle, and every Cayley digraph on Z/(2^L − q^k) is arc-even | PROVED (THM-4524 C1's involution, without the tournament hypothesis) |
| 12 | For a Mersenne prime M = 2^L − 1 (L ≥ 3), every unit cycle code induces the same circulant T_L inside QR_M; the only per-cycle data are Legendre signs, which do not separate integral from rational cycles (both signs occur on both sides, for 3x+1 and for 5x+1) | PROVED + FINITE-EXACT |
| 13 | q-interpolation obstruction: the w-cycle of x ↦ x/2, (qx+1)/2 is c_w(q)/(2^L − q^k) with deg c_w ≤ k − 1, so its integrality is not a function of (w, q mod N), for any N. The word 100 is the integral cycle {1,4,2} of 3x+1 and {−1,−4,−2} of 5x+1, and is integral for no other q | PROVED |
| 14 | Collatz via the Paley parity bridge | NO PROOF. The bridge supplies only q-independent (or residue) data, while integrality of cycles is archimedean. Collatz OPEN |
| 15 | NEW (byproduct): QR_127 restricted to μ14 = ±⟨2⟩ is a 14-vertex tournament with every arc on an odd number of HPs (H = 24540117). This settles the first open existence case (N = 14) of HYP-9167(b) | FINITE-EXACT (two independent engines) |
| 16 | Anti-circulant tournaments (a cyclic anti-automorphism through all vertices) exist only for N ≡ 2 (mod 4); in all of them the antipodal arcs are odd. They contain every QR_q − 0. All-odd isomorphism classes among them: 1, 1, 1, 2, 4, 2 for N = 6, 10, 14, 18, 22, 26. QR_p restricted to the sixth (tenth) roots of unity is all-odd iff p ≡ 7 (resp. 3) (mod 8) | PROVED + FINITE-EXACT; proposed conjectures in §6.5 |
| 17 | Snark, forbidden-minor and minimal-counterexample parallels | ANALOGY (no transfer) |

## 0. Plain-language answer to the owner

1. **The family of 7 Petersen graphs is real and verified.** It consists of K6, K3,3,1, the Petersen graph and four
   intermediate graphs, all with 15 edges, linked by Delta-Y moves (a triangle becomes a three-armed star and back). They
   are exactly the minimal graphs that cannot be drawn in space without two linked cycles (every other such graph
   contains one of them as a minor). Their root K6 is the underlying graph
   of the Paley tournament minus a vertex. K7, the underlying graph of the full Paley tournament, roots the 20-member
   Heawood family.
2. **In Paley coordinates everything is built from the 7 translates of {1, 2, 4}.**
   * Doing Delta-Y on all 7 translates turns K7 into the Heawood graph, the Fano incidence graph.
   * Delete the point 0. The 4 translates that avoid 0 (shifts by 0, 1, 2, 4: the residues and 0) are turned into the 4
     star centres of the Petersen graph. The 3 translates through 0 (shifts by 3, 5, 6: the non-residues) become its 3
     remaining "spokes".
   * K4,4 − e comes from the two codes of the trivial cycle: the residues {1, 2, 4} and the non-residues {3, 5, 6}.
3. **Your symmetry split is exactly what deleting a vertex does.**
   * Removing 0 kills all 7 translations and keeps exactly the 3 Frobenius symmetries, which are the trivial cycle's own
     dynamics.
   * The tournament that remains is the Paley tournament on the six codes of the cycle 1 → 4 → 2, and it is precisely the
     first tournament in which every arc lies on an odd number of Hamiltonian paths (THM-4524).
   * But the 7 graphs of the family do not correspond to the 7 translations: no member has any symmetry of order 7. The
     coincidence of the two sevens is numerology. The honest link is the split 7 = 4 + 3 above.
4. **The parity the family really carries is linking, and it is a different parity from Paley's.**
   * In every spatial drawing of any of the 7 graphs, the linking numbers of disjoint cycle pairs add up to an odd number
     (Conway–Gordon–Sachs).
   * For straight-line drawings of K6 we show that this is literally the parity of the crossings of a ±1 walk read off a
     tournament (the Gale dual).
   * The maximal case, 3 linked pairs, corresponds to C3[TT2, TT2, TT2], one of the two 6-vertex tournaments with the
     most Hamiltonian paths (45).
   * The other tournament with 45 paths is Paley minus a vertex. It carries THM-4524's parity but has no place in the
     linking picture. Hamiltonian-path parities are not preserved along the family at all.
5. **Collatz: no proof.** The parity bridge produces only objects that cannot see the multiplier 3:
   * Cayley digraphs of odd groups, which are always arc-even;
   * punctured Paley tournaments, which are the same for every cycle of a given length;
   * Legendre bits, which depend on the parity word alone.

   Whether a cycle consists of integers depends on the size of 2^L − 3^k, which no parity or residue invariant can see. For
   example, the same word 100 gives an integer cycle for 3x+1 ({1, 4, 2}) and for 5x+1 ({−1, −4, −2}), and for no other
   multiplier.
6. **A byproduct.** Testing these objects produced a tournament on 14 vertices in which every arc lies on an odd number of
   Hamiltonian paths. This settles the first open existence case of HYP-9167(b). It belongs to a new family
   (anti-circulant tournaments) that also contains all the Paley-minus-a-vertex examples.

## 1. The Petersen family, verified (FINITE-EXACT; CITED for the theory)

A **Delta-Y move** replaces a triangle by a new degree-3 vertex joined to its three corners. **Y-Delta** is the reverse: a
degree-3 vertex is deleted and its neighbours are joined. The **Petersen family** is the set of graphs reachable from K6
by these moves. It is the list of forbidden minors for linkless embedding (Robertson–Seymour–Thomas 1995; Sachs 1983).

Method:
* Closures were computed with nauty canonical forms and re-checked with networkx isomorphism (two independent methods).
* No Y-Delta move in any of the three families below creates a multiple edge, so the simple-graph convention does not
  matter.

| member | n | degrees | triangles | girth | bipartite | Aut order | undirected HPs | HCs |
|---|---|---|---|---|---|---|---|---|
| K6 | 6 | 5^6 | 20 | 3 | no | 720 | 360 | 60 |
| P7 = Delta-Y(K6) | 7 | 5^3 4^3 3 | 10 | 3 | no | 36 | 324 | 36 |
| K3,3,1 | 7 | 6 4^6 | 9 | 3 | no | 72 | 324 | 36 |
| P8 | 8 | 5 4^4 3^3 | 4 | 3 | no | 8 | 264 | 20 |
| K4,4 − e | 8 | 4^6 3^2 | 0 | 4 | yes | 72 | 324 | 36 |
| P9 | 9 | 4^3 3^6 | 1 | 3 | no | 12 | 192 | 8 |
| P10 = Petersen | 10 | 3^10 | 0 | 5 | no | 120 | 120 | 0 |

**Delta-Y Hasse diagram** (6 arrows, a tree): K6 → P7; P7 → P8; P7 → K4,4−e; K3,3,1 → P8; P8 → P9; P9 → P10.
* K6 and K3,3,1 have no degree-3 vertex, so no Y-Delta move applies to them.
* K4,4−e and P10 are triangle-free.
* The Delta-Y descendants of K6 are the 6 members other than K3,3,1.

**Code checks against the literature.**
* The closure of K7 (the **Heawood family**) has 20 graphs, all with 21 edges. Of these, 14 descend from K7 by Delta-Y
  alone. They are the 14 intrinsically knotted graphs with 21 edges (Kohara–Suzuki 1992; Lee–Kim–Lee–Oh, AGT 15 (2015)
  3305). The other 6 are not intrinsically knotted (Goldberg–Mattman–Naimi 2014; Hanaki–Nikkuni–Taniyama–Yamazaki).
* The closure of K3,3,1,1 has 58 graphs with 22 edges, 26 of them Delta-Y descendants (Goldberg–Mattman–Naimi, AGT 14
  (2014) 1801).

All three counts are reproduced.

## 2. Paley coordinates for both families (PROVED)

Work in Z/7 with D = {1, 2, 4} = QR7 and the 7 lines L_t = D + t (the Fano plane in Singer coordinates, S15 nineteenth
note Prop. 4). P7 is the Paley tournament: x → y iff y − x ∈ D.

**Lemma 2.1.** Each line L_t = {t+1, t+2, t+4} is a directed 3-cycle of P7, namely t+1 → t+2 → t+4 → t+1 (the
differences are 1, 2, 4). Moreover L_t = N^+(t). The 7 lines partition the 21 edges of K7.

**Proposition 2.2 (Heawood).** Delta-Y on all 7 lines of K7 gives the Heawood graph. Identify the centre of L_t with
t_out and the old vertex p with p_in; then the graph is the split (bipartite double) of P7: t_out ~ p_in iff t → p.
For any set S of lines, Delta-Y_S(K7) depends only on the GL(3,2)-orbit of S. The 10 orbits (sizes 1, 7, 21, 7, 28, 28,
7, 21, 7, 1) give 10 pairwise non-isomorphic members, all among the 14 Delta-Y descendants.
*Proof.* The lines are edge-disjoint triangles, so the moves commute. After all of them, every old vertex is joined
exactly to the centres of the three lines through it. Collineations are automorphisms of K7 that permute the lines. The 10
graphs are distinct because their degree sequences differ (`.out`, section B). ∎

**Proposition 2.3 (the Petersen family in Paley coordinates).** Delete the vertex 0. The underlying graph of P7 − 0 is K6.
* **Which lines survive.** 0 ∈ L_t iff t ∈ −D = NQR7 = {3, 5, 6}. So the 4 lines **avoiding** 0 are L_t with
  t ∈ {0} ∪ QR7 = {0, 1, 2, 4}. Any two of them meet in exactly one point (a Pasch configuration). The 3 lines **through**
  0 leave the pairs {1, 3}, {2, 6}, {4, 5} = {u, 3u} (u ∈ D), a perfect matching of K6.
* **The chain.** Delta-Y on k of the 4 Pasch lines gives K6, P7, P8, P9, P10 for k = 0, 1, 2, 3, 4, whichever k lines are
  chosen.
* **P10 is the Petersen graph.** It is cubic on 10 vertices with girth 5, and the Petersen graph is the unique such graph
  (the (3,5)-cage).
* **K4,4 − e** = Delta-Y on D = N^+(0) = {1, 2, 4} and on −D = N^−(0) = {3, 5, 6}. Both are directed 3-cycles of P7 − 0,
  and they are the real and 2-adic codes of the trivial cycle (S15 Prop. 2). The edge missing from K4,4 joins their two
  centres.
* **K3,3,1** = Y-Delta at the common point of two Pasch lines. For Delta-Y_{L_0, L_1}(K6) the common point 2 has neighbours
  y_0, y_1 and its matching partner 6. Y-Delta at 2 gives K3,3,1 with apex 6 and sides {1, 4, y_1} and {3, 5, y_0}.

*Proof.*
* **Edge count.** Fano lines meet in at most one point, so the 4 Pasch triangles are edge-disjoint. Every point p ≠ 0 lies
  on 3 lines, exactly one of which passes through 0. So p lies on 2 Pasch lines (4 of its K6 edges) and keeps one more
  edge, to its partner on the line through 0. Hence Delta-Y on all four gives a cubic graph on 10 vertices.
* **Girth 5.**
  * Two centres share at most one point, so there is no 4-cycle y–a–y'–b.
  * Point–point edges form a perfect matching, so there is no triangle y–a–b (a, b would lie on two common lines) and
    no 4-cycle through two point–point edges.
* **Independence of the choice.** The stabiliser of a point in GL(3,2) is S4, acting as S4 on the 4 lines missing it.
* The remaining identifications are direct (and `.out`, section B). ∎

**Proposition 2.4 (Petersen = Heawood minus a point).** Delete the point-vertex 0_in from the Heawood graph (= split(P7)).
The 3 lines through 0 drop to degree 2. Suppressing them gives the edges {u, 3u}, and what remains is
Delta-Y_{Pasch}(K6) = P10. The same holds for a line-vertex, by duality. Also Delta-Y_{Pasch}(K7) − 0 = P10.

FINITE-EXACT companion: for every Heawood-family graph G and every vertex v, the graph G − v with degree-2 vertices
suppressed was classified. Every Petersen-family member occurs except K3,3,1, which is instead K3,3,1,1 minus an apex
(`.out`, section B).

**Proposition 2.5 (the symmetry split, made exact).** Aut(P7) = {x ↦ ax + b : a ∈ QR7} has order 21.
* **Two stabilisers coincide.** Stab(0) = {x ↦ ax : a ∈ QR7} = ⟨x ↦ 2x⟩ ≅ C3, and Stab(D) is the same group (S15
  Prop. 5), because the bijection t ↦ L_t = N^+(t) between points and lines is Aut-equivariant.
* **The surviving symmetry.** Aut(P7 − 0) has order 3 (THM-4524) and contains Stab(0). So
  **Aut(P7 − 0) = Stab(0) = Stab(D) = the Frobenius.** Deleting 0 removes exactly the translations and keeps exactly the
  multiplier ×2. On D this multiplier is the dynamics of T on the trivial cycle (×4 = ×2^(−1) is T; S15 Prop. 2.4).
* **The six codes.** Under doubling, Z/7 splits into {0} ∪ QR7 ∪ NQR7: the code of the fixed point 0 (word 000), the real
  codes of {1, 4, 2}, and its 2-adic codes. So **P7 − 0 is the Paley tournament on the six codes of the trivial cycle**.
  It is THM-4524's first all-odd tournament: c = 13 on 12 arcs and 23 on the 3 antipodal arcs −u → u.
* **Arc classes.** The Frobenius splits the 15 arcs into five classes of 3:
  * the arcs inside D;
  * the arcs inside −D;
  * u → 3u (the Fano matching);
  * u → 5u;
  * the antipodal arcs −u → u.

  The antipodal arcs (c = 23) lie in Pasch triangles; the matching (c = 13) is what Delta-Y leaves of K6.

The Collatz reading is an EXPLAINED COINCIDENCE, as in S15: period 3 forces Paley codes, and nothing extends.

**Proposition 2.6 (no order-7 symmetry in the Petersen family).**
* The automorphism groups of the 7 members have orders 720, 36, 72, 72, 8, 12 and 120. None is divisible by 7.
* In the Heawood family exactly two orders are divisible by 7: K7 (5040) and the Heawood graph (336). These are the two
  graphs Delta-Y_S(K7) with S translation-invariant (S = ∅ or all lines).

Consequently the 7 translations act on no member of the Petersen family, and a bijection "7 translations ↔ 7 family
members" is not induced by any symmetry or construction: **NUMEROLOGY**. What is true is the 4 + 3 split: relative to the
deleted point, the 7 translates of {1, 2, 4} appear in the Petersen graph as 4 vertices (the Pasch centres) and 3 edges
(the matching).

**Proposition 2.7 (the snark dictionary in Paley coordinates; extends S15 §3).** The Petersen graph of Prop. 2.3 has 28560
Fano colourings (S15's count, reproduced).
* Exactly 24 are equivariant under the Frobenius (edge e ↦ 2e, colour c ↦ 2c).
* 6 of these use only 4 lines, and those 4 lines are always the pencil through the deleted point 0 (L_3, L_5, L_6) plus
  D = L_0.
* The other 18 use all 7 lines.

*Proof that the 4 lines are forced.* The line set used by an equivariant colouring is ⟨×2⟩-invariant. The ⟨×2⟩-orbits on
lines are {L_0}, {L_1, L_2, L_4} and {L_3, L_5, L_6}, so an invariant 4-set is L_0 together with one 3-orbit. Now
L_0 ∪ {L_1, L_2, L_4} is the quadrilateral missing 0, which only yields 3-edge-colourings (S15 Prop. 5(iv)); that is
impossible for a snark. ∎

So the Petersen graph is **built** (by Delta-Y) from the quadrilateral of lines missing 0, and it is **coloured**,
minimally and equivariantly, by the pencil at 0 plus the trivial-cycle line. The two 4-line sets share exactly D.

## 3. Parity across the family

### 3.1 Arc-HP parity over all orientations (FINITE-EXACT)

All 2^15 = 32768 orientations of each member were examined. H is the number of directed Hamiltonian paths; c(e) is the
number through the arc e.

| member | orientations with H odd | all arcs odd | all arcs even | of which H odd | sum of H |
|---|---|---|---|---|---|
| K6 | 32768 (Rédei) | 240 (= P7 − v, H = 45) | 0 | 0 | 737280 |
| P7 | 12960 | 0 | 3752 | 0 | 331776 |
| K3,3,1 | 11808 | 0 | 4592 | 144 | 331776 |
| P8 | 12064 | 0 | 10628 | 0 | 135168 |
| K4,4 − e | 9504 | 0 | 11156 | 0 | 165888 |
| P9 | 7552 | 0 | 19872 | 4 | 49152 |
| P10 | 3968 | 0 | 26480 | 0 | 15360 |

The sum of H equals 2 · #undirected HPs · 2^(15 − (n − 1)).

The C engine agrees with an independent Python enumeration on 280 random orientations.

**Proposition 3.2 (Rédei graphs are the complete graphs; PROVED).** A simple graph G has H(O) odd for every orientation O
iff G is complete.
* If u and v are non-adjacent, orient every edge at u or v away from it. Then both have in-degree 0, and H = 0.
* Complete graphs have the property by Rédei (1934).

So Rédei's parity belongs to the root K6 alone and is never transported by Delta-Y. The same holds for the all-odd
property (only K6, by the table).

**Proposition 3.3 (transit lemma; PROVED).** Let G' = Delta-Y(G) on the triangle abc, with centre y.
* **Cycles.** The undirected Hamiltonian cycles of G' correspond bijectively to the Hamiltonian cycles of G that use
  exactly one triangle edge (a–y–b ↔ ab).
* **Paths.** The Hamiltonian paths of G' are of two kinds:
  * those with y inside, which correspond to the Hamiltonian paths of G using exactly one triangle edge;
  * those with y at an end: for every Hamiltonian path of G using no triangle edge, one path per end lying in {a, b, c}.

So Hamiltonian counts depend on the local structure, and no parity of them is invariant. The lemma was checked on all 44
triangles of the family (`.out`, section D).

**Oriented Delta-Y has no canonical form.**
* A directed 3-cycle has rotational symmetry C3. A C3-invariant orientation of a Y makes its centre a source or a sink,
  which forces the centre to be an end of every Hamiltonian path.
* In Paley coordinates, consider the Frobenius-symmetric orientations of the Petersen graph built from P7 − 0: the matching
  arcs u → 3u are kept, the centre of D is a source or a sink, and the other three centres are oriented compatibly. There
  are 16 such orientations. They have H ∈ {0, 3, 6} (10, 4 and 2 of them) and never have all arcs odd (`.out`, section D).
* So the all-odd property of P7 − 0 survives no orientation of any Delta-Y move; P7 has no all-odd orientation at all.

### 3.4 The preserved parity: Conway–Gordon–Sachs (CITED + VERIFIED)

For every embedding of any of the 7 graphs in R^3, the sum over unordered pairs of vertex-disjoint cycles of the linking
numbers is odd (Sachs 1983; Conway–Gordon 1983). Delta-Y preserves this property (Sachs; Motwani–Raghunathan–Saran 1988).

This was verified on 180 random PL embeddings per member (straight edges, one bend and two bends per edge); the sum was odd
every time. The numbers of disjoint cycle pairs are 10, 9, 9, 8, 9, 7, 6 for K6, P7, K3,3,1, P8, K4,4−e, P9, P10.

Rédei's theorem and Conway–Gordon–Sachs are both degree-mod-2 arguments: an invariant changes by an even amount under a
local move (an arc reversal, a crossing change), and one base case is computed. That is an ANALOGY; §4 makes it an exact
statement for linear embeddings of K6.

## 4. Linear K6 and the Gale tournament (Theorem L)

Let p_1, ..., p_6 ∈ R^3 be in general position (no 4 coplanar), and draw K6 with straight edges.
* Its affine dependences form a 2-dimensional space. Writing a basis as the columns of a 6 × 2 matrix gives **Gale
  vectors** g_1, ..., g_6 ∈ R^2, with Σ g_i = 0 and no two parallel.
* The **Gale tournament** T* has i → j iff det(g_i, g_j) > 0. It is well defined up to reversing every arc.

**Theorem L.**
1. **Which tournaments occur.** T* is a circular (locally transitive), non-transitive tournament, and every such
   6-tournament occurs. There are five isomorphism classes (four up to reversing all arcs), with H = 17, 23 (a chiral
   pair), 41 (RT7 − v) and 45 (C3[TT2, TT2, TT2]).
2. **Piercing rule.** Partition the points into triangles A and B. The number of edges of B that pierce the triangle A
   equals the number of b ∈ B such that either
   * b → every vertex of A and every vertex of B − b → b, or
   * every vertex of A → b and b → every vertex of B − b.

   A and B are linked iff this number is 1.
3. **Walk form.** Sort the 12 rays ±g_i by angle. Let ε_p = +1 if the p-th ray is some g_i and −1 if it is some −g_i (so
   ε_(p+6) = −ε_p), and let h(k) be the number of negative rays among the 6 consecutive rays starting at p = k.
   * The splits {A, B} with conv(A) ∩ conv(B) ≠ ∅ and |A| = 3 are exactly the windows with h(k) = 3.
   * Such a split is linked iff ε_(k−1) = ε_k, i.e. iff the walk h crosses the level 3 at k.
4. **Count.** The number of linked triangle pairs is odd, and it is 1 or 3 (Hughes 2006; Huh–Jeon 2007; the linear
   Conway–Gordon–Sachs theorem). It is **3 iff T* ≅ C3[TT2, TT2, TT2]**, i.e. iff the Gale vectors form three tight
   pairs at 120°.
5. **The two maximizers.** The two 6-vertex tournaments with the maximum H = 45 are:
   * C3[TT2, TT2, TT2]: circular, 9 odd arcs, the maximally linked linear K6;
   * P7 − v: all 15 arcs odd (THM-4524), and **not circular** (the out-neighbourhood of 1 is the cyclic triangle
     {2, 3, 5}). So P7 − v is the Gale tournament of no linear K6.

*Proof.*
* **(2).** The segment b1b2 meets the triangle A iff A | {b1, b2} is the Radon partition of the five points other than b3.
  The affine dependences vanishing at b3 are λ_i = ⟨g_i, u⟩ with u ⊥ g_(b3), i.e. λ_i ∝ det(g_(b3), g_i). So the Radon
  partition of P − b3 is the sign split of i ↦ det(g_(b3), g_i), which is the rule. A closed triangle meets the plane of A
  in 0 or 2 points, with opposite crossing signs; hence lk(A, B) ≠ 0 iff exactly one edge of B pierces A.
* **(3).** A and B are not separated by a plane iff no affine function is positive on A and negative on B. By Gordan's
  theorem this happens iff the switched vectors {g_a : a ∈ A} ∪ {−g_b : b ∈ B} lie in an open half-plane, i.e. iff
  switching T* at A gives a transitive tournament. These switchings are the half-plane windows. In window k the switched
  tournament is transitive in ray order, and b ∈ B is uniform-opposite in T* iff it is the first or the last ray of the
  window. The first ray lies in B iff ε_k = +1, and the last iff ε_(k+5) = −ε_(k−1) = +1. Exactly one of them lies in B
  iff ε_(k−1) = ε_k.
* **(4).**
  * h moves by ±1 and h(k + 6) = 6 − h(k). Start a half period at a k_0 with h(k_0) ≠ 3. The walk ends on the other side
    of 3, so it crosses 3 an odd number of times; touches (ε_(k−1) ≠ ε_k) do not count.
  * h(k) = 3 only at positions of one parity, so at most 3 times per half period.
  * Three crossings force h = 2, 3, 4, 3, 2, 3, 4 (or the mirror pattern), i.e. ε = ++−−++. Then the positive rays sit at
    positions {0, 1, 4, 5, 8, 9}: three tight pairs at 120°, and T* = C3[TT2, TT2, TT2].
* **(1), (5) and an independent confirmation of everything above.** All 6! · 2^5 = 23040 order/sign patterns were
  enumerated. Of these, 18720 are totally cyclic, and every order type of 6 points in general position arises among them.
  The dual points were given exact rational coordinates. Linking numbers, piercing counts and the walk criterion agree in
  every case (`.out`, section F). ∎

**Theorem L′ (all odd dimensions; PROVED + VERIFIED).** Let n be odd, r = (n + 3)/2, and p_1, ..., p_(n+3) ∈ R^n in
general position.
* Define the Gale vectors, the 2(n + 3) rays ±g_i, the windows of n + 3 consecutive rays and the height h(k) exactly as
  above.
* Split the points into two r-sets A and B. The boundary spheres of the simplices conv(A) and conv(B) have dimension r − 2
  each, so their linking number is defined in R^n.
* They are linked iff the split is a window with h(k) = r at which the walk crosses the level r. Every linked pair has
  |lk| = 1.
* The number of linked pairs is odd. It is at most r when r is odd and at most r − 1 when r is even, and every odd value
  up to this bound occurs.
* For n = 3 this gives {1, 3} (Hughes); for n = 5 it gives {1, 3}; for n = 7 it gives {1, 3, 5}.

*Proof.* The argument is the same as for Theorem L.
* A facet conv(B − b) meets conv(A) iff (A, B − b) is the Radon partition of the n + 2 points P − b. This is the
  uniform-opposite rule.
* In a window, only the first and the last ray can be uniform-opposite.
* The linking number equals the signed count of B-facets piercing conv(A), and also, up to sign, the signed count of
  A-facets piercing conv(B). So two piercings of the same kind cancel, |lk| ≤ 1, and the pair is linked iff exactly one
  of the two end rays lies in B.
* Oddness comes from the intermediate-value parity of the walk. The bound comes from the fact that h(k) = r only at
  positions of one parity.
* Every walk that never reaches 0 is realised by a totally cyclic Gale configuration, hence by points; this gives the
  values that occur. ∎

Verified on random configurations (runner section F2): n = 3, 5, 7, 9, with the walk criterion exact in every case. The
observed counts were {1, 3}, {1, 3}, {1, 3, 5} and {1, 3, 5}; the bounds are 3, 3, 5, 5.

The theorem itself (linear CGS, the 1-or-3 count) is classical; Bogdanov–Matushkin (arXiv:1508.03185) give an algebraic
Radon-type proof of the existence statement in every odd dimension. We did not find in print the Gale-tournament form,
the walk form, the exact counts in Theorem L′, or the identification of the 3-link class with C3[TT2, TT2, TT2], but the
literature search was not exhaustive.

**What this says for the owner's question.**
* The parity that defines the Petersen family (odd linking), read at the root K6 = the underlying graph of P7 − v, is
  literally an intermediate-value parity of a ±1 walk.
* Its extremal case is dual to a Hamiltonian-path maximizer, but to the *other* maximizer. P7 − v, the Paley one with all
  arcs odd, has no place in the linear picture.
* So the Paley parity (THM-4524) and the Petersen parity (Conway–Gordon–Sachs) are two different parities on K6, and
  they pick different extremal tournaments.

## 5. Trying to prove Collatz with the Paley parity bridge (a different angle from S15)

S15's twentieth note read the bridge as parity *codes* (Collatz ⟺ 7Φ_T(n) ∈ Z). Here "parity bridge" is read as THM-4524's
arc-HP parity (Paley all-even, Paley minus a vertex all-odd, Rédei), applied to digraphs and tournaments built from cycles.

**Theorem 5.1 (PROVED).** Let Γ be an abelian group of odd order and S ⊆ Γ − {0} arbitrary. In Cay(Γ, S) every arc lies on
an even number of directed Hamiltonian paths.

*Proof.* For an arc u → v put φ(x) = u + v − x. Then φ maps Cay(Γ, S) onto Cay(Γ, −S), the reversed digraph, and swaps
u and v. If P is a Hamiltonian path through u → v at position i, then reverse(φ(P)) is a Hamiltonian path of Cay(Γ, S)
through u → v at position n − i. This is an involution, and a fixed point would need i = n − i, which is impossible for
n odd. ∎

THM-4524 C1 is the tournament case. Checks (`.out` section G): 84 random Cayley digraphs on Z/n (n odd ≤ 13) and on
Z3 × Z3; all code digraphs with L ≤ 4, e.g. QR7 with c = 54 on every arc, and Cay(Z/15, {1, 2, 4, 8}) with H = 239160 and
all 60 arc counts even. The control Cay(Z/6, {1, 2}) (even order) has odd arcs.

**Corollaries (PROVED).**
* For every cycle word w of length L, the code digraph Cay(Z/(2^L − 1), C(w)) is arc-even, whatever the map.
* Every Cayley digraph on Z/(2^L − q^k) (q odd) is arc-even. These moduli, the natural home of the q-aware residues of a
  cycle, are always odd.
* Any tournament on the points of a cycle of odd length L that is invariant under the cycle's own dynamics is a circulant
  on Z/L, hence arc-even. On the code side the dynamics is the multiplier ×2 (the Frobenius); on the index side it is a
  translation; either way Theorem 5.1 makes every dynamics-invariant arc parity even.

**Proposition 5.2 (Legendre bits; PROVED + FINITE-EXACT).** Let M = 2^L − 1 be prime with L ≥ 3, so that 2 ∈ QR_M.
* For a unit code c, the cycle code c⟨2⟩ induces in the Paley tournament QR_M the circulant
  T_L = Cay(Z/L, {d : (2^d − 1 | M) = +1}) when c ∈ QR_M, and its reverse (isomorphic to it) when c ∉ QR_M; this does not
  depend on c.
* T_3 = C3; T_5 ≅ RT5 and T_7 ≅ RT7 (the rotational tournaments, up to a multiplier). For the Mersenne exponents
  L = 5, 7, 13, 17, 19, 31, T_L is not a Paley tournament: for L ≡ 1 (mod 4) none exists, and for L = 7, 19, 31 no
  multiplier maps the connection set to QR_L (for prime order, isomorphic circulants are multiplier-equivalent).
* So the induced code tournament carries no information. The only per-cycle data are the Legendre signs of the real code
  r(w) and of the 2-adic code −c(w). These signs are opposite for every word with L ≤ 7, but not in general: at L = 13 they
  agree for 280 of the 630 necklaces.
* The signs do not separate integral from rational cycles:
  * L = 5: 5x+1's {1, 3, 8, 4, 2} (word 11000) has real sign −1, and 3x+1's {−5, −14, −7, −20, −10} (word 10100) has
    real sign +1; non-integral words of both signs occur.
  * L = 7: 5x+1's three cycles have real signs −1, −1, +1.

**Proposition 5.3 (puncturing).** Deleting 0 from a code digraph kills the translations and keeps the multipliers that fix
C(w).
* At L = 3 this gives P7 − 0, which is all-odd.
* At L = 4 the punctured digraph Cay(Z/15, {1, 2, 4, 8}) − 0 has H = 47208 and **all 52 arcs even**, and
  Cay(Z/15, {3, 6, 9, 12}) is disconnected.

So puncturing is not an all-odd mechanism beyond the Paley case (FINITE-EXACT). For a Mersenne prime M, the punctured Paley
tournament QR_M − 0 (all-odd, conjecturally: HYP-9167(a)) is a single object for all cycles of length L at once, and it
carries nothing cycle-specific.

**Proposition 5.4 (the archimedean obstruction; PROVED).** Consider x ↦ x/2, (qx + 1)/2 and a word w of length L with
k ≥ 1 ones.
* The w-cycle sits at x_w(q) = c_w(q)/(2^L − q^k), where c_w(q) = Σ_(t: w_t = 1) q^(#ones after t) 2^t is a polynomial
  with non-negative coefficients and degree ≤ k − 1.
* So for fixed w and any residue class q ≡ q_0 (mod N), |x_w(q)| → 0, and x_w(q) ∉ Z for all large q.
  **Integrality of a cycle is not a function of (w, q mod N), for any N.** Every invariant built from codes is a function
  of w alone (the case N = 1).

Examples:
* **The word 100 (non-shortcut)** gives x = 1/(4 − q). It is integral exactly for q = 3 (the cycle {1, 4, 2}) and q = 5
  (the cycle {−1, −4, −2}): the same Paley code QR7, with opposite signs.
* **The shortcut word 11000** gives 5x+1's cycle {1, 3, 8, 4, 2} (x = 1) and, for q = 3, the rational cycle through 5/23.

Small-q search (odd q < 200, positive starts ≤ 20000):
* Every q = 2^a − 1 has its cycle through 1.
* Only q = 5 (cycles through 13 and 17) and q = 181 (through 27 and 35) have other positive cycles.
* Within the class q ≡ 3 (mod 4) no other cycle was found.

So the rigorous form of the obstruction is the cycle-level statement above, not "Collatz fails in every residue class of q".

**Switching-class data (PROVED, short).** They do not help either. The switching class of a code tournament is
determined by its code set, so it is a function of the word and Prop. 5.4 applies. And the switching-class sum of H is
the constant 2·N! for every tournament (THM-1467; THM-4524 A2).

**Summary of the objects tried.**

| object built from a cycle word w (length L) | arc-HP parity behaviour | sees q? | sees integrality? |
|---|---|---|---|
| code digraph Cay(Z/(2^L − 1), C(w)) | arc-even (Thm 5.1) | no | no |
| any Cayley digraph on Z/(2^L − q^k) | arc-even (Thm 5.1: the modulus is odd) | through the modulus | no (for an integral cycle every residue is 0) |
| punctured code digraph | all-odd at L = 3; all-even at L = 4 | no | no |
| code subtournament of QR_M | always ≅ T_L, a circulant, so arc-even | no | no |
| two-code tournament Θ(w) (§6.1) | all-odd at L = 3; 35/45 at L = 5; all-odd iff palindromic at L = 7; not all-odd at L = 13 (palindromic) | no | no |
| Legendre signs of the codes | one ±1 per code orbit | no | no |
| switching class of a code tournament | class sum of H = 2·N! | no | no |

**Verdict (Collatz): NO PROOF.** Every parity object the bridge produces falls into one of four kinds:
* a Cayley or circulant digraph, hence arc-even by Theorem 5.1;
* the punctured Paley tournament, which is the same for all cycles of a given length;
* Legendre bits, which are functions of the word;
* the two-code tournaments of §6, which are also functions of the word.

Integral cycles of 3x+1 are determined by the archimedean size of 2^L − 3^k against c_w(3) (Baker theory: Steiner 1977;
Simons–de Weger 2005; Hercher 2023). Divergence is governed by drift and measure (Terras; Tao 2019). The parity bridge
supplies neither.

What would be needed is a parity invariant that sees the size of q (not q mod N) and the sign of the cycle. Prop. 5.4
shows that no invariant of codes and residues can do this.

## 6. Byproduct: an all-odd tournament on 14 vertices, and anti-circulant tournaments

### 6.1 Two-code tournaments

For a cycle word w of length L with M = 2^L − 1 prime, let Θ(w) be the Paley tournament QR_M restricted to the union of the
real code orbit r(w)⟨2⟩ and the 2-adic code orbit −c(w)⟨2⟩. It has 2L vertices. FINITE-EXACT:
* **L = 3:** Θ(100) = P7 − 0, all arcs odd.
* **L = 5:** every Θ(w) has 35 of 45 arcs odd (H = 15505).
* **L = 7:** Θ(w) has **all 91 arcs odd** (H = 24540117) for the 14 palindromic necklaces, and 63 of 91 for the 4 chiral
  ones.

If w is palindromic, then c(w) = r(reverse w) lies in r(w)⟨2⟩, so the vertex set is a coset of the subgroup μ14 = ±⟨2⟩
of F_127^*, and Θ(w) ≅ QR_127[μ14].

### 6.2 HYP-9167(b) at N = 14

**A tournament on N = 14 vertices with every arc on an odd number of Hamiltonian paths exists.** An example is QR_127
restricted to μ14 = {±1, ±2, ±4, ±8, ±16, ±32, ±64}, with x → y iff y − x is a nonzero square mod 127.
* H = 24540117.
* The 91 arc counts take the 7 odd values 3085307, 3125731, 3361361, 3510057, 3604781, 4056367, 4087295.
* |Aut| = 7: the multipliers ⟨2⟩, and nothing else. Scores 6^7 7^7.
* The exact C engine, an independent pure-Python big-integer DP and the bitset parity engine agree.

This settles the first open existence case of HYP-9167(b). HYP-9167 lists N = 14 as open and records that none of the 128
circulants on Z/15 minus a vertex is all-odd.

### 6.3 Anti-circulant tournaments (PROVED framework)

Let m be odd and s: (Z/2m) − {0} → {±1} with s(−d) = −(−1)^d s(d). Define T_s on Z/2m by x → y iff s(y − x) = (−1)^x.
* **The two symmetries.** x ↦ x + 2 is an automorphism of T_s, and x ↦ x + 1 is an anti-automorphism (it maps T_s onto
  its reverse).
* **Characterisation.** Every tournament with an anti-automorphism that is a single cycle through all N vertices is of
  this form, and then N ≡ 2 (mod 4). (If N ≡ 0 mod 4, the N/2-th power would be an automorphism of order 2, and tournaments
  have none.)
* **Examples.** QR_q − 0 (q ≡ 3 mod 4 a prime power; F_q^* is cyclic of order 2m) and QR_p[μ_2m] (p ≡ 3 mod 4,
  2m | p − 1) are anti-circulant: label the vertices by exponents of a generator ζ, a non-residue, and put
  s(d) = χ(ζ^d − 1).

So the all-odd condition N ≡ 1, 2 (mod 4) of THM-4524 and the existence condition N ≡ 2 (mod 4) of anti-circulants
match.

**Theorem 6.1 (PROVED; generalizes the antipodal half of THM-4524 C2).** In every anti-circulant tournament on 2m vertices
(m odd), the m antipodal arcs between x and x + m lie on an odd number of Hamiltonian paths.
*Proof.*
* σ: P ↦ reverse(τ(P)), with τ(x) = x + 1, maps Hamiltonian paths to Hamiltonian paths and arcs x → y to
  y + 1 → x + 1. So ⟨σ⟩ ≅ Z/2m acts on the arcs, and c is constant on its orbits.
* An arc is fixed by σ^m exactly when y − x = m. So the antipodal arcs form one orbit of odd size m, and all other orbits
  have size 2m.
* Σ_e c(e) = (2m − 1)H is odd by Rédei, so m · c(antipodal) is odd. ∎

### 6.4 Census (FINITE-EXACT)

All 2^m sign patterns s were examined, reduced up to multipliers and negation, which give isomorphic tournaments. The
numbers of all-odd isomorphism classes are:

| N | patterns | orbit representatives | isomorphism classes | all-odd classes | identified |
|---|---|---|---|---|---|
| 6 | 8 | 2 | 2 | 1 | P7 − v |
| 10 | 32 | 4 | 4 | 1 | QR11 − v |
| 14 | 128 | 12 | 12 | 1 | QR127[μ14] |
| 18 | 512 | 48 | 44 | 2 | QR19 − v (H = 117266659317), and a class with H = 116670839805 not of the form QR_p[μ18] for p < 2500 |
| 22 | 2048 | 104 | ≤ 104 | 4 | all four are QR_p[μ22]: p = 23 (QR23 − v), p = 199, p = 727, p = 1783 (and more p) |
| 26 | 8192 | 344 | ≤ 344 | 2 | QR27 − 0, and QR_p[μ26] for p = 2003, 2939 (the only all-odd one of the 9 prime-derived classes with p < 3000) |

(The isomorphism-class counts for N ≤ 18 come from a separate exact survey of all patterns in the lane's scratch; the
runner re-checks the all-odd counts for N ≤ 22. The N = 26 row comes from the separate script
`procgen_petersen_20261001_census26.py`, output `procgen_petersen_20261001_census26.out`: 689 checks, 18 minutes. It
uses the one-array engine for anti-circulants, which reconstructs the backward table from the forward one through the
anti-automorphism.)

* Every all-odd class with N ≤ 18 has |Aut| = m: only the shifts by 2.
* Theorem 6.1 was confirmed on every pattern.
* N = 18 and 22 use the bitset parity engine. It agrees with the exact engine on 200 random digraphs, and the N ≤ 14 rows
  use exact counts.

**Prime-derived patterns (data for p < 2500; the criteria for μ6 and μ10 are proved for all p).**
* QR_p[μ6] is all-odd (≅ P7 − v) iff p ≡ 7 (mod 8). *Proof:* with ζ = −ω (ω a primitive cube root of unity),
  s(1) = χ(ω²) = +1, s(2) = χ(1 − ω) and s(3) = χ(−2) = −χ(2). The all-odd orbit of patterns at N = 6 is exactly
  {s(1)s(3) = −1}.
* QR_p[μ10] is all-odd (≅ QR11 − v) iff p ≡ 3 (mod 8). *Proof:* with ζ = −η (η of order 5), and using
  χ(1 − x^(−1)) = −χ(1 − x) for x of odd order, the signs are s = (−ab, −b, ab, a, −χ(2)) with a = χ(1 − η),
  b = χ(1 − η²). For each of the four choices of (a, b) the pattern lies in the all-odd orbit iff χ(2) = −1. The
  primes p < 2500 confirm it (11, 131, 211, 251, ... all-odd; 31, 71, 151, ... not).
* QR_p[μ14] is all-odd for p = 71, 127, 463, 631, 743, ...
* QR_p[μ18] is all-odd for p = 19, 163, 307, 811, ...

At each of N = 6, 10, 14, 18 exactly one of the prime-derived classes (p < 2500) is all-odd; at N = 22 four are;
at N = 26 one of nine is (p < 3000).

### 6.5 Proposed hypotheses (for the orchestrator; no HYP file created)

* **(P1) Existence via anti-circulants.** For every odd m ≥ 1 some anti-circulant tournament on 2m vertices is all-odd.
  This would prove the existence half of HYP-9167(b). Evidence: m = 1, 3, 5, 7, 9, 11, 13 (exhaustive; all-odd class
  counts 1, 1, 1, 2, 4, 2 for m = 3, ..., 13; at m = 13 they are QR27 − 0 and QR2003[μ26]).
* **(P2) Prime-derived form.** For every odd m there are infinitely many primes p ≡ 3 (mod 4) with QR_p[μ_2m] all-odd.
  By Chebotarev, each pattern class that occurs occurs for a positive density of primes, so (P2) for a given m is a finite
  question about Frobenius classes. It is PROVED for m = 3 (exactly the p ≡ 7 mod 8) and m = 5 (exactly the
  p ≡ 3 mod 8), by the sign algebra above together with the finite classification at N = 6, 10.

**Collatz reading of §6: none.** Θ(w) depends on w alone (q-blind). Being all-odd is a property of the word's reversal
symmetry at L = 7, and it holds for integral and non-integral words alike: 5x+1's three 7-cycles and 11 non-integral
palindromic words. The integral 3x+1 cycle {−5, −14, −7, −20, −10} (L = 5) is not all-odd.

## 7. Analogies (typed)

* **Minimal obstructions (ANALOGY).** The Petersen family is the complete finite list of minor-minimal graphs without a
  linkless embedding, and the Petersen graph is the smallest snark. A minimal Collatz counterexample would be the minimum
  of a non-converging component. There is no minor order on Collatz structures, and the finite-unavoidable-set route is
  closed (S15 §3.5). No transfer.
* **Rédei ↔ Conway–Gordon–Sachs (ANALOGY, made exact only for linear K6, Theorem L).**
* **"7 translations ↔ 7 graphs" (NUMEROLOGY).** Prop. 2.6.
* **Other Petersen readings in the repo.** S15's companion note on `origin/main`
  (`odd_zeta_parallels_petersen_20261001.md`) reads the Petersen graph as the boundary graph of M̄_{0,5} (Brown's
  ζ(2) cellular integral). That reading is unrelated to the Fano/Pasch construction here; no claim connects them.

## 8. Open and directions

* **D76.** Is there a Gale-type correspondence for non-linear (e.g. book) embeddings of K6 in which P7 − v, which is not
  circular, plays a role? Theorem L covers straight edges only.
* **D77.** (P1) and (P2) of §6.5. The first cases beyond this note are N = 30 (QR31 − 0, i.e. HYP-9167(a) at q = 31)
  and N = 34 (the anti-circulants on Z/34, where no Paley tournament minus a vertex exists). Both are beyond the 2^N
  subset DP used here.
* **D78.** Prime-derived criteria beyond m = 5: at N = 14, 18, 22 several prime-derived classes occur, so χ(2) alone
  cannot decide; which quadratic characters of cyclotomic units do?
* **D79.** Frobenius-equivariant Fano colourings of the Heawood-family graphs.
* **HYP-9167(b), non-existence side.** N = 13 remains open.

## 9. Reproduction and audit guide

```bash
python3 04-computation/experiments/procgen_petersen_20261001_run.py > 05-knowledge/results/procgen_petersen_20261001.out
```

The run takes about 3 minutes (902151 checks); it needs networkx, nauty (`labelg`, `gentourng`) and gcc. The C engines
are compiled into `scratch/procgen_petersen/build/`.

Suggested audit order (the claims most worth an independent check):
1. **N = 14 all-odd (§6.2).** Build QR_127 on μ14 and count Hamiltonian paths through each arc with any exact method.
   A 14-vertex subset DP in pure Python takes about a second (`exact_hp_python` in the lib is one such implementation).
2. **Theorem L (§4).** The proof is short (Radon/Gale duality plus a switching argument). The exhaustive check is
   runner section F.
3. **Family facts (§1–2).** These are standard and cheap to recheck with any isomorphism tool.
4. **Prop. 5.4 and Theorem 5.1.** One-line proofs.
5. **Census counts at N = 18, 22, 26 (§6.4).** These rely on the two bitset engines, which the runner cross-checks against
   the exact engine.

The N = 26 census is a separate script (about 20 minutes, 270 MB):

```bash
python3 04-computation/experiments/procgen_petersen_20261001_census26.py > 05-knowledge/results/procgen_petersen_20261001_census26.out
```
