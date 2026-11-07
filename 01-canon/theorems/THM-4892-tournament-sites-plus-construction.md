---
id: THM-4892
title: "Tournament sites and the plus construction. (A) On the pair site, where a covering sieve must contain every pair, linear orders form a separated presheaf. Its sheafification, reached in one plus step, is the presheaf of tournaments; cyclic triangles are the minimal loci where local orders fail to glue. (B) For a fixed tournament T, the strong-component site, where a covering sieve must contain every strong component of T[S], is a Grothendieck topology whose points are the strongly connected subsets. Its sheaves are exactly the presheaves that are trivial at the empty set and multiplicative over strong components. F^+(S) is the fibre product of the F(C_i) over F(empty), so the second plus step acts only on data that is nontrivial at the empty tournament. Global data, such as the converse bit, sheafify in exactly two steps to one bit per strong component. (C) The 2-coloured odd-cycle packings form a sheaf with exactly H(T) global sections (OCF). So the multiplicativity of H is the sheaf condition, Rédei's theorem is the unique fixed point of the colour swap, and the H-spectrum is the multiplicative closure of the strong (point) spectrum. (D) Session 210's projective carrier is not separated, and its separation defect is exactly the relative-colour torsors. (E) THM-4891's relaxations of parallel coherence are the sheafifications of coherent face orientations for coarser topologies on the permutohedron's face poset. (F) Padding adds two points to the site and the gate collapses everything to one; OCF gives a second derivation of Session 210's gate formula"
status: >
  PROVED (new formulation; elementary): (A)-(F).
  INHERITED and credited:
  - Grothendieck's plus construction. F^+ is separated; if F is separated then F^+ is a sheaf; a = (+)(+).
    (Mac Lane-Moerdijk, Sheaves in Geometry and Logic, III.5; SGA4 II.3.)
  - The odd-cycle formula H = I(Omega(T), 2) (Grinberg-Stanley; repository THM-002).
  - Session 210's gate formula H(Gamma T) = H(T) + P_T(2), its Lemma 3 (F_T(x)F_T(-x) = 1), the projective carrier and
    its join law, and the coherent permutohedron (session210/*/PROOF.md).
  - THM-4891, THM-1370 (h-spectrum-omits-7-21-all-n), THM-4094.
  FINITE-EXACT checks (`sheafification_plus_construction_20261107/*.out`):
  - (A) for n <= 5;
  - pullback stability of the strong-component topology on 300 random cases;
  - the constant presheaf's two steps on all subsets of 40 random tournaments (688 subsets need both steps);
  - OCF and multiplicativity on 1,559 tournaments;
  - (D) on 1,559 tournaments and on 80 order joins of 2-3 strong pieces;
  - the (F) gate identity on 1,399 tournaments.
  NOT CLAIMED:
  - any topos, Lawvere or "2 = |Omega_Set|" bridge (cf. MISTAKE-529);
  - any functor identifying the gate with the plus construction - the "plus twice versus +2" correspondences in the
    source note are typed as analogies there;
  - any new H-value or progress on Chan-Pak (OPEN).
session: opus-2026-11-07-S85
source: 05-knowledge/results/sheafification_plus_construction_20261107/README.md
related:
  - THM-002
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-4094
  - THM-4891
  - THM-4900
  - THM-796
  - THM-801
  - session210/coherent-permutohedron/PROOF.md
  - session210/graded-spectrum-gates/PROOF.md
  - session210/projective-cycle-carrier/PROOF.md
  - session210/structure-atlas/STRUCTURE_ATLAS.md
verification:
  - 05-knowledge/results/sheafification_plus_construction_20261107/sheaf_check.py
  - 05-knowledge/results/sheafification_plus_construction_20261107/gate_packing_check.py
  - 05-knowledge/results/sheafification_plus_construction_20261107/q_torsor_check.py
---

# THM-4892 - tournament sites and the plus construction

## Background (classical)

For a presheaf F of sets on a site, the plus construction is

    F^+(S) = colim over covering sieves R of S of Match(R, F),

where Match(R, F) is the set of matching families. Grothendieck's facts:
- F^+ is always separated.
- If F is separated, F^+ is a sheaf.
- Hence F^(++) = a(F) is the sheafification.

So two applications always suffice. The first can only make F separated; the second then glues.

**Lemma 0.** On each site below, every object S has a smallest covering sieve R_min(S). Hence
F^+(S) = Match(R_min(S), F).

## Statement

**(A) The pair site.** Let J_2 be the topology on the poset of finite subsets of a set V in which a sieve R on S
covers iff R contains every subset of S with at most two elements.
1. **One plus step.** For every presheaf F, F^+(S) = lim over U subset S with |U| <= 2 of F(U). F^+ is already a
   sheaf, so a(F) = F^+.
2. **Linear orders.** Let L(S) be the set of linear orders of S, with restriction to the induced order.
   - L is separated, and not a sheaf once |S| >= 3.
   - L^+ = Tour, the presheaf of tournaments, with restriction to the induced subtournament.
   - L -> a(L) = Tour is the inclusion of the transitive tournaments.
3. **Where gluing fails.** A section T of Tour over S lies in the image of L iff every triple of S does. So the minimal
   subsets on which local linear orders fail to glue are exactly the cyclic triangles of T.

**(B) The strong-component site.** Fix a tournament T on V. Let J_ss be the topology on the poset of subsets of V in
which a sieve R on S covers iff R contains every strong component of T[S].
1. **Topology and points.** J_ss is a Grothendieck topology. Its points (objects whose only covering sieve is the
   maximal one) are the nonempty S with T[S] strongly connected. The empty set is covered by the empty sieve.
2. **Plus formula.** Let C_1, ..., C_c be the strong components of T[S]. Then

       F^+(S) = F(C_1) x_F(empty) ... x_F(empty) F(C_c),   and F^+(empty) = pt.

3. **Sheaves.** F is a sheaf iff F(empty) = pt and F(S) -> prod_i F(C_i) is a bijection for every S. So J_ss-sheaves
   are exactly the presheaves that are multiplicative over strong components.
4. **One step or two.**
   - If F(empty) = pt, then a(F) = F^+.
   - The second step changes F^+(S) only when two components of T[S] restrict to different points of F(empty).
   - For a constant presheaf A (A(S) = A for every S, the empty set included):
     - A^+(S) = A, the diagonal, for nonempty S;
     - A^(++)(S) = A^(c(T[S])), where c counts strong components;
     - so exactly two steps are needed wherever c(T[S]) >= 2 and |A| >= 2.
5. **Global symmetry becomes local symmetry.** The constant presheaf of groups Z/2, generated by the converse
   T -> T^op, sheafifies in exactly two steps to (Z/2)^(c(T)). This is the group of componentwise converses (reverse
   the arcs inside any set of strong components), and H is invariant under all of it.

**(C) Hamiltonian paths as global sections.** Let Pk_2(S) be the set of pairs (p, col): p is a set of pairwise
vertex-disjoint directed odd cycles of T[S], and col: p -> {0, 1}. Restriction to S' subset S keeps the cycles contained
in S'.
1. **H counts global sections.** Pk_2 is a J_ss-sheaf, and |Pk_2(V)| = H(T) (OCF).
   - The law H(A => B) = H(A) H(B) for order joins is its sheaf condition for the cover {A, B}.
2. **Rédei's theorem.** The global colour swap is an involution of Pk_2(V) with exactly one fixed point, the empty
   packing. Hence H(T) is odd.
3. **Spectrum.** {H(T)} is the multiplicative monoid generated by the point spectrum {H(T) : T strong}. The exclusions
   7 and 21 (THM-1370, h-spectrum-omits-7-21-all-n) and THM-4094's strong-lane reduction are statements about points.

**(D) The projective carrier is not separated.** Let Q = Pk_2/(global swap), the presheaf quotient.
- Q(V) is Session 210's projective carrier together with one basepoint, the empty packing, so |Q(S)| = (H(T[S]) + 1)/2.
1. **Not separated.** Q is not J_ss-separated as soon as two strong components of T[S] carry nonempty packings.
2. **One plus step.** Since Q(empty) = pt, a(Q)(S) = Q^+(S) = prod_i Q(C_i), and |a(Q)(S)| = prod_i (H(C_i) + 1)/2.
3. **The defect is the torsor.** Q(S) -> Q^+(S) is surjective. The fibre over a point whose packing meets m components
   is a torsor under (Z/2)^m / diagonal = (Z/2)^(m-1). These are exactly Session 210's relative-colour torsors.
4. **Counting.** Put q = (H - 1)/2. For an order join A => B, the presheaf count is |Q| - 1 = q_A + q_B + 2 q_A q_B
   (Session 210's join law), and the sheafified count is |Q^+| - 1 = q_A + q_B + q_A q_B. The factor 2 on q_A q_B is
   the separation defect.

**(E) THM-4891 as a table of sheafifications.** On the face poset of the permutohedron (ordered set partitions, ordered
by inclusion of faces), let O(B_1 | ... | B_k) = prod_j Tour(B_j) be the presheaf of coherent orientations of a face.
For X a set of 2-face types, let J_X be the topology in which a sieve on a face F covers iff it contains every edge of F
and every 2-face of F of a type in X.
- **The finest topology.** For X = {squares, hexagons}, O is a sheaf. Its global sections are the tournaments
  (THM-4891 (i)).
- **The coarser topologies.** For X = {squares}, {hexagons} and the empty set, O is separated, so a(O) = O^+. The
  global sections are, respectively:
  - the set-context (square-coherent) orientations;
  - the hexagon-coherent orientations;
  - all orientations.
- **The gaps.** By THM-4891 (iv), each of these three sheafifications has global sections with 7 and with 21 sinks.

**(F) The +2 operations on the strong-component site.**
1. **Padding** Sigma(T) (a source b and a sink a with b -> a) adds the two singleton points {b} and {a}. So
   a(A)(Sigma(T)) = A^(c(T) + 2). Every J_ss-sheaf that is a point on single vertices, Pk_2 included, is unchanged, so
   H(Sigma T) = H(T).
2. **The gate** Gamma(T) (a -> b -> T -> a) is strong, hence a single point, and a(A)(Gamma(T)) = A.
3. **Second derivation of the gate formula.**
   - The odd cycles of Gamma(T) that use the gate are a -> b -> (a Hamiltonian path of T[S]) -> a with |S| odd, and at
     most one of them occurs in a packing. So OCF gives

         H(Gamma(T)) = H(T) + 2 sum over S subset V with |S| odd of H(T[S]) H(T[V - S]).

   - This agrees with Session 210's H(Gamma T) = H(T) + P_T(2), where P_T(2) = sum over all S of H(T[S]) H(T[V - S]).
     The two sums are equal because the alternating sum over S vanishes (Session 210's Lemma 3).

## Proof

**Lemma 0.**
- In J_2, every covering sieve on S contains the sieve R_2(S) of all subsets of size at most 2, which itself covers.
- In J_ss, every covering sieve on S contains the downward closure of {C_1, ..., C_c}, which covers.
- In J_X, every covering sieve contains the sieve generated by the edges and the X-type 2-faces of F.
- The colimit over covering sieves, ordered by refinement, is therefore attained at R_min.

**(A).**
- **One step.** By Lemma 0, F^+(S) = Match(R_2(S), F) = lim over |U| <= 2 of F(U). For |U| <= 2, R_2(U) is maximal, so
  F^+(U) = F(U). Hence F^(++)(S) = lim over |U| <= 2 of F^+(U) = F^+(S). F^(++) is a sheaf by Grothendieck's lemma, so
  F^+ is too.
- **Linear orders.** A linear order on a pair is an orientation of that pair, and singletons and the empty set carry one
  order. So a compatible family on R_2(S) is exactly a tournament on S, which gives L^+ = Tour.
  - An order is determined by its pairs, so L is separated.
  - L -> L^+ sends an order to its comparability tournament, which is transitive, and every transitive tournament
    arises this way.
- **Gluing failure.** A tournament is transitive iff it has no cyclic triangle (a standard fact). A cyclic triangle is
  not transitive while all its pairs are.

**(B).**
1. **Topology.**
   - *Maximal sieve.* It contains every strong component, so it covers.
   - *Stability.* Let S' subset S. A strong component C' of T[S'] is strongly connected inside T[S], so it lies in a
     strong component C of T[S]. A covering sieve on S contains C, hence C', because sieves are downward closed. So
     R|S' covers S'.
   - *Transitivity.* Let R cover S, and let R' be a sieve with R'|U covering U for every U in R. Each component C of
     T[S] lies in R. Since T[C] is strong, R'|C covering C forces C in R'. So R' covers S.
   - *Points.* A strong S is its own unique component, so only the maximal sieve covers it. A non-strong S has the
     proper covering sieve generated by its components.
2. **Plus formula.** A matching family on R_min(S) is determined by its values x_i in F(C_i). Distinct components share
   only the empty set, so compatibility says exactly x_i|empty = x_j|empty. For S = empty, R_min is the empty sieve, with
   exactly one matching family.
3. **Sheaves.** F is a sheaf iff F -> F^+ is a bijection, and (2) evaluates F^+.
4. **One or two steps.**
   - If F(empty) = pt, the fibre product in (2) is the full product prod_i F(C_i), and F^+ is a sheaf by (3).
   - For a constant presheaf A:
     - for S nonempty, A^+(S) is the fibre product of copies of A over A(empty) = A with identity maps, i.e. the
       diagonal, which is A;
     - A^+(empty) = pt;
     - A^+(C) = A on points C;
     - so A^(++)(S) = prod_i A^+(C_i) = A^c.
   - A^+(S) and A^(++)(S) differ iff c >= 2 and |A| >= 2.
5. **Converse.** Sheafification of a presheaf of groups is computed on underlying sets, so (4) gives (Z/2)^c. On each
   strong component it acts as the converse.
   - H(T') = prod H(C'_i) = prod H(C_i), since reversing a Hamiltonian path of C gives one of C^op.
   - Componentwise converses keep each C_i strong and keep all arcs between components.

**(C).**
1. **H counts global sections.**
   - A directed cycle is strongly connected, so it lies inside a single strong component of T[S].
   - Packings in different components are automatically disjoint, and colourings factor. Hence Pk_2(S) = prod_i Pk_2(C_i),
     and Pk_2(empty) = {empty packing} = pt.
   - By (B3), Pk_2 is a sheaf.
   - |Pk_2(V)| = sum over packings of 2^|p| = I(Omega(T), 2) = H(T), by THM-002.
2. **Rédei.** The swap (p, col) -> (p, 1 - col) fixes (p, col) iff p is empty.
3. **Spectrum.** By (1), H(T) = prod_i H(C_i), and every product of strong values is realised by an order join of the
   corresponding strong tournaments.

**(D).**
- The swap commutes with restriction, so Q is a presheaf. By Q(empty) = pt and (B2), Q^+(S) = prod_i Q(C_i).
- **Fibres.** Take local classes [p_i], and let m be the number of components with nonempty p_i.
  - The components with nonempty packings give 2^m representative tuples. The simultaneous swap acts freely on them
    when m >= 1.
  - So the fibre has 2^(m-1) elements: a torsor under (Z/2)^m / diagonal.
  - If m >= 2, two distinct sections of Q(S) have equal restrictions to every component, so Q is not separated.
- **Counting.** For two components, the fibres over the mixed part (m = 2) have size 2, which gives the two counts
  stated.

**(E).**
- Restriction of coherent orientations to subfaces restricts each block tournament to its sub-blocks. Every pair inside
  a block is an edge of the face, so O is separated for every J_X.
- A matching family on R_min(F) is an orientation of the edges of F that agrees with a section of O on each X-type
  2-face. On a square, that section is a pair of parallel edge classes; on a hexagon, a tournament on its three values.
  So a(O)(F) is the set of orientations coherent on all X-type 2-faces.
- For X = {squares, hexagons}, THM-4891 (i) identifies a(O)(F) with O(F). The spectra are THM-4891 (iv).

**(F).**
1. **Padding.** The components of Sigma(T) are {b}, those of T, and {a}, since b dominates and a is dominated.
2. **Gate.** Gamma(T) is strong (Session 210, Theorem 1).
3. **Gate formula.**
   - An odd cycle of Gamma(T) through a or b uses both, in the form a -> b -> P -> a, with P a Hamiltonian path of
     T[S] for a nonempty S. Its length is |S| + 2.
   - All gate cycles share a and b, so a packing contains at most one of them. The rest of the packing lies in
     T[V - S].
   - Evaluating the independence polynomial at 2 gives the displayed sum.
   - Equality with Session 210's formula is the vanishing of sum over S of (-1)^|S| H(T[S]) H(T[V - S]). That is the
     coefficient of x_V in F_T(x) F_T(-x) = 1 for nonempty V.

QED.
