---
id: THM-4892
title: "Tournament sites and the plus construction. (A) On the pair site, where a covering sieve must contain every pair, linear orders form a separated presheaf. Its sheafification, reached in one plus step, is the presheaf of tournaments; the minimal loci where local orders fail to glue are the cyclic triangles (Condorcet; THM-512). (B) For a fixed tournament T, the strong-component site (a covering sieve must contain every strong component of T[S]) is a Grothendieck topology. By the comparison lemma, its sheaves are the presheaves on the poset of nonempty strong subsets, so a(F)(S) = prod_i F(C_i). The second plus step can act only through the empty set: global data, such as the converse bit, sheafify in exactly two steps to one bit per strong component. (C) H(T) is the number of global sections of the sheaf of 2-coloured odd-cycle packings (OCF); no bijection with Hamiltonian paths is claimed. (D) The presheaf of swap-classes of 2-coloured packings, whose global sections index the cells of Session 210's projective carrier, is not separated. Its separation defect is exactly the (Z/2)^(m-1) relative-colour torsors. (E) THM-4891's relaxations of parallel coherence are the sheafifications of coherent face orientations for finer topologies on the permutohedron's face poset. (F) Padding adds two points to the site and the gate makes the whole vertex set a point; OCF gives an odd-subset form of Session 210's gate formula"
status: >
  PROVED (elementary; new as a sheaf-theoretic formulation):
  Lemma 0, the comparison form, (A1), (B1)-(B5) as site statements, (C1), (D1)-(D3), (D4) for strong factors and in
  general form, (E), and the odd-subset gate formula of (F3).
  KNOWN, re-derived in site language:
  - (A2)-(A3): a tournament is a compatible family of pair orientations, and it is transitive iff it has no cyclic
    triangle (classical). The cyclic triangle is the Condorcet paradox (THM-512, the Arrow-Condorcet bridge).
  - (B5): H is invariant under componentwise converse, since H is multiplicative over strong components and
    H(C^op) = H(C).
  - (C2): Rédei's theorem from the parity of OCF (Grinberg-Stanley, arXiv:2307.05569, who derive the theorems of
    Rédei and Berge; THM-002).
  - (C3): the spectrum is the multiplicative closure of the strong spectrum (THM-1862; THM-1370 (1)).
  - (D) for two components: the two-element relative-colour torsor is Session 210's (STRUCTURE_ATLAS section 2;
    projective-cycle-carrier section 3). New here: the (Z/2)^(m-1) torsor for m components, and the separation reading.
  - (F1), (F2): Session 210 Theorem 1. The gate-cycle bijection used in (F3) is Session 210 Theorem 2.
  INHERITED and credited:
  - Grothendieck's plus construction (Mac Lane-Moerdijk III.5; SGA4 II.3): F^+ is separated; if F is separated, F^+
    is a sheaf; a = (+)(+); F is a sheaf iff F -> F^+ is a bijection.
  - The comparison lemma (SGA4 III.4.1; Johnstone, Elephant C2.2.3).
  - OCF H = I(Omega(T), 2) (Grinberg-Stanley; THM-002).
  - From Session 210: the gate formula (5), Lemma 3, the projective carrier and its join law, Theorems 1, 2 and 4.
  - THM-4891; THM-1370 (h-spectrum-omits-7-21-all-n); THM-4094; THM-1862; THM-512.
  FINITE-EXACT:
  - The session's cardinality checks (`sheafification_plus_construction_20261107/*.out`):
    - (A) for n <= 5 (by the count formula);
    - J_ss pullback stability on 300 random cases;
    - the constant presheaf's two steps on all subsets of 40 random tournaments;
    - OCF and multiplicativity on 1,559 tournaments;
    - |Q| = (H + 1)/2 on 1,559 tournaments and the join laws on 80 order joins (cardinalities only);
    - the (F3) identity on 1,399 tournaments.
  - The independent audit's structural checks (`audit/*.out`):
    - the plus construction as the colimit over all covering sieves, so Lemma 0 itself;
    - the J_2, J_ss and J_X axioms;
    - the fibres of Q and the (Z/2)^M action on them;
    - the comparison formula on 8,410 (object, presheaf) pairs, 402 of the presheaves non-separated;
    - the finer/coarser order of the face topologies.
  NOT CLAIMED:
  - any topos, Lawvere or "2 = |Omega_Set|" bridge (cf. MISTAKE-529);
  - any functor identifying the gate with the plus construction. "Plus applied twice" and "+2" share only the number
    2; this is a verbal ANALOGY in the source note.
  - any bijection between Pk_2(V) and Hamiltonian paths;
  - any new H-value or progress on Chan-Pak (OPEN).
  AUDIT: PASS WITH CORRECTIONS (`audit/AUDIT.md`; MISTAKE-731). There were 4 majors, all applied:
  - (E) had finer and coarser inverted;
  - the sheafified join count needs strong factors;
  - "Hamiltonian paths as global sections" was an overclaim;
  - credits were missing.
  The 16 minors and 4 upgrades are also applied.
session: opus-2026-11-07-S85
source: 05-knowledge/results/sheafification_plus_construction_20261107/README.md
related:
  - THM-002
  - THM-512
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-1862
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
  - 05-knowledge/results/sheafification_plus_construction_20261107/audit/aud_sites_colimit.py
  - 05-knowledge/results/sheafification_plus_construction_20261107/audit/aud_ocf_q_gate.py
  - 05-knowledge/results/sheafification_plus_construction_20261107/audit/aud_perm_site.py
  - 05-knowledge/results/sheafification_plus_construction_20261107/audit/aud_gate_points.py
---

# THM-4892 - tournament sites and the plus construction

## Background (classical)

For a presheaf F of sets on a site, the plus construction is

    F^+(S) = colim over covering sieves R of S of Match(R, F).

Here Match(R, F) is the set of matching families, and the colimit runs over covering sieves ordered by reverse
inclusion. Grothendieck's facts:
- F^+ is always separated, and it may already be a sheaf.
- Applied to a separated presheaf, the construction yields a sheaf.
- So F^(++) = a(F), and F is a sheaf iff F -> F^+ is a bijection.

**Lemma 0.** On each site below, every object S has a smallest covering sieve R_min(S). It is terminal in the index
category, so F^+(S) = Match(R_min(S), F).

**Comparison form (the audit's UPGRADE-1).** On each site below, the irreducible objects (those whose only covering
sieve is the maximal one) form a dense full subposet, and the topology it inherits is trivial. By the comparison lemma,
sheaves are the same as presheaves on that subposet. For every presheaf F,

    a(F)(S) = lim over irreducible U contained in S of F(U).

**Where the second step acts (the audit's UPGRADE-2).** On these sites R_min(S) is generated by irreducible objects.
The map F^+(S) -> F^(++)(S) can fail to be onto only through an overlap W of two generators that is itself covered by a
smaller sieve. Such overlaps are:
- the empty set on the strong-component site;
- the vertices on the face sites;
- none on the pair site.

## Statement

**(A) The pair site.** Let J_2 be the topology on the poset of finite subsets of a set V in which a sieve R on S
covers iff R contains every subset of S with at most two elements.
1. **One step.** For every presheaf F, F^+(S) = lim over U subset S with |U| <= 2 of F(U). F^+ is already a sheaf, so
   a(F) = F^+.
2. **Linear orders.** Let L(S) be the set of linear orders of S, with restriction to the induced order.
   - L is separated, and not a sheaf once |S| >= 3.
   - L^+ = Tour, the presheaf of tournaments, with restriction to the induced subtournament.
   - L -> a(L) = Tour is the inclusion of the transitive tournaments.
3. **Gluing failure.** A section of Tour over S lies in the image of L iff every triple of S does. So the minimal
   subsets on which local linear orders fail to glue are the cyclic triangles (the Condorcet paradox, THM-512).

**(B) The strong-component site.** Fix a tournament T on V. Let J_ss be the topology on the poset of subsets of V in
which a sieve R on S covers iff R contains every strong component of T[S].
1. **Topology and irreducible objects.** J_ss is a Grothendieck topology.
   - Its irreducible objects are the nonempty S with T[S] strongly connected.
   - By the comparison lemma, Sh(P(V), J_ss) is equivalent to the presheaves on the poset Str(T) of these S, so they
     are exactly the points of the topos (evaluation at S).
   - The empty set is covered by the empty sieve.
2. **Plus formula.** If C_1, ..., C_c are the strong components of T[S], then

       F^+(S) = F(C_1) x_F(empty) ... x_F(empty) F(C_c),   and F^+(empty) = pt.

3. **Sheaves.** F is a sheaf iff F(empty) = pt and F(S) -> prod_i F(C_i) is a bijection for all S. So J_ss-sheaves are
   exactly the presheaves that are multiplicative over strong components.
4. **One step or two.**
   - For every presheaf F, a(F)(S) = F^(++)(S) = prod_i F(C_i), and F^+(S) -> F^(++)(S) is the inclusion of the
     fibre product over F(empty).
   - This inclusion is a bijection unless c(T[S]) >= 2, every F(C_i) is nonempty, and the images of the maps
     F(C_i) -> F(empty) are not all one and the same singleton. In particular a(F) = F^+ whenever F(empty) = pt.
   - For a constant presheaf A (A(S) = A for every S, the empty set included):
     - A^+(S) = A, the diagonal, for nonempty S;
     - A^+(empty) = pt;
     - A^(++)(S) = A^(c(T[S])).
     So when |A| >= 2, two steps are needed exactly at the S with c(T[S]) >= 2.
5. **Global symmetry becomes local symmetry.**
   - The constant presheaf of groups Z/2 sheafifies to a(Z/2)(S) = (Z/2)^(c(T[S])). Two steps are needed exactly at
     the S with c(T[S]) >= 2.
   - Let Z/2 act on the presheaf Tour of all tournaments by the converse. This action commutes with restriction, so
     a(Z/2) acts on a(Tour), where a(Tour)(V) = prod_i Tour(C_i) and (Z/2)^(c(T)) acts by componentwise converse.
   - Reversing the arcs inside any set of strong components of T keeps the components, the arcs between them, and H
     (classical). So (T_i) -> prod_i H(T_i) is invariant under a(Z/2).

**(C) H as the number of global sections.** Let Pk_2(S) be the set of pairs (p, col): p is a set of pairwise
vertex-disjoint directed odd cycles of T[S], and col: p -> {0, 1}. Restriction to S' subset S keeps the cycles contained
in S'.
1. **Count.** Pk_2 is a J_ss-sheaf, and |Pk_2(V)| = H(T) (OCF, THM-002). No bijection with Hamiltonian paths is
   claimed. The law H(A => B) = H(A) H(B) for order joins (classical; THM-1862) is, via OCF, the cardinality of the
   sheaf condition for the covering sieve generated by {A, B}.
2. **Rédei (KNOWN).** The global colour swap is an involution of Pk_2(V) whose only fixed point is the empty packing.
   Hence H(T) is odd (Grinberg-Stanley's parity consequence of OCF).
3. **Spectrum (KNOWN).** {H(T)} is the multiplicative monoid generated by {H(T) : T strong} (THM-1862; THM-1370 (1)).
   Pk_2 gives the packing-level form of this reduction.

**(D) The swap quotient Q of Pk_2 is not separated.** Let Q = Pk_2/(global swap), the presheaf quotient. Q(V) is the
set of cells of Session 210's projective carrier Q_T (STRUCTURE_ATLAS section 1's orbit set Q(T)) together with one
basepoint, the empty packing. So |Q(S)| = (H(T[S]) + 1)/2.
1. **Separation.** Q is not J_ss-separated at S as soon as two strong components of T[S] carry nonempty packings.
   Equivalently, Q is separated at S iff at most one strong component of T[S] has at least 3 vertices.
2. **One step.** Since Q(empty) = pt, a(Q)(S) = Q^+(S) = prod_i Q(C_i), of size prod_i (H(C_i) + 1)/2.
3. **Torsors.** Q(S) -> Q^+(S) is surjective. The fibre over an element whose packing meets m components is a torsor
   under (Z/2)^m / diagonal = (Z/2)^(m-1). For m = 2 these are Session 210's relative-colour torsors.
4. **Counting.** Put q = (H - 1)/2.
   - **Presheaf count.** For every order join A => B, |Q| - 1 = q_A + q_B + 2 q_A q_B (Session 210's join law).
   - **General count.** |Q^+(A => B)| = |Q^+(A)| |Q^+(B)| = prod over the strong components C of A and B of
     (H(C) + 1)/2.
   - **Strong factors.** If A and B are strong (more generally, if each has at most one strong component carrying an
     odd cycle), then |Q^+| - 1 = q_A + q_B + q_A q_B. The factor 2 on q_A q_B is the separation defect.
   - **Example.** A = C3 => C3 and B = C3 give |Q^+| - 1 = 7, not 9.

**(E) THM-4891 as a table of sheafifications.** On the face poset of the permutohedron (ordered set partitions, ordered
by inclusion of faces), let O(B_1 | ... | B_k) = prod_j Tour(B_j) be the presheaf of coherent face orientations. For a
set X of 2-face types, let J_X be the topology in which a sieve on a face F covers iff it contains every edge of F and
every 2-face of F of a type in X. A topology is finer when it has more covering sieves. So J_{squares, hexagons} is the
coarsest of the four, and J_{} (edges only) is the finest.
- **The coarsest.** For X = {squares, hexagons}, O is a sheaf. Its global sections are the tournaments (THM-4891 (i)).
- **The finer topologies.** Every J_{squares, hexagons}-covering sieve is also J_{squares}-, J_{hexagons}- and
  J_{}-covering. For X = {squares}, {hexagons} and the empty set, O is separated but not a sheaf, so a(O) = O^+. The
  global sections are, respectively:
  - the set-context (square-coherent) orientations;
  - the hexagon-coherent orientations;
  - all orientations.
- **The gaps.**
  - By THM-4891 (iv)-(v), the square-coherent and the hexagon-coherent sheafifications have acyclic global sections
    with 7 and with 21 sinks (n = 5).
  - These sections are also global sections for X empty. Session 210's Theorem 4 shows that the acyclic orientations
    alone already realise every sink count 1, ..., n!/2.
  - The exclusions of 7 and 21 hold for O itself, the sheaf for the coarsest of the four topologies (THM-1370).

**(F) The +2 operations on the strong-component site.**
1. **Padding.** Sigma(T) (a source b and a sink a with b -> a) adds the two singleton points {b} and {a}, and no
   others, so a(A)(Sigma(T)) = A^(c(T) + 2).
   - For a family of sheaves defined uniformly for all tournaments, with one-element values on single vertices (Pk_2 is
     one), the value at Sigma(T) equals the value at T. So H(Sigma T) = H(T).
2. **Gate.** For nonempty T, Gamma(T) (a -> b -> T -> a) is strong (Session 210, Theorem 1). So its whole vertex set is
   a point, and a(A)(Gamma(T)) = A.
   - The site of Gamma(T) has 2^n + 1 more points than the site of T: {a}, {b}, and {a, b} u S for every nonempty S.
3. **Odd-subset gate formula.** The odd cycles of Gamma(T) that use the gate are exactly a -> b -> (a Hamiltonian path
   of T[S]) -> a with |S| odd. This description is Session 210's Theorem 2 restricted to odd cycles. Since at most one
   such cycle occurs in a packing, OCF gives

       H(Gamma(T)) = H(T) + 2 sum over S subset V, |S| odd, of H(T[S]) H(T[V - S]),   H(T[empty]) = 1.

   This agrees with Session 210's H(Gamma T) = H(T) + P_T(2), where P_T(2) = sum over all S of H(T[S]) H(T[V - S]),
   because sum over S of (-1)^|S| H(T[S]) H(T[V - S]) = 0 (Session 210's Lemma 3).

## Proof

**Lemma 0, the comparison form, and the second step.**
- R_min(S) is:
  - for J_2, the sieve of subsets of size at most 2;
  - for J_ss, the downward closure of the strong components;
  - for J_X, the sieve generated by the edges and the X-type 2-faces.
- Each is covering and contained in every covering sieve, so it is terminal.
- The irreducible objects are:
  - for J_2, the subsets of size at most 2;
  - for J_ss, the nonempty strong subsets;
  - for J_X, the edges and the X-type 2-faces.
- They are dense, since every object is covered by them, and the induced topology on them is trivial. The comparison
  lemma applies, and a(F)(S) is the limit over the irreducible objects below S. The audit verified this formula as a
  colimit computation over all covering sieves.
- **The second step.** F^(++)(S) = Match(R_min(S), F^+) consists of families y_G in F^+(G) = F(G), one for each
  generator G. They are compared in F^+(W) on overlaps W.
  - F(W) -> F^+(W) is the identity when W is irreducible.
  - So the compatibility conditions of F^+ and F^(++) differ only at overlaps W that are covered by smaller sieves.
  - On J_2 the overlaps (singletons and the empty set) are irreducible.
  - On J_ss the overlap of two components is the empty set.
  - On J_X the overlaps that are not irreducible are the vertices.

**(A).**
- **One step.** By Lemma 0, F^+(S) = lim over |U| <= 2 of F(U), and F^+(U) = F(U) for |U| <= 2. Under these
  identifications the canonical map F^+ -> F^(++) is the identity. Since F^(++) is a sheaf, so is F^+.
- **Linear orders.** A linear order on a pair is an orientation of that pair. So a compatible family on R_min(S) is
  exactly a tournament on S, which gives L^+ = Tour.
- **Gluing failure.** A tournament is transitive iff it has no cyclic triangle (classical). A cyclic triangle is
  intransitive while all its pairs are linearly ordered.

**(B).**
1. **Axioms.**
   - The maximal sieve covers.
   - Stability: a strong component of T[S'] is strongly connected in T[S], so it lies in a component of T[S]. Sieves
     are downward closed.
   - Transitivity: each component C of T[S] lies in R, and since T[C] is strong this forces C to lie in R'.
   - The irreducible objects and the equivalence are as in the comparison form.
2. **Plus formula.** A matching family on R_min(S) is determined by its values x_i in F(C_i). Distinct components meet
   only in the empty set, so compatibility says x_i|empty = x_j|empty.
3. **Sheaves.** This is the comparison form.
4. **One step or two.**
   - F^+(C_i) = F(C_i), since C_i is irreducible.
   - By the comparison form, a(F)(S) = prod_i F(C_i).
   - The fibre product in (2) equals the full product exactly under the stated condition.
   - For the constant presheaf, every map A -> A(empty) = A is the identity.
5. **Groups and the action.**
   - Sheafification of a presheaf of groups is computed on underlying sets.
   - The converse commutes with restriction, so the action passes to a(Tour) = prod_i Tour(C_i) componentwise.
   - Reversing a component C keeps it strong, and H(C^op) = H(C).

**(C).**
- A directed cycle is strongly connected, so it lies in one component.
- Packings in different components are disjoint, and colourings factor. So Pk_2(S) = prod_i Pk_2(C_i), and
  Pk_2(empty) = {empty packing}. By (B3), Pk_2 is a sheaf.
- |Pk_2(V)| = sum over packings of 2^|p| = I(Omega(T), 2) = H(T) by THM-002.
- The swap (p, col) -> (p, 1 - col) fixes (p, col) iff p is empty.
- The spectrum statement is THM-1862 and THM-1370 (1).

**(D).**
- The swap commutes with restriction, so Q is a presheaf. Q(empty) = pt, and (B2) gives Q^+(S) = prod_i Q(C_i).
- **Fibres.** Take local classes [p_i], and let m be the number of components with nonempty p_i.
  - These give 2^m representative tuples.
  - The simultaneous swap acts freely on them when m >= 1, so the fibre has 2^(m-1) elements: a torsor under
    (Z/2)^m / diagonal.
- **Separation.** If m >= 2, two distinct sections have equal restrictions to every component.
- **Strong factors.** For A and B strong, the components of A => B are A and B. The fibres over the mixed part (m = 2)
  have size 2, which gives the two counts. The general product formula is (B4) applied to all components of A and B.

**(E).**
- **Axioms.**
  - Each edge and each X-type 2-face is required of itself, which gives transitivity.
  - Edges and X-faces of a subface are edges and X-faces of the face, which gives stability.
  - A vertex has no edges, so every sieve on a vertex covers it, the empty sieve included, and O(vertex) = pt.
- **Order of the topologies.** Containing more required faces is a stronger requirement on a sieve, so the J_X with
  larger X have fewer covering sieves.
- **Separation.** Every pair inside a block is swapped by an edge of the face, so O is separated.
- **Sheafification.** A matching family on R_min(F) is an orientation of the edges of F that agrees with a section of O
  on each X-type 2-face. On a square that section is a parallel pair of edge classes; on a hexagon it is a tournament
  on its three values. So a(O)(F) is the set of orientations coherent on all X-type 2-faces.
- **The sheaf case.** For X = {squares, hexagons}, apply THM-4891 (i) to each factor of F ~= prod_j Perm(B_j). The
  squares that mix two blocks make the factors independent. This identifies a(O)(F) with O(F).
- **Spectra.** The global-section sink spectra are THM-4891 (iv). The case of X empty is covered by inclusion and by
  Session 210's Theorem 4.

**(F).**
- The components of Sigma(T) are {b}, those of T, and {a}.
- Gamma(T) is strong by Session 210 Theorem 1. Its additional strong subsets are as listed (audit,
  `aud_gate_points.out`).
- **Gate formula.**
  - An odd cycle through a or b uses both, and has the form a -> b -> P -> a, with P a Hamiltonian path of T[S] and
    |S| = length - 2 odd (Session 210 Theorem 2).
  - All such cycles share a and b, so a packing contains at most one of them, and the rest of the packing lies in
    T[V - S].
  - Evaluating the independence polynomial at 2 gives the formula.
  - The vanishing of the alternating sum is the coefficient of x_V in F_T(x) F_T(-x) = 1 for nonempty V.

QED.
