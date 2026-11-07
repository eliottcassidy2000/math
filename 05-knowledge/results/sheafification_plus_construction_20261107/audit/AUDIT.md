# Independent audit: THM-4892 and the S85 sheafification note

Audited deliverables (opus-2026-11-07-S85, checkpoint `4d3b0f90f4`):
- `01-canon/theorems/THM-4892-tournament-sites-plus-construction.md`
- `05-knowledge/results/sheafification_plus_construction_20261107/README.md`, `sheaf_check.py/.out`,
  `gate_packing_check.py/.out`, `q_torsor_check.py/.out`
- the two INDEX rows for THM-4892 (`05-knowledge/hypotheses/INDEX.md` line 2, `05-knowledge/results/INDEX.md` line 6)

Background read: THM-002, THM-4891, THM-4900, THM-796 (sections 2-3), THM-801 (sections 1-2), THM-1370, THM-4094,
THM-1862, THM-512, MISTAKE-529, `session210/{coherent-permutohedron, graded-spectrum-gates, projective-cycle-carrier}/PROOF.md`,
`session210/structure-atlas/STRUCTURE_ATLAS.md`.

The auditor wrote only into this `audit/` folder. No canon, results or index file was edited, and nothing was committed.
All code here is the auditor's own; nothing imports the session's modules.

## Verdict: PASS WITH CORRECTIONS

The sheaf theory is correct. All three topologies are Grothendieck topologies. Lemma 0 holds: the plus construction
computed as a colimit over every covering sieve agrees with the smallest-sieve formula. The plus formula, the sheaf
criterion, the one-step and two-step statements, the constant-presheaf computation, the Pk_2 sheaf, the torsor
fibres of Q and the gate formula are all correct, and all were re-derived and re-computed independently.

Four corrections are needed before the theorem is final:
1. **(E) inverts finer and coarser.** With the standard convention, the relaxations are sheafifications for *finer*
   topologies, not coarser ones.
2. **The sheafified count in (D)4 is false for a general order join.** It needs strong factors. Counterexample:
   A = C3 => C3 and B = C3 give 7, not 9.
3. **"Hamiltonian paths as global sections" overclaims.** The sections are 2-coloured packings; OCF makes them
   equinumerous with Hamiltonian paths, but no bijection is given.
4. **Several parts are typed as new but are known.** Rédei via the swap is Grinberg-Stanley's parity consequence of
   OCF. Spectrum = multiplicative closure is THM-1862 / THM-1370 (1). The gate-cycle bijection is Session 210's
   Theorem 2. The two-component relative-colour torsor is Session 210's. The Condorcet reading is THM-512.

The minor items are wording: "points", the "gate collapses everything to one" phrasing, and README sentences that
overstate. They also cover verification scripts whose torsor and count checks are tautological.

---

## 1. Line-by-line check of the proofs (task 1)

Conventions used for checking: sites are posets of subsets (resp. faces) with inclusions; a sieve on S is a down-set
of objects below S; J is a Grothendieck topology iff (i) the maximal sieve covers, (ii) R covers S and S' <= S imply
R|S' = {U in R : U <= S'} covers S', (iii) R covers S and R'|U covers U for every U in R imply R' covers S.
F^+(S) = colim_{R in J(S)} Match(R, F), the index category being J(S) ordered by reverse inclusion.

**Background and Lemma 0.**
- Correct: if R_min(S) is a covering sieve contained in every covering sieve, it is terminal in the index
  category, so the colimit is Match(R_min(S), F).
- Lines 57-58 need a wording fix (MINOR-9). "The first can only make F separated" is imprecise: F^+ is always
  separated, and it is often already a sheaf.

**(A) The pair site.**
- **Axioms (i)-(iii).** Correct. For (iii), each U with |U| <= 2 lies in R, and R'|U covers U forces U in R'.
- **(A1).** Correct. Match(R_2(S), F) = lim_{|U|<=2} F(U), and F^+(U) = F(U) for |U| <= 2, since these objects
  have only the maximal covering sieve.
  - Gap (MINOR-16): the proof should say that, under these identifications, the canonical map F^+ -> F^++ is the
    identity. Only then does "F^++ is a sheaf, so F^+ is too" follow.
- **(A2), (A3).** Correct and classical. A tournament is a choice of orientation on each pair, and a tournament is
  transitive iff every triple is.
- **Computational confirmation.** The colimit was computed over every covering sieve (17 on a 4-set) for 14 test
  presheaves, including non-separated ones. It confirms (A1) for every one of them.

**(B) The strong-component site.**
- **Maximal sieve and stability.** Correct, as written at lines 166-169. A strong component of T[S'] is strongly
  connected inside T[S], so it lies in a component of T[S].
- **Transitivity.** Correct (lines 170-171).
- **"Points".** Correct under the theorem's own definition, but that definition is not the standard meaning of the
  word (MINOR-1).
  - The standard meaning agrees here, and a short argument proves it: by the comparison lemma,
    Sh(P(V), J_ss) is equivalent to [Str(T)^op, Set], where Str(T) is the poset of nonempty strong subsets (UPGRADE-1).
  - The points of a presheaf topos on a finite poset are exactly the evaluations at its elements.
- **The empty set** is covered by the empty sieve. Correct. It forces F(empty) = pt for every sheaf.
- **(B2).** Correct. Distinct components meet only in the empty set, which lies in R_min(S).
- **(B3).** Correct. It uses "F is a sheaf iff F -> F^+ is a bijection". This is standard (if F -> F^+ is
  bijective, F is separated because F^+ is, so F^+ and hence F is a sheaf), but it is not in the INHERITED list
  (MINOR-16).
- **(B4).**
  - Bullet 1 is correct. The proof should add F^+(C_i) = F(C_i), since C_i is irreducible.
  - Bullet 2 is imprecise (MINOR-2). For every presheaf, a(F)(S) = prod_i F(C_i).
    - The map F^+(S) -> a(F)(S) is the inclusion of the fibre product.
    - It fails to be onto exactly when c >= 2, every F(C_i) is nonempty, and the images F(C_i) -> F(empty) are not
      all one common singleton.
    - Verified on 8,410 (object, presheaf) pairs.
  - The constant-presheaf computation A^+(S) = A (S nonempty), A^+(empty) = pt, A^++(S) = A^c is correct.
- **(B5).** The group computation is correct, and the sheafification of a presheaf of groups is indeed computed on
  underlying sets.
  - "Generated by the converse" and "this is the group of componentwise converses" are interpretations, not yet
    statements (MINOR-3). They become one when the constant group acts on the presheaf Tour by converse (UPGRADE-3).
  - The H-invariance is classical: H is multiplicative over strong components and H(C^op) = H(C).
  - "Exactly two steps" refers to the presheaf as a whole. At the global object it holds iff c(T) >= 2.

**(C) Pk_2.**
- **Presheaf and sheaf.** It is a presheaf: keeping the cycles inside S' is functorial. It is a J_ss-sheaf: a
  directed cycle lies in one strong component, and the restriction map Pk_2(S) -> prod_i Pk_2(C_i) is a bijection.
- **(C1).** |Pk_2(V)| = sum_p 2^|p| = I(Omega, 2) = H. This is THM-002 verbatim.
  - The heading "Hamiltonian paths as global sections" is not what is proved (MAJOR-3).
- **(C2).** Correct, and KNOWN. Grinberg-Stanley (arXiv:2307.05569, abstract) recover Rédei's theorem from their
  formula, and THM-002 notes that I(Omega, 2) is odd. The swap is the bijective form of "2^k is even for k >= 1".
- **(C3).** Correct, and KNOWN: THM-1862 (H is multiplicative under order join, with the strong tournaments as its
  atoms) and THM-1370 (1). Neither is credited (MAJOR-4).

**(D) Q = Pk_2 / global swap.**
- **Presheaf, values and fibres.** The swap commutes with restriction, so Q is a presheaf, and Q(empty) = pt.
  |Q(S)| = (H + 1)/2. Q^+ = a(Q) = prod_i Q(C_i). The fibres are torsors under (Z/2)^m / diagonal. All correct;
  the m = 2 case is Session 210's relative-colour torsor (STRUCTURE_ATLAS section 2, projective-cycle-carrier
  section 3).
- **Separatedness.** Q is not separated, and in fact Q is separated at S iff at most one strong component of T[S]
  has size >= 3 (verified).
- **(D)4.** The presheaf count holds for every order join. The sheafified count q_A + q_B + q_A q_B holds only when
  Q(A) = Q^+(A) and Q(B) = Q^+(B), for instance when A and B are strong. As stated it is false (MAJOR-2). The proof
  sentence "For two components ..." silently assumes A and B are the components.
- **Terminology.** "Q(V) is Session 210's projective carrier" conflates the CW complex Q_T with its cell set
  (MINOR-5).

**(E) The permutohedron face poset.**
- **Axioms of J_X, checked carefully.**
  - *Maximal sieve.* It covers.
  - *Stability.* Edges and X-type 2-faces of a subface F' <= F are edges and X-type 2-faces of F, so R|F' contains
    all of them.
  - *Transitivity.* Every required piece is required of itself: an edge e is the only edge of e, and an X-face Phi
    is an X-face of Phi. So R'|e covering e forces e in R', and R'|Phi covering Phi forces Phi in R'.
  - *Vertices.* A vertex has no edges, so every sieve on a vertex covers it, including the empty sieve. Every
    J_X-sheaf is therefore a point on vertices, and O(vertex) = pt.
  - *Checks.* Exhaustive for n = 3 and on all faces of dimension <= 2 for n = 4. Randomised, with the R_min
    candidates included, on the top face for n = 4.
- **O.**
  - O is a presheaf.
  - O is separated for every X: every pair inside a block is swapped by some edge of the face.
  - O is a sheaf exactly for X = {squares, hexagons}.
- **Global sections.** Those of a(O) = O^+ are the X-coherent orientations, with 2^(classes) of them. The class
  counts are 6, 24, 12, 36 at n = 4 and 10, 80, 60, 240 at n = 5, matching THM-4891's dictionaries.
- **Gaps in the proof.**
  - THM-4891 (i) is stated for the whole permutohedron. The proof applies it to every face F, which is
    isomorphic to prod_j Perm(B_j); this needs one sentence about the squares that mix two blocks (MINOR-6).
  - **Terminology error:** the "finest/coarser" labels are inverted (MAJOR-1).

**(F) Padding and gate.**
- **(F1).** Correct. The site of Sigma(T) has exactly two more points.
- **(F2).** Γ(T) is strong for nonempty T, so its top object is a point. But the site of Γ(T) has 2^n + 1 *more*
  points than the site of T, so "collapses everything to one" is wrong as a statement about the site (MINOR-7).
  Verified on 195 tournaments. The extra points are {a}, {b} and {a, b} u S for every nonempty S.
- **(F3).** Correct. It uses H(T[empty]) = 1, which should be stated.
  - The description of the gate cycles is Session 210's Theorem 2, restricted to odd cycles, and is not credited
    (MINOR-8).
  - The equivalence with Session 210's (5) is the coefficient of x_V in F_T(x) F_T(-x) = 1 (Session 210, Lemma 3).
    The sign (-1)^|V - S| versus (-1)^|S| is immaterial for vanishing.

---

## 2. Independent computations (task 2)

All counts below are from the auditor's scripts. Every reported condition is True; there are no failures.

| task | check | universe | result |
|---|---|---|---|
| (e) | J_2 axioms (exhaustive transitivity); L^+ and F^+, F^++ as the colimit over **all** covering sieves | n = 0..4, 14 presheaves including 2 non-separated ones and 6 random quotient presheaves | Lemma 0 holds. \|L^+(S)\| = 2^C(\|S\|,2). L is separated, and is a sheaf exactly when \|S\| <= 2. The image of L is the transitive tournaments (1, 1, 2, 6, 24). Membership is decided on triples. F^+ is a sheaf for every presheaf (A1). |
| (B) | J_ss axioms; irreducible objects; R_min; F^+, F^++ by the full colimit | all 75 labelled tournaments with n <= 4, plus 3 random with n = 5; 565 presheaf instances, 402 of them non-separated | 0 axiom violations. Irreducible objects = nonempty strong S (1,194 objects). Lemma 0, (B2), (B3), (B4) and the exact second-step criterion hold, and a(F)(S) has size prod \|F(C_i)\| (8,410 pairs each). |
| (a) | constant presheaf A (\|A\| = 2, 3; A(empty) = A), one and two steps, all subsets | 60 random tournaments, n = 3..6 | 3,216 subset checks of A^+ and A^++. All 2,158 subsets with c >= 2 need both steps. |
| (b) | OCF and the explicit bijection Pk_2(S) -> prod Pk_2(C_i) | all n <= 5 (every subset); random n = 6 (300), 7 (120), 8 (25) | OCF holds on 1,544 tournaments. The bijection holds on 34,311 (T, S) pairs. The swap fixes only the empty packing. H is multiplicative. Spectrum(n) = products of strong values for n <= 5. |
| (c) | Q explicitly: orbits, the map to prod Q(C_i), fibres, the (Z/2)^M action, separatedness, and both plus steps | all n <= 5, random n = 6, 7, 60 order joins of strong pieces, 31 random non-separated cases | 39,302 (T, S) pairs and 55,613 fibres (143 with m = 2, 24 with m = 3). Each fibre has size 2^(m-1), and the action is transitive with stabiliser exactly the diagonal. Q^+ -> Q^++ is bijective (6,208 objects). The (D)4 counterexample is described below this table. |
| (d) | gate formula; Session 210 (5); alternating identity; enumeration of the gate cycles; OCF on Γ(T); padding | all n <= 5 plus random n = 6, 7, 8 (1,339 tournaments) | All hold. Every gate odd cycle has the form a -> b -> P -> a, with P a Hamiltonian path of T[S] and \|S\| odd (1,179 tournaments with n <= 6). |
| (E) | J_X axioms; the order of the four topologies; O separated and sheaf; colimit plus; global sections; sink spectra; THM-4891 witnesses | face posets n = 3, 4 (axioms, colimit), n = 5 (sections) | J_{sq,hex} <= J_{sq}, J_{hex} <= J_{} on every face, strictly somewhere. Hexagon-coherent sink spectrum at n = 4: {0..9, 12}, with 48 orientations at 7 sinks. Square-coherent: {0..6} (all 2^24). Tournaments: {1, 3, 5}. The n = 5 tables give 7 and 21 sinks. |
| (B5) | componentwise converses | 300 random T, n = 2..8 | 1,468 cases: the strong components and H are unchanged. |

**Counterexample to (D)4.** Take A = C3 => C3 (not strong, q_A = 4) and B = C3 (q_B = 1).
- |Q(A => B)| - 1 = 13 = q_A + q_B + 2 q_A q_B, as stated.
- |Q^+(A => B)| - 1 = 7, not the stated q_A + q_B + q_A q_B = 9.

---

## 3. Prior work and literature (task 3)

**Spot-check of README section 1.** All rows are accurate unless noted.
- **category_deep2_S71o.py.** The "Sheaf condition verified" counter is incremented under `if True:` at line 296,
  inside the loop at lines 288-298. Confirmed.
- **The HP presheaf has no restriction maps.** Confirmed in both files. symbolic_abstract_S71q.py (lines 1040-1075)
  states F(T) = {HPs of T} on embeddings, with no contravariant action. tau_projective_duality_S71s.py
  (lines 761-766) states F(S) = {HPs of T|_S}.
- **beyond_ouroboros_S71u.py PART 6** is a "sheaf" on {Comb, Alg, Anal}, a metaphor. beyond_equation_S71x.py and
  new_object.py are metaphors too.
- **HYP-1479 and HYP-1511** exist only as rows of `INDEX-HISTORICAL-THROUGH-2026-07-21.md`, and those rows still say
  CONFIRMED (MINOR-13).
- **MISTAKE-529** retracts HYP-1556 and refutes HYP-1557. Accurate.
- **THM-796.** (2.2) X_n ~= (X_(n-1) x_[X_(n-2)] X_(n-1)) x F_2(apex), and (3.3) is a uniformly two-to-one C_2
  endpoint-phase torsor for n >= 5. Accurate.
- **THM-801.** At n = 5, 32 of the 64 pairwise-compatible triples have nonzero holonomy and do not glue. Accurate.
- **THM-813 and THM-818** exist; THM-818 is a three-way Čech join.
- **Remaining rows.** All present and accurately described:
  - HYP-6825, and HYP-2245 (OPEN);
  - the LRC sheaf HYPs and THM-2542, 2591, 2657, 2751;
  - session202 overlap-gluing (1), and session195 pathcover-algebra Theorem 3 (no fixed-dimensional bilinear
    algebra);
  - session165 fi_interfaces ([FI, Set]);
  - THM-4900 statement 4 (y -> (2y+1)/(y+1), the square of y -> 1 + 1/y);
  - HYP-1193 and HYP-1231;
  - `operads_derived_23.py:207`, the only Quillen-plus mention in tracked files, with no tournament link.

**Omissions from section 1** (MINOR-13):
- `04-computation/tournament_category_23.py`: "SHEAVES on the arc-flip poset ... Čech cohomology of this sheaf =
  obstructions to global HP count";
- `04-computation/euler_identity_S71w.py` and `functors_S71y.py`: the presheaf topos Set^{TOURN^op};
- THM-512, the Arrow-Condorcet bridge. It identifies the tiling cube with the social-choice cube, transitive with
  rational, and the cyclic triangle with the Condorcet paradox. It is the in-repository precursor of the pair-site
  reading in (A) and README section 3.
- No precise earlier sheafification or plus construction on tournaments was found. The search used git grep over
  every tracked file for sheafif|plus construction|grothendieck topolog|covering sieve|matching famil|presheaf
  (tournament files), and git log for same-day work. The README's main prior-work claim stands.

**External citations.**
- **Abramsky-Brandenburger.** New J. Phys. 13 (2011) 113036, arXiv:1102.0264: contextuality = obstruction to global
  sections. Accurate.
- **Sargsyan.** "Localizing Preference Aggregation Conflicts: A Graph-Theoretic Approach Using Sheaves",
  arXiv:2512.02416 (v1 2 Dec 2025, v2 15 Jul 2026). The paper exists, but the README's one-line description is
  inaccurate (MINOR-13).
  - The base is a graph of *voters*. Vertex stalks are total orders on each voter's visible alternatives, and edge
    stalks are total orders on the overlaps.
  - The Condorcet triangle appears as the smallest instance, where merging gives an empty stalk (a "local
    impossibility").
  - The paper does not use pairwise-preference stalks, Grothendieck topologies, sheafification, tournaments or
    Hamiltonian paths.
- **Literature search.** Neither "tournaments are the sheafification of linear orders on the pair site" nor "H as
  the global sections of an odd-cycle packing sheaf" was found.
  - Nearest context: Abramsky et al., "Contextuality, Cohomology and Paradox" (arXiv:1502.03097): cyclic paradoxes
    as local consistency with global inconsistency. Also arXiv:2601.07283 ("Condorcet's Paradox as
    Non-Orientability", topological).
  - Typing for (A): the content is folklore, and only the J_2 formulation is new.

---

## 4. Findings

### MAJOR

**MAJOR-1 (E): finer and coarser are inverted.**
- **The convention.** J' is finer than J iff J(c) is contained in J'(c) for every c, i.e. iff J' has more covering
  sieves (nLab "Grothendieck topology"; same convention in SGA4 Exp. II; the canonical topology is the *finest*
  subcanonical one). A finer topology has fewer sheaves.
- **The four topologies.** Every J_{sq,hex}-covering sieve is J_{sq}- and J_{hex}-covering, and these are J_{}-covering
  (verified on every face for n <= 4; strict). So J_{sq,hex} is the *coarsest* of the four and J_{} the *finest*.
  This matches O being a sheaf only for J_{sq,hex}.
- **Where the error appears.** The theorem's title (E) and lines 121 and 123, README lines 23, 94 and 96, and both
  INDEX rows all say the opposite.
- **Replacements.**
  - *Title (E):* "(E) THM-4891's relaxations of parallel coherence are the sheafifications of coherent face
    orientations for finer topologies on the permutohedron's face poset (topologies requiring fewer 2-faces have more
    covering sieves)".
  - *Line 121:* "- **The coarsest of the four topologies.** For X = {squares, hexagons}, O is a sheaf. Its global
    sections are the tournaments (THM-4891 (i))."
  - *Line 123:* "- **The finer topologies.** Every J_{sq,hex}-covering sieve is J_{sq}-, J_{hex}- and J_{}-covering
    (a topology is finer when it has more covering sieves). For X = {squares}, {hexagons} and the empty set, O is
    separated but not a sheaf, so a(O) = O^+. The global sections are, respectively:"
  - *README line 23:* "THM-4891's relaxations are the sheafifications for finer topologies (fewer required 2-faces,
    more covering sieves)."
  - *README line 94:* "Sheafifying for finer topologies gives exactly THM-4891's relaxations".
  - *README line 96:* "So the gaps 7 and 21 hold for O itself, the sheaf for the coarsest of the four topologies.
    Each finer sheafification has sections with 7 and 21 sinks."
  - *INDEX rows:* "THM-4891 = sheafifications for finer face topologies".

**MAJOR-2 (D)4: the sheafified join count is false for non-strong factors.**
- **Counterexample.** A = C3 => C3 and B = C3 give |Q^+(A => B)| - 1 = 7, but the stated value is 9.
- **Why.** Q^+(A => B) = prod over all strong components of A and B. It equals (q_A + 1)(q_B + 1) iff each of A and
  B has at most one strong component of size >= 3 (equivalently, carrying an odd cycle).
- **Replace line 113.** "4. **Counting.** Put q = (H - 1)/2. For every order join A => B the presheaf count is
  |Q| - 1 = q_A + q_B + 2 q_A q_B (Session 210's join law). If A and B are strong (more generally, if each has at
  most one strong component carrying an odd cycle), the sheafified count is |Q^+| - 1 = q_A + q_B + q_A q_B, and the
  factor 2 on q_A q_B is the separation defect. In general |Q^+(A => B)| = |Q^+(A)| |Q^+(B)| = prod over the strong
  components C of A and B of (H(C) + 1)/2. For example A = C3 => C3, B = C3 gives |Q^+| - 1 = 7, not 9."
- **Proof lines 210-211.** Replace them with "For A, B strong, the components of A => B are A and B; the fibres over the
  mixed part (m = 2) have size 2, which gives the two counts."
- **Elsewhere.** In README line 90 and the results INDEX row, append "(A, B strong)".

**MAJOR-3 (C): "Hamiltonian paths as global sections" overclaims.**
- **What is proved.** The global sections are 2-coloured odd-cycle packings. OCF gives |Pk_2(V)| = H(T), but no
  bijection with Hamiltonian paths is constructed.
- **Why it matters.** Session 210's STRUCTURE_ATLAS section 1 explicitly lists such a bijection as an open
  obligation, and the README's own section 7 asks for "a site on which Hamiltonian paths themselves ... form a
  presheaf".
- **Replacements.**
  - *Line 97:* "**(C) H as the number of global sections.**"
  - *README line 22:* "H is the number of global sections of the sheaf of 2-coloured odd-cycle packings (OCF); no
    bijection with Hamiltonian paths is claimed."
  - *Hypotheses INDEX row:* "H = number of global sections of the 2-coloured odd-cycle packing sheaf".
  - *Title clause:* "So the multiplicativity of H is the sheaf condition" becomes "So the (classical)
    multiplicativity of H over strong components is the cardinality of the sheaf condition".

**MAJOR-4: typing and credits.** The status line types (A)-(F) as "PROVED (new formulation; elementary)", but
several parts are known:
- (C2) Rédei from OCF parity: Grinberg-Stanley 2023, whose abstract says their results recover the theorems of Rédei
  and Berge.
- (C3): THM-1862 and THM-1370 (1).
- (A): folklore. A tournament is a compatible family of pair orientations, transitive iff it has no cyclic triangle.
  The Condorcet reading is THM-512.
- (D) for m = 2: Session 210's relative-colour torsor.
- (F1)-(F2): Session 210 Theorem 1. The gate-cycle bijection in (F3) is Session 210 Theorem 2.

Replace status lines 4-12 with:
```
status: >
  PROVED (elementary; new as a sheaf-theoretic formulation): Lemma 0, (A1), (B1)-(B5) as site statements,
  (C1), (D1)-(D3), (D4) for strong factors, (E), the odd-subset gate formula of (F3).
  KNOWN, re-derived in site language:
  - (A2)-(A3) content: a tournament is a compatible family of pair orientations, transitive iff no cyclic triangle
    (classical); cyclic triangle = Condorcet paradox (THM-512);
  - (B5) H-invariance under componentwise converse (H multiplicative over strong components, H(C^op) = H(C));
  - (C2) Rédei's theorem from the parity of OCF (Grinberg-Stanley, arXiv:2307.05569, who derive Rédei and Berge;
    THM-002 notes);
  - (C3) the spectrum is the multiplicative closure of the strong spectrum (THM-1862; THM-1370 (1));
  - (D) for two components: Session 210 STRUCTURE_ATLAS section 2 and projective-cycle-carrier section 3 (the
    two-element relative-colour torsor); new here: the (Z/2)^(m-1) torsor for m components and the separation reading;
  - (F1), (F2): Session 210 Theorem 1; (F3)'s gate-cycle bijection: Session 210 Theorem 2.
  INHERITED and credited:
  - Grothendieck's plus construction: F^+ separated; F separated => F^+ a sheaf; a = (+)(+); F is a sheaf iff
    F -> F^+ is a bijection. (Mac Lane-Moerdijk III.5; SGA4 II.3.)
  - The comparison lemma, if UPGRADE-1 is adopted (SGA4 III.4.1; Johnstone, Elephant C2.2.3).
  - OCF H = I(Omega(T), 2) (Grinberg-Stanley; THM-002); Session 210's gate formula, Lemma 3, projective carrier,
    join law and coherent permutohedron; THM-4891, THM-1370, THM-4094, THM-1862, THM-512.
```

### MINOR (with replacement wording)

**MINOR-1. (B1) "points".** "Objects whose only covering sieve is the maximal one" is not the standard meaning of
"point"; a point of a topos is a geometric morphism Set -> E. Replace line 78 with: "**Topology and irreducible
objects.** J_ss is a Grothendieck topology. Its irreducible objects, those whose only covering sieve is the maximal
one, are the nonempty S with T[S] strongly connected. By the comparison lemma, Sh(P(V), J_ss) is equivalent to the
presheaves on the poset Str(T) of these S, so they are exactly the points of the topos (evaluation at S). The empty
set is covered by the empty sieve."
- Use "point" only in this sense.
  - Line 111 ("over a point whose packing"): use "over an element".
  - Line 132 ("a point on single vertices"): use "one-element on single vertices".

**MINOR-2. (B4) second bullet (line 88).** Replace with: "For every presheaf F, a(F)(S) = F^++(S) =
prod_i F(C_i), and F^+(S) -> F^++(S) is the inclusion of the fibre product over F(empty). It is a bijection unless
c(T[S]) >= 2, every F(C_i) is nonempty, and the images of F(C_i) -> F(empty) are not all one and the same
singleton."

**MINOR-3. (B5) (lines 93-95).** Replace with: "**Global symmetry becomes local symmetry.** The constant presheaf of
groups Z/2 sheafifies to a(Z/2)(S) = (Z/2)^(c(T[S])); two steps are needed exactly at the S with c(T[S]) >= 2. Let
Z/2 act on the presheaf Tour of all tournaments by the converse. The action commutes with restriction, so a(Z/2)
acts on a(Tour), where a(Tour)(V) = prod_i Tour(C_i) and (Z/2)^(c(T)) acts by componentwise converse. Reversing the
arcs inside any set of strong components of T keeps the components, the arcs between them, and H (classical)."

**MINOR-4. (C1) line 101.** Replace with: "The law H(A => B) = H(A) H(B) for order joins (classical; THM-1862) is,
via OCF, the cardinality of the sheaf condition for the covering sieve generated by {A, B}."

**MINOR-5. (D) wording.**
- Line 107: "**(D) The swap quotient Q of Pk_2 is not separated.**"
- Line 108: "- Q(V) is the set of cells of Session 210's projective carrier Q_T (STRUCTURE_ATLAS section 1's orbit set
  Q(T)) together with one basepoint, the empty packing, so |Q(S)| = (H(T[S]) + 1)/2."
- Title clause (D): "Session 210's projective carrier is not separated" becomes "the presheaf of swap-classes of
  2-coloured packings, whose global sections index the cells of Session 210's projective carrier, is not separated".

**MINOR-6. (E) proof (lines 214-219).**
- "Every pair inside a block is an edge of the face" should read "every pair inside a block is swapped by an edge of
  the face".
- Line 219 should read: "For X = {squares, hexagons}, THM-4891 (i), applied to each factor of F ~= prod_j Perm(B_j)
  (the squares mixing two blocks make the factors independent), identifies a(O)(F) with O(F)."
- For X = empty, cite Session 210 Theorem 4, or note that the square-coherent orientations are among all
  orientations.

**MINOR-7. (F) title and lines 131-134.**
- Title clause: "Padding adds two points to the site; the gate makes the whole vertex set a point (the site of
  Γ(T) has 2^n + 1 more points: {a}, {b} and {a, b} u S for S nonempty); OCF gives a second derivation of
  Session 210's gate formula".
- (F1): "... adds the two singleton points {b} and {a}, and no others ... For a family of sheaves defined uniformly
  for all tournaments with one-element values on single vertices (Pk_2 is one), the value at Σ(T) equals the value
  at T; so H(ΣT) = H(T)."
- (F2): "The gate Γ(T) (a -> b -> T -> a, T nonempty) is strong (Session 210, Theorem 1), so its whole vertex set
  is a point and a(A)(Γ(T)) = A."
- README lines 31 and 106 and the results INDEX row ("the gate collapses everything into one point" / "collapses to
  one"): "the gate makes the whole tournament a point".
- README line 120: "The gate makes every tournament the top point of a larger site: Γ(T) is strong and
  HC(ΓT) = H(T)."

**MINOR-8. (F3) credits and conventions.**
- Append to the first bullet of (F3): "This description of the gate cycles is Session 210's Theorem 2 (restricted
  to odd cycles)."
- In the displayed formula, state H(T[empty]) = 1 (the S = V term).

**MINOR-9. Background, line 58, and README line 64.**
- Line 58: "So two applications always suffice. F^+ is always separated (and may already be a sheaf); applied to a
  separated presheaf, the construction yields a sheaf."
- README line 64: "... but a compatible family of F^+-sections need not come from a single F-family: gluings may
  fail to exist one level further (uniqueness already holds, since F^+ is separated)."

**MINOR-10. README lines 24 and 77 (where "twice" lives).**
- The stated reason for one step on J_2 is incomplete. On J_ss the components also have no proper covers, yet two
  steps are needed.
- Line 77: "Here one step always suffices, for every presheaf, because the covering pieces (pairs) and their
  overlaps (singletons and the empty set) have only the trivial cover. On the strong-component site the overlap of
  two components is the empty set, which is covered by the empty sieve; that is exactly where a second step can act."
- Line 24: "On the strong-component site the second plus step can act only through the empty set, the one object
  covered by the empty sieve. On the face site the vertices play this role; on the pair site the second step never
  acts."

**MINOR-11. README section 0.4 and section 4 (the "+2" analogy).**
- Line 32: "mac-mini's THM-4900 shows that the gate acts on the marked-resolvent invariant y = (s(T) - 1)/2,
  s(T) = F_T(1), as the square of the continued-fraction map y -> 1 + 1/y."
- Line 34: "'Plus applied twice' and the '+2' operations share only the word and the number 2. No mechanism or
  functor relates them (verbal ANALOGY; cf. MISTAKE-529: equal numerals do not identify mechanisms)."
- Line 107 (THM-796 row). The last sentence is unsupported, since no coverage is defined and no step count is
  proved. Replace with: "PROVED there. As iterated gluing it has one level per added vertex; no site is defined,
  so it is not a plus construction (ANALOGY at most)."
  - Also, "plus one torsor bit" conflates two levels. At the tiling level the apex is a free bit (x F_2(apex)); the
    C_2 torsor arises at the line level.

**MINOR-12. README line 119.** It says "This is THM-4094's strong-lane reduction, now as the sheaf condition for
Pk_2", which overclaims. Replace with: "This multiplicative reduction is classical (THM-1862; THM-1370 (1)); the
sheaf condition for Pk_2 is its packing-level form. THM-4094's two-lane reduction uses more (the finite carry atoms
49, 63, 343) and is not re-derived here."

**MINOR-13. README section 1, completeness and accuracy.**
- Add tournament_category_23.py, euler_identity_S71w.py and functors_S71y.py to the S71 row, and THM-512 as the
  in-repository Condorcet precursor of (A).
- Replace the Sargsyan line with: "Sargsyan, *Localizing preference aggregation conflicts: a graph-theoretic
  approach using sheaves* (arXiv:2512.02416). Cellular sheaves on a graph of voters (stalks: total orders on each
  voter's visible alternatives; edge stalks on overlaps); Condorcet-type configurations give empty merged stalks
  ('local impossibilities'), i.e. no global section."
- HYP-1479 ("H = number of global sections of HP sheaf", CONFIRMED) and HYP-1511 ("H as sheaf global section",
  CONFIRMED) still carry CONFIRMED in `INDEX-HISTORICAL-THROUGH-2026-07-21.md`, although the README now shows they
  are ill-typed. Per MISTAKE-529's reusable rule (re-audit restriction maps), record the retype in a short MISTAKES
  entry, or annotate those rows.

**MINOR-14. Verification scripts overstate what they test.**
- `sheaf_check.py` check (1) computes |L^+| as `2 ** len(pairs)`, a formula, not matching families.
- The "torsor fibre" tests in `sheaf_check.py` and `q_torsor_check.py` compare `glob == loc * 2**(m-1)` with
  glob = 2^(k-1) and loc = 2^(k-m). That is an identity, true for every input.
- `q_sheaf` is defined by the formula prod (H(C_i)+1)/2, so the sheaf join law is checked against itself.
- Transitivity of J_ss and Lemma 0 are not tested at all.
- Only 2 of the 1,559 tournaments have two strong components carrying cycles (2 packings with m = 2 in all).
- Status line 18 ("(D) on 1,559 tournaments and on 80 order joins") should therefore read "(D): cardinalities on
  1,559 tournaments and 80 order joins; fibres, the (Z/2)^M action and both plus steps in the independent audit".
  Add the audit scripts to `verification:`.

**MINOR-15. INDEX rows.**
- Both rows inherit the wording of MAJOR-1/2/3 and MINOR-7.
- Hypotheses INDEX line 2: "H = number of global sections ..."; "THM-4891 = sheafifications for finer face
  topologies".
- Results INDEX line 6: "(join law ... -> q_A + q_B + q_A q_B for strong A, B)"; "finer topologies"; "the gate makes
  the whole tournament a point"; "Rédei = the unique fixed point of the colour swap (Grinberg-Stanley's parity
  consequence of OCF)".

**MINOR-16. Small proof gaps.**
- (A1): state that the canonical map F^+ -> F^++ is the identity under the identifications F^+(U) = F(U),
  |U| <= 2.
- (B3): add "F is a sheaf iff F -> F^+ is a bijection" to the inherited facts.
- (B4), bullet 1: note F^+(C_i) = F(C_i), since C_i is irreducible.

### UPGRADES (optional)

**UPGRADE-1. Comparison lemma.** Str(T), the poset of nonempty strong subsets, is J_ss-dense, and the induced
topology on it is trivial. So Sh(P(V), J_ss) is equivalent to [Str(T)^op, Set], and for every presheaf F:

    a(F)(S) = lim over strong C contained in S of F(C) = prod_i F(C_i).

- Verified on 8,410 (object, presheaf) pairs, 402 of the presheaf instances non-separated.
- Together with (B2), this formula gives (B3), (B4) and the constant-presheaf computation at once.
- It makes "points" literally correct: the points of a presheaf topos on a finite poset are the evaluations at its
  elements.
- Likewise Sh(J_2) is equivalent to presheaves on the subsets of size <= 2, which explains (A1). Sh(J_X) is
  equivalent to presheaves on {edges} u {X-faces}.

**UPGRADE-2. A precise form of "where twice lives".** Suppose every object has a smallest covering sieve generated
by irreducible objects. Then F^+ -> F^++ can fail to be onto only through overlaps of generators that are themselves
covered by smaller sieves:
- on J_ss, only the empty set (covered by the empty sieve);
- on J_X, the vertices;
- on J_2, none.

This replaces README section 0.3 and the first row of section 4, and yields the exact criterion of MINOR-2.

**UPGRADE-3. (B5) as an action.** Write it as in MINOR-3: a(Z/2) acts on a(Tour) = prod Tour(C_i) componentwise, and
(T_i) -> prod H(T_i) is invariant.

**UPGRADE-4. (D) in general.** |Q^+(A => B)| = |Q^+(A)| |Q^+(B)| for all A and B, and the separation defect over an
element whose packing meets m components is (Z/2)^(m-1). The torsor action itself, not only its size, is now verified
(`aud_ocf_q_gate.out`). Cite it in place of the tautological fibre test.

---

## 5. Recommended typing, part by part

| part | type |
|---|---|
| Lemma 0; (A1) | PROVED (elementary) |
| (A2), (A3) | PROVED, elementary/folklore. Only the J_2 formulation is new. Condorcet reading: THM-512. Sheaf-theoretic context: Abramsky-Brandenburger 2011, Sargsyan 2025. |
| (B1)-(B4) | PROVED (new formulation, elementary; an instance of the comparison lemma) |
| (B5) | PROVED (computation); the H-invariance is KNOWN (folklore) |
| (C1) | PROVED (OCF + condensation; the count is THM-002) |
| (C2) | KNOWN (Grinberg-Stanley) |
| (C3) | KNOWN (THM-1862, THM-1370 (1)) |
| (D1)-(D3) | PROVED; m = 2 is Session 210 |
| (D4) | presheaf count KNOWN (Session 210 (7)); sheafified count PROVED for strong factors, FALSE as stated |
| (E) | PROVED (reformulation of THM-4891 (i)-(iii); terminology corrected) |
| (F1), (F2) | KNOWN (Session 210 Theorem 1) plus PROVED site bookkeeping |
| (F3) | PROVED, new derivation. The odd-subset form is found nowhere earlier in the repository; it is equivalent to Session 210 (5) via Lemma 3, and the cycle bijection is Session 210 Theorem 2. |
| "plus twice" against "+2"; the THM-4900 golden map | ANALOGY (verbal only) |
| NOT CLAIMED block (topos/Lawvere, a functor gate = plus, Chan-Pak) | correct and adequate; keep |

---

## 6. Scripts and outputs (this folder)

- `aud_sites_colimit.py` and `.out`: J_2 and J_ss axioms, exhaustive for n <= 4.
  - The plus construction as the colimit over all covering sieves, with F^+, F^++, units and restriction maps.
  - Lemma 0; (A1)-(A3); (B1)-(B4); the constant presheaf, one and two steps, for n <= 6. About 145 s.
- `aud_ocf_q_gate.py` and `.out`:
  - OCF and the explicit Pk_2 sheaf bijection; the Rédei swap; the spectrum;
  - Q: orbits, fibres, the (Z/2)^M action, separatedness, and both plus steps for non-separated cases;
  - the (D)4 laws and the counterexample;
  - the gate formula, Session 210 (5), Lemma 3's alternating identity and the gate-cycle enumeration;
  - padding; (B5). About 30 s.
- `aud_perm_site.py` and `.out`:
  - J_X axioms and the order of the four topologies (the finer/coarser evidence);
  - O separated and sheaf, by colimit, on all faces;
  - global sections at n = 3, 4, 5; sink spectra at n = 4 (numpy for 2^24);
  - recount of THM-4891's n = 4 and n = 5 witness tables. About 150 s.
- `aud_gate_points.py` and `.out`: the number of points of the sites of T, Σ(T) and Γ(T), giving +2 and +(2^n + 1).

Run, from this folder: `python aud_sites_colimit.py > aud_sites_colimit.out` (likewise for the other three).
Requirements: Python 3, standard library, plus numpy for `aud_perm_site.py`.
