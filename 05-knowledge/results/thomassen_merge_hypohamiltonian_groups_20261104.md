# The Thomassen merge: the 8 Thomassen graphs and {6,4,8,9,12}; linear reducible-or-obstruction for uniquely Hamiltonian tournaments; hypohamiltonian graphs against non-critical runs and the square-sum graphs; the trivial group, C_2 and {3,7,21}; Maschke's planar groups and the fives; even graphs

**Session:** kind-pasteur-2026-11-04-S10 (continuation). Worktree `.claude/worktrees/uh-one-arc-join-7beb92`.

**Owner prompt (paraphrased).**
- Merge in {6,4,8,9,12}, the 8 Thomassen graphs and the prior work on 105.
- Consider the 2026 theorem that every planar triangulation contains linearly many pairwise non-touching reducible
  configurations or obstructing cycles of length at most 5.
- Consider tournament/even-graph equinumerosity and the bijection attempts.
- Consider how planar hypohamiltonian and hypotraceable graphs relate to the square-sum work.
- Consider planar groups and Maschke, the five rotation-group families {Z_n, D_n, A_4, S_4, A_5} against the repo's
  other fives, and the trivial group and C_2 against {7, 21} and other pairs.

**Status.**
- **PROVED.**
  - THM-4870 (A): the reducible-or-obstruction dichotomy and its exact mean.
  - THM-4870 (B): 3 is forced symmetric.
  - THM-4870 (C): no square-sum graph is hypohamiltonian or hypotraceable.
- **FINITE-EXACT.**
  - THM-4870 (B): the asymmetric spectrum for H <= 1125.
  - THM-4870 (C): the square-sum vertex-deletion data.
  - THM-4870 (D): the census of uniquely Hamiltonian graphs.
- **PRE-REGISTERED** (HYP-9462 in `ffedbba714`): (A1) CONFIRMED and (B1) CONFIRMED, but (B2) REFUTED, so conjecture (B)
  is false.
- **INTERPRETATION / NUMEROLOGY** (typed as such): the Maschke and fives readings, the coincidence ledger (section 7),
  and the 4CT analogy itself.
- **Novelty.** Not claimed against the literature. Thomassen 1980 (Discrete Math. 31) and the 2026 four-colour paper
  were read only through their abstracts.

Independently audited: section 10.

## 0. What the prompt's objects are

**The 8 Thomassen graphs** (MathWorld, "Thomassen Graphs").
- They are the graphs on 20, 24, 32, 34, 41, 60, 94 and 105 vertices.
- Their roles:
  - 34 is "the" hypotraceable Thomassen graph (Thomassen 1974);
  - 94 and 105 are planar;
  - 105 is Thomassen's first planar hypohamiltonian graph (1976).

  Today's records: the smallest planar hypohamiltonian graph known has 40 vertices (Jooyandeh, McKay, Östergård,
  Pettersson, Zamfirescu), and the smallest planar hypotraceable graph known has 154.
- **{6,4,8,9,12}** is the set of distinct MathWorld crossing numbers of the six non-planar Thomassen graphs: 20 has 6,
  24 has 4, 32 and 34 both have 8, 41 has 9, 60 has 12 (so 8 occurs twice). In MathWorld's listing order the values
  read 6, 4, 8, (8), 9, 12.

**Coincidence in the repository** (inheritance search; NUMEROLOGY).
- The LRC radius-7 handoff `05-knowledge/results/radius7_multicluster_handoff_opus_S314.txt` lists the speed prefix
  [4, 6, 8, 9, 12] with sixth-speed lifts including 105 = 14 + 7·13. The audit found that the file lists all 792
  five-element subsets of {1..12}, that 105 is a lift for 462 of them, and that 14 + 13j is only one of seven lift
  residue classes. So the co-occurrence is forced by the file's design.
- {4,6,8,9,12} is exactly the set of composite {2,3}-smooth numbers <= 12. THM-806 (16) has
  {ab mod 13 : a, b in {2,3,4}} = {3,4,6,8,9,12}.
- The joint appearance with 105 therefore carries no evidential weight. No mechanism connects crossing numbers to
  lonely-runner speeds.

**The triangulation theorem.** Inoue, Kawarabayashi, Miyashita, Mohar, Thomassen, Thorup, *The Four Color Theorem with
Linearly Many Reducible Configurations and Near-Linear Time Coloring*, arXiv:2603.24880 (v3, 1 Oct 2026).
- Main theorem: every planar triangulation contains linearly many pairwise non-touching reducible configurations, or
  pairwise non-crossing obstructing cycles of length <= 5.
- Consequence: an O(n log n) 4-colouring.
- Read from the abstract and listings only (CITED). The figure of more than 8,200 D-reducible configurations comes from
  the paper's text as quoted in secondary summaries, not from the abstract.

**Thomassen's tournament paper.** *On the number of Hamiltonian cycles in tournaments*, Discrete Math. 31 (1980)
315-323. It is not cited anywhere in the repository, and its body could not be read here (CAPTCHA).

**Equinumerosity.** Royle, Praeger, Glasby, Freedman, Devillers (J. Algebraic Combin. 2023; arXiv:2204.01947) prove
that tournaments and even graphs are equinumerous for every n. The repository's programme on this is extensive:
- THM-4720 / THM-4721 (the formula);
- HYP-9340 (refined by odd symmetry: proved for odd p-groups, refuted at six mixed-prime symmetries up to n = 13, with
  the double-triple 5+3+3 boundary);
- seven bijection attempts with proved obstructions (THM-4531, THM-4722, THM-4725);
- the merge packet `even_graph_response_merge_20261104`.

## 1. Linear reducible-or-obstruction for uniquely Hamiltonian tournaments (THM-4870 (A); PROVED)

**The dictionary.**
- In a triangulation, a *reducible configuration* lets one pass to a smaller instance, and an *obstructing cycle* is a
  short separating cycle.
- In a uniquely Hamiltonian tournament T, a non-critical vertex v is a reduction, because T - v is again uniquely
  Hamiltonian (Lemma 1 of THM-4868).
- A cut arc of the cycle is an obstruction: the one-arc cut it crosses splits T as a one-arc join (THM-4867). These
  cuts are nested (S74 audit, section A).

**Theorem.** Let r be the number of non-critical vertices and k the number of cut arcs.
- (i) r + k = n + 1 (THM-4868 (D)). So every uniquely Hamiltonian tournament has at least (n + 1)/2 reductions or at
  least (n + 1)/2 obstructions. The reductions form one run of the cycle and the obstructions form the complementary
  path.
- (ii) Both extremes occur. A_n has r = 2 and k = n - 1. Classes with r = n - 1 and k = 2 number 2^(n-5).
- (iii) Over the F(2n - 6) classes of order n, E[r] = ((5 + √5) n + 12 - 8√5) / 10 + O(n φ^(-4n)) and
  E[k] = n + 1 - E[r]. The error is (72 + 161√5/5)(n - 8) φ^(-4n) + ..., per the audit.
  The typical fractions are therefore φ/√5 ≈ 0.7236 reductions and (5 - √5)/10 ≈ 0.2764 obstructions.

**Proof of (iii).**
- The generating function of THM-4868 (H) gives E[r] = [x^n] ∂_t G(x, 1) / [x^n] G(x, 1), with:
  - G(x, 1) = x^4 / (1 - 3x + x^2), which has a simple pole at φ^(-2);
  - ∂_t G(x, 1) = x^4 (x^3 - 2x^2 - 3x + 2) / (1 - 3x + x^2)^2, which has a double pole there.
- Partial fractions give a mean that is linear in n up to O(n φ^(-4n)). The slope is
  -ρ'(1)/ρ(1) = (5 + √5)/10, where ρ(t) is the dominant root of 1 - (1 + 2t)x + t x^2.
- The constant (6 - 4√5)/5 is checked exactly against the rows for n = 10-70, where the deviation is below 10^(-13)
  from n = 20 on.

**The analogy (INTERPRETATION).** Both theorems say the same kind of thing: linearly many local reductions, or
linearly many nested obstructions. Two differences:
- the 2026 configurations are pairwise non-touching, while here the reductions form one run of adjacent vertices;
- here the trade-off is an identity (r + k = n + 1), not a lower bound.

## 2. Hypohamiltonian graphs against non-critical vertices

**Tournaments.** No tournament is hypohamiltonian or hypotraceable: every tournament has a Hamiltonian path (Rédei),
and a non-strong tournament has non-strong vertex deletions. The right analogue is the non-critical set:
- in a uniquely Hamiltonian tournament it is ONE run of length 2..n-1, never all of V (THM-4868);
- a hypohamiltonian graph has every vertex non-critical while the graph itself is non-Hamiltonian.

**Uniquely Hamiltonian graphs (HYP-9462 (B); REFUTED).** The conjecture was the graph-side opposite of the tournament
run: never two consecutive non-critical vertices on the cycle. Results:
- it holds exhaustively at n <= 9 (labelled counts 3, 11, 69, 561, 6,903, 113,107);
- it fails first at n = 10, in 1,300 of 2,689,657 labelled graphs;
- the pre-registered random test at n = 12-15 found 5 witnesses (HYP-9462 Results).

So neither direction transfers. The non-critical set of a graph is usually scattered but can contain consecutive
vertices. The proved partial case: two "shortcut" deletions at consecutive vertices force a second Hamiltonian cycle
(HYP-9462 (B)).

**Square-sum graphs (THM-4870 (C)).** The inheritance comes from THM-4505, THM-4851 and the planarity fact recorded in
`collatz_procgen_20260922_synthesis.md` ("Q_N is planar iff N <= 24", with the first K_(3,3) through the ear 11-25-24).
New here:
- **(C1)** No Q_N with N >= 3 is hypohamiltonian or hypotraceable. (Q_2, two isolated vertices, is trivially
  hypotraceable.)
  - For 3 <= N <= 30, Q_N has a vertex of degree <= 1: degree 0 for N <= 6, minimum degree 1 for 7 <= N <= 30
    (THM-4505).
  - Q_31 has minimum degree 2 and is traceable. For N >= 32, Q_N is Hamiltonian (THM-4851 (a)). A hypohamiltonian graph cannot
  have a vertex of degree <= 2, and a hypotraceable graph cannot have a leaf (delete its neighbour).
- **(C2)** The vertex-deleted subgraphs of the non-traceable orders, which the repository had not studied:

  | N | 18 | 19 | 20 | 21 | 22 | 24 |
  |---|---:|---:|---:|---:|---:|---:|
  | vertices v with Q_N - v traceable | 1 | 1 | 1 | 1 | 2 | 4 |

  For Q_23, Q_25 and Q_26 (traceable, not Hamiltonian), Q_N - v is Hamiltonian only for v = 18, the leaf.
  Q_31 - v is Hamiltonian exactly for v in {2, 4, 13, 15, 26, 28, 29}.
- **(C3)** In the uniquely Hamiltonian Q_32 and Q_33, the non-critical vertices are:
  - Q_32: {2, 13, 15, 26, 29, 30} (6);
  - Q_33: {2, 13, 15, 18, 26, 29, 30, 32, 33} (9).

  Each set is pairwise non-adjacent in the graph, and no two of its vertices are consecutive on the unique cycle. This
  is the "pairwise non-touching" pattern of the triangulation theorem, and the opposite of the tournament run. It is
  not a law for uniquely Hamiltonian graphs in general (the refutation above).
- **Planar and hypohamiltonian.** The planar square-sum graphs (N <= 24) lie inside the regime with vertices of degree
  <= 1 (N <= 30), and the first Hamiltonian one (Q_32) is non-planar. So the square-sum family meets neither
  Thomassen's planar hypohamiltonian graphs (40-105 vertices and up) nor planar uniquely Hamiltonian ones.

## 3. The trivial group, C_2 and the pairs

**(B) of THM-4870: the asymmetric spectrum.**
- **PROVED.** H(T) = 3 forces a 3-cycle component, since strong tournaments of order >= 4 have H >= 5. So Aut(T)
  contains Z_3.
- **FINITE-EXACT.** For every odd H <= 1125, some tournament with trivial automorphism group has H Hamiltonian paths iff
  H is not in {3, 7, 21}. HYP-9462 (A1) was pre-registered on [225, 1125].
  - Below 225 the census is complete: strong order >= 11 has H >= f(11) = 225 (Busch 2006; THM-4744). Every value is
    realised by products of asymmetric strong classes of order <= 10 and of uniquely Hamiltonian tournaments, which
    are all asymmetric (THM-4868 (E)).
  - The audit independently exhibited an explicit asymmetric tournament of order <= 9 for every odd H <= 1125 except
    3, 7 and 21 (560 witnesses).
  - Across orders <= 10 the only strong value needing symmetry is 3. Within smaller order caps more do (171, 175, 189
    within orders <= 7). The |Aut| values seen are 1, 3, 5, 7, 9, 15 and 21.
  - The products of uniquely Hamiltonian H values alone miss 24 odd numbers up to 1125 (HYP-9462).
- **Reading.** The trivial group sees exactly the two crevasses plus one forced value, 3 = |Z_3|. 21 = 3·7 is the
  product of the other two.

**C_2 never acts on a tournament, but appears beside it in three places.**
- Aut(T) has odd order: an involution would reverse an arc. It acts freely on Hamiltonian paths, so |Aut| | H
  (THM-4747), and H / |Aut| is never 7 or 21 (THM-4812).
- *Converse.* The converse gives Aut±(T) / Aut(T) in {1, C_2}. The C_2 case is the self-converse "real points"
  (THM-4866 (E)), whose spectrum m·s^2 has the parity pair {3, 5}: 5 is impossible at odd order and 3 at even order.
  THM-4743 (2), Gluck's theorem with the converse added, shows that the primitive tournaments in which no vertex
  subset has trivial stabilizer under Aut± are exactly C_3, Z_5{1,2} and P_7.
- *Graphs.* The sign character χ_G : Aut(G) → C_2 (an automorphism's parity of reversed edges) has image 1 (even
  graph) or C_2 (odd graph). Even graphs and tournaments are equinumerous, and #graphs = #tournaments + #odd graphs.
- So on the tournament side the trivial group/C_2 dichotomy is the converse quotient, and on the graph side it is the
  sign character. The counting identity matches labelled tournaments with the C_2-trivial graphs under a
  Cauchy-Frobenius sum whose even permutations cancel (THM-4720).

**{7, 21} as group orders (KNOWN numerology).** P_7 has Aut = F_21, of order 21, containing the translations Z_7. It
acts freely on its 189 = 9·21 paths, and its orientation-plus-converse group is AGL(1,7)
(`monster_levels_square_sum_7_21_20261022.md` section 5). No mechanism relates this to the crevasses.

**Other pairs.**
- **{32, 33}:** the uniquely Hamiltonian square-sum graphs, one segment reversal apart (THM-4851).
- **{12875, 12877}:** klein's rearming targets (THM-4862).
- **Today's binary structures (THM-4868):**
  - 1 or 2 sole exits, with 1 or 2 sole entries;
  - 2 or 3 admissible Shelor parents;
  - the counts F(2n - 8) / F(2n - 7) for one or two sole exits.

## 4. Maschke's planar groups and the fives (INTERPRETATION unless stated)

**Maschke's theorem** (1896; not previously in the repository). The finite groups with a planar Cayley graph are
exactly the finite groups of isometries of the 2-sphere. The rotation groups among them are Z_n, D_n, A_4, S_4 and A_5.
The repository uses "Maschke" only for semisimplicity.

**PROVED (elementary).**
- A tournament automorphism group has odd order. The only odd-order groups among Maschke's rotation groups are the
  cyclic ones (D_n, A_4, S_4 and A_5 all have even order).
- So a tournament's automorphism group is "planar" iff it is cyclic.
  - The first non-cyclic one is Z_3 × Z_3, at order 6 (C_3 ⇒ C_3, H = 9).
  - The first strong non-cyclic ones are at order 7, and the census shows exactly two (FINITE-EXACT):
    - Z_3 × Z_3 (H = 99, hc = 9, self-converse). This is C_3 with two of its vertices replaced by 3-cycles: a group of
      order 9 on 7 points cannot be Z_9;
    - F_21, for P_7 (H = 189).
  - At order 6 every strong class has |Aut| 1 or 3.

**Converse groups of the three primitive exceptional tournaments of THM-4743 (2).**
- Aut±(C_3) = S_3 = D_3 and Aut±(Z_5{1,2}) = D_5: Maschke's dihedral family, the rotation groups of the triangular and
  pentagonal prisms.
- Aut±(P_7) = AGL(1,7), of order 42, is not a planar group.

**The fives.**
- The rotation families split as 2 infinite (Z_n, D_n) plus 3 exceptional (A_4, S_4, A_5).
- Shelor's five letters split as 2 (type A: a, b) plus 3 (type B: c, d, e). These are the admissible parents of
  THM-4868 (G), with |adm| = 2 or 3.
- ADE splits the same way (A_n, D_n; E_6, E_7, E_8). This one is NOT a coincidence: it is the McKay correspondence
  between the finite subgroups of SU(2) (the binary versions of the rotation groups) and the simply-laced Dynkin
  diagrams.
- The Fermat five and the Platonic five share only {3, 5}
  (`gilbreath_fermat_platonic_20260926.md`).
- `fermat_five_pepin_collatz_20261025/FIVE_LADDER.md` types the repository's fives as cardinality coincidences.
  Nothing here changes that typing for the repository's fives. The rotation-group/ADE match is the one genuine
  correspondence (McKay).

**One structural echo** (INTERPRETATION). The icosahedral group A_5 is the rotation group whose character table needs
√5. φ^(±2) are the eigenvalues of Shelor's transfer matrix [[1,1],[1,2]], and φ/√5 is the reduction fraction of
section 1. No map is claimed.

## 5. Equinumerosity: what today's objects add

- **Rigidity breaks equivariance on the tournament side.** Uniquely Hamiltonian tournaments are rigid (Aut = 1) from
  n = 4 on (C_3 has Aut = Z_3), and there are F(2n - 6) of them. There are no asymmetric graphs on <= 5 vertices. So no bijection
  preserves trivial automorphism groups. This is the repository's attempt-1 obstruction (|Aut| multisets) seen from
  the rigid side.
- **The {1, C_2} dichotomy on both sides** (section 3) suggests one invariant to carry: the converse quotient on
  tournaments against the sign character on graphs. This is a question, not a construction.
- **The canonical Shelor tree** of THM-4868 (G) gives the uniquely Hamiltonian classes a canonical recursive name, a
  path of letters. Any even-graph family counted by F(2n - 6) with a 2/3-branching grammar would be a candidate target.
  None was found in this session; this was not searched systematically.

## 6. Answers to the prompt, in one line each

- **{6,4,8,9,12}:** the distinct crossing numbers of the six non-planar Thomassen graphs (8 occurs twice). Their
  co-occurrence with 105 in the LRC handoff is forced by that file's design.
- **The 8 Thomassen graphs and 105:** the hypo-family whose tournament analogue is the non-critical run, and whose
  square-sum analogue fails (THM-4870 (C)).
- **Linearly many reducible configurations:** the tournament analogue (INTERPRETATION) is r + k = n + 1, with mean
  fraction φ/√5 (THM-4870 (A)).
- **Equinumerosity:** rigid uniquely Hamiltonian tournaments and the {1, C_2} dichotomies (section 5). No bijection.
- **Hypohamiltonian against square-sum:** no square-sum graph Q_N with N >= 3 is hypo-anything. The vertex-deletion data are new, and
  the non-touching non-critical sets in Q_32 and Q_33 are an instance only (section 2).
- **Maschke and the fives:** tournament symmetry is planar iff cyclic. The fives are cardinality coincidences
  (section 4).
- **The trivial group, C_2 and {7, 21}:** the asymmetric spectrum misses exactly {3, 7, 21} up to 1125. C_2 enters as
  the converse and the sign character (section 3).

## 7. Coincidence ledger (NUMEROLOGY; hostile probe)

The Thomassen orders {20, 24, 32, 34, 41, 60, 94, 105} against the square-sum landmarks:

| landmark | meaning |
|---:|---|
| 24 | last planar Q_N |
| 32 | first Hamiltonian Q_N, unique cycle |
| 34 | first terminal cycle |
| 41 | Gerbicz pairs start; the Q_41 gate |

Four of the five Thomassen orders in [20, 50] are landmarks.
- Only 12 of the 15 landmarks {15, 16, 17, 23, 24, 25, 31, 32, 33, 34, 38, 39, 40, 41, 42} lie in [20, 50] (31
  integers).
- So the expected number of chance hits is about 1.94, and P(>= 4 hits) is about 0.06 (corrected by the audit).
- The landmark set was chosen after seeing the orders, so even 0.06 overstates the evidence. This stays NUMEROLOGY,
  and no mechanism is proposed.

## 8. Typing

| item | type |
|---|---|
| (A) r + k = n + 1 dichotomy; mean E[r] with error O(n φ^(-4n)) | PROVED (THM-4868 (D), (H), and singularity analysis); exact check n <= 70; error term corrected by the audit |
| (B) 3 forces Z_3 | PROVED |
| (B) asymmetric spectrum = odd numbers minus {3, 7, 21}, for H <= 1125 | FINITE-EXACT (A1 pre-registered); all H: HYP-9462 (A), OPEN |
| (C1) no Q_N with N >= 3 is hypohamiltonian or hypotraceable | PROVED (from THM-4505, THM-4851) |
| (C2), (C3) square-sum vertex deletions | FINITE-EXACT |
| (D) uniquely Hamiltonian graphs: consecutive non-critical pairs absent at n <= 9, present at n = 10 (1,300) | FINITE-EXACT (exhaustive); conjecture (B) REFUTED |
| Maschke: tournament symmetry planar iff cyclic | PROVED (elementary) |
| 4CT analogy, fives, A_5/φ echo, Thomassen/square-sum coincidences, {4,6,8,9,12}+105 | INTERPRETATION / NUMEROLOGY |

## 9. Pre-registration record

- HYP-9462 was pre-registered in `ffedbba714`, with the programs and the exploratory record
  (`thomassen_merge_20261104_explore.out`).
- Confirmatory outputs: `thomassen_merge_20261104_confirm.out`.
- The down-set census (`uh_graph_downset_20261104.c`) was written after the refutation, to find the minimal order.

## 10. Audit

An independent adversarial audit ([report](../../04-computation/experiments/thomassen_merge_20261104_audit_S10/AUDIT.md))
wrote its own code. Verdict: **SOUND WITH CORRECTIONS**, with no blocking finding.

**Reproduced.** Every finite claim in (A)-(D), including:
- its own uniquely Hamiltonian census (n <= 11), with exact means;
- 560 explicit asymmetric witnesses for all odd H <= 1125 except 3, 7, 21;
- its own strong-class censuses (2,643 classes with H <= 224; 107,773 with H <= 1125 at orders <= 10);
- the square-sum data, and the uniquely Hamiltonian graph counts by two methods;
- the two non-cyclic groups at order 7, and the citations.

**Corrected (MISTAKE-708).**
- the error term (O(n φ^(-4n)), not O(φ^(-4n)));
- "exact analogue", retyped as interpretation, with the differences stated;
- six non-planar Thomassen graphs, not five (32 and 34 both have crossing number 8);
- the coincidence-ledger arithmetic;
- HYP-9462: chord-adjacent non-critical pairs first appear at n = 8 (16 graphs), not n = 9;
- (C1) for small N;
- rigid from n = 4, not n = 3;
- the THM-4743 sentence;
- the LRC file description;
- the Busch citation for f(11) = 225;
- smaller wording: the McKay correspondence named, "lies inside" for the planar regime, the source of the 8,200 figure,
  and the order-cap qualification of "only 3 needs symmetry".

## 11. Reproduction

```
gcc -O2 -fopenmp -o as  04-computation/experiments/asym_spectrum_autspec_20261104.c
./as autspec 10 1125 > autspec.out                 # strong classes of order <= 10, H <= 1125, with |Aut|
gcc -O2 -fopenmp -o uhh 04-computation/experiments/uh_hpath_values_20261104.c
./uhh 14 0 0 > uhH.out                             # H of every uniquely Hamiltonian class, orders 4-14
python3 04-computation/experiments/asym_spectrum_closure_20261104.py autspec.out uhH.out 1125
gcc -O2 -o g 04-computation/experiments/uh_graph_noncritical_20261104.c
./g 8 8 0 1 0.5 ; ./g 12 15 750 20261104 0.12 0.18 0.24 0.30
gcc -O2 -fopenmp -o d 04-computation/experiments/uh_graph_downset_20261104.c ; ./d 4 10
gcc -O2 -o q 04-computation/experiments/square_sum_vertex_deleted_20261104.c ; ./q 14 40 100000
```
