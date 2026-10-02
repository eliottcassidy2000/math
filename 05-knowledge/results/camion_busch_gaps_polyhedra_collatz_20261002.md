# Camion, Busch and the two tournament gaps; tournaments as polyhedra (A₃, the tetrahedron, P7 on the torus); why "{7, 21} are the only gaps" cannot be isomorphic to Collatz

**opus S15, 2026-10-02.**

- Scripts: `04-computation/experiments/camion_busch_gaps_polyhedra_collatz_20261002.py` (+ `.out`, ALL CHECKS PASSED) and
  the order-7 census `camion_busch_strong7_20261002.c` (+ `.out`).
- Related: [`noble_polyhedra_hill_fissary_20261002.md`](noble_polyhedra_hill_fissary_20261002.md) (gap 7 = Camion +
  Moon; the A₅ dictionary); [`procgen_tcpc_20261001_tournament_clock_prime_collatz.md`](procgen_tcpc_20261001_tournament_clock_prime_collatz.md)
  (the triplet principle, `F = U + S`, the 3-4-3 split of `p²qr`).
- Canon: THM-075, THM-079, THM-115, THM-200 and THM-1370 (the two gaps); THM-466 (2-adic digits of H); THM-4094, with
  THM-4097, THM-4102 and THM-4104 (completeness).

**Status.**
- **The gap at 21, as a theorem of Camion's type.** The two gaps sit in the crevasses of the strong floor `f(m)`, the
  minimum H over strong m-tournaments (THM-1370): 7 ∈ (f(4), f(5)) = (5, 9) and 21 ∈ (f(6), f(7)) = (15, 25).
  - PROVED + COMPUTED: Camion's theorem with Moon's cycle counts gives the direct count `f(m) ≥ B*(m)`. Moon's counts
    are c_k ≥ m − k + 1 (Moon 1966, Cor. 2.1, sharp), and B*(m) = 1 + 2⌊(m−1)²/4⌋ = 3, 5, 9, 13, 19, 25, 33 for
    m = 3, …, 9.
    - **B*(5) = f(5) = 9, so this count certifies the gap 7 exactly.**
    - **B*(7) = 19 does not exceed 21, so it does not certify the gap 21 at m = 7.** B*(8) = 25 already does.
    - A floor proof at m = 7 needs f(7) ≥ 23. The true floor is f(7) = 25, which is Busch's theorem (2006, Cor. 1) at
      m = 7. It sharpens the bound H ≥ m that Camion's theorem gives. The census here re-derives f(7) = 25.
    - What the count misses at m = 7 is the disjoint pairs. The tournament with the fewest odd cycles (9, exactly
      Moon's bound) is Moon's T′₇, which has three disjoint pairs and H = 31. The floor H = 25 is one isomorphism class
      with 12 pairwise-intersecting odd cycles.
  - COMPUTED (exhaustive): the birth of 21 at n = 6 is a **mod-4 lock**.
    - Every strong 6-tournament with exactly five cyclic triangles has an even number of 5-cycles (4 or 6) and at most
      one disjoint pair of triangles. So its H is 19, 23 or 27, all 3 (mod 4).
    - Four triangles give H ≤ 17; six or more give H ≥ 23.
    - The window around 7 at n = 5 is locked to the other residue: H ∈ {5, 9}, both 1 (mod 4).
  - PROVED (elementary; canon THM-075, THM-079, THM-115 and THM-200 have the case lists).
    - H = 7 iff the odd-cycle conflict graph is K₃.
    - H = 21 iff it is one of K₁₀, K₈ − e, K₆ − 2K₂, K₆ − P₃, P₄, K₃ + K₁.
    - Since I(P₄, x) = (1 + x)(1 + 3x), the 7-polynomial 1 + 3x is a factor of one of the four 21-polynomials: the
      21-problem contains the 7-problem.
- **Collatz versus "{7, 21} are the only gaps": no uniform isomorphism is possible.**
  - PROVED (elementary, from Camion's theorem): an odd m is an H-value iff some tournament on at most m vertices has
    H = m. So the completeness conjecture is a **Π₁** sentence: a counterexample would be certified by a finite search.
  - CITED: Collatz is Π₂, and the generalized Collatz problem is Π₂-complete (Kurtz–Simon 2007). Only its cycle half
    ("no cycle but 1 → 4 → 2 → 1") is Π₁. The divergence half has no tournament counterpart, because Camion's theorem
    bounds every witness.
  - NUMEROLOGY: {1, 2, 4} is both the trivial Collatz cycle and the set of quadratic residues mod 7 that defines P7.
    7 = 1 + 2 + 4 and 21 = 1 + 4 + 16 = 4² + 4 + 1 are the point counts of PG(2,2) and PG(2,4). In Syracuse terms 7 ≡ 3
    (mod 4) rises (7 → 11), while 21 ≡ 1 (mod 4) drops straight to 1.
- **Tournaments as polyhedra: transitivity, duality and the medial graph** (PROVED or COMPUTED; classical where cited).
  - Score vectors are the weights of `V_ρ`: `Σ_T x^score = Π (x_i + x_j) = s_δ`.
  - At n = 4 the score vectors form three shells around the centre:
    - 24 transitive (the truncated octahedron, multiplicity 1);
    - 8 with one cyclic triangle (a cube, multiplicity 2);
    - 6 strong (an octahedron at the cube's face centres, multiplicity 4).
  - Arc reversals are the 12 roots, i.e. the **cuboctahedron**, the medial graph of cube and octahedron. The score
    lattice is A₃ (face-centred cubic), whose Voronoi cell is the **rhombic dodecahedron**.
  - On the tetrahedron, planar duality of orientations sends cyclic triangles to sources and sinks. It swaps the 24
    transitive 4-tournaments with the 24 strong ones and fixes the 16 others as a class.
  - **P7 is the face-cyclic orientation of the 7-vertex torus triangulation.** Its 14 cyclic triangles are exactly the
    14 faces, which form two Fano planes. The only two such orientations are P7 and its converse.
    - Under the subgroup Aut(P7) (order 21) the map has the quasiregular orbit pattern: one orbit on vertices, one on
      edges, two alternating orbits on faces.
    - Its dual, the Heawood graph oriented points → lines, has the quasiregular-dual pattern, as the rhombic
      dodecahedron does.
    - Its medial graph is a Cayley graph of the Frobenius group of order 21.
    - (V, E, F) = (7, 21, 14): the two gap values are V and E. That is NUMEROLOGY.
  - Beyond 7 (COMPUTED; independent check pending, §6.6): P31 is the face-cyclic orientation of 1,280 Z₃₁-invariant
    triangulations of the genus-63 surface (112 up to multipliers). P19 and P43 have none of this kind.
- The triplet principle (TCPC: an involution with one fixed point plus a factor-2 cocycle). §7 records exact facts that
  carry the involution half. None of them supplies the factor-2 cocycle; for 60, TCPC §4.3 shows there is none between
  3 and 5. **As instances of the principle they are ANALOGY.**
  - The sharpest is **60 = 2²·3·5 = |A₅|**: its ten proper divisors are exactly the stabilizer orders and orbit sizes of
    the icosahedral rotation group (exact).
  - Reading the square 2² as the doubling of stabilizers from points to axes is ANALOGY.
- Independent audit: DONE (2026-10-02) for §1–§7. No item was unsound. The corrections are applied (§11; MISTAKE-560).
  The §6.6 counts are under a separate check.

## 0. The owner's prompt

The gap at 7 is Camion's theorem; find the theorem that plays this role for 21, and how the two are linked, possibly as
axioms with many consequences. Use the medial graph and the graph dual, and how transitivity takes us from tournaments
to polyhedra and to equations in three variables. Consider attempting to prove that the Collatz conjecture is
isomorphic to "{7, 21} are the only forbidden H-values". Apply the triplet principle (one member primarily
differentiated, the other two subtly, by halving/doubling) and the earlier `F = U + S` / `p²qr` work (3-4-3 split; "`p²qr`
is `p` and `p³` combined"). Draw inspiration from the rhombic dodecahedron. The owner attached a text on
quasiregular duals: the triality G, G*, M(G), whose vertices are the vertices, faces and edges of a polyhedron.

## 1. The two gaps as forbidden conflict graphs (§A)

By the odd-cycle formula (OCF, THM-002), `H(T) = I(Ω(T), 2) = 1 + 2α₁ + 4α₂ + …`. Here `Ω(T)` is the conflict graph of
the directed odd cycles and `α_k` counts k pairwise disjoint ones. Solving `α₁ + 2α₂ + 4α₃ = (H − 1)/2` over graphs
gives:

| H | α-vector | Ω | I(Ω, x) |
|---|---|---|---|
| 7 | (3, 0, 0) | K₃ | 1 + 3x |
| 21 | (10, 0, 0) | K₁₀ | 1 + 10x |
| 21 | (8, 1, 0) | K₈ − e | 1 + 8x + x² |
| 21 | (6, 2, 0) | K₆ − 2K₂, K₆ − P₃ | 1 + 6x + 2x² |
| 21 | (4, 3, 0) | P₄, K₃ + K₁ | (1 + x)(1 + 3x) |

Every other α-vector has no graph. This is canon's starting point (THM-200 for 7; THM-075, THM-079 and THM-115 for 21).
The new remark is the factorisation in the last row: the 7-polynomial `1 + 3x` divides a 21-polynomial. P₄ and K₃ + K₁
have the same independence polynomial, and K₃ + K₁ is literally "a 3-cycle next to an impossible K₃". So the
21-problem contains the 7-problem. THM-079 excludes K₃ + K₁ by the 7-gap, and P₄ by a direct argument on its middle
pair of cycles. (The audit notes that this argument, as written, does not treat two triangles sharing an edge. The
canon result stands through THM-1370.)

## 2. Camion for 7, Busch for 21 (§B, C census)

Write `f(m)` for the strong floor, the minimum H over strong m-tournaments. It is monotone (THM-1370, from Moon's
subtournament theorem and the insertion lemma). Two classical theorems bound it from below:
- **Camion** (1959): a strong tournament has a Hamiltonian cycle. Cutting it at each of its m arcs gives m distinct
  Hamiltonian paths, so `f(m) ≥ m`.
- **Moon** (1966, Thm 1, Thm 2, Cor. 2.1): every vertex of a strong m-tournament lies on a cycle of every length
  3, …, m, and there are at least `m − k + 1` directed k-cycles. This count is sharp: Moon's tournament T′_m
  (p_i → p_j iff i = j − 1 or i ≥ j + 2) has exactly m − k + 1 k-cycles for every k (checked for m = 5, 6, 7).

With the OCF and α₂ ≥ 0 these give the direct count

`f(m) ≥ B*(m) = 1 + 2 Σ_{odd k ≤ m} (m − k + 1) = 1 + 2⌊(m − 1)²/4⌋`.

The weaker pancyclicity count of THM-115 (c₃ ≥ m − 2, c_k ≥ ⌈m/k⌉) is listed as B(m).

| m | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| direct count B*(m) (Moon's c_k ≥ m − k + 1) | 3 | 5 | **9** | 13 | **19** | 25 | 33 |
| pancyclicity count B(m) (THM-115) | 3 | 5 | 9 | 13 | 17 | 21 | 25 |
| strong floor f(m) | 3 | 5 | **9** | 15 | **25** | 45 | 75 |

The floors are exhaustive for m ≤ 6 (§B) and for m = 7 (the C census: 1,677,488 strong labelled 7-tournaments, with the
OCF checked on each). f(8) and f(9) are from canon and Busch.

- **Gap 7.** f(5) = B*(5) = 9: the direct count pins the floor exactly at m = 5. Below 9 only the strong values 3 and 5
  and the product 9 exist, so 7 is skipped for every n. On five vertices this is the familiar form: strong ⟺ c₃ ≥ 3 ⟺
  Hamiltonian cycle, and H jumps from 5 to 9.
- **Gap 21.** It needs f(m) > 21 for every m ≥ 7. The direct count gives B*(7) = 19 at m = 7 and B*(8) = 25 at m = 8,
  so **m = 7 is the one size where it falls short**. (THM-115's weaker count gives B(7) = 17 and B(8) = 21, and first
  exceeds 21 at m = 9.) This is a statement about the direct count, not a proof that no argument from Camion's and
  Moon's theorems could reach f(7) = 25.
- **The theorem that plays Camion's role for 21** is the next rung of the same ladder: **Busch's theorem** (Electron.
  J. Combin. 13 (2006) #N3, Cor. 1). It shows that Moon's 1972 upper-bound construction is optimal, so
  `f(m) = min{3^a·5^b : 2a + 3b = m − 1}`. At m = 7 that is f(7) = 5² = 25 > 21, and monotonicity carries it to every
  m ≥ 7.
- **The two rungs.** In Busch's formula f(5) = 3² and f(7) = 5², the squares of the two atoms 3 = H(C₃) and
  5 = H(strong T₄). The direct Camion–Moon count reaches the first square (B*(5) = 9) and not the second
  (B*(7) = 19 < 25). The gaps sit just below the squares: 7 < 3², and 21 < 5².
- **What the count misses: disjoint pairs** (computed).
  - At m = 5, 6, 7 the floor is attained only with no two disjoint odd cycles (Ω complete), with α₁ = 4, 7 and 12.
  - At m = 7 the floor is a single isomorphism class, with scores (1,1,2,3,4,5,5).
  - The fewest odd cycles at m = 7 is 9, exactly Moon's bound 5 + 3 + 1. It is attained only by Moon's T′₇ (unique up to
    isomorphism), and only together with three disjoint pairs of triangles, so H = 31.
  - So the count is exact in α₁. What it misses at m = 7 is α₂: the true constraint is α₁ + 2α₂ ≥ 12. Few odd cycles
    force disjoint pairs; no disjoint pairs force many odd cycles.

## 3. The mod-4 locks (§C)

By THM-466, `H ≡ 1 + 2α₁ (mod 4)`, so H mod 4 is the parity of the number of odd cycles.
- **n = 5, around 7:** the tournaments with 5 ≤ H ≤ 9 are exactly (c₃, c₅, H) = (2, 0, 5) and (3, 1, 9), 240 each. Both
  are 1 (mod 4), while 7 is 3 (mod 4). Three triangles force a fourth odd cycle, a 5-cycle (Camion), which flips the
  parity of α₁.
- **Strong n = 6, around 21:**

  | c₃ | c₅ | d₃₃ | H |
  |---|---|---|---|
  | 4 | 2, 3, 4 | 0, 1 | 15, 17 |
  | **5** | **4, 6** | **0, 1** | **19, 23, 27** |
  | 6 | 4–8 | 0–2 | 23, 25, 29, 31, 33, 37 |
  | 7 | 7, 9 | 0–2 | 33, 37 |
  | 8 | 6, 8, 11, 12 | 1, 2, 4 | 41, 43, 45 |

  Only c₃ = 5 reaches the window around 21. There c₅ is always even and d₃₃ ≤ 1, so H = 11 + 2c₅ + 4d₃₃ is 3 (mod 4),
  while 21 is 1 (mod 4). **The gap 21 is born as a parity lock: five cyclic triangles force an even number of 5-cycles.**
  (For c₃ = 7, c₅ is always odd. For even c₃ there is no lock, consistent with THM-466(iv): H mod 4 is not a function
  of the scores.)
- So the two gaps take the two odd residues mod 4, 7 ≡ 3 and 21 ≡ 1. Each sits in a window locked to the other residue.
  In the owner's language this is the "subtle even/odd" differentiation: both are odd (Rédei), and they differ in the
  next binary digit.
- The lock is a finite computed fact: 3360 labelled tournaments in six isomorphism classes (three converse pairs). A
  hand proof is listed in §9.

## 4. Camion's theorem as an axiom: what it buys

Take as axioms Rédei/OCF (H = I(Ω, 2)), Camion (strong ⟹ Hamiltonian cycle), Moon (cycle counts and the
subtournament theorem), Busch (the floor) and multiplicativity over strong components. They generate:
- **Rédei and the 2-adic digit tower** (THM-466): H is odd; H mod 2^k reads the counts of fewer than k disjoint odd cycles.
- **The gap 7** (Camion + Moon + multiplicativity), and **its permanence** (the monotone floor).
- **The gap 21 and its permanence**: Busch for m ≥ 7; the exhaustive strong spectrum for m ≤ 6, whose delicate part is
  the computed n = 6 lock; and the factorisation 21 = 3·7.
- **Decidability: completeness is a Π₁ sentence.**

  **Lemma (PROVED).** For odd m, m is an H-value iff some tournament on at most m vertices has H = m.

  *Proof.* Take T with H(T) = m. Deleting the trivial (one-vertex) strong components changes no other strong component,
  so it leaves H unchanged. Each remaining component C has |C| ≤ H(C), by Camion (a Hamiltonian cycle gives |C|
  distinct Hamiltonian paths), and H(C) ≥ 3. A sum of integers that are each at least 2 is at most their product, so
  Σ|C| ≤ ΣH(C) ≤ ΠH(C) = m. ∎

  With Busch's floor, |C| ≤ 1 + 3 log₅ H(C), so at most log₃ m + 3 log₅ m = O(log m) vertices are needed.
- **The reduction of completeness to two prime lanes** (THM-4094), and the solid-interval constructions that have
  realized every allowed value up to 80,405 (THM-4097, THM-4102 and THM-4104).

The link between the two gaps: each lies just below a square rung of Busch's floor, in the crevasse that rung opens
(7 < 3² = f(5), 21 < 5² = f(7)). Camion's theorem alone gives the linear floor H ≥ m. With Moon's counts, the direct
count reaches f(5) but not f(7) (B*(7) = 19). So the gap 7 is within reach of the direct count, and the gap 21 is the
first gap beyond it: the first one where disjoint odd cycles matter.

## 5. Collatz versus "{7, 21} are the only gaps": the attempt and the verdict

**What an isomorphism would have to do.** It would have to send a counterexample of one problem to a counterexample of
the other: a divergent or nontrivially cyclic Collatz orbit to an odd m ∉ {7, 21} that is not an H-value, or the
reverse.
- **Logical form.** By the Lemma of §4, the {7,21} conjecture is Π₁: a counterexample is certified by a finite search.
  Collatz ("∀n ∃k: T^k(n) = 1") is Π₂, and its natural generalization is Π₂-complete (S. A. Kurtz, J. Simon, *The
  undecidability of the generalized Collatz problem*, TAMC 2007, LNCS 4484). A Π₂-complete problem cannot be reduced to
  a Π₁ one, not even by a Turing reduction. So **no uniform reduction (one that works for the generalized problems) can
  exist**.
- **The specific sentences.** Two true sentences are trivially equivalent, so an equivalence proof is not ruled out
  logically. But a proof of the equivalence together with a proof of completeness would prove Collatz. Once
  completeness is proved (the solid-interval route of THM-4097/4102/4104 is a route of exactly that kind), the
  "isomorphism" is as hard as Collatz. If completeness were false, the equivalence would refute Collatz.
- **The correct structural match** is the **cycle half** of Collatz: "the only cycle is 1 → 4 → 2 → 1". It is Π₁, and
  like the {7,21} conjecture it asserts that a known finite exceptional set is the whole exceptional set. The divergence
  half has no tournament counterpart, because Camion's theorem bounds the size of every witness.
- **Objects examined in the attempt** (§G, §I):
  - *Tournaments with T − v transitive.* Here H = 1 + 2K(w), where w is the binary word recording which vertices v
    beats, and K(w) = Σ over "10-inversions" b < a of 2^max(0, a − b − 2).
    - *Proof.* T − v is transitive, so every cycle passes through v and Ω is complete. A cycle is
      v → b → (an increasing path) → a → v, with w_b = 1, w_a = 0 and b < a. It is odd iff it uses an even number of
      the a − b − 1 vertices strictly between b and a, which gives 2^max(0, a − b − 2) odd cycles per pair. (Checked
      for every word of length ≤ 7.)
    - This is the most Collatz-like (2-adic) H-formula we found.
    - The family misses K = 3 and K = 10 (H = 7, 21) as it must, but it also misses 111 of the first 256 values of K.
      (Leading 0s and trailing 1s of w do not affect K, and a stripped word of length L has K ≥ 2^(L−3), so words of
      length ≤ 10 already give every K < 256.)
    - It is not a model of the spectrum, let alone of Syracuse dynamics.
  - *Numerology.*
    - {1, 2, 4} = QR(7) = the trivial cycle.
    - 7 = 1 + 2 + 4 and 21 = 1² + 2² + 4², i.e. q² + q + 1 at q = 2, 4: the point counts of PG(2,2) and PG(2,4).
      Equivalently (H − 1)/2 = 3 and 10 are the triangular numbers T₂ and T₄.
    - 7 ≡ 3 (mod 4) rises under Syracuse (7 → 11 → 17 → 13 → 5 → 1), while 21 ≡ 1 (mod 4) drops to 1 at once.
    - None of these maps Syracuse dynamics to Hamiltonian paths.
- **Verdict.** The isomorphism claim fails in every uniform sense. As a specific equivalence it is as hard as Collatz
  once completeness is proved. The true kinship is narrower and exact: the {7,21} conjecture has the logical shape of
  the Collatz *cycle* conjecture, and it has that shape *because of* Camion's theorem.

## 6. Tournaments as polyhedra (§D, §E, §F)

**6.1 Transitivity and the root system A_{n−1} (classical).**
- A tournament chooses one root from each pair ±(e_i − e_j) of A_{n−1}.
- Transitive tournaments are the Weyl chambers, i.e. the n! orderings.
- The score vector is the half-sum of the chosen roots plus a constant, and `Σ_T x^score(T) = Π_{i<j}(x_i + x_j) = s_δ`
  (the Schur function of the staircase δ). *Proof:* the bialternant is V(x²)/V(x), where V is the Vandermonde
  determinant; checked for n ≤ 5.
- So the number of tournaments with a given score vector is the weight multiplicity of `V_ρ`. Landau's theorem says the
  score vectors are the lattice points of the permutohedron conv(S_n·ρ).
- The transitive tournaments, joined by single arc reversals, form the permutohedron's graph. At n = 4 this is the
  truncated octahedron, which is also the flag graph of the tetrahedron (§E): **a transitive 4-tournament is a flag
  (vertex, edge, face) of the tetrahedron.**

**6.2 The n = 4 shells, the cuboctahedron and the rhombic dodecahedron.** Around the centre (3/2, 3/2, 3/2, 3/2), the 38
score vectors of 4-tournaments form three shells:

| squared distance | score vectors | polyhedron | tournaments per vector | kind |
|---|---|---|---|---|
| 1 | 6, type (1,1,2,2) | octahedron | 4 | strong |
| 3 | 8, types (3,1,1,1), (0,2,2,2) | cube | 2 | a 3-cycle plus a source or sink |
| 5 | 24, type (0,1,2,3) | truncated octahedron | 1 | transitive |

- The octahedron points are the face centres of the cube.
- Reversing an arc moves the score by a root e_j − e_i. The 12 roots form a **cuboctahedron**, the medial graph of the
  cube and of the octahedron: the two cyclic shells.
- The score lattice is A₃, the face-centred cubic lattice. Its Voronoi cell is the **rhombic dodecahedron**, whose 14
  vertices are 6 octahedral and 8 tetrahedral holes.
- The centre, the would-be "regular" score that no 4-tournament has, is an octahedral hole, surrounded by the 6 strong
  score vectors.
- At n = 5 the centre (2,2,2,2,2) is a lattice point of multiplicity 24: the 24 regular tournaments, which the noble
  note identified with the faces of the dodecahedron and the great stellated dodecahedron.

**6.3 Duality on the tetrahedron (PROVED; checked on all 64).** Give K₄ its embedding as the tetrahedron and orient the
dual edge of each arc from its left face to its right face. A face is a directed triangle iff the dual vertex is a
source or a sink (planar duality of cycles and cuts). So `c₃(T) = #sources + #sinks of T*`.
- This is an involution (up to converse): transitive (24) ↔ strong (24), and the middle class (16) maps to itself.
- The arcs are fixed. They are the vertices of the medial graph M(K₄), the octahedron. This is the pasted triality
  (V ↔ F, E fixed) acting on tournaments.

**6.4 P7 on the torus (PROVED; computed).** K₇ triangulates the torus with faces {i, i+1, i+3} and {i, i+2, i+3}
(mod 7): the A- and B-faces, the lines of the two Fano planes {0,1,3} + i and {0,2,3} + i.
- **The 14 faces are exactly the 14 cyclic triangles of P7** (i → j iff j − i ∈ {1, 2, 4}), and exactly two
  orientations of the 21 edges make every face a directed triangle: P7 and its converse.
- Every arc runs along its A-face and against its B-face. So in the dual every A-face is a source and every B-face a
  sink. The dual is the Heawood graph (cubic, bipartite, girth 6), oriented points → lines.
- Under Aut(P7) = Z₇ ⋊ Z₃ (order 21) the map has:
  - one orbit on vertices (7);
  - one regular orbit on edges (21);
  - two orbits on faces (7 + 7), alternating around every vertex.

  **This is the orbit pattern of a quasiregular polyhedron** (cuboctahedron: V, E transitive, faces 8 + 6). The dual
  has the pattern of a **quasiregular dual** (rhombic dodecahedron: F, E transitive, vertices 8 + 6; Heawood: vertices
  7 + 7).
- These are patterns of the subgroup Aut(P7). The map itself is the chiral regular map {3,6}_(2,1). Its group
  AGL(1,7), of order 42, is face-transitive; for example x ↦ 3x swaps A ↔ B. The octahedron shows the same pattern under
  the tetrahedral rotation group.
- The medial graph (21 vertices, 4-regular) carries a regular action of Aut(P7), so it is a Cayley graph of the
  Frobenius group of order 21. It is the fixed point of the duality here, as the cuboctahedron is for cube and
  octahedron.
- H(P7) = 189 = 9·21. P7 has 24 Hamiltonian cycles: one free Aut-orbit of 21, plus the 3 circulant cycles x → x + d,
  d ∈ {1, 2, 4}.
- (V, E, F) = (7, 21, 14) = 7·(1, 3, 2), with V + F = E. The two forbidden H-values are V and E of this map.
  **NUMEROLOGY**: no mechanism connects the torus map to H = 7 or H = 21. One difference is real: doubly regular
  tournaments (the Paley type) need n ≡ 3 (mod 4). P7 exists, and there is no "P21", since 21 ≡ 1 (mod 4).

**6.5 Equations in three variables (classical; the reading is ANALOGY).**
- Euler's V − E + F = χ (2 on the sphere, 0 on the torus).
- The triangle groups (p, q, r), classified by the sign of 1/p + 1/q + 1/r − 1:
  - (2,3,5) is spherical, the icosahedral group of order 2/(31/30 − 1) = 60, home of the noble polyhedra;
  - (2,3,6) is flat, and the K₇ torus map is a {3,6} map;
  - (2,3,7) is hyperbolic: the Klein quartic, whose group PSL(2,7) contains Aut(P7) and is the collineation group of
    the Fano plane.
- Hill's orbit parameters (a, b, c) are distances to the three mirrors of a (2,3,5) or (2,3,4) triangle.

**6.6 Beyond P7: Paley face-cyclic triangulations (COMPUTED; independent check pending).**
- Every face 2-colourable orientable triangulation of K_n induces a tournament in which every face is a directed
  triangle: orient each edge along its colour-1 face. Such a triangulation is a biembedding of two Steiner triple
  systems.
- Search, for primes p ≡ 1 (mod 6) and p ≡ 3 (mod 4):
  - take Z_p-invariant Steiner triple systems whose blocks are all directed triangles of the Paley tournament P_p;
  - take pairs of such systems with no common block, orienting one along P_p and the other against it;
  - keep those whose link at a vertex is a single cycle (a genuine surface).

| p | systems | disjoint pairs | genuine surfaces (genus) | up to multipliers |
|---|---|---|---|---|
| 7 | 2 | 1 | 1 (torus) | 1 |
| 19 | 8 | 4 | 0 (every link splits into 3 cycles) | 0 |
| 31 | 192 | 9,056 | 1,280 (genus 63) | 112 |
| 43 | 1,024 | 93,696 | 0 (links split into 3, 5 or 7 cycles) | 0 |

- So P31 plays P7's role on a genus-63 surface, while P19 and P43 have no Z_p-invariant surface of this kind. With two
  data points on each side, the split p ≡ 7 versus p ≡ 19 (mod 24) is a pattern, not a conjecture yet.
- On P31's surfaces the faces are 310 of P31's 1,240 cyclic triangles. P7 is the only case where the faces are all the
  cyclic triangles.

## 7. The triplet principle: exact facts, analogical reading

TCPC (§4) made the principle exact as an involution with one fixed point plus a factor-2 cocycle on the swapped pair.
It is exact for Z/3 with doubling (2 = −1) and for Pythagorean legs. New instances (exact facts; the triplet reading is
ANALOGY):

| triplet | primarily differentiated | the pair and its halving/doubling | status |
|---|---|---|---|
| (V, E, F) of a polyhedron | E (fixed by duality; the medial graph's vertices) | V ↔ F swapped by duality | exact facts (classical); triplet reading ANALOGY |
| 4-tournaments under tetrahedral duality | the 16 "middle" ones (fixed class) | transitive ↔ strong (24 ↔ 24) | exact facts (§6.3); triplet reading ANALOGY |
| n = 4 score shells | transitive (24 vectors, multiplicity 1) | cube ↔ octahedron, a dual pair, multiplicities 2 and 4 | exact facts; triplet reading ANALOGY |
| the gaps {7, 21} | both odd (Rédei) | 7 ≡ 3, 21 ≡ 1 (mod 4): the next binary digit (α₁ mod 2) | exact facts (§3); triplet reading ANALOGY |
| 60 = 2²·3·5 = \|A₅\| | 2: the edge orbit 30 = 60/2, fixed by duality | 3 ↔ 5: the orbits 20 ↔ 12 (dodecahedron ↔ icosahedron vertices), swapped by duality; from points to axes the stabilizers double (2, 3, 5 → 4, 6, 10) and the orbits halve (30, 20, 12 → 15, 10, 6) | exact facts (§H); triplet reading ANALOGY |
| P7 / K₇ on the torus | E = 21 = \|Aut(P7)\| (regular orbit) | V = 7 and F = 14 = 2V (two chiral face classes) | exact numbers; triplet reading ANALOGY |

None of these rows supplies a factor-2 cocycle telling the swapped pair apart. For 60, TCPC §4.3 shows the divisor
lattice has none between 3 and 5.

**7.1 The 3-4-3 split of `p²qr` read at 60.** The ten proper divisors of 60 split 3-4-3 as primes {2,3,5},
squarefree composites {6,10,15,30} and square-containing {4,12,20} (TCPC §4.3). At N = 60 they are **exactly** the
stabilizer orders {2,3,5,4,6,10} and orbit sizes {30,20,12,15,10,6} of the icosahedral rotation group.
- The five complementary pairs (d, 60/d) are the three point orbits (2, 30), (3, 20), (5, 12) (edges, dodecahedron
  vertices, icosahedron vertices), and the axis orbits (4, 15) and (6, 10) = (10, 6). The last pair serves both the
  3-fold axes (D₃, 10 axes) and the 5-fold axes (D₅, 6 axes), and duality swaps them inside it.
- Reading the square 2² as the doubling C_k → D_k when oriented points become unoriented axes (the antipodal quotient)
  is ANALOGY. This is the same quotient that folds the dodecahedron onto the Petersen graph in the noble note.
- For general `p²qr` this has no polyhedral meaning: it is exact at 60 only. "`p²qr` is `p` and `p³` combined" stays
  ANALOGY, as in TCPC.

## 8. Verdicts

- The theorem that plays Camion's role for 21 is **Busch's theorem at m = 7** (f(7) = 25), plus the **n = 6 parity
  lock** for its birth. The direct Camion + Moon count certifies the gap 7 exactly (B*(5) = f(5)) and falls short of
  the gap 21 at m = 7 only (B*(7) = 19). What it misses there is disjoint odd cycles. PROVED + COMPUTED.
- The 21-problem contains the 7-problem: (1 + x)(1 + 3x). PROVED.
- "Collatz is isomorphic to {7,21}-completeness": **no** uniform isomorphism (Π₁ versus Π₂-complete). A specific
  equivalence is as hard as Collatz once completeness is proved. The real kinship is with the Collatz cycle half,
  through Camion's theorem. PROVED (logic) + CITED (Kurtz–Simon).
- Tournaments → polyhedra: the root-system and duality dictionary of §6 is exact.
  - The rhombic dodecahedron is the Voronoi cell of the n = 4 score lattice, and the cuboctahedron its arc-reversal
    roots.
  - P7 is the face-cyclic orientation of the minimal torus, with the quasiregular / quasiregular-dual orbit patterns
    under Aut(P7). P31 does the same for 1,280 genus-63 surfaces.
  - The links from these objects to the values 7 and 21 are NUMEROLOGY.

## 9. Next steps

- **A hand proof of the n = 6 lock**: a strong 6-tournament with five cyclic triangles has an even number of 5-cycles.
  This would be the elementary "Camion lemma" for the birth of 21.
- **A short OCF proof of f(7) ≥ 23** (equivalently α₁ + 2α₂ ≥ 11 for strong 7-tournaments), trading odd cycles against
  disjoint pairs as T′₇ and the floor tournament do. Busch's proof covers all m; an OCF proof at m = 7 would show
  exactly where the disjoint pairs enter.
- **Face-cyclic orientations of neighbourly triangulations** (§6.6). K_n triangulates an orientable surface iff
  n ≡ 0, 3, 4, 7 (mod 12) (Ringel–Youngs). Which primes p admit a Z_p-invariant Paley face-cyclic triangulation? The
  data so far: yes for p = 7, 31 and no for p = 19, 43. Compare with the literature on cyclic biembeddings of Steiner
  triple systems.
- **Completeness** is to be pursued by the solid-interval constructions (THM-4097, THM-4102, THM-4104), not through
  Collatz.

## 10. Reproduction

- `python 04-computation/experiments/camion_busch_gaps_polyhedra_collatz_20261002.py` takes about 25 seconds and runs
  sections A–I.
- `gcc -O2 camion_busch_strong7_20261002.c && ./a.out` takes about 40 seconds. It is the exhaustive strong 7-tournament
  census:
  - f(7) = 25, with the floor a single isomorphism class;
  - minimum α₁ = 9, attained only by Moon's T′₇, with α₂ = 3;
  - the OCF checked on all 1,677,488;
  - H ≠ 21 at n = 7.
- `python 04-computation/experiments/camion_busch_paley_triangulations_20261002.py <p>` reproduces the §6.6 table.
- Not computed here: f(8) = 45 and f(9) = 75 (canon THM-1370; Busch 2006), and the Kurtz–Simon theorem (cited).

## 11. Audit record

An independent audit on 2026-10-02 used its own code (an exhaustive census of all 2,097,152 7-tournaments, graph
enumeration, and polyhedra checks). It read the primary sources (Moon 1966, Busch 2006, Kurtz–Simon 2007) and
reproduced both committed outputs exactly.
- **Verdicts.**
  - SOUND: §1, §3, §9 (reproduction).
  - SOUND with minor corrections: §4, §5, §6.
  - SOUND WITH CORRECTION: §2, §7, typing.
  - UNSOUND: nothing.
- **Main correction (§2).** The draft used THM-115's pancyclicity count (c_k ≥ ⌈m/k⌉) and called it "Moon's bound".
  Moon 1966 Cor. 2.1 proves the sharp c_k ≥ m − k + 1. The corrected direct count B*(m) = 1 + 2⌊(m−1)²/4⌋ already
  exceeds 21 at m = 8. So:
  - "Moon's bound of 8" became 9 (attained by T′₇);
  - "B(8) = 21; first exceeds 21 at m = 9" became the m = 7-only shortfall;
  - "proving Moon's 1972 conjecture" became Busch's actual statement (Moon's construction is optimal).
  The verdict that the direct count does not certify 21 survives, now located exactly at m = 7.
- **Other corrections.**
  - The Π₁ lemma bound is m vertices (sum ≤ product), and O(log m) with Busch, not O(log² m).
  - The "as hard as Collatz" argument is now stated relative to a proof of completeness.
  - The K(w) formula now has its proof.
  - §6.4 now says that the quasiregular patterns are those of the subgroup Aut(P7) (the full map group AGL(1,7) is
    face-transitive).
  - The §7 "exact instances" are retyped: exact facts, but the triplet reading is ANALOGY (no factor-2 cocycle).
  - Canon citations: the case lists are in THM-075, THM-079, THM-115 and THM-200. The completeness theorems are THM-4097,
    THM-4102 and THM-4104 (THM-4098–4101 and THM-4103 are unrelated).
  - The script's "Camion + Moon cannot certify it" print is corrected.
- **Canon aside from the audit.** THM-079 Part B's P₄ argument, as written, does not treat two triangles sharing an
  edge. The all-n result stands through THM-1370.
- §6.6 was added after the audit and is under a separate independent check.
