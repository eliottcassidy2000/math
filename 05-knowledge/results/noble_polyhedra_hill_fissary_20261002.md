# Noble polyhedra (Hill, arXiv:2607.28711) through this repo: four fissary objects, one golden quartic field, genus 21, Petersen folding, and an A₅ dictionary for 5-tournaments

**opus S15, 2026-10-02.**

- Script: `04-computation/experiments/noble_polyhedra_hill_20261002.py` (+ `.out`, ALL CHECKS PASSED).
- Data, fetched at run time and not vendored: Hill's model library (GPL-3.0, commit `a801da7`) and the paper's
  appendix tables (parsed from the arXiv HTML).
- Related: [`petersen_product_4polytope_attempt_20261001.md`](petersen_product_4polytope_attempt_20261001.md) (THM-4535,
  Petersen and the hemi-dodecahedron); the tournament H-spectrum canon (THM-338, THM-343, THM-1370, THM-4094);
  [`collatz_golden_holonomy_20261001.md`](collatz_golden_holonomy_20261001.md) (base φ).

**Status.**
- CITED: Hill's theorem. Besides the stephanoids and disphenoids there are exactly 146 noble polyhedra. The proof is
  computer-assisted. Checked here against Hill's model files:
  - the library holds 146 nobles plus two fissary models;
  - exactly four nobles have coplanar faces, as the paper says;
  - seven printed table entries contradict the model files (§2): four rows of Tables 5–9 and three decimal orbit
    locations in Tables 12–13. The paper's own data refutes each one (a printed dual row, or the printed minimal
    polynomial). Neither the count nor any minimal polynomial is affected.
- **Correction to the prompt.** The paper has **four** fissary objects, not two of the 146. They are the duals of four of
  the 146 (gD-19.1, gD-28.1, D-4, D-5) and are not themselves counted. Two of them, tI-F and rD-F, are in Hill's model
  library; they are the two whose faces are pentagons in either reading.
- COMPUTED (exact, from the model files and the paper's tables):
  - **The orientable genus spectrum of the 146 is {0, 4, 6, 7, 11, 13, 16, 21, 31, 61}.** Genus 21 occurs exactly once,
    for D-4, the parent of the fissary D-F1. Genus 7 occurs five times.
  - **All four fissary-related orbits lie in one quartic field**, `Q(φ, √(4φ−3))`, of discriminant −5²·19. In it the tI-F
    location `x` satisfies `x(x+1) = 1/φ`, and rD-F is `x + 1`. The field is not fissary-specific: of the 100 orbits in
    Hill's location tables, exactly one more lies in it, sD-9.
  - The dodecahedral nobles fold under the antipodal map onto the **Petersen graph** or onto `K₁₀` = Petersen ∪ `J(5,2)`.
- PROVED (elementary; checked exhaustively):
  - **The cyclic triangles of the 24 regular tournaments on 5 points are exactly the faces of the dodecahedron and of the
    great stellated dodecahedron**, read in A₅. The two A₅-orbits of regular tournaments are these two noble polyhedra,
    which the Galois conjugation √5 ↦ −√5 swaps.
  - **The gap 7 at n = 5 is Camion's theorem plus Moon's bound.** A 5-tournament is strong ⟺ it has at least 3 cyclic
    triangles ⟺ it has a Hamiltonian cycle. Strong ones have H ≥ 9; the others have H ≤ 5.
- NUMEROLOGY: the two known permanent gaps of the Hamiltonian-path spectrum first appear at n = 5 (gap 7) and n = 6
  (gap 21). These happen to be the degrees of the two smallest nontrivial transitive actions of A₅: the natural one, and
  PSL(2,5) on the 6 axes of the icosahedron. First appearance is forced by max H(4) = 5 and max H(5) = 15, not by A₅.
- **The {7, 21} link proposed in the prompt: NUMEROLOGY.**
  - The fissary arithmetic involves 5 and 19, never 7.
  - Genus 21 (D-4) and genus 7 do occur, but no mechanism connects surface genus to Hamiltonian path counts.
  - What is real is an exact but elementary A₅ dictionary at n = 5: regular 5-tournaments correspond to the faces of D-1
    and D-6. It does not involve the values 7 or 21.
- Independent audit: DONE (2026-10-02). No item was unsound, and its corrections are applied (§11; MISTAKE-559). Found
  after the audit, computed by the script but not re-audited: the three location errors of Tables 12–13, the Camion–Moon
  reading of the gap 7, and the isomorphism-class reading of the `c₃ = 3` tournaments (§6).

## 0. The owner's prompt

"Consider deeply the attached groundbreaking paper and how it deeply intertwines with the themes of this repo. Connor
Hill was able to prove the complete classification of noble polyhedra and it involves two infinite families plus 146
individual shapes, but of those 146, 2 are fissary, and thus those should be most investigated regarding connections to
our repo's work, perhaps with {7,21} being the tournament hamiltonian path forbidden values."

## 1. What Hill proves, and the shape of the hard direction

A **noble polyhedron** is vertex-transitive and face-transitive; faces are planar and may self-intersect. The precise
conditions are Hill's Definition 2.3:
- the realization is faithful;
- adjacent faces are never coplanar.

**Theorem (Hill).** Besides the stephanoids (crown polyhedra) and disphenoids, there are exactly 146 noble polyhedra up
to similarity.
- Nonprismatic counts by orbit type: T, O, C, tO, tC, rC: 1 each; I: 4; ID: 6; D: 7; tI: 17; tD: 6; rD: 19; sC: 7;
  gC: 3; sD: 33; gD: 38.
- 33 of them were found in 2020 (Mikloweit and others). Two (sD-10.1 and sD-12.1) were found by Ben Klein while the paper
  was being written.

**The hard direction is completeness.** The proof has four steps:
1. **Orbits.** Vertex sets of nobles are orbits of point groups. Every point group is normal in a reflection group, so
   the orbits fall into **orbit types** with 0, 1 or 2 real parameters `(a, b)`: 23 nonprismatic types and 5 prismatic
   classes. The parameters are the distances of a generating point to the three mirrors of a fundamental triangle.
2. **Criticality.** Four vertices are coplanar iff a 4×4 volume determinant vanishes. That determinant is a polynomial
   of degree ≤ 3 in the parameters. Orbits with the same coplanarity pattern (`≡_c`) have the same noble facetings, so it
   suffices to facet one orbit per class, the typical class included. The search finds nobles in 1- and 2-parameter
   types only at **critical orbits**, where coplanarities hold that do not hold for the whole type. These lie on the
   zero sets of the volume polynomials, which are curves in the 2-parameter case.
3. **Finiteness.** Distinct irreducible cubic curves meet in finitely many points, so each type has finitely many
   classes.
4. **Search.** A computer facets one representative of each class, using sympy, numpy and the Wolfram Language.

Prismatic symmetry is handled by hand: faces have 3 or 4 sides, which forces disphenoids or stephanoids.

## 2. The fissary objects

A noble polyhedron can have two faces in one plane provided they do not share an edge. Its polar dual then has distinct
vertices that coincide. Hill calls such a dual **fissary**. Read with coinciding vertices it is a degenerate polyhedron;
read with a compound vertex figure it is not a polyhedron at all. Here is Hill's Table 4, re-derived here (§A, §B):

| fissary | parent | parent's coplanar faces | points | abstract vertices | edges | faces | Schläfli | orientable genus |
|---|---|---|---|---|---|---|---|---|
| tI-F | gD-19.1 | 60 planes × 2 faces | 60 | 120 | 300 | 120 | {5,5} | 31 |
| rD-F | gD-28.1 | 60 planes × 2 faces | 60 | 120 | 300 | 120 | {5,5} | 31 |
| D-F1 | D-4 | 20 planes × 3 faces | 20 | 60 | 120 | 20 | {12,4} | 21 |
| D-F2 | D-5 (chiral, symmetry 532) | 20 planes × 3 faces | 20 | 60 | 90 | 20 | {9,3} | 6 |

- The genus is that of the abstract surface, which is shared with the parent.
- The coplanarity check over all 148 model files finds exactly these four parents.
- tI-F and rD-F are the two in Hill's library. Their faces are pentagons in either reading.
- For D-F1 and D-F2 a face visits coinciding points, so with compound vertex figures there is no Schläfli type.

**Table inconsistencies (§H, §H2).** Seven printed entries disagree with the model files, and the paper's own data
refutes each.
- Tables 5–9 (Schläfli type and counts). The printed dual row agrees with the model each time:
  - D-2 is printed {9,3}; the model is {3,9} (dual row rD-5.1: {9,3}, 60/90/20).
  - D-7 is printed {5,3}; the model is {3,9} (dual row tI-5.7: {9,3}, 60/90/20).
  - tI-5.6 is printed {5,5} with 150 edges; the model is {6,6} with 180 edges (dual row rD-5.7: {6,6}, 180 edges).
  - rD-5.2 is printed with 180 edges and 60 faces; the model has 120 edges and 30 faces (dual row ID-5: {4,8},
    30/120/60).
  - The D-2, D-7 and rD-5.2 rows also fail `2E = pF = qV` on their own.
  - (Audit: tI-5.6 and tI-5.2 have the same 60 face planes and face vertex-sets, but different hexagons; tI-5.6 is
    chiral (532), tI-5.2 is not (*532). The same holds for rD-5.7 and rD-5.3.)
- Tables 12–13 (orbit locations). Each wrong decimal is copied from a neighbouring row. The model value is a root of the
  printed minimal polynomial each time:
  - sD-3 b is printed 1.27703280443828 (sD-1's b); the model gives 0.17390637642569.
  - sD-16 b is printed 0.86777225107780 (sD-17's a); the model gives 0.78083000327643.
  - gD-15 a is printed 0.39234579229990 (gD-14's a); the model gives 0.55463631656377.
  - The other 173 printed locations agree with the parameters recovered from the models (the distances to the mirrors
    of the Möbius triangle, normalised as in Hill's Table 1).

## 3. The number fields of the orbit locations (§D)

Hill's minimal polynomials for the orbit locations of tI-F, rD-F, gD-19 and gD-28 all have **discriminant −475 = −5²·19**.
Each factors over `Q(√5)` into two quadratics. The closed forms below, together with `√(5φ+1)·√(4φ−3) = φ + 4`, show that
all four lie in one quartic field, `K = Q(φ)(√(4φ−3)) = Q(φ)(√(5φ+1))`.
- **tI-F** sits at `x = (−1 + √(4φ−3))/2 = 0.43168…`, so `x(x+1) = 1/φ`.
- **rD-F** sits at `y = x + 1 = 1.43168…`, so `y(y−1) = 1/φ`. This is the golden equation `φ(φ−1) = 1` with `1` replaced by
  `1/φ`.
- **gD-28** sits at `a = b = (φ + √(5φ+1))/2 = 2.31651…`. **gD-19** sits at `a = (−φ + √(5φ+1))/2 = 0.69848…` with `b = φ`.
- Both radicands have norm −19 in `Z[φ]`.
- **No 7 and no 21.** In `Q(√5)` the primes 3 and 7 are inert, so no element of `Z[φ]` has norm ±3, ±7 or ±21.

**Census of the location fields (all 100 orbits of Tables 10–13).** Every quartic location in Hill's 1-parameter table
(Table 10) is quadratic over `Q(φ)`. The radicands, up to squares and units, are:

| radicand | norm | Table 10 | Tables 11–13 | field discriminant |
|---|---|---|---|---|
| φ | −1 | tI-2, tD-3, rD-4 | — | −400 |
| 2 | 4 | tI-1, rD-1 | sD-18 (b), sD-20, sD-21, sD-29 | 1600 (`Q(√2, √5)`) |
| 1+4φ or its conjugate | −11 | tI-4, tI-7, tD-1, rD-6, rD-7 | sD-10, sD-12 | −275 |
| 4φ−3 | −19 | tI-F, rD-F | sD-9 (a), gD-19 (a), gD-28 | −475 |

- tI-7's polynomial has discriminant −3²·275; its radicand is `5 + 9φ = φ²(1+4φ)`.
- The other Table 10 locations are rational (tC-1, rC-1, tD-2, rD-5), in `Q(φ)` (tI-3, tI-5, tI-6, tD-4, rD-2), in
  `Q(√2)` (rD-3, rD-8), or cubic (tO-1, discriminant −31).
- −475 and −275 are forced, since each is 25 times a squarefree number and the discriminant of a quartic field
  containing `Q(√5)` is divisible by 25. −400 follows from the tower formula (audit), and 1600 is the discriminant of the
  biquadratic field `Q(√2, √5)`.
- **The orbits all of whose location parameters lie in K are exactly tI-F, rD-F, gD-19, gD-28 and sD-9.** So K is not
  fissary-specific. sD-9 has `a` = the tI-F location and `b = 1/φ`. Its noble sD-9.1 is a self-dual chiral `{6,6}` with
  no coplanar faces, and its 60 points are not similar to the tI-F points.

## 4. Surfaces, automorphisms, and where 7 and 21 occur (§A, §C)

Each noble polyhedron is a map on a closed surface. Fissary points are split into one abstract vertex per link component.

**Genus spectrum of the 146, plus the 2 fissary models.** The 146 alone have 50 at genus 31.
- Orientable:

  | genus | 0 | 4 | 6 | 7 | 11 | 13 | 16 | 21 | 31 | 61 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | count | 7 | 2 | 5 | 5 | 1 | 4 | 34 | 1 | 52 | 1 |

- Non-orientable (crosscap number): 14 four times, 32 twenty-two times, 62 ten times.
- **Genus 21: D-4 only.** It is `{4,12}` with 20 vertices, 120 edges and 60 faces, and its dual is the fissary D-F1.
- **Genus 7:** tC-1.1, rC-1.1, sC-3.1, sC-5.1, sC-6.2. These are the orientable 24-vertex `{5,5}` nobles. The other four
  sC are non-orientable with crosscap number 14 = 2·7.
- No noble polyhedron has 7 or 21 points, vertices, edges or faces.

**Abstract map automorphisms (§C).** These are counted on flags. The geometric symmetry orders are Hill's (*532: 120,
532: 60), which the audit confirmed against the models.
- D-1, D-6, I-2 and I-3 are regular, with 120 automorphisms.
- **D-3 (`{6,6}`, genus 11) is an abstractly regular map with 240 automorphisms**, twice its geometric symmetry.
- **D-4 has 240 automorphisms in 2 flag orbits**, and **D-5 has 120 in 3 flag orbits**. Each is twice the geometric
  symmetry (120 and 60).
- The fissary pair tI-F and rD-F, and their parents, have exactly the 120 geometric automorphisms. Their Petrie
  (`r₀r₁r₂`) orbits all have length 50.

## 5. Petersen: folding the dodecahedral nobles (§E)

All seven nobles on the 20 dodecahedron vertices (the D orbit) are drawn from the dodecahedral distance classes 1, 3
and 4. The antipodal map turns the dodecahedron into the hemi-dodecahedron, whose graph is the Petersen graph `K(5,2)`.
Two antipodal classes are joined in `K(5,2)` when their lifts are at distances 1 and 4, and in the complementary Johnson
graph `J(5,2)` when they are at distances 2 and 3. Each quotient edge has four lifts, two at each of its distances.
`K(5,2)` has 15 edges and `J(5,2)` has 30.

| noble | edges at dodecahedral distance | centrally symmetric | quotient graph | Petersen-edge multiplicity | `J(5,2)`-edge multiplicity |
|---|---|---|---|---|---|
| D-1 (dodecahedron) | 1 | yes | **Petersen** | 2 | — |
| D-6 (great stellated dodecahedron) | 4 | yes | **Petersen** | 2 | — |
| D-3 (`{6,6}`, regular map, genus 11) | 1, 4 | yes | **Petersen** | 4 | — |
| D-2 (`{3,9}`) | 1, 3 | yes | `K₁₀` | 2 | 2 |
| D-7 (`{3,9}`) | 3, 4 | yes | `K₁₀` | 2 | 2 |
| **D-4** (fissary parent, genus 21) | 1, 3, 4 | yes | `K₁₀` | 4 | 2 |
| **D-5** (fissary parent, chiral) | 1, 4, and 30 of the 60 distance-3 pairs | no | `K₁₀` | 4 | 1 |

`J(5,2)` is the graph of the rectified 5-cell, a calibration graph in the 4-polytope pipeline. The two dodecahedral
fissary parents sit over the complete graph on the Petersen vertex set. D-4 covers every Petersen edge 4 times and every
`J(5,2)` edge twice. The chiral D-5 takes exactly one of the two distance-3 lifts of each `J(5,2)` edge. The central
inversion maps D-5 to its mirror image, which takes the other lift.

## 6. Tournaments and A₅

**Repo canon.**
- The number `H(T)` of Hamiltonian paths of a tournament is odd, and never 7 (THM-343; THM-338 at n = 5) or 21
  (THM-1370), for every n. Every other odd value up to 609 occurs by n = 8 (THM-1370). That 7 and 21 are the *only* gaps
  is the open spectrum-completeness conjecture (THM-1370, THM-4094).
- The mechanism of the two gaps: `H` is multiplicative over strong components, and the strong minima are 3, 5, 9, 15,
  25, 45, 75, … (Busch).
- The gap 7 first appears at n = 5, where `H ∈ {1,3,5,9,11,13,15}`. The gap 21 first appears at n = 6, where the missing
  odd values up to 45 are 7, 21, 35, 39 (§F2). 35 and 39 occur at n = 7 (explicit witnesses, §F3).

**n = 5 and the dodecahedron (§F).** The 20 vertices of the dodecahedron can be identified with the 20 3-cycles of A₅:
A₅ acts by conjugation as the rotations, and the antipode of a vertex is the inverse 3-cycle. A cyclic triangle
`a→b→c→a` of a tournament on {1,…,5} is the 3-cycle `(abc)`. So **the cyclic triangles of a 5-tournament form an
antipode-free set of dodecahedron vertices**, and its directed 5-cycles are elements of order 5 of A₅. A₅ carries two
invariant dodecahedral graphs on its 3-cycles: the edge graphs of D-1 and of D-6, swapped by the outer automorphism.
Which one is called D-1 is a convention, and every statement below is symmetric under the swap.

**Proposition (proved; checked on all 1024 labelled 5-tournaments).**
1. On 5 points `H = 1 + 2(c₃ + c₅)` (odd-cycle formula). The pairs `(c₃, c₅)` that occur are (0,0), (1,0), (2,0), (3,1),
   (4,1), (4,2), (4,3), (5,2), so `c₃ + c₅ = 3`, i.e. `H = 7`, never occurs. **This is Camion's theorem plus Moon's
   bound.** A 5-tournament is strong ⟺ `c₃ ≥ 3` ⟺ it has a Hamiltonian cycle.
   - A strong tournament has `c₃ ≥ n − 2 = 3` (Moon) and `c₅ ≥ 1` (Camion), so `H ≥ 9`.
   - A non-strong one has no Hamiltonian cycle, and its strong components have at most 4 vertices, so `c₃ ≤ 2` and
     `H ≤ 5`.
2. **The cyclic triangles of each of the 24 regular 5-tournaments form a "layer"**: an orbit of a 5-fold rotation on the
   20 vertices. The 24 layers are exactly the 12 faces of the dodecahedron D-1 and the 12 pentagram faces of the great
   stellated dodecahedron D-6.
3. **The two A₅-orbits of regular tournaments (12 + 12) are exactly D-1's faces and D-6's faces.** Odd permutations swap
   the orbits. In icosahedral geometry the outer automorphism of A₅ is the Galois conjugation √5 ↦ −√5. On exact
   coordinates it sends the edges of the dodecahedron to the edges of a great stellated dodecahedron (§F), and those of
   the icosahedron to the edges of a great icosahedron (§F2).

*Proof of 2–3.* A regular 5-tournament is a relabelled circulant `i → i+1, i+2`. Its cyclic triangles are the translates
of `{0,1,3}`, i.e. the orbit of one 3-cycle under `⟨ρ⟩`, where `ρ` is a Hamiltonian 5-cycle of the tournament. That is a
layer. The two layer types (pairwise dodecahedral distances 1, 2 versus 2, 4) are the two A₅-classes. The script checks
every case against the geometric D-1 and D-6 faces. ∎

**Three cyclic triangles (§F).** The 240 tournaments with `c₃ = 3` all have score sequence (1,1,2,3,3) and `c₅ = 1`.
- In 120 of them the three cyclic triangles lie in one layer (60 in a D-1 face, 60 in a D-6 face), and that layer is
  always an orbit of the forced 5-cycle ρ.
- In the other 120 the three triangles are pairwise at dodecahedral distance 2. They are the three neighbours of one
  vertex, which is a 3-cycle on a transitive triple of the tournament, and they lie in no layer.
- These two sets of 120 are exactly the two isomorphism classes of 5-tournaments with `c₃ = 3`.

So the layer picture tells the two classes apart, but it does not explain why three cyclic triangles force a 5-cycle.
Camion's theorem does.

**n = 6 and the icosahedron (§F2).** A₅ also acts transitively on 6 points: PSL(2,5) on `P¹(F₅)`, which is A₅ acting on
the 6 axes of the icosahedron.
- The 20 triples of axes split into the 10 face triples of the icosahedron I-1 and the 10 face triples of the great
  icosahedron I-4. These are the two hemi-icosahedra, and each is a 2-(6,3,2) design.
- Complementation swaps the two classes.
- Hence every pair of complementary triples is {an I-1 face triple, the complementary I-4 face triple}. This includes
  every vertex-disjoint pair of cyclic triangles counted by `d₃₃` in `H = 1 + 2(c₃ + c₅) + 4d₃₃` (checked on all 32768
  labelled 6-tournaments). It holds for every tournament and every labelling of the 6 vertices by the axes, so it
  carries no tournament information.

**Typing.**
- NUMEROLOGY: the gaps 7 and 21 first appear at n = 5 and n = 6, the degrees of the two smallest nontrivial transitive
  actions of A₅. First appearance is forced by max H(4) = 5 and max H(5) = 15.
- DICTIONARY (exact, elementary): at n = 5 the cyclic-triangle sets of the 24 regular tournaments (H = 15) are the faces
  of D-1 and D-6, one A₅-orbit each, swapped by the outer automorphism (Galois √5 ↦ −√5).
- TAUTOLOGY: at n = 6 the 20 triples split into I-1 and I-4 face triples, swapped by complementation, for any
  tournament.

None of these statements concerns the values 7 or 21. Their permanence (THM-343, THM-1370) comes from multiplicativity
and Busch's strong minima, which owe nothing to A₅.

## 7. Structural analogies: the shape of the hard direction

These are typed ANALOGY. They are what the paper offers this repo beyond the shared objects.
- **Existence only on a coincidence locus.** Nobles live only at critical orbits, where cubic volume polynomials vanish.
  Repo analogues:
  - Collatz cycles occur only at resonances `2^A ≈ 3^p` (cf. the S19 digest: negative cycles are 3-adic resonances).
  - Lonely-runner extremal configurations are rigid coincidences of fractional parts.
  - In Hill's case finiteness comes from counting intersections of finitely many cubic plane curves (Bézout). The
    analogy stops there. Finiteness of Collatz cycles is open, and Baker-type bounds on `|2^A − 3^p|` only give lower
    bounds on cycle length. Lonely-runner tight configurations have no such finiteness theorem either.
- **Infinite families plus finitely many sporadics.** The prismatic families (stephanoids, disphenoids) come from a
  degenerate structure, an `n`-fold axis that can be scaled, and the 146 sporadics come from the seven rigid finite
  groups. The `H`-spectrum is generated by products over strong components, and its two known permanent gaps live in
  the small rigid range `m ≤ 6`, below Busch's minimum 25 at `m = 7`.
- **Fissary = locally fine, globally unfaithful.** A fissary dual is abstractly a good polyhedron whose realization
  identifies vertices. Two repo objects have a similar shape:
  - the hemi-dodecahedron, a good abstract polyhedron with no faithful *noble* (symmetric) realization in R³. An orbit of
    10 vertices and one of 6 faces would need a point group of order divisible by 30, hence an icosahedral one, and
    these have no 10-point orbit. Non-symmetric faithful realizations do exist, for example from 6 generic planes
    (audit);
  - PPS's cellular `S¹×RP²` with graph `P×C₃`, which satisfies every local 4-polytope condition and fails only globally
    (THM-4535).

  Hill's nondegeneracy conditions (a faithful realization; adjacent faces not coplanar) play the role that (I) and the
  facet-intersection cuts play in the 4-polytope pipeline. They exclude realizations whose combinatorics is fine but
  whose geometry collapses. PPS's `S¹×RP²` complex passes those cuts and is excluded only by the homology cut (H), so the
  fissary ↔ `S¹×RP²` pairing goes through (H).
- **Method.** Both Hill's theorem and THM-4535 are exhaustive computer-assisted classifications of realizations under
  necessary conditions. Both are reproducible from a script, and neither has a proof certificate.

## 8. Verdict on the prompt's {7, 21} suggestion

- The fissary objects carry the primes 5 and 19 and the golden ratio, not 7.
- 21 does occur, as the genus of D-4, the unique genus-21 noble and the parent of the fissary D-F1. 7 occurs as the genus
  of five octahedral nobles.
- Nothing connects these genera to Hamiltonian path counts. **NUMEROLOGY.**
- No link between noble polyhedra and the gaps {7, 21} was found. The one exact connection to tournaments is the n = 5
  dictionary of §6: regular 5-tournaments, which have H = 15, correspond to the faces of D-1 and D-6. It does not
  involve 7 or 21. The match of n = 5 and 6 with A₅'s two smallest nontrivial actions is numerology. At n = 5 the gap 7
  is explained by Camion's theorem and Moon's bound, not by geometry.

## 9. Next steps

- **D-3, the regular `{6,6}` map of genus 11 with 240 automorphisms.** Identify it in Conder's census. Its antipodal
  quotient covers every Petersen edge 4 times. Is it the regular map of the "Petersen double cover" type?
- **The A₅ dictionary at n = 5.** Extend the `c₃ = 3` picture (a layer of ρ, or the neighbours of a transitive-triple
  vertex) to all twelve isomorphism classes of 5-tournaments.
- **Faithful versus fissary as a SAT condition.** Add Hill's coplanarity data as constraints to the 4-polytope pipeline,
  for star-polytope versions of the realization question.
- **Report the seven table inconsistencies (§2) to the author.**
- **Canon (flagged separately).** OPEN-Q-094's withdrawal note calls the spectrum "co-finite … with exactly 2 gaps" and
  says canon proves it. Busch's lower bound proves only that 7 and 21 are omitted; completeness is open per THM-1370 and
  THM-4094.

## 10. Reproduction

`python 04-computation/experiments/noble_polyhedra_hill_20261002.py` takes about 15 seconds once the data are cached. The
first run downloads Hill's model library (148 files) and the paper's HTML. Sections A–H2 are in the `.out`. Everything
in this note is computed by the script, except the following, which the independent audit confirmed with its own code:
- the validity and pairwise non-similarity of the 146 models, and their geometric symmetry orders;
- the flag-isomorphism of the split fissary models with the abstract duals of their parents;
- the field discriminant −400;
- the tI-5.6 / tI-5.2 remark in §2;
- the non-symmetric realization of the hemi-dodecahedron in §7.

## 11. Audit record

An independent audit on 2026-10-02 checked every claim with its own code and reproduced the script's `.out` exactly.
- Verdicts:
  - SOUND: the paper summary, surfaces and automorphisms, Petersen folding, the `(c₃, c₅, H)` table, and the
    regular-tournament proposition.
  - SOUND WITH CORRECTION: the number-field section, the `c₃ = 3` statistic, the n = 6 statement, the canon citation,
    and the typing in §6–§8.
  - UNSOUND: nothing.
- Corrections applied:
  - §3: the 1-parameter table does not use only the norms −11 and −19. √φ and √2 also occur, tI-7 had been omitted,
    and the −475 field also contains sD-9.
  - §6, canon: "takes every odd value except exactly 7 and 21" overstated canon. Completeness is the open conjecture of
    THM-1370 and THM-4094, and THM-1370 (21 is never attained, for all n) is now cited.
  - §6, `c₃ = 3`: the statistic misdescribed the 120 non-layer cases.
  - §6, n = 6: the disjoint-triangle statement was a tautology, and is now typed so.
  - §6–§8, typing: "the two smallest transitive actions of A₅" was typed DICTIONARY; it is now NUMEROLOGY, with
    "nontrivial". "Locates where the gaps first appear" and "the genuine link" are withdrawn.
  - §7: Collatz-cycle finiteness is open, and Baker bounds give only lower bounds on cycle length. The hemi-dodecahedron
    has non-symmetric faithful realizations. The `S¹×RP²` pairing runs through the homology cut (H).
  - Script: section G's canon wording is fixed. The odd-cycle formula is now checked on all 32768 6-tournaments.
    `abstract_structure` asserts two face-sides per segment. Every claim the audit found uncomputed is now computed:
    the `c₃ = 3` statistic, the field census, the cover multiplicities, the n = 7 witnesses, and the dual-row
    confirmation. The exceptions are listed in §10.
- Found after the audit while making the script compute every claim, and not re-audited:
  - the three location errors in Tables 12–13, checked against parameters recovered from all 100 models;
  - the identification of the layer / no-layer split of the `c₃ = 3` tournaments with their two isomorphism classes;
  - the Camion–Moon reading of the gap 7.
- MISTAKE-559 records the overclaims.
