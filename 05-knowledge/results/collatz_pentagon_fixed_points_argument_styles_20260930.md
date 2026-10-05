# Two answers to the pentagon operator's open problems and a re-pricing of Collatz argument styles: `P_13` is pentagon-expanding (an induced icosahedron with 108 tadpole hats inside `C_5^2(P_13)`), the fixed points of `C_k` on vertex-transitive graphs come from three mechanisms (face self-duality for `k = 4`, the Eisenstein rotation of vertex links for `k = 6`, winding cycles in powers of cycles for `k = 4, 5`; the Paley graph `P_13` is a fixed point of the hexagon operator and the honeycomb torus falls into it), the carry's weighted-mediant law generates every rational Collatz cycle, and the orbit series obeys a Hadamard identity

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, thirteenth note), 2026-09-30.
**Owner's directive:** "please spend a long session trying to find the
doubling gadget inside `C_5(P_13)` and to classify the `C_k`-fixed
vertex-transitive graphs. consider other possible argument styles for
collatz, and be unique in generating new possible proposals that could
eventually close the proof, look for surprising connections to other fields
of math that would allow the proof to become simpler than it seems."

**Status.** Part A: PROVED conditional on Theorem 3.5 of Gervacio–Maehara–Ramos
(arXiv:2604.18984, CITED; its edge-removal step verified here at level one)
with a FINITE-EXACT certificate; the monotonicity lemma PROVED. Part B:
Theorems B1–B3 PROVED (elementary), the census FINITE-EXACT (circulants
`n ≤ 24`, degree `≤ 6`; the Cayley-graph part of the census is recorded in
§2.6), the honeycomb and antiprism identifications FINITE-EXACT; the
[twelfth note](collatz_pentagon_operator_coherence_20260930.md)'s conjecture
D36 REFUTED (MISTAKES entry). Part C: the Hadamard identity, the
weighted-mediant theorem and the one-run inequality PROVED (elementary); the
censuses FINITE-EXACT; every proposal typed and priced; Collatz OPEN.
Independent audit OWED (weekly subagent limit); self-checks: two counters for
the icosahedron search (pole search and generic induced-subgraph search
agree), the mediant law on 500 random triples and on the `-17` word, the
Hadamard identity to order 40 for six seeds, isomorphisms by colour
refinement plus backtracking.
Scripts (`04-computation/experiments/`): `collatz_pentagon_gadget_20260930.py`
(+ `_selfdouble.py`), `collatz_ck_fixed_points_20260930.py` (+ `_maps.py`;
outputs `_part1.out`, `_part1_large.out`, `_census.out`, `_maps.out`),
`collatz_argument_styles_20260930.py` (+ `_p4.out`).

## 1. Part A: the doubling gadget inside `C_5(P_13)`

**Lemma 1 (monotonicity and additivity; PROVED).** If `H` is an induced
subgraph of `G`, then `C_5(H)` is an induced subgraph of `C_5(G)`; and
`C_5(G ⊔ G') = C_5(G) ⊔ C_5(G')`. *Proof.* Induced pentagons of `H` are induced
pentagons of `G`, and "sharing an edge" is intrinsic; a pentagon lies in one
component. ∎ Consequently, if some iterate `C_5^k(G)` contains an induced
copy of a pentagon-expanding graph, `G` is pentagon-expanding.

**What `C_5(P_13)` is.** `P_13 = C_13(1,3,4)` has 39 induced pentagons (all
non-contractible on its torus, §2), so `C_5(P_13)` has 39 vertices; it is
16-regular; `Aut(P_13) = Z_13 ⋊ Z_6` (order 78) is transitive on the 39
pentagons, and its index-2 subgroup `Z_13 ⋊ Z_3` acts freely, so `C_5(P_13)`
is a Cayley graph of the non-abelian group of order 39. It contains no
induced `P_13`, but it contains induced icosahedra (156 found by a pole
search; 50 by a generic search), none of which carries a tadpole hat.
`C_5^2(P_13)` has 4563 vertices and 758862 edges (degrees 266–387), in 62
orbits of `Aut(P_13)` (7 of size 39, 55 of size 78).

**The certificate (FINITE-EXACT).** For each induced icosahedron `I` of
`C_5(P_13)`, its twelve vertex-link pentagons are induced pentagons of
`C_5(P_13)`, hence vertices of `C_5^2(P_13)`, and by Lemma 1 they induce
`C_5(I) ≅ I` there. The 156 icosahedra lift to 13 distinct induced
icosahedra of `C_5^2(P_13)`, and **each of them carries 108 tadpole hats**
(vertices outside the icosahedron whose neighbourhood inside it is a triangle
with a pendant edge, `T_{3,1}`). So the hatted icosahedron `I_1` of the paper
is an induced subgraph of `C_5^2(P_13)`.

**Theorem A.** `P_13` is pentagon-expanding, given Theorem 3.5 of the paper
(an icosahedron with `h ≥ 1` tadpole hats is pentagon-expanding). *Proof.*
`I_1 ≤_ind C_5^2(P_13)`, so `C_5^k(I_1) ≤_ind C_5^{k+2}(P_13)` by Lemma 1, and
`|V(C_5^k(I_1))| → ∞`. ∎

**On Theorem 3.5.** Its proof passes through `I_2 := C_5(I_1)` minus the edge
`[p, q]` joining the two pentagons through the hat and asserts that
`C_5(I_2)` is an induced subgraph of `C_5^2(I_1)`; edge removal is not
covered by Lemma 1, so this was checked: `C_5(I_1)` has 14 vertices with the
two hat pentagons adjacent, `C_5(I_2)` has 16 vertices, it *is* an induced
subgraph of `C_5^2(I_1)` (17 vertices), and it contains an induced icosahedron
with four hats, as the paper says. The direct trajectory of `I_1` is
`13, 14, 17, 31, 408` vertices, with icosahedra carrying `2, 4, 10, 16` hats at
the four levels. A self-contained alternative (some `C_5^m(I_1)` containing
two vertex-disjoint, mutually non-adjacent copies of `I_1`, which would give
`|C_5^{mk}(I_1)| ≥ 13·2^k` by Lemma 1) fails through level 4: the copies
overlap (24, 48, 120, 4527 hatted copies at levels 1–4, no two disjoint and
non-adjacent). So the conclusion rests on the paper's induction, whose first
step is verified; it is typed accordingly.

**A reusable certificate, and `P_17`.** Any graph whose second pentagon
iterate contains an induced hatted icosahedron is expanding. The test does
not need `C_5^2` in full: for an induced icosahedron `I` of `C_5(G)`, the
candidate hats are the induced pentagons of `C_5(G)` through the edges of
the twelve link pentagons, and a hat is one that meets exactly four links
forming a tadpole in the lifted icosahedron (adjacency there being "share an
edge", the distance-2 relation of the twelfth note, not adjacency in `I`).
Validated on `P_13` (108 hats, as in the full computation:
`_p17_validation.out`), the test gives **`P_17` an induced icosahedron in
`C_5(P_17)` (272 vertices) whose lift carries 402 tadpole hats**
(`_p17.out`): `P_17` is pentagon-expanding under the same hypothesis as
Theorem A. For the Clebsch graph the pole search is unusable (a link in
`C_5(Clebsch)`, 192 vertices of degree 90, has more than `300000` induced
pentagons), but a generic induced-subgraph search finds an icosahedron whose
lift carries **583 hats** (`_clebsch.out`): the Clebsch graph is
pentagon-expanding too. So the three growing graphs of the twelfth note's
table (`P_13`, `P_17`, Clebsch) all expand, each by a hatted icosahedron two
levels down.

## 2. Part B: the fixed points of `C_k` on vertex-transitive graphs

**The counting constraint.** If `G` is vertex-transitive on `n` vertices and
`C_k(G) ≅ G`, then `G` has exactly `n` induced `k`-cycles and every vertex
lies in exactly `k` of them. The census below uses this as its filter.

### 2.1 Locally cyclic graphs: the induced `k`-cycles are the links (PROVED)

**Theorem B1.** Let `G` be a connected graph in which every vertex link is an
induced `k`-cycle with `k ≥ 5` (so `G` triangulates a closed surface with all
degrees `k`). Then (i) every contractible induced `k`-cycle of `G` is a vertex
link; (ii) distinct vertices have distinct links; (iii) if `G` has no
non-contractible induced `k`-cycle, then `C_k(G) ≅ R(G)`, the *rhombus graph*
on `V(G)` in which `x ~ y` iff `x, y` are the apexes of two triangles on a
common edge.

*Proof.* (i) A contractible induced `k`-cycle bounds a disk `D` triangulated
by `G`, with `V_int` interior vertices (all of degree `k`) and `k` boundary
vertices, `d(v) ≥ 1` triangles of `D` at each. Counting faces and corners,
`F = 2V_int + k - 2` and `3F = k V_int + Σ_∂ d(v)`, so
`Σ_∂ d(v) = (6 - k) V_int + 3k - 6`. Since the cycle is induced and `k ≥ 4`,
`d(v) ≥ 2` for every boundary vertex (a vertex with `d(v) = 1` makes its two
boundary neighbours adjacent, a chord). For `k > 6` this gives
`(k - 6) V_int ≤ k - 6`, so `V_int ≤ 1`; `V_int = 0` is impossible (a
triangulated polygon without interior vertices has a triangle with two
boundary edges, again a chord); so `V_int = 1` and all `d(v) = 2`. For
`k = 6`, `Σ d(v) = 12` with six terms `≥ 2` forces all `d(v) = 2`. In both cases
the triangle on each boundary edge has an interior apex (a boundary apex
would be a chord), consecutive boundary triangles share the edge from their
common boundary vertex to their apex, so all apexes coincide in one interior
vertex `u` of degree `k` adjacent to the whole cycle: the cycle is the link
of `u`. For `k = 5` the only locally-`C_5` graph is the icosahedron (CITED,
classical), whose induced pentagons are its twelve links (computed).
(ii) Two vertices with the same link `C` would both be apexes over every edge
of `C`, so the surface would be the double cone over `C`, whose equatorial
vertices have degree `4 ≠ k`. (iii) With links as the only induced
`k`-cycles and the link map injective, `C_k(G)` has vertex set `V(G)`, and
`link(x)`, `link(y)` share an edge `uv` iff `u, v ∈ N(x) ∩ N(y)` with `u ~ v`;
`x ~ y` would put a triangle in a link, impossible for `k ≥ 4`; so `x, y` are
the two apexes over `uv`. ∎

### 2.2 The Eisenstein tori are fixed points of the hexagon operator (PROVED)

**Theorem B2.** Let `I` be an ideal of `Z[ω]` (`ω^2 + ω + 1 = 0`) of index `n`,
and `G = Cay(Z[ω]/I, {±1, ±ω, ±ω^2})` the triangular torus, assumed locally
`C_6` with no non-contractible induced hexagon. Then `C_6(G) ≅ G` iff
`3 ∤ n`; the isomorphism is `x ↦ (1 - ω) x`.

*Proof.* By B1, `C_6(G) ≅ R(G)`; over the edge `{u, u+1}` the two apexes are
`u + 1 + ω` and `u - ω`, whose difference `1 + 2ω = (1 - ω) ω` has norm 3, so
`R(G) = Cay(Z[ω]/I, (1 - ω)·{units})`. Multiplication by `1 - ω` is a
well-defined bijection of `Z[ω]/I` iff `1 - ω` is invertible modulo `I`, i.e.
iff `3 ∤ n` (its norm is 3); then it carries the unit set to the rhombus set
and is a graph isomorphism. If `3 | n`, `(1 - ω)` maps onto a subgroup of
index 3 and `R(G)` is disconnected while `G` is connected. ∎

**Examples (FINITE-EXACT, `_part1.out`, `_part1_large.out`).** `P_13 =
C_13(1,3,4)` is locally `C_6`, its 13 induced hexagons are the vertex links,
`C_6(P_13) ≅ R(P_13) = Cay(Z_13, {non-squares}) ≅ P_13` (the multiplier
`1 - ω = 1 - 3 = -2` is a non-square, and Paley graphs are
self-complementary): **the Paley graph `P_13` is a fixed point of the hexagon
operator.** So are `C_37(1,10,11)`, `C_43(1,6,7)`, `C_49(1,18,19)`,
`C_61(1,13,14)`, `C_67(1,29,30)` (the circulants `C_n(1, a, a+1)` with
`a^2 + a + 1 ≡ 0 mod n`) and the square tori `Z_7 × Z_7`, `Z_8 × Z_8`; `C_19(1,7,8)`
and `C_31(1,5,6)` have non-contractible induced hexagons (`> 3n` hexagons)
and are not fixed; `Z_6 × Z_6` (`3 | 36`, 54 hexagons) and `Z_4 × Z_8` (not an
ideal quotient, 64 hexagons) are not fixed. `P_13` also has 39 induced
squares and 39 induced pentagons, all non-contractible (its torus has
systole 4), which do not disturb the hexagon count.

**The honeycomb falls into it.** The honeycomb torus over `Z[ω]/(13)` (26
vertices, cubic, 13 induced hexagons = its faces) has `C_6 ≅ P_13`, then fixed:
the `{6,3}` map is `C_6`-preperiodic into the `{3,6}` fixed point, exactly as
the dodecahedron is `C_5`-preperiodic into the icosahedron (`_maps.out`).

### 2.3 Square tori are fixed points of the quadrangle operator (PROVED)

**Theorem B3.** Let `G = Cay(Z^2/L, {±e_1, ±e_2})` be a square torus whose
induced 4-cycles are exactly its `n` faces `{x, x+e_1, x+e_1+e_2, x+e_2}`. Then
`C_4(G) = G` (literally, on the face labels `x`). *Proof.* Two faces share an
edge iff their labels differ by `±e_1` or `±e_2`: the dual of `{4,4}` is a
translate of itself. ∎ (No Eisenstein-type condition: the self-duality of
the square tiling needs no rotation.) In the census every `C_4`-fixed
circulant with `n ≥ 12` is such a torus (`C_n(a, b)` with faces the only
induced squares, e.g. `C_12(1,4)`, `C_13(1,5)`, `C_14(1,4)`, `C_15(1,4)`,
`C_16(1,6)`, `C_17(1,4)`, `C_20(1,4)`, `C_25(1,7)`), while `C_8(1,3) = K_{4,4}`,
`C_9(1,3)`, `C_10(1,3)`, `C_11(1,3)` have extra induced squares and are not
fixed. The octahedron, the only locally-`C_4` graph, is not fixed because
antipodal vertices share their link (`C_4 = 3K_1`).

### 2.4 The winding mechanism (FINITE-EXACT)

Three fixed circulants are neither face nor link constructions:
`C_7(1,2)` (`k = 4`; its seven induced squares are the rotations of the step
pattern `(1,2,2,2)`, winding once around `Z_7`; the lattice faces of
`C_7(1,2)` as a square torus have chords), `C_8(1,2)` (the 4-antiprism,
`k = 5`; eight induced pentagons with steps `(2,1,2,1,2)`; `C_5(C_8(1,2)) ≅
C_8(1,2)`, realised by the multiplier 3 up to relabelling) and
`C_14(1,2,3)` (`k = 5`, degree 6; fourteen induced pentagons with steps
`(3,3,3,3,2)`). In each, the fixed cycles are the shortest induced cycles that
wrap the ring once, and the pentagon graph is the same power of a cycle. The
twelfth note's conjecture D36 ("fixed iff locally `C_k` with the distance-2
condition") is therefore false; MISTAKES entry.

### 2.5 The answer to the paper's question, and the mechanism behind every fixed point

For each fixed point found, `C_k(G) ≅ G` is realised by an explicit
self-map: the identity (`K_4`, `k = 3`; the triangles form a `{3,3}` map, and
for any vertex-transitive `C_3`-fixed graph whose `n` triangles form a map
Euler's formula gives `χ = n/2`, so `n = 4`), the translation by
`(1/2, 1/2)` (square tori, `k = 4`), the antipodal map (icosahedron,
`k = 5`), multiplication by `1 - ω` (Eisenstein tori, `k = 6`), and a unit
multiplier of `Z_n` (the winding circulants). Periodic polyhedral fixed points
exist for `k = 3` (tetrahedron) and `k = 5` (icosahedron) on the sphere, for
`k = 4` and `k = 6` on the torus, and never for `k = 4` on the sphere. Dual
pairs of maps appear as operator steps (`C_3(oct) = cube`, `C_4(cube) = oct`,
`C_3(I) = D`, `C_5(D) = I`, `C_6(honeycomb) = P_13`).

### 2.6 The census (FINITE-EXACT)

Circulants `C_n(S)`, `n ≤ 24`, degree `≤ 6`, `k = 3..6`: `k = 3`: `K_4` and
its disjoint unions (`C_{4m}(m, 2m)`); `k = 4`: `C_7(1,2)` and the square tori
listed in §2.3 (every `n` from 12 to 24 has some, `n = 8..11` none);
`k = 5`: `C_8(1,2)`, `C_14(1,2,3)`, and `2 C_8(1,2)`, `3 C_8(1,2)`; `k = 6`:
`P_13` only. **Cayley graphs** (abelian groups of rank 2 with `ab ≤ 36`;
`A_4`, `S_4`, `SL(2,3)`, `D_3..D_8`, `Q_8`, `Z_7 ⋊ Z_3`, `Z_13 ⋊ Z_3`; degree
`≤ 6`; `_census_groups.out`, `_dihedral.out`): every `k = 3` fixed graph is
`K_4` or a disjoint union of `K_4`'s; the `k = 4` fixed graphs are square tori
on `Z_a × Z_b` (dozens of connection sets per group), the dihedral ones being
the circulants again (`Cay(D_6, S) ≅ C_12(1,4)`, `Cay(D_8, S) ≅ C_16(1,6)`),
and the `Z_7 ⋊ Z_3`, `Z_13 ⋊ Z_3` ones three disjoint copies of `C_7(1,2)`,
`C_13(1,5)`; the `k = 5` fixed graphs are the icosahedron (`A_4`, degree 5),
`C_8(1,2)` (as a Cayley graph of `D_4`) and its disjoint unions (`Z_2 × Z_8`,
`Z_4 × Z_8`, `Z_2 × Z_16`, `D_8`, `S_4`), `C_14(1,2,3)` (`D_7`, degree 6) and
its double (`Z_2 × Z_14`), two icosahedra (`S_4`, degree 5), and **one new
sporadic example: a connected 6-regular Cayley graph of `D_8` on 16 vertices
with 16 induced pentagons, six triangles at each vertex, and common-neighbour
counts `0, 2, 4` on edges, isomorphic to no 6-regular circulant on 16
vertices and neither to the Shrikhande graph (the `4 × 4` triangular torus,
96 induced pentagons, not fixed) nor to the rook graph**; the `k = 6` fixed
graphs are three copies of `P_13` in `Z_13 ⋊ Z_3` and nothing new. Caveat: the
isomorphism test is budgeted after colour refinement, spectrum and triangle
counts agree, and it gave up on 39 pairs (36-vertex abelian and 21/39-vertex
non-abelian candidates), which are reported as non-fixed; those could hide
further square tori, not new mechanisms.

## 3. Part C: argument styles for Collatz, re-priced, with three exact reformulations

### 3.1 The Hadamard identity (PROVED)

Let `F_n(z) = Σ_k 2^{d_k} z^k` be the valuation series of the Syracuse orbit
`m_0 = n, m_1, ...` of an odd `n` (`d_k` the halvings before the `k`-th odd
step) and `G_n(z) = Σ_k m_k z^k` the orbit series. Then, with `⊙` the
termwise product,
`(1 - 3z) (F_n ⊙ G_n)(z) = n + z F_n(z)`.
*Proof.* `2^{d_k} m_k = 3^k n + Σ_{i<k} 3^{k-1-i} 2^{d_i}` (the Terras identity)
summed against `z^k`. ∎ (Checked to order 40 for six seeds.) The orbit series
is the termwise quotient of a Möbius transform of `F_n` by `F_n`; Collatz
over `N` is "`G_n` is rational for every `n`", and `F_n` is rational iff `G_n`
is (Pólya; Prop 10 of the [first note](collatz_posets_dags_zeta5_20260927.md)).
The Hadamard-quotient theorem (Pourchet, van der Poorten, Rumely) says a
termwise quotient of two rational series that is integral is rational, so
nothing is gained from it: the identity organises the Pólya/Bézivin
statements, it does not close them.

### 3.2 The weighted-mediant law generates every rational cycle (PROVED)

With the cocycle of the twelfth note (`S_{uv} = 3^{p_v} S_u + 2^{A_u} S_v`,
`D_{uv} = 2^{A_u} D_v + 3^{p_v} D_u`, `x_w = S_w/D_w`):

**Theorem C1.** `x_{uv} = (α x_u + β x_v)/(α + β)` with `α = 3^{p_v} D_u`,
`β = 2^{A_u} D_v`, and the Farey determinant `S_u D_v - S_v D_u = -β(u, v)` is
the twelfth note's commutation defect. Hence every fixed point `x_w` is an
iterated weighted mediant of the letter points `x_{(v)} = 1/(2^v - 3)`
(`-1, 1, 1/5, 1/13, 1/29, ...`, the one-step rational cycles of the Collatz map
on `Z_(2)`), and the rational cycles of Lagarias (1990, CITED) are exactly
the values of this tree. *Proof.* Divide the cocycle identities. ∎

*Corollary (convexity).* If a word has no 1-step, all its sub-clocks are
positive (`2^A ≥ 4^p > 3^p`) and its fixed point lies in `(0, 1]`; so a
positive integer cycle other than `{1}` contains a 1-step (an element
`≡ 3 mod 4`), the classical fact, here for free.

*Corollary (the denominator law; D40 in the coprime case).* Always
`D_u S_{uv} ≡ 2^{A_u} β(u, v) (mod D_{uv})` (substitute `3^{p_v} D_u ≡
-2^{A_u} D_v` into the cocycle identity), so when `gcd(D_u, D_v) = 1` the
denominator of `x_{uv}` in lowest terms is `|D_{uv}| / gcd(β(u,v), D_{uv})`:
the mediant of two words with coprime clocks is integral iff the clock of
the whole divides the Farey determinant of the parts. (Checked on all
`16129` pairs of words with `A ≤ 7`, `13780` of them with coprime clocks;
`_klein_farey.out`.)

*Corollary (the one-run inequality).* If `w = (1)^k v` with `v` free of
1-steps and `x_w = N ≥ 2` an integer, then
`1 + 1/N < 2^k D_v / (3^{p_v} (3^k - 2^k)) ≤ 1 + 2/(N - 1)`:
Steiner's 1-cycle condition in mediant form (his theorem needs Baker's
method to finish; the inequality alone does not).

**Census (FINITE-EXACT, all `262143` words with `A ≤ 18`).** Denominators of
`x_w` (rational cycles) in lowest terms: `1` (46 words: the four integer
cycles and their repeats), `5` (24), `7` (8), `11` (26), `13` (74), `17`, `19`,
`23`, `25`, ...; e.g. `(1,1,3) ↦ 19/5`, `(1,3) ↦ 5/7`, `(1,1,2) ↦ -19/11`,
`(1,1,1,1,4) ↦ 211/13`. The closest positive near-misses of an integer
`N ≥ 2`: the word `(3,3,3,1,1,1,2,2,2)` has `x = 484921/242461 = 2 - 1/242461`
and `(4,1,1,4,1,1,3,2,1)` has `2 - 9/242461`, both of shape `(18, 9)` with
clock `2^18 - 3^9 = 242461`: **carries come within 1 of a multiple of the
clock**, so no metric ("the carry stays away from `0 mod D`") argument can
exclude cycles; only the exact residue decides. And the Farey determinant is
always even with minimum `|det| = 2` (at `(1), (2)`): the mediant tree has no
unimodular pairs, so the Stern–Brocot uniqueness argument (determinant one,
Euclidean descent) has no analogue.

### 3.3 The never-descending set and the sign (FINITE-EXACT)

Let `K ⊂ Z_2` be the set of 2-adic integers whose word never descends
(`2^{d_j} ≤ 3^j` for every prefix), the sibling ladder's `Bad(3)`, of
Hausdorff dimension `h(log_3 2) = 0.95`. `K ∩ [-10^5, -1] = {-1, -5, -17}` and
`K ∩ [1, 10^5] = ∅` (the last positive integer to leave `K` below `10^5` is
`35655`, at odd step 85). Collatz over `Z` is the statement
`K ∩ Z = {-1, -5, -17}`. The sign is invisible 2-adically (`N` and `-N` are
both dense in `Z_2`); the barrier atlas's SHEET control (any argument must
distinguish `3x+1` from `3x-1`) is, in these coordinates, the statement that
`K` and the mirror set of `3x - 1` differ only by which sign the three
integers carry.

### 3.4 Local against global (FINITE-EXACT; the obstruction typed)

Over the words of shape `(18, 9)`, `(18, 11)`, `(18, 12)` the carries `S_w`
modulo `ℓ = 5, 7, 11, 13` are uniform to within `14%`, `2%`, `2%`, `5%` (max
relative deviation), and the deviation decays with `A` (the residue map is a
transfer operator on `Z/ℓ` with a spectral gap — PLAUSIBLE, not proved here).
So there is no local obstruction to cycles at any fixed prime: the class of
the carry in `H^1(<w>; Z) = ⊕_{ℓ^e || D} Z/ℓ^e` (twelfth note) vanishes at a
given `ℓ` for about `1/ℓ` of the words. The cycle condition is global, and
counting cannot prove non-existence: a main term below one is not a proof.
This is why the cycle problem is Diophantine (Baker for few runs) and open
for many runs.

### 3.5 Proposals, priced

| proposal | closing statement it would need | why it could be simpler | price / verdict |
|---|---|---|---|
| **P1 mediant-tree descent** | an operation on words that lowers an integer height while preserving integrality of `x_w`, as the Euclidean algorithm does for Stern–Brocot | the tree structure is exact (C1) and the cycle condition is "a mediant lands on an integer" | blocked: no unimodular pairs (3.2); integrality does not descend to splits (`-17 = (1,1,1)(2,1,1,4)` has `x_v = 103/175`) |
| **P2 the sign of `K`** | `K ∩ N = ∅` for an explicit closed 2-adic set of dimension `0.95` | one set, one lattice `Z ⊂ R × Z_2`; a Diophantine-approximation shape (Cantor sets avoiding integers) | the set is sign-blind; the content is SHEET; no known tool couples the archimedean sign to the 2-adic word except the size price of the [sixth note](collatz_generic_price_topology_20260927.md) |
| **P3 rationality of the orbit series** | `G_n` rational for all `n` | Hadamard products/quotients of rational series are classical | circular: `(n + zF)/(1 - 3z) ⊙^{-1} F` is rational iff `F` is |
| **P4 cohomological transfer** | a transfer `H^1(<u,v>) → H^1(<uv>)` forcing the class of `uv` from those of `u, v` | classes are local at clock primes; the cocycle is free on letters | no finiteness (the monoid is free); the classes of splits do not determine the class of the whole (`β ≡ 0 mod D_{uv}` is all one gets) |
| **P5 expanding/vanishing certificates** | guarded residue classes with pointwise bounded descent, plus a proof that their rules cover all inputs | in the operator world certificates are induced subgraphs (Lemma 1) | **CORRECTED 2026-10-04:** such descent classes do exist: every `n=5 mod8` has `U(n)<n`. Descent on selected classes does not establish convergence without closing their dependencies. Separately, every complete arithmetic progression is sufficient in the common-future sense ([Monks et al., introduction and Theorem4.1](https://arxiv.org/pdf/1204.3904), CITED and primary checked): proving convergence of all its members would settle Collatz. These are different predicates. |
| **P6 fixed points as self-maps** | a self-map of the Syracuse level operator playing the role of `1 - ω` | every `C_k` fixed point came with an explicit self-map (§2.5) | the level operator (THM-4520, [circulant note](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)) has spectrum on a half-circle, not a fixed graph; ANALOGY only |
| **P7 two-regime cycles** | Baker for few runs, an exact residue argument for many runs | the first regime is a theorem (Steiner; Simons–de Weger; Hercher) | the second has no method: equidistribution is the wrong direction (3.4) |

None closes; P1 and P2 are the two with an exact new object behind them (the
mediant tree, the set `K`), and D40–D41 say what to compute next.

## 4. Verdicts

| claim | status |
|---|---|
| Lemma 1 (monotonicity, additivity) | PROVED |
| `C_5(P_13)` Cayley on `Z_13 ⋊ Z_3`; no induced `P_13`; icosahedra without hats; `C_5^2(P_13)` contains 13 lifted induced icosahedra with 108 hats each | FINITE-EXACT |
| Theorem A: `P_13` pentagon-expanding | PROVED conditional on the paper's Theorem 3.5 (CITED; its first step verified) |
| Theorem B1 (links are the contractible induced `k`-cycles; `C_k = R`) | PROVED (elementary; `k = 5` via the icosahedron) |
| Theorem B2 (Eisenstein tori fixed iff `3 ∤ n`); `P_13`, `C_37, C_43, C_49, C_61, C_67`, `Z_7^2`, `Z_8^2` fixed; honeycomb preperiodic | PROVED + FINITE-EXACT |
| Theorem B3 (square tori `C_4`-fixed); census `k = 4` | PROVED + FINITE-EXACT |
| winding fixed points `C_7(1,2)`, `C_8(1,2)`, `C_14(1,2,3)`; D36 refuted | FINITE-EXACT; MISTAKES |
| Hadamard identity; weighted-mediant theorem; convexity and one-run corollaries | PROVED |
| rational-cycle census; near-misses at distance `1/D`; no unimodular pairs; `K ∩ Z` to `10^5`; local uniformity of carries | FINITE-EXACT |
| proposals P1–P7 | typed and priced; none closes |
| Collatz | OPEN |

## 5. Directions

* **D39.** Apply the hatted-icosahedron test to `C_5^2(P_17)` and
  `C_5^2(Clebsch)`; and find, for `I_1`, a level `m` at which two disjoint
  non-adjacent copies appear (or prove they never do), to make Theorem A
  independent of the paper's edge-removal induction.
* **D40.** The mediant tree: classify the words whose fixed point has
  denominator `≤ 13` (the small rational cycles) by their splits; test whether
  the denominator of `x_{uv}` is determined by those of `x_u`, `x_v` and the
  clocks (a "twisted Farey" law for denominators).
* **D41.** The set `K`: compute the 2-adic distance from `K` to the positive
  integers below `2^N` as a function of `N` (the exit-time distribution) and
  compare with the record `35655 ↦ 85`; a divergence proof must show this
  distance is bounded below by a function of the size, which is the size
  price in another coordinate.
* **D42.** Fixed points of `C_7`: the census stopped at `k = 6`. A first
  hyperbolic test is negative: among the `2997` degree-7 connection sets of
  `S_4`, `96` give locally-`C_7` Cayley graphs (Klein's `{3,7}` map on 24
  vertices is one, `PSL(2,7) ⊃ S_4` acting regularly), in four isomorphism
  classes, and every one has more than `120` induced 7-cycles (short
  non-contractible cycles), so none is a `C_7` fixed point
  (`_klein_farey.out`). Larger hyperbolic quotients with systole `> 7` and the
  winding family (`C_n(1..m)`) remain.
* **D43.** The Cayley-graph census beyond degree 6, and vertex-transitive
  non-Cayley graphs (Petersen-like) for `k = 5`: is the 4-antiprism the only
  sporadic `C_5` fixed point?
