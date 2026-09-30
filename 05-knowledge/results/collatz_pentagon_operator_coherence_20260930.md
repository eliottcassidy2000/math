# The pentagon graph operator and the Collatz trichotomy: the Platonic graphs under the induced-cycle operators (the tetrahedron and the icosahedron are the fixed points, the dual pairs are operator steps, Petersen's pentagon graph is `K_12` minus the co-polar matching and then vanishes), the carry as a 1-cocycle with two hexagon identities (the commutation defect of two words is the product of their clocks times the difference of their fixed points), integer cycles as coboundaries (`H^1(<w>; Z) = Z/(2^A - 3^p)`), and Mac Lane's coherence against Terras and Conway

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, twelfth note), 2026-09-30.
**Owner's directive:** "Pursue promising direction and new ones as they
emerge. Just as tournament classification evaluates how 3-vertex and 5-vertex
subgraphs structure a massive directed network, the triangle and pentagon
axioms act as the base cases for monoidal categories. Mac Lane's Coherence
Theorem states that if you enforce structural consistency at the level of 3
elements (the Triangle Axiom) and 4 elements mapped across 5 operations (the
Pentagon Axiom), all higher-order configurations are automatically consistent
and isomorphic. Think how the triplet of {vanishing, periodic, expanding}
relates to other abstract things of 3 we have been studying, and possibly
creatively draw connections between the hexagon identities in braided
monoidal categories relate deeply with the work we have done, especially with
the number 6 appearing in primes are related to collatz", with the attached
paper S. V. Gervacio, H. Maehara, P. C. F. Ramos, *The Pentagon Graph
Operator*, arXiv:2604.18984v1 (21 Apr 2026).
**The paper.** `C_5(G)` has the induced 5-cycles of `G` as vertices, two
adjacent iff they share an edge. Every graph is exactly one of
pentagon-vanishing (`C_5^k(G) = ∅` for some `k`), pentagon-periodic, or
pentagon-expanding (`|V(C_5^k(G))| → ∞`), a consequence of Gervacio's
fundamental theorem on graph operators (their [4]). `C_5(dodecahedron) ≅
icosahedron` and `C_5(icosahedron) ≅ icosahedron`; an icosahedron with `h`
"tadpole hats" (a new vertex joined to an induced triangle-plus-pendant) is
expanding, each hat begetting two, so `≥ h 2^k` hats after `k` steps. The
paper asks which `k` admit periodic polyhedral fixed points for the operators
`C_k`, and which graphs expand.

**Status: PROVED (elementary) for Propositions 1–4; FINITE-EXACT for every
trajectory, the Platonic table and the Petersen identification; CITED for the
paper, van Rooij–Wilf, Mac Lane, Conway, Terras and the repo results named by
path; ANALOGY and NUMEROLOGY typed where used; Collatz OPEN. Independent
audit OWED (weekly subagent limit; self-checks recorded: the cocycle and
hexagon identities on 400 random word triples and on every split of the `-17`
word, the isomorphisms by brute force).** Script
`04-computation/experiments/collatz_pentagon_operator_20260930.py`, output
alongside.

**Reconciliation.** The parallel note of the same day,
[`collatz_circulant_20260930_circulants_lucas_cubic_monotile.md`](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)
(THM-4520 and the Pillai–Lucas census), answers the eleventh note's directive
from the circulant side; an addendum to the
[eleventh note](collatz_lucas_monotile_discrepancy_20260930.md) now cites it.
This note does not touch it.

## 1. The trichotomy is the Collatz trichotomy; where it is a theorem, a Lyapunov function does it (ANALOGY, both sides exact)

For a self-map of a countable set an orbit is either eventually periodic or
leaves every finite set, so a Collatz orbit is *vanishing* (absorbed by the
cycle `{1, 2}`, the analogue of the empty graph, which `C_5` fixes),
*periodic* (a non-trivial cycle) or *expanding* (divergent: a non-periodic
orbit takes distinct values, hence `|T^k n| → ∞`). The conjecture over `N`
is "every orbit vanishes"; the paper's Theorem 4.1 is the same elementary
trichotomy for graphs. What the graph world has and Collatz lacks is a
*classified* instance: for the line graph operator, van Rooij and Wilf (1965,
the paper's [10]) proved that paths vanish, cycles and `K_{1,3}` are periodic
(`K_3` is fixed) and every other connected graph expands, the proof being the
monotone quantity `|E(L(G))| = Σ_v C(deg v, 2) ≥ |E(G)|` with an exact
equality case (all degrees `≤ 2`) — a Lyapunov function whose level sets are
the periodic graphs. The pentagon operator has no such quantity, and the
paper offers examples on both sides (periodic: the golden solids; expanding:
a doubling gadget) and two open classification problems. In the
[tenth note](collatz_two_carries_typology_20260927.md)'s typology,
provability tracked an exact invariant or monotone quantity; the line graph
operator is on the provable side, the pentagon operator and Collatz on the
other, for the same reason.

## 2. What the operator does to the repo's graphs (FINITE-EXACT)

| graph | `|V(C_5^k)|`, `k = 0, 1, ...` | fate |
|---|---|---|
| `C_5 = P_5`, and the `5x+1` trivial cycle `1, 3, 8, 4, 2` | `5, 1, 0` | vanishing (the paper's first example) |
| `K_5`, `K_{3,3}` (the Kuratowski pair), `K_6 = I/±1` | `5, 0`; `6, 0`; `6, 0` | vanishing at once (no induced pentagon) |
| Heawood, Desargues (girth 6) | `14, 0`; `20, 0` | vanishing at once |
| Petersen `= D/±1` | `10, 12, 0` | vanishing in two steps |
| Seidel double `D(P_5)` (eighth note's tower) | `10, 12, 1, 0` | vanishing in three steps |
| dodecahedron `D` | `20, 12, 12, 12, ...` | periodic from step 1 |
| icosahedron `I` | `12, 12, 12, ...` | fixed |
| Paley `P_13`, `P_17`; Clebsch; `D^2(P_5)` | `13, 39, 4563`; `17, 272`; `16, 192`; `20, 384` | growing (expansion not proved) |

**Proposition 1 (Petersen).** The Petersen graph has twelve induced
pentagons; two share `0, 1` or `2` edges in `6, 30, 30` of the `66` pairs, the
six edge-disjoint pairs being a perfect matching; so `C_5(Petersen) = K_12 -
6K_2` (the cocktail-party graph `K_{2,2,2,2,2,2}`), which has no induced
pentagon (its non-edges form a matching, and the complement of an induced
`C_5` is a `C_5`), so `C_5^2(Petersen) = ∅`. *Proof.* Computation for the
counts; the last two sentences are the argument. ∎

**The folds.** The [coincidence atlas](collatz_coincidence_atlas_20260927.md)
read the golden solids through their antipodal folds `D → Petersen`, `I →
K_6`. The operator sees the double cover, not the fold: `C_5(D) = I` is
periodic while `C_5(D/±1) = K_12 - 6K_2` vanishes; `C_5(I) = I` while
`C_5(I/±1) = ∅`. The six non-edges of `C_5(Petersen)` are the six co-polar
pairs of the paper (the vertices of `I` at distance 3, the antipodal face
pairs of `D`): the fold keeps exactly the co-polar structure and destroys
the rest.

**A correction to the paper's wording.** The map `x ↦ [N_I(x)]` (the link
pentagon of `x`) is *not* an isomorphism `I → C_5(I)`: the links of adjacent
vertices share two vertices (the apexes of the two triangles on the edge)
but no edge, and the links of vertices at distance 2 share the edge formed
by their two common neighbours. So `C_5(I)` is the distance-2 graph of `I`,
which is isomorphic to `I` through the co-polar involution (`d(x, y) = 2` iff
`d(x, y') = 1` for `y'` co-polar to `y`); the isomorphism is `x ↦ [N_I(x')]`.
The theorem `C_5(I) ≅ I` stands (verified by brute force); the stated map is
off by the antipodal twist — the same twist the folds exhibit.

**Proposition 2 (the Platonic graphs under `C_3, C_4, C_5`).**

| solid | `C_3` | `C_4` | `C_5` |
|---|---|---|---|
| tetrahedron `K_4` | `K_4` (fixed) | `∅` | `∅` |
| octahedron | cube | `3K_1` then `∅` | `∅` |
| cube | `∅` | octahedron | `∅` |
| icosahedron | dodecahedron | `∅` | icosahedron (fixed) |
| dodecahedron | `∅` | `∅` | icosahedron |

When the faces of a polyhedral graph are induced `k`-gons and are all its
induced `k`-cycles, `C_k(P) = P^*` (the dual: faces adjacent iff they share
an edge); so the two dual pairs appear as operator steps, the self-dual
tetrahedron is the fixed point of `C_3`, and the icosahedron is the fixed
point of `C_5` because its twelve induced pentagons are the vertex links
rather than faces. The octahedron/cube pair vanishes under both `C_3` and
`C_4` (the octahedron's induced squares are its three equators, pairwise
edge-disjoint). *Proof.* Computation; the dual statement is the definition
of the dual. ∎ In the [eighth note](collatz_doubling_tower_brocard_20260927.md)'s
triple dictionary KT-b ("dual pair + self-dual"), the Platonic graphs sort
under the operators as `{tetrahedron}` self-dual fixed, `{octahedron, cube}`
dual pair vanishing, `{icosahedron, dodecahedron}` dual pair periodic. For
the paper's question, among the Platonic graphs the operators `C_k` have
polyhedral fixed points exactly at `k = 3` (`K_4`) and `k = 5` (`I`), and none
at `k = 4`.

**A remark on multipliers (NUMEROLOGY, exact).** The trivial cycle of the
shortcut map for `qx + 1` through `1` is a digon for `q = 3` (`1, 2`), a
pentagon for `q = 5` (`1, 3, 8, 4, 2`) and a triangle for `q = 7` (`1, 4, 2`).
Every `3x + 1` functional graph (over `N`, or over `Z` with its `-5` triangle
and `-17` eleven-cycle) has no induced pentagon and is pentagon-vanishing at
once; only the `q = 5` sibling's trivial cycle is a pentagon.

## 3. The carry is a 1-cocycle; the two hexagons; cycles are coboundaries (PROVED, elementary)

Write a word `w = (v_1, ..., v_p)` for `p` odd Syracuse steps with valuations
`v_i ≥ 1`, `A_w = Σ v_i`, `p_w = p`, carry `S_w = Σ_{i<p} 3^{p-1-i} 2^{d_i}`
(`d_i = v_1 + ... + v_i`, `d_0 = 0`), clock `D_w = 2^{A_w} - 3^{p_w}`, and fixed
point `x_w = S_w/D_w` of the affine map `f_w(x) = (3^{p_w} x + S_w)/2^{A_w}`.

**Proposition 3 (cocycles).** For all words `u, v`:
`S_{uv} = 3^{p_v} S_u + 2^{A_u} S_v` and `D_{uv} = 2^{A_u} D_v + 3^{p_v} D_u`.
Hence `S` is a derivation (1-cocycle) on the free monoid of words with
values in `Z` made a bimodule by `w · c = 3^{p_w} c` on the left and
`c · w = 2^{A_w} c` on the right; it is the unique such derivation with
`S(letter) = 1`. *Proof.* `f_{uv} = f_v ∘ f_u` and multiply out; the clock
identity is `2^{A_u}(2^{A_v} - 3^{p_v}) + 3^{p_v}(2^{A_u} - 3^{p_u})`. ∎ (The
one-letter case is the rotation formula of the
[necklace note](collatz_necklace_20260929_fair_splits_power_clocks_basins.md);
word reversal gives S20's carry reciprocity.)

**Proposition 4 (the commutation defect and the two hexagons).** Let
`β(u, v) = S_{uv} - S_{vu}`. Then
`β(u, v) = D_u D_v (x_v - x_u)`, `β(v, u) = -β(u, v)`, and
`β(u, vw) = 3^{p_w} β(u, v) + 2^{A_v} β(u, w)`,
`β(uv, w) = 3^{p_v} β(u, w) + 2^{A_u} β(v, w)`.
So two words commute (their two orders give the same map) iff they have the
same fixed point; the defect with a product is the twisted sum of the
defects with the factors, one hexagon for the `3`-side (left action) and one
for the `2`-side (right action). *Proof.* From Proposition 3,
`S_{uv} - S_{vu} = (3^{p_v} - 2^{A_v}) S_u - (3^{p_u} - 2^{A_u}) S_v = D_u S_v
- D_v S_u`, and `S_w = D_w x_w`. For the first hexagon, expand `D_u D_{vw}
(x_{vw} - x_u)` with `D_{vw} x_{vw} = S_{vw} = 3^{p_w} S_v + 2^{A_v} S_w` and
`D_{vw} = 2^{A_v} D_w + 3^{p_w} D_v`; the second is symmetric. ∎ (Checked on
400 random triples, script part 2.)

**What this says about braiding.** A braided monoidal structure would be a
commutation isomorphism `u ⊗ v → v ⊗ u` satisfying the two hexagons; here
the "commutation" of two words is the translation by `β(u, v)/2^{A_u + A_v}`,
antisymmetric, so it is a *symmetry* (order two), the Yang–Baxter equation
holds trivially, and the two hexagons are Proposition 4. There is no
braided-but-not-symmetric structure in the word monoid: the whole
non-commutativity of "halve" and "triple" is the scalar `D_u D_v (x_v - x_u)`.
The "6" the owner asks about is the bimodule: `3` acts on the left, `2` on
the right, one hexagon each; and every clock is a Syracuse-core residue,
`2^A - 3^p ≡ (-1)^A (mod 6)` (odd, and `≡ (-1)^A mod 3`), the sign being the
parity of the halving count (script part 4).

**Proposition 5 (cycles are coboundaries).** Restricted to the cyclic
submonoid `<w>`, the carry is a coboundary over a ring `R` iff `S_w ∈ D_w R`,
and `H^1(<w>; R) = R/D_w R`. Over `Z` this is `Z/(2^A - 3^p)`, the fixed-point
group of the [eleventh note](collatz_lucas_monotile_discrepancy_20260930.md)'s
Proposition 1, and the class of the carry is that note's torsion point: an
integer cycle of word `w` exists iff the class vanishes, the cobounding
element being the cycle point (`S_w = x_w D_w`; for `-17`, `2363 = (-17)(-139)`).
Over `Z_2` and `Z_3` the clock is a unit, `H^1 = 0`, and every word is a
cycle (the 2-adic and 3-adic periodic points); the obstruction lives at the
primes dividing the clock, which the eleventh note's Proposition 2 locates in
a lattice. *Proof.* Derivations on `<w> ≅ N` are free on `S(w)`; inner ones
are `w · c - c · w = (3^{p_w} - 2^{A_w}) c`. ∎ A hexagon-shaped necessary
condition follows: for every split `w = uv` of a cycle word,
`D_u D_v (x_v - x_u) ≡ 0 (mod D_w)`; all six splits of the `-17` word satisfy
it (script part 2, the defects `-1112, -2780, -5282, -3336, -6116, -10286`
are multiples of `139`). It is equivalent to rotation invariance of the cycle
condition and adds no power; it is recorded because it is the exact form the
owner's hexagon takes.

## 4. Mac Lane's coherence against Terras and Conway (ANALOGY; the bijection typed)

Coherence says: enforce the triangle and pentagon axioms and every diagram
of associators and unitors commutes (Mac Lane 1963; the pentagon is the
associahedron `K_4`, the `5 = C_3` bracketings of four letters, and
Stasheff's `K_n` carry the higher cases). Two local base cases propagate to
all sizes. The Collatz thread has the opposite structure, exactly:

* **No finite level propagates.** Terras's bijection says every parity word
  of length `N` occurs exactly once modulo `2^N`: the residue class of `n`
  modulo `2^N` fixes its first `N` steps and nothing beyond. Consistency at
  level `N` is total and says nothing at level `N + 1`; there is no
  "triangle and pentagon" whose commutation forces the rest.
* **No finite set of base cases can certify the family.** Conway (1972,
  *Unpredictable iterations*): generalized Collatz maps are Turing
  complete, so no finite local check decides termination across the family.
  Coherence is the statement that finitely many local axioms decide
  everything; Conway's theorem is the statement that for Collatz-like maps
  they cannot in general.
* **The associahedra do appear, as a bijection.** Their vertices are
  Catalan-counted (`1, 2, 5, 14, 42`), as are the critical spine blocks at the
  formal slope `2` of the [seventh note](collatz_catalan_ramsey_20260927.md)
  (a rise followed by a Dyck path): the five bracketings of the pentagon axiom
  correspond to the five critical blocks of length 3 through Dyck paths. A
  bijection of counts, not a mechanism; typed NUMEROLOGY.
* **The base cases `3` and `5` in the sibling family.** In the
  [fourth note](collatz_hedgehog_family_20260927.md)'s family `T_q`, the
  certification density tends to `1` iff `q ≤ 3`, and `q = 5` is the first
  multiplier with a positive-measure never-descending set (`f_∞(5) = 0.177`).
  The triangle and pentagon of the owner's analogy are the two odd multipliers
  on either side of the critical `4 = 2^2`: the last that (conjecturally)
  always vanishes and the first that (conjecturally) expands. Typed
  NUMEROLOGY.

**The owner's triple `{vanishing, periodic, expanding}` against the thread's
triples.** It is the orbit trichotomy itself (absorbed / cycle / divergent),
the hedgehog note's fates, and the axis "fate" of the tenth note's typology
(where provability tracked an exact invariant, as in section 1). Under the
operators the Platonic graphs give the same triple with a twist: the
tetrahedron is fixed, the octahedron/cube pair vanishes, the
icosahedron/dodecahedron pair is periodic — the KT-b dictionary ("dual pair +
self-dual") of the eighth note read as dynamics.

## 5. Verdicts

| claim | status |
|---|---|
| graph-operator trichotomy = Collatz orbit trichotomy; van Rooij–Wilf classified the line graph operator by a Lyapunov function; pentagon operator and Collatz lack one | ANALOGY (CITED theorems) |
| trajectories table; `C_5(Petersen) = K_12 - 6K_2` then `∅`; the folds; `C_5(I)` is the distance-2 graph, isomorphic to `I` through the co-polar involution (the paper's stated map is off by the antipodal twist) | FINITE-EXACT / PROVED |
| Platonic table: `C_k(P) = P^*` for induced `k`-gonal faces; fixed points `K_4` (`k = 3`), `I` (`k = 5`), none at `k = 4` | PROVED (elementary), FINITE-EXACT |
| carry cocycle, clock cocycle, `β(u,v) = D_u D_v (x_v - x_u)`, antisymmetry, two hexagons; the structure is symmetric, not braided | PROVED |
| cycles are coboundaries; `H^1(<w>; Z) = Z/(2^A - 3^p)` = the eleventh note's fixed-point group; `H^1 = 0` over `Z_2`, `Z_3` | PROVED |
| clocks `≡ (-1)^A mod 6`; `5x+1`'s trivial cycle is a pentagon | EXACT, NUMEROLOGY |
| coherence vs Terras/Conway; associahedra vs critical blocks; `3, 5` around `4` | ANALOGY / NUMEROLOGY |
| Collatz | OPEN |

## 6. Directions

* **D35 (expansion of the Paley and Clebsch graphs).** `P_13 → 39 → 4563`,
  `P_17 → 272`, Clebsch `→ 192`: find the doubling gadget (the paper's
  "local attachment to a cycle-rich core") inside `C_5(P_13)`, or a monotone
  quantity, to decide expansion.
* **D36 (fixed points of `C_k`).** Conjecture from the Platonic table: a
  vertex-transitive graph is `C_k`-fixed iff its vertex links are induced
  `k`-cycles that are all its induced `k`-cycles and meet in an edge exactly
  at distance 2; classify for `k = 5` beyond the icosahedron (the paper's
  first open problem).
* **D37 (the cohomology of the word monoid).** `H^1` of the free monoid with
  the `(3, 2)`-bimodule coefficients is large; the cycle question is the
  vanishing of the restriction classes on cyclic submonoids. Is there a
  cohomological operation (transfer along `<w> ⊂ <u, v>`) that relates the
  classes of `w`, of its splits and of its rotations beyond Proposition 4?
* **D38 (the operator on Collatz residue graphs).** Define the graph on
  `Z/2^N` with edges between residues whose words differ in one letter
  (Terras's cube) and study `C_4` and `C_5` on it; the words of a fixed shape
  are the vertices of a Johnson-type graph, whose induced cycles are exactly
  the pentagon-free structures — a place where the paper's operator and the
  thread's word combinatorics coincide by construction.
