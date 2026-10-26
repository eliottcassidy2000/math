# Cyclic triples, decorated prime composition, and necklace metagraphs

2026-10-26. **PROVED relative to the cited modular reconstruction theorem:** the
first failure of cyclic-triple determination of automorphism order is exactly
12 vertices; a complete infinite family of cyclic-triple fibres is classified
below, including its groups, move multiplicities, spectrum and pointed growth.
**FINITE-EXACT:** the explicitly bounded checks in section 7. No claim of a
complete census of all prime tournaments or all tournament metagraphs is made.
No external novelty claim is made. Local theorem labels belong to this file.

The user's emblem is also a correction to
[THM-4803, intransitivity/modular-tree census](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md).
Its original assertion that the cyclic-triple hypergraph determines `|Aut(T)|`
is false. Its labelled fibre formula and finite census remain intact.
The concurrent [composition classification](../tournament_composition_classification_20261026/README.md)
has been integrated: it proves the stronger prime-specific `4p` boundary,
uniqueness of the bad order-12 fibre, and full path-cover invariance. This
packet contributes complementary complete carousel fibres, quotient spectra
and pointed growth; its order-12 calculation is an independent recovery.

## 1. Inheritance and the right classification object

The closest proved mechanism is the unique strong-module decomposition and
the substitution stabilizer formula, recovered in
[session56's tournament hierarchy](../session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md).
The canonical hostile is `C3[A+,A+,A+]` versus `C3[A-,A+,A+]` below. The corrected
near miss is THM-4803's inference from `Aut(Q*)=Aut(Q)` to invariance under
independent moves in descendants. The least-used sidecar is the actual
isomorphism-type colouring of the quotient's children, together with the
attachment orbit when size changes.

**CITED inputs.** Boussaïri–Chergui–Ille–Zaidi,
[arXiv:1805.04607](https://arxiv.org/html/1805.04607), Theorem 10 gives the shared
strong modules of a realizable hypergraph and its realizations. Theorem 24
(attributed there to Boussaïri–Ille–Lopez–Thomassé, 2004) gives prime
reconstruction up to converse. Proposition 51 describes independent quotient
choices; Corollary 57 counts them. These are imported results, not new proofs.

Write `C3(T)` for the **labelled set** of cyclic three-element subsets; a
module is a set seen uniformly by every vertex outside it. At a prime
strong-module node, write the tournament as `Q[T_1,...,T_r]`. Then

    1 -> product_i Aut(T_i) -> Aut(T) -> Stab_Aut(Q)(i -> [T_i]) -> 1,       (1)
    |Aut(T)| = product_i |Aut(T_i)| * |Stab_Aut(Q)(i -> [T_i])|.

The sequence splits after choosing identifications between isomorphic children.
Its splitting is not canonical. Every automorphism permutes maximal strong
modules, inducing the displayed type-preserving quotient automorphism;
conversely such a quotient map lifts using child isomorphisms. The kernel is
exactly the product. At a linear node the children have a fixed transitive
order and the automorphism group is just their product.

This supplies a recursive classification: keep the ordered child signatures
at linear nodes, and the `Aut(Q)`-orbit of actual child signatures at prime
nodes (with quotient isomorphisms when comparing different representatives).
Leaves are singletons. Induction proves two signatures agree exactly when
the tournaments are isomorphic. Forgetting order and prime-converse choices
gives the coarser `C3` observer; it does **not** preserve the stabilizers in (1).

This is the useful part of the simple-group analogy: composition has canonical
indecomposable quotients and a recursive isomorphism test. A prime tournament
need not have a simple automorphism group: the prime carousel `Q_9` below has
group `C_9`. A list of prime factors or child sizes omits essential action data.

## 2. The emblem, its entire fibre, and the sharp order boundary

Put `A+ = C3 -> 1` and `A- = 1 -> C3`. On fixed labels `0,1,2,3`, keep
`0->1->2->0`; the singleton is respectively a sink or a source. Both blocks
have exactly the cyclic triple `{0,1,2}` and automorphism group `C_3`, but are
nonisomorphic. Use three consecutive four-vertex blocks, with quotient
`0->1->2->0`. Let `T_w = C3[A_(w_0),A_(w_1),A_(w_2)]`, where `0=+`, `1=-`.

**Theorem A (PROVED; literal arc verification supplied).** All these tournaments
have the same 67 labelled cyclic triples: three internal triangles and all
`4^3=64` triples with one vertex in each block. A triple meeting exactly two
blocks is transitive. Nevertheless:

| Word representative | Automorphism group | Order | Labelled realizations of this class | `H` | `H/|Aut|` |
|---|---|---:|---:|---:|---:|
| `000` | `C_3^3 semidirect C_3` | 81 | 16 | 462915 | 5715 |
| `100` | `C_3^3` | 27 | 48 | 462915 | 17145 |
| `110` | `C_3^3` | 27 | 48 | 462915 | 17145 |
| `111` | `C_3^3 semidirect C_3` | 81 | 16 | 462915 | 5715 |

The semidirect product cyclically permutes the three factors; it has 36
elements of order 9, whereas `C_3^3` has none. Thus even the group exponent
changes. The complete labelled fibre has 128 tournaments and exactly the four
ordinary isomorphism classes shown, not merely the two examples. The `H`
entries are **FINITE-EXACT**, independently checked by subset path DP and
inclusion-exclusion of walks. The other entries also follow from section 3.

**Theorem B (PROVED).** If tournaments `S,T` have equal labelled cyclic triples
and at most 11 vertices, then `|Aut(S)|=|Aut(T)|`. The bound is sharp by A.

*Proof.* A tournament automorphism has no even vertex cycle: if it had a cycle
of length `2a`, its `a`th power would swap an antipodal pair, reversing its
arc. Hence every nontrivial vertex cycle has length at least three.

Induct on the common order. The cited reconstruction theorem gives a common
strong-module tree. At a linear root, (1) is a product and induction applies.
At a prime root, the two quotients are identical or converse on the same
child slots; these have the same permutation automorphism group `G`. Restrict
to its subgroup `G_s` preserving child sizes. If an element of `G_s` moves a
child of size at least four, its cycle uses at least three such children, so
the total order is at least twelve. Thus for order at most eleven every large
child slot is fixed by **every** element of `G_s`. Children of size at most
three have their isomorphism type fixed by their `C3` structure (the only
three-vertex choices are transitive or cyclic). Consequently the subgroup
preserving all child types is the same in both realizations. Induction gives
equal child automorphism orders, and (1) finishes the proof. Singleton roots
are the base case. ∎

**Stronger current result (PROVED in the integrated source):**
[fibre F2](../tournament_composition_classification_20261026/fibre/RESULT.md)
proves `v_p(|Aut|)` invariant for `n<4p`, sharply for every odd prime `p`.
Its refinement of this argument restricts a Sylow `p`-subgroup to fewer than
`p` mutable child ports; all of them must be fixed. The construction at
`4p` is the constant word versus a word marking one slot in our family at
`m=p`. Thus (2) realizes every sharp prime threshold, while (3) classifies
all of that carousel fibre, including composite-period strata.

The first failed implication in THM-4803 was treating a child-type stabilizer
as unchanged while the child types themselves change. The strong-component
vertex partition survives; its actual tournament types need not. Both
tournaments in A are strong and nonisomorphic, also refuting that assertion.
The strongest order-only survivor is B; for arbitrary size (1) is the repair.

## 3. Complete classification of an infinite family

For every odd `m>=3`, define the carousel `Q_m` on `Z/m` by

    i -> j iff 1 <= (j-i mod m) <= (m-1)/2.

For `w in {0,1}^m`, put `T_w=Q_m[A_(w_i)]`, and denote their common labelled
cyclic-triple hypergraph by `H_m`. This `H_m` is a hypergraph, not a path count.

**Lemma (PROVED).** `Q_m` is prime and `Aut(Q_m)=C_m`.
Indeed each in- and out-neighbourhood is transitive. A proper module lies in
one of these neighbourhoods of an outside vertex, so is transitive. Global
regularity makes all its internal outdegrees equal, forcing it to be a
singleton. Rotations act transitively; an automorphism fixing a vertex fixes
its two linearly ordered neighbourhoods pointwise. The stabilizer is trivial.

**Theorem C (PROVED, completeness uses the cited reconstruction theorem).**
The ordinary isomorphism classes of *all* realizations of `H_m` are exactly
the binary necklaces of length `m`, i.e. words modulo rotation, not reflection
and not complementation. If `d|m` is a word's least period, then

    Aut(T_w) = C_3^m semidirect C_(m/d),
    |Aut(T_w)| = 3^m m/d.                                                (2)

The cyclic factor acts by rotation through `d` child slots. Every divisor
`d|m` occurs: repeat any primitive length-`d` binary word (a single 1 works
for `d>1`). The possible automorphism orders are precisely `3^m s`, `s|m`.

*Completeness and boundary data.* The root and each internal triangle are
prime; each block has a linear node with two children. Thus the complete
labelled fibre has `2^(m+1) 2^m = 2^(2m+1)` elements. The internal triangle
directions are irrelevant to unlabelled block type. Reflection converts the
converse root back to `Q_m`. Once that root is fixed, only its rotations can
identify words. Formula (1) now gives (2). The hypergraph automorphism group
is `S_3^m semidirect D_(2m)`, where `D_(2m)` has order `2m`: the shared strong
modules force this action, and prime reconstruction restricts the root to
automorphisms or anti-automorphisms of `Q_m`. All these permutations work.
Hence a period-`d` class has exactly `2^(m+1)d` labelled realizations.

Global converse sends a word to its complemented reflection. No necklace of
odd length is fixed: a reflection of an odd polygon has a fixed position,
where the bit would have to equal its complement. Therefore every `T_w` is
non-self-converse. Converse pairs are distinct vertices until explicitly
quotiented; that quotient is not being taken here.

The number of classes of least period `d` and the total number are

    P_2(d) = (1/d) sum_(e|d) mu(e) 2^(d/e),
    N_2(m) = (1/m) sum_(e|m) phi(e) 2^(m/e).                              (3)

These follow respectively by Möbius inversion of the primitive-word
decomposition and by counting fixed words of rotations. They imply the
weighted mass identity

    sum_[w] 1/|Aut(T_w)| = 2^m/(m 3^m).                                  (4)

| `m` | Tournament order | Classes | Classes by least period |
|---:|---:|---:|---|
| 3 | 12 | 4 | `1:2, 3:2` |
| 5 | 20 | 8 | `1:2, 5:6` |
| 7 | 28 | 20 | `1:2, 7:18` |
| 9 | 36 | 60 | `1:2, 3:2, 9:56` |
| 11 | 44 | 188 | `1:2, 11:186` |
| 13 | 52 | 632 | `1:2, 13:630` |

The arithmetic is structural: divisors measure stabilizers, primitive periods
measure symmetry breaking, and totients count rotation types. Size alone is
insufficient; a composite `m` introduces additional symmetry strata.

## 4. A fully specified metagraph and its spectrum

**Theorem G (PROVED; general prime-quotient reduction).** Fix any prime
tournament `Q` of order `r>=3`, and let `G=Aut(Q)` act on its vertex set.
The fixed-root family `Q[A_(w_i)]` is classified by `G`-orbits on `{0,1}^r`.
Its automorphism orders are `3^r |Stab_G(w)|`. Its block-cut flip metagraph
is the `G`-quotient of the `r`-cube, retaining all physical move multiplicities.
Its eigenvalue `r-2k` has multiplicity

    b_G(k) = [x^k] (1/|G|) sum_(g in G) product_(cycles c of g)(1+x^|c|).  (8)

Its total class count is `(1/|G|) sum_g 2^(number of cycles of g)`.
Proof: (1) proves the isomorphism and group assertions; orbit counting and
the Walsh proof of E give the rest, without requiring that `G` be cyclic.
Orbit-size edge balance gives stationary mass proportional to `1/|Aut|`.
For the ordinary classes of the *whole* labelled `C3` fibre, the
converse-root choice supplies one such slice if `Q` is self-converse,
and two disjoint ordinary-class slices
otherwise. A root-converse move is an additional specified move; the
block-cut graph alone does not connect the two slices. These statements
classify this decoration family for every prime quotient, not all possible
child tournaments or all possible metagraph operations.

**FINITE-EXACT noncyclic control.** On `Z/7` take arcs with differences
`{1,2,4}`. Exhausting all 7! permutations and all subsets verifies primeness,
self-converseness and `|G|=21`. Its 28-vertex decorated fibre has 12 ordinary
classes, with automorphism orders `45927` (2 classes), `6561` (6), and
`2187` (4). The multiplicities for weights `0..7` are `1,1,1,3,3,1,1,1`.
The full 12-dimensional Walsh eigenbasis and cycle-index formula agree.

Define `M_m` to have the necklaces of Theorem C as vertices. A physical move
flips **one block cut**, reversing its three triangle-to-singleton arcs. Keep
the directed multiplicity `c(u,v)`: how many of the `m` positions of any
representative of `u` yield `v`. This is the rotation quotient of the binary
`m`-cube. It is a connected, loop-free, bipartite multigraph with weighted
degree `m`. There is no intrinsic tournament orientation here; forcing one
would discard reversible moves. This is a modular-flip metagraph, not an
identification with the repository's Hamilton-path-complement black/blue graph.

**Theorem D (PROVED).** If `d(u)` is least period, then

    d(u)c(u,v) = d(v)c(v,u),
    c(u,v)/|Aut(T_u)| = c(v,u)/|Aut(T_v)|.                                 (5)

Count cube edges between the two rotation orbits to prove the first identity;
(2) gives the second. The uniform-position walk has stationary probability
`d(u)/2^m`. It has period two, so stationarity does not imply convergence
from a single starting vertex without making the walk lazy.

For `m=3`, order the vertices by weights `0,1,2,3`. The matrix `C=(c(u,v))` is

    [0 3 0 0]
    [1 0 2 0]
    [0 2 0 1]
    [0 0 3 0]

and stationary masses are `(1,3,3,1)/8`. The unweighted graph is a path;
discarding the multiplicities changes its operator and spectrum.

**Theorem E (PROVED).** The eigenvalues of `C` are `m-2k`, `0<=k<=m`, with
multiplicities

    a_m(k) = (1/m) sum_(e|gcd(m,k)) phi(e) binom(m/e,k/e).                 (6)

Here `gcd(m,0)=m`. The characteristic polynomial is the product of
`(x-(m-2k))^a_m(k)`. In particular the two weight-one/weight-`m-1` eigenspaces
have multiplicity one, and the random-walk eigenvalues are `1-2k/m`.

*Proof.* The cube Walsh characters indexed by subsets `S` have eigenvalue
`m-2|S|`: flipping a coordinate changes the sign exactly for coordinates in
`S`. Rotations permute these characters. Sums over subset-rotation orbits
form a basis of the invariant subspace, which is exactly the quotient
function space. Thus the multiplicity is the number of weight-`k` necklaces.
Burnside's lemma gives (6): a rotation of cycle length `e` fixes such a subset
only if `e|k`, then chooses `k/e` of its `m/e` cycles. ∎

This exact spectrum is specific to the multiplicity-retaining modular flip.
It says nothing by itself about spectra or connectivity of a different
metagraph. In particular, the common path count in A does not determine
`H/|Aut|`: the Hamilton-path presentation orbit counts differ by a factor of
three. An automorphism fixing a full ordered Hamilton path fixes all vertices,
so that division is exact, not a heuristic symmetry normalization.

## 5. Growth needs a marked attachment

**Theorem F (PROVED).** Write `m=2k+1`. There is an induced embedding
`Q_m -> Q_(m+2)` by

    f(i)=i+1 for 0<=i<k;  f(i)=i+2 for k<=i<m.

The omitted positions are `0,k+1`. For a pointed word define

    grow(w) = 0, w_0,...,w_(k-1), 0, w_k,...,w_(m-1).                     (7)

Then `T_w` induced-embeds in `T_grow(w)`, adding eight vertices.
For `i<j`, their difference either stays `j-i` or becomes `j-i+1`.
If it stays unchanged, the same-side maximum difference is `k`; if it
crosses the inserted gap, the old threshold `j-i<=k` becomes the new
threshold `j-i+1<=k+1`. Thus every old arc is preserved. Whole blocks then
give the stated embedding of substitutions.

The point is necessary: `011` and `101` are the same necklace, but (7) gives
`00011` and `01001`, distinct length-five necklaces. The growth map therefore
does not descend to unpointed isomorphism classes. A marked cyclic gap or its
stabilizer orbit is a sufficient attachment sidecar.

Symmetry need not increase with size. Start from the all-ones word. Its group
order is `m 3^m`. The child has exactly two zeros and odd length `m+2`; any
nontrivial repetition count would divide both two and `m+2`, so its least
period is `m+2`. The child group order is `3^(m+2)`, a ratio `9/m`. At `m=11`
this gives an induced embedding on `44 -> 52` vertices with group orders
`1,948,617 -> 1,594,323`. No extension of all parent symmetries is implied by
an induced embedding.

## 6. Research board and recovered connections

The maintained board, the seeded random sample, explicit connection maps,
hostiles and stopping conditions are in [RESEARCH.md](RESEARCH.md). The anchor
is decorated prime composition; the niche is arithmetic orbit strata and
the modular-flip operator; the wildcard is pointed growth and transport of
finite-state boundary data. LRC(14) remains **OPEN**. The sampled clock result
supplies a comparison of cyclic actions, not a transfer of an LRC inequality.

The next classification layer is **OPEN**: arbitrary prime quotients with
several child types, their stabilizer strata and compatible growth ports.
Formula (1) gives the exact reduction, but a prime census and a specified
move functor are still required. Transfer to Hamilton-path complements needs
the endpoint/path-cover response at each block; the cyclic-triple tree or
group order alone cannot replace it. No new universal meta-pattern is promoted
from this single session.
The integrated [metagraph M4-M5](../tournament_composition_classification_20261026/metagraph/RESULT.md)
gives a sharper existing hostile for that different operation: two half-edges
of one order-four blue loop acquire singleton extensions with `H=11` and
`H=13`. Our gap-marked growth is a separate operation with the same explicit
need to retain the attachment orbit, not an identification of these graphs.

## 7. Exact universe, controls and reproduction

Run from the repository root with Python 3, standard library only:

    python3 05-knowledge/results/cyclictriple_composition_20261026/verify.py

`check.py` exhausts all 128 labelled realizations supplied by the cited fibre
theorem, assigns each by actual arc isomorphism to exactly one of four
representatives, and counts their actual automorphisms without assuming a
module decomposition. It also exhausts **all labelled tournaments** of orders
1 through 5 and their permutation automorphisms; all binary words at odd
lengths 3 through 13; and all physical flip ports in their necklace quotients.
It checks every candidate module of `Q_3,Q_5,Q_7,Q_9`, all quotient arcs in
the growth map for odd orders through 25, and all substituted growth arcs
for every binary word of lengths 3 and 5. These finite tests do not replace
the all-order proofs or claim an enumeration through order 11.

`audit.py` imports no producer code or receipt. It builds arcs by independent
coordinate rules, checks all `6*6^3` module-compatible permutations on each
of the emblem's two representatives, and counts paths by an independent
walk inclusion-exclusion formula. For `m=3,5,7,9` it constructs the quotient
matrix and every Walsh orbit-sum eigenvector, verifies the eigen-equations
exactly, and proves full column rank modulo 1,000,003. It checks the period
counts by a divisor recurrence without the producer's Möbius routine.
It separately exhausts the 7! permutations of the difference-`{1,2,4}`
quotient, all its subsets and binary decorations, and the complete
12-dimensional quotient eigensystem for the noncyclic control of G.

Positive controls include constant words (full root rotations), aperiodic
words (no root rotations), and the direct small-order tournament census.
Hostiles include the 81/27 split, absence of order-nine elements in the second
group, non-self-converse necklaces, loss of the operator when port counts
are erased, and the two different growth children of one parent class.
Both scripts pass normally and under `-O`; each pair of receipts is byte-identical.
This is an independent computational path by the same session, not an
independent human or agent proof audit. `SHA256SUMS` pins the delivered files.
