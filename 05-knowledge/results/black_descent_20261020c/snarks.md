# Cyclic-five Petersen reductions with seven-port support and Kempe defects

**PROVED:** the support-domination implication and its stability under
gluing; the explicit NN substitution identity and rooted minor construction,
with its finite Petersen-negator base case checked exhaustively; and the
cyclic-five substitution criterion in section 7.
**FINITE-EXACT:** the named-graph certificates, complete three-vertex
deletion censuses, 112 independent closure checks, and the four larger
cyclic-five examples. **CITED:** the standard all-order connectivity of
the underlying flower snarks. **OPEN:** a universal reduction rule for all
cyclically five-connected snarks.

## 1. Inheritance and the change in the reduction target

The closest mechanism is the ordered-boundary gluing identity in
[THM-4116](../../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md),
the rooted minor lemmas in the
[previous snark note](../blue_quotient_20261020b/snarks.md), and the complete
seven-port coordinates in
[junctions.md](../blue_quotient_20261020b/junctions.md).
The inherited hostile is Petersen's four-pole compression to the looped
dumbbell: retaining a colouring obstruction does not ensure that the
closed graph remains simple and bridgeless.

The concurrent [rooted Kempe packet](../defect_transitions_20261020/petersen.md)
and its [integration note](../defect_transitions_20261020/incoming.md) were
read from `origin/main` after commits `cb2874b2e7` and `c75e354bdd`.
Their four-pole and dot-product work overlaps the previous lane.
Their general J-support seven-pole argument isolates a precise additional
obligation: connect the third claw terminal while avoiding the two wire
paths. We do not assume that simultaneous linkage. The new NN carrier
below supplies actual paths, its cut predicate supplies the graph-class
certificate, and double-star uses strict support inclusion beyond the
earlier exact-equivalence grammar.

The anchor here is a reduction between actual cyclically five-connected
snarks. The niche is the distinction between equality and containment of
labelled boundary support. The wildcard is whether a nonnegative tensor
remainder can be realized by another cubic pole.

The concept board is: ordered support, the seven Kempe defects, a rooted
minor model, the unchanged exterior, and the class of the closed graph.
The double-star example below needs all five coordinates. In particular,
its boundary counts cannot merely be scaled to match Petersen.

Let `f_Q(s)` count proper three-edge-colourings of a cubic pole Q with
ordered boundary word s. Write `S_Q={s:f_Q(s)>0}`.
If Q' is a rooted minor of Q and

```text
S_Q' ⊆ S_Q,                                                   (1)
```

then for every exterior R glued in that same order,

```text
Q∪R uncolourable  =>  Q'∪R uncolourable.                       (2)
```

Indeed, a colouring of Q'∪R would give a boundary word in both S_Q' and
S_R, hence in S_Q and S_R, and would extend to Q∪R. This is the useful
direction of support containment. Equality gives equivalence; containment
is enough for descent of an existing obstruction. The rooted minor and
the graph-type check remain separate requirements.

This weaker target is decisive. For the double-star seven-poles considered
below, exact support equality with Petersen is impossible, but (1) holds
for an explicitly certified Petersen minor.

## 2. Named graphs and the exact cyclic-connectivity audit

The two Loupekine inputs are constructed from three copies of Petersen
with a three-vertex path removed. Each resulting five-pole has two
two-port connectors and one residual port. Join them cyclically by their
two-port connectors and attach all residual ports to one new vertex.
Straight and crossed closing connectors give nonisomorphic graphs.

For background and the known NN substitution, see Mazák, Rajník and
Škoviera, [*Morphology of small snarks*, EJC 29(4), P4.30 (2022),
§5.6](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v29i4p30/pdf/).
The paper identifies the two order-22 Loupekine snarks and discusses the
order-30 double-star in its order-by-order classification. The NN
colour-equivalence used below is established there; no novelty claim is
made for that substitution. The double-star adjacency data is transcribed
from the official
[Sage `DoubleStarSnark` constructor](https://raw.githubusercontent.com/sagemath/sage/develop/src/sage/graphs/generators/smallgraphs.py).
Its order is **30**, not 20 or 22.

Every input is independently checked for cubicity and zero colouring
count. Colourings are counted by enumerating perfect matchings in one
colour and even complementary two-factors in the other two colours.

| Graph | Vertices | Girth | Perfect matchings | Tait colourings | Four-edge deletion sets checked |
|---|---:|---:|---:|---:|---:|
| Loupekine 0 | 22 | 5 | 36 | 0 | 40,920 |
| Loupekine 1 | 22 | 5 | 36 | 0 | 40,920 |
| Double-star | 30 | 6 | 124 | 0 | 148,995 |

For every deletion set of zero, one, two, three, or four edges, the script
finds all connected components and tests which contain a cycle. No
deletion leaves two cyclic components. The table reports just the
four-edge cases; every smaller size is tested too. All three graphs also
have explicit cyclic five-cuts: isolate a Loupekine negator, or separate double-star vertices
0–14 from vertices 15–29. Thus their cyclic edge-connectivity is exactly
five.
The [graph audit receipt](snarks_graph_audit.json) records every edge set,
the five-cut witnesses, and the deletion counts separately for sizes 0–4.

## 3. A strict seven-port reduction of double-star to Petersen

Use the Sage labels and retain exterior vertices

```text
X={0,1,3},  E(X)={(0,1)}.
```

The ordered crossing edges are

```text
(0,14), (0,15), (1,2), (1,11), (2,3), (3,4), (3,18).
```

Delete X and put labelled ports 0,…,6 at the surviving endpoints. This
produces a connected pole Q with 27 internal vertices. Delete the following
ten internal edges:

```text
(4,5), (7,8), (9,24), (12,13), (15,29),
(17,18), (19,23), (20,28), (21,22), (26,27).                  (3)
```

They are a matching. Suppress their twenty degree-two endpoints. The
result Q' has seven internal vertices. Its complete ordered edge set,
using port labels -1,…,-7 in the code, is

```text
(-7,25), (-6,14), (-5,2), (-4,11), (-3,2), (-2,16), (-1,14),
(2,6), (6,10), (6,16), (10,11), (10,14), (11,25), (16,25).
```

After reattaching the unchanged three-vertex exterior, the resulting
ten-vertex cubic graph is isomorphic to Petersen. The
[certificate](snarks_double_star_minor.json) records the original graph,
unchanged exterior, all deleted edges, fourteen internally disjoint
subdivision paths, the reduced pole, and an explicit Petersen vertex map.
Thus this is an actual rooted topological minor, not an inference from
an unlabelled Petersen-minor theorem.

Full enumeration gives

```text
|S_Q|=192,       sum f_Q=1008,       nonzero entries of f_Q ∈ {4,6};
|S_Q'|=120,      sum f_Q'=120,       nonzero entries of f_Q' = 1;
S_Q' ⊊ S_Q.                                                   (4)
```

Among the 120 surviving words, 72 have old multiplicity four and 48 have
old multiplicity six. The 72 lost words all have old multiplicity six.
In particular, (1) holds, and the closed colouring obstruction is
preserved. Both the input and output are simple, bridgeless, cyclically
five-connected snarks.

The seven orbit-normalized Kempe defects `d_i=<f,w_i>` are

```text
d(Q)  = (-16,-8,0,-24,0,0,0),
d(Q') = (  0,-4,0, -4,0,0,0).                                (5)
```

An independent traversal of every physical colouring's two-colour
components gives

```text
N(Q)  = (96,48,0,144,0,0,0),
N(Q') = ( 0,24,0, 24,0,0,0),
N_i=-6d_i.
```

These are the counts of colourings with three mixed-endpoint boundary
Kempe paths at the unique-colour port, as defined in the inherited note.
The seven coordinates therefore change in a controlled, measured way;
they do not have to remain equal for an admissible descent.

As an independent verification of both full tensors, each pole was closed
against all 105 junction-and-two-wire kernels and the seven additional
P5 kernels. A separate perfect-matching/two-factor counter, preserving
parallel edges and treating loops as uncolourable, agrees with all
112 tensor contractions for each pole. This is 224 independent closure
counts, and the completed bank determines the full 91-coordinate tensor.

## 4. The seven defects detect a concrete realizability failure

Entrywise, (4) gives `f_Q ≥ 4 f_Q'`. It is tempting to interpret the
nonnegative remainder

```text
h = f_Q - 4 f_Q'
```

as an additional cubic-pole contribution. This interpretation is false:

```text
|supp h|=120,  sum h=528,  nonzero entries ∈ {2,6},
<h,w_i>=(-16,8,0,-8,0,0,0).                                 (6)
```

Every actual cubic seven-pole has every defect nonpositive. The positive
coordinate at port 1 rules out h as an actual pole tensor, or even as a
nonnegative linear combination of actual pole tensors.

More generally, if `f_Q=λ f_Q'+h` with `h` in that physical cone and
`λ≥0`, the port-1 inequality forces

```text
-8+4λ≤0, hence λ≤2.
```

Pointwise nonnegativity alone permits λ≤4. The seven extra coordinates
therefore provide a strictly stronger realizability constraint in this
specific cyclic-five example. This does not prove realizability at λ=2.
It also does not invalidate support domination: realizability of the
algebraic remainder is unnecessary for (2).

## 5. Why full-support equality cannot be the universal seven-port target

The finite universe here is explicit: all three-vertex subsets inducing
exactly one edge, with every crossing edge retained as an ordered port.
Automorphisms are used only to avoid duplicate computations; all orbit
sizes are recorded.

| Graph | Labelled triples | Automorphism orbits | Seven-pole support range | Triples with support 120 |
|---|---:|---:|---:|---:|
| Petersen | 60 | 1 | 120–120 | 60 |
| Loupekine 0 | 528 | 46 | 120–384 | 54 |
| Loupekine 1 | 528 | 47 | 120–378 | 54 |
| Double-star | 1080 | 25 | 192–468 | 0 |

For a cubic graph of girth at least five the labelled triple count is
`|E|(|V|-6)`, which the census checks. The three complete receipt files are
[Loupekine 0](snarks_Loupekine0_triples.json),
[Loupekine 1](snarks_Loupekine1_triples.json), and
[double-star](snarks_DoubleStar_triples.json).

Support cardinality is invariant under port permutation. Consequently,
none of the 1080 double-star poles in this universe can be replaced by a
Petersen seven-pole with exact support equality while preserving that
three-vertex exterior. This is a scoped impossibility result, not a claim
that all seven-cuts or all larger intermediate reductions fail.

There is also a stronger obstruction to choosing the interface
arbitrarily. All 7! port permutations of a Petersen seven-pole yield
exactly 630 distinct supports. Exhaustively compare every one with each
double-star pole in the census. No candidate support is contained in the
old support for these three exterior orbits:

| Representative exterior triple | Orbit size |
|---|---:|
| (0,1,5) | 40 |
| (0,5,10) | 20 |
| (1,2,5) | 40 |

Thus **100 of the 1080 labelled interfaces forbid even strict
support-dominated direct replacement by a Petersen seven-pole**.
No search for minor branch sets can overcome that algebraic obstruction
at the fixed interface. The successful triple (0,1,3) admits exactly two
candidate supports; section 3 realizes one as an actual rooted minor.
The [candidate audit](snarks_double_star_support_candidates.json) records
all 25 orbits. The hostile pole for (0,5,10) has all seven Kempe defects
equal to zero, so vanishing of these coordinates alone does not certify
that an interface admits the desired descent.

## 6. Exact support equality still gives a second useful cyclic-five move

Let N be the Petersen negator, obtained by deleting the path 0–1–2.
In port order `(I0,I1,O0,O1,r)`, its support consists exactly of words
where one two-port connector is equal-coloured, the other is
different-coloured, and the residual colour is the sum through the
different-coloured connector in `F2^2`. There are 36 such words and each
has one extension; the script checks all of them.

Form NN from two copies: join the output connector of the first to the
input connector of the second, and join their residuals to a new cubic
vertex with one new port. In any colouring, the joined connector must
have equal colours: otherwise both residual colours would coincide at
the new vertex. The two external connectors therefore have unequal
colours. Conversely every colouring of a three-vertex path P3 with its
five dangling edges extends to NN in exactly three ways, one for each
common colour on the joined connector. Thus

```text
f_NN = 3 f_P3.                                                (7)
```

In the script's labels NN has fifteen internal vertices. Delete

```text
(0,7), (4,5), (6,9), (8,106), (104,105), (107,109).
```

The resulting graph is a subdivision of the required five-port P3.
Its three branch vertices are 10,110,108, in that path order; the
five-port incidences after suppression are

```text
port0 at110,  ports1,2 at10,  ports3,4 at108.
```

The [rooted certificate](snarks_nn_minor.json) explicitly checks that the
minor has exactly the same ordered support. Replacing this NN block in
either Loupekine graph restores Petersen.

This equality lifts to every retained exterior boundary, including all
seven-port interfaces. For a concrete nonzero seven-port test, delete
vertices `(204,205,206)` in the third negator. Then

```text
f_old=3 f_new,  |supp f_old|=|supp f_new|=120,
d_old=(-12,0,0,0,-12,0,0),  d_new=(-4,0,0,0,-4,0,0).
```

All eighteen one-edge triples lying in the third negator are checked for
each Loupekine variant, with
[two](snarks_loupekine0_lift.json)
[receipts](snarks_loupekine1_lift.json).
Permuting the three negators explains the 54 minimal-support triples in
each whole-graph census.

## 7. A cut certificate gives an infinite cyclic-five nonflower family

The following sufficient criterion upgrades the finite larger examples
to an all-order result.

**Cyclic-five replacement lemma (PROVED).** Let G be a simple cubic graph
of girth at least five and cyclic edge-connectivity at least five, with at
least ten vertices. Replace an induced three-vertex path by a simple
cubic five-pole Q in the corresponding port order. Suppose:

1. Q contains a rooted subdivision of that five-port path.
2. Its five port-owning vertices are distinct and pairwise nonadjacent.
3. Every subset A of at least three internal vertices of Q satisfies
   `b_Q(A)≥5`, where `b_Q(A)` counts edges leaving A, including dangling
   edges whose internal endpoint lies in A.

Then the resulting cubic graph G* is simple and cyclically
five-connected.

**Proof.** In a cubic graph, a vertex set of size at least three and
boundary at most four necessarily contains a cycle: the degree sum gives
`e(A)=(3|A|−|δ(A)|)/2≥|A|`, using integrality when `|A|=3`.
Thus the hypotheses on G exclude a cut of size at most four with at
least three vertices on both sides.

Suppose G* has such a cut. If each side contains at least three exterior
vertices, restrict the cut to the rooted path subdivision and read the
side of each of its three branch vertices. Every edge of the restored
path whose endpoints are on different sides requires a crossing along
its corresponding subdivision path. Since these paths are internally
disjoint, the induced cut in G is no larger and has at least three
vertices on each side. This is impossible.

Otherwise choose the side containing j≤2 exterior vertices, call their
set X, and let A be its internal Q vertices. Its total size is at least
three. Each exterior vertex is adjacent to at most one of the five
old path ports: two neighbours on a three-vertex path would give a
triangle or four-cycle in G. Write e_X for the number of edges within X
and t_A for the number of attachments from X to A. Then `t_A≤j` and

```text
|δ_G*(A∪X)| = b_Q(A)+3j−2e_X−2t_A.                         (8)
```

If `|A|≥3`, condition 3 gives b≥5, and the correction term is nonnegative
for j=0,1,2, since `e_X≤1` and e_X=0 for j≤1.
If `|A|=2`, then j≥1. Nonadjacent A has b=6, so (8) is at least five.
Adjacent A has b=4 and, by condition 2, at most one port owner;
thus `t_A≤1`, again giving at least five.
Finally, if `|A|=1`, then j=2; that vertex owns at most one port, so
`b=3`, `t_A≤1`, and (8) is at least five.
The case A empty cannot have total size at least three.
All cases contradict the chosen cut. Distinct port owners and simplicity
of G and Q ensure simplicity after substitution. QED.

The Petersen NN pole satisfies all three conditions. Section 6 supplies
the rooted subdivision. Its five port owners are
`(0,5,6,104,108)` and form an independent set. An exhaustive check of
all 32,647 internal vertex subsets of size at least three gives the
following minima of b_Q(A), indexed by cardinalities 3,…,15:

```text
5,6,5,6,5,6,7,8,7,8,7,6,5.
```

This fixed finite gadget certificate is recorded in
[snarks_nn_cyclic_gate.json](snarks_nn_cyclic_gate.json).
The lemma applies to every host satisfying its hypotheses; it does not
extrapolate a finite host census.

Let F_k be the graph obtained from the flower J_k by replacing the
three-vertex path `(b0,a0,c0)` with the Petersen NN pole above. For every
odd k≥5 this construction adds twelve vertices, so

```text
|V(F_k)|=4k+12.
```

Equation (7) proves uncolourability for every such k. The explicit rooted
minor restores exactly the original labelled J_k; the inherited flower
reductions then give Petersen. These orders are divisible by eight,
whereas every odd flower snark has order four modulo eight, so the family
is not a renaming of the flower family.

The standard cyclic connectivity of flowers supplies the host hypothesis;
see the discussion preceding Proposition 2.5 in Chladný and Škoviera,
[*Factorisation of Snarks*, EJC 17 (2010), R32,
§2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v17i1r32/pdf/).
The smallest hosts J5 and J7 are also checked directly by the script.
Consequently **every F_k for odd k≥5 is cyclically five-connected**.
The cut surrounding NN has five edges and cycles on both sides, so its
cyclic edge-connectivity is exactly five.

The following are additional exhaustive finite checks of whole graphs:

| k | Vertices of F_k | Perfect matchings | Tait colourings | Four-edge deletion sets checked |
|---:|---:|---:|---:|---:|
| 5 | 32 | 178 | 0 | 194,580 |
| 7 | 40 | 658 | 0 | 487,635 |
| 9 | 48 | 2,482 | 0 | 1,028,790 |
| 11 | 56 | 9,490 | 0 | 1,929,501 |

The [family receipt](snarks_nn_flower_family.json) retains every edge set.
All tested members have an explicit five-cut around NN and no smaller
cyclic cut. Their reduction chains consist of genuine graph minors.
The first move has full boundary equality (7), and the subsequent flower
strip moves have the previously proved full support equality.

The all-order conclusion comes from the three-condition substitution
lemma and its fixed NN certificate. The table independently checks its
consequences on larger graphs. The uncolourability and rooted minor
chain hold separately by (7).

## 8. Pythagorean interpretation and what remains

At a retained interface let f be the pole tensor and g the exterior
tensor. The equal-norm vectors `U=f⊗g` and `V=g⊗f` give orthogonal
`U+V` and `U−V`, with squared lengths

```text
2(M+C²), 2(M−C²),  M=||f||²||g||²,  C=<f,g>.
```

For nonzero f and g, the right triangle is isosceles precisely when
C=0. The NN move scales f by three and therefore preserves the normalized
triangle for every exterior. The double-star move changes the triangle
for general exteriors; support containment preserves its isosceles
condition whenever the original closure is obstructed. This is the exact
surviving geometric connection.

The larger open step is to find an admissible rooted support-dominated
minor in every cyclically five-connected snark. Seven-port tomography
identifies the tensor, and the Kempe defects constrain realizability, but
neither supplies a structural minor by itself. Double-star demonstrates
that strict containment can be the correct target when equality is
provably unavailable.

## Reproduction

From the repository root:

```text
python 04-computation/black_descent_20261020c_snarks.py
python -O 04-computation/black_descent_20261020c_snarks.py
```

The script uses NetworkX and the previous committed pole enumerator.
Every check is explicit and survives `-O`. The matching counter used for
the 224 closure checks is independent of that pole recursion. The
double-star discovery used a fixed random search, but the recorded
certificate is verified directly; no random search is needed to reproduce
the result.
