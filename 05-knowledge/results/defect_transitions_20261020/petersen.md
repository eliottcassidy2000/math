# Rooted Kempe reductions and the cubic-vertex repair

**Status: PROVED elementary reductions + FINITE-EXACT examples.** This packet
extends the constructive Petersen reduction beyond the odd flower ring, to
iterated dot products and explicit vertex substitutions. It does not show that
every snark has a reduction of one of these types. The boundary support and
the graph-minor certificate are both retained.

Companion: `petersen_check.py`, with deterministic transcript
`petersen_check.out`. Run from the repository root:

```text
python 05-knowledge/results/defect_transitions_20261020/petersen_check.py
python -O 05-knowledge/results/defect_transitions_20261020/petersen_check.py
```

## 1. Inheritance, precise objects, and prior art

The closest proved mechanism is **THM-4116**,
`01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md`,
especially its adjacent-edge parity and Kempe compression. It already proves
the full uniform four-port profile used below. The new local use here is to
keep the actual Kempe paths as a rooted minor certificate.

The flower-family support/minor reduction and the hostile bridge example are
in `05-knowledge/results/boundary_support_20261020_snarks.md`. The latter
shows that a four-pole's Petersen-like support does **not** force a Petersen
minor. Our target minor for that pole is only two wires. The corrected near
miss is the seven-port equality merger in
`boundary_support_20261020_incoming.md`: three meeting edges require
all-distinct colours, not a common colour. The least-used sidecar is the
actual coloured path routing, discarded by a support table.

Dot products and these support types are established mathematics.
Chladný and Škoviera, *Factorisation of Snarks* (2010), §3, define
isochromatic and heterochromatic four-poles and prove the corresponding
colour-equivalence statements; see in particular Propositions 3.1–3.6.
Their Theorem 3.7 supplies a dot-product decomposition at a four-edge cut in
a critical snark. The paper's global factorisation theorem has
cyclically five-connected irreducible factors, not just Petersen factors.
[Primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v17i1r32/pdf/).
The repository already considered Petersen dot products in
`05-knowledge/results/collatz_label_map_paley_snarks_20261001.md` §3.
No new-family or literature-priority claim is made here.

A pole has cubic interior vertices and ordered, individually distinguished
boundary semiedges. Distinct semiedges may meet the same interior vertex.
Three labelled colours mean the three nonzero elements of `F_2^2`, not
`Z/3Z`. Its extension vector `f(s)` counts all interior colourings with
ordered boundary word `s`. Gluing to an exterior vector `q` gives
`sum_s f(s)q(s)`. In graph-minor statements, temporarily append four or
three distinct terminal leaves and keep their labels fixed.

The live concept board is: boundary support; Kempe paths; rooted minors;
cubic-vertex repair; dot-product factorisation. MISTAKE-507 forbids replacing
the exact colouring predicate by prime analogies or the wrong colour group;
MISTAKE-643 forbids transferring geometric line-incidence claims without the
correct map. Neither shortcut is used.

## 2. The support theorem and the actual rooted two-wire minor

Put

```text
B = {(a,a,b,b) : a,b in {0,1,2}},
W = the two boundary wires 0--1 and 2--3.
```

**Lemma 1 (inherited support, explicit proof).** Let `H` be an uncolourable
finite simple cubic graph, delete the endpoints `u,v` of an edge, and retain
their other four edges as ports, ordered by deleted endpoint. Its profile is

```text
f(s) = m 1_B(s),  m >= 0.
```

When nonempty, `m` is even and `supp(f)=B`.

For each colour `c`, the number `k_c` of boundary occurrences satisfies
`|V(H)|-2 = 2m_c+k_c`, so all `k_c` are even. Up to a global colour
permutation the only patterns are `0000,0011,0101,0110`. Either crossed
pattern would restore `u,v`: each endpoint receives two different colours,
both omit the same third colour, and that third colour can colour `uv`.
Uncolourability excludes both crossed patterns.

In an extension of `0000`, take the two-colour `0/1` subgraph including
terminal leaves. Every interior vertex has degree two in it. It consists of
cycles and two boundary paths. A path between the different port pairs would,
after switching its two colours, create a forbidden crossed word. Therefore
the paths pair `0--1` and `2--3`. Switch the latter path to obtain `0011`.
The same no-cross argument starting at `0011` gives the inverse: the
underlying two-colour subgraph is unchanged by switching. Thus the numbers
of extensions of these two words are equal. Global colour permutations give
the asserted constant profile. In particular, existence of either type
forces both types and all nine labelled words.

There are `3m` extensions with constant boundary word. The global `S_3`
action on that set is free because there is an interior cubic vertex and all
three colours occur there. Hence `6 | 3m` and `m` is even. Nonemptiness is
essential; it is not supplied merely by `H` being a snark.

**Lemma 2 (rooted linkage).** Every realizable cubic four-pole with support
exactly `B` has `W` as a rooted minor, and replacing it by `W` preserves
colourability in every exterior with the same port order.

Choose any extension of `0011`. Its `0/1` paths cannot cross the port
pairing, by the switch argument above. They give two vertex-disjoint paths
between the designated terminal pairs. Delete all other edges and vertices,
then contract path interiors while preserving their terminal endpoints.
These are legal graph-minor operations and yield `W`. Its extension profile
is `1_B`, so both old and new gluings are colourable exactly when the
exterior supports some word in `B`. If the old profile is `m1_B`, its
colouring count is exactly `m` times that of the reduced closure.

This also works for a nonempty support contained in `B`: the preceding
Kempe switches and palette permutations force equality with `B`. There is
no need to infer global structure from nine numbers; a single colouring
extracts the two required routes.

**Minimal hostile to the converse.** Take two cubic vertices joined by one
edge, with two boundary ports at each vertex. Deleting their connecting edge
gives the rooted wire minor `W`. But their support has the twelve crossed
words `a!=b, c!=d, {a,b}={c,d}`, disjoint from `B`. A rooted wire minor alone
does not certify support equivalence. Also, changing `W`'s pairing changes
its support: `0011` and `0101` distinguish the original and crossed pairings.

## 3. A concrete reduction beyond flower snarks

Use the following fixed Petersen graph:

```text
01 04 05 12 16 23 27 34 38 49 57 58 68 69 79.
```

Delete vertices `0,5`. Its ordered ports attach at `(1,4,7,8)`, paired
`(1,4)|(7,8)`, and its extension profile is `2 1_B`. The two wire paths are

```text
1--2--3--4,       7--9--6--8.
```

They cover all eight remaining vertices. Delete the other four pole edges
`16,27,38,49` and suppress the two paths. This is an explicit rooted
minor certificate; no minor solver is involved.

Given distinct edges `ab,cd` of a finite simple cubic graph `G`, delete them,
insert a fresh copy of this eight-vertex pole, and join its four ports to
`a,b,c,d` in that order. Call the result `I(G;ab,cd)`. Then

```text
#Tait(I(G;ab,cd)) = 2 #Tait(G),
G is a minor of I(G;ab,cd).
```

More strongly, the displayed paths make a **spanning subdivision of `G`**
inside the new graph. The four unused pole edges are the only extra edges.
Consequently a connected bridgeless `G` gives a connected bridgeless result:
every edge of the subdivision remains on a cycle, and any additional edge
has a path between its endpoints in that connected subdivision.

If the selected edges are independent and `girth(G)>=5`, the result has
girth at least five. A cycle inside the pole already has length at least
five. A cycle meeting both sides uses either two or four attachment edges.
In the two-attachment case, the distinct pole ports are at distance at
least two, and the distinct exterior endpoints need at least one exterior
edge, giving length at least five. In the four-attachment case there are
already four attachment edges and at least two nonempty pole paths.
An exterior-only cycle was already a cycle of `G`. Independence matters:
sharing an exterior endpoint can introduce a four-cycle.

For independent edges this is a Petersen right-factor dot product. Starting
with Petersen and iterating it gives graphs of order `10+8r`, all
uncolourable, connected, bridgeless, and of girth at least five. A strict
definition of “snark” sometimes also requires cyclic edge-connectivity at
least four; that additional property is not asserted for arbitrary
substitution schemes here. The two named first-step graphs satisfy it,
as independently checked in the finite audit.

**Certificate composition.** Maintain fifteen paths labelled by the fifteen
edges of the original Petersen graph. When a source edge is replaced, replace
its occurrence in any certificate path by the corresponding fresh wire
path. New interiors are disjoint from the old graph and from each other.
The fifteen paths remain a Petersen subdivision after every insertion,
including an insertion that cuts an edge currently unused by the model.
This proves all-depth transport; the depth-100 run is a check of the
implementation, not the proof of an infinite claim.

**Other donors.** Lemma 1 and Lemma 2 apply identically to any donor
`H-{u,v}` with `m>0`; its local Kempe paths replace the two source edges.
The count multiplier is `m`, and every previously certified Petersen
subdivision transports through those paths. The checker gives explicit
paths for Blanuša I (`m=4`), Blanuša II (`m=6`), and two cuts of `J_5`
(`m=12,4`). These allow arbitrary mixtures of certified donors.
Unvisited donor vertices may be deleted in the minor certificate.
Connectivity, girth, and cyclic connectivity for a general donor require
their own checks; support equivalence alone does not supply them.

## 4. The cubic-vertex repair as an actual seven-port gadget

**Lemma 3 (three-port substitution).** Every nonempty colourable cubic
three-pole has precisely the six all-distinct boundary words, each with
the same positive multiplicity `r`. Its three boundary ports lie in one
connected component, which has the rooted cubic claw as a minor.

Boundary parity makes all three colour counts odd, hence each is one.
Palette permutations give every all-distinct word with the same count.
A connected component with a single boundary port would have nonzero
`F_2^2` boundary sum and is impossible in a proper colouring. A distribution
of three ports among several components would necessarily contain such a
component. Thus all three ports belong to one component. Delete any closed
components, retain a minimal tree joining the three terminal leaves, and
suppress its degree-two vertices. This gives the claw, with labels fixed.
Equivalently, contracting the interior component gives the required cubic
vertex. The tree argument also gives a topological claw routing, sufficient
to transport subcubic subdivision certificates.

Replacing a cubic vertex by this pole therefore preserves colourability
for every exterior, multiplies the count by `r`, and has the original graph
as a minor. A triangle with one port at each vertex has `r=1`.
`K_{2,3}` with ports at its three degree-two vertices has `r=2`.
The Petersen graph with one vertex deleted has empty support and is a
hostile control: the nonempty hypothesis cannot be dropped.

Now form a seven-pole as the disjoint union of a triangle three-pole and
the Petersen four-pole above. Keep the ordered boundary

```text
(triangle port 0,1,2; Petersen port 0,1,2,3).
```

Its exact profile is

```text
2 J(s), where
J(s) = [s0,s1,s2 all distinct] [s3=s4] [s5=s6].
```

It has `6*9=54` supported words, each with multiplicity two. All three
colour counts are odd, as required at a seven-edge cubic boundary.
Contract the triangle to the cubic vertex and use the two Petersen
wire paths. This reduces to **one cubic vertex plus two wires**, with all
seven port labels retained. For every exterior `Q`,

```text
#Tait(Q glued to this seven-pole)
    = 2 #Tait(Q glued to the cubic-vertex-and-two-wire repair).
```

The supported boundary words agree exactly, and the reduced graph is an
actual minor. Any pre-existing Petersen model in the repaired closure
transports through the local claw and two-wire routes. Replacing the
triangle by any colourable three-pole gives multiplier `2r`.

This completes a concrete repair construction. It does **not** establish
that every arbitrary seven-pole whose support happens to be `J` contains
this rooted forest minor. Nor does it prove that every unresolved cubic
graph exposes such a seven-port decomposition. Those are separate
structural questions.

### A forced partial linkage in every realizable J-support seven-pole

There is a useful general statement even without the product construction.
Let an arbitrary cubic seven-pole have support exactly `J`, and choose an
extension of the supported word

```text
(0,1,2;0,0;1,1).
```

Its `0/1` subgraph has six terminal endpoints. Its three boundary paths
**must** pair `(0,1)`, `(3,4)`, and `(5,6)`. Switching a path that touches
only one endpoint of one of these groups would break either the
all-distinct condition at the first triple or an equality condition at a
wire pair. That would produce a colouring whose boundary word is outside
`J`, a contradiction. Thus this pole always contains two disjoint desired
wire paths and a disjoint `0--1` path.

The other two colour pairs give, separately,

```text
0/2 paths: (0,2), (3,4);
1/2 paths: (1,2), (5,6).
```

The same forbidden-switch argument proves each pairing. Paths within any
one line are vertex-disjoint; paths from different lines need not be.
In particular these statements cannot be concatenated without checking
intersections.

The exact remaining forest obligation is to join terminal 2 to the
`0--1` path while avoiding the two chosen wire paths. If such a joining
path exists, truncate it at its first intersection with `0--1`; their
union is a three-terminal tree and yields the rooted claw, disjoint from
both wires. This is a sufficient explicit additional certificate.
The support argument above does not yet supply that avoiding path or
prove that a different choice of the two wire routes always allows it.
The product construction supplies it because the third terminal is in
the separate three-pole component.

## 5. Exact audit and the next obstruction

The finite universe and controls are explicit:

* Six adjacent-vertex deletion poles: Petersen, the two authenticated
  Blanuša edge banks, two `J_5` edges, and a `J_3` edge with empty support.
  The two Blanuša banks and their naming are inherited from
  `04-computation/experiments/snark_boundary_transfer_20261018a.py`;
  they are copied as literal edge data, not imported with its algorithms.
* All 75 unordered independent-edge pairs in the labelled Petersen graph,
  using increasing endpoint order. Exactly 15 outputs are isomorphic to
  the inherited Blanuša I bank and 60 to Blanuša II. Every graph has zero
  Tait colourings, a verified Petersen subdivision, and girth at least five.
  The two banks have no cyclic cut of size at most three. This is an
  identification with those explicit named banks; no census completeness
  assertion for all 18-vertex snarks is needed.
* All three independent-edge insertions in `K_4`: six Tait colourings become
  twelve. Pole profiles use edge-assignment backtracking; whole-graph
  colouring counts use an independent perfect-matching/even-cycle method.
* The 54-state seven-port repair and its actual closures over `K_4` and
  Petersen. The former has twelve colourings; the latter has none.
  Triangle inflation introduces triangles, so these closures are not
  claimed to satisfy the stricter girth convention.
* A deterministic depth-100 insertion sequence has 810 vertices and 1215
  edges. All hundred successive Petersen subdivision certificates pass.
* Empty support, changed boundary pairing, and the two-vertex hostile
  separate the failed implications from the surviving theorem.

The strongest next target is now precise: find a certified local
isochromatic four-pole or colourable three-pole in a graph outside the
current recursive grammar, or exhibit a smallest graph with neither useful
site. At cyclically five-connected irreducible factors, the four-edge
decomposition is unavailable. The factorisation theorem does not say these
factors are Petersen. That is the structural obstacle left after the
boundary algebra and certificate transport have been handled.

No implication for the Collatz conjecture follows from these graph
reductions. A transfer there would still need an explicit arithmetic
encoding that preserves the fixed integer and its unresolved predicate.
