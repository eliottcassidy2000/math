# Kempe-certified minor reductions beyond flower snarks

**PROVED:** the elementary rooted-pole minor lemmas, support transport,
the explicit dot-product reduction chain, and the connected seven-pole
construction below. **FINITE-EXACT:** the stated graph census and receipts.
**CITED:** the classical dot-product connectivity theorem. A universal
reduction of arbitrary snarks to Petersen remains **OPEN in this work**.

## 1. Inheritance, source audit, and the new target

The anchor is an actual minor operation preserving the **full labelled
boundary support**, rather than just preserving a zero colouring count.
The niche is the repaired seven-port closure

```text
J(s) = [s0,s1,s2 all distinct] [s3=s4] [s5=s6].
```

The wildcard is whether a connected cubic tree can realize that repair.
The answer for trees is decisively negative; cyclic poles succeed.

The closest proved mechanism is
[THM-4116, ordered gluing and four-port signatures](../../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md).
The immediately inherited result is the
[flower strip reduction](../boundary_support_20261020_snarks.md).
The [incoming addendum](../boundary_support_20261020_incoming.md), section 2,
corrects the seven-port equality mistake: an equal triple plus two equal
pairs violates cubic boundary parity, while a distinct triple plus those
pairs is an actual cubic junction. The underused sidecar here is a
**rooted path model** extracted from a colouring, together with checks on
the graph obtained after the ports are glued.

**Primary-source audit, 2026-10-20.** Chladny and Skoviera,
[*Factorisation of Snarks*, EJC 17 (2010), R32](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v17i1r32/pdf/),
already treats colour-open four-poles, their isochromatic/heterochromatic
classification (Proposition 3.5), and substitution by colour-equivalent
poles (Lemma 4.1 and Remark 4.2). Printed page 12 states preservation of
cyclic four-edge-connectivity and girth at least five under the dot product
with independent left edges. Example 5.5, printed pages 23–24, already
uses the flower-strip equivalence `N4 ~ N2`. These are antecedents, not new
claims of this session. The proofs below spell out rooted minor witnesses
and the necessary graph-type guard; no priority claim is made.

The concept board is: boundary support; Kempe pairing; rooted minor model;
seven-port cubic parity; and the graph class after closure. Every positive
result below retains the first three, while the two hostile tests show
why neither parity nor support alone handles the last coordinate.

## 2. Isochromatic four-poles contain the required two-wire minor

Represent ports by labelled degree-one vertices. All other vertices have
degree three. Isolated edges with two ports are permitted. Colours are
`0,1,2`, identified when useful with the nonzero elements of `F2^2`.
Let

```text
B = {(a,a,b,b): a,b in {0,1,2}},
A = {(a,b,c,d): a!=b, c!=d, {a,b}={c,d}}.
```

Here `|B|=9` and `|A|=12`. B is the support of two independent wires joining
ports 0–1 and 2–3. A is the support of two adjacent cubic vertices, with
ports 0,1 at the first and 2,3 at the second.

**Lemma 1 (PROVED).** Any colourable cubic four-pole with support exactly B
contains two vertex-disjoint paths joining ports 0–1 and 2–3. Deleting
unused edges and vertices and suppressing the internal vertices of those
paths gives the two-wire pole as a rooted minor. This minor has exactly
the original full labelled support.

**Proof.** Choose an extension of boundary word `0011`, which exists by
the support assumption. Retain only edges of colours 0 and 1. Every
internal vertex has degree two in this subgraph; its components are cycles
and two port-to-port paths. If one path joined the two different couples,
switching 0 and 1 on that path would produce a boundary word outside B.
That is forbidden. Hence the paths pair the ports within their prescribed
couples. The two path components are vertex-disjoint. Their union is the
required rooted subdivision of two wires. QED.

The actual witness from the Petersen eight-vertex four-pole `P-{0,1}` has
ordered port neighbors `(4,5,2,6)` and internal paths

```text
4-9-7-5,       2-3-8-6.
```

Both paths are independently authenticated from a complete colouring, not
inferred from the fact that Petersen has a minor. This four-pole's extension
vector is exactly `2 B`, where B now denotes its indicator vector.

The same procedure succeeds on the marked adjacent-vertex poles of both
Blanusa graphs and J5. The retained paths need not cover every vertex:
closed Kempe cycles and other unused material are simply deleted.

## 3. Heterochromatic four-poles contain the adjacent-junction minor

**Lemma 2 (PROVED).** Any colourable cubic four-pole with support exactly A
has the two-adjacent-vertex four-pole as a rooted minor, with the same
ordered endpoint couples and full labelled support.

**Proof.** Choose boundary word `0101` and inspect the 0/1 Kempe subgraph.
If a path connected ports 0–2, switching it would create `1111`, outside
A. If a path connected ports 0–3, switching it would create `1100`,
also outside A. Thus the two Kempe paths necessarily connect 0–1 and 2–3.
They are vertex-disjoint.

All four ports belong to one connected component of the original pole.
Indeed, a colourable component cannot have exactly one port, by parity.
If the four ports split into two two-port components, each pair would
always have equal colours; after independent colour relabelling both
components would admit the all-equal word, contrary to A. Components with
no ports do not affect this argument.

Therefore some internal path connects the two Kempe paths. Choose one
shortest between their vertex sets. Its interior meets neither path.
The union is a tree with four labelled leaves and two cubic branch
vertices, one separating ports 0,1 and the other separating ports 2,3.
Delete all remaining material and suppress its degree-two vertices.
The result is exactly the adjacent-junction pole, whose support is A.
QED.

The script authenticates this construction for **all 75 pairs of
independent edges in Petersen**: disconnecting either pair produces a
heterochromatic four-pole, and the Kempe paths plus connector produce the
claimed rooted minor in every case.

## 4. The general four-cut consequence, with its type guard

**Proposition 3 (PROVED).** Suppose an uncolourable cubic graph is split
by a four-edge cut into two colourable poles. After ordering the ports,
one pole has support B and the other A. Each can therefore be replaced by
its corresponding rooted minor from Lemma 1 or Lemma 2 without changing
the labelled extension predicate in any retained exterior context.

**Proof.** Four-port cubic parity allows only the four colour orbits with
representatives `0000,0011,0101,0110`. Every colourable four-pole admits at
least two orbits: a Kempe switch on a port-to-port path changes the orbit,
starting either with an all-equal word or with a two-colour word. Since
the two poles glue to an uncolourable graph, their supports are disjoint.
They must therefore occupy two orbits each. Exactly one contains the
all-equal orbit. Relabel the ports to make its other orbit `0011`; this
is B, and the complementary pair of orbits is A. Apply the two lemmas.
Finally, any exterior colouring induces a fixed port word, whose
extendability is unchanged. QED.

This explains the established dot-product factorization in operational
terms: Kempe paths supply the structural minor certificate in addition to
the Boolean substitution rule.

**The graph-type guard is essential.** A rooted pole minor is not
automatically a useful reduction inside the class of simple bridgeless
snarks. Glue the Petersen eight-vertex four-pole to the two adjacent vertices that
complete it to Petersen. Compressing its B-pole to two wires now creates
a loop at each of those two vertices, joined by their remaining edge:
the looped dumbbell. The colour obstruction survives, but simplicity and
bridgelessness do not. Deleting the loops to take an ordinary simple minor
would lose the cubic obstruction.

Consequently a global step toward Petersen must additionally verify that
the completed factor remains in the required graph class. A cycle-separating
minimum four-edge cut in a cyclically four-connected cubic graph has
distinct endpoints on each side, preventing this particular loop failure;
the remaining graph properties must still be checked or established by a
factorization theorem. We do not infer a universal Petersen reduction
from Proposition 3 alone.

For the explicit family below, the type guard is automatic: every reverse
step restores a previously specified simple cyclically four-connected
graph, with exactly its original labelled edge set.

## 5. An infinite nonflower family with certified Petersen reductions

Start with `G_0=Petersen`. Given `G_r`, select two independent edges, sever
them into four endpoints, and insert a fresh Petersen eight-vertex four-pole, with
each endpoint pair attached to its corresponding equal-colour couple.
Call the result `G_(r+1)`.

The extension vector of the inserted pole is exactly `2B`. Thus, for an
arbitrary retained labelled boundary as well as for a closed graph,

```text
f_(G_(r+1)) = 2 f_(G_r),                                  (1)
```

where the boundary lies outside the inserted pole. In particular all
members are uncolourable. Lemma 1 gives a rooted minor reduction
`G_(r+1) -> G_r`, so reversing the insertions gives an explicit
obstruction-preserving chain to Petersen.

Each insertion adds eight vertices, hence

```text
|V(G_r)| = 10+8r.
```

For r>0 these orders differ from every flower-snark order, which is
divisible by four. The construction is the classical iterated Petersen
dot product. The cited dot-product connectivity statement guarantees
cyclic four-edge-connectivity and girth at least five for the entire
family when the two chosen edges are independent. It is used only for
that all-order graph-class assertion; the exact support identity and
minor model have been proved above.

The deterministic implementation picks the lexicographically first pair
of independent edges at each step. This fixes actual graphs, not just
their orders. It checks the following members independently:

| r | vertices | perfect matchings | labelled Tait colourings |
|---:|---:|---:|---:|
| 1 | 18 | 19 | 0 |
| 2 | 26 | 57 | 0 |
| 3 | 34 | 163 | 0 |
| 4 | 42 | 458 | 0 |
| 5 | 50 | 1288 | 0 |
| 6 | 58 | 3655 | 0 |
| 7 | 66 | 10108 | 0 |
| 8 | 74 | 28425 | 0 |

For each graph the code verifies cubicity, connectedness, absence of
bridges, girth five, absence of **every** cycle-separating cut of size at
most three, and exact restoration of its predecessor by the rooted
minor. The cut search is exhaustive: any minimal cut of size at most
three appears as a bridge after deleting at most two other edges.

Colour counts in the table are independently obtained by enumerating a
colour-0 perfect matching and checking whether every complementary
two-factor cycle is even. This is independent of the pole colouring
recursion and support proof. As a positive control, the same insertions
starting with K4 give exactly 12, 24 and 48 colourings after one, two and
three steps, as (1) predicts.

## 6. A connected seven-pole realizes the repaired cubic junction exactly

Begin with the disconnected seven-pole consisting of one cubic vertex
with ports 0,1,2 and the two independent wires 3–4 and 5–6. Its vector is
the indicator J from section 1.

First replace the edge leading to port 0 and wire 3–4 by a Petersen
eight-vertex four-pole, preserving their endpoint pairings. This joins the first two
components and multiplies the boundary vector by two. Then do the same
with the edge leading to port 1 and the remaining wire 5–6. The resulting
pole Q is connected and has

```text
17 internal cubic vertices, 7 labelled ports,
f_Q = 4J,    |supp(f_Q)| = 54,    sum_s f_Q(s) = 216.        (2)
```

Both inserted poles have the rooted two-wire minors from Lemma 1.
Undoing the insertions gives an actual rooted minor from Q to the cubic
junction plus two wires. Formula (2) is checked at all `3^7` labelled
boundary words by a direct exhaustive edge-colouring engine.

This is a realizable connected patch for the parity repair. It is more
than an abstract kernel: its complete edge set, port order, colour table,
and minor contractions are generated and checked by the script. It does
not by itself imply that a general seven-pole can be replaced by Q.

The relation remains valid inside any larger pole or cubic graph, with
the same graph-type caveat from section 4. Thus Q can be inserted into a
Petersen-derived construction without changing the colouring obstruction,
provided the final graph passes the stated loop, simplicity and
connectivity checks.

## 7. Cheap hostile probes explain why a tree is different

A connected cubic tree with seven ports has five internal vertices.
There are two possible unlabelled internal-tree shapes with maximum
degree three. Every such pole has exactly 96 extensions and 96 distinct
boundary words: choose a colouring at one vertex in six ways, then two
choices at each of the remaining four vertices. Internal edge colours
are uniquely determined by the boundary through `F2^2` cut sums.

In particular its support cannot equal the 54-word support of J.

There is a stronger obstruction to the tempting minor operation.
Delete two disjoint internal tree edges and suppress the resulting
degree-two vertices. The output is one cubic junction plus two wires,
but its support is **disjoint** from the original tree support.
For any deleted edge, the original colour is the XOR sum of boundary
colours on one side and must be nonzero. After the deletion, that side
is a union of complete components of the reduced pole, so the same XOR
sum must be zero. No boundary word extends to both.

The script checks both tree shapes and all five allowable choices of two
disjoint internal edges. Every original support has size 96; every
reduced support has size 54; every intersection is empty. This is a
concrete failure of obstruction-preserving reduction despite correct
seven-port parity.

Another cheap probe removes three vertices inducing exactly one edge
from Petersen. All 60 choices give seven-poles whose support and total
mass are both 120. None is equivalent to any 54-word repaired junction
kernel. The connected seventeen-vertex construction works because its
cycles enforce the required equal-colour couples; vertex count and
boundary parity alone do not.

## 8. Pythagorean predicate and what remains

Ordered support equality preserves zero overlap with every exterior
extension vector. If f,g are nonzero and

```text
U=f tensor g,  V=g tensor f,
M=||f||^2||g||^2,  C=<f,g>,
```

then `U+V` and `U-V` are perpendicular and their squared lengths are
`2(M+C^2)` and `2(M-C^2)`. The isosceles endpoint is exactly `C=0`.
Thus every admissible reduction above retains the same Pythagorean
predicate used in the prior session. General support equivalence need
not preserve nonzero counts or angles. The special Petersen insertion
has the stronger scalar identity (1): at a fixed retained interface,
replacing f by 2f scales U,V and every triangle length by two. Hence its
normalized Pythagorean point and angle are preserved exactly, even for
colourable exteriors. The connected seven-pole replacement has the same
property with scale factor four. This refers to the marked interface
geometry, not to a separately chosen maximum-cut invariant of the whole
graph.

The remaining structural task is narrower now: find applicable,
graph-class-preserving rooted replacements in cyclically
five-connected cores, or in cuts where one pole is already uncolourable.
Neither four-port factorization nor the connected J gadget proves that
such a replacement always exists. The parent session's extra seven-port
Kempe coordinates are relevant precisely because J responses alone do
not determine a general seven-pole.

### Independent check of the seven missing Kempe coordinates

The [seven-port reconstruction proof](junctions.md) defines `w_i` to be
+1 when port i is a singleton colour in pattern `(5,1,1)`, -1 when it is
the singleton in `(3,3,1)`, and zero otherwise, using one coordinate per
global-colour orbit. It proves `6<f,w_i>=-N_i`, where N_i counts actual
colourings having a unique colour at i and three other-colour Kempe paths
whose endpoints all have different colours.

The independent engine in this packet enumerates every full colouring,
constructs the actual two-colour components, and obtains:

| Pole | `<f,w_i>` for i=0,...,6 | `N_i` |
|---|---|---|
| Connected 17-vertex Q | `(0,0,0,0,0,0,0)` | `(0,0,0,0,0,0,0)` |
| Five-vertex path tree | `(0,0,0,-4,0,0,0)` | `(0,0,0,24,0,0,0)` |
| Petersen minus `{0,1,3}` | `(0,-4,0,-4,0,0,0)` | `(0,24,0,24,0,0,0)` |

The path's port incidences are `(0,0,1,2,3,4,4)`. For the Petersen deletion,
the live endpoint of each ordered cut edge is `(4,5,2,6,2,4,8)`. These
receipts verify the identity with genuine multiplicities and physical
paths; no rank calculation is used in this check. The switch-orbit proof,
its factor of six, the tree pairings, and the modular-to-rational rank
argument in the companion were also read and audited.

## 9. Reproduction and exact universe

Run

```text
python3 -B 04-computation/blue_quotient_20261020b_snarks.py
python3 -B -O 04-computation/blue_quotient_20261020b_snarks.py
```

The [script](../../../04-computation/blue_quotient_20261020b_snarks.py)
and [receipt](../../../04-computation/blue_quotient_20261020b_snarks.out)
use Python and NetworkX 3.5. The universe is: the four marked isochromatic
poles; all 75 independent-edge openings of Petersen; eight deterministic
dot-product graphs through 74 vertices; three colourable insertion
controls; the connected seven-pole at every boundary word; both
five-vertex tree shapes and their five smoothing choices; and all 60
specified Petersen vertex triples. No unlabelled-snark census or
random-sampling completeness is claimed. Three additional full-colouring
enumerations independently check the singleton/Kempe-path identity above.
