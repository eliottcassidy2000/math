# Flower strips: obstruction-preserving reductions to Petersen

**Status: PROVED elementary statements below; FINITE-EXACT the specified
computations.** This is a research result packet, not a canon promotion or
a proof of the general cubic Petersen-minor theorem. The larger flower
family has a complete constructive answer: a local replacement preserves
the full Boolean boundary relation and produces a smaller graph minor.

## Inheritance and the actual target

The closest proved mechanism is
[THM-4116, ordered boundary gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md):
two nonnegative extension vectors glue by an inner product. Its adjacent
edge signature `(m,m,0,0)` is the inherited four-port mechanism. The
canonical hostile is the Petersen four-pole; K4 is the positive control.
The corrected near miss is MISTAKE-507 in
[MISTAKES.md](../../01-canon/MISTAKES.md): cubic Tait colourings use the
three nonzero elements of `F2^2`, not a `Z/3Z` flow.

The preceding [graph/Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md)
identified zero boundary overlap with an exact real isosceles right
triangle, but did not give an obstruction-preserving minor reduction.
The incoming [snark boundary transfer packet](snark_boundary_transfer_20261018a.md)
already supplies authenticated Petersen subdivisions in both Blanusa
graphs, J5 and J7, and the colourable Desargues graph. Its warning that a
minor alone does not determine an extension relation is retained here.
The [source-context packet](source_orientation_boundary_contexts_20261018k.md)
provides a second explicit warning about losing boundary relations under
subdivision. These packets were read before the present derivation.

The live concept board is:

| Lane | Object and operation | Retained predicate | Essential sidecar |
|---|---|---|---|
| Anchor | Flower strips; delete two centres and suppress strand vertices | Extendability in every exterior | Six ordered ports and their strand identities |
| Niche | Adjacent-edge four-poles; three smoothing closures | Exact change in Tait count | Four-port colour orbits and multiplicity |
| Wildcard | Pythagorean polarization of boundary vectors | Zero overlap, hence isosceles endpoint | Both nonzero vectors; lengths need not be preserved |
| Hostile | Two subdivided K4s joined by a bridge | Same four-pole support as Petersen | Bridgelessness and a rooted minor model |

The successful move is to enlarge the local replacement from two modules
to four. Two modules cannot be replaced by three bare wires. Four modules
can be replaced by two, uniformly in the exterior.

## 1. The 27-state transfer is completely explicit

A **flower module** is a claw with centre `a` and leaves `b,c,d`. Each leaf
has an incoming and an outgoing strand edge. Order the three strands
`(b,c,d)`. For incoming word `x` and outgoing word `y` in `{0,1,2}^3`, let
`T[x,y]` count proper colourings of the three spokes. Then

```text
T[x,y] = 1
  iff x_i != y_i for each i
  and the three omitted colours {0,1,2}\{x_i,y_i} are all distinct;
T[x,y] = 0 otherwise.                                      (1)
```

The spoke colours are forced by `x,y`, so all entries are zero or one.
Concatenating `r` modules gives exactly `T^r`. There is no quotient by
colour permutations here.

**Proposition 1 (PROVED).** As the adjacency matrix of a graph on the 27
words,

```text
T = 2 K3 disjoint_union 3 C6 disjoint_union 3 isolated vertices. (2)
```

**Proof.** The three constant words have no extension: all three spokes
would have to avoid the same colour. On the six all-distinct words, a
transition cyclically permutes the colours. The even and odd permutation
classes are therefore two triangles.

For the remaining words write `q` for the singleton colour, `r` for the
repeated colour, and `i` for the singleton position. A transition preserves
`q`, changes `r` to the third colour, and changes `i` to either of the other
two positions. For each of the three choices of `q` this is
`K_(3,3)` minus a perfect matching, hence `C6`. These exhaust 27 words.
Equivalently, identifying colours with the nonzero elements of `F2^2`,
`q` is the nonzero boundary charge `x_b+x_c+x_d`. QED.

Squaring separates each six-cycle into its two bipartition classes. Thus
there are eight active blocks of size three, and on every block

```text
T^2 = I3 + J3,             T^4 = I3 + 5 J3.                 (3)
```

The three constant words remain isolated. More generally, on each active
block and for `h>=1`,

```text
T^(2h) = I3 + ((4^h-1)/3) J3.                              (4)
```

On each permitted odd-step block, with rows and columns indexed by the
singleton position or permutation coordinate,

```text
T^(2h+1) = ((2^(2h+1)+1)/3) J3 - I3.                      (5)
```

For a triangle the odd block maps the block to itself; for a six-cycle it
maps one bipartition class to the other. Formula (5) also holds for `h=0`.

## 2. A six-port replacement that preserves every exterior

**Theorem 2 (PROVED).** Replace four consecutive flower modules by two,
retaining the three incoming and three outgoing ports in the same order.
The two poles have identical Boolean boundary relations. The smaller pole
is a graph minor of the larger pole. Consequently, gluing either pole to
the same exterior yields a three-edge-colourable graph if and only if
gluing the other does.

**Proof.** By (3), `supp(T^4)=supp(T^2)`: both have 72 nonzero ordered
boundary pairs. In a four-module strip numbered `0,1,2,3`, delete the
centres of modules 1 and 2. Their six leaves now have degree two. Suppress
those six leaves along the three strands. The surviving modules 0 and 3
are joined directly, and the six outer ports are unchanged. Deletion and
degree-two suppression are minor operations.

For any fixed exterior colouring the six port colours extend across one
pole exactly when they extend across the other. Taking the union over all
exterior colourings proves the equivalence. Parallel edges can be retained
as distinct edges when the exterior creates them; the flower reductions
below stay simple. QED.

In fact `supp(T^(r+2))=supp(T^r)` for every `r>=2`; the four-to-two move is
the smallest such uniform deletion. Multiplicities are deliberately not
claimed equal: diagonal entries change from 6 to 2, off-diagonal entries
from 5 to 1.

**Minimal failed replacement.** Replacing two modules by three wires would
replace `T^2` by the identity. The constant boundary pair `(000,000)` is
allowed by the wires and forbidden by the modules. Conversely, `T^2`
allows unequal words in an active three-state block, whereas three wires
require equality. This records both failed inclusion directions.

## 3. All odd flowers reduce to Petersen with the obstruction retained

Define `J_k` on vertices `a_i,b_i,c_i,d_i`, `i modulo k`, with each `a_i`
joined to its three leaves. The b-strand is a k-cycle. The other two
strands join into the 2k-cycle

```text
c_0,c_1,...,c_(k-1),d_0,d_1,...,d_(k-1),c_0.
```

This is the standard odd flower family for `k>=5`. The same graph
definition is simple and cubic for every integer `k>=3`, including the
even values; J3 has a triangle and is not being called a girth-five snark.

Let `R` swap the c and d positions in a boundary word. Then

```text
C(J_k) = sum_x (T^k)[x,Rx].                                (6)
```

**Proposition 3 (PROVED).**

```text
C(J_k) = 0                       if k is odd,
C(J_k) = 3*2^(k+1)               if k is even.             (7)
```

For the untwisted closure of three separate k-cycles, the counts are
`2^(k+1)-4` for odd k and `2^(k+3)+16` for even k.

**Proof.** On all-distinct words, R interchanges the two K3 components, so
the twisted closure has no contribution. On each C6, R preserves the
repeated-colour class, while an odd number of transitions reverses it.
This proves the odd zero. At even k put `alpha=(2^k-1)/3`. In each of the
two three-state classes for a fixed charge, R fixes one position and
interchanges two. Formula (4) contributes
`(1+alpha)+2alpha=2^k` per class. There are six classes, giving
`6*2^k=3*2^(k+1)`. Taking ordinary traces in (2), or using (4)-(5), gives
the untwisted counts. QED.

The same calculation classifies every permutation closing the three
strands. A three-cycle of strand positions gives `2^(k+1)+2` colourings
for odd k and `2^(k+3)-8` for even k. Indeed it acts without fixed points
on each active three-state block, whereas a transposition has one fixed
point. Thus, for `k>=3`, a permutation-closed strip is uncolourable
**exactly when k is odd and the closing strand permutation is odd**.
The three conjugacy classes have different nonzero counts, so the actual
permutation remains meaningful after retaining its parity.

**All bijective strand wirings around a claw ring (PROVED).** More
generally, allow an arbitrary permutation `sigma_i in S3` to match the
three leaves of module i to those of module i+1, including the final
join. Relabel the three leaves independently in successive modules.
This makes every matching except the final one the identity; the last
matching is the product `pi` of the original permutations, up to
conjugacy. The order of composition is the order in which a strand
travels around the ring. This is also visible from (1): T commutes with
every permutation of strand positions.

Thus the three conjugacy classes of the resulting permutation classify
these marked ring wirings up to independent leaf relabelling. For
`k>=3`, their colouring obstruction is precisely
`k odd and pi odd`. Every odd permutation in S3 is a transposition, so
every obstructed graph in this entire class is isomorphic to J_k and
inherits (8). In the raw labelled universe of `6^k` choices of junction
permutations, exactly half obstruct when k is odd, and none obstruct
when k is even. This counts wiring choices, not graph isomorphism classes.
It applies to a cyclic arrangement of claws with three wires per
junction, not to arbitrary complete metagraph wiring.

**Corollary 3.1 (PROVED; explicit minor chain).** For every odd `k>=5`,

```text
J_k -> J_(k-2) -> ... -> J_3 -> Petersen                    (8)
```

is a sequence of obstruction-preserving pole replacements and minor
operations. Each initial step uses Theorem 2. In J3 contract the triangle
`b_0,b_1,b_2` to one vertex. A triangle with three boundary edges extends
exactly when those three colours are distinct, and then uniquely, so this
last contraction also preserves the full boundary relation.

For a completely explicit final isomorphism to Petersen's standard labels,
send the contracted triangle to 0, send `(a_0,a_1,a_2)` to `(1,4,5)`, send
`(c_0,c_1,c_2)` to `(2,3,8)`, and send `(d_0,d_1,d_2)` to `(6,9,7)`.
The resulting edge set is the inherited fifteen-edge Petersen graph.

This supplies the requested implication for an infinite larger-snark
family through actual reductions, rather than importing a general
Petersen-minor theorem.

## 4. Closed formulas for every flower adjacent-edge signature

For odd `k>=3`, use the four-port convention of THM-4116: delete the
endpoints of an edge, and list the two remaining semiedges at one endpoint
followed by the two at the other. The four colour-orbit representatives
are `0000,0011,0101,0110`. The signature is `(m_e,m_e,0,0)`.

**Proposition 4 (PROVED).**

| Edge type | Number of such edges | `m_e` |
|---|---:|---:|
| `a_i b_i` | k | `(2^k+4)/3` |
| `a_i c_i`, `a_i d_i` | 2k | `(2^k-2)/3` |
| b-cycle edge | k | `(2^(k-1)-4)/3` |
| c/d-cycle edge | 2k | `(2^(k-1)+2)/3` |

The graph's rotations, including the c/d switch at the twisted join, and
the global c/d interchange act transitively on each listed class. Thus
these are all edges, not four arbitrarily selected samples.

**Proof.** It suffices to count the all-zero boundary word, because
Proposition 3 and the parity/Kempe proof in THM-4116 give the other three
orbit counts.

For a deleted spoke there remain `k-1` full modules. Put
`alpha=(2^(k-1)-1)/3` in (4). At the damaged module the two intact strand
leaves have their spoke semiedges coloured 0, so the incoming and outgoing
colours at each must be 1 and 2 in opposite order. The strand of the
deleted leaf has both outer colours 0. If its strand is b, exactly two
local assignments connect equal coordinates in the appropriate block
after the closing c/d swap; the other two assignments connect different
blocks. The count is `2(1+alpha)`. If its strand is c or d, exactly two
assignments connect distinct coordinates in the same active block; the
others connect different blocks. The count is `2alpha`. These are the
first two rows.

For a deleted horizontal edge, its two endpoint leaves disappear. The two
adjacent centres each have a boundary spoke coloured 0. The other two
spokes at each centre must be 1 and 2 in opposite order. There remain
`r=k-2` full modules, with r odd; put `beta=(2^r+1)/3` in (5).
The deleted strand has colour 0 at both outer ends.

The complete local two-module possibilities on the other two strands
are: `(0,0)` returns to itself in two ways; a word with exactly one zero
and one of 1,2 returns to itself in one way; `(1,2)` or `(2,1)` can return
to itself or reverse in one way each; a repeated nonzero word has no
extension. Constant triples and the one-zero pair cases cannot close
through the permitted odd-step blocks. For a deleted b-strand edge, the
two reversed all-distinct cases contribute `beta-1` each. For a deleted c
or d edge they contribute `beta` each. Thus the last two rows are
`2(beta-1)` and `2beta`, as claimed. QED.

For k=3 the b-triangle edges have `m_e=0`; all other entries are positive.
For every odd `k>=5`, every edge has `m_e>0`. This identifies the exact
failure boundary of the next reduction attempt.

## 5. Why deleting just two adjacent vertices destroys the obstruction

**Proposition 5 (PROVED).** Let a cubic graph's complementary four-pole at
an edge have signature `(m,m,0,0)`, with `m>0`. Remove the two endpoints
and join the four ports in the three possible perfect pairings. The three
labelled Tait counts are

```text
within-endpoint pairing: 9m;
either cross pairing:   3m.                                (9)
```

**Proof.** In the within-endpoint closure all nine words `(a,a,b,b)` are
allowed. In either cross closure only the three all-equal words in that
support are allowed. Each word has multiplicity m. QED.

Thus all three elementary smoothings are colourable on every edge of an
odd flower with `k>=5`. The four-module reduction works by coordinating
more vertices while retaining a six-port relation; no sequence whose
first step is one of these three smoothings can preserve uncolourability.

The exact Blanusa examples behave the same way. At edge `(7,8)` of
Blanusa I, `m=4` and the smoothing counts are `(36,12,12)`. At edge `(9,15)`
of Blanusa II, `m=6` and they are `(54,18,18)`. Each selected edge is an
actual branch edge of the inherited authenticated Petersen subdivision.
Keeping its endpoints, deleting the other unused edges, and suppressing
the recorded path interiors gives a **rooted** Petersen four-pole minor.
The original complement and Petersen complement both have exactly
`B={(a,a,b,b)}` as their support. Replacing one by the other therefore
preserves extendability in every four-port exterior. This is a certified
reduction for those examples, not a new derivation of their minor models.

## 6. Support alone cannot force the minor

Take two copies of K4, subdivide one edge of each, and join the two new
subdivision vertices by an edge. The result is a simple cubic graph on ten
vertices and fifteen edges with a bridge. Deleting the endpoints of the
bridge leaves two independent K4-minus-edge two-poles. Each requires its
two boundary colours to agree, and has two extensions per common colour.
Consequently the four-pole signature is

```text
(4,4,0,0), with support B={(a,a,b,b)}.                      (10)
```

This is exactly Petersen's support. Yet this graph has no Petersen minor:
a minor retaining ten vertices cannot contract or delete a vertex; a
fifteen-edge Petersen minor would then have to retain all fifteen edges
as well. That would make the graph itself Petersen, contradicted by its
bridge (or its triangles).

The first failed implication is therefore “Petersen boundary support
implies Petersen minor.” The strongest survivor is interchangeability in
colouring contexts. The missing additional coordinates are the structural
graph assumptions and an authenticated minor operation/model.

## 7. The Pythagorean connection survives as a predicate

For ordered boundary extension vectors f and g, the parallelogram
identities give

```text
||f-g||^2 + ||f+g||^2 = 2(||f||^2+||g||^2),
||f+g||^2 - ||f-g||^2 = 4<f,g>.                           (11)
```

The vectors f+g and f-g themselves are perpendicular only when f and g
have equal norm. To obtain an actual orthogonal pair for arbitrary
nonzero f,g, retain the tensor construction from the preceding synthesis:
set `U=f tensor g`, `V=g tensor f`, `M=||f||^2||g||^2`, and `C=<f,g>`.
Then `||U||^2=||V||^2=M`, `<U,V>=C^2`, and

```text
(U+V) perpendicular (U-V),
||U+V||^2 = 2(M+C^2),   ||U-V||^2 = 2(M-C^2).             (12)
```

These form a possibly degenerate real right triangle: proportional f,g
give a zero difference leg. For nonzero f,g the triangle is nondegenerate
and isosceles exactly when the glued graph is uncolourable.
Support-equivalent pole replacements preserve
that zero-overlap predicate in every fixed exterior. The positive entries
in (3) change, so lengths and nonzero overlap angles need not be preserved.
What the flower reduction carries down to Petersen is the isosceles
endpoint, not an integer Pythagorean triple or a literal isometry. The
repeated square-root/integer distinction in the previous session remains
necessary.

There is a separate primitive-integer chart, obtained in the root's
independent audit. The block multiplicity
`alpha_h=(4^h-1)/3` in (4) gives the two distinct block eigenvalues
`1` and `1+3alpha_h=4^h`. Euclid's map with parameters `(2^h,1)` therefore
gives

```text
(3alpha_h, 2^(h+1), 3alpha_h+2)
    = (4^h-1, 2^(h+1), 4^h+1),             h>=1.           (13)
```

The parameters are coprime and of opposite parity, so these triples are
primitive. They begin `(3,4,5)`, `(15,8,17)`, `(63,16,65)`. Thus the
arithmetic parameter is an actual extension multiplicity in the reduction
matrix. This family flattens as h grows; it does not approach the
isosceles endpoint in (12). No Pell recurrence or graph isomorphism with
an arithmetic object is asserted.

## 8. Exact checks and reproduction

The companion
[script](../../04-computation/boundary_support_20261020_snarks.py)
and [saved output](../../04-computation/boundary_support_20261020_snarks.out)
use Python's standard library. They explicitly enumerate all 27 incoming
and 27 outgoing local words, check the entire T-squared block identity,
check closure formulas through k=41, and authenticate explicit minor
chains through J101 (404 vertices).

An independent engine enumerates perfect matchings and counts even
complementary two-factors, without recursive edge colouring. It agrees
with the transfer count for every k from 3 through 15, up to 60 vertices,
and checks untwisted controls through k=8. A separate fixed-boundary
edge-colouring recursion checks all four edge types for odd k through 15.
The inherited perfect-matching four-pole engine independently checks both
marked Blanusa profiles. Positive controls are K4 and the even/untwisted
flower closures; hostile controls are Petersen, constant boundary words,
all three elementary smoothings, and the bridged pair of K4s.

```text
python3 -B 04-computation/boundary_support_20261020_snarks.py
python3 -B -O 04-computation/boundary_support_20261020_snarks.py
```

The root agent independently reconstructed T from two local definitions,
checked the block algebra, and independently checked twisted closure
counts in
[the independent audit](../../04-computation/boundary_support_20261020_independent.py).
That audit is separate from the larger-graph engine here.

**Remaining open obligation.** A general reduction proof for arbitrary
bridgeless cubic snarks would need a theorem guaranteeing an applicable
support-preserving minor replacement or an explicit Petersen model when
none applies. The flower strip identity and the marked Blanusa models do
not supply that global structural dichotomy.
