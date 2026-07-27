---
id: THM-2606
title: "Affine V4 parity channels, partial cubes, and the Feuerbach origin boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On an affine V4
  torsor the three nonzero difference
  parities are exactly the three translation-invariant C4 partial cubes,
  the three omitted V4 directions, and the three 2+2 quartic-resolvent
  channels.  GL(2,F2)=S3 acts faithfully on this three-set: an order-three
  element cycles the binary structures and an involution fixes one and
  swaps two.  Full AGL(2,F2)=S4 invariance leaves only the empty graph and
  K4, so no connected invariant partial cube exists; after choosing an
  origin, the unique connected GL-invariant partial cube is K1,3.  The six
  unordered pairs carry the octahedral incident/disjoint association
  scheme and no S4-invariant tournament; their four opposition-compatible
  C6 partial cubes are themselves equivalent to choosing an origin.
  Feuerbach tangency supplies precisely such an origin on the in/excircle
  sign labels, while four-body Jacobi squares realize the three resolvent
  channels.  These are exact label/kinematic correspondences, not N-body,
  supergroup, Keller, or LRC dynamics.
source: codex-2026-07-27-affine-v4-geometry
depends_on:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
script: 04-computation/v4_parity_partial_cube_feuerbach_thm2606.py
output: 05-knowledge/results/v4_parity_partial_cube_feuerbach_thm2606.out
script_sha256: 59f01db7ef84c572ef14876b4a85979a4a499c1ce1c90bc6e64245b1dfa6f0f7
output_sha256: a9c925aa8a3c98d908f01be5a35a4c237dfe3ba9e5bff887c02e2034ca72c1d7
hash_basis: working-tree bytes (LF)
---

# THM-2606 -- three binary faces on one affine four-state object

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The group/graph
census and all displayed finite controls
are reproduced by the dependency-free exact companion.  The quartic
application uses THM-2598, and the marked modular action uses THM-2595.
The Feuerbach paragraph uses the classical tangency theorem only to identify
one sign label; it does not turn the four circles into a geometric group.

## 1. Parity channel = omitted direction = `C4`

Let `X` be an affine torsor for

```text
V=F_2^2={0,a,b,c},                 c=a+b.                       (1)
```

A **parity channel** on `X` means a nonzero covector
`ell in V^*`.  It supplies an unordered affine bipartition: an absolute
function `p:X->F_2` with

```text
p(x+v)-p(x)=ell(v)                                             (2)
```

exists, and its only other choice is `p+1`.  Thus a channel does not by
itself choose which class is even.  Its kernel is

```text
ker(ell)={0,c_ell}                                             (3)
```

for a unique nonzero direction `c_ell`.  Define a graph on `X` by

```text
x ~_ell y       iff       ell(y-x)=1.                          (4)
```

The connection set in (4) is the two-element complement
`V\ker(ell)`.  Those two directions are a basis of `V`, so their coefficient
map embeds (4) isometrically as the square `Q_2`.  Hence `Gamma_ell` is a
`C4` partial cube.  Its two missing diagonals are exactly

```text
{ {x,x+c_ell}: x in X }/duplication,                           (5)
```

the two orbits of translation by `c_ell`.

Conversely, a translation-invariant simple graph on `X` is a Cayley graph
with connection set `D subset V\{0}`.  It is connected and a partial cube
exactly when `|D|=2`: sizes `0,1` are disconnected, size `2` is (4), and
size `3` is `K4`, which contains triangles.  Therefore the following
three-element sets are canonically equivariantly identified:

```text
nonzero ell in V^*
  <--> nonzero omitted directions c_ell in V
  <--> unordered 2+2 affine partitions of X
  <--> connected translation-invariant partial cubes C4 on X. (6)
```

This is the exact binary structure carried by the affine four-state object.
There are three such structures, not one canonical binary tree.

## 2. Full affine symmetry forbids the choice; an origin repairs it

The affine automorphism group is

```text
AGL(V)=V semidirect GL(V)=V4 semidirect S3 isomorphic to S4.    (7)
```

For an `AGL(V)`-invariant graph the connection set `D` above must be
`GL(V)`-invariant.  Since `GL(V)` is transitive on the three nonzero
directions, `D` is empty or all of `V\{0}`.  Thus the complete classification
is

```text
AGL-invariant graphs on X:                 empty, K4.           (8)
```

The only connected one is `K4`, hence no connected `AGL(V)`-invariant
partial cube exists.  More strongly, no translation-invariant tournament
exists on `X`: translation by `y-x` interchanges any prescribed pair
`x != y` and reverses its putative edge.

Now choose an origin `o in X`.  The remaining linear group `GL(V)=S3`
has two edge orbits:

```text
{ {o,x}: x!=o }                       (the three star edges),
{ {x,y}: x,y!=o, x!=y }               (the opposite triangle). (9)
```

Consequently

```text
GL-invariant graphs after choosing o:
empty, K1,3, K3 plus an isolated o, K4.                         (10)
```

The unique connected partial cube in (10) is the star `K1,3`, embedded as
`000` with leaves `100,010,001`.  Thus a distinguished origin is an exact
repair, not decorative metadata.

Existential graceful labeling does **not** select the lost sidecar.  All
three relevant graphs are graceful:

```text
K4:    vertex labels {0,1,4,6},       differences {1,...,6};
C4:    cyclic labels 0,4,1,2,         differences {1,...,4};
K1,3:  centre 0, leaves 1,2,3,        differences {1,2,3}.      (11)
```

Hence gracefulness distinguishes neither affine symmetry from a chosen
origin nor one parity channel from the other two.

## 3. The literal `C2/C3` co-occurrence frame

`GL(V)=S3` acts faithfully on each of the three sets in (6).  An element of
order three cycles all three parity/`C4` channels.  An involution fixes one
channel and swaps the other two.  For the marked matrices of THM-2595 this
is exactly

```text
C2*C3 --> GL(2,F2)=S3 --> Sym(V^*\{0}).                        (12)
```

Thus three binary faces and one ternary permutation really do live on one
object.  The statement is finite and sharp: (12) factors through `S3` and
forgets the free-product word/Bass--Serre coordinate.  An order-three cycle
of three binary structures is not a rooted ternary tree; THM-2596 separately
types the latter as a `PGL_2(Z)` reduction cross-section.

A nontrivial `Z_2` grading of the **pointed** group `V` is likewise a
nonzero `ell:V->F_2`.  There are three, permuted transitively by `S3`, and no
nonzero grading is `S3`-invariant.  This is the exact content of the
super-parity analogy.  It does not construct a Lie supergroup, a `Z_3`
grading, or any physical supersymmetry.

## 4. Quartic resolvent channels are exactly the three diagonals

On the separable depressed open of THM-2598, the four reconstruction
sections form an affine `V4` torsor `X`.  Its three nonidentity translations
are

```text
(12)(34),                 (13)(24),                 (14)(23), (13)
```

the three perfect matchings/`2+2` partitions of four root labels.  They also
index the three roots of the cubic resolvent.  For the channel represented
by `c` in (13), the unique covector with `ker(ell)={0,c}` gives the `C4` in
(4); its nonedges are exactly the two `c`-translation pairs.  Therefore

```text
one resolvent root/channel
  <--> one nonzero V4 direction
  <--> one 2+2 partition of the four sections
  <--> one C4 partial-cube structure.                           (14)
```

The `S3=S4/V4` resolvent action permutes the three choices.  Full `S4`
monodromy consequently cannot choose one.  The stabilizer of a chosen
channel is `V4 semidirect C2`, the order-eight automorphism group of its
square.  A chosen quartic section instead supplies an origin and hence the
star of (10).  These are different repairs: a channel gives a square but no
absolute even class; a section gives an origin but no preferred channel.

This refines THM-2598's lost-origin boundary.  It does not supply a rational
section of the quartic source, make the cubic resolvent a Keller map, or
exclude a degree-four point-cap map.

## 5. Six pair vertices: an association scheme, not a tournament

Let

```text
E=binom(X,2),                         |E|=6.                    (15)
```

On unordered pairs of distinct vertices of `E`, the `S4` action has exactly
two orbits: the original `K4` edges either intersect or are disjoint.  The
intersection graph is the octahedral graph `L(K4)` and the disjointness
graph is its complement, three independent edges.  Both orbitals are
self-paired.  Indeed, a transposition swaps two incident edges, and a double
transposition swaps two disjoint edges.  Therefore no `S4`-invariant
tournament exists on these six vertices.  The intrinsic binary carrier is
the two-relation incident/disjoint association scheme.

There is nevertheless a precise six-cycle repair.  For an origin `o`, put
an edge in `C_o` between the two edge-vertices

```text
ov  --  vw             (o,v,w pairwise distinct).              (16)
```

Every vertex in `E` has degree two in (16), and the graph is connected, so
`C_o` is `C6`.  Its antipodal vertices are exactly disjoint edges of `K4`,
and it is a partial cube (cyclic embedding
`000,100,110,111,011,001`).  Conversely, an antipode-compatible cyclic
order chooses one element from each of the three disjoint edge pairs in its
first half; there are

```text
3! * 2^3 / 12 = 4                                             (17)
```

unoriented cycles.  The four graphs in (17) are precisely the `C_o`, and
`o -> C_o` is `S4`-equivariant.  Thus identifying an abstract six-vertex
bicycle/partial cube with the six-edge quartic carrier spends exactly one
quartic-origin choice.

This boundary is also invisible to gracefulness.  `C6` is not graceful:
modulo two, the sum of its required edge differences would be
`1+...+6=21`, while it is also the sum of degree-two vertex labels and hence
even.  The graceful `K4/C4/K1,3` controls in (11) and the nongraceful `C6`
show that graceful and partial-cube predicates do not recover each other.

## 6. Feuerbach gives the origin, not a V4 motion

For a nondegenerate, non-equilateral labelled real triangle with positive side lengths
`a,b,c`, the incenter and excenters have homogeneous barycentric sign labels

```text
I   =( a: b: c),          I_A=(-a: b: c),
I_B =( a:-b: c),          I_C=( a: b:-c).                     (18)
```

Modulo simultaneous sign reversal, the four sign classes form

```text
{+/-1}^3 / <(-1,-1,-1)>  isomorphic to V4,                    (19)
```

with `I` as the all-positive identity and coordinate permutation `S3`
fixing `I` while permuting the three excenter labels.

The classical Feuerbach theorem says that the nine-point circle is
internally tangent to the incircle and externally tangent to the three
excircles.  Hence the internal/external tangency-type label singles out
exactly the identity in (19).  After that origin is retained, the unique
connected `S3`-invariant partial cube is the star in (10).

What transfers is therefore exact and limited:

```text
Feuerbach tangency type  --> distinguished V4 sign origin
                         --> canonical abstract K1,3 carrier. (20)
```

Multiplication of sign labels in (19) is not a motion taking one circle to
another, and the star in (20) is not the circle-tangency graph.  No Euclidean
metric, radius identity, or triangle dynamics transfers to the quartic.

## 7. Four-body Jacobi squares are the resolvent cubic

There is also a literal scalar four-body kinematic model.  Let

```text
f(T)=product_i(T-x_i)=T^4+pT^2+qT+r,       sum_i x_i=0,        (21)
```

and define the three pairing coordinates

```text
s_1=x_1+x_2,       s_2=x_1+x_3,       s_3=x_1+x_4,
u_i=s_i^2.                                                       (22)
```

Then

```text
(U-u_1)(U-u_2)(U-u_3)
 =U^3+2pU^2+(p^2-4r)U-q^2,        s_1s_2s_3=-q.               (23)
```

The centred Hadamard/Jacobi coordinates

```text
X=s_1*2=x_1+x_2-x_3-x_4,
Y=s_2*2=x_1-x_2+x_3-x_4,
Z=s_3*2=x_1-x_2-x_3+x_4                              (24)
```

satisfy

```text
X^2+Y^2+Z^2=4 sum_i x_i^2,       u_1+u_2+u_3=sum_i x_i^2.     (25)
```

The three double transpositions in (13) act on `(X,Y,Z)` by the three
even sign flips, so their squares descend to the `S3` channel permutation.
This is exactly the cubic resolvent, not merely a cluster-count analogy.

On the real chamber `x_1<x_2<x_3<x_4`, write the positive consecutive gaps
as `A,B,C`.  Direct centering gives

```text
s_1=-(A+2B+C)/2,       s_2=-(A+C)/2,       s_3=(C-A)/2.        (26)
```

Thus when `q!=0`,

```text
u_1>u_2>u_3>0,                                             (27)
```

and the ordered real-root chamber supplies a local origin/sign section.
The repair is order-theoretic and cannot be analytically continued through
the complex discriminant complement without monodromy.  Equations
(21)--(27) describe collinear centred coordinates only; they imply no
four-body Hamiltonian, force law, collision regularization, integrability,
or general `N`-body theorem.

## 8. Consequence and stopping rule

The user's `2/3/4/6` co-occurrence has one exact finite core:

```text
four affine states X
  --three nonzero differences/parities--> three binary C4 faces
  --GL2(F2)=S3--------------------------> C2 fixes one, C3 cycles three
  --unordered pairs---------------------> six octahedral edge states. (28)
```

But every promoted structure spends a sidecar:

| desired carrier | required choice | information still lost |
|---|---|---|
| one `C4` parity face | one resolvent channel/covector | absolute even class, section |
| invariant `K1,3` | one origin/section | preferred channel |
| one edge-state `C6` | one origin | tournament orientation |
| a six-vertex tournament | non-invariant ordered gauge | affine symmetry |

Therefore neither the modular finite quotient, Feuerbach labels, a graceful
labeling, nor the Jacobi square cubic supplies the missing algebraic quartic
section.  The cheapest valid G1 follow-up remains THM-2598's normalized
resolvent-order/branch test; the cheapest modular follow-up remains
THM-2596's Gram-owner cocycle.  No LRC row, Keller stratum, or knot family is
changed by this theorem.

## 9. Reproduction

Run

```bash
python 04-computation/v4_parity_partial_cube_feuerbach_thm2606.py
python -O 04-computation/v4_parity_partial_cube_feuerbach_thm2606.py
```

The dependency-free companion enumerates `GL(2,F2)`, all affine
permutations, the three parity squares, both invariant-graph classifications,
all tournaments on four and six vertices, the six-edge association scheme,
all four opposition-compatible hexagons, the graceful controls, the
Feuerbach sign quotient, and exact integer four-body/resolvent controls.
Normal and optimized transcripts must byte-match the stored output and end
in `FAILED CHECKS: NONE`.
