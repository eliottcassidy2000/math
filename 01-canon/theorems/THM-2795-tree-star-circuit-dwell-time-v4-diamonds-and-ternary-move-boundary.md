---
id: THM-2795
title: "Tree star-circuit dwell time, V4 diamonds, and ternary move boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The vertex-star vectors of an m-edge tree form one minimal binary
  circuit.  Every vertex order is its simple closed prefix walk, and each
  edge difference is exactly that cut coordinate's dwell time; gracefulness
  means the m dwell times are 1,...,m.  An adjacent swap changes one prefix
  through an affine V4 diamond.  On three consecutive vertices the binary
  swap and ternary rotation generate the local S3 quotient of C2*C3, not a
  free modular action.  Graceful-preserving move graphs are disconnected
  already on P4, and stars give an all-size two-component boundary.
source: root/tree-star-circuit-v4-c3-2026-07-28
depends_on:
  - THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion
  - THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction
related:
  - THM-2768-modular-c2-c3-a4-s4-bass-serre-quotient
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
  - THM-2777-marked-d3-six-root-determinant-tournament-and-cokernel-spectrum
script: 04-computation/tree_star_circuit_v4_c3_moves_thm2795.py
output: 05-knowledge/results/tree_star_circuit_v4_c3_moves_thm2795.out
script_sha256: fceb90a2246d60c5595cedc3cb3881446d8a04c3bc5d076d55f8229a665c3289
output_sha256: f90fd061d5b0cc930a6d11cf5ab0c40c9ec7b24bffed6df6108244a5328d9549
hash_basis: LF-normalized bytes
---

# THM-2795 -- tree star-circuit dwell time, V4 diamonds, and ternary move boundary

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2787 turns graceful labeling into bounded signed path integration.
THM-2793 turns its gap matrix back into a rooted vertex order.  There is a
third exact representation between them: order the vertex-star vectors in
the binary cut space and record how long each edge coordinate remains on.

This representation realizes the suggested binary/ternary grammar on one
object.  The binary move is literally an affine `V4` diamond.  A ternary
rotation joins it into the six-state `S3` permutohedron.  The local action is
a quotient of `C2*C3`, however, and the graceful locus is highly
disconnected.

## 1. The vertex stars are one minimal binary circuit

Let `T=(V,E)` be a tree with

```text
|E|=m,                       |V|=m+1.
```

For each vertex `v`, let

```text
s_v in F_2^E
```

be its edge-incidence star:

```text
(s_v)_e=1 iff e is incident with v.                     (1)
```

Every edge has two endpoints, so

```text
sum_(v in V) s_v=0.                                     (2)
```

For `A subset V`,

```text
sum_(v in A)s_v=delta_T(A),                              (3)
```

the binary cut vector.  Since `T` is connected, `(3)` is zero only for
`A=empty` or `A=V`.  Hence the labelled family

```text
{s_v:v in V}                                             (4)
```

is a minimal binary circuit: it has rank `m`, and `(2)` is its unique
nonzero dependence.

For `m>=2`, all `s_v` are distinct and nonzero.  Equality of two vertex
stars in a simple connected tree would force the tree to consist of the
single edge joining those vertices.  The one-edge tree is the labelled
two-element circuit boundary; its two generator labels represent the same
nonzero vector.

## 2. Vertex orders are closed cut-prefix walks

Choose a vertex order

```text
v_0,v_1,...,v_m
```

and define

```text
c_t=sum_(j<t)s_(v_j)=delta_T({v_0,...,v_(t-1)})
        in F_2^E,                  0<=t<=m+1.             (5)
```

Then

```text
c_0=c_(m+1)=0,                                          (6)
```

while `c_0,...,c_m` are pairwise distinct.  Otherwise the star sum over
one nonempty proper consecutive block of the order would vanish, contrary
to minimality.  For `m>=2`, `(5)` is therefore a simple length-`m+1`
closed walk in the Cayley graph of the cut space with labelled generators
`s_v`, using every generator once.

Give `v_j` the label

```text
f(v_j)=j.                                                (7)
```

For an edge `e={v_a,v_b}` with `a<b`, its coordinate in `(5)` obeys

```text
(c_t)_e=1 iff a<t<=b.                                   (8)
```

Its binary dwell time is consequently

```text
D_e=sum_(t=1)^m(c_t)_e=b-a=|f(v_a)-f(v_b)|.              (9)
```

Thus:

> **Star-circuit dwell theorem.** A vertex order is graceful if and only if
> its `m` coordinate dwell times are exactly `{1,...,m}`.

This is equivalent to THM-2787, but it exposes a different operation:
graceful labeling is a sequencing problem for one minimal binary circuit.
The consecutive-ones gap matrix is precisely the coordinate history of
the prefix walk `(5)`.

## 3. Adjacent swaps are affine V4 diamonds

Suppose consecutive positions in the order contain distinct vertices
`a,b`, and let `x` be the prefix before them.  The two local paths are

```text
x -> x+s_a -> x+s_a+s_b,
x -> x+s_b -> x+s_a+s_b.                                (10)
```

For `m>=2`, `s_a,s_b` are linearly independent, so

```text
{x,x+s_a,x+s_b,x+s_a+s_b}                               (11)
```

is an affine `V4=F_2^2` torsor.  Swapping `a,b` changes exactly one
intermediate prefix state.

The resulting integer dwell-time law is equally local.  For each edge `e`,

```text
D'_e-D_e in {-1,0,1},                                   (12)
```

and it is nonzero exactly when `e` is incident with exactly one of `a,b`:

```text
D'_e!=D_e iff (s_a+s_b)_e=1.                             (13)
```

If `e={a,b}`, both star coordinates are one and its dwell time is
unchanged.  Equations `(10)--(13)` are the exact `V4` tree move; no
metaphorical torsor is being inferred from the prime two alone.

## 4. Three-place rotations give S3, not the free modular group

On three consecutive vertices `(a,b,c)`, let

```text
s:(a,b,c)->(b,a,c),
r:(a,b,c)->(b,c,a).                                     (14)
```

The six local orders form the `A2` permutohedron.  Directly,

```text
s^2=r^3=(sr)^2=1,                                      (15)
```

and `s,r` act transitively on all six orders.  Hence their local action is

```text
S_3 = (C_2*C_3)/<(sr)^2>.                               (16)
```

This is the precise relation to THM-2768's modular grammar.  The binary and
ternary generators co-occur on the star-circuit order, but the extra face
relation in `(15)` means this is not a faithful action of
`PSL_2(Z)=C_2*C_3`.  Nor does the local `S3` action preserve gracefulness
automatically.

## 5. The graceful move locus is disconnected

Define two finite move graphs for a fixed tree:

```text
vertices: graceful vertex orders;
binary edges: one adjacent swap preserving gracefulness;
ternary edges: one consecutive triple rotation preserving gracefulness.
                                                               (17)
```

The exact companion obtains the following aggregate census over all
nonisomorphic trees:

```text
n                         2    3    4    5     6      7       8
graceful orders           2    4   16   68   376   2184   15096

binary components         1    2    6   18    94    402    2304
binary isolated           0    0    4   12    72    260    1554

binary+ternary components 1    1    6   14    62    334    2000
binary+ternary isolated   0    0    4    8    40    202    1244.       (18)
```

The first obstruction is already the four-vertex path with edges

```text
{(1,0),(1,2),(0,3)}.                                   (19)
```

The graceful order

```text
(0,2,1,3)                                               (20)
```

has no graceful adjacent swap and no graceful consecutive triple rotation.
Thus `V4` diamonds plus local `C3` faces do not connect even the smallest
nonstar graceful locus.

There is also an all-size boundary.  For the star `K_(1,m)`, a graceful
order must put its centre at label zero or `m`; an interior centre cannot
produce edge difference `m`.  When `m>=3`, no graceful-preserving move in
`(17)` can carry the centre through an interior position.  The two extreme
classes are each connected by leaf swaps, so the mixed move graph has
exactly

```text
two components of size m!                               (21)
```

for every such star.  The ternary move joins the two ends only in the
three-vertex boundary `m=2`.

This proves a stopping rule: any induction or Markov argument based solely
on graceful-preserving local `C2/C3` moves is incomplete.  It needs a
controlled excursion through a defect stratum, a nonlocal block move, or an
additional reference that can cross `(21)`.

## 6. Exact verification and scope

Run

```bash
python 04-computation/tree_star_circuit_v4_c3_moves_thm2795.py
python -O 04-computation/tree_star_circuit_v4_c3_moves_thm2795.py
```

The assertion-free companion checks the circuit and dwell theorem on every
vertex order of every nonisomorphic tree through seven vertices.  It checks
`355,836` affine `V4` diamonds, the complete six-state local `S3` orbit,
and both move graphs through all 23 eight-vertex tree shapes.  It verifies
the isolated `P4` order and the star theorem through eight vertices.

```text
PROVED HERE (candidate):
  minimal binary vertex-star circuit;
  simple closed cut-prefix walk;
  edge difference equals coordinate dwell time;
  graceful iff dwell spectrum is 1,...,m;
  exact adjacent-swap V4 diamond and signed unit dwell law;
  local C2/C3 action is S3 with the face relation;
  exact graceful move-graph census through n=8;
  isolated P4 and all-size star disconnection boundaries.

NOT PROVED:
  connectivity after allowing controlled nongraceful defects;
  the marked gap-tail extension conjecture;
  the Graceful Tree Conjecture;
  a faithful PSL2(Z) action;
  an LRC, Keller, JC(2), DC(2), knot, or tournament consequence.       (22)
```

QED, conditional only on candidate status promotion after independent
hostile audit.
