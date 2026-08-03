---
id: THM-3260
title: "Bispanning reset-link holotopy atlas and C12/C7 carrier boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  THM-3254's 23 first-link-blocked row pairs form a connected bipartite
  graph with one leaf bridge.  Its leafless core has 12 vertices, 22 edges
  and first Betti number 11.  It is nonplanar, has intrinsic automorphism
  group C2, and has no faithful C12 vertex action.  Adjoining the eight
  delayed pairs gives a relative integral homology layer of rank six, but
  the full graph has no order-seven automorphism.  Positively, the core is bispanning:
  exactly 4,960 unordered complementary spanning-tree pairs form one
  connected 43,408-edge symmetric-exchange atlas.  Every chart gives an
  integral unimodular 11-dimensional polarization.  The C2 symmetry fixes
  no chart.  Hence the dimension-11 match with THM-3255's C12 augmentation
  defect is noncanonical, but a supplied C12 vertex label plus one tree-pair
  chart gives the exact missing integral sidecar.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins THM-3254's theorem, script and
  output plus THM-3255; parses only the two literal pair banks; reconstructs
  the graph; verifies connectivity, bridge, bipartition, Betti number, an
  explicit K3,3 subdivision and all degree-compatible automorphisms; exhausts
  the full delayed relative boundary and both uncoloured and edge-coloured
  automorphism groups; exhausts
  all 4,960 complementary tree pairs and every symmetric exchange; checks
  connectedness and the free C2 action; and builds an explicit pair of
  unimodular reduced incidence matrices with integral GL11 transition.
  Normal, optimized and stored replay and the declared LF hashes are required.
  An independent NetworkX/SymPy reconstruction rederived the original and
  delayed graph censuses, automorphism groups, relative boundary, tree-pair
  atlas and integral cyclic-difference typing; it replayed normal, optimized
  and stored bytes exactly.  Accepted after the faithful-action, integral
  lattice and loop-holonomy scope repairs.
depends_on:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
related:
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
script: 04-computation/gmc_bispanning_reset_link_holotopy_thm3260.py
output: 05-knowledge/results/gmc_bispanning_reset_link_holotopy_thm3260.out
script_sha256: 7adf81a692bf477d493860483f50b291c01dd267ed578caa12dd31bca9128616
output_sha256: eb0b2d2de4808ae1aabdede94640390aeb35cc7a7deec55b193da44dcd043e37
hash_basis: LF-normalized bytes
---

# THM-3260 -- bispanning reset-link holotopy atlas and C12/C7 carrier boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3254 proves that every lawful two-row covering pair has a sharp two-state
positive clutch, and that 23 of the 31 pairs already fail on the first link of
the reset.  THM-3255 independently finds an eleven-dimensional multiplicative
phase defect.  The equality of dimensions suggests a transport, but dimension
alone is not structure.

The reset-link graph supplies a surprisingly rich answer.  It has no planar
or faithful cyclic symmetry capable of identifying the two spaces canonically,
yet it admits a complete atlas of integral tree/cotree polarizations.  The
delayed edge layer repeats the phenomenon in rank six.  This is the precise
sense in which the bridge is a holotopy: locally trivial and globally connected
by exchanges, but with no distinguished chart or equivariant carrier.

## 1. The reset-link graph and its core

Let the vertices be the response-row labels occurring in THM-3254's 31
covering pairs.  Delete its eight delayed pairs and retain the 23 pairs whose
positive ratio line is already covered by reset-link traps.  Call the resulting
graph `G`.

Exact reconstruction gives

```text
|V(G)|=13, |E(G)|=23, beta_1(G)=11.                    (1)
```

The graph is connected.  It has one bridge, `(13,14)`, and row 14 is a leaf.
Delete that leaf and bridge.  The resulting core `G_0` has

```text
|V(G_0)|=12, |E(G_0)|=22, beta_1(G_0)=22-12+1=11.      (2)
```

Its bipartition is

```text
{2,11,16,18,22} | {3,7,10,13,17,19,21}.               (3)
```

Thus the core has exactly twice a spanning-tree edge count.

## 2. Two canonical routes fail

First, `G_0` is nonplanar.  An explicit subdivision of `K_(3,3)` has branch
parts

```text
{10,3,18} | {11,16,22}.                                (4)
```

Six branch paths are direct.  The remaining three are

```text
18-13-11,
18-19-2-21-16,
18-7-22.                                               (5)
```

Their interiors are pairwise disjoint and avoid the six branch vertices.
Consequently there is no planar-dual cycle/cut Hodge star.

Second, exhaustive degree-compatible permutation gives

```text
Aut(G_0)=C_2,                                           (6)
```

where the nontrivial element swaps rows 17 and 21 and fixes every other
vertex.  In particular the unlabeled graph has no order-12 automorphism, hence
no faithful `C_12` action and no `C_12`-torsor vertex action.  (It does admit
the nonfaithful parity action obtained by mapping `C_12` onto `(6)`.)  The equality

```text
beta_1(G_0)=dim Aug(Q[C_12])=11                         (7)
```

is therefore not an equivariant identification with THM-3255's missing
character module.

## 3. The delayed relative layer and the second dimension trap

Restore the eight delayed pairs and call the full 31-edge covering graph `H`.
It has

```text
|V(H)|=15, |E(H)|=31, beta_1(H)=17,                    (7a)
```

and its unique bridge is `(3,9)`.  Relative to `G`, the new vertices are
exactly 9 and 12.  The relative cellular complex has eight edge generators,
two vertex generators and boundary rank two.  Since `G` and `H` are connected,
the graph-pair exact sequence therefore gives

```text
H_1(H,G;Z) = H_1(H;Z)/H_1(G;Z) = Z^6.                 (7b)
```

Concretely, five delayed edges already have both endpoints in `G`; the two
edges incident to new row 12 contribute one further relative cycle, while
`(3,9)` and one row-12 edge attach the two new vertices.

This produces a second exact numerical resonance,

```text
rank H_1(H,G;Z)=dim Aug(Q[C_7])=6.                     (7c)
```

It is again not equivariant.  Exact degree-compatible enumeration gives

```text
Aut(H)=S_3({10,17,21}) x C_2({14,18}),
element orders = {1,2,3,6}.                            (7d)
```

Thus `H` has neither an order-seven nor an order-twelve automorphism.  If the
reset/delayed edge type is retained as a colour, its automorphism group drops
to the single `C_2` swapping rows 17 and 21.  Consequently even the combined
dimension identity

```text
beta_1(H)=11+6=dim Aug(Q[C_12])+dim Aug(Q[C_7])        (7e)
```

does not supply a canonical `C_12 x C_7` module, a splitting, or labels.

## 4. The bispanning atlas

Despite `(4)--(7)`, the 22 edges of `G_0` split into two edge-disjoint spanning
trees.  One exact ordered example is

```text
T = {(2,10),(2,13),(2,17),(2,19),(3,11),(3,16),
     (7,18),(10,11),(10,22),(13,18),(16,21)},

T'={(2,21),(3,22),(7,22),(10,16),(11,13),(13,22),
     (16,17),(17,22),(18,19),(19,22),(21,22)}.          (8)
```

Exhausting all 11-edge subsets and identifying a tree with its complement
gives exactly

```text
4,960 unordered decompositions E(G_0)=T disjoint_union T'. (9)
```

Their ordered edge-index digest is

```text
28285b3925c94703da8a0267642ffe8d20c30f5d76b443c52c146aa525c2b822. (10)
```

Join two decompositions when one can be obtained from the other by exchanging
one edge of `T` with one edge of `T'` while both sides remain trees.  On
unordered decompositions this symmetric-exchange graph has

```text
4,960 vertices,
43,408 edges,
minimum/maximum degree 13/27,
one connected component.                               (11)
```

Thus every tree-pair chart is reachable from every other by lawful local
exchanges.

## 5. Integral polarization in each chart

Choose a root vertex and orientations.  For a tree pair `(T,T')`, let `B_T`
and `B_T'` be the two reduced `11 by 11` incidence matrices.  Each is
unimodular.  For the chart `(8)`, rooted at row 22,

```text
det B_T=det B_T'=-1,
U=B_T^(-1)B_T' in GL_11(Z),
det U=1.                                                (12)
```

The exact matrix digest is

```text
4d692df61815d6f39e0dbe40259ab9fe6e0dc19ef26722497834ec25ad6019cc. (13)
```

The columns `(-Uz,z)` are the fundamental integral cycles obtained by adding
the `T'` edges to `T`.  Dually, either reduced incidence matrix identifies
vertex potentials modulo constants with integral tree-edge differences.
Every symmetric exchange changes these bases by an integral unimodular
transition.  Equation `(11)` therefore gives one connected `GL_11(Z)` atlas
of local cycle/cut frames.

This atlas has genuine chart ambiguity rather than a preferred origin.  The
nontrivial automorphism in `(6)` fixes none of the 4,960 unordered tree pairs;
it acts freely in 2,480 pairs of charts.  No tree decomposition is intrinsic
to the unlabeled response graph.  This statement does not assert nontrivial
loop holonomy: the direct basis-change matrices telescope around every closed
exchange loop unless an additional transport rule is supplied.

## 6. Conditional bridge to the missing `C_12` modes

Suppose two additional choices are supplied:

1. a bijection from `V(G_0)` to `C_12`, including a selected zero/root; and
2. one ordered tree-pair chart `(T,T')`.

Write `R=Z[C_12]=Z[x]/(x^12-1)` and `N=1+x+...+x^11`.  Multiplication by
`x-1` gives the integral cyclic-difference isomorphism

```text
R/(Z N) --(x-1)--> Aug(R).                             (14a)
```

The vertex label identifies `R/(Z N)` with integral vertex potentials modulo
constants.  Given `a in Aug(R)`, invert `(14a)`, choose the unique root-gauge
representative `p` whose selected-root coordinate is zero, and put

```text
z=B_T'^T p.
```

These are the `T'` edge differences.  Since `B_T'` is unimodular, and the
cycle boundary equation is

```text
B_T(-Uz)+B_T'z=0,
```

the assignment `a |--> (-Uz,z)` constructs an explicit integral isomorphism

```text
Aug(Z[C_12])  -->  H_1(G_0;Z).                         (14)
```

Tensoring `(14)` with `Q` gives the rational character-space bridge.  The
cyclic-difference step is load-bearing: directly identifying the integral
augmentation lattice with root-gauged potentials would hide an index-12
lattice ambiguity.

Changing the tree pair transports `(14)` through the connected exchange atlas.
Thus the numerical resonance in `(7)` is not empty: it becomes an exact bridge
once an external cyclic vertex label and a polarization chart are provided.

Neither sidecar is presently canonical.  The response rows do not come with
owner phases, and `(6)` cannot generate them.  The tree-pair atlas supplies
transport between chosen charts, not a distinguished choice.

## 7. Scope

The theorem uses only which two-row pairs are first-link blocked.  It does not
transport the signs, rational clutch thresholds, physical state labels or
positive cones attached to those edges.  It constructs no row-to-owner map,
no faithful `C_12` torsor action, no faithful/order-seven graph action on the
`C_7`-sized relative layer, no physical Gaussian-moment response, and no positive LRC
current.  The conditional isomorphism `(14)` is not a consequence of
THM-3254 or THM-3255 without its two explicit sidecars.

In particular, the shared number 11 is a proved source of a possible chart,
not a proof that the two eleven-dimensional spaces are the same.  The exact
gain is the replacement of a vague dimensional analogy by a complete
bispanning exchange atlas and a precise statement of the missing data.

No row exclusion, factorial staircase closure, arbitrary-radial NC2 theorem,
Gaussian Moment Conjecture, or `LRC(14)` decrement follows.

## 8. Exact companion

Run

```text
python 04-computation/gmc_bispanning_reset_link_holotopy_thm3260.py
python -O 04-computation/gmc_bispanning_reset_link_holotopy_thm3260.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
only exact integer and rational arithmetic.  It has no floating point,
randomness, discovery cache, graph-library dependency, or
optimization-sensitive assertion.

QED.
