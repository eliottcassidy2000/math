# High girth forces geometric crowding; fixed spacing forces three colours

2026-10-11. **PROVED:** the elementary growth/packing bounds and fixed-spacing
degeneracy theorem below. **CITED:** Tomon's high-girth preservation theorem,
taken as a preprint input rather than independently audited here.
**CONDITIONAL:** six/seven cascades using that input and the accepted H6
lower-bound premise. **OPEN:** the unrestricted six-versus-seven decision.

## Inheritance and the decisive change of scale

The closest mechanisms are the
[critical algebraic witness reduction](plane_algebraic_witness_20261011.md)
and [odd-cycle channel certificates](plane_odd_cycle_channels_20261011b.md).
The hostile is an arbitrarily long odd cycle, invisible to a fixed short-cycle
catalogue. The corrected near miss is treating local tree structure as a
global colouring theorem without retaining the geometry's resolution.
The least-used sidecar is the minimum separation between all distinct points,
not just the fixed distance one on edges.

**CITED input.** Istvan Tomon, *Geometric triangle-free graphs of large chromatic
number*, arXiv:2610.03517v1, Theorem 1.3: every chromatic number realized by a
unit-distance graph in dimension d>=2 is realized by an induced unit-distance
graph in the same dimension with arbitrarily large girth. Theorem 1.5 supplies
the combinatorial mechanism: induced subgraphs of suitable Cartesian powers
retain the chromatic number and have prescribed large girth.
[Primary source](https://arxiv.org/html/2610.03517v1).

This is not the categorical-product construction in the knight-torus work.
The retained dimension does not mean retained vertex separation, coordinate
size, or number of vertices. Those are the variables the next proof measures.

## An elementary Moore-growth / disk-packing bound

Let G be a finite simple graph realized by distinct points in the Euclidean
plane, with every edge of length one. Additional distance-one pairs may be
omitted: the argument works for either induced or non-induced unit graphs.
Assume minimum graph degree at least d>=2 and girth at least 2r+2, where r>=1
is an integer. Then a radius-r breadth-first ball has at least

\[
 M(d,r)=1+d\sum_{i=0}^{r-1}(d-1)^i
 =\begin{cases}1+2r,&d=2,\\
 1+d((d-1)^r-1)/(d-2),&d>2.
 \end{cases}                                                    \tag{1}
\]

Indeed, two distinct paths in the breadth-first tree cannot identify a
vertex at these depths, or add an internal non-tree edge: that would give
a cycle of length at most 2r+1. Each non-root vertex before the last level
has at least d-1 distinct children. All ball vertices lie in a Euclidean
disk of radius r about the root, by the triangle inequality along unit edges.

If the point set has pairwise separation at least h>0, place disjoint open
disks of radius h/2 around those vertices. They fit inside a disk of radius
r+h/2. Comparing areas gives

\[
 M(d,r)\le(1+2r/h)^2.                                           \tag{2}
\]

For the actual global minimum separation delta, therefore,

\[
 \boxed{\delta\le {2r\over\sqrt{M(d,r)}-1}.}                     \tag{3}
\]

No equality classification or sharpness is asserted. Overlap of graph edges
is irrelevant; only the disjoint auxiliary disks are used.

For a k-vertex-critical graph, every vertex has degree at least k-1: otherwise
a (k-1)-colouring after deletion extends back. Hence a k-critical unit graph,
k>=4, of girth at least 2r+2 satisfies (3) with

\[
 M_k(r)=1+(k-1)((k-2)^r-1)/(k-3),\qquad
 \delta=O_k\bigl(r(k-2)^{-r/2}\bigr).                             \tag{4}
\]

Vertex-critical means every vertex deletion reduces the chromatic number;
it does not assume an unproved geometric minimality principle. Induced
vertex deletion preserves both the geometry and the lower girth bound.

## Fixed separation plus high girth gives an actual upper bound

**Theorem.** Fix h>0. Choose any integer r>=1 such that

\[
 1+3(2^r-1)>(1+2r/h)^2.                                        \tag{5}
\]

Every finite h-separated unit-distance graph in the plane with girth at least
2r+2 is 2-degenerate, and thus three-colourable.

**Proof.** If a nonempty subgraph had minimum degree at least three, it
would inherit the separation and girth hypotheses. Equations (1)--(2) with
d=3 contradict (5). Thus every nonempty subgraph has a vertex of degree at
most two. Repeatedly delete one and colour in reverse order, always choosing
one of three colours absent from its at most two already coloured neighbours.
This also produces a full colouring certificate. QED.

An r satisfying (5) exists because the left side grows exponentially and
the right side quadratically. For a supplied positive rational h it can be
found using only exact rational inequalities. The method does not require
computing exact irrational closest-pair distances: any certified positive
rational lower bound on separation suffices.

Some sufficient, nonoptimal thresholds are:

| Minimum point separation h | r from (5) | Girth at least | Guaranteed colours |
|---:|---:|---:|---:|
| 1 | 6 | 14 | 3 |
| 1/2 | 9 | 20 | 3 |
| 1/10 | 15 | 32 | 3 |
| 1/100 | 23 | 48 | 3 |
| 1/1000 | 31 | 64 | 3 |

At h=1/1000 the sharper critical bounds already exclude a 6-critical
witness at girth>=32 and a 7-critical witness at girth>=28. These do not
require assuming that the original graph was critical: any graph of chromatic
number at least k contains an induced k-critical subgraph. Hence every such
h-separated unit graph is six-colourable at girth>=28 and five-colourable at
girth>=32. These are sufficient bounds, not claimed optimal constructions or
coordinate-bit lower bounds.

## The six/seven cascade, retaining its quantifiers

Using the cited high-girth theorem and compactness gives

\[
 \chi(\mathbb R^2)\le6
 \quad\Longleftrightarrow\quad
 \exists g\ \forall\text{ finite unit graphs }G\text{ of girth}\ge g,
 \ \chi(G)\le6.                                                 \tag{6}
\]

The forward direction is restriction. For the reverse, a finite graph of
chromatic number greater than six would have a high-girth counterpart of the
same chromatic number, contradicting the right side. Under H6, (6) is an
exact-six criterion. If the plane instead requires seven colours, there
are finite induced seven-chromatic witnesses of every prescribed girth;
critical extraction and the inherited real-algebraic transfer retain this
property. Their separations must obey (4).

The quantifier order is crucial. Our fixed-spacing result proves

\[
 \forall h>0\ \exists g(h)\ \forall G\text{ with separation}\ge h,
 \quad\operatorname{girth}(G)\ge g(h)\Rightarrow\chi(G)\le3.
\]

It does not prove one girth threshold uniformly over all positive h.
Substituting each finite graph's own positive separation and then silently
interchanging those quantifiers would be a false proof of a plane bound.
This identifies the missing uniformity coordinate precisely.

## Connection web and next experiment

* **Channels:** the odd-cycle hypergraph reduction remains exact. If the plane
  requires seven colours, no fixed cutoff on cycle lengths captures all its
  seven-colour witnesses.
* **Geometry:** high-girth critical growth forces near-collisions. A fixed
  lattice mesh h implies separation>=h, hence eventually enters the proved
  three-colour regime. Refining the mesh is mathematically necessary for
  representing the preserved high-girth chromatic obstructions.
* **Pants and boundary transport:** local pieces can look completely harmless
  while global cycle transport carries the obstruction. The parallel is an
  operation on constraints; no pants decomposition of the plane graph is
  presumed and no surface genus is identified with chromatic number.
* **Compression:** deleting short cycles or compressing a boundary does not
  authorize deletion of its long-cycle parity or its geometric source data.
  The surviving coordinate is what the external context can still observe.

The next concrete target is a multiscale boundary construction retaining exact
unit-edge equations, ternary labels and binary seam parity while refining
clusters of close points. It must either lift three labels globally or export
a finite algebraic obstruction. Merging nearby points without an exact
neighbourhood congruence is not valid: the unit-distance relation is equality,
not a continuous tolerance test. This is a target, not a completed algorithm.

## Reproduction and independent audit

```text
python -B -X utf8 04-computation/experiments/plane_girth_crowding_20261011b.py
python -O -B -X utf8 04-computation/experiments/plane_girth_crowding_20261011b.py
```

The 2,309 exact checks cover degrees 2..11 and radii 1..40 by independently
accumulated levels; 500 rational spacing/degree pairs with both the original
and denominator-cleared strict inequalities; the degree-two growth boundary;
equality, girth and typed-input hostiles. The output certifies arithmetic
bounds, not geometric existence. The proof of (1)--(5) and its hereditary
degeneracy consequence received a separate agent audit. Normal and optimized
runs agree. Saved output: [plane_girth_crowding_20261011b.out](plane_girth_crowding_20261011b.out).
