# Borsuk four: small geometric blockers and the order-realizability boundary

2026-10-14. **PROVED:** the finite direction/blocker alternative, the
eight-point four-dimensional case, and the marked order-transfer statements
below. **FINITE-EXACT:** the declared small configurations and complement
profiles. **CITED / ACCEPTED INPUT:** the supplied Hexagon asymptotic theorem.
**OPEN:** the unrestricted five-piece Borsuk assertion in dimension four.
No priority claim is made for the classical geometric or polyhedral ingredients.

## 1. Target, inheritance, and the useful paper interface

The [Conjectures.io target](https://conjectures.io/problems/borsuk-borsuk-conjecture-four)
is `Borsuk.BorsukConjecture 4`: every bounded, nontrivial subset of Euclidean
four-space has a cover by five sets of strictly smaller diameter. This is the
literal [formal source](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/Wikipedia/BorsukConjecture.lean),
not a theorem merely about a finite graph. This session makes no submission.

The closest repository bridge was the [operator/tensor prize route](collatz_operator_tensor_bridges_20261007.md).
Its maintained boundary is essential: the accepted projector construction
starts from lines in $\mathbb R^4$ but its trace-one symmetric-matrix carrier
has dimension nine. It does not refute Borsuk four. Other repository searches
for Borsuk primarily led to Borsuk--Ulam/topological results, which have
different hypotheses and conclusions. The useful older sidecar is the
[full boundary-colour relation](boundary_palette_compression_20261011.md),
rather than a count of colourings without their geometric ports.

The selected Hexagon input is [*An exponential improvement for Borsuk's
problem*, 2610.00163v1](https://hexagonmath.org/2610.00163v1),
[PDF](https://hexagonmath.org/pdf/2610.00163v1). Theorem 1.1, p.1, gives an
asymptotic exponential base $\sqrt{9536/6375}$. Section 2, p.3, reduces
partitioning to directions selecting strict local maxima in a near-diameter
graph. Sections 3--4, pp.4--7, use an antipodal pairing, a one-coordinate
radial lift and reflection in a self-dual cone to estimate the available
directions. We accept its stated theorem as requested; we do not re-audit it.
Its sufficiently-large-dimension hypothesis cannot be replaced by $n=4$.
The finite directional interface below is proved directly and uses no
asymptotic estimate. Relevant full pages 1--3 and the geometric construction
were read; pp.3 and 6 were also visually inspected.

The live board is: diameter versus near-diameter graphs / strict margins /
five direction classes / convex-dependence blockers / linear extensions /
the realized projection-order measure. The anchor is Borsuk four; the niche
is its small forbidden configuration; the wildcard is balanced witness
splitting and the golden-ratio threshold.

## 2. An actual finite Borsuk result: at most eight points

Let $X$ be a finite set of distinct points in $\mathbb R^d$, of diameter
$D>0$. Its diameter graph has an edge precisely at distance $D$. Because the
set is finite, a proper $k$-colouring is exactly a partition into at most
$k$ parts of diameter strictly below $D$. This equivalence uses finiteness.

**Lemma 1 (a forbidden near-clique).** A diameter graph in $\mathbb R^d$
cannot contain $K_{d+2}$ with one edge removed, even as a non-induced subgraph.

*Proof.* Its $d$ shared vertices form a regular $(d-1)$-simplex of edge $D$.
Let $c$ be its circumcentre. Its circumradius satisfies
$R^2=(d-1)D^2/(2d)$. A point at distance $D$ from all these vertices has
projection $c$ onto their affine span, and perpendicular displacement of
length

$$h=D\sqrt{\frac{d+1}{2d}}.$$

The perpendicular space has dimension one. The two tips therefore coincide,
or their distance is $2h=D\sqrt{2(d+1)/d}>D$. Both possibilities contradict
the marked distinct diameter configuration. $\square$

For $d=4$, the forbidden graph is $K_6-e$. Equivalently, the common
neighbourhood of any diameter triangle has maximum degree at most one:
a two-edge path in that neighbourhood would form the forbidden graph with
the triangle. This is a real geometric restriction, not a statement about
all unit-distance graphs.

**Theorem 2.** Every nontrivial set of at most eight points in $\mathbb R^4$
has a partition into five sets of strictly smaller diameter.

*Proof.* Suppose its diameter graph is not five-colourable, and choose a
vertex-minimal non-five-colourable induced subgraph $G$ with $n\le8$.
Deleting any vertex gives a five-colouring, so $\chi(G)=6$, and every vertex
has degree at least five. Thus its complement $H$ has maximum degree
$n-6\le2$. Each component of $H$ is an isolated vertex, a path, or a cycle.

Colouring $G$ is the same as partitioning $H$ into cliques. For an $s$-vertex
path, or an $s$-cycle other than $C_3$, this saves $\lfloor s/2\rfloor$
colours relative to singletons; a triangle saves two. Componentwise comparison
with the independence number gives the exact identity

$$\chi(G)=\alpha(H)+c_{\rm odd\ge5}(H), \tag{1}$$

where the last term counts odd-cycle components of length at least five.
By Lemma 1, $G$ has no $K_6$, so $\alpha(H)\le5$. There can be at most one
such odd component. To obtain $\chi(G)=6$ therefore requires
$\alpha(H)=5$ and one odd component. A $C_7$ leaves at most one vertex and
gives independence number at most four. A $C_5$ leaves at most three vertices,
and equality five forces all three to be isolated. Hence necessarily

$$n=8,\quad H=C_5\sqcup3K_1,\quad G=C_5\vee K_3.$$

Three consecutive vertices of the $C_5$, together with $K_3$, contain
$K_6-e$, contradicting Lemma 1. $\square$

The checker independently enumerates every path/cycle component profile
through eight vertices and verifies its complement's chromatic number by
exact backtracking. There are 46 profiles at eight vertices; the displayed
one is the sole non-five-colourable, $K_6$-free profile. No assertion about
arbitrary graphs on eight vertices is inferred from a sample.

The geometric hostile is explicit. Take the four tetrahedral points

$$ (1,1,1,0),(1,-1,-1,0),(-1,1,-1,0),(-1,-1,1,0) $$

and tips $(0,0,0,\pm1)$, with positive metric matrix
$\operatorname{diag}(1,1,1,5)$. Fourteen distances squared are eight, but
the tip distance squared is twenty. Thus $K_6-e$ really is realizable as
an equal-distance graph in four-space; the missing **upper diameter bound**
is what fails. Its threshold-at-eight graph is $K_6$, while its actual
diameter graph has just the tip edge. The production `diameter_partition`
requires the exact diameter threshold and rejects the former invocation.

This is a scoped proof, not a claimed best known cardinality theorem.
Beyond eight points, the complement of a minimal obstruction can have
degree three, so (1)'s path/cycle classification no longer applies.

## 3. Exact direction certificates and five-vector blockers

Fix a finite configuration $X\subset\mathbb R^4$ and a distance threshold
$\tau>0$. Join $x,y$ when $\|x-y\|\ge\tau$. A linear functional $v$ selects

$$B_v=\{x\in X:\langle v,x-y\rangle>0
                  \text{ for every neighbour }y\}. \tag{2}$$

Every $B_v$ is independent: two neighbours cannot both lie strictly ahead
of one another. A cover by five such classes therefore supplies a
five-colouring, with all within-class distances below $\tau$.

For a proposed class $S$, retain the **marked** collection

$$\mathcal D_S=\{(x,y;x-y):x\in S,\;xy\text{ is a threshold edge}\}. \tag{3}$$

**Theorem 3 (finite alternative).** If (3) is nonempty, exactly one holds:

1. A direction has positive product with every vector in (3), and hence
   selects all of $S$.
2. There are at most five vectors $d_1,\ldots,d_r$ from (3), with their
   original ordered endpoints, and positive weights satisfying
   $\sum_i\lambda_i=1$, $\sum_i\lambda_i d_i=0$.

*Proof.* A direction excludes a zero convex combination. Conversely, if zero
is outside the compact convex hull, let $z$ be its nearest point to zero.
Then $z\cdot d\ge\|z\|^2>0$ for every vector in the hull, giving the
direction. If zero is in the hull, Carathéodory reduces a convex representation
to at most five vectors; discard zero coefficients. $\square$

For rational coordinates both certificates can be rational. The executable
finds a feasible solution to $d_i\cdot v\ge1$ by enumerating up to four
linearly independent active normals of its minimum-norm solution. If that
fails, it enumerates affinely independent convex supports of size at most
five. Every output is independently checked by rational arithmetic. The
metric matrix is retained to authenticate which pairs are threshold edges;
the linear-functional feasibility itself uses the coordinate differences.
A positive rational Gram matrix is an exact Euclidean realization through
an invertible real linear change of coordinates, with no rank substitution.

**A sharp loss example.** Set
$u_0=(2,-1,-1)$, $u_1=(-1,2,-1)$, $u_2=(-1,-1,2)$ and
$X=\{\pm u_0,\pm u_1,\pm u_2\}$, embedded in four-space. Its diameter graph
is three disjoint antipodal edges. The three positive $u_i$ form an independent
set, but their marked difference vectors are $2u_i$ and sum to zero. Thus this
colour class cannot be selected by one direction. The entire graph still has
a two-direction cover. Directional selectability is a stronger predicate
than independence, not an equivalence silently supplied by the graph.

The five-direction route now has a precise target: cover all vertices by
five classes passing Theorem 3. Local failures have tiny exact certificates;
the theorem does not assert a globally compatible five-class choice exists.

## 4. What balanced order partitions preserve

For a finite poset $P$ on $n$ labelled elements, its order polytope is

$$\mathcal O(P)=\{z\in[0,1]^n:z_a\le z_b\text{ whenever }a<_P b\}.$$

Up to measure-zero ties, its chambers are exactly the linear extensions,
each of volume $1/n!$. Therefore

$$\operatorname{vol}\mathcal O(P)=e(P)/n!,\qquad
 \frac{\operatorname{vol}(\mathcal O(P)\cap\{z_a<z_b\})}
      {\operatorname{vol}\mathcal O(P)}
 =\frac{e(P+a<b)}{e(P)}. \tag{4}$$

This is the classical [Stanley order-polytope interface](https://math.mit.edu/~rstan/pubs/pubfiles/66.pdf),
also proved by the chamber decomposition just given. The
[1/3--2/3 conjecture](https://arxiv.org/abs/1706.04985) asks for a comparison
whose two extension counts are both at least one third whenever $P$ is not
a total order; we do not assume or prove that general conjecture here.

If every internal node of a specified decision tree has an
$\alpha$-balanced comparison, its surviving number of extensions is at most
$(1-\alpha)^h e(P)$ after $h$ comparisons. Thus a bound on comparison depth
follows once all those splitting premises have been authenticated. The
choice $\alpha=1/3$ gives the factor $2/3$; a golden-ratio interpretation
$\alpha=\varphi^{-2}=(3-\sqrt5)/2$ would give factor $\varphi^{-1}$.
That stronger universal pair-splitting assertion is already false for the
three-element poset with the single relation $0<1$: each incomparable pair
splits its three extensions as one versus two. A continuous golden cut of a
volume is a different operation from comparing two labelled elements.
This one-comparison golden-ratio interpretation is also distinct from
Peczarski's two-comparison gold condition; its exact context-weighted
inequalities are treated in the [gold transport packet](gold_partition_context_transport_20261014j.md).

**Positive exact geometric transfer.** For the five vertices of a regular
four-simplex, every projection order is realizable. More generally arbitrary
distinct scores on any affinely independent configuration can be interpolated
by an affine functional. For the regular simplex, its orthogonal permutation
symmetries additionally show that a uniformly random direction gives every
one of the $5!=120$ orders equal probability. Conditional on comparisons $P$,
equation (4) is therefore also the exact geometric direction probability.
All 120 realized orders are checked explicitly in the rational Gram model.

**The general transfer fails.** Take collinear points $x_0=0,x_1=1,x_2=2$ and
retain the comparison $x_0<x_1$. Abstract extensions are

$$012,\quad021,\quad201.$$

The next comparison $1<2$ is perfectly $1/3$--$2/3$ balanced. But the only
geometric projection order consistent with the first comparison is $012$;
the larger branch is empty in direction space. Projection scores live in the
image of a dimension-at-most-four linear map, not freely throughout
$\mathcal O(P)$. Its realizable cones and their angular measure are the
missing sidecars. An abstract order count alone cannot certify the five
geometric classes in (2).

## 5. Strict margins: the remaining infinite-set obligation

For a compact set $K$ of diameter one, an exact equivalent formulation of
a five-piece partition is

$$\exists\,0<\tau<1\quad\text{the graph on }K
       \text{ with edges at distance }\ge\tau\text{ is five-colourable}.
       \tag{5}$$

One direction chooses $\tau$ strictly above the maximum of the five part
diameters; the other bounds each part's diameter by $\tau$. By finite graph
compactness, (5) is equivalent to
$\exists\,0<\tau<1\;\forall\text{ finite }F\subset K$, the same threshold
graph on $F$ being five-colourable. The threshold is fixed before the finite
subset. A collection of finite partitions with arbitrarily vanishing margins
does not supply this quantifier order. A bounded set can first be replaced
by its compact closure, without changing its diameter.

There is an elementary robust consumer: if $F$ is an $\varepsilon$-net of
$K$, and a five-colouring of $F$ has each colour class of diameter at most
$1-3\varepsilon$, assign every point of $K$ to a nearest net point with
retained tie ownership. The resulting parts have diameter at most
$1-\varepsilon<1$. This is a sufficient certificate with the net and margin
as explicit premises. The eight-point theorem neither constructs an
arbitrarily fine eight-point net nor supplies such a uniform margin.

The concrete next obstruction is therefore twofold: larger critical diameter
configurations escape the degree-two complement classification, and finite
direction/order decisions need compatible strict margins to cover a compact
body. The accepted Hexagon asymptotic estimate attacks the latter through
many directions in large dimension, not a proof of five directions in
dimension four.

## 6. Reproduction and source contracts

[Checker](../../04-computation/experiments/borsuk_order_interfaces_20261014j.py)
and [saved JSON](borsuk_order_interfaces_20261014j.json):

```text
python -B -X utf8 04-computation/experiments/borsuk_order_interfaces_20261014j.py
python -B -O -X utf8 04-computation/experiments/borsuk_order_interfaces_20261014j.py
```

The finite production domain is two through eight distinct labelled rational
points, a positive-definite rational Gram matrix of dimension one through
four, and a positive rational squared threshold. `diameter_partition` handles
the exact diameter case. `separate`, `audit_separation`, `compile_cover` and
`audit_cover` retain the original point/metric/threshold source, edge endpoints,
strict inequalities and colour witnesses. Malformed numeric aliases, replaced
sources and altered weights are rejected. No floating-point feasibility test,
unrestricted graph oracle or external theorem is hidden inside these consumers.

Normal and optimized runs agree exactly at **1,092 checks**. The saved JSON,
normalized to LF, has SHA256
`367757403024572fb8233fa3607a67af892e3a451b445422b6c45e3a5749160c`.
Independent proof/API review and optimized replay passed. The general proofs
are separate from the declared finite controls.
