# Elliptic contacts, a cyclotomic ear, and a seventeen-vertex extension theorem

2026-10-16. **PROVED:** the finite exterior extension theorem, the exact
eleven-cycle attachment, and its complete boundary-colouring relation.
**FINITE-EXACT:** the rational field certificates and the marked elliptic
test patch. **INHERITED:** the infinite elliptic contact family and the old
field's five-colour ceiling. **OPEN:** a six-colour forcing context and the
Hadwiger--Nelson six/seven alternative. No priority claim is made for the
classical field, list-colouring, or path-count ingredients.

## 1. Inheritance and the next predicate

The incoming [recursive interfaces](recursive_interfaces_20261016b.md),
§2, proves infinitely many rational points on

$$C:\quad y^2=-275r^4+150r^2-3,$$

and turns them into unit contacts. Its independently checked sample of
49 points is a tree. The two retained contact parameters are

$$ (r,y)=(1/5,8/5),\qquad (91/135,10552/3645). $$

Writing $j=\sqrt{-11}$, $k=\sqrt{-55}$, the coordinates and rotations in
that construction lie in

$$K_0=\mathbb Q(\sqrt5,\sqrt{-11})
\subset K=\mathbb Q(\zeta_5,\sqrt{-11})
\subset L=\mathbb Q(\zeta_{55}). \tag{1}$$

The [rhombus/elliptic packet](rhombus_overlap_elliptic_20261011.md),
§G2, already supplies a five-colouring of the entire unit graph of $K$,
using the conjugation-fixed prime over 19 and an explicit residue colouring.
Its dependency is the maintained
[THM-4558 residue-colouring theorem](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md),
not an inference from the field's degree. Thus rational elliptic additions,
translations and field arithmetic staying in $K_0$ cannot yield a
six-chromatic graph. They remain useful sources of exact contact geometry.

The same older note proves that two-circle intersections in the odd-degree
extension $L/K$ cannot leave $K$. We retain that premise rather than rename
it as a new theorem. The new question is operational: how many exterior
vertices are needed before the **full six-colour boundary relation** can
change, and what does an actual field-escaping attachment do to it?

The live board is rational elliptic points / actual unit contacts / field
extension / full labelled boundary support / extension multiplicity /
coupled interior cores. The canonical hostile is a new unit contact whose
addition changes only the number of colourings, not which boundary
assignments extend. The least-used sidecar is the selected old field,
including all permitted denominators.

## 2. The inherited field gate, made into an extension theorem

The fields in (1) have $[K:\mathbb Q]=8$, $[L:\mathbb Q]=40$ and
$[L:K]=5$. The quadratic subfields of $L$ are precisely
$\mathbb Q(\sqrt5),\mathbb Q(\sqrt{-11}),\mathbb Q(\sqrt{-55})$,
as recovered in the incoming note. In particular $\sqrt{-3}\notin L$.

**Lemma 1 (at most one old neighbour).** A point of $L\setminus K$ can be
at unit distance from at most one point of $K$.

*Proof.* Suppose $a,b\in K$ are distinct and $|z-a|=|z-b|=1$.
Put $t=(z-a)/(b-a)$. The norm equations imply

$$t+\bar t=1,\qquad t\bar t=1/|b-a|^2.$$

Consequently $t$ is quadratic over $K$. Its degree also divides
$[L:K]=5$, so $t,z\in K$, a contradiction. This is the older two-circle
argument with all old points quantified, not merely a chosen pair of ports.
$\square$

There are two further exact geometric constraints. Unit graphs in $L$ are
triangle-free: an equilateral triangle would put
$(1\pm\sqrt{-3})/2$ in $L$ after dividing two unit differences. Also any
two distinct plane points have at most two common unit neighbours, since
two distinct unit circles intersect in at most two points.

**Theorem 2 (full-boundary extension).** Let $A\subset K$ be any set with
a specified proper six-colouring of its unit graph. Let
$B\subset L\setminus K$ contain at most seventeen distinct points. Every
such old colouring extends to the full induced unit graph of $A\cup B$.
All old and new cross-edges are included.

*Proof.* Each $b\in B$ loses at most one available colour by Lemma 1, so
its remaining list has at least five colours. Every induced subgraph of
the unit graph on $B$ has a vertex of degree at most four. Otherwise some
nonempty induced subgraph has minimum degree $\delta\ge5$. Choose a vertex
$v$ of degree $d\ge5$. Its $d$ neighbours are independent, and each has
at least four other neighbours outside $\{v\}\cup N(v)$. This gives at
least $4d$ incidences to the second neighbourhood. Each second-neighbour
vertex can account for at most two, by the common-neighbour bound. There
are therefore at least

$$1+d+2d\ge16$$

vertices. This excludes all possible cores of size at most fifteen.

For a sixteen-vertex core all the inequalities are equalities, at every
choice of root. It is five-regular. Translate one root to zero and write
its five distinct unit neighbours as $u_1,\ldots,u_5$. Its ten second
neighbours exhaust the ten unordered pairs: each is the second unit-circle
intersection $w_{ij}=u_i+u_j$. An antipodal pair would have no second
intersection, so cannot occur. A vertex $w_{ij}$ has the two first-layer
neighbours $u_i,u_j$ and needs three further neighbours. Pairs $w_{ij},w_{ik}$
cannot be adjacent because their difference is $u_j-u_k$, and the first
layer is independent. Therefore all three disjoint-pair choices must be
adjacent to $w_{ij}$.

Put $S=\sum_i u_i$. For each four-subset, the three signed pair-partition
vectors all have norm one. Their squared norms sum to
$16-|S-u_m|^2$, where $m$ is the omitted index. Hence
$|S-u_m|^2=13$ for every $m$. Summing gives
$3|S|^2+5=65$, so $|S|^2=20$ and $S\cdot u_m=4$ for all five indices.
A fixed line meets the unit circle in at most two points, contradicting
their distinctness. The checker independently obtains the stronger formal
Gram consequence: all ten pairwise inner products would equal $3/4$, giving
a rank-five Gram matrix (determinant $1/64$), not a planar one.

The incidence pattern itself is consistent: in $\mathbb R^5$ set
$u_i=\tfrac12 e_i+\tfrac3{10}(1,1,1,1,1)$. The sixteen points
$0,u_i,u_i+u_j$ have exactly the stated unit edges, forming a triangle-free
five-regular graph with forty edges. Their first-layer Gram matrix has
diagonal one and off-diagonal $3/4$. The checker retains all coordinates
and checks every pair. This is a positive higher-dimensional control, not
a minimal-dimension theorem: the planar argument used the unique second
intersection of two circles, which need not hold in higher dimensions.

For a seventeen-vertex core, a vertex of degree six or more would require
at least $1+6+\lceil6\cdot4/2\rceil=19$ vertices by the same incidence
bound. Thus the core would be five-regular on seventeen vertices, contradicting
the even degree sum. All core sizes at most seventeen have now been excluded.
Repeatedly remove a vertex of degree at most four and reinsert in reverse
order, choosing from its list of at least five. This constructs the extension
while retaining every old colour. $\square$

This is stronger than saying that the enlarged graph happens to be
six-colourable: **every** proper old assignment survives. No count or
uniformity of extension fibres is asserted in this theorem. The old set
may be infinite; only the exterior set is bounded.

**Corollary 2a (odd-degree transfer).** Theorem 2 holds with the same bound
seventeen for any conjugation-stable field $K\subset\mathbb C$ with
$\sqrt{-3}\notin K$ and any finite odd-degree extension $M/K$, replacing
$L\setminus K$ by $M\setminus K$. The extension field itself need not be
conjugation-stable. The quadratic equation in Lemma 1 has coefficients in
$K$, so its degree divides both two and $[M:K]$; likewise a unit triangle
would give a quadratic subextension by $\sqrt{-3}$. The entire remaining
argument is planar geometry and list extension. The common ambient field
and its odd degree are hypotheses; this is not a statement about unrelated
pointwise degree labels.

The same counting proof gives a general palette version. For $q\ge2$, at
most $q(q-1)/2$ exterior vertices preserve every proper $q$-colouring.
A failed extension has a vertex-minimal list obstruction of minimum degree
at least $q-1$, and its two-neighbourhood has at least
$1+(q-1)+(q-1)(q-2)/2=1+q(q-1)/2$ vertices. This triangular number counts
incidence capacity; no tournament orientation or golden-number identification
is involved. These are sufficient cardinality bounds, not claimed sharp ones.

The strongest useful necessary condition for a failed six-colour extension
is source-sensitive. In a vertex-minimal failed exterior set $H$,

$$\deg_H(b)\ge 6-t(b),\qquad
t(b)=\begin{cases}1&\text{if }b\text{ has an old unit neighbour},\\
0&\text{otherwise}.\end{cases} \tag{2}$$

Thus a prospective coupled connector needs a five-core, and its vertices
without an old contact need degree at least six. Eighteen vertices are the
first size not excluded here; no existence at that size is asserted.

The exterior field marker is necessary. Six distinct unit points
$\gamma^i\in K$, $0\le i\le5$, where $\gamma=(5+\sqrt{-11})/6$,
have no edges between them (an edge would form an equilateral triangle with
zero). Colour them with all six labels. Adding the single point $0\in K$
blocks this particular boundary assignment. The entire graph is only a
star and is two-colourable. This is failure of boundary preservation,
not a seven-chromatic graph, and it lies outside the theorem's
$B\subset L\setminus K$ premise.

### Why rotating the old contact family does not supply coupling

If $a,b\in K$ are nonzero and $|\gamma|=1$, a unit edge between $a$ and
$\gamma b$ forces

$$\bar a b\gamma^2-(|a|^2+|b|^2-1)\gamma+a\bar b=0. \tag{3}$$

For $\gamma\in L\setminus K$ this is impossible by degree five.
Nonzero coincidences are likewise impossible. Therefore the eleven copies
$\zeta_{11}^i A$, $0\le i<11$, of any old $K$ patch meet only at their
retained common centre, if present, and have no cross-edges away from that
centre. Their boundary relation is the corresponding one-point amalgam of
the individual relations. An infinite arithmetic supply and a larger
coordinate field have added no noncentral constraint. The finite controls
retain both the disjoint-copy and the exceptional shared-centre cases.

## 3. An exact field-escaping attachment with controlled boundary relation

Set $\zeta=\zeta_{11}$ and

$$p_j=1+\zeta+\cdots+\zeta^{j-1}\quad(0\le j\le10),\qquad p_0=0.$$

These are the vertices of a regular unit eleven-gon, with $p_1=1$.
Its rotation centre $1/(1-\zeta)$ lies outside $K$, so this operation
escapes the common-$K$-centre restriction above.

**Lemma 3 (complete old-field contact classification).** For $2\le j\le10$,
$p_j\notin K$. The set of points of $K$ at unit distance from $p_j$ is

$$\{1\}\ (j=2),\qquad\{0\}\ (j=10),\qquad
\varnothing\ (3\le j\le9). \tag{4}$$

This statement quantifies over the whole field, not a sample of old points.
Here is a small exact linear certificate proof. Put
$F=\mathbb Q(\zeta)$ and $E=\mathbb Q(\sqrt{-11})=F\cap K$.
The degrees $[F:E]=5$ and $[K:E]=4$ are coprime, so these fields are
linearly disjoint over $E$. In the rational basis
$1,\zeta,\ldots,\zeta^9$ use

$$J=\sum_{i=1}^{10}\left(\frac{i}{11}\right)\zeta^i,
\qquad J^2=-11.$$

The unit equation for an old point $a$ implies the linear equation

$$-p_j\bar a-\bar p_j a+z=1-p_j\bar p_j,
\qquad z=a\bar a. \tag{5}$$

First allow the three unknowns $\bar a,a,z$ to range independently over
$E$. The six rational columns are
$(-p_j,-Jp_j,-\bar p_j,-J\bar p_j,1,J)$.
Appending the right side and its $J$ multiple gives the following ranks
over $E$; the checker obtains them by exact rational elimination and
retains annihilating row witnesses for each inconsistent case.

| $j$ | coefficient rank | augmented rank | compatible old point |
|---|---:|---:|---|
| 2 | 3 | 3 | $a=1$ |
| 3,4,5 | 3 | 4 | none |
| 6 | 2 | 3 | none |
| 7,8,9 | 3 | 4 | none |
| 10 | 3 | 3 | $a=0$ |

These ranks and solutions remain the same after extending scalars from
$E$ to $K$, because $F\otimes_E K=L$. Thus the linear impossibilities
already exclude every $a\in K$; the two surviving solutions also satisfy
conjugation and $z=a\bar a$. Membership in $E$ is separately rejected for
every $p_2,\ldots,p_{10}$ by the same exact field coordinates, so none lies
in $K$. This completes the finite algebraic proof of (4).

**Theorem 4 (eleven-cycle ear).** Let an arbitrary old configuration
$A\subset K$ contain a marked unit edge $a b$. Adjoin

$$a+(b-a)p_j,\qquad j=2,\ldots,10. \tag{6}$$

The full induced unit graph is exactly the old unit graph plus an eleven-cycle
sharing just the marked edge and its endpoints. There are no additional
unit cross-edges and no coincidences with any old point.

Indeed normalize by the unit affine map $x\mapsto(x-a)/(b-a)$ and apply
(4). Internal incidence is exactly the regular polygon cycle: equal chord
length one occurs only at adjacent vertices. All coordinates in (6) lie
in $L\setminus K$.

For any $q\ge3$, restriction from proper labelled $q$-colourings of the
enlarged graph onto old colourings is surjective, with the constant fibre

$$c_q=\frac{(q-1)^{10}-1}{q}. \tag{7}$$

With the different endpoint colours fixed, the new path has ten edges.
The adjacency matrix of the complete loopless palette graph has eigenvalues
$q-1$ on constants and $-1$ on the zero-sum space. Its off-diagonal tenth
power entry is (7). The checker also computes this by elementary path DP.
For six colours the exact factor is **1,627,604**. For two colours it is
zero: the attachment creates an odd cycle. Thus this is not a claim that
every palette is preserved.

For a nonempty old unit graph containing the marked edge, its chromatic
number after attachment is exactly $\max(\chi(G),3)$. In particular this
single operation cannot turn an old $K$ graph into a six-chromatic graph,
even though its new coordinates leave $K$.

Every chosen subset of old labels inherits the same extendable boundary
assignments before and after the operation. For a finite old graph the
number of extensions of each such assignment is multiplied by the same
factor (7). This is the exact equal-fibre mechanism also isolated in the
[gold context transport packet](gold_partition_context_transport_20261014j.md):
support and multiplicity are different observables. No transfer of the
gold conjecture or of its comparison operation is asserted.

## 4. Source-owned elliptic control and executable consumers

For each of the two rational quartic points retained in §1, the incoming
formulas give

$$\alpha=(1+J)/2,\quad s=(1+55r^2)/4,\quad
b=\alpha(s+2+(y/4)J)/6,$$
$$\delta=(55r^2-1+2r\sqrt{-55})/(1+55r^2).$$

Use the seven points consisting of $\alpha$ and
$\delta^{-1}b,b,\delta b$ for each parameter. All equations are verified
in $\mathbb Q[X]/\Phi_{55}(X)$, including the actual selected embedding's
complex conjugation. The full old graph has seven vertices and six edges;
it is a tree. Attach (6) to the edge from $\alpha$ to the first $b$.
The result has sixteen vertices and sixteen edges and is exactly
three-chromatic. Its complete six-colour count is

$$6\cdot5^6\,c_6=5^6(5^{10}-1),$$

and all old six-colour boundary assignments survive. The source is small,
but Theorems 2 and 4 are not restricted to this particular elliptic patch.

The production APIs retain exact rational coordinates in the degree-40
cyclotomic basis. `compile_exterior` authenticates old-field membership,
all unit pairs, distinctness and a four-degeneracy removal order.
`extend_boundary` consumes an actual old six-colouring and returns a
complete colouring without changing any old label. `compile_cap` retains
the marked edge and authenticates the entire induced incidence graph;
`audit_cap` rejects altered coordinates, ports or incidence data.
The implementation's explicit materialization cap is 128 total vertices;
the theorem itself has no old-source cardinality restriction.

## 5. What remains productive

Incoming commit `f2988e75c2`, received during close-out, supplies a precise
alternative carrier in [elliptic cross-edges, G7](elliptic_crossedges_20261016c.md).
For the real root $\theta\in(4/5,1)$ of
$2\theta^3+2\theta^2-2\theta-1=0$, it constructs a mixed norm-one direction
over $K'=K(\theta)$ inside $L'=L(\theta)$, with
$[L':K]=15$ and $[L':K']=5$. Corollary 2a applies both relative to the
original $K$ and relative to $K'$: the latter is conjugation-stable and
cannot acquire $\sqrt{-3}$ in its odd-degree extension of $K$.
Thus this positive arithmetic escape still needs at least eighteen exterior
vertices to delete an old six-colour assignment. No five-colour ceiling
for $K'$ is inferred. The current executable compiler retains its original
degree-forty $L$ basis; the general corollary is a mathematical transfer,
not an unimplemented generic-field API claim. This corollary received a
separate independent proof audit after reading the incoming construction.

The old-field elliptic family is an infinite source of geometry, not an
escape from its five-colour ceiling. Simple cyclotomic rotations leave
different copies uncoupled. The deliberate eleven-cycle attachment does
escape the field and creates a genuine odd cycle, but its complete
six-colour boundary relation is unchanged, not merely its existence count.

A candidate that changes this relation must leave these uniform-ear and
small-exterior mechanisms: it needs at least eighteen exterior vertices
with a surviving five-core and the source-sensitive degree constraints
(2). Multiple interacting ears are not automatically covered by Theorem 4:
after the first attachment, the new source is no longer entirely in $K$.
Their cross-edges must be authenticated again. A coupled connector passing
these tests remains a valid next target; neither its existence nor its
impossibility is established here.

## 6. Reproduction

[Script](../../04-computation/experiments/hn_elliptic_bridge_20261016c.py)
and [saved JSON](hn_elliptic_bridge_20261016c.json).

```text
python -B -X utf8 04-computation/experiments/hn_elliptic_bridge_20261016c.py
python -B -O -X utf8 04-computation/experiments/hn_elliptic_bridge_20261016c.py
```

Exact controls include the nine full-field contact systems, every unit pair
of the elliptic patch and its attachment, all 192 proper old three-colourings
consumed at palette six, explicit old assignments using all six labels,
labelled path DP for palettes two through seven, eleven rotated segment
copies, the retained-centre star, the old-field six-label hostile, the exact
sixteen-core Gram system and its explicit five-dimensional unit realization,
and malformed/source-replacement controls.
No floating distance comparison or unrecorded SAT verdict is used.
Normal and optimized runs agree exactly at **474 checks**. The saved JSON,
normalized to LF, has SHA256
`6bb4e9ef7fa72b14464c26485e4f33fc32f68416f1d96e85e888424afbd52262`.
Independent full proof/API review passed, including the seventeen-vertex
bound and the complete old-field contact certificate.
