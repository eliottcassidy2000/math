# Native seam cycles, a Hodge determinant, and the unit-gap selector

**Status: PROVED.** The results are finite combinatorial analogues of
periods, regulators and orders of vanishing. They make no claim about
algebraic Hodge classes, elliptic-curve L-functions or BSD. They do
not prove that the unit-gap selector works at every count.

## 1. The native cubical complex

Label n vertices, and coordinatize the tournament flip cube by the
orientations of the binom(n,2) arcs. Let X_n be the union of the
coordinate faces on which there is a **single directed Hamiltonian
spine valid everywhere on that face**. In particular it contains
every labelled tournament by Rédei. Its edges are exactly single-arc
flips with a common Hamiltonian path before and after. This is the
common-spine seam geometry related to THM-5271. It is **not** the
tile-complement-line metagraph of THM-4739; that theorem's Hodge
decomposition is on yet another carrier (the arc flow on K_n).

**Which single flips are missing?** Precisely the arcs between two
consecutive singleton strong components in the ordered strong
decomposition. Indeed absence means that every Hamiltonian path
uses the arc: any path avoiding it survives the change. Within a
nonsingleton strong component, cyclically cutting a spanning
directed cycle supplies a path omitting any prescribed arc, and
supplies a path with any starting vertex or any ending vertex.
(Existence follows by maximal-cycle insertion: any outside vertex
not insertable is uniform to the cycle; strongness forces an arc
from a vertex losing to the cycle to one beating it, and those two
can be inserted consecutively.) Hamiltonian paths of the full
tournament concatenate component paths, so a cross arc is
unavoidable only at a boundary of two adjacent singleton
components. In particular **every** single flip from a strong
tournament is native, whether or not its other endpoint is strong.

**Theorem.** X_n is d-connected whenever n≥4d+5. Its ordinary
isomorphism quotient and its converse-merged quotient are simply
connected for n≥9. For either quotient, reduced homology vanishes
through dimension d over any coefficient ring where the group order
is invertible when n≥4d+5.

**Proof, native complex.** For each directed path word P on all labels
take the cube face where P is forced. Finite intersections are empty
or a single coordinate face, so their nerve is homotopy equivalent to
X_n by the nerve lemma. The nerve is flag: pairwise compatibility
means no variable is forced to opposite values, and so guarantees
simultaneous compatibility.

Given at most r path words, prohibit the reverses of their edges
in the complete bidirected graph. Each vertex loses at most r
outgoing and r incoming arcs. If n≥2r+1, the total degree is at least
n−1; [Ghouila-Houri, Theorem 6](https://www.numdam.org/item/ASENS_1964_3_81_3_267_0.pdf)
gives a Hamiltonian path. This new spine is in a common closed star
with all the given words. Thus the flag nerve is r-conic. The
[(2d+2)-conic connectivity theorem, Theorem 4](https://link.springer.com/article/10.1007/s41468-023-00134-9)
gives the result with r=2d+2. In degree one this also has an
elementary proof: route each edge in a loop via a vertex adjacent
to the edge's endpoints and to a fixed base vertex, applying a
triangle homotopy. For each original loop vertex, its two neighboring
new vertices and the base form a quadrilateral with it; cone the
quadrilateral from a common star of those four vertices. This fills
the loop. A purely elementary bound n≥12 follows because a longest
allowed path has at least n−4 vertices (otherwise an outside vertex
can be appended), whereas an absent vertex can block at most eight
of its insertion positions.

**Proof, quotients.** Subdivide faces barycentrically so group
elements preserving a simplex preserve it pointwise. Any quotient
edge walk lifts from a chosen representative vertex, by translating
an edge representative to the current lift. A lifted closed quotient
walk goes from x to g x. A generator h of the quotient action
having a fixed point p can be joined by a path from x to h x whose
image retraces: take a path x→p, then follow the reverse of its
h-image. Concatenate translates of these paths for a word in the
generators. Then a given lifted walk x→g x and this constructed walk
have the same endpoints, so are homotopic relative endpoints in the
simply connected native complex. The original quotient loop is
therefore null.

Every labeling transposition (a b) fixes the midpoint of a native
edge. Make a,b twins outside their mutual arc, placed as one block
in a directed cyclic sequence of blocks {a,b}→{c}→{d}→{a,b}.
The path a→c→d→b avoids the single free arc a,b and is a common spine;
prepend the remaining vertices as a transitive block. The
transposition exchanges the endpoints of this cube edge. When
converse is also included, converse composed with reversal of a
fixed label order fixes its transitive tournament vertex, so this
and the transpositions generate the whole action with fixed points.
Finally, chains on the subdivided quotient are coinvariants of
chains on X_n. When the group order is invertible the averaging
idempotent makes this functor exact, proving the homology assertion.
∎

In particular every integer edge cochain of the labelled native
complex whose sums vanish on all oriented native squares is already
a gradient for n≥9. This is integral exactness: integrate along a
path and use the square boundaries, not merely a real eigenvalue
calculation. In the quotient the analogous statement holds on the
subdivision, with the subdivided faces retained. It need not hold
if one forgets the faces or discards all but a chosen set of seams.

## 2. A genuine finite rank/regulator identity

Let A be any nonempty finite subcomplex selected from the native
geometry, e.g. keeping only seams within a target height window
and their available 2-cells. For a quotient work with its
subdivided cells. Give the integral cellular bases the standard
Euclidean metric. Put

```text
Δ₁ = ∂₁ᵀ∂₁ + ∂₂∂₂ᵀ,          D_A(s)=det(sI+Δ₁).
```

Then the following are exact, including torsion:

```text
ord_(s=0) D_A(s) = b₁(A),
leading coefficient = (product_i v_i) τ² t₁² D₂ / D_H.
```

Here v_i are vertex counts of the connected components of the
1-skeleton, τ is the product of their ordinary spanning-tree counts,
t₁=|Tor H₁(A;Z)|, D₂ is the Gram discriminant of the integer lattice
ker(∂₂), and D_H is the Gram discriminant of the lattice of
orthogonal harmonic projections of integral 1-cycles. An empty
lattice has discriminant one. The latter projection identifies
this period lattice with H₁(A;Z) modulo torsion.

**Proof.** Hodge decomposition gives orthogonal summands
im(∂₁ᵀ), harmonic 1-chains, and im(∂₂). Thus the zero order is
the harmonic rank, and the nonzero eigenvalue product is the
product of the two boundary-map pseudodeterminants. For any integral
matrix B on a standard Euclidean lattice, the squared singular
value product is

```text
disc(ker_Z B) · disc(im_Z B).
```

Indeed, choose a lattice basis extending a basis of the saturated
kernel. Projection of the other basis vectors to its orthogonal
complement has covolume reciprocal to that of the kernel, by
block-triangular volume. The restricted B multiplies this covolume
by its nonzero singular value product and gives the image lattice.

The integral graph cycle lattice Z=ker_Z(∂₁) has discriminant τ:
in a fundamental cycle basis, Cauchy–Binet sums squared minors;
the minor for an edge set of size rank Z is ±1 when its complement
is a spanning forest and zero otherwise. The image of ∂₁ is the
zero-component-sum lattice and has discriminant product_i v_i.
This gives the first pseudodeterminant. Write B=im_Z(∂₂) and
B_sat=(span_R B)∩Z. Then [B_sat:B]=t₁, and block-triangular volume
inside Z gives disc(Z)=disc(B_sat)·D_H. Thus
disc(B)=t₁² τ / D_H; the matrix identity gives the second
pseudodeterminant D₂·disc(B). Multiplying proves the result. ∎

For n≥13, in the labelled or either quotient geometry X, the
relative boundary map gives a rational isomorphism

```text
H₂(X,A;Q) → H₁(A;Q).
```

Indeed H₂(X;Q)=H₁(X;Q)=0 by §1. In this scope the order of the
Hodge zero thus counts independent native loop classes that fill
in the ambient geometry but do **not** already fill within the
selected seam window. The discriminant formula separates this
rank from hidden torsion and metric period size. This is a
precise finite rank/leading-term counterpart to the vocabulary
of BSD; it is not analytic continuation or the BSD conjecture.

For the literal boundary of one square the Hodge zero has order 1,
with period-lattice discriminant 4 and leading coefficient 16.
Filling its native face removes the zero, and the determinant at
zero is 64. This illustrates why graph circulations before and after
retaining native 2-cells are not the same invariant.

## 3. A determinant valuation that certifies a physical target rail

Take the native seam graph, with or without identifying vertex
isomorphisms/converse, and two anchors u,v in the same component.
Loops may be discarded; keep parallel edges if desired. Label
each edge e=xy with its nonnegative integer cost

```text
w_e = |H(x)−H(y)| / 2.
```

For z>0 give that edge conductance z^w_e. Let T(z) be the weighted
spanning-tree polynomial of this component and F_uv(z) its
weighted spanning-two-forest polynomial separating u from v.
By the tree/forest formula T/F_uv is effective conductance.
Then

```text
ord_(z=0) [T(z)/F_uv(z)] = min_(u-v paths P) max_(e in P) w_e = β.
```

**Proof.** All coefficients of both polynomials are nonnegative,
so their orders are their minimum tree and separating forest costs.
Add to a minimum separating forest an edge of an optimal bottleneck
path that crosses its cut; it costs at most β. Conversely, remove
the heaviest edge on the u–v path in a minimum spanning tree; it
costs at least β and the remaining forest separates u and v. The
two inequalities give precisely the displayed valuation. In
particular it is a physical path certificate and can also be
computed by sorting edges and finding when the anchors connect. ∎

This valuation is also the exact disappearance threshold of a
relative connector. Let A_t contain the native vertices and the
edges of cost at most t, with any compatible 2-cells. In a simply
connected ambient X, the class of an ambient u–v path in H₁(X,A_t)
is independent of the path, maps to [v]−[u] in reduced H₀(A_t),
and vanishes exactly when u,v belong to the same component of
A_t. For quotient coordinates, take the image of enabled edges
as a subcomplex after subdivision; a folded edge can only add a
dangling part and does not affect anchor connectivity.

**Chan–Pak use.** If β≤1 and the anchor counts surround an odd
target h, a path of actual flips with adjacent gaps at most two
must visit h. A valuation above one says that this particular
native unit-step route cannot certify that interval; it **does not**
say h is unattainable. No first homology or period calculation
on the ambient complex alone fills a missing integer in the
Hamiltonian-count image.

## Finite controls

The companion script enumerates the shared-spine native vertices,
edges and squares from coordinate constraints through order five;
then it uses exact F₂ cycle elimination and independent subset-DP
Hamiltonian counts. Labelled native squares span the cycle space
over F₂ at orders 4 and 5 (93 and 3817 dimensions). This is not
an integral simple-connectivity claim for these small orders.
For its displayed anchored extremes 1 and 5 at order four the
bottleneck is 1. For the extremes 1 and 15 at order five it is 2;
the exact spectrum there is 1,3,5,9,11,13,15, illustrating a real
unit-rail obstruction without inventing a global obstruction.

Exact finite replay: `python3 check_native_hodge.py`.
