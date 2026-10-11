# Simplex, sphere, and cube: polarity, marked rays, and a height boundary

**Status: PROVED elementary representation statements; FINITE-EXACT controls.**
Started October 9 and completed October 11, 2026. The polarity facts below are
classical and are proved here to fix their meaning. The certificate interface
is a marked homogeneous ray: changing its norm section preserves an
authenticated arithmetic state, its original source, and its guards.
This supplies no new ROOT coverage or uniform Collatz height bound.

## 1. Inheritance and distinct meanings

The closest proved topological mechanism is
[THM-2688, affine facet normalization, vertical sphere, diagonal simplex, and
lens quotient](../../01-canon/theorems/THM-2688-affine-facet-normalization-vertical-sphere-diagonal-simplex-and-lens-quotient.md).
Its boundary simplex becomes a sphere by radial projection from the barycenter.
It explicitly separates a simplicial image from an orbit quotient and from a
physical lift. We inherit this distinction, not an arithmetic orbit conclusion.

The owner's historical
[Eureka, Zeckendorf, and the Simplex–Cuboid Duality](../../07-reflections/eureka-zeckendorf-simplex-cuboid.md)
uses another meaning. Its relevant elementary identity is
\[
k^2=\binom{k}{2}+\binom{k+1}{2},\qquad
2^{k^2}=2^{\binom{k}{2}}2^{\binom{k+1}{2}}.
\]
Here “simplex” and “cuboid” refer to edge/bit sets and their counts, not convex
polarity. The path independence polynomial at fugacities 1 and 2 changes a
counting weight, not a polar body. Other historical tournament assertions in
that reflection are not dependencies of this note.

[Small port refinement](small_port_refinement_20261004.md) already retains
ports and frames under changes of chart.
[Tetrahedral pair refinement](tetrahedral_pair_refinement_20261004.md) retains
lattice congruences despite an invertible Hadamard transform.
The corrected absolute-coordinate adjacency rule in
[edgewise cumulative coordinates](edgewise_cumulative_coordinates_20261004.md)
is a hostile precedent: an invertible coordinate map does not by itself
authorize new transitions.

The concept board is polarity, radial sections, marked affine guards, lattice
membership, compactification, and grounded terminals.

## 2. Classical polarity and dimension exceptions

For a full-dimensional closed convex body $K\subset\mathbb R^d$ with the
origin in its interior, fix the Euclidean pairing and set
\[
K^\circ=\{y:\langle x,y\rangle\le1\text{ for every }x\in K\}.
\]
The pairing, origin, and scale are part of this definition.

| Body | Polar in this pairing | Preserved sense of shape |
|---|---|---|
| Euclidean unit ball | The same unit ball | Literal equality |
| Centered regular simplex | A reflected and rescaled simplex | Similarity and combinatorial duality |
| Cuboid with half-widths $a_i$ | $\sum_i a_i|y_i|\le1$ | A weighted crosspolytope |
| Ellipsoid $x^TQx\le1$, $Q>0$ | $y^TQ^{-1}y\le1$ | The same family under linear changes of coordinates |

For the simplex, in $H=\{x\in\mathbb R^{d+1}:\sum x_i=0\}$ let
\[
v_i=e_i-\frac1{d+1}\mathbf1,\qquad \Delta=\operatorname{conv}\{v_i\}.
\]
Then
\[
\Delta=\{x\in H:x_i\ge-1/(d+1)\},\qquad
\Delta^\circ=-(d+1)\Delta.
\]
Indeed $a_i=-(d+1)v_i$ satisfies
$\langle a_i,x\rangle=-(d+1)x_i$ on $H$, giving precisely the facets.
Moreover $\langle a_i,v_j\rangle=1$ for $i\ne j$ and $-d$ for $i=j$.
If the vertices instead have unit length, the polar scale is $-d$.

For the cuboid, signs can be optimized independently:
\[
\max_{|x_i|\le a_i}\langle x,y\rangle=\sum_i a_i|y_i|.
\]
Its polar vertices are $\pm e_i/a_i$. A $d$-cube has $2^d$ vertices, whereas
its polar has $2d$; these differ for $d\ge3$, excluding combinatorial
selfduality. A self-similar product subdivision of the cube is a different
operation.

In dimension one all three body families are intervals up to scale, with
$[-a,a]^\circ=[-1/a,1/a]$. In dimension two the square polar is a rotated,
rescaled square: $(x,y)\mapsto((x-y)/2,(x+y)/2)$ sends square vertices to
diamond vertices. A general rectangle instead has a rhombus as polar.
For ellipsoids the displayed formula follows from Cauchy–Schwarz after
applying $Q^{1/2}$; a spheroid has reciprocal polar semiaxes.

**Literal selfpolarity lemma.** In the same Euclidean pairing, $K=K^\circ$
implies $K$ is the Euclidean unit ball.
For $x\in K=K^\circ$, $\langle x,x\rangle\le1$, hence $K\subset B_2$.
Reversing inclusion under polarity gives $B_2\subset K^\circ=K$.
Thus simplex selfduality allows a reflected/scaled identification rather
than literal equality.

Classical definitions and vertex/halfspace descriptions also appear in
Jesús A. De Loera's
[Polyhedral Convexity class notes, §1.1](https://www.math.ucdavis.edu/~deloera/TEACHING/MATH168/notesconvexity.pdf).
The polarity identities needed here are proved directly above.

## 3. Boundary, filled body, and three ray sections

A triangle boundary is homeomorphic to a circle; a filled triangle is
homeomorphic to a closed disk. A triangle is not affinely equivalent to a
disk, since affine maps preserve boundary segments. This is the
boundary/body distinction behind the radial map in THM-2688.

For $z\in\mathbb R_{\ge0}^3\setminus\{0\}$ and
$g\in\{\|\cdot\|_1,\|\cdot\|_2,\|\cdot\|_\infty\}$, put $N_g(z)=z/g(z)$.
The three sections are
\[
\Sigma_1=\{q\ge0:\sum q_i=1\},\quad
\Sigma_2=\{q\ge0:\sum q_i^2=1\},\quad
\Sigma_\infty=\{q\ge0:\max q_i=1\}.
\]
They are a filled triangle, a spherical octant, and three square patches
on the positive cube boundary. The octant is a curved triangle, not the
entire sphere. Each ray meets each section once, so
\[
N_h(N_g(z))=N_h(z).
\tag{1}
\]
Section changes are mutually inverse homeomorphisms. Zero-coordinate face
labels are preserved. On the cube, retain every coordinate attaining the
maximum: a corner has several incident facet labels.

For rational rays the script stores a section as a primitive nonnegative
integer vector divided by $\sqrt D$, where $D$ is the squared sum, sum of
squares, or squared maximum. Thus $(3,5,1)/\sqrt{35}$ is represented
exactly, without pretending it is a rational lattice point. The API covers
rational rays; the topological statement covers all real rays in this cone.

## 4. Marked source, composition, and native guards

For positive integers $x$ (current) and $n$ (immutable original source),
encode $z=(x,n,1)$. The marked last coordinate of each section is positive,
and
\[
x=q_1/q_3,\qquad n=q_2/q_3.
\tag{2}
\]
The coordinate order and marker convention make this map lossless.
An affine guard $ax+bn+c\ge0$ has exactly the sign of
$a q_1+b q_2+c q_3$. Original-source payment $x<n$ becomes $q_1<q_2$.

For a valuation word $w=(a_1,\ldots,a_r)$ with formal map
$F_w(x)=(Px+B)/Q$, $P=3^r$, $Q=2^{\sum a_i}$, use
\[
M_w=\begin{pmatrix}P&0&B\\0&Q&0\\0&0&Q\end{pmatrix}.
\tag{3}
\]
Positive homogeneity gives
$N_g(M_wN_h(z))=N_g(M_wz)$ for a positive image. Matrix multiplication
follows the original chronological word. Changing sections or reassociating
the product therefore preserves the full affine record.

Actual arithmetic still needs authentication. At each step the implementation
decodes (2), checks $x,n$ odd and $x>1$, and checks $a=v_2(3x+1)$.
ROOT can be the final result; another step from ROOT is rejected. Neither
the matrix nor a geometric halfspace proves these congruence and first-hit
guards by itself.

**Lost-marker hostile.** Unmarked $(3,5)$ and $(9,15)$ define the same ray
and have identical section images, but
\[
v_2(3\cdot3+1)=1,\qquad v_2(3\cdot9+1)=2.
\]
The marked triples $(3,5,1)$ and $(9,15,1)$ are distinct. Dropping the unit
marker destroys even the next branch. Keeping it repairs exact decoding.
After the supplied ROOT words $(1,1,10)$ and $(1,2,8)$, the interface also
keeps $(1,151,1)$ distinct from $(1,75,1)$.

Source, target, map, predicate, loss, and repair are now explicit: ordinary
source monitors map to norm sections; exact affine guards and authenticated
words survive; unmarked rays lose scale and integer branches; the fixed unit
marker repairs that loss. No new route to ROOT follows.

## 5. Compactification and a quantitative height threshold

Write $m_g$ for the last coordinate of $N_g(x,n,1)$ and $M=\max(x,n)$.
For $x,n\ge1$,
\[
m_1=\frac1{x+n+1},\quad
m_2=\frac1{\sqrt{x^2+n^2+1}},\quad m_\infty=\frac1M,
\tag{4}
\]
and
\[
\frac1{2M+1}\le m_1\le m_2\le m_\infty=\frac1M,\qquad
\frac1{\sqrt{2M^2+1}}\le m_2.
\tag{5}
\]
Consequently $m_g\to0$ if and only if $M\to\infty$. An unbounded sequence
has an escaping subsequence; it need not escape at every index.
Compactness then gives subsequential limits on the last-coordinate-zero
boundary. Such points are directions at infinity for this chart, not finite
ordinary integer monitors.

For the explicit integer sequence $(x_j,n)=(2^j+1,3)$ every finite marker
is positive but the markers tend to zero. This is a geometric hostile to
uniform boundedness, not a claim that this sequence is a Collatz orbit.
The chart condition $q_3>0$ does not supply a uniform positive floor.

**Quantitative marker theorem.** If an independently justified rational
$0<\eta\le1$ satisfies $m_g^2\ge\eta$, then
\[
\max(x,n)\le B(\eta):=
\operatorname{isqrt}\!\left(\left\lfloor1/\eta\right\rfloor\right).
\tag{6}
\]
This follows from $m_g\le1/M$. The common bound is sharp for the cube.
The implementation uses rational square comparisons and integer square root.
For positive integer $B$ there are also the exact gauge-specific equivalences
\[
m_1\ge1/(B+1)\iff x+n\le B,\qquad
m_2^2\ge1/(B+1)\iff x^2+n^2\le B,\qquad
m_\infty\ge1/B\iff M\le B.
\]
These are conditional state bounds, not assertions that an orbit satisfies
their premises.

**Conditional finite-core consequence.** If the same floor $\eta$ has been
proved at every preterminal state of one actual trajectory, every visited
integer is at most $B(\eta)$. If the finite directed graph on these
admissible states is additionally authenticated and has no reachable
non-ROOT cycle, the trajectory reaches ROOT. An infinite deterministic
trajectory in that box would repeat a state. A height bound alone forces
eventual periodicity, not the distinguished terminal. Both a whole-orbit
floor and the finite cycle/grounding check remain separate obligations.
The production function marked_height_bound checks a floor at one supplied
monitor only.

## 6. Dual markers and lattice loss

For invertible real $A$,
\[
(AK)^\circ=A^{-T}K^\circ,\qquad
\langle Ax,A^{-T}y\rangle=\langle x,y\rangle.
\tag{7}
\]
Primal/facet incidence and pairing signs survive when the dual marker is
transformed contragrediently. Standard integer-lattice membership need not:
$A=\operatorname{diag}(2,1)$ sends the first dual basis vector to $(1/2,0)$.
Preserving the same standard lattice in both directions requires an integer
unimodular map. Otherwise the transformed lattice or its congruence data
must remain in the record. A real geometric equivalence can preserve one
predicate while losing an arithmetic one.

## 7. Exact controls and scope

Script:
[simplex_sphere_cube_geometry_20261009.py](../../04-computation/experiments/simplex_sphere_cube_geometry_20261009.py).
Transcript:
[simplex_sphere_cube_geometry_20261009.out](simplex_sphere_cube_geometry_20261009.out).

Declared universes:

- Nonzero triples in $\{0,\ldots,4\}^3$ in all sections, including boundaries
  and cube-facet ties: 372 section cases.
- All odd current/source pairs from 1 through 63, actual and deliberately
  wrong valuations, and two supplied ROOT words.
- All positive current/source pairs from 1 through 31 for marker formulas
  and floors: 2,883 section cases; 32 escape stages in each chart.
- Simplices of dimensions 1 through 6 with 19,530 rational facet-membership
  cases; cuboid polar support in dimensions 2 through 6; square and cube
  dimension controls.
- Irrational-coordinate and nonunimodular controls, plus 17 malformed
  type/guard/marker inputs.

Reproduction:

~~~text
python -B -X utf8 04-computation/experiments/simplex_sphere_cube_geometry_20261009.py
python -B -O -X utf8 04-computation/experiments/simplex_sphere_cube_geometry_20261009.py
~~~

Both runs perform **43,642 exact checks**. These supplement the proofs;
they establish no uniform orbit-height bound by sampling, no new ROOT
coverage, and no identification of the historical meanings of duality.
