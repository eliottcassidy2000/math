# Fivefold characters, retained geometry and graph boundary guards

**PROVED:** the elementary weighted-character, lattice-window and graph
interfaces below. **FINITE-EXACT:** the declared rational/integer controls.
Classical finite-conic and Platonic correspondences are **CITED / inherited**
from the corrected sources, not newly claimed classifications. No result here
establishes Collatz convergence, six-colour equality, or a new tiling theorem.

The useful links are maps with a retained coordinate. Golden conjugation
selects four exact weight differences. A four-dimensional fivefold lattice
has a planar projection whose missing conjugate coordinate controls
discreteness. Petersen contraction preserves integer circulations with an
explicit lift, whereas preserving electrical energy requires a different,
weighted reduction. The common integer five alone supplies none of these maps.

## 1. Inheritance and correction boundaries

The current source order for this recovery is:

- [105 research atlas, section3](105_research_atlas_20261018f.md): the
  arithmetic character cube of $\mathbb Q(\sqrt3,\sqrt5,\sqrt7)$, its
  golden marked pairs, and the 105 optimal matching atlases.
- [Quadratic-form/conic packet](pythagorean_form_conic_graphs_rp_margin_20261018.md):
  the same form $x^2+y^2-z^2$ over the integers and finite fields; corrected
  attribution to Coxeter for the finite-conic graphs (MISTAKE-641).
- [Petersen/Pythagorean synthesis, sections3,6](graph_pythagorean_20261018_synthesis.md):
  literal graph maps and complete five-pole colouring supports;
  [THM-4116, ordered boundary gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
  is the inherited count-gluing mechanism.
- [Forest tilings](forest_20260926_tilings.md): a complete flat-star census,
  forced boundary-word obstructions, and valid periodic/fractal selections.
- [Simplex/sphere/cube geometry](simplex_sphere_cube_geometry_20261009.md):
  polarity versus marked-ray normalization, with source and height retained;
  [THM-2688](../../01-canon/theorems/THM-2688-affine-facet-normalization-vertical-sphere-diagonal-simplex-and-lens-quotient.md)
  supplies its older facet-normalization route.

The historical [five Platonic tournaments](../../07-reflections/the-five-platonic-tournaments.md)
is an idea source, not an exact five-way correspondence. For example, the
polar of the six-dimensional tournament cube is a crosspolytope with twelve
vertices, not the three-dimensional octahedron with six vertices. The exact
octahedron on the six edges of a tetrahedron is the line graph $L(K_4)$.
Duality of a polytope also does not identify reversal of all tournament arcs.
The [corrected angle-defect reflection](../../07-reflections/the-angle-defect-spine-five-solids-five-lattices-and-the-seventh-prime.md)
and [cardinality atlas](cardinality_bridges_20261018m.md) retain this distinction.

Anchor: authenticated boundary constraints. Niche: a weight-preserving
fivefold transform. Wildcard: lattice projection and electrical reduction.
The live concepts are source labels, cyclic order, conjugation, paired weight,
removed-edge flow, and energy. The hostile probes below change exactly one of
these while retaining a tempting coarser statistic.

## 2. Golden conjugation gives an exact weighted-cancellation criterion

Let $F=\mathbb Q(\sqrt3,\sqrt5,\sqrt7)$ and identify
$G=\operatorname{Gal}(F/\mathbb Q)$ with $\mathbb F_2^3$ by its signs on
the three displayed square roots. For

\[
 D=\{1,3,5,7,15,21,35,105\},\qquad
 g(\sqrt d)=\chi_d(g)\sqrt d,
\]

the eight $\chi_d$ are all real characters of $G$. Multiplication of radical
labels is modulo squares. Thus $d\mapsto5d$ means the fixed-point-free
involution with pairs

\[
 \{1,5\},\quad\{3,15\},\quad\{7,35\},\quad\{21,105\}.
\]

For arbitrary real or complex weights $w_d$, put
$f_w(g)=\sum_{d\in D}w_d\chi_d(g)$. With representatives
$R=\{1,3,7,21\}$, the exact identity is

\[
 f_w(g)=\sum_{d\in R}\chi_d(g)
             \bigl(w_d+\chi_5(g)w_{5d}\bigr).                 \tag{1}
\]

**Weighted cancellation theorem.** The observer $f_w$ vanishes at **every**
$g$ that conjugates $\varphi=(1+\sqrt5)/2$ to $1-\varphi$ if and only if
$w_d=w_{5d}$ for all four pairs.

Indeed that Galois half is $\chi_5=-1$. The restrictions of the four
$\chi_d$, $d\in R$, to this affine half form a four-by-four invertible
Walsh matrix. Equation(1) therefore vanishes there exactly when its four
coefficients $w_d-w_{5d}$ vanish. No positivity premise is needed. Cancellation
at one selected Galois element is weaker and cannot replace all four tests.

Equivalently the label involution $J$ gives rank-four projectors
$(I+J)/2$ and $(I-J)/2$. The golden-preserving half reads the pair sums;
the golden-conjugating half reads the pair differences. This is a finite
Fourier statement about a marked character, not an identification of the
character group with physical Platonic vertices.

**Hostile.** Take $w_1=2,w_5=1$ and all other weights zero. At the pure
golden flip, $f_w=1$, although the two labels are in the same golden pair.
Grouping them and assuming cancellation would discard the weight defect.

The atlas subsequently maps the eight affine characters to $H_8$ and
matching partitions. It does not make the golden pairing one of the five
Platonic solids. Swapping the radical3 and5 address bits is a code symmetry;
it is not a $\mathbb Q$-automorphism that exchanges those quadratic subfields.
The source field and the chosen golden character remain marked.

## 3. Fivefold symmetry in four dimensions: what planar projection loses

Put $\zeta=e^{2\pi i/5}$ and

\[
 A_4=\{a\in\mathbb Z^5:\sum_i a_i=0\},\qquad
 z_j(a)=\sum_{i=0}^4 a_i\zeta^{ji}\quad(j=1,2).
\]

Cyclic permutation of the five coordinates is an integral order-five
isometry. It intertwines with multiplication by $\zeta$ and $\zeta^2$ on
the two complex coordinates. The discrete Fourier identity gives

\[
 |z_1(a)|^2+|z_2(a)|^2=\frac52\sum_i a_i^2.                \tag{2}
\]

Hence $a\mapsto(z_1,z_2)$ is a scaled real-linear isometry on the
four-dimensional sum-zero space, and its integral image is a discrete
four-dimensional lattice. Projection to $z_1$ alone is injective on $A_4$:
a rational polynomial of degree at most4 vanishing at $\zeta$ is a multiple
of $\Phi_5=1+X+X^2+X^3+X^4$; the sum-zero condition kills that multiple.

Nevertheless that planar projection is **not discrete**. Here is an exact
sequence, not a dimension heuristic. Since
$\varphi=1+\zeta+\zeta^{-1}$ and $\varphi^{-1}=\varphi-1$ are integral
units, every

\[
 z_n=(1-\zeta)\varphi^{-n}\ne0
\]

is in the image of $A_4$. To check the last assertion, represent it by a
degree-at-most3 integral polynomial $b$. Its coefficient sum is divisible
by5, because evaluation at1 modulo $\Phi_5(1)=5$ preserves divisibility
by $1-\zeta$. Adding a multiple of $\Phi_5$ gives the unique sum-zero
five-coefficient source. Physically $z_n\to0$. In the other embedding,
$\varphi$ becomes $1-\varphi=-\varphi^{-1}$, so

\[
 |z_{2,n}|=|1-\zeta^2|\varphi^n\longrightarrow\infty.
\]

The link to section2 is restriction to the same quadratic subfield
$\mathbb Q(\sqrt5)$: $\zeta\mapsto\zeta^2$ restricts to golden conjugation.
The full field $\mathbb Q(\zeta_5)$ has cyclic Galois group of order4,
whereas the earlier multiquadratic field has group $\mathbb F_2^3$.
Only the specified restriction is shared; their entire actions are not
identified.

**Positive surviving target.** For every bounded internal window
$W\subset\mathbb C$, the selected physical set
$\{z_1(a):a\in A_4,\ z_2(a)\in W\}$ is uniformly discrete. Differences
with physical length at most1 have bounded internal coordinate in $W-W$;
equation(2) leaves only finitely many lattice differences. Their nonzero
physical lengths have a positive minimum by injectivity. A window, edges,
and covering conditions still have to be supplied before this is a tiling.
No Penrose classification or relative-density theorem is inferred here.

The dimension distinction is exact. Any nonidentity order-five integer
matrix has $\Phi_5$ dividing its characteristic polynomial, hence dimension
at least4. Irreducibility follows, for example, by applying Eisenstein at5
to $\Phi_5(X+1)$. Thus a periodic lattice in dimensions2 or3 cannot admit
such an order-five linear symmetry, while $A_4$ supplies it in dimension4.
This does not exclude finite icosahedral solids or nonperiodic patterns in
ordinary three-space.

## 4. Petersen to K5: what does and does not survive contraction

Use literal vertices $0,\ldots,9$, outer edges $i(i+1)$ modulo5, inner edges
$(5+i)(5+i+2)$ modulo5, and spokes
$M=\{i(5+i):0\le i<5\}$. Contracting these five disjoint forest edges by
$\pi(v)=v\bmod5$ gives $K_5$: each of its ten edges has one non-spoke lift.
The JSON retains every edge, orientation, and the six basis cycles.

### Integral circulations: an exact isomorphism with a compulsory lift

Orient each edge from its smaller label to its larger label, with divergence
positive at the tail. Given a $K_5$ circulation $x$, pull each edge flow back
with its induced orientation sign. Let $d_i$ be the resulting divergence at
outer vertex $i$. The unique missing spoke flow is

\[
 j_{i,5+i}=-d_i.                                           \tag{3}
\]

The quotient conservation equation says the inner divergence is $-d_i$,
so (3) cancels both. Restriction inverts this construction, and integrality
is preserved in both directions. This proves an isomorphism of the two
rank-six integral cycle lattices; equality of ranks alone would not prove it.
More generally, contraction of a forest preserves integer circulations by
the unique tree-flow solution to the retained vertex divergences.

For unit edge resistances the retained source energy is

\[
 \mathcal E_P(j)=\sum_{e\in E(K_5)}x_e^2+\sum_{i=0}^4d_i^2. \tag{4}
\]

The last term is the removed-edge bill. A quotient triangle has energy3;
its lifted source cycle has energy5. The receipt records both six-by-six
Gram matrices, rather than silently identifying them. Zero-cost contraction
discards the second term; it is not a physical unit-resistance reduction.

### Electrical boundary response: the correct Schur complement

Keep the five outer potentials $v$, and minimize over the five inner ones.
With $L_o,L_i$ the outer/inner cycle Laplacians, the source block Laplacian is

\[
 \begin{pmatrix}L_o+I&-I\\-I&L_i+I\end{pmatrix}.
\]

The minimizing inner potential is $(L_i+I)^{-1}v$. Thus the exact response is

\[
 L_{\rm eff}=L_o+I-(L_i+I)^{-1}.
\]

This is a complete graph with conductance $12/11$ on the outer pentagon
and $2/11$ on its five diagonals, diagonal entry $28/11$. It is not the
unit-conductance $K_5$ obtained by graph contraction. For example, resistance
between original outer vertices0 and1 is $3/5$, while in the contracted
unit $K_5$ it is $2/5$. The Schur response retains the former value.

### Proper vertex colouring: the removed-edge predicate changes

A proper colouring of $K_5$ pulls back to a proper colouring of $P-M$
that is constant on each spoke pair. Conversely every such colouring
descends. It cannot be a proper colouring of $P$, because each spoke was
an actual edge. With five named colours there are120 compatible colourings
of $P-M$, and zero proper colourings of $P$ satisfying the same equality
predicate. This is the exact guard needed before using a minor as a
colour-extension interface.

## 5. Petersen to K3,3, and its genuine fivefold gluing law

Delete vertex0 and suppress its three resulting degree-two neighbours.
The six retained branch vertices are $\{2,3,6,7,8,9\}$. Besides the six
direct edges, retain the paths

\[
 2-1-6,\qquad9-4-3,\qquad8-5-7.
\]

Suppression gives $K_{3,3}$ with shores $\{2,8,9\}$ and $\{3,6,7\}$.
The path lengths are part of the interface. For two colours, a path
requires equal endpoints at even length and different endpoints at odd
length. Therefore $K_{3,3}$ has two proper two-colourings, while this actual
subdivision has none. For three colours, each length-two path allows any
endpoint pair, with two fillings for equal endpoints and one otherwise.
The actual branch support has66 named assignments; direct $K_{3,3}$ has42.
Suppression discards24 supported assignments. These are complete finite
support counts, not only a comparison of chromatic numbers.

For cubic **edge** colourings a different, valid fivefold compression is
already proved in the [Petersen synthesis](graph_pythagorean_20261018_synthesis.md#6-petersens-actual-isosceles-triangle).
Cut its five spokes, leaving two pentagon poles on the same five marked
boundary labels. Each pole has30 supported boundary words, each with one
extension. Its word has colour multiplicities $(3,1,1)$; the two singleton
positions form an edge of that pole's cyclic order $C$. Consequently

\[
 \#\text{named Tait colourings after gluing }C,D
       =6|E(C)\cap E(D)|.                                  \tag{5}
\]

Petersen uses complementary cycles and has zero support intersection;
the pentagonal prism uses equal cycles and has30. All twelve unoriented
five-cycles on the marked labels give count histogram
$\{0:1,12:5,18:5,30:1\}$. The sizes and local pieces agree; the gluing
permutation is decisive.

For arbitrary nonnegative extension weights, the source-preserving formula
is the full sum $\sum_b f(b)g(b)$ over matching marked boundary words.
It vanishes exactly when positive supports are disjoint. Signed weights can
cancel without empty support; changing to a signed observer requires keeping
that distinction. This is the direct connection to weighted cancellation,
whereas a Seidel radius or a shared Platonic count loses these predicates.

## 6. Which solid and tiling statements are actually recovered?

The finite-conic packet gives a particularly concrete Platonic bridge. Over
$\mathbb F_5$, polar conjugacy for $x^2+y^2-z^2$ gives ten internal points
forming Petersen and fifteen external points forming five disjoint triangles.
Its six conic points label six Petersen perfect matchings. Berggren's
integer branches act modulo5 through $\operatorname{PSL}_2(5)\cong A_5$;
the larger $\operatorname{PGL}_2(5)\cong S_5$ supplies additional graph
symmetries. These are finite-field incidence maps, not planar unit-distance
realizations. In particular the conic three-plus-three split at an internal
point is not automatically the branch-vertex labelling of the $K_{3,3}$
subdivision; the corrected packet gives the stabilizer obstruction.

For regular maps, $(p-2)(q-2)<4$ gives the five spherical Schläfli pairs;
equality gives exactly the three regular Euclidean plane tilings.
The five planar Bravais lattice types are a different classification.
The [classical rotation-group recovery](gilbreath_fermat_platonic_20260926.md#4-the-five-platonic-solids-and-the-fields-of-size-3-4-5-classical-verified)
retains the actual actions: $A_4$ on tetrahedron vertices, $S_4$ on the four
cube diagonals, and $A_5$ on icosahedral/dodecahedral axes and compounds.
The [cardinality atlas, section3](cardinality_bridges_20261018m.md#3-an-explicit-fifteen-bridge-including-its-action-boundary)
marks a sixth point to map a six-label perfect matching to an $A_5$
involution/icosahedral half-turn axis. Without that mark, the full $S_6$
point and matching actions differ.

The exact regular-polygon tiling recovery is stronger than matching counts:

- The flat-star equation gives17 unordered multisets or21 cyclic types.
  The locally flat $3.8.24$ star cannot extend because neighbours around a
  triangle would alternate8 and24 for an odd number of edges.
- Replacing separated six-triangle clusters in a triangular lattice by
  hexagons gives genuine whole-plane tilings. Periodic centre selections
  have finite symmetry-orbit counts; Sierpinski centre selection
  $Q=2Q+\{(0,0),(1,0),(0,1)\}$ has unbounded addresses and infinitely many
  vertex orbits despite the same two local stars. The tile boundaries
  themselves are not fractals.
- Product or projection operations need their own metric and covering
  checks. The four-dimensional lattice in section3 offers a valid
  fivefold carrier and a window criterion, not an automatic three-dimensional
  aperiodic solid tiling or an identification with the five Platonic solids.

Likewise simplex, sphere and cuboid have several distinct maps. Polarity
exchanges cube and crosspolytope; a simplex is self-polar up to its stated
reflection/scale, and the Euclidean unit ball is literally self-polar.
For the source-owned ray $(x,n,1)$, normalization by $\ell^1,\ell^2,\ell^\infty$
instead preserves coordinate ratios, affine guards and the payment $x<n$.
Dropping the final marker identifies sources such as $(3,5)$ and $(9,15)$
whose native arithmetic predicates differ. These are useful exact maps,
but they are not the same duality.

## 7. Consequences and the next obligations

| Source | Map | Preserved predicate | Necessary retained data |
|---|---|---|---|
| Eight radical characters with weights | Golden pair sum/difference | Vanishing on the full golden-flip half | Four signed weight defects and the marked character |
| Integral $A_4$ fivefold carrier | Two complex Fourier coordinates | Rotation and norm identity | Internal coordinate/window for planar discreteness |
| Petersen circulation | Five-spoke contraction plus unique lift | Integral conservation | Spoke flows and the energy correction |
| Petersen electrical network | Interior Schur complement | Outer Dirichlet energy/response | Rational edge conductances |
| Two pentagon poles | Ordered boundary inner product | Exact colouring count/support | Cyclic labels and the gluing permutation |
| Flat local polygon star | Compatible boundary construction | Actual geometric extension | Neighbour word, not merely angle sum |

These mechanisms suggest concrete tasks rather than global conclusions:
check pair-weight equality in an actual arithmetic observer; carry the full
energy matrix through a proposed graph compression; or compute the exact
boundary relation of a new geometric piece before gluing. The cycle rank6
here is not identified with another flow space merely because dimensions
agree. For example the separate weighted CRT $K_{3,7}$ construction has
its own marked row/column incidence and rank12 before further residue
layers; its source map must be supplied independently.

## 8. Exact reproduction

[Script](../../04-computation/experiments/fivefold_boundary_20261018g.py) and
[saved receipt](fivefold_boundary_20261018g.json) use only the Python standard
library and active exceptions, not assertions. No repository files are
written by running the checker.

```text
python -B -X utf8 04-computation/experiments/fivefold_boundary_20261018g.py
python -B -O -X utf8 04-computation/experiments/fivefold_boundary_20261018g.py
```

Normal and optimized output agree: **63,385 exact checks**; saved LF SHA256
`c2c484b9326373f5a2b3df375f3ffb1f05718bf50715ee8c56b2c9e93cdaefb9`.
The finite universes are all $3^8$ weights, all $3^6$ signed basis-cycle
combinations, $3^5$ outer voltage assignments, complete two/three-colour
branch support, all12 pentagon cyclic orders,625 bounded $A_4$ source
vectors, and13 exact shrinking projections. The displayed theorems have
proofs independent of these finite ranges. Typed source tuples reject
boolean/float aliases; non-circulating quotient flows are rejected.

**Independent peer audit PASS.** The reviewer checked the complete proof and
core API, including the golden-half iff, source-image/window argument,
integral forest lift and energy bill, Schur response, and the two colouring
quotient boundaries. Its independent optimized serialization exactly matches
the saved63,385-check receipt and hash above. The three owned files are frozen.
