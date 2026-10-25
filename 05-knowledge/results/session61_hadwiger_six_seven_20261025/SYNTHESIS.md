# Session 61: Hadwiger–Nelson, six-root K3s, and the seven-root hourglass

## 1. The common hexagon is an actual metric object

The hourglass branch arrangement is E,E±x,E±z,E±(x+z). Its six outer
*dual* labels identify with the six units of Z[omega], omega²+omega+1=0;
the central dual label is zero. There is also a primal identification:
with x=2a−b and z=−a+2b, the region |x|,|z|,|x+z|≤1 is the Voronoi
hexagon of the Eisenstein lattice, with norm (x²+xz+z²)/3.
The primal branch E=0 is the line at infinity, not the cell center.

Exactly two index-seven sublattices make the center-plus-six-neighbors
patch a complete residue system: the conjugate ideals (2−omega) and
(3+omega). Their color maps are a+2b and a+4b modulo 7. The proof covers
all eight index-seven kernels, not a finite collection of sampled shapes.
They are classical perfect one-codes in the triangular graph T.

This explains two different sevens without confusing them. T has chromatic
number 3, but its graph square T² has chromatic number 7: the seven-neighbor
flower is a clique only after adding graph-distance-two edges. To color
the entire Euclidean plane with seven colors, one instead retains the
Voronoi cells and their metric separation. If h is nearest-center spacing,
cell diameter is 2h/√3 and minimum like-color intercell separation is
h√(7/3). The strict admissible interval is √(3/7)<h<√3/2; h=3/4 works,
including arbitrary assignment of shared boundaries. These classical
coloring mechanisms do not give a lower bound of seven for the plane.

Rotations preserve a chosen norm-seven ideal; reflections exchange the
two. A single color quotient therefore loses half the marked geometric
symmetry. Modulo 7, the factorization
Q(a,b)=(a+2b)(a+4b) identifies the six roots with the six points of the
split norm-one torus {(u,u^−1):u≠0}. This is a representation-level
connection, not merely six objects on either side.

## 2. Keeping triangular faces selects a specific hourglass quotient

The strongest new-to-route bridge appears when we do not discard the
triangular two-cells. Quotient the triangular tiling by (2−omega).
The resulting simplicial torus has 7 vertices, 21 edges and 14 faces, with
face types {t,t+1,t+3} and {t,t+2,t+3}, t∈F7. Each type separately is a
Fano plane: every pair appears in exactly one of its seven triples.
The torus/Fano ingredients are classical; explicit coordinates and
independent exhaustive checks make their application here unambiguous.

Pull these Fano structures back to the fixed seven hourglass branch labels.
Neither contains an actual branch-crossing triple. Taking both conjugate
ideals and both face types gives four distinct codes, exactly the special
four-element D12 orbit among Session 60's sixteen free Fano kernels.
Consequently the construction selects the Campedelli isomorphism class
with p_g=q=0, K²=2 and full automorphism group C2×A4. The other free-code
orbit of size 12 is not produced by this particular construction.

Thus the seven-color lattice quotient supplies actual free quotient codes
of the K²=16 hourglass. The map is: marked lattice and triangles → torus
face design → binary inertia marking → free deck kernel. It is not a
map preserving plane colorability. The graph covering T→K7 is itself an
explicit hostile: its source is three-colorable despite its seven-chromatic
target. Source colorings need not descend through the chosen period lattice.

This also repairs old unqualified 'hexagon is Fano' prose. Euclidean
collinearity is not Fano incidence. The ideal, triangle type, marking and
binary code are essential added data. Central inversion exchanges the two
face designs instead of preserving either, consistent with Session 60's
no-invariant-Fano-plane result.

## 3. Six colors reconnect to genus-two pairings through a missing bit

The earlier six-color method writes each color as a ternary label and a bit.
For a fixed ternary assignment q, retain only same-label edges H_q.
The bit lift exists precisely when the constant edge cochain 1 is a
coboundary in H_q, equivalently when each component is bipartite.
The graph H_q varies with q; a fixed cycle basis of the original graph
cannot replace all these tests.

The fifteen pairings of six palette labels are exactly the fifteen
maximal isotropic subgroups of genus-two two-torsion used for Richelot
correspondences, and exactly the fifteen possible central pencils of the
thirty Fano structures. Each pencil has TWO Fano completions. Swapping one
pair's endpoints fixes the pairing but exchanges completions, so there is
no S6-equivariant rule choosing a completion from a pairing alone.
This is an explicit missing bit, not a numerical identification of fifteens.

With the three pairs ordered, frame transport takes values in F2³.
Projection to the completion bit is parity p(x)=x₁+x₂+x₃. Full frame
gluing requires the F2³-valued cycle obstruction to vanish; completion
bits test only its parity image. The forgotten obstruction space has
dimension 2b₁ of the gluing graph. A triangle with gain (1,1,0) on every
edge passes parity but fails two channels. Unordered pair transport needs
a local coefficient system and is not silently treated by this formula.

A chosen one-root-per-pair embedding makes this parity character the
restriction of the six-root K3 canonical character. Label permutations,
deck flips and graph lifts remain distinct actions. The comparison
transports their explicitly marked character response, not colorability.

## 4. Pairwise data cannot see which quotient surface is present

Across all thirty Fano kernels, the first and second moments of the seven
canonical deck characters are identical. The third moment on a triple is 1
exactly when it is a Fano line, and 0 otherwise. Summing these third moments
on the three actual branch-crossing triples recovers d=0,1,3 and hence the
inherited 16/12/2 split of labeled kernels into Campedelli, elliptic and
Enriques cases.

This supplies a geometric counterpart of the earlier palette warning:
pairwise boundary summaries can erase the obstruction. Here every observer depending only on these uniform one-/two-character
marginals misses distinctions that change the birational type of a surface.
The triple product is a valid projective character, represented by the
complementary weight-four subset; no odd-support character was inserted
into the projective sign group.

## 5. The CM connection has opposite good primes for two purposes

The central K3 quotient's associated genus-two curve is geometrically
y²=u⁶−1, with j=0 elliptic factors and the same Eisenstein arithmetic.
Its rational-model twist is not discarded. For p≥5, those elliptic factors are ordinary
when p≡1 (mod 3) and supersingular when p≡2 (mod 3). These statements concern the
elliptic factors; no unaudited assertion about reductions of every root
model or new Picard rank is needed.

For field-plane residue colorings the relevant condition runs differently:
all norm-one edge vectors are integral at a prime P iff conjugation fixes
P. Split primes such as 7 fail this condition. A norm-seven element divided
by its conjugate gives a unit-length vector with a pole at one prime.
This obstruction persists in every larger field containing Q(√−3).
For K=Q(√−3), the prime above 3 gives the familiar three-color field-plane
bound, while the displayed degree-six curve is a bad model. This prime
need not remain usable in a larger field. Geometry, arithmetic
reduction and coloring therefore require different prime sidecars.

The older Q(ζ₅₅) residue obstruction sharpens after adjoining the actual
equilateral carrier: in Q(ζ₁₆₅), the least fixed-prime residue parameter
is 659, and Hoffman–Weil forces at least 14 colors for every coloring of the
specified single-prime residue type. The field plane nevertheless has a
seven-coloring by restricting the classical Euclidean construction.
This is an obstruction to a METHOD, not a lower bound of 14, or 7, for the
original graph. Even vertex/edge-surjective reduction of the triangular
lattice modulo 2 produces K4 from a three-colorable graph.

## 6. A historical gap closes with a different proof

THM-419's all-rank-two-lattice cap survives, but its old proof contains
false similarity and shell-size assertions. The repair is structural:
three distinct unoriented equal-norm integer directions force the
normalized Gram matrix to be rational. With fewer directions the
components are isolated vertices, paths or grids. In the rational case the norm graph embeds
into an imaginary quadratic field after normalization, and a ramified
prime supplies a two- or three-color residue lift. Thus every single
positive-distance graph on a discrete rank-two Euclidean lattice is
at most three-colorable.

The formerly problematic Q=x²+xy+6y², D=36 layer has exactly three opposite
step pairs and odd girth 23; an explicit residue map to C23 gives its exact
chromatic number 3. This is a hostile control rather than another bounded
census standing in for the general theorem. Stronger constructive periodic
and field-discriminant refinements are documented in the recovery packet.
The same conjugate norm-seven ideals also reconstruct a prior 21-point,
57-edge distance graph with exact chromatic number 3. Their relative
rotation has cosine 11/14 and infinite order, yet the two lattices remain
commensurate. Irrational rotation angle and high edge density therefore
do not by themselves escape the low-chromatic field. No rotation-rank
formula or bound for arbitrary higher-rank planar additive groups follows.

## What remains open

The plane six/seven decision is not resolved. The current lower bound 6 is
an external cited theorem, not a Lean build performed here; none of this
session's new deductions depends on independently certifying that proof.
No integer-square hourglass solution, new Picard rank, or new six-color
unit-distance equality gadget is claimed.

The useful next problem is to retain the full geometric boundary relation
while exploiting the explicit three-channel lift, rather than discard it
for pairwise statistics or a finite residue quotient. This is a research
direction, not automatically running work after this session ends.
