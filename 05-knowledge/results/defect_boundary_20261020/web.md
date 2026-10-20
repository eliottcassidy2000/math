# Three-edge tetrahedral graphs, signed interactions, and the lost orientation bit

20 October 2026. **PROVED below:** an all-order adjacent-edge Hessian sign
law and reconstruction theorem; elementary complement, opposition, and
switching identities. **FINITE-EXACT:** independent exhaustive controls
at orders three through six for the reconstruction, and the complete four-vertex
graph controls. No external novelty claim, Collatz proof, or scalar
enumeration closure is asserted.

The strongest recovered route is an old June question about the
Hamiltonian-path interaction graph. Its unsigned form forgets orientations.
Its signed adjacent-edge entries recover every tournament up to converse,
at every order. Only a spanning tree of the line graph and one original
edge direction are needed. This is an exact encoding, not a faster method
of computing Hamiltonian-path counts or A000568.

## 1. Inheritance and the five live concepts

The live board is **complement / incidence / cycle-cut parity / signed
interaction / stabilizer weight**.

- Closest historical question:
  [the endpoint-star floor, June S683c/d](../../../07-reflections/character-ratio-is-blind-to-dichromatic-the-OCF-governs-both-it-and-the-delta-field-s683c.md),
  Part B and its final "Honest notes", explicitly call the never-zero
  adjacent-edge coupling an observation through six vertices and request a
  proof. Only that question and its bounded data are inherited here; its
  separate spectral/dichromatic discussion is not a dependency.
- [HYP-2268's S699 reflection](../../../07-reflections/H-delta-field-gradient-hessian-and-the-7-21-holes-s699.md)
  defines the mixed flip difference. [HYP-2271's S686 reflection](../../../07-reflections/arc-interaction-graph-refines-H-s686.md)
  identifies the octahedral support at order four and asks whether an
  interaction graph reconstructs a tournament. Its warning that available
  Fourier support need not be realized nonzero is retained.
- Closest proved parity mechanism:
  [THM-2753, six-edge parity erasure](../../../01-canon/theorems/THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md)
  and [THM-2756, opposite-edge projectors](../../../01-canon/theorems/THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch.md).
  These distinguish the faithful six-edge action from its coarse cycle type
  and from the three-matching quotient.
- Least-used sidecar:
  [the August tetrahedral Haar/cycle bridge](../../../07-reflections/tetrahedral-k4-haar-xor-cycle-bridge-codex-20260816.md)
  and its [representation atlas](../../../07-reflections/tetrahedral-k4-representation-atlas-tournaments-xor-berggren-h1-codex-20260816.md).
  Their opposite-face and marked-square coordinates retain incidence and
  signed orientation that a matching histogram loses. Their physical
  ancestry-map gap remains a gap.
- Canonical current hostile:
  [the sharp weighted-defect boundary](../twin_quotient_defect_weights_20261020/SHARP_DEFECT_BOUNDARY.md).
  The independent triple blowup of P4 has defect +2. An operation fixing
  that quotient cannot be a sign-reversing pairing of those graph weights.

The targeted search found the historical never-zero question and its
computations, but no later canon proof of the adjacent sign/reconstruction
law. This records an inherited question answered below, not a claim of
priority over the mathematical literature.

## 2. All-order signed adjacent interactions reconstruct the tournament

Let T be a tournament on n>=3 labelled vertices, and let H(T) count its
directed Hamiltonian paths. For an unordered pair e, write T^e for the
tournament obtained by reversing that edge. Define

    M_ef(T) = H(T)-H(T^e)-H(T^f)+H(T^{e,f}).             (1)

Fix reference directions on the edges of K_n. Let s_e=+1 when T agrees
with the reference direction, and -1 otherwise. Let b_e be its incidence
column, with -1 at the reference tail and +1 at the reference head.
For distinct edges sharing a vertex, put

    q_ef = -<b_e,b_f>,                    q_ef in {+1,-1}.

**Theorem.** For every such adjacent pair,

    M_ef(T) = q_ef s_e s_f N_ef,          N_ef>0.        (2)

N_ef is the number of Hamiltonian vertex orders that contain both edges
and agree with T on all the other consecutive edges. The directions on e
and f themselves are not imposed in this definition.

Intrinsically, M_ef is positive when the two actual arcs form a directed
two-edge path through their common vertex, and negative when they form a
common-head or common-tail fork.

**Proof of the formula.** Expand H as a sum of the path indicators over
all vertex orders. A summand not using both edges cancels in (1). Write
e={u,v}, f={v,w}. A path using both must contain the consecutive block
u,v,w or w,v,u. Its other edges either all agree with T, or that summand
vanishes in all four terms of (1). For a compatible block, the mixed
difference of its two direction indicators is q_ef s_e s_f. Both possible
block orders give the same sign: they reverse both reference comparisons.
Summing proves the equality with the stated nonnegative count.

**Proof of strict positivity.** Choose a Hamiltonian path
a_1,...,a_k on the remaining vertices; for k=0 either block works.
Existence needs only the elementary insertion proof for tournaments:
insert a new vertex before the first old path vertex it beats, or append
it if it beats none.

If u or w beats a_1, prepend the block whose last endpoint does so.
Otherwise a_1 beats both u and w. If every a_i beats both, append either
block. In the remaining case take the first a_j that does not beat both.
Then j>=2, its predecessor beats both endpoints, and one endpoint z beats
a_j. Insert the block ending at z between a_(j-1) and a_j. The preceding
vertex beats the other endpoint, so every required external arc exists.
This constructs at least one compatible block path. QED.

The center v's other incident directions and the edge uw play no role in
N_ef. Positivity did not use a signed Fourier sum or a cancellation
estimate. A compatible block path need not leave a Hamiltonian path after
the block is deleted: the new adjacency across the gap may be backward.
Thus the existence construction is not an exact enumeration by remaining
Hamiltonian paths.

**Reconstruction corollary.** Let G=L(K_n) have one vertex for each
unordered edge of K_n, with two vertices adjacent when their edges meet.
Equation (2) gives

    q_ef sign(M_ef) = s_e s_f.                           (3)

G is connected. Along any chosen spanning tree of G, the products in (3)
recover every s_e after one root sign is specified. Without that sign,
there are exactly two solutions, s and -s, corresponding to T and T^op.
Consequently:

- adjacent Hessian signs with their incidence addresses determine every
  labelled tournament up to converse;
- only binom(n,2)-1 such product bits are needed;
- one original edge direction supplies the last bit and gives a lossless
  encoding of all binom(n,2) tournament bits.

The number of bits has not decreased. The useful change is a local
consistency carrier: every cycle in the line graph must have product +1
after dividing out q. It supports reconstruction and error detection.
The expensive quantities H(T^e) are not thereby easier to compute.

For n>=5 the full interaction support can include disjoint edge pairs.
Equation (3) uses the specified incidence subgraph L(K_n); it does not
identify that subgraph by nonzero support alone. The checker exhausts all
1024 tournaments at n=5 and all 32768 at n=6, verifies the adjacent sign
law, and confirms that the full labelled Hessian fibers are exactly
converse pairs. At n=5 it independently counts every formal block to check
the coefficient N_ef, not merely its sign.

## 3. Why four vertices recover their own incidence sidecar

At n=4, H(T)=1+2c_3(T). Equivalently a direct expansion over its four
triples gives

    M_ef(T) = -2<b_e,b_f>s_e s_f.                       (4)

Thus every shared-endpoint pair has magnitude two, and every disjoint
pair has value zero. The unsigned interaction support is the same
octahedron L(K4) for all 64 labelled tournaments: it contains no
orientation information.

The octahedron has eight triangular cliques: four are vertex stars of
K4, and four are faces of K4. For a clique with edge-vertices e,f,g, the
product of the three entries of M/2 is

    -1 for a vertex star;          +1 for a face.        (5)

Indeed the factors s_e s_f s_f s_g s_g s_e cancel. For a star, each
incidence sign at its common vertex occurs twice, leaving (-1)^3.
For a triangular face the two incidences at each of its three vertices
have opposite signs, supplying the additional (-1)^3.

The unsigned octahedron has 48 automorphisms. Its three pairs of opposite
vertices may be permuted and flipped independently, giving 3!*2^3=48.
Only 24 preserve the four vertex stars; these are precisely the
permutations of the four original vertices. Opposite-edge duality swaps
the star and face families.

Equation (5) recovers the star family from the signed Hessian. Hence any
edge permutation identifying two signed Hessians at order four must come
from a vertex permutation. Together with (3), this proves that the
unaddressed signed Hessian determines the unlabelled tournament up to
converse at order four. A single orientation bit then restores the
labelled lift once the vertex identification has been chosen.

This supplies a concrete answer to the old interaction-graph question:
unsigned support is too coarse already at four, while signed local
interactions retain the needed orientation products.

## 4. The three graphs and three different meanings of duality

There are exactly three unlabelled simple graphs with four vertices and
three edges: P4, the claw K_(1,3), and K3 disjoint union K1. Their labelled
multiplicities are 12,4,4. The degree sequences distinguish all three.

Let E=E(K4), C(D)=E\D, and J replace every edge {i,j} by the opposite
pair V\{i,j}. C is graph complementation; J is tetrahedral opposition.
They are different involutions. C exchanges the claw and triangle plus
isolate, and preserves the P4 type. J also exchanges stars and opposite
triangles, but need not equal C on an individual P4.

The line graphs of the claw and the triangle are both K3: their three
edges meet pairwise. The isolated vertex disappears. Thus line-graph
incidence alone cannot distinguish a cut star from a triangular cycle.
No general line-graph reconstruction theorem is needed for this explicit
ambiguity. Rank, a vertex-star cover, or the signed clique distinction (5)
restores the missing coordinate in this small example.

For the graphic matroid of K4, J exchanges circuits and cocircuits:
triangular faces go to three-edge stars, and four-cycles go to the
corresponding four-edge cuts. It identifies this matroid with its dual.
There are 16 spanning trees, of which 12 are paths and four are stars.
Only the paths have spanning-tree complements, giving six unordered
complementary tree pairs. But J(C(B)) is a spanning tree for every
spanning tree B. Matroid self-duality uses J; graph complementation alone
does not preserve bases.

This makes the limitation of
[THM-3260's bispanning atlas](../../../01-canon/theorems/THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary.md)
concrete: a tetrahedral planar duality is available here, whereas that
nonplanar carrier needs an explicitly chosen tree-pair chart. Equal
dimensions do not supply the missing chart.

## 5. An exact switching law on every odd tetrahedral support

Regard D as a binary vector on the six edges. Write d(D) for its four
degree parities, h_v(D) for the edge parity on the face opposite v, and
m=|D| mod 2. Directly,

    h(D)=d(D)+m*1,              h(JD)=d(D).              (6)

A face is the complement, within E, of a vertex star. Opposition carries
that face to the star, which proves both equations.

Triangle parity identifies binary edge fields modulo cuts. Its kernel
is exactly the cut space: if all triangular sums vanish, choose root 0
and write D_ij=D_0i+D_0j. Therefore when m=1,

    h(C(D))=h(JD),
    C(D) xor J(D) is a cut.                            (7)

More specifically, the difference is zero or a balanced four-cycle.
For a three-edge support, h has weight zero on a claw, four on a
triangle plus isolate, and two on a P4. The apparent three-way family
therefore becomes a complete local parity classification.

There is a characteristic-two boundary. The cut and cycle spaces of
K4 each have dimension three, but their intersection has dimension two:
its nonzero elements are the three balanced four-cycles. They do not
form a direct sum over F2. J exchanges the two spaces, without producing
a binary orthogonal decomposition.

Over the reals, the signed Hodge star on exterior two-vectors instead
does exchange the orthogonal cut and cycle spaces. With an oriented
four-vertex basis it obeys

    *rho(g)=sgn(g)rho(g)*.                              (8)

Unsigned J commutes with S4; signed * has the orientation twist (8).
They must not be substituted for one another. This is consistent with
THM-2756's different, unsigned opposition eigenspaces and its index-eight
integral gluing defect.

For tournament defect bits relative to a reference order, vertex
relabelling is affine: D -> P_g D xor epsilon_g. Equation (7) is a
fixed-gauge binary switching statement, not a tournament-isomorphism
criterion. The gauge correction in MISTAKES.md dated 18 October explicitly
rules out treating tournament bit relabelling as a pure permutation.

## 6. Complement transports the sign; it does not reverse every defect

For any graph X and g in Aut(X)=Aut(complement(X)), let chi_X(g) be its
edge-reversal character relative to an arbitrary reference orientation.
The edges of X and its complement partition those of K_n, so

    chi_X(g) chi_complement(X)(g)=sgn(g).                (9)

The right side counts all inversions of g. Thus complementation twists
the character by the vertex sign, and flips it only on odd permutations.
For the three-edge trio:

| Graph | Automorphism | Vertex sign | Edge-reversal sign |
|---|---|---:|---:|
| claw | a leaf transposition | -1 | +1 |
| triangle plus isolate | the corresponding transposition | -1 | -1 |
| P4 | reversal, a double transposition | +1 | -1 |

P4 and its complement are both sign-odd. This is the first immediate
hostile to calling graph complementation a uniform sign-reversing rule.
For every odd-order permutation, sgn(g)=+1, so complement preserves the
edge character on that sector. It does not remove the missing odd-cycle
terms of a fixed-point-free tournament/core cycle-index boundary.

There is a second loss. Complementation changes independent false twins
into adjacent true twins. It leaves the independent-triple packing
carrier used by the current defect formula. For X=P4 with each vertex
replaced by three independent vertices:

- its defect is +2 by the current weighted formula;
- complementing the four-vertex quotient and again taking independent
  triples produces an isomorphic graph, still with defect +2;
- taking the actual graph complement produces clique modules and no
  independent twin triple, so its defect is zero.

Neither version pairs +2 with -2. The claw and triangle-plus-isolate
independent triple blowups both have zero defect; the claw's three leaf
blocks merge into one false-twin class of size nine. Treating its four
drawn blocks as four canonical quotient vertices would be another error.

## 7. Precise continuation targets and reproduction

The positive next construction is an incidence-addressed signed local
compiler: retain a spanning tree of L(K_n), the q-corrected Hessian signs,
and one edge orientation. It is exact and compatible with local defect
edits. For example, flipping s_e changes exactly the products incident to
the line-graph vertex e. The remaining question is whether the defect
family provides these products cheaply from its structural parameters.
Re-evaluating full Hamiltonian-path counts to obtain them is not a
complexity gain.

For the decorated defect cancellation, the useful new constraint is
(9), together with the false/true-twin swap. A candidate that uses
complementation must transport a packing stabilizer and its vertex-sign
twist, and must represent clique modules as well as independent modules.
The order-twelve P4 blowup is the first mandatory control. No
sign-reversing involution is supplied here.

The independent standard-library checker is [web_check.py](web_check.py);
its stored output is [web_check.out](web_check.out). It uses explicit
exceptions, includes positive and hostile controls, and is run normally
and with `-O`, requiring identical output. The finite universes are all
64 four-vertex graphs, their 264 graph/automorphism pairs, all 24 vertex
permutations, all 720 edge permutations for the octahedron audit, and all
labelled tournaments at orders three through six for the stated Hessian controls.

    python 05-knowledge/results/defect_boundary_20261020/web_check.py
    python -O 05-knowledge/results/defect_boundary_20261020/web_check.py
