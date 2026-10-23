# The actual path forces joint source-target refinement to be discrete

2026-10-23. **PROVED:** the all-size refinement bound, its fractional-symmetry
consequence, and the path-clock identity below. **FINITE-EXACT:** independent
controls in `verification.json`. **OPEN:** finding the correct target and
vertex alignment efficiently from compressed source records. No theorem here
classifies the unmarked tournament metagraph from its margins.

## Inheritance and scope

The closest mechanism is the existing directed-agreement certificate in
[Session25, Theorem 1](../session25_spectrum_wiring_20261021/session25-wiring/RESULT.md)
and the [ordered-graph trijection, section 6](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).
The [signed-curl decoder](../metagraph_port_rigidity_20261023/PROOF.md),
Theorems 1--3, already supplies complete port acceptance with retained path
coordinates. We do not rediscover those tests as a new classification.

The canonical hostile is a directed triangle plus an isolated vertex: it has
the right arc count and degree caps but is not a path. The corrected near miss
is accepting a path-plus-cycles factor on endpoint degrees alone. With distinct
named endpoints the smallest hostile becomes a two-vertex path plus a triangle,
on five vertices. The least-used sidecar is the *jointly aligned target*: the
two separate tournament partitions can be coarse while their common refinement
is forced to be discrete.

Anchor: actual endpoint/shared-incidence refinement. Niche: the square-sum
path-cover boundary. Wildcard: the Pascal-parity evolution on a retained path.
Concept board: owned endpoints / agreement arcs / common equitable partition /
detached cycles / integer clock / phase transport. Every map below retains the
vertex ownership needed to compare its source and target.

## 1. An all-size common-refinement theorem

Let A be the zero-one adjacency matrix of a tournament on a fixed vertex set
V, |V|=n>=2. Let P=(v_0,...,v_(n-1)) be an actual directed Hamiltonian path.
Let B be the adjacency matrix of its complement partner K_P(A), which keeps
the path arcs and reverses every other arc. Let F be the directed path matrix.
Then, over the integers,

    C = A+B-(J-I) = F-F^T.                              (1)

Start with one vertex color. At each synchronous refinement round replace a
vertex's color by its old color together with its two vectors of outgoing
counts into the current cells, one vector from A and one from B. The counts
for A and B must use the SAME cells on the SAME vertices; independently sorted
histograms are not this operation.

**Theorem 1.** This joint refinement is discrete after at most floor(n/2)
rounds. In particular, the only partition equitable for both A and B is the
singleton partition.

**Proof.** The matrix J-I is equitable for every partition. Thus the two count
vectors determine all C-counts into each current cell. In round one the total
C-count is +1 at v_0, -1 at v_(n-1), and zero elsewhere, so the two actual
endpoints become distinct singleton cells. This also handles n=2 and n=3.

If v_i is already singleton, its C-column has exactly two possible nonzero
entries, at its adjacent path vertices, with opposite signs. Once the vertex
on the already exposed side is known, its other neighbor is uniquely picked
out from the remaining vertices by the C-count into {v_i}. Hence round r
exposes at least the first r and last r vertices of P individually. For odd n,
the one vertex left in the middle is automatically singleton. This proves the
bound. A common equitable partition cannot split under this refinement, so it
must already be discrete. QED.

The upper bound is attained in the exhaustive tests for each n=2,...,7. We
make no claim here that it is sharp for every n. The tested alternating-distance
tournaments, for example, stabilize in only two or three rounds for large n.

**Consequence for the seven/triple-cell frontier.** An equitable source
partition containing non-singleton cells cannot also be an equitable partition
of its *aligned actual partner*. Source-cell affine holonomy remains a valid
source invariant, but transporting the same partition unchanged to the target
would erase precisely the path incidences that force it to split. This is a
reason to query target incidence at the owned source cells.

**Limit.** Both the target and its alignment are inputs to this theorem. It
does not find either from the source alone. Joint discreteness is necessary,
not sufficient, for a candidate to be a complement partner. At n=4, for
example, the exact census has 3,840 nonnative aligned pairs with discrete joint
refinement, alongside 192 native pairs.

## 2. There is no residual common fractional symmetry

**Corollary 2.** If a nonnegative doubly stochastic matrix X commutes with
both A and B of an actual port, then X=I.

**Proof.** Such X commutes with J and I, hence with C. Since X1=1,

    X(e_start-e_end) = XC1 = CX1 = C1 = e_start-e_end.

A convex combination of the entries of e_start-e_end can equal its unique
maximum +1 only by putting all its mass at the start. Therefore the start row
of X is the start basis row, and the column-sum condition fixes its column as
well. The unique minimum fixes the end row and column in the same way. If the
first i vertices are fixed, commuting with C applied to e_(v_(i-1)) fixes the
next one, because C e_(v_(i-1)) is the difference of its two neighboring basis
vectors (with the appropriate one-sided endpoint case). Induction gives X=I.
QED.

Nonnegativity and both stochastic conditions are essential to this argument.
This is a rigidity statement for the already aligned two-matrix carrier, not
a claim that source-only linear programming solves tournament isomorphism.

## 3. Actual endpoint degrees still need one global guard

For any aligned tournaments A,B define the directed agreement matrix

    D_ij = A_ij B_ij.

Suppose distinct owned vertices s,t satisfy

    (in_D(s),out_D(s))=(0,1),
    (in_D(t),out_D(t))=(1,0),
    (in_D(v),out_D(v))=(1,1) for every other v.          (2)

These separate in/out conditions are stronger than the net divergence in
(1). They imply that D consists of one directed s-to-t path and zero or more
disjoint directed cycles. No cycle can have length one or two, since D is an
oriented subgraph of a tournament. The following are equivalent:

- the existing Session25 directed-agreement certificate accepts;
- the s-to-t path visits all n vertices;
- D is nilpotent;
- det(I-D)=1;
- (D^(n-1))_(s,t)=1.

Indeed, a path block is nilpotent with determinant one for I-D; a directed
cycle block has eigenvalue one and contributes determinant zero. With (2),
the matrix entry counts the unique s-to-t walk of length n-1, which exists
exactly when all vertices lie on that path. These equivalences require (2);
they are not arbitrary matrix substitutes for native port validation.

At n=5 retain the arcs 0->1, 2->3, 3->4, 4->2 of any tournament A that contains
them, and reverse every other arc to form B. The endpoint degrees (2) hold,
but det(I-D)=0 and only 0,1 are visited from the start. This is minimal with
distinct start/end: a detached oriented cycle needs at least three vertices
and their path needs at least two. The claim concerns this specified
alignment; it does not exclude another alignment between the isomorphism
classes from being native.

Converse sends (A,B,D,s,t) to (A^T,B^T,D^T,t,s). It preserves the acceptance
test and refinement partitions, and reverses the recovered path. No unsigned
identification is substituted for this operation.

## 4. Sierpinski's triangle appears on the actual retained path clock

Work over F2, and now suppose D=F is the authenticated path matrix. The
stay-or-advance operator U=I+F satisfies, for every integer time t>=0,

    U^t = sum_k binom(t,k) F^k
        = product_(j: bit_j(t)=1) (I+F^(2^j)).           (3)

This follows by the binomial theorem and the characteristic-two identity
(I+F)^(2^j)=I+F^(2^j). Acting on the start row, the coordinate at v_k equals
one exactly when every set bit of k is a set bit of t. Thus the space-time
diagram is precisely the Pascal-parity Sierpinski pattern, truncated at n
vertices. Here its positions are actual path vertices, and its clock is the
integer t; neither is an unowned motif count.

This is an exact directed-graph connection, not an additional proof of port
existence or a new compression algorithm: U at time one already contains F.
If one keeps only coarser dyadic powers, intermediate owner information can
be lost. Likewise a cycle component does not vanish under powers of F, so the
nilpotence guard in section 3 is the stopping condition needed before this
finite-path interpretation is justified.

There is a useful clock guard here. For n>=2, U has exact order
`2^ceil(log2(n))` over F2: its power at that exponent is I because F^n=0,
whereas the preceding power of two leaves a nonzero F-power. Thus an operator
snapshot aliases times separated by this period. The integer clock must be
retained separately; a periodic Pascal snapshot is not a cycle of the original
directed path dynamics.

### Incoming carry lift: exact periods modulo every power of two

Incoming commit `69648c025d`, [Session40's carry identity, section 3](../session40_marked_ports_square_sum_20261023/hopf-incidence/HOPF_SIERPINSKI_SHARED_INCIDENCE.md),
proves `v2(binomial(2^r,j))=r-v2(j)` for 0<j<2^r and its five-tap
modulo-eight consequence. Substituting the actual path matrix F into that
polynomial is valid; F^k vanishes for k>=n.

**Corollary 3 (clock lift).** For n>=2 and s>=1, the exact order of I+F over
Z/(2^s) is

    2^(s-1+ceil(log2(n))).                              (4)

**Proof.** Write m=ceil(log2(n)); then floor(log2(n-1))=m-1. For r>=m,
all nonconstant diagonals of (I+F)^(2^r) have coefficients with valuations
r-v2(j), 1<=j<n. Their smallest valuation is r-(m-1), achieved at
j=2^(m-1). Hence this power is I modulo 2^s exactly when r>=s+m-1.
For r<m, its F^(2^r) coefficient is one and that diagonal exists, so it
is not I. Since some power of two is the identity, the group order is a
power of two; the stated smallest such exponent gives (4). QED.

For N=2^r, r>=2, the inherited five-tap formula becomes

    (I+F)^N = I+4F^(N/4)+6F^(N/2)+4F^(3N/4)+F^N (mod 8).

This records the missing carries in the actual finite path carrier. It
does not identify this operator with the native metagraph Laplacian. For
example, n=31 gives snapshot periods 32,64,128 modulo 2,4,8 respectively;
these are clock-alias periods, not periods of a directed Hamiltonian path
and not Bott periodicity.

This Pascal recursion is distinct from the [three-piece affine boundary
gluing](../hopf_sierpinski_port_recovery_20261023/REPORT.md), whose state records
compatible owned boundary assignments. Both retain the gluing coordinates;
their carrier and composition law differ. Neither identity implies Bott
periodicity or the native electrical mod-8 obstruction.

## 5. Guard and verification ledger

The previous [marked mod-8 potential guard](../metagraph_five_port_phase_frontier_20261023.md)
is retained as an independent necessary test on candidate merged wiring.
Nothing here weakens its ownership or asserts that a passing electrical
signature is sufficient. Joint incidence and mod-8 data answer different
questions; they should be intersected, not identified.

Run `python -B 05-knowledge/results/aligned_port_refinement_20261023/check.py`
from the repository root. It exhausts all fixed-path tournaments for n=2..7
(33,867 total), with exact direct complement, signed cancellation, actual
agreement-path tracing, and converse controls. A separate enumeration checks
every aligned tournament pair for n=2..4 against a permutation-path recognizer.
The five-vertex detached-cycle hostile passes the endpoint degrees and fails
the determinant guard. Pascal parity is independently checked against both
chronological updates and binomial coefficients for n=2..64 and t=0..127.

The reproduction records finite evidence for the stated universe. The
all-size claims rest on the proofs, not extrapolation from those tests.

**Next obligation.** Starting with a real source partition and only partial
candidate-target records, choose the smallest owned incidence queries that
force the common partition to split, while retaining cycle rejection. That
is a meaningful compression problem beyond the complete-data theorem above.
