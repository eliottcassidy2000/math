# Session 192: weighted-cut geometry of tournament birth

**Status: PROVED elementary reductions; exact certificates supplied. Independently audited.**
The global Hamiltonian-path spectrum conjecture remains OPEN. This packet
connects the inherited ear response to graphical zonotopes, integral network
flows, braid-fan coarea and toric normalization. It supplies exact operations
and counterexamples, not geometric analogies offered as proofs.

Write H(U) for the number of directed Hamiltonian vertex words, H(empty)=1.
For a new vertex x let S={i:x->i}, and put Delta_U(S)=H(U+x)-H(U).
Let Q_ij count old-vertex words having i immediately followed by j, all other
adjacencies valid, while the exposed i,j step is unrestricted. All weights
below are recomputed from this same parent U; they are not fixed under births.

## 1. A positive capacity on each old arc

Set c_ij=Q_ij+Q_ji. This is symmetric in i,j. Then

    Delta_U(S) = sum_{i in S, j outside S, i->j in U} c_ij.       (1)

Thus a root triangle x->i->j->x corresponds exactly to an old outgoing cut
arc, and its intrinsic weight is c_ij. Every cut response is a nonnegative
weighted triangle sum, including the contribution of all longer odd cycles.

### Proof of the cancellation

Inherited THM-4097, `01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`,
and THM-4104, `01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md`, give

    Delta_U(S)=sum_{i in S}(start_i-end_i)
                 +sum_{a outside S,b in S}Q_ab.              (2)

For an actual Hamiltonian path, outgoing adjacent steps minus incoming ones
are +1 at its start, -1 at its end, and zero elsewhere. Summing paths gives

    start_i-end_i = sum_{i->j}Q_ij - sum_{j->i}Q_ji.           (3)

Summing (3) on S cancels all internal edges. At a crossing edge directed
outside-to-inside, its negative valid Q is canceled by (2). At an edge
i->j directed inside-to-outside, (2) supplies Q_ji in addition to Q_ij.
This is (1). No approximation, selected orientation or strongness is needed.

The inherited parity/second-difference argument gives c_ij even. There is
also a useful strictly positive bound:

    c_ij >= 2 H(U-{i,j}) >= 2.                              (4)

To prove it, take a Hamiltonian word z_1,...,z_k of U-{i,j}. Insert the
adjacent block (i,j) or (j,i) into every valid gap; its internal step is
exposed. For k>=1 write a_l=1_(z_l->i), b_l=1_(z_l->j). The number of valid
insertions of the two block orientations is exactly

    2 + sum_{l=1}^{k-1}[(1-b_l)a_(l+1)+(1-a_l)b_(l+1)] >= 2.

For k=0 both words are valid. Deleting the block recovers the original
word, so all counted outputs are distinct. Q may contain further words
whose deletion creates an invalid adjacency; we claim a lower bound, not
equality. This proves (4).

The ordinary number of root triangles is the same directed cut with every
capacity replaced by 1. That replacement need not preserve even the minimum
cut; Section 5 gives an explicit strong-parent counterexample.

## 2. Exact spectrum reformulation and integral certificates

For strong U on n>=3 vertices, every nonempty proper S produces a strong
child. Conversely Moon's strong-deletion theorem gives a strong-deletion
parent for every strong tournament of order >=4, as inherited in THM-4097.
Consequently, with S_n^str denoting the full strong H-spectrum,

    S_(n+1)^str = union over strong U of order n of
      {H(U)+capacity_U(delta^+(S)): empty!=S!=V(U)}.           (5)

This equality retains the **discrete weighted directed-cut support** of
actual parents. Neither an interval between extrema nor arbitrary weights
can replace that support.

For arbitrary U of order at least two, let A and B be its source and sink SCCs. The child U+x_S
is strong exactly when S meets A and its complement meets B. Necessity
follows because x must reach A and be reached from B. Sufficiency follows
by walking through the ordered old SCCs. Hence the minimum strong-child
increment is

    min_{a in A, b in B, a!=b} minimum a-to-b cut capacity.   (6)

For each pair use the directed network U with arc capacity c_ij. Force
x_a=1 and x_b=0 by the usual source/sink convention; no auxiliary negative
capacities or offset is required. Integral max-flow/min-cut provides an
integer flow and cut having equal value. The verifier checks every capacity,
every conservation equation, and equality of flow value and cut capacity.
There are |A||B|-|A∩B|<=n(n-1) eligible terminal pairs. Computing Q still requires
its exact path-boundary information (the supplied subset DP costs
O(n^2 2^n)); (6) removes exhaustive 2^n cut scanning after that information
has been computed. It does not make Hamilton-path counting polynomial-time.

For fixed a,b, minimizing cuts are closed under union and intersection:
submodularity gives f(X)+f(Y)>=f(X union Y)+f(X intersection Y), and both
new sets retain a and omit b. Thus the minimizers form a distributive
sublattice. This pinning is indispensable. For U=C3, all six nontrivial
cuts minimize Delta=2, but {0} intersection {1}=empty is inadmissible.

A single optimal flow certifies only its pinned pair. Certifying the
unrestricted strong-child minimum requires (6), or an independently proved
coverage reduction. We retain all pair certificates, not just the best cut.

## 3. SCCs are precisely the face strata of a weighted permutohedron

Set w_ij=c_ij/2 and, in R^n modulo the constant direction, form the
centrally symmetric graphical zonotope

    Z_U = sum_{i<j}[-w_ij,w_ij](e_i-e_j),
    h_U = sum_{i->j} w_ij(e_i-e_j).                         (7)

All weights are positive by (4), so Z_U has dimension n-1 and its normal
fan is the braid fan, exactly as for a permutohedron. This uses standard
zonotope geometry, not a new general theorem about permutohedra.

The tournament point h_U lies in Z_U because it chooses one endpoint of
each generating segment. Importantly, such an endpoint choice need not be
a vertex: cyclic orientations put it on a larger face or in the interior.
For every nonempty proper S, the slack of the inequality with normal -1_S is

    support_Z(-1_S)-h_U dot(-1_S)
      =sum_{i in S,j outside S} w_ij + sum_{i in S}h_i
      =Delta_U(S).                                         (8)

Hence the minimal ear increment of a strong parent is its smallest integer
facet slack. This is an arithmetic slack, not a Euclidean distance: the
facet normals have different Euclidean lengths.

Let U have k SCCs, ordered C1->...->Ck. For any y,

    support_Z(y)-h_U dot y
      =2 sum_{i->j} w_ij (y_j-y_i)_+.                      (9)

Equality holds exactly when y_i>=y_j on every old arc. This forces y to be
constant on each SCC and its k component values to decrease weakly in the
SCC order. Modulo constants this normal cone has dimension k-1. Therefore

    dimension of the smallest face of Z_U containing h_U = n-k.  (10)

In particular U is strong iff h_U is interior. If U is transitive, h_U is
a vertex. More generally the active cut-facet normals correspond exactly
to the k-1 nonempty proper SCC suffixes; including the trivial cuts gives
exactly k+1 zero-increment attachments. This recovers and geometrizes S191's
neutral-padding law without claiming that the scalar H itself determines
any of this face data.

**Map and loss ledger.** U maps to the labelled pair (Z_U,h_U). It preserves
every one-birth response and its SCC face partition. The unlabelled shape
alone loses the orientation-score point; h alone loses the facet weights.
After a birth one must recompute the next weighted zonotope. A fixed
polytope is not an iterated growth carrier.

## 4. Coarea and the convex/toric normalization firewall

The convex piecewise-linear extension of (1) is

    E_U(z)=sum_{i->j}c_ij max(z_i-z_j,0),  z in [0,1]^n.

It is linear on every braid simplex. Writing S_t={i:z_i>t} gives the exact
coarea identity

    E_U(z)=integral_0^1 Delta_U(S_t) dt.                    (11)

On a pinned face z_a=1,z_b=0, all thresholds contain a and omit b.
Therefore some physical threshold cut has value at most E_U(z), so the
convex minimum equals the integral minimum. This is the geometric reason
min-cut relaxation is exact. Equality E_U(z)=m, however, only makes m an
average of physical values; it gives no cut with value exactly m.

The S191 unique H=13 core U=(6,8,26,17,3) makes this failure concrete.
Its complete strong one-birth H support is

    {23,25,29,31,33,37,43,45}.

Cut 4 has H=37; cut 5 has H=43. These cuts are nested and share the pin
2 in S, 1 outside S. Give coordinate 2 value 1, coordinate 0 value 1/3,
and all others value 0. Its thresholds are cut 5 for 0<t<1/3 and cut 4
for 1/3<t<1. Thus H(U)+E_U(z)=39 although no physical one-birth cut has H39.
S191 proves the stronger all-order exclusion of strong H39 descendants
containing this core; the present convex hostile does not replace that proof.

There is an exact toric normalization version. Put

    A={Delta_U(S)/2:empty!=S!=V}={5,6,8,9,10,12,15,16},
    M=< (1,a): a in A > subset Z^2.

Its group is Z^2 (consecutive values 5,6 suffice), its cone is
5d<=k<=16d, and its saturation fills every lattice point in this cone.
Indeed 11(d,k)=(16d-k)(1,5)+(k-5d)(1,16), with nonnegative
integer coefficients for every lattice point in that cone.
The degree-one holes are k=7,11,13,14; the only degree-two hole is k=29;
there are no holes at degree>=3. For the last assertion, direct addition
gives 3A=[15,48], and adding A to an interval of that length creates no gap,
so induction gives dA=[5d,16d] for d>=3. The full semigroup therefore has
exactly five holes. The H39 attempt is the degree-one hole (1,13).

This normalization is an algebraic relaxation, not a tournament operation.
Adding two semigroup columns sums two independent ear energies; an ordered
join of the actual child tournaments multiplies their H values. Dropping
that operation label would turn a true toric saturation into a false
spectrum construction. It is a precise warning against using normality,
convexity or lattice gcd alone to certify physical landing.

## 5. A strong parent whose minimum ear uses MORE triangles

Take

    U=(24,121,235,48,4,145,185,27), H(U)=241.                 (12)

U is strong. Among all nontrivial attachment cuts, the minimum root-triangle
count is 1. The best such cut has Delta=224 (for example S=16), hence H465.
But the minimum increment over ALL cuts is 220, attained precisely by
S=8 and S=253, each having two root triangles; those children have H461.
Thus an exact weighted-cut minimizer need not be an unweighted-cut minimizer.

The discovery scan exhausts all strong labelled order-six parents and then
samples orders seven and eight deterministically. It found (12) at order
eight. It does NOT prove order eight is the smallest possible order: the
order-seven phase was sampled. The separate verifier checks every one of
this parent's 256 children by direct subset DP and both winning children
plus the one-triangle hostile by literal permutation counting.

This is not a counterexample to S191's global chain-minimum conjecture.
At child order nine the global minimum increment is only 26, much smaller
than 220. What it refutes is a potentially tempting fixed-parent reduction:
'reorient the root optimally, and some optimal root lies in one triangle.'
S191's separate greedy arc-smoothing hostile changes the parent; (12)
instead keeps EVERY old arc fixed and globally optimizes all root arcs.

## 6. Sources of inspiration and surviving approaches

1. **Permutohedral order geometry.** The repository reflection
`07-reflections/lrc-permutohedron-geometry-and-n14-attempt-s521o.md` separates
order from metric walls. Here the faithful weights are c_ij and the target
point is h_U; actual facet slack, unlike a bare order skeleton, is exactly
Hamilton growth. The reflection is inspiration, not a theorem dependency.
2. **Owner-transition min-cuts.** THM-2267,
`01-canon/theorems/THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel.md`,
uses cut certificates to retain gluing data. Formula (1) transports that
method to exact ear growth with no multi-owner relaxation gap.
3. **Polytope selector warnings.** THM-3137,
`01-canon/theorems/THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall.md`,
separates averaged selectors from physical objects. Section 4 supplies the
matching exact-growth hostile, even on one pinned braid simplex.
4. **Toric quotient/normalization.** The algebraic-geometry work including
THM-2615, `01-canon/theorems/THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice.md`,
requires retaining the physical map rather than quotient support. The
five-hole ear semigroup gives a small exact model of that same loss.

Promising next precise question: which additional constraints on the
intrinsic capacity matrix C_U force its minimum outgoing cut to match a
low-triangle cut? Generic nonnegative capacities cannot answer this, and
(12) rules out the unrestricted tournament version. For spectrum coverage,
study exact integer cut supports across parent banks, not convexified
supports. For S191 growth, use flow certificates to replace repeated
exhaustive root-cut scans, then seek inequalities on Q that survive changing
the parent. None of these statements settles the global {7,21} conjecture.
