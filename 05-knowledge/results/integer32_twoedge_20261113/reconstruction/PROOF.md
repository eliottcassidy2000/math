# Exact quadratic-jet image and sharp repair from relative orientations

**Status: PROVED IN NOTES, R1–R3, after independent review.** Finite controls below are **FINITE-EXACT**. No theorem ID
is reserved. The inherited reconstruction and M−1 query results are not claimed
as new. The general binary cut-space mechanism is classical; the contribution
here is its exact path-observer realization guard and explicit optimal-radius
repair, with the distinction between one and two slots retained.

## Inheritance, live board, and scope

The closest proved mechanism is Session320's
[WEDGE_RECONSTRUCTION.md](../../session320_second_moments_and_native_contexts_20261113/WEDGE_RECONSTRUCTION.md):
positive one-slot wedges recover the positive directed support up to converse,
while two optional slots make each local mixed-orientation test observable.
Its fixed spanning-tree bank of `M−1` binary queries is already optimal for
labelled tournaments modulo converse. We do not rebrand that query count.

The canonical hostile is the strong order-24 module-converse pair from
[STRONG_SHADOW.md](../../session320_second_moments_and_native_contexts_20261113/STRONG_SHADOW.md),
which preserves all full-tournament per-slot first marginals and every induced
unmarked Q. The corrected near miss is promoting separate endpoint and edge
marginals to their joint observation. The least-used sidecar is the complete
set of consistency equations among relative edge orientations, including the
negative answers in the two-slot observer.

The board has five entries: module-converse fibres; actual two-contact tests;
local image constraints; erroneous relative signs; marked native incidence.
Anchor: exact realization of the reconstruction observer. Niche: robust repair
of inconsistent finite data. Wildcard: the induced action of native path moves,
owned separately by the root. The source is a labelled tournament; the target
is a Boolean field on labelled wedges. Global converse is lost, and labels and
the declared slot protocol must be retained. Ordinary-class native matrix C,
vertex weighted shadow W, and the Boolean wedge field are different objects.

## Definitions and orientation gauge

Fix `n>=3`, vertices `0,...,n−1`, and `E=binom(V,2)` of size `M`.
A wedge `(v;u,w)` consists of edges `{v,u}`, `{v,w}`, with `u<w` and all
three vertices distinct. There are

    L = n binom(n−1,2) = 3 binom(n,3)

wedges. Put `beta(v;u,w)=1` exactly when one incident arc points into v
and the other points out. This is the support of their joint occurrence
in a two-slot spanning path cover, allowing empty slots: in the positive
case put the directed triple in one slot and a Rédei path of its complement
in the other. Conversely a cover using both edges must use them consecutively
in one slot. This is inherited, with the protocol made explicit.

Equivalently, in the labelled unordered-cover polynomial

    J_T(t;y)=sum_C t^(n−|C|) product_(e used by C) y_e,

these are the coefficients of `t² y_e y_f` for incident edges. The constant
term is 1, every linear edge coefficient is 1, disjoint-pair coefficients are
1, and repeated-edge coefficients vanish. These fixed coefficients plus beta
are the entire quadratic truncation. Here t is a formal edge-count variable,
not the number of optional slots or a literal field of size one.

Reference orientation is `i->j` for `i<j`. Write `x_e=1` for an arc reversed
from that reference. Its wedge field beta0 is 1 iff `u<v<w`. Then

    beta(v;u,w) + beta0(v;u,w) = x_{uv}+x_{vw} mod 2.       (1)

Thus relative to beta0 the observer is the binary edge coboundary of x on
the line graph `L(K_n)`, whose vertices are the tournament edges.

## R1. The exact image has local parity checks

**Candidate theorem.** A proposed Boolean field beta is realized by a
labelled tournament iff both families below hold:

1. At each v, for three distinct other vertices a,b,c,

       beta(v;a,b)+beta(v;a,c)+beta(v;b,c)=0 mod 2.         (2)

2. For each three distinct vertices a,b,c,

       beta(a;b,c)+beta(b;a,c)+beta(c;a,b)=1 mod 2.        (3)

When realized, exactly two labelled tournaments realize it, and they are
mutual converses. The image has `2^(M−1)` elements and is an affine binary
code of length L and dimension `M−1`. Its codimension is `L−M+1`.
Together with the fixed coefficients above, (2)–(3) are an exact realization
criterion for the quadratic path-cover truncation. They do not certify a
proposed arbitrary higher-degree polynomial or a prescribed native matrix C.

*Proof.* Let `a_vu=1` mean `v->u`. Then
`a_vu+a_uv=1` and `beta(v;u,w)=a_vu+a_vw`. Cancellation proves (2), and
summing around a vertex triangle proves (3).

Conversely, (2) says that at each v there exist incidence bits `a_vu`
realizing all beta values on that star. Choose one neighbour as a local
reference and set the others using its pairwise beta values; the triangle
checks ensure all pairs agree. Each such local choice is defined up to one
constant bit. Put `c_uv=1+a_uv+a_vu`. Equation (3) implies that c sums to
zero around every triangle of the complete graph. Consequently there are
vertex bits g_v with `c_uv=g_u+g_v`: fix g_0=0 and set g_v=c_0v. Adding g_v
to every incidence at v now gives `a_uv+a_vu=1` on each edge and leaves
beta unchanged. These are the arcs of a tournament.

Equation (1) also shows uniqueness modulo global converse. If two reversal
masks give the same beta, their difference is constant across every adjacent
pair of edges. The line graph is connected, so the difference is constant
everywhere: either no arcs or all arcs change. The linear kernel has dimension
one, proving the image size and dimension. The local equations are therefore
a possibly redundant complete set of parity checks of the stated codimension.
QED.

Both check families matter. At n=3 all-zero beta passes every star condition
(vacuously) but violates (3). At n=4, start from beta0 and toggle the two
wedges at vertices 0 and 1 of triangle `{0,1,2}`. Every equation (3) still
holds, but the corresponding star conditions fail. These are exact
nonrealizable formal observer fields, not hypothetical tournaments.

The cases n=1,2 have one observer value up to converse and no wedge queries;
the distance and repair statements below are restricted to n>=3.

## R2. Exact distance: one arc changes the fewest wedge bits

If two tournaments differ on a subset F of unordered pairs, their wedge
fields differ exactly where one incident edge is in F and the other is not.
Writing d_F(v) for its degree in the graph F, their Hamming distance is

    sum_v d_F(v)(n−1−d_F(v)) = |cut_(L(K_n))(F)|.          (4)

**Candidate theorem.** The minimum nonzero distance in the image is
`d_min=2n−4`. This is a metric on fixed labelled wedge coordinates, not a
metric lower bound after arbitrary relabelling.

*Proof.* Flipping exactly one tournament arc changes its comparisons with
all `2(n−2)` incident other edges, giving the upper bound. For the lower
bound, any nontrivial line-graph cut has adjacent vertices
`e={v,u}`, `f={v,w}` on opposite sides. The following are `2n−4`
edge-disjoint line-graph paths from e to f:

    e--f;
    e--{u,w}--f;
    e--{v,z}--f                       for z outside {v,u,w};
    e--{u,z}--{w,z}--f                for z outside {v,u,w}.

Their distinct edge pairs show disjointness directly. Each crosses the cut
at least once, so the cut contains at least `2n−4` edges. A change of all
arcs is the excluded zero-distance global converse. QED.

This proves exact distance without assuming an unproved expansion property
of tournament metagraphs. The auxiliary line graph is explicitly `L(K_n)`.

## R3. A constructive sharp error radius

**Candidate theorem.** Any full wedge field within `n−3` bit errors of a
realizable beta can be uniquely repaired in `O(n^4)` elementary bit operations.
No uniform unique repair is possible at radius `n−2`.

*Decoder and proof.* Subtract beta0 from the observed bits. For each adjacent
pair e,f in the line graph, compute the XOR along each of R2's `2n−4`
paths. On an error-free field every such path has value `x_e+x_f` by
telescoping. Since the paths are edge-disjoint, each erroneous observed wedge
can affect at most one vote. At most `n−3` errors leave a strict majority
of correct votes. Repeat for all L pairs to recover the entire true beta,
then integrate (1) along a spanning tree with `x_first_edge=0`.
There are `O(n^3)` wedges, `O(n)` constant-length votes for each, and linear
propagation afterward. The same procedure rejects tied votes or an output
that fails R1; finally verify that its distance to the original input is at
most `n−3`. This last check makes success a certificate of the stated radius,
including for inputs outside the guarantee.

For sharpness, take two image words differing by one tournament arc. Their
distance is `2n−4`. A field agreeing with one on half of those differing
coordinates and with the other on the rest has distance `n−2` from both.
No unique decoder can resolve that midpoint without further information.
QED.

The guarantee concerns errors in the **full two-slot local beta field**, or
in the corresponding labelled quadratic jet. It does not apply to zeros in
the one-slot HP co-usage field: absence there can mean failure of a spanning
extension rather than a local orientation clash. For example, the transitive
order-four tournament has just two positive HP wedges, whereas its two-slot
field has four. Treating the HP field as beta violates R1.

The decoder is equivariant under addition of a valid linear cut word after
subtracting beta0: path votes shift by the same endpoint difference. This
justifies testing all bounded error patterns at one codeword plus separate
nonzero-word controls; it is not a claim of full arbitrary-field enumeration
at larger n.

## Exact controls and the retained limits

`check.py` exposes `edges(n)`, `wedges(n)`, `encode(n,x)`, `decode(n,beta)`,
and `repair(n,noisy)`. Wedges are ordered by v, then the sorted other pair.
Exact decoding returns a reversal mask with bit zero equal to zero, or None.
Bounded repair returns the same type and certifies the stated radius.

The checker independently enumerates all image codewords through n=6 and
all possible formal wedge fields through n=4. It computes the rank of the
local check matrix, the full distance distribution, and the explicit short
path packing. Every error pattern of size at most n−3 at the reference
codeword is checked through n=6, together with translation controls and the
ambiguous midpoint. It separately reconstructs one-slot HP supports and
literal two-slot covers for every inherited ordinary representative through
six, using immutable metadata only. No old producer implementation is imported.
A labelled six-vertex cyclic-module control checks all induced H values and
full per-slot first marginals at q=1,2,3 by an independent assignment-of-vertices
algorithm; it is not a second census of the 24-vertex witness.

The only nontrivial class-bank input is pinned in `source-pins.json`. All
new receipts are small; no new large tournament bank is written. These results
provide an exact image and repair for a declared observer, not a scalar
A000568 recurrence, a proof that arbitrary prescribed incidence is realizable,
or a closure of the unrestricted Hamiltonian spectrum.

Reproduce from the worktree root:

```sh
python3 05-knowledge/results/integer32_twoedge_20261113/reconstruction/check.py > 05-knowledge/results/integer32_twoedge_20261113/reconstruction/normal.json
python3 -O 05-knowledge/results/integer32_twoedge_20261113/reconstruction/check.py > 05-knowledge/results/integer32_twoedge_20261113/reconstruction/optimized.json
```

Both receipts pass 57,186 explicit gates and agree byte for byte. The image
sizes at n=3,4,5,6 are 4,32,512,16384; minimum distances are 2,4,6,8.
The exhaustive reference-error pattern counts are 1,13,466,36051. Actual
one-slot and two-slot controls cover 2,4,12,56 ordinary representatives,
respectively, without extending the inherited class bank.

| Artifact | SHA256 |
|---|---|
| `check.py` | `2f6991f470496b2b802fe8a20760fd459562148678b9015be2fab1272cc8ab05` |
| `normal.json` | `bdd9fbc57328a637c1153cde5bea2c7c6f29d650d5744b6eb1d758055dbfcab1` |
| `optimized.json` | `bdd9fbc57328a637c1153cde5bea2c7c6f29d650d5744b6eb1d758055dbfcab1` |
