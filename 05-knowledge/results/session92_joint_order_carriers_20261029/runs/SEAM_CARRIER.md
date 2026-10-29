# Exact native-complement seams and a strong-source infinite hostile family

**PROVED (elementary, scoped); FINITE-EXACT controls.** This is an independent
mechanism, not a recovered replay of the user's four quoted Hamilton counts.
The exact ordered run lists and definition of each individual cursor response
have not yet been located. In particular, no equality of cursor responses is
claimed for the independent example below.

## Inheritance and the crucial distinction

Source pin: `b07ec51b204e5c0e09e52eed206cb2bd3cc29ebc`.
The inherited `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`
proves that ordinary substitution `Q[S_i]` has Hamilton count determined by
all child unordered path-cover counts. Reordering path runs while leaving the
children and quotient unchanged does **not** change that tournament or its H.
An obstruction using different assembled partner counts must retain the
operation that produces the partner. Native complement is path-dependent;
its definition appears in `strong_attachment_tail_compiler_20261028/native.py`
and the presentation conventions routed there. The theorem below treats that
operation explicitly, without presuming it is the user's missing construction.

Concept board: ordinary substitution / path-dependent complement / seam
relation / shared endpoint / exact reconstruction. The corrected near miss
is confusing a failure of a transformed observable with a failure of the
ordinary path-cover substitution law. The least-used sidecar is the exact
set of restored Hamilton path edges, rather than a per-run statistic.

## 1. Seam identity

Let T be a named tournament, A its skew adjacency matrix (+1 on arcs), and P
a directed Hamilton path. Define its native partner C_P(T) by reversing every
arc except the consecutive arcs of P. Let D_P be the skew adjacency matrix
of the sparse directed path P. Then

    B_P = -A + 2 D_P,
    B_P - B_Q = 2(D_P-D_Q).

This follows entry by entry. On a path edge D_P=A; elsewhere D_P=0.
Consequently the exact unordered edge set on which the partners differ is
`E(P) symmetric_difference E(Q)`. No numerical approximation or genericity
assumption enters.

If P and Q are assembled from the same directed run components, each used
once and internally unchanged, their internal path edges cancel. The entire
partner discrepancy lies at the seams connecting the runs. Thus an exact
sufficient carrier for this fixed underlying T and fixed run family is the
**joint directed seam relation**. Its isolated marginals need not specify it.
This is a carrier of the labelled partner, not a proved minimum for H alone.

For fixed T, the labelled partner actually reconstructs P: its agreement
edges with T are precisely the n-1 path edges, whose unique source and
successor relation recover the ordered Hamilton path. Hence a summary closed
under exact labelled-partner output cannot compress distinct Hamilton paths.
This does not prevent small certificates for a particular H comparison, nor
compression for unlabelled outputs or weaker operation languages.

## 2. Strong-source padding theorem

Suppose P and Q start at the same vertex r of a strong tournament T of order
at least three. Form T+ by adding x with the sole outgoing old arc x->r;
every other old vertex beats x. Prefix x to both paths.

* Both extended paths are Hamilton paths of T+.
* T+ is strong: x reaches all old vertices via r, and every old vertex reaches
  an old vertex different from r (trivial unless starting at r), then x.
* In either native partner, x beats **all** old vertices: x->r is preserved,
  and every other incident arc reverses. Therefore

      C_(xP)(T+) = {x} -> C_P(T),
      H(C_(xP)(T+)) = H(C_P(T)).

Any Hamilton path of a source join must start at its source, giving the last
identity. Iterate, using the new shared first vertex at each step. Thus every
shared-start counterexample gives an infinite family at all larger orders,
with strong original tournaments and exactly the same two partner H counts.
The partners after padding are not strong; this is explicitly not an
all-children-strong or all-partners-strong extension theorem. It preserves
any property determined solely by the unchanged induced old pair, but does
not automatically preserve a cursor signature whose language changes with n.

## 3. Independent exact seed

Rows encode outgoing neighbours by binary bitmasks. Let the regular T5 have
rows `(6,12,24,17,3)`, equivalently i->i+1,i+2 modulo5.
Take

    P=(0,1,2,3,4), Q=(0,2,3,4,1).

Their partners have rows

    (26,20,9,18,4), (28,1,10,18,6)

and H counts **11 and 15**. The changed unordered edges are
`{01,02,12,14}`. The included checker verifies strong sources and preserved
counts at orders5 through9, along with the exact seam symmetric difference.
The all-order padding claim is justified by the proof, not the finite tests.

Reproduce (standard-library Python):

    python3 check_seam_carrier.py
    python3 -O check_seam_carrier.py

The subset dynamic program counts all directed Hamilton paths, with arbitrary
initial vertex, and uses explicit failures rather than optimized-away asserts.
No comparison here divides H by an automorphism order.

## 4. What is still needed for the quoted obstruction

For each quoted pair, recover: all child adjacency rows; the quotient and
its precise walk; both ordered run lists; the definition and range of each
cursor query; the assembly and partner operation; and whether H counts paths
or automorphism-orbits. Only then can 2483993/2491877 and2899675/2906091 be
independently replayed or the padding theorem specialized to them. A numerical
match alone would not certify the claimed equality of all cursor responses.

## 5. Arbitrarily high correlation order under independent restorations

Here the operation language is deliberately enlarged: independently restore
any subset of the edges of a fixed spanning undirected path F in the reverse
of T. Such a subset is not required to be the edge set of a Hamilton path,
so this theorem is not a legal-seam or cursor theorem.

Let F have consecutive vertices0,...,n-1; let x_i be its restoration bits,
and let s_i=+1 if T has i->i+1, and -1 otherwise. The Hamilton count of the
restored tournament is a multilinear polynomial f(x) of degree at most n-1.
Indeed each candidate vertex permutation contributes the product of its
n-1 adjacency indicators, each constant or affine in one distinct x_i.
The coefficient of all n-1 bits is exactly

    [x_0 ... x_(n-2)] f = (1+(-1)^(n-1)) product_i s_i.

Only two vertex permutations use every edge of F as a consecutive pair:
the forward and backward traversals. Their leading coefficients are
product s_i and (-1)^(n-1) product s_i. Every other permutation omits some
variable. This proves the formula for **every** T, not merely a random core.

For odd n, the coefficient is ±2. Take the uniform distributions over bit
strings of even and odd parity. Their projections on every proper subset
of coordinates are identical: the omitted bit uniquely completes either
parity, and additional omitted bits give equal completion multiplicities.
Nevertheless their expected Hamilton counts differ by

    E_even f - E_odd f = 2^(3-n) product_i s_i.

Thus arbitrarily high-order dependence is visible to this specific tournament
observable even when every proper marginal is supplied. The expected gap is
rational because these are ensembles; individual Hamilton counts remain
integers. For even n the top coefficient cancels, so this particular parity
construction has zero gap. This is an exact oddness mechanism, not a claim
that all even-order joint dependence disappears.

The companion `check_high_order.py` verifies two explicit orientations at each
n=3,...,8, every restoration mask, every maximal proper marginal, and the
coefficient and expectation identities, under ordinary and optimized Python.
No claim is made that these restoration ensembles obey native Hamilton-path
guards. That stricter realizability problem remains separate.

### 5a. Even-order repair and the exact global polynomial degree

Let all `binomial(n,2)` reversal/restoration bits vary independently. For odd
n>=3 the total Hamilton-count polynomial has degree exactly n-1. For even
n>=2 its degree is exactly n-2.

The odd claim follows above. For the even upper bound, any degree n-1 monomial
arising from a candidate Hamilton path is supported on its spanning path.
The two orientations cancel as above; no other edge support can contribute
at that degree. For the lower bound at even n>=4, use a path on n-1 vertices
and one isolated vertex, and examine the coefficient of its n-2 edges. Each
contributing permutation treats the path as a contiguous block, in one of
its two orientations, with the isolated vertex either before or after. Across the two orientations and the two placements, the four bridge
indicators are z->a, b->z, z->b, a->z, where a,b are the path endpoints.
They sum to2 by pairing opposite arcs. The internal coefficient is the same
for both orientations because n-2 is even. The coefficient is therefore
`2 product s_i`, independent of all the other arcs. At n=2 H=1, degree0.

Consequently even n>=4 also supports an exact parity-ensemble obstruction,
using n-2 bits rather than n-1: all proper marginals agree while the expected
H gap is `2^(4-n) product s_i`. Oddness changes the highest detectable order
by one; it does not remove the joint-information obstruction.

## Independent review

The Session92 audit lane independently checked the H11/15 seed by all120
vertex permutations per partner and source-join padding through order13.
Its independent parity audit covered all1098 labelled tournaments at
orders2 through5 using permutation constraints rather than this packet’s
subset DP. The auditor also reviewed the all-order proofs and the separate
n=2 degree-zero boundary. Audit artifacts live in the parent `audit/` lane;
this statement reports that review rather than copying its implementation.
