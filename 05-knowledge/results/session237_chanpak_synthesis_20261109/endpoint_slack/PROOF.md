# Endpoint extremality and the tail pair-energy bridge

**Status: PROVED; independent audit PASS.** Session237, 2026-11-09.
No full Hamilton-path spectrum theorem is asserted.

## Inheritance and exact advance

The nonsink endpoint bound `2 e_a <= H-1` and its insertion injection are
[THM-4876, statement7](../../../../01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md).
Its [audit D](../../chanpak_endcount_audit_D_20261104.md) observed the equality
classification only through order9. Session231 [descent ROUND2](../../session231_overnight_completion_20261109/descent/ROUND2.md)
proves related zero-defect and minimum-join faces, including the endpoint
extremal slice with extra cover defect1; it does not classify endpoint equality
without that extra hypothesis. THM-5001/session233 endpoint universality realizes
every positive endpoint count, but does not minimize the total H.

The pair-distance energy Q and its digit recurrences are inherited from
Sessions197/198/205/206 and THM-4888/5085. The new synthesis identifies that
same energy with exact endpoint slack on a physical apex-chain family and
classifies equality for arbitrary tournaments. No claim of external priority
is made.

## 1. Complete equality classification

Let T be a tournament, a a nonsink vertex with e_a(T)>0, and h=H(T).
Then `h=2e_a+1` iff T is a transitive prefix followed by a tournament R_n
(n>=3), rooted at its penultimate vertex. Here R_n is the transitive tournament
`v1 -> ... -> vn` with only `vn -> v1` reversed. At n=3 any vertex is a
valid root by cyclic symmetry. Consequently

    (h,e_a)=(2^(n-2)+1, 2^(n-3)).

The positive endpoint hypothesis is essential: an arbitrary nonsink in a
transitive tournament has (h,e)=(1,0) and also satisfies the displayed scalar
equality. Such roots are excluded here.

### Proof

Hamilton paths respect the ordered strong components. The root lies in the
last component C. Write M for the product of the H-counts of all earlier
components. Then `h=M H(C)` and `e_a=M e_a(C)`. Equality gives
`M(H(C)-2e_a(C))=1`, so M=1 and the earlier components form a transitive
prefix. The last component has at least3 vertices because a is not a sink.

We may therefore suppose T strongly connected. Put O=N+(a) and let e_O be
the number of paths ending in O. The inherited injection is explicit: in a
path ending a, move a immediately before its first outneighbor. Its image
ends in N-(a), not a or O. Hence `2e_a+e_O<=h`. Every vertex of a strong
tournament ends at least one Hamilton path (cut a Hamilton cycle), so equality
`h=2e_a+1` implies O={b} and e_b=1.

Removing b leaves a as a sink. Every Hamilton path of T-b thus ends a and
can be followed by b; conversely every path ending b arises in this way.
Thus `H(T-b)=e_b=1`. A tournament has one Hamilton path iff it is transitive:
a cyclic triangle already has3, and vertex insertion injects paths into any
larger tournament. Write T-b as `x1 -> ... -> x_l -> a`. Strongness forces
b->x1. Set S={i:b->x_i}. The exact calculation in Section2 gives

    h-2e_a=1+Q(S),

where Q(S)>=0, and Q(S)=0 iff |S|<=1. Since 1 belongs to S, equality forces
S={1}. This is exactly R_(l+2), rooted at a. Conversely Section2 or the
known R_n formula gives equality for that family. This proves both directions.

## 2. The physical pair-energy identity

Start with `x1 -> ... -> x_l -> a`, append b with a->b, and choose
b->x_i exactly at positions S. All other x_i beat b. Empty S and l=0 are
allowed. Define

    m(S)=sum_(i in S) 2^(l-i),
    Q(S)=sum_(i<j in S) c(j-i),
    c(1)=2,  c(d)=2^(d-1) for d>=2.

Then exactly

    e_a=m-Q,
    H=1+2m-Q,
    H-2e_a=1+Q.                                      (1)

Distances are unchanged by reversing the bit order, so Q(S) is exactly the
inherited tail energy Q(m) of the binary integer m, not merely an analogy.

### Direct counting proof

Every path not ending a must contain a,b consecutively because b is a's
unique outneighbor. There is exactly one path ending b, namely the full
transitive chain followed by b. All other paths not ending a have shape
`P,a,b,Q`, where P,Q partition x1,...,x_l, each in increasing order,
Q is nonempty, and min(Q) belongs to S. The number of these partitions is
m: fixing min(Q)=x_i leaves the l-i later positions free, and earlier positions
must lie in P.

Paths ending a instead have shape `P,b,Q,a`, under the same conditions,
with the additional condition that P is empty or max(P) does not belong to S.
Thus e_a=m-B, where B counts partitions with both max(P) and min(Q) in S.
For a pair i<j in S, the case max(P)=x_j, min(Q)=x_i gives
`2^(j-i-1)` partitions. The reverse case max(P)=x_i,min(Q)=x_j is possible
only if j=i+1, and then gives one extra partition. Thus an adjacent pair
contributes2 and every other pair contributes2^(j-i-1). Hence B=Q,
proving all of (1). This also proves positivity and the Q=0 boundary without
using a numerical recurrence or a full-count conjecture.

In particular the arithmetic root-tail sequence `1+2m-Q(m)` of the directed
edge seed is the full path count of this explicitly labelled apex family.
Its second coordinate `m-Q(m)` is a genuine endpoint count, while Q(m) is
the excess over the optimal endpoint half bound. This is an exact map,
retaining the bit order and chosen endpoint. It destroys the general rooted
cover response and supplies no arbitrary-target selector.

## 3. Sharp consequence for endpoint universality

For N>=1 let f(N) be the least H of a strongly connected tournament having
an endpoint with exactly N spanning paths. THM-5001 makes f(N) finite.
The preceding theorem implies

    f(N)>=2N+1,
    f(N)=2N+1 iff N is a power of2.

For N not a power of2, parity gives f(N)>=2N+3. The R_(k+3) witness has
N=2^k and exactly2N+1 paths. No upper formula for f(N) off powers of2 is
claimed. Endpoint universality therefore cannot be converted into complete
H universality merely by saturating the half bound.

This is the same limitation seen in the symmetry-graft packet: its four-vertex
precursor has (h,R,e,d)=(5,6,2,3), so it lies on the extremal endpoint face;
one bridge leaves that face and creates the rigid arithmetic trap of H=9.
The actual full state, rather than the endpoint alone, controls the operation.

## 4. Exact tests and controls

`04-computation/experiments/session237_endpoint_slack/verify.py` independently
counts paths by subset DP for every labelled tournament of order<=6. It
recognizes the equality shape intrinsically from SCCs and score-ordered
transitive deletions. It checks every positive nonsink root. Separately it
checks every apex mask with l<=9 against all three formulas in(1), including
empty S, singleton S, adjacent and separated pairs. Normal/optimized outputs
and explicit --fail controls are recorded in the final replay.

The finite census tests the implementation; the proofs above provide the
unbounded quantifiers. Chan–Pak completeness remains OPEN.
