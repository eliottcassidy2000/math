# Tiny physical experiments classify every same-root fusion continuation

**Session281. PROVED from the inherited rooted-fusion law.** The linear algebra
and recognizable-series principle are elementary; no claim of novelty for
Hankel-rank theory is made. The contribution is an explicit small physical
observation basis, exact repair iff tests, and a concrete boundary when the
operation language is enlarged to rerooting. Chan–Pak remains OPEN.

## 1. Change the question: what information can any continuation consume?

A successful receiver repair must preserve future behavior, not just today's
Hamiltonian count. Session280 exhibits the distinction concretely: its repaired
three-count bank has common `(h,V,e)`, whereas the widened six-count union loses
V. Instead of guessing more sidecars, ask for a finite **complete experiment
basis** for the actual operation language.

The distant source is
[THM-2810, factorial Hankel faithfulness](../../../../01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md):
there, multiplication contexts detect every polynomial and prohibit finite
memory. Here, rooted fusion has an inherited finite matrix carrier, so the
opposite outcome is possible. The common map is `object -> all continuation
readouts`. The target predicate is equality of those readouts; arbitrary graph
operations are not among the contexts. A second source is the
[integral-lifting discussion](../../automata_arithmetic_geometry_session_20261019/INTEGRAL_LIFTING_CONNECTIONS.md):
keep the integer normalization when recovering coordinates from observations.
Our inverse uses halves of even differences, not division modulo2.

For a pointed tournament T write H for its Hamiltonian path count, h for the
count after removing its root, Q for the ordered two-path-cover count of that
deletion (empty slots allowed), U,V for the weighted root-start/end masses,
and d,e for their odd-minus-even path-length masses. The inherited
[Session264 fusion law](../../session264_constructive_20261111/fusion/PROOF.md)
is

    M(T) = [ h  U-h  d-h  H-h ]
           [ 0   Q    0   V-h ]
           [ 0   0    k   e-h ]
           [ 0   0    0    h  ].                         (1)

Physical fusion `L*T` identifies roots and orients every nonroot L vertex to
every nonroot T vertex. Matrices multiply. For a nontrivial tournament k=0;
for the singleton I, M(I) is the identity and k=1. No new derivation of (1)
is claimed. Let `lambda(M)=M11+M14`, so H(T)=lambda(M(T)).

## 2. Seven Hamiltonian counts recover the full nontrivial carrier

All row masks below use root0. Define these actual probes:

    I=(0)           singleton,
    S=(2,0)         root-source edge,
    E=(0,1)         root-sink edge,
    EE=(0,1,3)      transitive triple with root sink,
    C=(4,1,2)       directed triangle.

For any nontrivial T measure

    p0=H(T),        p1=H(T*E),     p2=H(T*EE),
    p3=H(T*C),      p4=H(S*T),     p5=H(S*T*E),
    p6=H(C*T).                                            (2)

Every experiment adds at most TWO vertices. Direct substitution in (1) gives

    p0=H,
    p1=H+U-d,
    p2=H+3U-2h-d,
    p3=H+2U,
    p4=H+V-e,
    p5=H+Q+U+V-d-e,
    p6=H+2V.                                              (3)

Therefore the inverse is explicit:

    H=p0,            U=(p3-p0)/2,       V=(p6-p0)/2,
    d=p0+U-p1,       e=p0+V-p4,
    h=(p1+p3-p2-p0)/2,
    Q=p5+p0-p1-p4.                                        (4)

All divisions in (4) are exact over the integers for actual observations.
The coefficient matrix of (3), on coordinates `(H,h,Q,U,V,d,e)`, has determinant
of absolute value8. One must not reduce observations mod2 and then try to
invert: every physical p_i is odd, and their parity shadow has lost the data.
The theorem makes no assertion that an arbitrary integer tuple passing the
congruences is physically realizable.

**Theorem1 (complete finite experiment basis).** For two nontrivial pointed
tournaments T,T', the following are equivalent:

1. the seven observations (2) agree;
2. all seven rooted coordinates in (1) agree;
3. `H(L*T*R)=H(L*T'*R)` for every actual pointed L,R.

Proof: (1)=> (2) by (4); (2)=> (3) by the inherited matrix law; (3)=> (1)
because the seven experiments are such contexts. This is an exact behavioral
congruence, not an isomorphism classification. Arbitrarily many common left
and right fusions collapse to L and R by associativity. ∎

For singleton T, p5 contains the additional k=1. The order-one object is
recognized separately. The seven-probe formula is NOT applied blindly there.
If order is not supplied, an eighth experiment `H(S*T*EE)` permits recovery
of k as well, using the full eight-entry algebra. Its observation determinant
has absolute value16. The verifier includes this boundary.

## 3. Exact iff tests for repair and for central count gates

Let `delta=H(T')-H(T)`, with no restriction on its sign.

**Theorem2 (universal two-sided count translation).** The following are
equivalent for nontrivial T,T':

* all seven differences `p_i(T')-p_i(T)` equal delta;
* every rooted coordinate other than H agrees;
* for all actual pointed L,R,

    H(L*T'*R)-H(L*T*R)=delta h_L h_R.                    (5)

These are also equivalent to the same difference being delta for ALL contexts
with h_L=h_R=1. All probes in (2) are of that type. For the forward implication,
(4) recovers zero coordinate differences except H. Conversely
`M(T')-M(T)=delta E14` and `M(L)E14M(R)=h_Lh_R E14` prove (5).
The scaling law is inherited; the finite iff experiment test is the point here.
No subtraction of matrices is implemented as a tournament operation.

For receiver-only repair use SS=(6,4,0), a transitive triple with root source,
and the four LEFT observations

    q0=H(T), q1=H(S*T), q2=H(SS*T), q3=H(C*T).

Then

    q1=H+V-e, q2=H+3V-2h-e, q3=H+2V,
    V=(q3-q0)/2, e=q0+V-q1,
    h=(q1+q3-q2-q0)/2.                                  (6)

**Theorem3 (complete renewable receiver test).** Pairwise differences of all
four q_i equal delta iff T,T' share `(h,V,e)`. This is equivalent to count
difference delta under every left context of deletion count1. With arbitrary
left L, the corresponding difference is h_L delta. The equivalence follows
from (6) and multiplication in (1). This exactly characterizes the continuing
receiver used in Session280, not merely a sufficient list of coordinates.

For a finite bank, check against one representative. If (6) passes, every
sequence of common left unit-deletion fusions preserves all its H spacings.
If it fails, one of the three explicit nonidentity probes already exposes
failure. Thus a supposed infinitely renewable repair always has a decisive
counterexperiment adding at most two vertices.

## 4. Sharp rational-linear memory and the two different dimensions

The seven transition coordinates and the four continuation-state coordinates
are different objects. The row and column bases from actual contexts are

    L = [1 0  0 0]      R = [1  1  1 3]
        [1 1 -1 0]          [0  1  3 2]
        [1 3 -1 0]          [0 -1 -1 0]
        [1 2  0 2],         [1  1  1 1].

L consists of first rows for I,S,SS,C; R consists of `M(R)(e1+e4)` for
I,E,EE,C. Their determinants are4 and−4. The physical count table is

    LR = [1 1  1 3]
         [1 3  5 5]
         [1 5 11 9]
         [3 5  9 9],       det(LR)=−16.                  (7)

Any rational linear word representation of the fusion-count series has
Hankel matrices of rank at most its state dimension. Hence (7) proves that
no such representation of dimension<4 works for all physical letters and
words. The inherited four-state carrier attains this bound. This is a precise
minimal linear-memory theorem for this operation language, not an all-purpose
lower bound for algorithms or nonlinear encodings.

Similarly the seven-coordinate nontrivial transition space really has
rational dimension7. Seven actual pointed tournaments of orders2,3,4 have
coordinate rows

    (1,1,2,1,2,1,0), (1,1,2,2,1,0,1),
    (1,1,4,1,4,1,0), (1,1,4,2,2,0,0),
    (3,1,4,3,3,1,1), (1,1,4,4,1,0,1),
    (1,1,8,2,4,0,0).

Their determinant is16. Thus fewer than seven scalar linear observations
cannot recover the arbitrary transition operator by a rational LINEAR
inverse on its full span. This does not rule out fewer nonlinear queries on
restricted physical families, nor contextual compression after fixing a
particular seed or alphabet. The exact matrices and witnesses are in the
verifier.

## 5. The entire fusion orbit has an exact seven-count content certificate

Define `c(T)=gcd(H,h,Q,U,V,d,e)>0` for nontrivial T. Since H and h are odd,
c(T) is odd. Then

    gcd { H(L*T*R) : actual pointed L,R }
       = gcd(p0,p1,...,p6) = c(T).                       (8)

The divisibility from c(T) follows from the integer matrix law. Conversely,
if g divides every p_i then g is odd (it divides p0=H), so2 is invertible
modulo g. The exact reconstruction (4) shows that g divides all seven
coordinates. This proves (8), including composite and prime-power content.
No asymptotic distribution or finite-horizon extrapolation is involved.

Consequently any finite seed bank with all c(T)>1 cannot supply any sufficiently
large prime Hamilton counts through arbitrary common-root fusion contexts:
a prime output must itself be one of those finitely many prime contents.
This does not apply to interspersed root changes or core surgery, and c(T)=1
is only a necessary gate, not a target-coverage theorem.

The prior [Session237 fixed-port content theorem](../../session237_chanpak_synthesis_20261109/symmetry_grafts/PROOF.md)
already established the corresponding four-test statement for O/Z-style
fixed-port extensions. Here the operation language is arbitrary left AND right
rooted fusion, and the exact full-carrier content is read by seven tiny probes.
The elementary gcd/integer-reconstruction method is inherited. The larger
language and its sufficient and necessary observation list are stated explicitly.

## 6. The exact quotient is unsafe as soon as rerooting is allowed

There are two literal six-vertex pointed tournaments

    A=(0,41,3,5,15,29), B=(0,33,3,7,15,29)

with the same entire packet

    (H,h,Q,U,V,d,e)=(9,9,84,9,84,9,0).                   (9)

Thus no same-root fusion experiment can distinguish them, even with arbitrary
large common physical contexts, or policies whose choice of common contexts
depends only on earlier fusion-count observations. Nevertheless their multisets of deletion
counts over ALL possible roots are

    A: {1,3,3,5,5,9},       B: {1,1,5,5,5,9}.             (10)

The first has a root of deletion count3; the second has none. They cannot be
isomorphic as unrooted tournaments, since the displayed deletion menu is an invariant. Their current
root is a sink; adding this sink to distinct count9 cores is what makes the
fusion observation quotient forget the alternative root resources.

The prior [Session237 symmetry-graft work](../../session237_chanpak_synthesis_20261109/symmetry_grafts/PROOF.md) already warned that fixed-port response twins can separate after a root change. Our claim is the explicit full-fusion-packet witness (9), not novelty for that general warning.

This is a concrete counterexample to extending Theorem1 to a grammar with
root changes. It is not a failure of Theorem1. In particular, it is unsafe to
merge search candidates solely by the full seven-coordinate packet and then
allow arbitrary rerooting, even if every current fusion continuation was
proved equivalent. A root-menu sidecar or finer actual-object storage is
necessary. No sufficiency of the deletion-count multiset alone is claimed.

## 7. Consequences for constructive search, and limits

The practical object is an operation-relative repair test:

* same-root fusion only: seven counts give a complete quotient;
* future unit left contexts only: four counts give the exact receiver test;
* root changes or substitutions: neither quotient is promised complete.

This separates two difficulties previously easy to conflate. Finding a
candidate that passes the complete receiver test is a finite constructive
problem; proving that enough such candidates exist at arbitrary absolute
counts is still an infinite supply problem. The tests do not solve that
supply problem, do not create new Hamiltonian counts by taking differences,
and do not prove Chan–Pak. They make every proposed repair's continuing-state
obligation exact and locally falsifiable, while the hostile (9)–(10) prevents
an otherwise legitimate compression from erasing the route out of a stalled
fusion language.

The provenance search found the central transport and C41/C43 switch already
in Session276 and earlier sources; neither is promoted as new here. The new
packet is a physical observation/repair interface and explicit reroot boundary.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. Literal Hamiltonian counts
use a backward subset recurrence separate from the inherited forward profile
helper. The finite controls include all labelled pointed tournaments of
orders2 through5, the two witnesses in (9), the singleton boundary, the exact
rank determinants, and the three/six-bank receiver distinction. These checks
support the explicit inputs; the all-context conclusions follow from (1)–(7).

## Prepared experiment basis (restored extension)

Let C1 be the rooted word10 block, with row masks(2,4,1). The carrier law gives

    H(C1*T)=H+2V, h(C1*T)=h,
    V(C1*T)=4V-h, e(C1*T)=h.

The three counts r0=H(T), r1=H(C1*T), r2=H(C1*C1*T) therefore recover

    V=(r1-r0)/2, h=(5r1-4r0-r2)/2.

For two inputs whose count difference is delta, all three experiment differences equal delta if and only if their h,V agree. In that case their C1-prepared receivers agree in h,V,e, and all subsequent common left fusion contexts preserve that offset. This is an exact iff for this preparation, not for arbitrary graph operations.

The raw row-mask pair(14,8,2,4,15) and(12,5,8,2) has(H,h,V,e)=(3,3,6,0) and(5,3,6,2). The two prepared objects have(15,3,21,3) and(17,3,21,3). The gap was not initially renewable under every left context, but becomes renewable after preparation.

An adaptive-context statement requires a common context or a policy driven only by the prior count observations. A policy that inspects the full hidden graph may distinguish equal carriers; no equivalence for that stronger policy is asserted.

## Exponentially many fusion aliases (restored extension)

Remove root0 from the six-vertex witnesses in section6 to obtain strong five-vertex cores A0,B0. Each has H=9 and ordered two-cover count Q=84. For any binary word of length m>=1, take its ordered SCC concatenation of these cores and append a single sink root. Strong-component order is canonical, and A0,B0 are nonisomorphic, so these give2^m nonisomorphic rooted tournaments.

All have the same carrier

    H=h=9^m, Q=84^m, U=d=9^m, V=84^m, e=0.

Indeed path counts and two-cover counts multiply across ordered SCCs. At the sink, root-start paths have only the singleton root, while root-end cover sums range over every ordered two-cover of the deletion. Its alternating cover sum is zero for a nonempty deletion. These facts give the displayed coordinates directly.

Deleting a vertex from one core multiplies that core's deletion count by9^(m-1); deleting the sink gives9^m. In A0 the deletion counts are1,3,3,5,5, while B0 gives1,1,5,5,5. Thus the number of occurrences of3*9^(m-1) distinguishes m+1 menus according to how many A0 factors occur. Equal fusion carriers can conceal unboundedly many different rerooting resources. This is an all-order construction, not an extrapolation of the six-vertex pair.
