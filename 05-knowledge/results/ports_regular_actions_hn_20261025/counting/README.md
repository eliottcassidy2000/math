# A strictly smaller fixed-target recurrence for unrestricted A000568

2026-10-25. **PROVED:** the fixed-degree gauge quotient, its exact recurrence,
and its minimality for preserving every normalized individual continuation.
**FINITE-EXACT:** independent controls and response-graph sizes. This is a
new quotient of the inherited Burnside computation, not a Burnside-free
formula, a bounded-state scalar recurrence, or a tournament-isomorphism
classification. No historical-priority claim.

## 1. Inheritance and the coordinate that can be removed

The closest mechanism is [the future-response transfer](../../defect_boundary_20261020/boundary.md),
with fixed-permutation formula from
[THM-516, Burnside core kernel phi reframe](../../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md).
The [packed engine](../../tournament_enumeration_20261021a/README.md) and
[degree-local normalization](../../degree_local_tournament_enumeration_20261021/session22-enumeration/README.md)
accelerate arithmetic while retaining the same response graph. We change
that graph for the task of computing one fixed coefficient A000568(N).

Canonical hostile: equal degree and cycle count do not determine future
gcd interactions. Corrected near miss: a full divisor profile is invertible,
so is not by itself compression. Least-used sidecar: every completed
continuation has the same **actual permutation degree**, not merely a
degree bound. The concept board is exact degree / future responses /
affine gauge / signed substitution excess / literal ports lost by counting.

Process odd cycle lengths in decreasing order. Before processing ell,
all past cycle lengths c exceed ell. Let

    u=sum(past lengths), m=number of past cycles, R=N-u,
    f(d)=sum_(past c) gcd(c,d).

The old key retains u and all affordable future f(d), for odd
1<=d<=min(ell,R). Instead retain

    D(d)=d*m-f(d),                 D(1)=0.             (1)

If the old accumulated prefix weight is w, store

    w_tilde = 2^(R*m) w.                               (2)

Different histories with the same u,D now have their w_tilde weights
summed. Actual cycle count m is not part of the new key.

## 2. Why the quotient is exact, and its precise minimality

For any future odd partition mu of exactly R, its interaction with the
past is

    sum_(d in mu) f(d) = R*m - sum_(d in mu) D(d).      (3)

Thus (2) absorbs everything lost by removing m. The remaining multiplier
depends only on D, the future partition, and its own intrinsic weight.
The equality is for every individual continuation, so it remains true
after summing all of them. This is an exact weighted-state quotient.

It is also the coarsest such quotient after normalization by (2), at a
fixed R and stage. Every future uses odd lengths at most ell, including1.
The all-ones continuation has normalized interaction1. Comparing it with
the continuation consisting of one d-cycle and R-d ones reads exactly
2^(-D(d)), for each affordable odd d. Therefore two histories have the
same normalized multiplier for every continuation iff every retained
D(d) agrees. This proof uses decreasing stages: if1 had already been
processed, those probe continuations need not exist.

This is minimality for the stated observable of every individual future
partition. It does not forbid a different representation, cancellation
method, or aggregation preserving only a final sum from compressing more.

Two strict mergers illustrate that the reduction is real:

| Past histories | Degree | Old (f(1),f(3)) | New (D(1),D(3)) |
|---|---:|---|---|
| (23), (5,9,9) | 23 | (1,1), (3,7) | both(0,2) |
| (27), (9,9,9) | 27 | (1,3), (3,9) | both(0,0) |

These histories are legal before stage3, and the first merger is present
at target N=26. Their old responses differ by a multiple of the degree
function d. Conservation of the completed degree makes that difference
a scalar weight rather than a state distinction.

More generally, continuation count vectors x obey sum d*x_d=R.
Their difference vectors span the hyperplane orthogonal to (d).
The discarded response direction is exactly this annihilator. This is
a concrete affine quotient, not an analogy based on matching counts.

## 3. The unrestricted recurrence

For odd ell define F_ell(R,D) by summing the normalized weights of all
odd partitions mu of R with parts at most ell:

    F_ell(R,D)=sum_mu 2^[e(mu)-sum_(d in mu)D(d)] / z_mu,

where e(mu)=sum(d-1)/2+sum_(unordered cycle pairs)gcd(d,e), and
z_mu=product_d d^(m_d) m_d!. These are the inherited fixed-permutation
quantities. Choose k copies of the next length ell, put R'=R-ell*k, and
update each affordable smaller odd coordinate by

    D'(d)=D(d)+k*(d-gcd(ell,d)).                       (4)

Then exactly

    F_ell(R,D)
      =sum_(k=0)^floor(R/ell)
         2^[ k*(R-D(ell)) - (ell*k*k+k)/2 ]
         /(ell^k*k!) * F_(ell-2)(R',D').              (5)

If ell>R, only k=0 is possible; simply skip to the next stage without
reading or storing the unaffordable coordinate D(ell).

The exponent is integral because ell is odd and k(k+1) is even. It can
be negative, so the displayed recurrence uses exact rational arithmetic.
At ell=1 one may directly use

    F_1(R,D)=2^binom(R,2)/R!,                         (6)

since D(1)=0. Set F_ell(0,D)=1 and F_ell(R,D)=0 for R>0,ell<1. With
L the largest odd integer at most N,

    A000568(N)=F_L(N,0),     A000568(0)=1.             (7)

**Proof of the transition.** The k new ell-cycles contribute
k(ell-1)/2+ell*binom(k,2) internally and -kD(ell). Their interaction
with the smaller future is absorbed by (4), at the price of adding kR'
to the exponent. Their denominator is ell^k*k!. Summing gives (5).
Equivalently, substitute (1)--(2) into the inherited transfer; every
occurrence of m cancels. Equations(6)--(7) follow from the single possible
all-ones partition and the original Burnside sum. QED.

Only affordable coordinates are stored. Because past c>ell,

    0<=D(d)<=m(d-1)<=u*(d-1)/(ell+2)<N

for d<=ell. Thus the old N-sized packed-coordinate radix is still
sufficient; the gauge does not require wider response digits. The supplied
reference implementation uses tuples and Fractions to make the proof
transparent. It is not a new integer or packed performance benchmark.

## 4. Exact graph reduction, not an unqualified speed claim

`check.py` compares two otherwise identical fixed-target tuple/Fraction
engines: inherited f responses, and the new D gauge. Both force completion
at stage1, use the same affordable-coordinate mask, and return only A(N).
They do not compare a one-coefficient run with an all-coefficients engine.
The measured quantities below are exact graph counts, not timings.

| N | Old input states, summed over stages | New | Old transitions | New |
|---:|---:|---:|---:|---:|
| 30 | 267 | 212 | 531 | 472 |
| 40 | 598 | 470 | 1,403 | 1,229 |
| 60 | 2,178 | 1,720 | 6,563 | 5,499 |
| 80 | 6,349 | 5,119 | 22,355 | 18,109 |
| 100 | 16,230 | 13,444 | 62,766 | 50,044 |

At100 this removes about17% of stage states and20% of contributions.
The transition graph is smaller, unlike the inherited arithmetic-only
normalization. No asymptotic speedup or polynomial state bound follows.
The new rational weights may cost more individually. Reusing computations
between targets N and N+1 is also a separate optimization problem.

## 5. Why clone-allocation excess cannot be silently included

The pair-free-core substitution uses excess degree that is not actual
permutation-cycle degree. If the total remaining substitution degree is R
but only a of it is actual cycle degree, (3) instead contains a*m. Scaling
by2^(R*m) leaves2^((a-R)m), so m has reappeared.

For the merger (27) versus(9,9,9), take R=4, one actual future fixed point
and three substitution-excess units. Its normalized interaction is1/8
for the first prefix and1/512 for the second, even though their D values
agree. This is a minimal-purpose hostile to extending (5) unchanged to
core mode. A repair must retain actual remaining degree or a grading that
tracks the excess. No clone merger or canonical-owner identity is inferred
from the ordinary counting quotient.

The source histories already forget individual tournament arcs, labelled
ports, and automorphism identities. Our map preserves their exact ordinary
counting contribution, including classes with clones; it does not recover
those forgotten structures. The [canonical owner forest](../../ports_owners_f1_20261025/anchor.md)
continues to supply a separate lossless structural encoding when they are
needed. Its growing rank is not made bounded by (5).

## 6. An exact 30030-offset horizon example

The two histories

    A=(30030-17,30030+17)=(30013,30047),
    B=(30030-1,30030+1)=(30029,30031)

both have degree60060 and two cycles. The first three displayed numbers
are prime, while30031=59*509. Direct trial division controls these facts.
Hence for every odd future length d<=57, both histories have f(d)=2
and D(d)=2(d-1). Their responses therefore agree at every actual remainder
R<=58. At d=59 they separate:

    f_A(59)=2, f_B(59)=60;
    D_A(59)=116, D_B(59)=58.

Thus this particular offset coincidence has a precise future-response
interpretation, with a sharp first differing odd coordinate. It is an
inherited-style finite-horizon merger, not an additional example of the
new affine gauge: their old responses already agree before59. The degree
and aggregate weights still must be retained. In particular the single
cycles30013 and30047 cannot be merged as full states merely because their
small responses agree; their degrees differ34.

This provides an actual source/target map, observable and failure boundary
for the user's offset theme. No plane-coloring or geometric quotient
theorem follows from these integer factors alone.

## 7. Verification and strongest next target

Run `python check.py` and `python -O check.py`. They use explicit exceptions.
Independent unmerged odd partitions verify every N0..27; literal permutation
actions on unordered edges verify N0..7 without using the gcd exponent.
All individual normalized future multipliers are compared on the declared
small complete prefix/future universe. The two larger state engines agree
at30,40,60,80,100. The strict mergers, invalid-input controls, core-excess
hostile and30030 horizon boundary are checked directly.

The next useful target is to combine this genuine quotient with the existing
odd-factorial integer transfers, then retain one explicit excess grading
for pair-free allocation. It is not to drop a still-observable port or
degree. Formula(5) is unrestricted in N but grows in state; the earlier
non-P-recursiveness obstruction to a fixed-order polynomial-coefficient
linear scalar recurrence remains valid.
