# A physical Farey equation, with a finite partner bound

**PROVED:** the arithmetic equivalence, entry guards and finite partner
bound below. The two-module surgery and its exact count formulas are
proved in the parent packet; this note extracts their ordered Farey
interface. The literal order-eight witness was found by the root agent's
response-bank search. No infinite family of physical partners or complete
ternary compiler is asserted.

## 1. The old operation and its new carrier

[THM-2596,§1](../../../../01-canon/theorems/THM-2596-modular-free-factor-farey-gram-owner-cocycle.md)
retains an ordered unimodular flank U=[u v]. Its two Stern children are
[u,u+v] and[u+v,v]. These are actual column operations preserving the
determinant. Its Gram/owner sidecar is essential when the objective
depends on more than the address. The new carrier below is a pair of
actual tournament observations; it is not identified with every integer
flank or every Stern word.

For nonempty tournaments M,N put

    h=H(M), p=P_M(2), j=H(N), q=P_N(2).

The parent packet's same-name three-apex-arc surgery changes the count
from

    A=12pq+2pj+6hq−hj

to Z=A+4D, where D=pj−hq. The prepared root has readiness defect hj.
Its actual NR exit has count

    Y=3Z+2hj=3A+2e,
    e=6D+hj=j(6p+h)−6hq.                              (1)

Thus an exact ternary±2 target relative to the ORIGINAL source count A
occurs precisely when

    det[(6p+h,h)^T,(6q,j)^T]=e=±1.                    (2)

Unlike asking the surgery itself to be count-neutral, (2) permits a
controlled nonzero count change which cancels the later readiness term.
This is the new useful role of the determinant. No intermediate count
is asserted equal to A.

## 2. Entry guards and arithmetic addresses

Set a=6p+h and b=6h. For a fixed physical M, equation(2) is

    a j−b q=e, e∈{−1,1}.                              (3)

It has integer solutions iff gcd(a,b)=1. Since h is odd, this condition
is equivalent to

    gcd(h,p)=1 and 3∤h.                                (4)

Indeed a is coprime to2 automatically; a≡h modulo3; and the common
prime factors of a and h are exactly those of6p and h. These observations
give both directions. Any actual partner also necessarily has
gcd(j,q)=1 and3∤j. Reducing(3) modulo3 gives the signed guard hj≡e(mod3).

Given one integer solution(q0,j0) for a chosen sign, every solution is

    (q,j)=(q0+a t,j0+b t), t∈Z.                       (5)

The equation already forces j odd because a is odd and b even. The
condition q even selects one parity of t, since a is odd. It therefore
leaves an arithmetic ray of formally valid positive observations, with
increment(2a,2b), once t is large enough. These are only observation
pairs; they are not tournaments unless an actual module realizes them.

The order of the modules matters. Swapping M,N sends D to−D and

    e ↦ 2hj−e.                                         (6)

For hj>1 a successful±1 pair cannot remain a successful±1 pair in the
opposite order. Unordered profile sets lose the exact output digit.

## 3. The infinite arithmetic line has only finitely many physical partners

Every actual nonempty N satisfies q≥(|N|+1)j, by the HP-cut injection.
Equation(3) instead gives

    q/j=p/h+1/6−e/(6hj).

Since j≥1, any physical partner has the explicit bound

    |N|≤floor(p/h−5/6+1/(6h)).                         (7)

For e=+1 the slightly smaller upper bound floor(p/h−5/6) is valid.
Consequently a FIXED M has only finitely many actual partners N, up to
isomorphism, even though(5) contains infinitely many positive integer
pairs. One can decide its whole partner question by an ordinary tournament
inventory only through the explicit order bound(7), retaining (H,P2).

No source-preserving padding escapes this bound: it leaves H(N) fixed but
multiplies P_N(2) by a power of2, whereas(3) fixes their allowed ratio to
a bounded interval. A family attempting to make this exact-target
compensation recur at unbounded scale must vary M as well, or change the
response interface. This is a finite source statement, not a negative
claim about those enlarged constructions.

### Exact selector for independent transitive padding

There is also a complete finite test when BOTH modules vary only by
transitive ordinal padding. Such padding of lengths t,s≥0 leaves h,j
unchanged and replaces p,q by2^t p,2^s q. Fix the ordered core pair,
assume hj>1, and choose the desired sign e∈{−1,1}. Set

    P=pj, Q=hq, C=(e−hj)/6.

If C is not an integer, there are no solutions. Otherwise C is nonzero,
and the exact target condition is

    2^t P−2^s Q=C.                                    (8)

For m=min(t,s), necessarily0≤m≤v2(|C|). For each such m test precisely
these two alternatives:

    t=m, s=m+d:  2^d=(P−C/2^m)/Q;
    s=m, t=m+d:  2^d=(Q+C/2^m)/P.                     (9)

Accept an alternative iff its right side is a positive integral power
of2, including1; d is then its unique nonnegative exponent. Deduplicate
the t=s case. This algorithm finds ALL actual padding pairs: necessity
is division by2^m in(8), and sufficiency is direct substitution. Every
accepted pair has literal modules from the two supplied cores and their
transitive padding. At most2(v2(|C|)+1) candidates survive before
deduplication; no infinite dyadic search is needed.

The exceptional case hj=1 has h=j=1, so both modules are transitive.
For e=−1, C is not integral. For e=1, C=0 and the condition is
2^t p=2^s q. Here p,q are powers of2, and infinitely many equalized
paddings exist. This is the trivial transitive family, separately scoped
from the nontrivial hj>1 result.

Thus any finite bank of nontrivial ordered core pairs, even allowing
unbounded independent transitive padding, supplies only finitely many
exact compensated±2 sources. The conclusion does not cover new strong
cores or operations changing their two-slot observations.

## 4. A concrete nontrivial physical pair

The root agent found, in the complete inherited response bank through
order8,

    M=(4,9,82,21,35,143,187,31),      (h,p)=(181,3424),
    N=(6,12,48,21,131,219,159,15),    (j,q)=(179,3416).

The rows are outgoing-neighbor bit masks. Exact arithmetic gives

    D=−5400, hj=32399, e=−1,
    A=145259777,
    Z=145238177,
    Y=435779329=3A−2,
    gcd(Y,Z)=1.                                        (10)

The move reverses three old apex arcs; NR then adds two vertices. It
retains the module interiors and the designated root star, but it is not
an induced extension on all retained old vertex names. Its balanced
output is primitive; no prime value is assumed in(8).

The partner bound for this M is18. Thus an unlimited search along its
integer Bezout ray is unnecessary: physical partners have bounded order.
The current order-eight bank is not asserted complete through18.
The exact padding selector(9) returns only(t,s)=(0,0) for the negative
digit of this displayed pair, and no positive-digit pair. So independent
transitive padding does not renew this particular compensation.

The root search found one nontrivial±1 profile pair among its1922
distinct profiles, besides transitive trivial cases. That is a scoped
finite report, not a uniqueness theorem for all modules. The independent
checker in this lane recounts the two displayed (H,P2) profiles and
checks(10); it does not repeat the broad profile search.

## 5. What has been transferred

Source: ordered unimodular Farey columns. Target: two physical module
observations in a three-arc surgery followed by NR. Map: columns
(6p+h,h) and(6q,j). Preserved predicate: determinant±1 is exactly the
desired final correction±2. Destroyed information under a raw rational
slope: module existence, path-cover multiplicity, labels and allowable
old-arc edits. Required sidecar: actual rows and the all-module response
formula. Cheapest hostile: the infinitely many integer solutions(5)
cannot be treated as an infinite physical family, by(7).

This differs from a Stern-tree analogy based only on binary branching.
It supplies an exact determinant equation, an actual successful pair,
and a finite physical bound for each fixed source. It does not yet
supply the changing-source rule needed for a complete compiler.
