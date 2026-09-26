# Exact precision collisions and returning arithmetic cores

**PROVED elementary identities / FINITE-EXACT controls / OPEN global
descent, 2026-09-26.** No novelty claim. This develops the specific
precision failure in section7 of the
[carry-defect boundary note](nextforest_20260926_boundary.md), rather
than treating that failure as an unexplained breakdown.

## 1. Inheritance and a precise change of representation

The closest mechanism is the boundary note's exact lift
`q*2^r+u`, valid until the small core uses the available dyadic precision.
The hostile is27, whose smaller core9 does converge but whose division
budget exceeds the available five bits. The corrected near miss is
equating core convergence with a source-preserving contractive lift.
The least-used sidecar is which summand has the smaller2-adic valuation
when their valuations cross.

The board is the actual odd source, affine coefficient q, precision r,
core u, collision valuation, and canonical-reset clock. The anchor is
preserving one source through a failed shadow; the niche is exact
valuation comparison; the wildcard is a recursive core returning to27.
The [corrected Q2 endgame](collatz_procgen_20260922_q2_endgame.md),
section3.2, already warns that consumed precision can be recreated by
the next arithmetic value. Its different dynamics do not imply any of
the identities below; the warning supplies the hostile question.

Write U for the odd Collatz map. A state is an exact positive odd integer

    y=q*2^r+u,   q,u positive odd, r>=1.                (P1)

The core u need not be below2^r. That distinction is load-bearing.

## 2. The three cases are an exact arithmetic controller

Put `3u+1=2^a v`, with a>=1 and v positive odd. Direct substitution gives

    3y+1=3q*2^r+2^a v.

If a<r, the parenthesis after dividing by2^a is odd, so

    (q,r,u) -> (3q, r-a, v).                         (P2)

If a>r, divide instead by2^r; the remaining low term3q is odd, giving

    (q,r,u) -> (v, a-r, 3q).                         (P3)

This swaps the roles of coefficient and core. If a=r, the two odd
terms meet and there is extra cancellation:

    U(y)=oddpart(3q+v).                              (P4)

Every one of(P2)--(P4) consumes exactly one U step. Both triples in
(P2)--(P3) again satisfy(P1). The failure of the previous shadow is
therefore replaced by an exact transition, not by another chosen source.

At a collision, retain the exact endpoint x in(P4). One possible reset
rule is: if x>1, compute z=U(x), and if z>1 write

    z=2^R+w, R=floor(log2 z), 0<w<2^R,
    new state=(1,R,w).                              (P5)

The arithmetic x->z costs one additional U step; splitting z costs none.
Thus collision followed by this fresh canonical reset costs two U steps
in total. The terminal1 is recorded explicitly. One could instead split
x immediately, but that is a different reset convention and does not
give the same core labels used below.

## 3. A positive boundary lemma, and why its hypothesis can disappear

Call(P1) canonical at the boundary when `0<u<2^r`, with no restriction
on the positive odd coefficient q. If a>=r then

    3u+1 < 3*2^r,

so the only possibilities are:

- **a=r.** Then v is odd and less than3, hence v=1. The equality
  `u=(2^r-1)/3` is integral only for even r. Thus `U(y)=U(q)<y`.
- **a>r.** Necessarily a=r+1 and v=1, with
  `u=(2^(r+1)-1)/3`; here r is odd. Thus `U(y)=3q+2`.
  This is strictly below y for r>=3, but strictly above y for r=1.

The inequalities are immediate: in the first case r>=2 and
`U(q)<2q<q*2^r+u`; in the second, r>=3 gives
`3q+2<8q+5<=q*2^r+u`, while r=1 gives u=1 and
`3q+2>2q+1`.

For q=1, the canonical collision directly reaches1; a canonical swap
reaches5. These are the familiar inverse fibres, now obtained as precise
boundary cases of the controller. They do not classify noncanonical
states.

The hypothesis is not invariant under(P2). In the actual27 example,
`(3,3,7)` is canonical, but its successor `(9,2,11)` is not. Applying
the preceding v=1 conclusion to that successor would be false. The
full exact triple survives; its canonical boundary property does not.

## 4. The27 core reappears after the first collision

The initial canonical reset27->41 gives `(1,5,9)`. The exact controller
then gives the following table; every row uses the same original orbit.

| current state | value y | a=v2(3u+1) | next operation |
|---|---:|---:|---|
| (1,5,9) |41|2| consume ->(3,3,7) |
| (3,3,7) |31|1| consume ->(9,2,11) |
| (9,2,11) |47|1| consume ->(27,1,17) |
| (27,1,17) |71|2| swap ->(13,1,81) |
| (13,1,81) |107|2| swap ->(61,1,39) |
| (61,1,39) |161|1| collision ->121 |

The next reset is121->91, and

    91=2^6+27.                                      (P6)

Thus the original source27 has returned as the *core* of a larger
actual state. It has not returned as an orbit vertex and no integer
cycle has been exhibited. A recursive claim that smaller cores suffice
must track this distinction.

Independently start at91. Its first reset gives

    91->137=2^7+9->103->155->233.                     (P7)

The same core9 now has seven bits of initial precision instead of five.
In particular233 is the actual value of `(27,3,17)` on this orbit.
Equations(P6)--(P7) give a specific arithmetic connection among27,91,
and233; a numerical coincidence is not substituted for it.

The boundary note's section6 already proves the infinite constant-defect
family `n_k=(4^k+17)/3`, k>=3, beginning27,91,347,1371,... . Its
large members have first descent exactly6, despite the37 and28 odd
steps needed by27 and91. Our script independently checks its six-step
identities using r=2k-1; this is a control of that inherited theorem,
not a second claim of discovery. Extending the algebraic source formula
to r=1,3 gives7,11, but those are not the same canonical d=18 boundary
because2^(r+1)<18. Canonical interpretations here start at r>=5.

## 5. Why subtraction of valuations is not a Euclidean rank

The consumed precision in(P2) decreases, but the swap precision(P3)
can increase arbitrarily. Let

    u_k=(4^k-1)/3, q=1, r=1, k>=2.

Then a=2k,v=1, so

    (1,1,u_k) -> (1,2k-1,3),
    U(u_k+2)=2^(2k-1)+3.                             (P8)

The precision grows from1 to2k-1 in one actual step. For k>=3 these
are exactly the E=12 hostile sources in section4 of the boundary note.
That note repairs them with a selected descent, not with monotonicity
of precision. Thus the new controller cannot be declared a Euclidean
algorithm merely because its branches contain differences of exponents.

This is a useful interface to the excursion forest: retain each actual
triple, which branch occurred, and the exact collision endpoint; a
collision or swap cannot be silently charged as a completed descent.
An indefinitely open controller run remains possible unless an additional
inequality excludes it. Its impossibility is not assumed here.

| Source -> target | Preserved | Lost without sidecar | Required hostile |
|---|---|---|---|
| Actual y ->(q,r,u) | exact value and next U step | representation is not unique | keep its chosen chronological decomposition |
| Small-core shadow -> three-case controller | exact arithmetic through precision failure | simple monotone consumption |(P8) recreates unbounded precision |
| Canonical boundary -> v=1 lemma | selected descent in its stated cases | canonicality after a consume step |(3,3,7)->(9,2,11) |
| Collision -> fresh core | actual endpoint and one explicitly charged reset | source/core roles and clock |161->121->91, whose core is27 |

The next question is an amortized estimate on *actual source-preserving
resets*, with both coefficient/core roles and the carry-defect boundary
retained. The current identities supply the exact transitions on which
to test it. They supply no well-founded rank by themselves.

## 6. Reproduction

Run the [independent script](../../04-computation/experiments/nextforest_20260926_precision.py)
normally and with `-O`. Its [output](nextforest_20260926_precision.out)
records all20480 triples with q odd1..63,u odd1..127,r1..10;256
canonical boundary cases; capped96-step runs from all1023 odd sources
3..2047; the displayed27/91 paths;127 precision-recreation examples;
and the inherited family's six-step controls through r257. The capped
runs do not presume that a collision or termination must occur.
All checks raise explicit exceptions and remain active under `-O`.
