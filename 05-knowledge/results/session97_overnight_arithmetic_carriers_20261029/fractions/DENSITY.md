# Arithmetic index versus the positive real window: an exact one-third law

Session97 third phase, 2026-10-29. **PROVED given the CITED rank/generator
input of `FULL_IMAGE.md`**, plus elementary real elliptic-curve topology.
This note has no new numerical density estimate and no effective error
bound. It does not infer a density from a finite census.

## 1. The counted objects and the exact claim

Use the complete carrier from `FULL_IMAGE.md`:

    Q_(n,e)=nH+eU on E',  e∈{0,1},
    q(Q)=−Z(Q)/196−1,
    I=(1/6,(29+5sqrt(65))/392).

For a fixed e, let A_e(N) count integers `1<=n<=N` for which `q(Q_(n,e))`
lies in I. These are rational quotient points with three distinct
positive real cubic roots; rationality of those roots is not assumed.
Let L_e(N) count the same points that actually lift to rational fruit
triples. Then

    lim_(N→∞) L_e(N)/A_e(N) = 1/3.              (1)

Equivalently, two thirds of the positive-real quotient candidates fail
rational reconstruction. The exact mod3 classifier of `FULL_IMAGE.md`
rejects them without recovering the integer source index. Nevertheless
that index and the root labels remain necessary for reconstruction.
The limit is about this explicitly specified ordering, not uniform
sampling of rational numbers q by numerator/denominator bounds.

The same limit holds after combining both values of e, or allowing both
signs of n with duplicate conventions fixed consistently. It also holds
in the natural canonical-height ordering on E':
`hhat(nH+eU)=n² hhat(H)`, so a height cutoff is just an index cutoff
up to integer rounding. Canonical height, not naive q height, is the
ordering being asserted. The latter generally needs a separate proof.

## 2. Real components and the positive window

The real roots of E's target cubic are

    0, (-421−5sqrt(65))/2, (-421+5sqrt(65))/2.

Thus E'(R) has two components. H=(-196,196) lies on the bounded component;
U=(0,0) lies on the identity component. The positive window q∈I lies
on the bounded component. Consequently only odd n can contribute to
A_e(N).

Translate that bounded component to the identity component, a circle.
The positive window becomes a finite union of open arcs, with positive
Haar measure μ and finitely many boundary points. Its measure does not
depend on translating by the fixed torsion point eU. Positivity of μ
is immediate from the nonempty real interval and continuity; the exact
35H certificate is also an explicit rational witness in the window.
We do not evaluate μ numerically.

H has infinite order, so6H has infinite order in the identity component.
In a circle group an infinite-order element is an irrational rotation.
For each r∈{1,3,5}, the progression

    Q_(6k+r,e) = rH+eU + k(6H)

is equidistributed on the bounded component. All three progressions
therefore meet the positive window with the same asymptotic frequency μ.
In particular

    #{n<=N: n=r mod6, q(Q_(n,e))∈I} = μN/6+o(N).

Summing gives `A_e(N)=μN/2+o(N)`. The exact rational image criterion is
`3|n`, so among odd residues only r=3 lifts. Hence
`L_e(N)=μN/6+o(N)` and(1) follows. This does not give independence at a
finite cutoff: the three remainder terms can differ substantially.

For clarity, the equidistribution step does not need a deep analytic
estimate. For irrational rotation α and every nonzero integer k,

    (1/M) sum_(j=0)^(M−1) exp(2πikjα)
      = (1−exp(2πikMα))/(M(1−exp(2πikα))) → 0.

Trigonometric approximation gives convergence against continuous
functions; upper/lower continuous approximations of the finitely many
arc boundaries give convergence for the window indicator. This proves
the required limit, but does not make its convergence rate uniform or
effective from the present data.

## 3. The general clock/coset compatibility rule

The mechanism extends beyond this curve. Suppose an indexed orbit nP
visits components labelled by `n mod c`, and every progression with
fixed component is an irrational rotation on a circle. Select an open
window of Haar measure μ>0 in component a mod c, with finite boundary.
Suppose the arithmetic acceptance rule is `d|n`. Put `g=gcd(c,d)`.
Then the asymptotic accepted share within the window is

    0,       if g does not divide a;
    g/d,     if g divides a.                     (2)

Indeed the simultaneous congruences are inconsistent in the first
case. In the second they occupy one residue modulo lcm(c,d), and
irrational rotation supplies the same μ on that progression as on
the full component progression. The counts are respectively
`μN/lcm(c,d)+o(N)` and `μN/c+o(N)`.

Here c=2, a=1, d=3, so g=1 and the share is1/3. If d were even, an odd
component would instead contain no accepted points. Thus arithmetic
cosets are not automatically independent of real-component parity;
compatibility of their clocks is exactly the missing guard. This is a
legitimate connection to the repository's changing-modulus/source-clock
work, not a claim that elliptic dynamics is Collatz dynamics.

## 4. Boundaries retained

- Rank1 and the specified generator are cited inputs, not re-proved by
  this density argument.
- The finite torsion label e and the three cyclic kernel owners are
  not erased by an asymptotic proportion.
- Excluded denominator points occur only at n=0 and remain excluded.
- No least positive height, effective next-return bound, uniform error
  term, or density under arbitrary naive-height enumeration is proved.
- The count concerns rational points on one fixed isogeny target; it is
  unrelated to a proposed upper bound21 for elliptic-curve ranks.

## 5. Exact Haar-arc integral

The proportion μ above can be expressed without choosing an elliptic
logarithm or numerically locating a first positive point. Set

    q_-=(29−5sqrt(65))/392,
    q_+=(29+5sqrt(65))/392,
    D(q)=(q+1)(1+29q−196q²).

The bounded oval has q∈[q_-,q_+] with two W branches. Its invariant
real measure is proportional to `|dq/W|`; the factor from the two
branches cancels between numerator and denominator. Therefore

    μ = [ integral_(1/6)^(q_+) dq/sqrt(D(q)) ]
        /[ integral_(q_-)^(q_+) dq/sqrt(D(q)) ].  (3)

Both improper endpoint integrals converge, and the denominator is
positive. For a nonsingular integration chart put

    Delta=q_+−q_-,
    theta_0=arcsin sqrt((1/6−q_-)/Delta),
    q=q_-+Delta sin²(theta).

Then the common factor1/7 cancels and(3) becomes

    μ = [ integral_(theta_0)^(pi/2)
                 (1+q_-+Delta sin²(theta))^(-1/2) dtheta ]
        /[ integral_0^(pi/2)
                 (1+q_-+Delta sin²(theta))^(-1/2) dtheta ].

This is an exact incomplete/complete elliptic-integral ratio, not an
asserted algebraic number. The two torsion partners have translated,
disjoint positive windows by `FULL_IMAGE.md` §5; hence2μ<=1. The
one-third reconstruction ratio cancels μ completely, so neither a
numerical quadrature nor a conjecture about μ is required for(1).
