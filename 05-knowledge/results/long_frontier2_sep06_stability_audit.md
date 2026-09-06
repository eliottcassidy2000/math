# Independent audit of regular-root global-modulus barriers

**Status: INDEPENDENT ANALYTIC AND EXACT SOURCE AUDIT — PASS.**
This audits [the proof](long_frontier2_sep06_stability.md),
[source](../../04-computation/long_frontier2_sep06_stability.py), and
[frozen output](long_frontier2_sep06_stability.out). No mathematical
correction was required. Root performed this audit independently of the
author and independently derived the normalization and radical identities.

## Scope and actual normalization

The domain is the inherited finite real list with p1=p2=1 and E>0.
R, Delta and M agree with THM-4455 / three-atom minimizing-sequence
rigidity and THM-4456 / sharp finite-length signed-root stability
asymptotics at their full canonical slugs. The local-modulus note supplies
the local lower limits; neither global right endpoint is advertised as
a proved optimal coefficient.

For any fixed square-normalized q, its sum S is finite. The displayed
alpha is positive for n>=2 because

    n(n+S^2-1)-S^2=(n-1)(n+S^2)>0.

Direct substitution verifies the chosen quadratic root. With
beta=(1-alpha S)/n, the first moment is exactly one and the quadratic
is equivalent to alpha^2+n beta^2=1. Alpha tends to one; beta is O(1/n),
and the dust square mass tends to zero. Positive ordered coordinates and
third/fourth moments therefore converge. This conclusion does not require
the limiting list's first sum to be one: the dust carries its signed defect.
The quotients used for m>=4 have strictly positive limiting denominators.

## Regular families and both barriers

Independently setting p3=x,p4=x^2 and the leading coordinates equal to x
in the definitions gives exactly

    R=4(sqrt2-1)/(3(1+x)),
    (R-K3)/Delta=2(sqrt2-1)/(3(1+sqrt3)(1+x)).

This computation used direct radical simplification, not the producer's
formal two-square remainder path. For integer m>=4, x=1/sqrt(m) decreases
with m, so the distance quotient increases strictly and its family minimum
is at m4. Substituting x=1/2 gives

    kappa4=2(sqrt2-1)(sqrt3-1)/9.

For the moment quotient, M=(1/sqrt3-x)^2. Its denominator after cancellation
contains (1+x)(1/sqrt3-x), whose derivative is negative in this range.
Thus the quotient decreases with m, with diffuse limit
4(sqrt2-1)/(1+sqrt3). The diagonal choice of actual lift parameters is
valid; no uniform two-parameter error estimate is asserted. The excluded
m3 endpoint is singular and cannot supply the local constant by this
cancellation.

The positivity of both global infima is a valid sequential contradiction.
Delta<=2 and M<=4/3+2/sqrt3 are uniform on the actual domain. A vanishing
global ratio would force R to K3 and then violate the inherited strictly
positive local liminf bound. This proves existence of positive global
coefficients without claiming either upper barrier as a lower bound.

## Finite hostile and multiplicity loss

For the rational family, literal sums give first numerator2n+1 and square
numerator(2n+1)^2. All entries have magnitude less than one, so p4<p2=1
and the energy is positive for every n>=1. At n1 the leading positive tie
causes no change to the sorted top-three sum. The degree-four ordinary
product reconstruction is an independent route to the duplication ratio.
The exact n1000 interval lies below1/10 and the finite length is1005.
The result is a counterexample to that global distance coefficient, not
to the different finite-length coefficient inequality.

The x4/7, formal multiplicity49/16 example is correctly marked nonactual.
Equality in p3^2<=p2 p4 forces r_i^2=p3 r_i at every coordinate, hence every
nonzero coordinate must equal4/7. Integer multiplicity then contradicts
p2=1. The lost coordinate is therefore concrete, and the example cannot
be used as an actual counterexample to the candidate kappa4 bound.

## Replay and frozen identities

Root independently replayed the source in normal and optimized Python,
captured raw bytes, and compared both with the frozen transcript. All40
always-active checks pass and every byte agrees. A separate direct-radical
calculation also verified the general quadratic root, regular distance
identity, four-atom constant and diffuse moment constant.

    source SHA256 bf0734ba632810b307b0a3dde8b788fcf6dafb4f226ed0c27308d86b7d086f8e
    output SHA256 cd9662c07adeca90b24086fd029186f1db13987cf240e20480d651d5a10c588c

Reproduction:

    python3 -B 04-computation/long_frontier2_sep06_stability.py
    python3 -B -O 04-computation/long_frontier2_sep06_stability.py

The global matching lower bounds, exact finite-N optima, and the proposed
N>=6 coefficient bound remain OPEN. No external priority or proof-assistant
verification is claimed.
