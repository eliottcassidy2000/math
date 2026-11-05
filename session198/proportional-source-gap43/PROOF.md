# A remaining target already requires leaving singleton-word transport

Session 198. **PROVED** M-word rigidity and source-ray obstruction at 43;
**FINITE-EXACT** all-order-bounded C1 control through 43. This is a precise
failure of a proposed coverage architecture, not a negative spectrum theorem.
Every source pair below consists of actual nontrivial start-one sources.

## 1. The M form has the same three-test rigidity as N

For an actual nontrivial source v=(D,E,U,d), let J(v)=(D-E+U,D-E,D+U,D). Put

    delta_k=M(v,J^k w)-M(J^k v,w),

where M is the complete physical C2 form. Then delta_1=delta_2=delta_3=0 if
and only if v,w are proportional as full four-coordinate numerical profiles.
This is also equivalent to delta_k=0 for every k>=0.

Here is an exact rank proof. As linear rows acting on w, the first three defects
are

    [-E+3U-2d,  D+d,      -3D-d,       2D-E+U],
    [-E+7U-2d,  D-2U+d,   -7D+2E-3d,   2D-E+3U],
    [-E+11U-4d, D-2U+d,   -11D+2E-5d,  4D-E+5U].

Their signed maximal minors equal -4Q(v)(D,E,U,d), where

    Q=2D^2-DE-DU+3Dd+EU-Ed+U^2-Ud+d^2.

The physical nonnegative coordinates

    a=D-E-d, f=U-D+E-d, e=E, p=d-E

satisfy d=e+p>=1, and substitution gives

    Q=2a^2+10ae+af+7ap+12e^2+2ef+17ep+f^2+2fp+6p^2
      >=6(e+p)^2>0.

Thus the row rank is three. The source vector v is in its nullspace by
symmetry of M, so that nullspace is exactly the ray through v. Proportionality
then gives equality for every k by bilinearity, linearity of J, and symmetry.
This proves every implication.

The separate N_WORD_RIGIDITY.md proves the identical three-test equivalence for
N by two simple difference identities. Hence, in either physical C2 family,
unrestricted shifting of a singleton word across a fixed source pair is
possible exactly when the full source profiles are proportional.

## 2. Prime counts in a proportional pair must be diagonal

For an actual nontrivial source, define its integral content as gcd(D,E,U,d). Since each
C2 expression has integral bilinear coefficients, the contents of either input
divide its output. If the output is a positive prime p, each content is one or
p. The latter is impossible: it would force U>=p, while the other nontrivial
source has U>=2 and the reviewed product bound gives H>=3U_1U_2+1>p for N;
M has an even stronger lower bound. Thus both source profiles are primitive.

Two positive primitive integer vectors that are rationally proportional are
equal. Consequently every prime N or M output from proportional actual source
profiles must be a diagonal value N(v,v) or M(v,v). This conclusion does not
assert that primitive vectors satisfying the inequalities are physical.

## 3. The exact obstruction at 43

No proportional actual source pair realizes the prime 43 in either C2 family.

For diagonal N, the source-product bound is

    N(v,v)>=3U^2+1.

Thus U>=4 gives at least 49. At U=2 and U=3, the only physical source responses
are respectively t=(1,0,2,1) and q=(3,1,3,1). They give

    N(t,t)=13,       N(q,q)=41.

There is no other nontrivial source with U<4. For diagonal M,

    M(v,v)>=6U^2+2U+1.

Every U>=3 gives at least 61, and the unique U=2 source gives M(t,t)=29.
The regular-five exception has count 15. This proves the exclusion at 43.

But 43 is genuinely attained by the nonproportional physical pair

    v=(1,0,2,1),        w=(3,0,6,1),
    N(v,w)=43.

The left source is transitive on two vertices. Build the right source as a
triangle with a universal sink appended, retaining a triangle root. The N graft
has seven vertices and exactly two Hamilton cycles. The verifier reconstructs
and counts this literal graph.

For this particular pair the first three N singleton-shift defects are

    (Delta_1,Delta_2,Delta_3)=(-8,-16,-28).

Thus its needed nonproportional data cannot be replaced by a scalar-preserving
singleton-word shift. By the theorem, every N or M pair realizing 43 must fail
at least one of the first three shift tests.

## 4. Why 43 is relevant to the remaining-target problem

The complete strong C1 class through 43 omits exactly

    7,13,21,41,43.

This is an all-graph-order finite numerical statement, not an order census.
The complete physical source inventory through U=86 is sufficient because every
strong start-one source with D<=43 has U<=2(D-E)<=86. The reviewed C1 normal
form identifies its H image with the strong C1 spectrum; the actual generator
is used, not the inequality cone. Among those omitted values, 13 and 41 have
the diagonal N witnesses just shown. After separately setting aside 7 and 21,
43 is the first target where restricting to proportional C2 sources already
fails. No unrestricted impossibility of 7 or 21 is proved in this packet.

The conclusion is constructive as well as restrictive: an explicit physical
source pair reaches 43, but the simplest whole-singleton-word reassociation
architecture cannot account for that witness or any other witness of 43.
Other target-dependent operations remain available and exact coverage remains
open.

## 5. Verification

Run `python3 verify.py` and `python3 -O verify.py`. Normal and optimized receipts
must agree. The standard-library checker verifies the M defect matrices on all
basis pairs, the polynomial maximal-minor identity and positivity expansion on
finite unisolvent grids, the three-test criterion on a complete bounded physical
source bank, the complete C1 control through43 using its proved source cap, the
small-source diagonal restrictions, and the literal H43,C2 witness. The grids
are polynomial-identity checks with explicit degree bounds, not empirical
extrapolation of a spectrum. Mathematical dependencies are frozen below.
