# A positive square-root-phase repair, with an exact field and branch boundary

**Status: PROVED for the named actual row, with independently audited exact
rational certificates.** The [independent referee](continuing1_20260906_laurent_amplitude_audit.md)
accepts the complete proof, reconstructing the literal coefficients and
isolating the square-root variable independently; its 55 always-active
gates pass in normal and optimized Python. This is a separate sidecar to the
[actual integer-Laurent cone obstruction](continuing1_20260906_laurent_cone_separator.md).
The main producer is unchanged. The result concerns its fixed genuine
support (-21,1,12), not a new all-parameter Laurent separation theorem.

The fixed separator fails already at the real exponent 1/2. More positively,
the actual response admits a strictly positive combination of J_0,
sqrt(s)J_0, and sJ_0 on the three positive first roots, with explicit real
algebraic coefficients. This is impossible with rational coefficients,
and impossible for a nonnegative half-Laurent identity on all six roots
of the squared-variable equation. The field and selected branch are
essential parts of the repair.

## 1. A certified crossing of the fixed separator

Keep the main note's actual quotient p(s)=s^3-28s^2+14s-1/3,
ordered positive roots r_1,r_2,r_3, window responses J_r, target T,
and rational functional ell. For an arbitrary real exponent theta,
define its action by the actual positive-root evaluation:

    ell(s^theta J_r)=sum_(i=1)^3 N(r_i) r_i^theta J_r(r_i)/p'(r_i).

For theta=1/2 there is no logarithm or floating-point evaluation in the
certificate. The main note's rational root brackets and integer square
roots give enclosures

| Positive square root | Lower bound | Upper bound |
|---|---:|---:|
| sqrt(r_1) | 39579743193/250000000000 | 15831897293/100000000000 |
| sqrt(r_2) | 347760017139/500000000000 | 69552094839/100000000000 |
| sqrt(r_3) | 1048640776277/200000000000 | 2621601940693/500000000000 |

Each lower bound squared is below the corresponding original root bracket,
and each upper bound squared is above it. Rational interval arithmetic then
proves

    -19914989775031856088736
       <= ell(sqrt(s)J_0) <= -6371590793283941536065 < 0,

    -88591556071110678488277235694
       <= ell(sqrt(s)J_6) <= -31956183471936284505432288886 < 0.       (1)

Thus nonnegativity of the fixed dual does not extend from integer to real
exponents. This matches the main proof: strict convexity forces the
integer minimum to occur at zero or one, but does not prevent a negative
dip between them. Since sqrt(s)>0 and J_r<0 at every first root,
sqrt(s)J_0 is a new sign-certified generator that actually crosses the wall.
Wall crossing alone is not a representation; the next section provides one.

## 2. Positive algebraic coefficients on the selected branch

The main certificate gives T==R(s)J_0 modulo p, where R(s)=a+bs+cs^2 and

    a=142115562175391338833022911/115962939903341750549938130,
    b=137850584919079100136401223/115962939903341750549938130,
    c=-825111792094668079242879/23192587980668350109987626.

Let z_i=sqrt(r_i)>0, and let E_1,E_2,E_3 be their elementary symmetric
functions. In particular E_3=1/sqrt(3). Set

    A=a+cE_1E_3,
    B=c(E_3-E_1E_2),
    C=b+c(E_1^2-E_2).                                           (2)

The exact identity in independent symbols is

    R(z^2)-(A+Bz+Cz^2)
      =c(z+E_1)(z^3-E_1z^2+E_2z-E_3).                          (3)

The right side vanishes at all three positive z_i. Applying rational
interval arithmetic to the preceding square-root bounds gives

    1.100291 <= A <= 1.100292,
    0.974426 <= B <= 0.974428,
    0.029415 <= C <= 0.029416.                                  (4)

Every endpoint in (4) is an exact rational with denominator one million.
Hence all three constants are strictly positive, and at each actual first
root

    T(s)=(A+B sqrt(s)+Cs)J_0(s)<0.                              (5)

This is a real algebraic, positive-branch observer certificate. It recovers
the known actual negative response through a mechanism which the main
integer-Laurent cone cannot express. The positive square-root amplitude
is already present when the original coefficient is written
gamma=i sqrt(s); using only t=gamma^2=-s had removed that square-root
branch coordinate. This observation does not change the
original coefficient law or grant an all-parameter amplitude theorem.

## 3. Why rational coefficients and all six roots remain obstructed

Let z^2=s. Suppose a finite nonnegative combination of
z^m J_r(z^2), m in Z, represented T(z^2) modulo p(z^2), on all six roots.
The odd part of this identity must vanish. At any positive z_i, every
nonzero odd-power summand has strictly negative value, because z_i>0
and J_r(r_i)<0. Thus all odd-power coefficients must be zero. The remaining
even terms give exactly an integer-Laurent representation excluded by the
main theorem. This argument allows arbitrary real nonnegative coefficients;
the all-six-root identity is impossible even over the reals.

There is also an obstruction to rational-coefficient identities restricted
only to the positive branch. The primitive integer cubic

    3s^3-84s^2+42s-1

has no rational root: the only candidates are +/-1 and +/-1/3, and direct
substitution excludes all four. Thus p is irreducible. In K=Q(r_i),
Norm_(K/Q)(r_i)=1/3. If r_i were a square in K, its norm would be a
rational square; 1/3 is not. Therefore adjoining sqrt(r_i) has degree
two over K, and p(z^2) is irreducible of degree six over Q.

A rational Laurent identity valid at even one positive z_i must therefore,
after clearing its nonzero monomial denominator, hold modulo p(z^2).
The preceding parity argument excludes it. Consequently the algebraic
constants and the positive-root factor in (2)--(3) are not cosmetic.
Equation (5) is not a rational coefficient identity modulo p(z^2).

## 4. Exact reproduction and scope

[Source](../../04-computation/continuing1_20260906_laurent_real_power_boundary.py)
and [output](continuing1_20260906_laurent_real_power_boundary.out)
use the frozen primary JSON by a pinned SHA256. There are 35 always-active
gates: original polynomial brackets, exact rational square-root bounds,
negative half-power responses, positive coefficients, a formal polynomial
identity in E_1,E_2,E_3,z, and the rational-root/norm controls. No numerical
root or logarithm is used.

```
python -B 04-computation/continuing1_20260906_laurent_real_power_boundary.py
python -B -O 04-computation/continuing1_20260906_laurent_real_power_boundary.py
```

Normal and optimized transcripts have identical LF bytes. Frozen SHA256:

    source d5b03a6a0034c378e944185f192e3af4d0e4826100d2839c58fd615a663250e6
    output f6d26209df03fe3d0c81eb54c8eef889304e409a708691a2c7fc9db602bd3075
    input  4b1ee5770b484e4164e692fbf2934f4099800b0b85d379ce01d8afef71040cc0

The exact cubic sidecar is complete. Any extrapolation to larger first-row
degree, rational coefficients, or an unselected square-root branch is a
different claim and is not used here.
