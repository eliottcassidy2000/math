# Independent audit: the positive square-root phase repairs the named Laurent cone

**PASS: the exact positive-branch repair and its field boundary are accepted.**
This audits
`05-knowledge/results/continuing1_20260906_laurent_real_power_boundary.md`
and its source, together with the actual cubic definitions and integer-cone
obstruction in the frozen main report
`continuing1_20260906_laurent_cone_separator.md`. The result concerns the
already-safe support(-21,1,12). It is a new observer certificate for this
fixed row, not a general Laurent noncancellation theorem.

No mathematical repair is requested. The independent source is
`04-computation/continuing1_20260906_laurent_amplitude_audit.py`.
It imports no producer code or certificate JSON. It uses only exact integer
and rational arithmetic, including direct bisection of the sextic in the
square-root variable. Its55 always-active gates pass normally and under
optimization with identical LF output bytes.

## Exact source and reduction

Literal enumeration of counts x+y+z=n and−21x+y+12z=0 gives the first
gamma-exponents1,3,5,7 with coefficients110,4620,9240,330. At n=22 the
gamma-exponents0,2,...,14 have coefficients

    22,43890,5969040,135795660,640179540,597500904,77597520,319770.

Thus the first phase equation is exactly
`p(s)=s^3−28s^2+14s−1/3`, and the carried target T=s Q_raw(−s) has constant
term−22. The lower carry has not been omitted or independently normalized.

The verifier reconstructs H=(1+u)^11(s^4−20s^3u^2+28s^2u^4−10su^6+u^8)
coefficient by coefficient, differentiates each coefficient, and sums all
ordered source degrees adding to14. This gives all seven original J_r,
r=0,...,6, by a separate coefficient route. Reducing these polynomials
modulo p verifies exactly that T=R J_0 for the reported rational
R(s)=a+bs+cs^2. It also verifies a,b>0 and c<0. No value of R was obtained
by numerical interpolation at the roots.

Let z_i be the three **positive** square roots and E_1,E_2,E_3 their
elementary symmetric functions. Their selected cubic is

    q_+(z)=z^3−E_1z^2+E_2z−E_3,

with E_3=1/sqrt(3)>0. Multiplying(z+E_1)q_+(z) gives

    z^4−(E_1^2−E_2)z^2−(E_3−E_1E_2)z−E_1E_3.

Consequently the producer's exact identity is correct:

    R(z^2)−(A+Bz+Cz^2)=c(z+E_1)q_+(z),
    A=a+cE_1E_3,
    B=c(E_3−E_1E_2),
    C=b+c(E_1^2−E_2).

This is an identity in independent symbols before imposing any root
relations. The independent sparse-polynomial check verifies every
coefficient. At the three selected z_i it yields
`T(s)=(A+B sqrt(s)+Cs)J_0(s)`.

## Exact signs and the half-power separator

The independent root engine starts from rational z-brackets
(158318,158320)/10^6, (695519,695522)/10^6, and
(5243203,5243205)/10^6. Opposite sextic signs and48 rational bisections
give disjoint positive intervals; their squares isolate all three roots
of the cubic. These brackets are constructed directly in z, rather than
using the producer's square roots of its s-brackets.

Rational interval arithmetic verifies all21 signs J_r(s_i)<0, all three
signs R(s_i)>0, and the following exact millionth enclosures:

    1.100291 <= A <= 1.100292,
    0.974426 <= B <= 0.974427,
    0.029415 <= C <= 0.029416.

In particular C's smaller positive margin is retained. All three
coefficients are real algebraic and strictly positive. This proves the
positive-combination certificate and hence T(s_i)<0 on the selected branch,
relative to the already certified original-window signs.

The Lagrange numerator was also reconstructed directly from the main
functional y=(y_0,y_1,y_2):

    N(s)=(y_2−28y_1+14y_0)+(y_1−28y_0)s+y_0s^2.

This agrees with the producer's polynomial. Evaluating
`sum_i N(s_i) sqrt(s_i) J_r(s_i)/p'(s_i)` through the independent intervals
gives strictly negative bounds for r=0 and6. Both independent enclosures
lie inside the producer's published intervals. Thus the fixed separator
really ceases to be nonnegative at exponent1/2. The integer-only convexity
proof does not control the minimum between integer exponents. Separator
crossing alone would not imply a representation; the formal identity and
positive A,B,C above provide that additional step.

## The field and branch restrictions

Suppose a finite nonnegative real combination of
`z^m J_r(z^2)`, m an integer and r=0,...,6, equals T(z^2) at all six roots
of p(z^2). Subtract its values at z_i and−z_i. The even target disappears.
Every nonzero odd-power summand on the positive side is strictly negative:
its coefficient is nonnegative, z_i^m>0 even for negative m, and J_r(s_i)<0.
Their sum can vanish only if all odd-power coefficients are zero. The
remaining even powers are exactly the nonnegative integer-Laurent window
combination excluded by the main cone theorem. The full-six-root obstruction
therefore holds even with arbitrary real coefficients.

The rational positive-branch obstruction is independently confirmed by a
different irreducibility certificate. Reducing the primitive integer cubic
3p modulo7 gives3s^3−1. Cubes modulo7 are0,1,6, while its required cube is5;
there is no root, so the cubic is irreducible over F_7 and over Q. For any
cubic root s_i, its norm from Q(s_i) is1/3. A square in that field would
have rational-square norm;1/3 is not a rational square, since its3-adic
valuation is−1. Thus sqrt(s_i) has degree two over Q(s_i), and degree six
over Q. The degree-six polynomial p(z^2) is irreducible.

A rational Laurent identity valid at even one positive z_i, after its
nonzero monomial denominator is cleared, must therefore vanish modulo this
irreducible sextic and hold at all six conjugates. The parity argument then
excludes a nonnegative representation. This also applies to the entire
seven-window family, not just J_0. It does not exclude signed rational
identities, which already exist, and does not contradict the real algebraic
positive-branch coefficients A,B,C. Algebraic conjugation changes the
selected positive factor and its coefficient data.

The original Laurent coefficient can be written gamma=i sqrt(s), so the
new coordinate has a concrete source meaning. This observation supplies
no uniform higher-degree amplitude theorem. Selecting the positive square
roots and allowing their real algebraic coefficient field are essential
parts of the proved repair, not removable presentation choices.

## Reproduction and promotion scope

```text
python -B 04-computation/continuing1_20260906_laurent_amplitude_audit.py
python -B -O 04-computation/continuing1_20260906_laurent_amplitude_audit.py
```

The final audited sidecar source has SHA256
`d5b03a6a0034c378e944185f192e3af4d0e4826100d2839c58fd615a663250e6`;
its output hash was
`f6d26209df03fe3d0c81eb54c8eef889304e409a708691a2c7fc9db602bd3075`.
Its original source docstring's phrase denying a "full nonnegative
representation" was flagged only as a wording ambiguity. The producer
clarified it before final pinning to distinguish the selected-branch repair
from rational/all-six identities. The proof and frozen output are unchanged.
No producer file was edited by this audit.

Promotion is supported by a complete independent proof read, literal
constant-term source reconstruction, a separate differentiated-window
coefficient engine, exact quotient equality, independent sextic root
isolation and signs, formal algebra, and the separate mod7/norm degree
argument. The inherited all-integer cone obstruction remains the named
dependency for the two impossibility statements. No all-parameter Laurent
closure or rational positive-coefficient identity is promoted.
