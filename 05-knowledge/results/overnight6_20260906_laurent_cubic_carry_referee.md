# Independent mathematical referee: endpoint-21 cubic carry family

**Verdict: PASS, analytic all-parameter proof and exact certificate.**
Reviewed `overnight6_20260906_laurent_cubic_carry.md/.py/.out` and the
root's separate raw-fiber/resultant interpolation audit. This review uses
an additional standard-library polynomial path, importing neither script.
No numerical roots or enlarged census are used. The only prose correction
found is that Appendix A's displayed coefficients are positive rationals,
not all integers; the parent has acknowledged that correction. It does not
alter a coefficient or mathematical claim.

## 1. Source, completeness, factorials, and detection scope

The original charge equation is exactly `g(2y+3z)=21m`, with
`x+y+z=m`. For integral `g>=11`, `gcd(g,21)=1` gives `g|m`.
At m=g, solving `2y+3z=21` gives `z=1+2j`, `y=9-3j`,
`x=g-10+j`, with exactly `0<=j<=3`. At m=2g the corresponding solution
is `z=2j`, `y=21-3j`, `x=2g-21+j`, exactly `0<=j<=7`.
All entries needed for the factorial denominators are nonnegative, and
the first x entries are positive on the stated parameter domain.

The phase ratio is `tau=alpha*gamma^2/beta^3`, because the primitive
fiber step is `(1,-3,2)`. The first anchor is
`X=alpha^(g-10)*beta^9*gamma`; the second anchor is **X^2/tau**.
This independently verifies the lower carry, including its direction.
Dividing the first raw coefficients by `(g)_7/9!` gives, in increasing
phase order,

```text
(g-7)(g-8)(g-9), 84(g-7)(g-8), 504(g-7), 72.
```

Dividing the second raw coefficient j by `(2g)_14` gives exactly
`(2g-14)_(7-j)/((21-3j)!(2j)!)`. Both contents are strictly positive.
Thus the two literal identities in the producer are correctly normalized;
they are identities of complete coefficient monomials, not coefficient
specializations. At nonzero coefficients, X and tau never vanish.

No positive return mass exists before g or strictly between g and 2g.
The strict negative anchored second value at each first root therefore
proves the stated **first detection exactly g or 2g**, for every nonzero
complex coefficient triple. Positive coefficients attain g; setting
`beta=gamma=1` and alpha equal to a first root attains 2g. The normalized
negative real sign does not assert that the unnormalized complex moment
is itself negative. The gcd-dropped g=12 control correctly has first
support return four and so prevents an overbroad detection interpretation.

## 2. Real spectrum, all characteristic signs, and degree bounds

The monic cubic has coefficients

```text
p2=7(g-7), p1=7(g-7)(g-8)/6, p0=(g-7)(g-8)(g-9)/72.
```

An independent five-term cubic discriminant formula reproduces the
producer's positive factorization. Positive coefficients exclude every
nonnegative root; positive discriminant gives three distinct real roots.
Consequently evaluation diagonalizes the quotient multiplication operator
over the reals. All three positive characteristic coefficients then
exclude every nonnegative eigenvalue. The middle coefficient is essential
in dimension three: `(-5,1,1)` has the misleading trace/product signs.
Real spectrum is equally essential; the producer's positive-coefficient
cubic `z^3+z^2+z+2` has its sole real root in `(-2,-1)` and hence a
nonreal conjugate pair with positive real parts.

The interpolation degree justification is complete. Assign g and tau
weight one. Each monic cubic coefficient has degree at most `3-j`, so
replacing `tau^3` by `-p2*tau^2-p1*tau-p0` never increases weight.
Every noninverse term of `Qbar/tau` has weight at most six. Its sole
inverse contribution has the exact cancellation

```text
q0/p0=1152(g-10)(2g-15)(2g-17)(2g-19)/21!,
q0/tau=-(q0/p0)*(tau^2+p2*tau+p1) mod p.
```

It too is polynomial of weight at most six; no rational pole is left.
Thus `deg_g R_j<=6-j`. Newton recursion gives
`deg_g Tr(C^ell)<=ell` for every ell>=0. Expanding `R(C)^a` then gives
`deg_g Tr(R(C)^a)<=6a`. Newton identities for the characteristic
coefficients give `deg B_i<=6i`. This proves the degree bounds 6,12,18
without using finite samples. Nineteen distinct rational parameter values
suffice for exact identities. The root's samples include nonprimitive g;
that is valid because the mass-g and mass-2g fibers and polynomial identities
remain correct there. Their first-return interpretation is not used in
interpolation. No unproved extrapolation is present.

## 3. Independent full symbolic certificate and the SOS boundary

The companion `overnight6_20260906_laurent_cubic_carry_referee.py` implements
exact arithmetic in `Q[g]` using only `fractions.Fraction`. It derives
the two rows directly from falling factorials, performs exact content
division, reduces the carried response modulo the monic cubic, and builds
its three-column multiplication matrix. It obtains B2 from its three
principal 2x2 minors and B3 from the six determinant terms, separately
from the producer's trace formula and the root's resultant interpolation.
All **39** rational coefficients of the three shifted polynomials agree
exactly with the published arrays and are strictly positive. This is a
direct full symbolic identity check, not another interpolation.

Additional controls use literal two-loop multiplicity enumeration at five
primitive parameters to recover all ten complete rows, multinomial
weights, and every phase-step exponent. The shifted discriminant is
derived independently, and both the gcd and real-spectrum hostiles are
retained. All **130** optimization-live gates pass; normal and optimized
LF outputs are identical.

The ordinary-core exclusion is correctly mapped to the source. Write
`f(u)=u^-21 h(u^g)` with `h(s)=alpha+beta*s^2+gamma*s^3`. After a scalar
and variable gauge, `h` becomes `H(s)=1+b*s^2+s^3`, with
`b^3=1/tau`. At the displayed g=11 root in `(-28,-27)`, its discriminant
is `-4/tau-27<0`; H and its powers have nonreal roots. Hence the
real-rooted-core hypothesis of the current
`THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md` does
not cover this cancellation phase. This is a valid nonsubsumption
example, not a counterexample to that theorem. The first-row polynomial's
own real-rootedness is a different object and cannot repair the core
hypothesis.

The proof may be accepted at its declared fixed endpoint-21, unbounded-g
scope. It proves neither a general endpoint-21 strip nor the next quartic
family, and makes no external priority claim. The lower carry, complete
characteristic vector, and real-spectrum condition are all indispensable
parts of the preserved negative-root-response predicate.

```text
python -B 04-computation/overnight6_20260906_laurent_cubic_carry_referee.py
python -B -O 04-computation/overnight6_20260906_laurent_cubic_carry_referee.py
published coefficient record SHA-256
88d4c05a95fa48066151b72063c5f9fffb0cbecb18f79ca91c975da63f935315
referee source SHA-256
2f8fd8cb6eb0540ba3ee31551a0e33eac18fdc7845495f6c8d66fe160876a9d1
referee output SHA-256
13c53421a1db90f431dc3eabdcab9f9f0e720cb1bdbb52472987b854313bd061
```

**Filing:** root integrated this audited report after `f5f0f7f75`;
portable reproduction paths are shown above. The exact producer and
transcript bytes remain pinned by the sixth manifest.
