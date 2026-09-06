# Independent endpoint-27 and complete midpoint-transport audit

**PASS: analytic proof plus independent exact reconstruction.** The
[endpoint-27 family](overnight7_20260906_laurent_quartic_carry.md) and the
all-parameter identities and degree bounds in
[the midpoint report](overnight7_20260906_laurent_midpoint_transport.md)
are accepted. The 90 grouped signs are **FINITE-EXACT** only. Universal
group negativity and general trinomial two-rung separation remain **OPEN**.

## Endpoint-27 source and degree audit

The charge equation is `g(2y+3z)=27m`. For `gcd(g,27)=1` all support
returns have mass divisible by g. At g and 2g its complete nonnegative
fibres are respectively

```text
(g-13+j,12-3j,1+2j), j=0,...,4,
(2g-27+j,27-3j,2j), j=0,...,9.
```

All counts are valid at g>=14. The first and second anchor ratio retains
the factor `tau^-1`; deleting it changes the rootwise sign. Nonzero complex
coefficient monomials may be divided out, but a sign assertion applies to
the normalized scalar response, not to an arbitrary complex raw moment.

The first monic polynomial p has coefficients

```text
p_j=(g-9)_(4-j) * 9! / [(12-3j)!(1+2j)!], j=0,...,4.
```

Its four roots are simple and strictly negative by the complete factorial
row theorem [THM-4436 / complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
with the checked positive-count tuple `(A,B,h,r,z,x)=(2,3,4,0,1,g-13)`.
A positive quartic discriminant alone would not establish this reality.

Divide the doubled raw row by `K=(2g)_18` to obtain Qbar. The only apparent
parameter denominator in reducing Qbar/tau modulo p comes from its constant
term times tau^-1. It cancels exactly:

```text
Qbar_0/p_0 = 42240(g-13)(2g-19)(2g-21)(2g-23)(2g-25)/27!.
```

Assign weight one to g and tau. Every noninverse term of Qbar/tau has
weight at most eight; the inverse substitution has the same bound after
this cancellation. Reduction by the monic weighted polynomial p preserves
the bound, so the reduced coefficient R_j has degree at most `8-j`.
Multiplication by R has entry-degree bound `8+column-row`. Each size-k
principal minor therefore has degree at most 8k, because the row and column
indices cancel term by term. The four characteristic coefficient bounds
are consequently 8,16,24,32. This is an all-parameter degree proof.

The independent verifier imports no producer mathematics. At every integer
g=14,...,46 it solves the raw charge equation afresh, forms exact factorial
coefficients, and computes the resultant of `p(tau)` and
`tau*X-Qbar(tau)`. Its monic normalization has exactly the four response
values Qbar(rho)/rho as roots. All coefficients match the producer's shifted
polynomials; independent rational interpolation reconstructs the complete
identities. Thirty-three points exceed every proved degree bound.
Nonprimitive parameters are valid algebraic checks, not first-return claims.

All 84 entries of the shifted numerator arrays are strictly positive, with
positive denominators. Thus the characteristic polynomial is positive at
every nonnegative real X, including zero. Every response is real by the
first-root theorem, so all four are strictly negative. This proves first
detection at g or 2g for every admissible parameter and coefficient triple.
Both alternatives occur, since the first polynomial is nonconstant with
four nonzero roots and is nonzero at generic phase. Width is 3g; no sharp
width claim follows. The gcd-dropped g=15 hostile has first support mass5.

## Complete all-h midpoint identities

The audit checks the proof before evaluating any sign. The alpha parity
split is the even-index part of the binomial convolution:
`A_double=O^2+t^-1 E^2`. The beta split cuts a lattice path at the weighted
level `y+3x`. A hit gives B^2. A vertical step may skip that level from
either of its two preceding levels; the two resulting products are both
tCD after the stated index shift. These cases partition the path set,
proving `B_double=B^2+2tCD` without a sign argument. Complete supports
include the negative beta indices. Truncating them would destroy the identity.

Expanding the two factors before the Hadamard product yields exactly the
three displayed correction groups. It preserves the full actual doubled
multinomial row, rather than merely its virtual subfamily. The Euler
identities follow coefficient by coefficient from binomial ratios on the
full supports, whose denominator bounds are strictly positive in the stated
domain. The general cancellation formula for Qbar_0/p_0 follows by splitting
the doubled falling factorial into its even and odd factors. The same
weighted argument gives reduced degree `2h-j` and characteristic bound
`2hk` at arbitrary h. It supplies no positivity of those polynomials.

An independent standard-library convolution implementation reconstructs
every coefficient of these Laurent identities at the 12 named `(h,g)`
controls. For the 30 first roots it tests each of the three groups using
Hermite trace forms, a different sign certificate from the producer's
isolating-interval evaluation. For a response R the form has entries
`Tr(R T^(i+j))`, where T is multiplication by tau in the first quotient.
With simple real roots it is Vandermonde-congruent to the diagonal of the
rootwise R values. All signed leading principal minors are positive, hence
each form is negative definite: exactly 90 strict rootwise signs.

The audit also retains `(h,g)=(1,5),rho=-2`, where
`O*B=0` but `O*C=15`, `O*D=10`, `E*B=-16` (stars denote Hadamard
products evaluated at rho). Replacing a factor does not preserve the
coefficient-zero hypothesis of signed duplication. Positive ordinary
coefficients likewise do not determine evaluation signs at negative rho.
These are boundaries on the proposed bridge, not exceptions to the proved
path identities or the named finite signs.

## Reproduction

```text
python -B 04-computation/overnight7_20260906_laurent_independent_audit.py
python -B -O 04-computation/overnight7_20260906_laurent_independent_audit.py
```

All **553 optimization-live gates PASS**; normal and optimized outputs
agree byte for byte. The producer's quartic and transport checks also pass
normally and optimized. LF hashes of the independent audit are

```text
source 4eac6948a3c850982d322f79f94aa296716fbdc80cf299ea9b1d88b0be5991c1
output db64cee82d2aee965840752ba9a07045bc22dabe833b318edf95ccc97787cf3f
quartic record 7f304a55990b177e13d063706bf8f2df7697ac9ed86f317ba3ab9b48a7af650c
```

Incoming independently audited endpoint-21 Hermite certificates provide
another lower-degree instance of the same carried sign question; this audit
does not count their concurrent proof as a new family. The exact common
object for the next step is the original coefficient-zero hyperplane with
all three coupled correction responses retained.
