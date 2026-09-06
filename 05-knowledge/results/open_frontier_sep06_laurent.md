# Fixing the genuine linear anchor does not restore full signed transport

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The linearly anchored full response is positive for every real
R>=384 in the family below, and an algebraic ratio in (256,384) gives exact
model double cancellation. Actual binomial path transport and
general trinomial two-rung separation remain **OPEN**.

## Inheritance and target

The incoming [coupled-window theorem](creative_20260906_laurent_bridge.md)
and its [independent audit](creative_20260906_laurent_audit.md) were read first.
They reduce all coupled derivative certificates at one vanished coefficient
to `min(k,n-k)` generators. Their conclusion is sufficient negativity for an
explicit weighted opposite-pair response; the actual grouped response still
needs an equality to such a kernel.

The closest joint mechanism and canonical hostile are the audited
[combined Euler/pencil model](combined_pencil_empty_core_morning_sep06.md).
It preserves both contiguous Euler identities, all-weight beta compatibility,
the common binomial carrier and the original first zero. Its free-scale
counterexample leaves the genuine beta linear coefficient unfixed. The
corrected near miss is its normalized ratio-64 control: restoring the first
coefficients `(13,12,11)` restores a negative full response, while its skip
and crossing classes stay positive.

This probe fixes those three coefficients throughout. The underused sidecar
is the beta root **shape after its sum has been fixed**. The six live concepts
are: the full factorial coefficient law; its linear anchor; the relative root
shape; simultaneous Euler/pencil admissibility; the original common zero;
and membership in the new coupled derivative cone. Shape changes preserve
the anchor and the first four model constraints, but need not preserve the
full coefficient law or cone membership.

Targeted repository searches for the linear anchor and its synonyms found
only the predecessor's open next test. No external priority is claimed.

## 1. A rational anchored witness

Fix `q=5,h=4,x=1,g=14,k=9`. For a positive ratio R define

```text
S_R=1+R+R^2+R^3+R^4,
a_i=13 R^i/S_R,                         0<=i<=4,
B(v)=product_i(v-a_i)=sum b_i v^i,
c_i=3(i+1)b_(i+1)/(7+2i),
d_i=(2+3i)c_i/(6+2i),                   0<=i<=4,
C(v)=sum c_i v^i,  D(v)=sum d_i v^i.
```

For **R=512**, the B-root tuple is exactly

```text
(13,6656,3407872,1744830464,893353197568)/68853957121.
```

B,C,D are monic. Reverse with alternating signs to obtain

```text
G_B(t)=product_i(1+a_i t),
G_C(t)=sum_(j=0)^4 (-1)^j c_(4-j)t^j,
G_D(t)=sum_(j=0)^4 (-1)^j d_(4-j)t^j.
```

Every coefficient is strictly positive, their constant coefficients are one,
and the **genuine linear anchors hold exactly**:

```text
([t]G_B,[t]G_C,[t]G_D)=(13,12,11).
```

Exact rational evaluations in the frozen output give
`C(a_i)/B'(a_i)>0` and `D(a_i)/B'(a_i)>0` for all five i. Thus both C and D
strictly interlace B. For every w>0, `wB^2-2vCD` has ten simple positive
roots: its signs alternate through 0, the five B roots, the four intervening
C roots, and positive infinity. At w=0 the degree is nine and the roots are
0 together with those of C,D. Reversal therefore proves all-weight
nonnegative-coefficient real-rootedness of `wG_B^2+2tG_CG_D` for every w>=0.
No finite sample of weights is used for this universal assertion.

## 2. The same original zero has a positive full response

Keep all lower Laurent exponents in the definitions

```text
beta=t^-1 G_B,  C_raw=t^-1 G_C,  D_raw=t^-1 G_D,
O=sum_j binom(14,2j+1)t^j,     E=sum_j binom(14,2j)t^j,
P=O star beta,
V=O^2 star beta^2,
W=(O^2+t^-1 E^2) star beta^2,
R_skip=(O^2+t^-1 E^2) star(2t C_raw D_raw),
Q_model=W+R_skip.
```

The doubled row has a genuine exponent -1 term. The verifier performs raw
Laurent convolution and coefficientwise products before clearing this carry.
At t=-s, put

```text
K_B=s^5 B(u^2/s), K_C=s^4 C(u^2/s), K_D=s^4 D(u^2/s),
H_i=(1+u)^14 K_i.
```

These satisfy the same two Euler identities as the actual contiguous rows:

```text
K_B'=(2/3)u(7K_C+uK_C'),
6K_D+uK_D'=2K_C+(3/2)uK_C'.
```

The exact extraction identities retain the common binomial factor and zero:

```text
P(-s)=(-s)^-1 [u^9]H_B,
W(-s)=s^-2 [u^18]H_B^2,
R_skip(-s)=-2s^-1 [u^16]H_C H_D.
```

There is a unique smallest positive root s_* of P(-s) in

```text
1518379/1000000 < s_* < 1518380/1000000.
```

The verifier prints the exact integer polynomial defining this root. Endpoint
signs and a strictly negative derivative on the whole interval from zero to
the upper endpoint prove both existence and the smallest-root claim.
Rational Horner interval bounds certify, on this entire bracket,

```text
sW(-s) in [-163465.760909,-163465.386864],
sR_skip(-s) in [237855.873369,237865.815274],
sQ_model(-s) in [74390.359366,74400.181498].
```

The terminating decimals here are exact rational outward bounds. In
particular `P(-s_*)=0` and `Q_model(-s_*)>0`, despite every retained model
predicate and all three genuine linear anchors. The virtual response V and
the completed hit W stay strictly negative. The crossing class K stays
positive. This is a full-response hostile, not merely a positive individual
summand.

## 3. An unbounded anchored family and exact double cancellation

The finite witness extends by an exact polynomial certificate. For clarity,
let B_base have roots `(1,R,R^2,R^3,R^4)` before the positive normalization
`lambda=13/S_R`. Put

```text
t=1/R,    Sbar(t)=1+t+t^2+t^3+t^4,
s_base=u t^3,   sigma=s_base/lambda=Sbar(t)u/(13t).
```

Here t is the reciprocal **shape parameter**, not the Laurent variable in
Sections 1--2. Using v as the ordinary carrier's coefficient variable, write

```text
p_base(s)=[v^9]H_(B,base)(s,v)/s,
f(t,u)=t^4 p_base(u t^3)=sum_(j=0)^4 f_j(t)u^j,
f_0=-14(1+t+t^2+t^3+t^4),
f_1=364(1+t+2t^2+2t^3+2t^4+t^5+t^6),
f_2=-2002(t+t^2+2t^3+2t^4+2t^5+t^6+t^7),
f_3=3432(t^3+t^4+t^5+t^6+t^7),
f_4=-2002t^6.
```

In particular `f(0,u)=-14+364u`, selecting the finite phase limit `u=1/26`.
For every `0<=t<=1/256`, the verifier proves the following rational polynomial
interval inequalities. Quotients by t extend polynomially at t=0:

```text
-3 < f(t,1/26)/t < -29/10,
2/3 < f(t,1/26+t/100)/t < 3/4,
363 < partial_u f(t,u) < 366
      for 0<=u<=1/26+1/25600.
```

Thus for every `0<t<=1/256` there is a unique smallest positive root

```text
u(t)=1/26+t z(t),     0<z(t)<1/100.
```

The nonzero derivative proves this branch is continuous (indeed real
analytic). The original first-row phase is exactly
`sigma(t)=Sbar(t)u(t)/(13t)`, throughout, rather than a root selected from
one of the contiguous replacements.

Admissibility also holds on this entire parameter interval. For each
`i=0,...,4` and F=C_base,D_base, normalize
`(-1)^i F(R^i)` by its highest power of R and substitute R=1/t. This gives
ten explicit univariate polynomials. Every one is strictly greater than
`1/14` on `[0,1/256]`, as certified by rational Horner intervals. The
positive sign is exactly the sign of `F(R^i)/B_base'(R^i)`, because the
ordered B roots give `sign B_base'(R^i)=(-1)^i`. Both strict interlacings,
the three linear anchors, the two Euler relations and all-weight pencils
therefore persist for **every real R>=256**.

For a compact full-response certificate, let
`w_base(s)=s W_base(-s)` and `r_base(s)=s R_skip,base(-s)`. Define

```text
T(t,u)=13t^6 w_base(u t^3)+t^2 Sbar(t) r_base(u t^3).
```

All negative powers cancel before evaluation: T is a polynomial of degree
18 in t and 9 in u, with 122 nonzero rational coefficients. Raw carry-aware
homogeneity gives the exact normalized response identity

```text
T(t,u)=t^2 Sbar(t) sigma Q_model,R(-sigma),
sigma=Sbar(t)u/(13t).                                  (12)
```

Its boundary polynomial and selected value are

```text
T(0,u)=-1512u^2(314847u^2-92391u+3025)/143,
T(0,1/26)=47439/48334 > 0.
```

Using the actual root tube before interval evaluation removes the otherwise
large cancellation. Direct polynomial substitution `u=1/26+t z` gives

```text
1/40 < T(t,1/26+t z) < 33/32
for 0<=t<=1/384 and 0<=z<=1/100.                        (13)
```

Equations (12)--(13) prove `Q_model,R(-sigma_R)>0` for **every real R>=384**.
The threshold is a sufficient certificate boundary, not a sharp transition
claim. The proof also gives

```text
sigma_R/R -> 1/338,
Q_model,R(-sigma_R)/R -> 47439/143.
```

The full response at R=256 is strictly negative by Section 4's independent
control, while it is positive at R=384. The same admissible first-root branch
is continuous on the entire closed interval. Therefore some
`R_0 in (256,384)` satisfies exactly

```text
P_(R0)(-sigma_(R0))=Q_model,(R0)(-sigma_(R0))=0.          (14)
```

R_0 can be specified algebraically, without any numerical zero inference.
Clear the rational denominators in R from P and `sigma Q_model`. Their
resultant in sigma is a polynomial in R that is not identically zero: the
verifier finds gcd degree zero at R=256, where both polynomial degrees remain
full. Every solution of (14) therefore has algebraic R_0, and then algebraic
sigma_(R0). One exact choice is the least such ratio in (256,384); there are
only finitely many because this resultant is nonzero. All three linear
anchors remain fixed at this algebraic model. No independent beta scale is
reintroduced.

## 4. Scope, mechanism and controls

Fixing one coefficient removes the predecessor's independent scale parameter;
it does not fix beta shape. The genuine row is the complete composition
polynomial F_15, with coefficients `(1,13,55,84,35,1)`. This model matches
its constant and linear coefficients only; higher coefficients vary with R.
Consequently this is **not** an actual path-count or Laurent-support
counterexample.

The exact ratio-256 control has the same three anchors and admissibility
certificates. At its original smallest root, in
`(760984/10^6,760985/10^6)`, its full Q is strictly negative. The sign change
under a shape operation is the new obstruction. Section 3 supplies the
additional continuous admissibility and original-root certificate needed to
turn it into exact double cancellation.

At the limiting shape, four B roots tend to zero while the fifth tends to
13. The full beta polynomial G_B approaches the degree-one row `1+13t`, and its
original vanished-coefficient phase escapes to infinity at rate R/338.
For each finite R all degrees and coefficients remain full. This degeneration
explains why fixing the coefficient sum leaves an unbounded obstruction;
it is not a legitimate replacement of the complete actual coefficient law.

The new derivative cone remains sound. At this witness, a proposed equality
of Q with a nonzero nonnegative combination of its same-zero windows is
impossible, because each generator is strictly negative while Q is positive.
Any successful actual-path certificate therefore has to use a further part
of the coefficient law; linear anchoring cannot supply it.

## 5. Verification status

The [standalone source](../../04-computation/open_frontier_sep06_laurent.py)
uses Fraction arithmetic, elementary-symmetric subset construction, literal
raw Laurent products and a separate ordinary-carrier convolution path.
Its finite universe is R=256 and R=512 at the displayed q,h, plus the actual
`(-9,1,6)` control, where `P(-2)=0,Q(-2)=-6305`. All stated sign and
root-isolation gates are exact. The floating-point discovery pass supplies
no proof input and no degree-minimality claim.

The symbolic family is constructed directly from subset exponents as
polynomials in R, then transformed by exact reciprocal substitution. No
interpolation or floating root finder is used. The certificate checks the
uniform interlacings, phase tube, derivative, transformed response identity
and its strict rectangle bound. The latter is an inequality on a continuous
rectangle, not a sample of ratios. Direct specializations of both symbolic
row identities are compared with the separate raw Laurent engine at four
rational phases, and an exact Euclidean algorithm proves the coprime
specialization used to make the resultant nonzero.

Current source/output: **203 exact gates PASS**. The
[output](open_frontier_sep06_laurent.out) records both original-zero
polynomials, all strict-interlacer residues and response intervals.
Normal and optimized runs are byte-identical to that frozen output.
Reproduce from the repository root:

```bash
python3 -B 04-computation/open_frontier_sep06_laurent.py
python3 -B -O 04-computation/open_frontier_sep06_laurent.py
```

SHA256:

```text
source   37fdc13870f7cf3122f7eebd413cbbfcf8df6a3451744027a8172b1ecdf5f8f0
output   c9d40353e33585d6c1c86d465e09f2862345c8e5af83744595d4357af015e8dd
semantic 096cdf5d1d82f5d5dd28d2805b365732312d66a6e6ca8b57c288dacdc6cdb247
```

Root independently read the complete proof and source and replayed both
203-gate runs: **PASS**. The separate
[analytic and source referee](open_frontier_sep06_laurent_audit.md) likewise
accepted the homogeneity, full carry, phase tube, all-weight proof,
continuous-parameter certificate, degree persistence, algebraic cancellation
and asymptotic constants, and independently reproduced the frozen output:
**PASS**. All source/output claims are restricted to this abstract model and
its exact declared parameter range.

The applied research cards were
[re-evaluate after a fibre-changing operation](../../00-navigation/META-PATTERNS.md#re-evaluate-a-certificate-after-every-fibre-changing-operation)
and [attack by degeneration](../../00-navigation/META-PATTERNS.md#attack-a-proposed-bound-before-extending-it).
The result supplies a counterindication to treating a repaired scalar anchor
as a full coefficient law. No shared meta-pattern file was edited.
