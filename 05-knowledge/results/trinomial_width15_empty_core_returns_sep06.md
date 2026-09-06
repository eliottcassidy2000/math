# An unbounded three-channel trinomial family beyond the endpoint-eight strip

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For every integer `g>=8` with `gcd(g,15)=1`, and every nonzero complex
coefficient triple, a Laurent polynomial supported on

```text
(-15, 2g-15, 3g-15)
```

has first nonzero constant-term moment exactly `g` or `2g`. Both values
are attained on every such support. At either scalar zero of the first
moment, the second moment divided by the square of the first channel's
coefficient monomial is strictly negative. The proof keeps the full
six-channel second row and its lower carry.

This is one fixed `(a,A,B)=(15,2,3)` family, not a theorem for every
support of smaller endpoint degree at most fifteen. The filename's
`width15` is shorthand for its sole-negative charge fifteen; its total
width is `3g`, which is unbounded.

## 1. Inheritance, prior coverage and the research move

The closest current mechanism is the
[carry-corrected Hermite trace formulation](trinomial_trace_sign_empty_core_certificate_sep06.md).
The [preceding root-sign probe](trinomial_root_sign_empty_core_returns_sep06.md)
proved normalized negativity for an unbounded endpoint-four family, then
identified the present family as its precise next symbolic test. The
canonical hostile `(-15,1,9)` has **positive** canonical second values,
while the essential factor `tau^-1` makes the first-anchor-normalized
values negative. Dropping that factor is the corrected near miss.
The least-used sidecar is the quotient norm's sign, which distinguishes
same-sign root values from bare coprimality.

Individual simple negative roots are supplied in much greater generality by
[THM-4436, complete factorial-row simple negative roots and trinomial phase wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
Here the first row is quadratic, so its exact discriminant proves the
needed real-root statement directly. The new cross-mass signs do not
follow merely from that individual-root theorem.

The prior [endpoint-eight theorem](overnight_20260906_moments_width8.md)
already introduced complete symbolic return rows, exact resultants, and
positivity after translating the unbounded parameter. Its family bank
has sole-negative charge at most eight, together with the finite
opposite-endpoint cases of positive endpoint at most eight. Every
support here has `min(15,3g-15)>=9`, so none belongs to that strip.
[THM-4432, two-channel two-rung noncancellation](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md)
has exactly two first channels, while this family has three. The earlier
[height-sixty coprimality census](synthesis_20260905_moments_trinomial.md)
already includes some specializations; those individual finite cases
are controls, not new discoveries.

Exact family/support strings, normalized tuple `(15,2,3)`, and the
distinctive certificate coefficients were checked against theorem and
relevant trinomial/moment result and producer paths on `origin/main`.
No prior unbounded proof for this family was recovered. This is an
addition to the recovered repository coverage, not a literature-priority
claim.

During the proof audit, incoming `e5a45df05f`, read on `origin/main` at
`533bf3d6e9`, supplied an independently audited
[221-support, 1,015-root sign bank](overnight3_20260906_moments_root_geometry.md).
Its Section 6 already retains the raw Laurent factor `tau^-epsilon_z`
and explicitly leaves the general negative-sign candidate open. Its
declared parameter bank includes this family's `g=8,11` controls.
Those signs and the carry correction are concurrent shared evidence,
not separate novelty here. The present contribution is the unbounded
analytic slice proved below. Incoming
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md)
is **RESERVED / UNPROVED EMPTY STUB** at that revision, with no proved
dependency; it is not used in this proof.

The live board is: full affine fibres; lower-carry normalization;
quadratic quotient algebra; trace versus determinant; unbounded parameter
positivity. The proof combines all five. Its decisive test is an exact
degree-seven polynomial sign, not a larger specialized gcd census.

## 2. Primitive support, full rows and the exact consumer

Let

```text
f(u)=alpha*u^-15+beta*u^(2g-15)+gamma*u^(3g-15),
alpha*beta*gamma!=0, g>=8 integral, gcd(g,15)=1.
```

The positive charges are distinct, and

```text
gcd(15,2g-15,3g-15)=gcd(15,g)=1,
gcd(15+(2g-15),15+(3g-15))=g.
```

A balanced multiplicity vector `(x,y,z)` of mass `m` satisfies

```text
g(2y+3z)=15m.
```

Thus `g|m`. At mass `g`, all nonnegative positive-charge counts solve
`2y+3z=15`; at mass `2g`, they solve `2y+3z=30`. Consequently the
**complete** rows are respectively

```text
(g-7+j, 6-3j, 1+2j),     j=0,1,2,
(2g-15+j, 15-3j, 2j),    j=0,...,5.
```

All first coordinates are positive for `g>=8`. There are three first
channels and six second channels, including the extra lower carry.
The first mass is exactly `g`, and there is no support-return mass
strictly between `g` and `2g`.

Set

```text
tau=alpha*gamma²/beta³,
X=alpha^(g-7)*beta⁶*gamma,
L(g)=(g)_5/720,
P_g(tau)=6tau²+20(g-5)tau+(g-5)(g-6),
Q_g(tau)=sum_(j=0)^5 (2g)_(15-j) tau^j/[(15-3j)! (2j)!].
```

The falling factorial notation is `(v)_k=v(v-1)...(v-k+1)`. The full
moment identities, with no suppressed zero-capable factor, are

```text
CT(f^g)=X*L(g)*P_g(tau),
CT(f^(2g))=X²*tau^-1*Q_g(tau).                    (1)
```

The canonical parameters are `(A,B,h,r,z,x)=(2,3,2,0,1,g-7)`.
Their affine step is `(1,-3,2)`, and the second left anchor equals
twice the first left anchor minus that step. This accounts exactly
for `tau^-1`; the upper carry is zero.

The leading coefficient of `P_g` is positive and its discriminant is

```text
8(g-5)(47g-232)>0.                                (2)
```

Its sum of roots is `-10(g-5)/3<0`, and its product is
`(g-5)(g-6)/6>0`; hence it has two distinct strictly negative roots.
Both are attainable on the coefficient torus: take `beta=gamma=1`
and `alpha=tau`. Also `L(g)>0` throughout the stated domain.

## 3. Exact anchored remainder and trace/norm factorization

Work modulo `P_g`. Since its constant coefficient is nonzero,

```text
tau^-1 = -[6tau+20(g-5)]/[(g-5)(g-6)].              (3)
```

Divide the full product in `(3)` times `Q_g` by `P_g`. The remainder
has degree one:

```text
R_g(tau)=tau^-1 Q_g(tau) mod P_g = C(g)+D(g)tau.
```

For a compact factorization define the positive factor

```text
K(g)=g(g-4)(g-3)(g-2)(g-1)
       *(2g-9)(2g-7)(2g-5)(2g-3)(2g-1),
U(g)=1106459g³-17466623g²+91449662g-158894736,
V(g)=27352253g³-404127357g²+1990079478g-3266244982.
```

Exact polynomial division gives

```text
C(g)=K(g)(g-5)U(g)/12259447200,
D(g)=K(g)V(g)/15324309000.                         (4)
```

The following computation is an identity over `Q(g)`, checked symbolically
and independently by exact interpolation with proved degree bounds;
it is not an empirical assertion from samples. Define

```text
J(g)=106089635g³-1564109559g²+7685968926g-12588295720,

H(g)=8051838961249g⁷-295063248863031g⁶
     +4630220682540559g⁵-40333487193583421g⁴
     +210638390555236696g³-659512619763975524g²
     +1146316900083491136g-853264618035013184.
```

Let `lambda_1,lambda_2` be the two first roots. Their sum and product
from Section 2 give exactly

```text
R_g(lambda_1)+R_g(lambda_2)
 =2C-(10/3)(g-5)D
 =-K(g)(g-5)J(g)/18389170800,                       (5)

R_g(lambda_1)*R_g(lambda_2)
 =C²-(10/3)(g-5)CD+(g-5)(g-6)D²/6
 =K(g)²(g-5)H(g)/3757351141239696000000.             (6)
```

These are the quotient multiplication trace and norm. Their only
nontrivial signs reduce to `J(g)` and `H(g)`.

## 4. Strict positivity closes every admissible parameter

Put `g=s+8`, `s>=0`. Expansion gives

```text
J(s+8)=106089635s³+982041681s²+3029425902s+3114337032,

H(s+8)=8051838961249s⁷+155839732966913s⁶
       +1288856301033727s⁵+5903575385111259s⁴
       +16172002313744184s³+26489396415060828s²
       +24017452577936640s+9296533901880000.
```

Every displayed coefficient is strictly positive. All linear factors
of `K(g)` and the factor `g-5` are positive for `g>=8`. Therefore
the trace in `(5)` is negative and the norm in `(6)` is positive.
The two remainder values are real. A negative sum and positive product
force both values to be strictly negative. Hence

```text
P_g(tau)=0  =>  tau^-1 Q_g(tau)<0.                 (7)
```

No equality or common-root case occurs in the domain. In particular
the complete first and second rows are coprime for every admissible
integral parameter. Substituting `(7)` in `(1)` shows that a vanishing
first moment is followed by a nonzero second one. Mass divisibility
then proves the exact detection law

```text
m_*(f)=g     if P_g(tau)!=0,
m_*(f)=2g    if P_g(tau)=0.                        (8)
```

The off-root case occurs, for example, at positive coefficients.
The two roots are attainable as noted above, so `2g` is the exact
worst detection order for each support. It is strictly below the
total exponent width `3g`; this is not an equality case for that width.

The same trace/norm sign identities hold for every real `g>=8` in
the falling-factorial extension. The first-return interpretation needs
the integer and coprimality hypotheses. Dropping the latter gives the
explicit hostile `g=9`: `(-15,3,12)` has a return at mass three,
not first mass nine. The displayed mass-nine row remains an identity;
its interpretation as the first row fails.

## 5. Controls, independent path and reproduction

The [standalone producer](../../04-computation/trinomial_width15_empty_core_returns_sep06.py)
uses exact rational symbolic arithmetic. It reconstructs both complete
falling-factorial rows, the inverse carry, remainder `(4)`, trace and
norm factorizations, discriminant, and every positive shift coefficient.
No repository producer or theorem is imported.

Seven named admissible parameters are `8,11,13,16,19,23,31`. For each,
an original charge-equation scan and independent repeated Laurent
multiplication agree on every coefficient of both rows. The latter
also checks every earlier mass. Specialized inverse/remainder and
quotient multiplication matrices independently check `(5)` and `(6)`,
as well as negative definiteness of the two-dimensional weighted trace
form. These are controls for the symbolic proof, not a census.

At `g=8`, the inherited hostile is recovered exactly:

```text
P_actual=56+560tau+56tau²,
Q_actual=16+10920tau+400400tau²+1681680tau³
         +720720tau⁴+8008tau⁵,
tau^-1 Q_actual mod P_actual=4743488+47087024tau.
```

Its Hermite form has negative leading entry `-461383264` and
positive determinant `587632760217600`, matching the earlier trace
companion. The nonprimitive `g=9` hostile is independently replayed.

```sh
python3 -B 04-computation/trinomial_width15_empty_core_returns_sep06.py
python3 -B -O 04-computation/trinomial_width15_empty_core_returns_sep06.py
```

The [frozen output](trinomial_width15_empty_core_returns_sep06.out)
records **113 explicit gates**. Normal and optimized outputs are
byte-identical. Semantic digest:
`155c07c958f6a6087266d09059d53ee7603d59b469fc8484795bac7e3c56e8ab`.

Raw SHA-256:

```text
source 28793dd2d8de3b45dfdb11dd8f3bf00888122a20f628b44bb9ad1245ca167b9c
output 58f8cbe33fbf671a0bac9e9a7f086f05c50e9a366cf715e47a3d44d9d397015a
```

The [independent audit](trinomial_width15_audit_empty_core_certificate_sep06.md)
uses a separate
[standard-library Fraction producer](../../04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py)
and [output](trinomial_width15_audit_empty_core_certificate_sep06.out),
importing neither SymPy nor the primary source. It first proves the
polynomial cancellation in the inverse-carry term after removing
`K=(2g)_10/32`. The normalized constant, linear, trace and norm have
degrees at most `4,3,4,8`. Nine exact recurrence evaluations at distinct
`g=8,...,16` therefore certify all four polynomial identities for every
parameter. Nonprimitive evaluation points are used only for those
identities. A separate binomial shift checks all twelve positive
coefficients, and native Laurent multiplication checks complete rows and
earlier emptiness at the four primitive parameters `8,11,13,14`.
All **69 gates** pass; normal and optimized bytes agree.

```sh
python3 -B 04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py
python3 -B -O 04-computation/trinomial_width15_audit_empty_core_certificate_sep06.py
```

Independent audit hashes:

```text
source ca0a5416a4e2c92d590fc3e701a9ac5d3a26ded1535d81d839eeffd2d834d83a
output 45c71b8626b1efc1aa3b97b58b18ba0482591db0c0ace7fd3bb7df1f6579c622
```

The certificate lane's independent analytic, identity and final-text
audits are **PASS**.
Root separately read the complete proof and primary source, checking mass
divisibility, complete rows, the monomial carry ratio, discriminant,
trace/norm signs and the exact detection consequence: **PASS**.

## 6. Connection contract and next boundary

The source is a complete weighted first/doubled affine fibre with one
lower carry. The target is a pair of scalar signs in a quadratic
quotient. The map first removes nonzero coefficient monomials, retains
`tau^-1`, and reduces the complete second row modulo the first.
It preserves every zero/nonzero moment and the anchored sign at each
first zero. Taking the norm alone loses the individual signs; the trace
is the needed second scalar. Individual real-rootedness supplies neither.

The polynomial translation `g=s+8` converts these signs into twelve
positive coefficients, recovering the earlier endpoint-eight symbolic
operation on a new fixed family. The full weighted fibres and the carry
are still required; pair compression would change the first equation.
The session's empty-hexagon source motivates retaining the validity
sidecar while compressing a certificate at the method level only.

The next structural boundary is **OPEN**: a cubic first row, for example the
already named `(-21,1,12)` extended to
`(-21,2g-21,3g-21)`, `g>=11`, `gcd(g,21)=1`. It has first degree three,
second degree seven, and one lower carry. A trace and norm no longer
determine all three signs; the full signed leading minors of its
Hermite form are needed. This is a precise next proof obligation,
not a claim that the present quadratic argument already extends.
