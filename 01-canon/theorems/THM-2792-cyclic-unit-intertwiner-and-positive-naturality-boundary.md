---
id: THM-2792
title: "Cyclic unit intertwiner and positive-naturality boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every oriented
  determinant cycle of the THM-2625 endpoint current and THM-2791's
  transferred semantic chain are cyclic generators of regular K[C13]
  modules.  They therefore determine one unique origin-free module
  isomorphism F(rA)=rJ.  In chosen origins its kernel is T=J A^(-1), with
  exact relative-origin covariance.  The normalized semantic inverse has
  Smith boundary Z/2, and the actual kernel has a nontrivial class modulo
  scalar monomials and admits no bi-positive circulant realization.  This
  is coefficient-module descent, not one-way positivity, frame naturality,
  pointwise ancestry, physical Cech descent, row exclusion, or LRC(14).
source: root/central-cycle-unit-intertwiner-2026-07-28
audit: >
  root/audit-2792-2026-07-28 independently rederived the Fourier unit
  inference, inverse/Smith/integrality boundary, coordinate-free map and
  origin covariance, scalar-monomial and bi-positive obstructions, formal
  frame hostile, lower-central quotient loss, and Cech scope; exact
  companion replay, hashes, and documentation checks were accepted.
depends_on:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2647-endpoint-anchored-two-point-deconvolution-and-thirteen-halves-signed-tax
  - THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
related:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
script: 04-computation/lrc14_cyclic_unit_intertwiner_thm2792.py
output: 05-knowledge/results/lrc14_cyclic_unit_intertwiner_thm2792.out
script_sha256: a8108344826409a43d615a2d9cb88330096af2b910f8ede1c8a0261d5633aba0
output_sha256: 9a6350b38b8d9f92367129f894357be15956ea1a6c67a240f32fa603bab86834
hash_basis: LF-normalized bytes
---

# THM-2792 -- cyclic unit intertwiner and positive-naturality boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The full endpoint-cycle spectrum from THM-2625 and THM-2790 and the
two-point semantic transfer from THM-2791 solve the split cyclic
coefficient-algebra problem exactly.  Each side is a cyclic generator, so
there is a unique equivariant coefficient-module isomorphism between them.
That isomorphism is independent of cycle origins even though its displayed
convolution row is not.

The theorem also isolates the exact failure boundary.  The isomorphism is
not induced by a pointwise torsor identification, it has no bi-positive
circulant realization, its natural integral model has an unavoidable
factor-two obstruction, and the quotient algebra forgets the lower-central
physical chord.  No physical ancestry or LRC(14) conclusion follows.

## 1. Exact coefficient field and cycle polynomial

Let

```text
K=Q(zeta_N),      N=50,334,435,734,703,120,
```

the characteristic-zero cyclotomic field of THM-2625.  It contains
`zeta=zeta_13` and all endpoint coefficients.

For every `s in F_13^2\{0}` and `Delta in F_13`, choose one origin `R_0`
with

```text
det(s,R_0)=Delta.
```

The associated endpoint cycle is

```text
R_j=R_0+j s,                    j in F_13,
```

and its coefficient polynomial is

```text
J_(s,Delta,R_0)(z)
 =sum_(j=0)^12 P_(R_j+s) Q_(R_j) z^j
 in K[C_13]=K[z]/(z^13-1).                               (1)
```

There are

```text
(13^2-1)*13=2,184                                      (2)
```

oriented determinant cycles before origins are chosen.

## 2. Every endpoint cycle is a group-algebra unit

Because `K` contains every thirteenth root and has characteristic zero,
Fourier evaluation is an algebra isomorphism

```text
K[C_13] -> product_(k=0)^12 K,
F(z) |-> (F(zeta^(-k)))_k.                              (3)
```

Consequently,

```text
F is a unit  <==>  F(zeta^(-k))!=0 for every k.          (4)
```

At `k=0`, equation `(1)` gives

```text
J_(s,Delta,R_0)(1)
 =sum_(R:det(s,R)=Delta) P_(R+s)Q_R
 =Sstar(s,Delta).                                        (5)
```

This is nonzero by THM-2625 because `s!=0`.  For each
`k=1,...,12`, equation `(1)` evaluated at `zeta^(-k)` is exactly
THM-2790's central mode

```text
Jhat_(s,Delta)(k).                                       (6)
```

THM-2790 proves all

```text
168*13*12=26,208                                        (7)
```

values in `(6)` nonzero in characteristic zero.  Combining the `2,184`
zero modes from `(5)` with `(7)`, all

```text
2,184*13=28,392                                         (8)
```

Fourier coordinates are nonzero.  Therefore every
`J_(s,Delta,R_0)` is a unit of `K[C_13]`.

Changing the named origin only multiplies `J` by a monomial unit, so this
unit statement is origin-independent.

## 3. Semantic unit and the exact integral boundary

THM-2647 supplies the general odd-cycle inverse, determinant-two index,
and unavoidable `13/2` signed tax.  Here that mechanism specializes to
THM-2791's semantic transfer, and the coefficientwise integrality boundary
can be stated exactly.

THM-2791's transferred semantic-arm chain is

```text
A=c z^6(1+z),
c=790,161,473,087,466,480!=0.                            (9)
```

Since

```text
(1+z)(1-z+z^2-z^3+...+z^12)=1+z^13=2,                  (10)
```

one has in `K[C_13]`

```text
A^(-1)
 =1/(2c) z^7(1-z+z^2-z^3+...+z^12).                    (11)
```

Thus `A` is a unit.  If the harmless nonzero scalar is removed and

```text
A_0=z^6(1+z),
```

then the only denominator in its inverse is exactly two.

This denominator is sharp integrally.  Multiplication by `1+z` on
`Z[C_13]` has

```text
det(I+shift)=2.                                          (12)
```

Deleting row zero and column zero leaves a determinant-one triangular
minor.  Hence the Smith form is

```text
diag(1,...,1,2),                                         (13)
```

with twelve ones, and the integral cokernel is exactly `Z/2`.

There is also a coefficientwise criterion for a particular normalized
intertwiner.  Put

```text
S=1-z+z^2-z^3+...+z^12.
```

For `J in O_K[C_13]`,

```text
J A_0^(-1)=1/2 J z^7 S
```

is integral if and only if `J(1)` lies in `2 O_K`.  Modulo `2`, all
coefficients of `S` equal one, so every coefficient of `JS` is congruent
to `J(1)`; this proves both directions.  The endpoint theorems establish
`J(1)!=0`, not divisibility by two, so no uniform endpoint-lattice claim
follows.  The Boolean cube splits sharply into `4,096` integral and
`4,096` genuinely half-integral normalized coefficient solutions.  These
Boolean controls need not themselves be units.

For raw `A`, the determinant is `2c^13` and `(11)` has denominator `2c`.
Passing to the normalized coefficient line divides out `c` over `K`; it
does not remove the factor two or create a positive inverse.

## 4. Unique coordinate-free coefficient intertwiner

Fix a nonzero endpoint direction `s`, a determinant value `Delta`, and the
orientation which identifies the abstract generator of

```text
G=C_13
```

with the semantic central successor on the source and translation by `s`
on the endpoint cycle.  Let `V_A,V_J` be the two resulting regular
`K[G]`-modules.  The semantic chain `A in V_A` and endpoint cycle
`J in V_J` are cyclic generators by Sections 2--3.  Therefore

```text
mu_A:K[G]->V_A,             r |-> r A,
mu_J:K[G]->V_J,             r |-> r J                    (14)
```

are isomorphisms.  Define

```text
F_(s,Delta)=mu_J mu_A^(-1),                              (15)
F_(s,Delta)(r A)=r J.
```

This is the unique `K[G]`-module isomorphism satisfying

```text
F_(s,Delta)(A)=J.                                        (16)
```

Equation `(15)` is coordinate-free: it requires the oriented group
actions but no choice of either cycle origin.

After origins are chosen, both modules are written as `K[C_13]` and the
same abstract map is represented by convolution with

```text
T=J A^(-1).                                              (17)
```

Its coordinate inverse is `A J^(-1)`.  Writing
`Jhat(k)=J(zeta^(-k))`, the endpoint inverse is

```text
J^(-1)
 =sum_(j=0)^12 [
    1/13 sum_(k=0)^12 Jhat(k)^(-1) zeta^(k j)
   ] z^j.                                               (17a)
```

Thus every one of the `2,184` oriented determinant cycles has a unique
origin-free `K[C_13]`-module isomorphism from the semantic regular module
to the endpoint regular module carrying `A` to `J`.  In chosen origins,
its matrix is the unique circulant kernel `(17)`.

Here uniqueness means uniqueness among `K[C_13]`-linear maps.  Many
non-equivariant linear automorphisms send one nonzero vector to another.

The coordinate coefficients of `T` lie in `(2c)^(-1) O_K` because `J`
has cyclotomic-integral coefficients.  This does not make `F` an
automorphism of the natural integral lattices: THM-2625 and THM-2790 prove
the Fourier coordinates of `J` nonzero, not algebraic units.  The inverse
may require additional cyclotomic denominators.

## 5. Exact origin and frame covariance

Let the endpoint origin advance by `h_E` and the semantic origin by `h_A`.
Their coordinate rows change by

```text
J -> z^(-h_E)J,                 A -> z^(-h_A)A,
T -> z^(h_A-h_E)T.                                      (18)
```

Equation `(18)` is the basis-change law for the same abstract map
`F_(s,Delta)`.

1. A simultaneous common shift `h_A=h_E` leaves the displayed circulant
   row unchanged.
2. An endpoint-only shift gives `T->z^(-h_E)T`.  Since `T` is a unit, the
   two coordinate rows differ for `h_E!=0`; this says only that one
   coefficient vector cannot be invariant under an independent target
   basis change.  It is not an obstruction to `(15)`.

There is a formal algebraic covariance.  If both coefficient coordinates
are simultaneously reindexed by

```text
alpha_(u^(-1)): z->z^(u^(-1)),                           (19)
```

then

```text
T->alpha_(u^(-1))(T).                                    (20)
```

It would be false to identify `(19)` with the physical frame change
`s->u s`.  The current itself contains the increment:

```text
J_s(R)=P_(R+s)Q_R,
J_(us)(R)=P_(R+us)Q_R.                                   (21)
```

Even after the same geometric line is re-enumerated, the left endpoint in
`(21)` changes from a one-step edge to a `u`-step edge.  The two arrays
obey no universal coefficient-permutation identity, even after an
arbitrary cycle-origin shift.  This is a formal boundary: it does not
claim that the canonical current has been exhaustively computed for every
accidental cross-frame equality.  Orientation reversal is included
because the separate `P,Q` factors are ordered and need not be symmetric.

Thus a THM-2779 frame change changes the geometric central edge, not just
its coordinates.  THM-2790 proves the unit predicate for every new
direction but supplies no transport of coefficient values between
directions.  Equation `(15)` is a `2,184`-member oriented-cycle bank, not
one canonical `SL_2(F_13)`-natural operator.

## 6. Sharp hostile boundaries

### 6.1 Fourier support does not imply positivity

Take the abstract endpoint control

```text
J=delta_0.
```

All thirteen Fourier coordinates equal one.  With normalized `A_0`, the
unique intertwiner is

```text
T=A_0^(-1).
```

Its coefficient vector has seven positive and six negative entries, all
of magnitude `1/2`.  Complete Fourier support and invertibility therefore
do not imply a positive, stochastic, Boolean, or mass-preserving
transport.  In particular, a nonnegative circulant inverse of `1+z`
cannot exist: the equations

```text
b_j+b_(j-1)=delta_(j,0)
```

force the unique alternating half-vector.

### 6.2 No pointwise torsor or bi-positive realization

Let

```text
M=K^x <z> subset K[C_13]^x.
```

A fixed-action equivariant bijection of two thirteen-point torsors
linearizes, after origins are chosen, to a monomial `z^h`; allowing one
common nonzero weight gives an element of `M`.

The semantic vector `A` has support two in every origin.  THM-2790 gives
support thirteen for the actual endpoint vector `J` on every determinant
cycle.  Hence `(15)` cannot be the linearization of a fixed-action torsor
bijection, even after scalar weighting: if `T=lambda z^h`, then `TA`
would still have support two.

Equivalently, the actual class

```text
[T] in K[C_13]^x / (K^x <z>)                            (22)
```

is nontrivial.  This statement is origin-independent.  The thirteen
abstract torsor identifications classified by THM-2611 still exist, but
none realizes the actual coefficient-module map.

There is a sharp ordered-cone strengthening.  Let `B,C` be real
coefficientwise-nonnegative kernels with

```text
B C=delta_0.
```

For every `x in supp(B)` and `y in supp(C)`, the convolution coefficient
at `x+y` contains the positive term `B_x C_y`.  With no cancellation and
a product supported only at zero, every such pair has `x+y=0`.  Fixing
either support point forces both supports to be singletons.  Therefore a
circulant isomorphism whose kernel and inverse preserve the standard
nonnegative cone is exactly a positive scalar monomial.

The actual `F` is not scalar-monomial, so it has no bi-positive circulant
realization in any origin gauge.  This does not exclude a one-way
nonnegative kernel with signed inverse.  The normalized chain
`A_0=z^6(1+z)` is itself the sharp control: it is nonnegative and
invertible, while its inverse has seven positive and six negative
coefficients.

### 6.3 Formal frame rescaling is not cycle reindexing

On an abstract thirteen-cycle independently tag endpoints by

```text
P_j=2^j,                         Q_j=3^j.
```

For every `u=2,...,12`, the physical `u`-step edge is

```text
J_u(j)=P_(u(j+1))Q_(uj).
```

After formally reindexing `J_1` by `j->uj` and shifting the cycle origin
by `h`, the comparison edge is

```text
P_(u(j+h)+1)Q_(u(j+h)).
```

Unique factorization forces `uh=0` from the powers of three, hence `h=0`,
and then `u=1` from the powers of two, a contradiction.  All
`11*13=143` nonidentity scaling/origin pairs fail.  This is the first
failed implication in a formal frame-natural identity; it is not an
exact canonical-current cross-frame census.

### 6.4 Quotient algebra forgets the lower-central class

In `Z[C_(13^5)]`, the two physical chains

```text
1+u,
1+u^53028
```

have the same pushforward `1+z` to `C_13`.  THM-2791 proves

```text
53028=1+13*4079,               4079=10 mod13,            (23)
```

so the second has first transgression `10 in Z_2/Z_3`, while the first is
the pure `Z_1` edge.  The quotient intertwiner `(15)` sees no difference.
It cannot reconstruct physical descent, same-high-digit ancestry, or the
lower-central chord.

### 6.5 Formal coefficient equivalence is not physical identification

The multiplication action of `K[C_13]^x` on itself is simply transitive.
Once `A` and `J` are units, `T=JA^(-1)` exists for formal algebraic
reasons.  Replacing `J` by any unrelated unit, such as `delta_0` or
`2+z`, produces the same kind of unique module map.  Coordinate-free
uniqueness is therefore not evidence of shared ancestry, endpoint phase,
or a physical carrier.  Such a conclusion needs an additional reduction
of structure group, naturality square, or allocation square.

## 7. Strongest valid holotopy/descent statement

Equation `(18)` gives an exact two-stage descent law.  The pair of origin
gauges `C_13 x C_13` acts on coordinate kernels through

```text
(h_E,h_A) |-> z^(h_A-h_E).                               (24)
```

The diagonal subgroup is the kernel, and the residual relative-origin
quotient is one `C_13`.  Coordinate rows descend under common
re-origining and transform covariantly under the remaining gauge.  They
are presentations of the single global module map `(15)`, not thirteen
distinct abstract intertwiners.

The origin-free fixed-action point-transport obstruction is the class
`(22)`.  Origin changes and scalar normalization do not alter it, and
Section 6.2 proves it is nontrivial.  In bundle language, the coefficient
isomorphism has structure group `K[C_13]^x`, but it does not reduce to
the scalar-monomial subgroup that linearizes pointwise torsor transport.
On a real ordered form, the positive part of that subgroup is exactly the
bi-positive circulant automorphisms.  The quotient class does not decide
whether a merely one-way nonnegative representative exists.

A genuine Cech statement requires physical data absent here.  If a future
base cover supplied local semantic and endpoint origins, their relative
shifts

```text
h_ij=(endpoint shift)_(ij)-(semantic shift)_(ij)
```

would satisfy

```text
h_ij+h_jk=h_ik,
T_j=z^(-h_ij)T_i.                                        (25)
```

These monomials are the coordinate transition cocycle of the already
global `F`; they do not establish a nontrivial physical Cech class.  If
the future cover instead supplied local pointwise torsor reductions,
their exponent transitions would form the `C_13` cocycle of THM-2611: a
cyclic global section would exist exactly at zero holonomy and then would
have thirteen choices.  No such cover, local reduction, overlap map, or
holonomy cycle is constructed by THM-2625, THM-2790, or THM-2791.

Thus coefficient descent succeeds, while reduction to scalar-weighted
fixed-action point transport fails locally.  Any later Cech obstruction
begins only after a different packet supplies local physical reductions.
More general noncirculant or one-way-positive physical maps remain outside
this no-go.

## 8. Source, target, map, and missing sidecar

```text
source:
  THM-2791's positive transferred central chain A, a cyclic generator over K;

target:
  one oriented THM-2625/2790 endpoint determinant-cycle generator J;

map:
  the unique coordinate-free K[C_13]-module map F(rA)=rJ,
  represented in chosen origins by T=J A^(-1);

preserved:
  the complete thirteen-character coefficient module, cyclic translation,
  independent-origin covariance, and formal simultaneous coefficient
  reindexing;

lost or absent:
  a set-level, bi-positive, or integral realization; any decision on
  one-way positivity; naturality across physical frame directions;
  physical high digits; lower-central class; semantic ancestry;
  root/Bockstein identification; and endpoint chronology;

needed sidecar:
  a same-ancestry physical allocation reducing the coefficient map to a
  lawful positive carrier while identifying the THM-2791 graded chord germ
  with one THM-2625 endpoint cycle before quotient transfer.
```

This closes the split cyclic coefficient-algebra problem.  It does not
close the physical naturality square.

## 9. Exact companion and audit

Run

```bash
python 04-computation/lrc14_cyclic_unit_intertwiner_thm2792.py
python -O 04-computation/lrc14_cyclic_unit_intertwiner_thm2792.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_cyclic_unit_intertwiner_thm2792.out.
```

The companion verifies `(11)`, the determinant and primitive minor in
`(12)--(13)`, the exact integrality criterion on all `8,192` Boolean
coefficient vectors, a nontrivial rational unit/intertwiner control, all
`169` independent source/target origin pairs, all `13`
scalar-monomial support controls, all `8,191` nonempty bi-positive support
patterns, every formal unit reindexing, all `143` tagged physical
frame-step hostiles, the alternating-sign hostile, the lower-central
transfer hostile, and the exact endpoint Fourier census.  The endpoint
nonvanishing inputs themselves are inherited from THM-2625 and THM-2790.

The stored transcript retains its original `SCRATCH AUDIT` banner so the
independently audited bytes remain immutable; the theorem status and
frontmatter are the canonical truth surface.

```text
script_sha256 = a8108344826409a43d615a2d9cb88330096af2b910f8ede1c8a0261d5633aba0
output_sha256 = 9a6350b38b8d9f92367129f894357be15956ea1a6c67a240f32fa603bab86834
hash_basis    = LF-normalized bytes
```

The independent audit rederived every load-bearing algebraic step,
corrected the origin-covariance interpretation, strengthened the
set-level boundary to the scalar-monomial quotient and bi-positive no-go,
repaired the formal frame hostile, and confirmed that a Cech obstruction
cannot be asserted without local physical reductions.

No positivity decision beyond the bi-positive boundary, physical endpoint
allocation, row exclusion, or LRC(14) conclusion is claimed.  QED.
