---
id: THM-3946
title: "Affine one-factor internal splits are reducible or have at least two infinity places"
status: >
  PROVED + VERIFIED-EXACT + CORRECTED AFTER MISTAKE-472; PENDING INDEPENDENT
  HOSTILE AUDIT.
  In the complete one-variable affine grammar, absorb the UFD split scalar
  into the actual opposite-side factors A and B.  If both are nonconstant,
  affine gauge gives A=Y and B=mY+c with a genuine slope ratio m.  For c=0,
  THM-3947 makes the full discriminant a reducible union of one-place
  line/parabola components.  For
  c!=0, every row is irreducible and its hidden plane cubic has one ordinary
  cusp, hence normalization P1.  It has exactly two infinity places at
  m=1 and m=-omega^2, and three otherwise; m=1 is the explicit quartic while
  the other rows are sextics.  The infinity polynomial has no quadratic term
  and constant -2, explaining the uniform no-one-place gate.  Exactly one
  constant factor reduces to THM-3942; two nonzero constant factors make the
  discriminant a reducible vertical cubic.  Thus no irreducible full
  discriminant in this entire affine internal-split grammar is one-place.
source: jc-degree6-place / full slope-ratio repair of the THM-3942--3947 internal-split lane, 2026-08-24
audit: >
  SELF-HOSTILE EXACT CANDIDATE.  The first equal-slope draft was rejected
  before promotion because the actual quotient-factor slope ratio is not a
  removable scalar.  The repaired companion retains m, verifies the complete
  affine gauge, the hidden cubic and its birational inverse, its unique cusp
  and two exact cusp-line normalizations, the corrected special tangent at
  m=-omega^2, the exact infinity polynomial and end count, the balanced
  quartic's basepoint-free normalization, and the corrected THM-3944
  order-regular character ledger.  Normal and optimized runs byte-match the
  frozen output at CHECKS=66.  The companion covers the complete
  nonconstant-factor affine gauge; MISTAKE-472 records the omitted
  constant--constant boundary, now closed directly in Section 5 without a
  changed exact certificate.  Independent audit remains required.
depends_on:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3947-scalar-weighted-repeated-square-split-trichotomy
related:
  - THM-3949-coprime-one-variable-internal-factor-splits-are-reducible-or-multi-ended
script: 04-computation/jc2_affine_internal_split_two_end_collision_thm3946.py
output: 05-knowledge/results/jc2_affine_internal_split_two_end_collision_thm3946.out
script_sha256: 37fe5a90e08c31ff71ba868e431542ab6a6dc2edd911eb9d8faaa8e76c5a5b21
output_sha256: 05289ab2b6520b7ef2ddf4d9caf69919ed80097141c9c8f9a7f75fa9498aac5d
semantic_sha256: 04eaeed97e485f1561d7e83a6462ae2147f6ec20bc96d96e621bb4fe5208924c
hash_basis: raw LF bytes
---

# THM-3946 -- affine internal splits never produce an irreducible one-place branch

**PROVED + VERIFIED-EXACT + CORRECTED AFTER MISTAKE-472; PENDING INDEPENDENT
HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
characteristic zero.  Fix

```text
omega^2+omega+1=0,          delta=omega-omega^2,          delta^2=-3. (1)
```

Let `P,Y` be independent affine coordinates and consider two torus rows

```text
H=q0^2-4p0^3=q1^2-4p1^3.                              (2)
```

Suppose one cube-difference factor is internally split between the two
complementary rows, while the other two factors are assigned whole.  After a
cube-root gauge on `p0`, exchanging the two sides if necessary, and absorbing
the UFD scalar into the factors actually appearing in the rows, this means

```text
p0=P,                         p1=P+A(Y)B(Y),
q1-q0=2A(Y)(p1-omega p0),
q1+q0=2B(Y)(p1-omega^2 p0),                         (3)
```

where `A,B` are affine polynomials in the same complementary coordinate.
Then the following cases are exhaustive.

1. If both `A,B` are nonzero constants, then `H` is a degree-three polynomial
   in `P` and hence is reducible over `k`.  A zero factor gives duplicate
   torus data.
2. If exactly one of `A,B` is a nonzero constant and the other is
   nonconstant, `(3)`
   is a whole-factor singleton split.  THM-3942 gives two infinity places.
3. If both are nonconstant and have the same zero, the full discriminant is
   reducible by THM-3947: generically it is three distinct parabolas; at
   `m=1` it is a doubled `p0` line plus a parabola; and at `m=-omega^2` it is
   a doubled `p1` parabola plus a distinct parabola.
4. If both are nonconstant and coprime, the full discriminant is irreducible.
   Its normalization is `P1` minus exactly two points for
   `m in {1,-omega^2}` and minus three points otherwise.  The `m=1` row is an
   irreducible quartic; every other row is an irreducible sextic.

Consequently no irreducible **full** discriminant in the affine one-factor
internal-split grammar is one-place.  The assertion deliberately distinguishes
the full discriminant from an individual component of a reducible divisor;
THM-3947 shows that such individual components can themselves be smooth
one-place lines or parabolas.

## 1. The slope ratio is a genuine coordinate

The scalar ambiguity in the UFD factorization should be absorbed before
normalizing coordinates.  Namely, define the actual quotient factors

```text
A=(q1-q0)/(2(p1-omega p0)),
B=(q1+q0)/(2(p1-omega^2 p0)).                         (4)
```

Write `A=aY+b` and `B=eY+f`, with `ae!=0`.  Replacing the complementary
coordinate by `Y'=A` gives

```text
A=Y',                 B=mY'+c,
m=e/a,                 c=f-eb/a.                       (5)
```

Thus `m` is not a removable factorization scalar; exchanging `A` and `B`
inverts it.  The first draft of this theorem incorrectly tried to force
`m=1`.  Formula `(5)` is the repaired exhaustive gauge.  Moreover, `c=0`
holds exactly when `A,B` have the same zero, while `c!=0` is exactly the
coprime affine-factor row.

Drop the prime from `Y'` and put

```text
A=Y,                 B=mY+c,              D=AB,
L1=p1-omega p0,      L2=p1-omega^2 p0.                   (6)
```

Equations `(3)` determine

```text
q0=B(D+(1-omega^2)P)-A(D+(1-omega)P),
q1=B(D+(1-omega^2)P)+A(D+(1-omega)P).                  (7)
```

Their product identity is

```text
(q1-q0)(q1+q0)=4D L1 L2=4(p1^3-p0^3),                 (8)
```

so `(2)` follows exactly.

## 2. The hidden branch equation is birational to the discriminant

On an irreducible component of `H=0` not contained in `P=0`, set

```text
h=q0/(2P).                                               (9)
```

Then `q0^2=4P^3` gives `h^2=P` and `q0=2h^3`.  Conversely these two equations
imply `H=0`.  Substituting `P=h^2` in `(7)` gives the exact hidden row

```text
C_{m,c}(h,Y)=
 m(m-1)Y^3+c(2m-1)Y^2
 +[c^2+{m(1-omega^2)-(1-omega)}h^2]Y
 +(1-omega^2)c h^2-2h^3=0.                            (10)
```

For `c!=0`, neither `P=0` nor `h=0` is a component, because

```text
C_{m,c}(0,Y)=Y(mY+c)((m-1)Y+c),
H(0,Y)=[Y(mY+c)((m-1)Y+c)]^2.                         (11)
```

Localizing at `P` and `h`, respectively, gives inverse ring maps

```text
k[P,Y]/(H)[P^-1]  <-->  k[h,Y]/(C_{m,c})[h^-1],
P |-> h^2,                    h |-> q0/(2P).            (12)
```

Since neither distinguished coordinate is a component, `(12)` proves the
equivalence

```text
H irreducible  <=>  C_{m,c} irreducible,                (12a)
```

as well as birationality.  This is the bridge that lets the exact
normalization and infinity points of `(10)` classify the full discriminant,
rather than merely an auxiliary cover.

## 3. Every coprime row is an irreducible cuspidal curve

Assume `c!=0`.  Scaling

```text
(h,Y,P)=(c hbar,c Ybar,c^2 Pbar)                       (13)
```

extracts `c^3` from `(10)`, so it is enough to set `c=1`.  Put

```text
alpha=delta(m+1)+3(m-1),       beta=delta+3,
M=m+omega^2.                                             (14)
```

Twice the hidden equation is the plane cubic

```text
F_m=-4h^3+alpha h^2Y+beta h^2
    +2m(m-1)Y^3+2(2m-1)Y^2+2Y.                         (15)
```

We now classify its singularity uniformly for every `m!=0`.  On `h=0`,

```text
F_m(0,Y)=2Y(mY+1)((m-1)Y+1),                           (16)
```

and each finite root is simple.  Away from `h=0`, the exact critical-line
identities are

```text
(F_m)_h=2h(-6h+alpha Y+beta),

F_m((alpha Y+beta)/6,Y)=2delta(MY+1)^3/9,
(F_m)_Y((alpha Y+beta)/6,Y)=2delta M(MY+1)^2/3.         (17)
```

If `M!=0`, these equations give the unique affine singular point

```text
S_m=(h,Y)=(-omega/M,-1/M).                              (18)
```

With `u=h+omega/M` and `v=Y+1/M`, the whole local cubic, not merely its first
terms, is

```text
F_m=Q_m+C_m,
Q_m=(6omega/M)(u-alpha v/6)^2,
C_m=-4u^3+alpha u^2v+2m(m-1)v^3.                       (19)
```

On the unique tangent `u=alpha v/6`,

```text
C_m=2delta M^3v^3/9 !=0.                               (20)
```

Thus `S_m` is a unibranch ordinary cusp and its tangent line is not a
component.

If `M=0`, the first restricted polynomial in `(17)` is the nonzero constant
`2delta/9`, so there is no affine singularity.  The cusp has moved to

```text
[h:Y:Z]=[1:omega^2:0].                                  (21)
```

In the chart `h=1`, put `Y=omega^2+v`.  The corrected exact expansion is

```text
F_{-omega^2}=
 3(1+delta)(v-delta Z/3)^2
 +2v(-v^2+delta vZ+Z^2).                               (22)
```

The cubic part on the tangent `v=delta Z/3` is
`2delta Z^3/9`, again nonzero.  The coefficient `delta/3` is important: an
earlier hostile calculation proposed `(1+delta)/3`, but direct substitution
rejects it.  Formula `(22)` is the repaired exact tangent.

These are the only projective singularities.  Indeed the infinity binary
cubic for `m notin {0,1}` is

```text
2phi_m(r)=2m(m-1)r^3+2[m(1-omega^2)-(1-omega)]r-4,     (23)
```

and

```text
disc_r(2phi_m)=-192delta m(m-1)M^3.                    (24)
```

Away from `m=-omega^2` its roots are simple.  At `m=1`, the remaining
infinity point `[0:1:0]` is smooth because `(F_m)_Z=2` there; at
`m=-omega^2`, `(22)` gives the sole cuspidal point plus one simple point.

Finally, every reducible plane cubic over an algebraically closed field has a
line component.  That line meets the residual conic at a projective singular
point.  Since the singular point above is unique and has one tangent, the line
would have to be that tangent; but `(20)` or `(22)` says the cubic restriction
to the tangent is nonzero.  This contradiction proves `F_m` irreducible for
every `m!=0`.  By `(12a)`, the full discriminant `H` is therefore irreducible
for every coprime affine row.  It is a sextic when `m!=1` and a quartic when
`m=1`.

The cusp-line construction also gives exact normalizations.  For `M!=0`, set

```text
K_m(s)=-4s^3+alpha s^2+2m(m-1),
v(s)=-(6omega/M)(s-alpha/6)^2/K_m(s),
h(s)=-omega/M+s v(s),            Y(s)=-1/M+v(s).       (25)
```

Substitution gives `F_m(h(s),Y(s))=0`, and

```text
K_m(alpha/6)=2delta M^3/9 !=0,                         (26)
```

so there is no hidden cancellation.  For `M=0`, lines `v=sZ` through `(21)`
give

```text
Z(s)=3(1+delta)(s-delta/3)^2/[2s(s^2-delta s-1)],
h(s)=1/Z(s),                     Y(s)=omega^2/Z(s)+s.   (27)
```

Both maps are birational, proving that the projective normalization is `P1`.
Their pole divisors are exact.  From

```text
disc_s(K_m)=-192delta m(m-1)M^3,                       (28)
```

`K_m` has three distinct roots for
`m notin {1,-omega^2}`.  At `m=1` it has one double and one simple root,
hence two distinct poles.  In `(27)`, the two infinity parameters are
`s=delta/3` and `s=infinity`.  Consequently

```text
# infinity places = 2,    m=1 or m=-omega^2;
# infinity places = 3,    otherwise.                   (29)
```

In particular, every coprime affine row is rational but multi-ended.  The
small no-triple-root mechanism remains visible directly in `(23)`: its
quadratic coefficient is zero and its constant is `-4`, so it cannot be one
cube of a linear factor.

## 4. Equal slopes give an explicit two-place quartic

Set `m=1` and retain `c!=0`.  Then

```text
q0=[Pc delta+3Pc+2P delta Y+2c^2Y+2cY^2]/2.            (30)
```

The hidden row becomes quadratic in `Y`:

```text
cY^2+(c^2+delta h^2)Y+c(delta+3)h^2/2-2h^3=0,          (31)
```

with exact discriminant

```text
disc_Y(31)=(c-h)^3(c+3h).                               (32)
```

Over `k(c)[h]`, the factors in `(32)` are coprime and have odd valuations.
Thus `(32)` is not a square, `(31)` is irreducible by Gauss, and `(12)` makes
the plane quartic `H` irreducible as well.

Removing the square `(c-h)^2` in `(32)` gives the smooth conic

```text
u^2=(c-h)(c+3h).                                       (33)
```

An exact normalization parameter is

```text
h= 2c(1-s)/(s^2+3),
u=-c(s-3)(s+1)/(s^2+3),
P= 4c^2(s-1)^2/(s^2+3)^2,
Y=-c(s-1)^2(s^2+2s+3+2delta)/(s^2+3)^2.               (34)
```

The rational inverse on a dense open is

```text
h=q0/(2P),
u=(2cY+c^2+delta h^2)/(c-h),
s=(u-c)/h.                                              (35)
```

In homogeneous normalization coordinates `[S:R]`, the target map is

```text
[P:Y:Z]=[
 4c^2(S-R)^2R^2,
 -c(S-R)^2(S^2+2SR+(3+2delta)R^2),
 (S^2+3R^2)^2].                                        (36)
```

These coordinates have no common zero for `c!=0`.  The target line at
infinity pulls back to

```text
S^2+3R^2=0,                                             (37)
```

which has two distinct roots, and the `P` numerator is nonzero at both.
Therefore the affine normalization is `P1` minus two points, hence `Gm`, and
the quartic has exactly two infinity places.  This is the equality boundary
of the general obstruction.

## 5. Collision and constant-factor boundaries

If both `A,B` are nonzero constants, formula `(7)` makes `q0` affine-linear
in `P`.  Therefore

```text
H=q0^2-4P^3 in k[P]
```

has degree three with leading coefficient `-4`.  Since `k` is algebraically
closed, it factors in `k[P]`, hence is reducible in `k[P,Y]`.  If either
factor is zero, the two torus rows are duplicates.  This is the
constant--constant boundary omitted from the first promoted case split.

If `c=0`, then `D=mY^2`.  Choose `G` with `G^2=mY^2` and put
`lambda=m^(-1/2)`.  The row is exactly THM-3947 with

```text
t=lambda^2=1/m.                                        (38)
```

That theorem proves the full discriminant always has three one-place
components counted with multiplicity.  Generically all three are parabolas;
the endpoint `m=1` has a doubled `p0` line and is THM-3944:

```text
H=-P^2(4P+3Y^2).                                       (39)
```

Its nonnormal quadratic order has regular locus `Gm^2`; the first Cardano
radicand gives the genuine boundary vector `(2,1)` there, while the second is
an exact cube.  The first vector does not extend across the conductor to the
full `A2` normalization.  The other exceptional value `m=-omega^2`
corresponds under `(38)` to `t=-omega` and is the swapped endpoint of
THM-3947.  Thus the collision boundary is reducible, not a hidden one-place
irreducible branch.

Finally, if exactly one of `A,B` is constant and the other is nonconstant,
then `p1-p0` is affine-linear and the internal assignment is actually a
whole-factor singleton split.  THM-3942 identifies its normalization with
`Gm` and again gives two ends.  This completes the affine case split.

## 6. Reproduction, scope and successor

Run

```bash
python3 04-computation/jc2_affine_internal_split_two_end_collision_thm3946.py
python3 -O 04-computation/jc2_affine_internal_split_two_end_collision_thm3946.py
```

The companion verifies the complete nonconstant-factor `m,c` coefficient
grammar, both torus rows,
the common discriminant, the hidden cubic and inverse, the unique cusp and
both cusp-line normalizations, the corrected exceptional tangent, the exact
two-versus-three end count, the balanced conic normalization and inverse, the
basepoint-free two-place projective map, the repeated-square seam, and the
corrected THM-3944 character ledger.  The constant--constant boundary is the
direct degree-three factorization at the start of Section 5 and does not
alter the 66-gate certificate or its hashes.

The theorem is exhaustive only for one internally split cube-difference
factor whose two assigned factors are affine polynomials in one complementary
coordinate.  It does not treat higher-degree coprime factors, simultaneous
internal splitting of two or three cube-difference factors, non-coordinate
`P,Y`, gcd overlap across distinct `Li`, or individual components of a
reducible full discriminant.  THM-3949 is the higher-degree one-variable
Newton boundary, where `(10)` persists with polynomial `A,B`.
