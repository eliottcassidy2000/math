---
id: THM-2848
title: "Whitened moving-plane multipoles and the strict Pearson boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Whitening turns binary quadratic divisibility into vanishing of one
  pure complex moment.  When the cubic and quartic pure moments vanish,
  the Pearson residual is an explicit nonnegative real conic, strictly
  positive for factorial polynomials.  General positive laws admit sharp
  fifth-root, two-radius, six-atom, exponential-circle, and spherical
  multipole hostiles.  On a factorial four-slot mean-zero space, harmonic
  projection reduces the bad locus to a vanishing quartic multipole or a
  shared cubic--quartic real multipole line; the exact resultant is a
  product of twelve squared cross products.  Integration by parts exposes
  a derivative witness, but the scalar moment sequence does not provide
  the required multiplier/lowering access.
source: root/gmc-four-slot-moving-plane-boundary-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
related:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
script: 04-computation/gmc_whitened_moving_plane_multipole_thm2848.py
output: 05-knowledge/results/gmc_whitened_moving_plane_multipole_thm2848.out
script_sha256: 08416715a0bfd66d66ad7e650c1aec0573fab457c4d70d29b1b17c3d13cd6bc4
output_sha256: 3855c15afe796a40997a135681a7ec56a03b4a9d9f658a476514b1e7ed5c1e35
hash_basis: LF-normalized bytes
---

# THM-2848 -- whitened moving-plane multipoles and the Pearson boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem has two deliberately separate layers.

1. For an arbitrary conjugation-compatible positive probability
   functional, whitening identifies binary divisibility exactly and gives
   a nonnegative Pearson conic.  Positivity alone does **not** exclude the
   bad plane; four explicit probability laws mark different failure
   boundaries.
2. For the factorial functional, faithful exponential-square positivity
   makes the Pearson conic strict.  In a four-slot mean-zero space, real
   harmonic multipoles give an exact global resultant criterion and
   integration by parts gives a transverse derivative witness.  The
   original scalar moments do not contain that multiplier, so this is a
   reduction and boundary theorem, not four-slot SFC or GMC closure.

## 1. Whitening makes divisibility a pure moment

Let `L` be a conjugation-compatible positive probability functional on a
commutative complex `*`-algebra:

```text
L(1)=1,              L(conj(f))=conj(L(f)),
L(f conj(f))>=0.                                             (1)
```

Let `E` be a centered real two-plane on which

```text
q(H)=L(H^2)                                                 (2)
```

is positive definite.  Choose an oriented conformal basis `U,V`:

```text
L(U^2)=L(V^2)=rho/2,       L(UV)=0,       rho>0.             (3)
```

Put

```text
Z=U+iV,             eta=x+iy,
H_eta=xU+yV=(conj(eta)Z+eta conj(Z))/2.                      (4)
```

Then

```text
L(Z)=L(Z^2)=0,       L(Z conj(Z))=rho,
q(eta)=rho |eta|^2/2.                                       (5)
```

For `m>=2`, write `p_m(eta)=L(H_eta^m)`.  Exact expansion gives

```text
2^m p_m
 =sum_(r=0)^m binom(m,r)
   conj(eta)^r eta^(m-r)
   L(Z^r conj(Z)^(m-r)).                                    (6)
```

The quadratic `q` is a nonzero scalar multiple of
`eta conj(eta)`.  Every interior term of `(6)` is divisible by that
product.  Its two endpoint coefficients are `L(Z^m)` and its conjugate.
Therefore

```text
q divides p_m
 iff
L(Z^m)=0.                                                   (7)
```

This equivalence is invariant under the remaining oriented gauge
`Z -> e^(i theta)Z`: the pure moment is merely multiplied by
`e^(im theta)`.

## 2. Cubic and quartic quotients

Assume from now through Section 4 that

```text
L(Z^3)=L(Z^4)=0.                                            (8)
```

Define

```text
a=L(Z^2 conj(Z)),
b=L(Z^3 conj(Z)),
s=L(|Z|^4).                                                 (9)
```

By `(7)`, both `p_3` and `p_4` are divisible by `q`.  Direct expansion
gives the exact quotients

```text
ell(eta):=p_3/q
 =3[a conj(eta)+conj(a)eta]/(4rho),                         (10)

R(eta):=p_4/q
 =3s|eta|^2/(4rho)
  +[b conj(eta)^2+conj(b)eta^2]/(2rho).                     (11)
```

For each real `H=H_eta`, the polynomial

```text
H^2-q(eta)-ell(eta)H                                       (12)
```

is orthogonal to both `1` and `H`.  Consequently

```text
L[(H^2-q-ell H)^2]
 =q [R-q-ell^2].                                           (13)
```

The quotient in brackets is the real conic

```text
T(eta)=A|eta|^2+C conj(eta)^2+conj(C)eta^2,                 (14)

A=[6rho s-4rho^3-9|a|^2]/(8rho^2),
C=[8rho b-9a^2]/(16rho^2).                                 (15)
```

Thus arbitrary positivity proves the sharp non-strict inequality

```text
6rho s-4rho^3-9|a|^2
 >=
|8rho b-9a^2|.                                             (16)
```

Equality is possible; Section 4 gives a six-atom witness.

The following conditions are equivalent:

```text
T(eta)>0 for every eta!=0;
6rho s-4rho^3-9|a|^2 > |8rho b-9a^2|;
{1,H,H^2} is L^2-independent for every nonzero H in E.     (17)
```

In particular, `(17)` holds for a real nonconstant polynomial `H` under
the factorial integral

```text
L(f)=integral_0^infinity f(t)e^(-t)dt.                     (18)
```

Indeed, equality in `(13)` would make `(12)` zero almost everywhere,
hence identically as a polynomial.  A nonconstant polynomial cannot
satisfy a quadratic equation over the constants.

Two useful but weaker Bessel bounds follow by projecting `|Z|^2` onto
the orthogonal triples `{1,Z,conj(Z)}` and
`{1,Z^2,conj(Z)^2}`:

```text
s >= rho^2+2|a|^2/rho,
s^2 >= rho^2 s+2|b|^2.                                    (19)
```

They do not imply the phase-sensitive conic `(16)`.

## 3. The exact information carried by the conic

On a fixed two-plane, simultaneous divisibility of the quadratic,
cubic, and quartic moment forms is exactly

```text
L(Z^3)=L(Z^4)=0.                                           (20)
```

The coefficients `a,b,s` do not obstruct `(20)`; they describe how
strongly the real directions surrounding its two complex null points
leave the null locus.  Formula `(16)` is therefore a transverse
positivity condition, not a proof that the complex null points are
absent.

For factorial polynomials the strict inequality in `(17)` says that a
bad plane, if it exists, has a strictly positive real Pearson conic
around it.  This is consistent with THM-2846: its exact moment-three
hostile has nonzero fourth moment, so it leaves the cubic null locus at
the next rung.  The conic does not by itself force that exit for an
arbitrary bad plane.

## 4. Four sharp positive-functional hostiles

These examples all have a positive definite real plane
`E=span_R(Re Z,Im Z)` and

```text
L(Z^m)=0,                      1<=m<=4.                  (21)
```

They are not factorial-polynomial counterexamples.

### 4.1. The support-minimal fifth-root law

Let `Z` be uniform on the fifth roots of unity.  Then `(21)` holds and

```text
rho=s=1,              a=b=0.                             (22)
```

Five support points are minimal among positive laws on the unit circle
satisfying `(21)`.  If `z_1,...,z_r`, `r<=4`, were the support, the
support polynomial

```text
P(z)=product_(j=1)^r (z-z_j)
```

would vanish almost surely.  But `(21)` makes `L(P(Z))` equal its
nonzero constant term `(-1)^r product_j z_j`, a contradiction.

### 4.2. A strictly positive two-radius law

Give `R=1,2` probability `1/2` each.  Conditional on `R`, give the angle
the density

```text
[1+epsilon(h_1(R)cos theta+h_2(R)cos 2theta)]/(2pi),
0<epsilon<=1/4,                                           (23)

h_1(1)=h_2(1)=1,
h_1(2)=-1/2,             h_2(2)=-1/4.
```

The density is strictly positive.  Its first two Fourier modes cancel
between the radii, and it has no third or fourth mode.  Hence `(21)`
holds, while

```text
rho=5/2,       s=17/2,       a=b=-3epsilon/4.            (24)
```

Thus nonconstant radius and strictly positive angular densities do not
exclude the bad plane.

### 4.3. The six-atom Pearson equality

Let `U` and `V` be independent, with

```text
P(U=1)=P(U=-1)=1/2,
P(V=0)=4/5,
P(V=sqrt(5))=P(V=-sqrt(5))=1/10,                         (25)
```

and put `Z=U+iV`.  This is a six-atom law and `(21)` holds.  Exact mixed
moments are

```text
rho=2,          s=8,          a=0,          b=-4.        (26)
```

Both sides of `(16)` equal `64`.  Along `H=U`,

```text
H^2=1
```

pointwise, so `{1,H,H^2}` is dependent and the Pearson residual vanishes.
This is the sharp reason strictness cannot be stated for an arbitrary
positive functional.

### 4.4. A positive exponential base with nonlinear observables

Let `S` have density `e^(-s)` on `(0,infinity)` and put

```text
Z=exp(2pi i exp(-S)).                                    (27)
```

The variable `T=exp(-S)` is uniform on `(0,1)`, so

```text
L(Z^m)=integral_0^1 exp(2pi i m t)dt=0
```

for every nonzero integer `m`.  Even a positive exponential base measure
does not exclude `(21)` when the chosen observables are not polynomials.
The factorial strictness in Section 2 uses both the measure and the
polynomial algebra.

## 5. Three-dimensional harmonic multipoles

Let `W` be a centered real three-space with positive quadratic, cubic,
and quartic moment forms

```text
Q(X)=L(H_X^2),       C(X)=L(H_X^3),       F(X)=L(H_X^4). (28)
```

After a real linear whitening, take

```text
Q=x^2+y^2+z^2.                                          (29)
```

The harmonic projections are

```text
C^o=C-Q Delta C/10,                                     (30)

F^o=F-Q Delta F/14+Q^2 Delta^2 F/280.                   (31)
```

They are respectively harmonic of degrees three and four, and restriction
to the null conic `Q=0` forgets exactly the `Q`-multiples removed in
`(30)--(31)`.

Use the conic parametrization

```text
nu(r,s)=(r^2-s^2, i(r^2+s^2), 2rs).                     (32)
```

A nonzero real harmonic polynomial of degree `d` restricts to a binary
form of degree `2d`.  Its roots are paired by the fixed-point-free real
structure

```text
[r:s] -> [-conj(s):conj(r)].
```

Every pair is the intersection with a unique real projective line.
Therefore the real Maxwell multipole factorization is

```text
C^o=gamma Harm product_(i=1)^3 (u_i dot X),
F^o=delta Harm product_(j=1)^4 (v_j dot X),             (33)
```

with real projective directions `u_i,v_j`; `gamma,delta` are nonzero
when their harmonic is nonzero.  This follows directly by pairing the
conic roots.  The difference between the two sides of either identity in
`(33)` is harmonic and divisible by `Q`, hence zero by uniqueness of the
Fischer decomposition.

For a real vector `u=(u_1,u_2,u_3)`, its conic quadratic is

```text
q_u=(u_1+i u_2)r^2+2u_3rs+(-u_1+i u_2)s^2.             (34)
```

Exact elimination gives

```text
Res(q_u,q_v)=-4 ||u cross v||^2.                        (35)
```

Resultant multiplicativity now yields

```text
Res(C^o o nu,F^o o nu)
 =2^24 gamma^8 delta^6
  product_(i=1)^3 product_(j=1)^4 ||u_i cross v_j||^2.
                                                                  (36)
```

Thus the complete real harmonic bad-locus criterion is

```text
Q=C=F=0 has a nonzero complex solution
iff
C^o=0
or F^o=0
or some cubic multipole line equals a quartic multipole line.     (37)
```

The zero-harmonic branches must be stated separately: multipole lines
are undefined there.  If `C^o=0`, then `C` vanishes on the whole conic
and the binary octic `F|_(Q=0)` has a root.  If `C^o!=0` and `F^o=0`,
any of the six roots of `C|_(Q=0)` is common.  When both are nonzero,
`(35)--(36)` give the last alternative in `(37)`.

With THM-2843's normalization

```text
c_6=(C o nu)/2,              f_8=(F o nu)/3,
```

the same formula is

```text
Res(c_6,f_8)
 =(2^16/3^6) gamma^8 delta^6
  product_(i,j)||u_i cross v_j||^2.                    (38)
```

## 6. Factorial four-slot specialization

For four normalized factorial slots, the real mean-zero space has
dimension three and `(28)` is exactly the first-window moment triple of
THM-2843.  In this factorial setting

```text
C^o!=0.                                                 (39)
```

Indeed, `C^o=0` would mean `Q` divides `C` on the entire three-space.
Restricting to any three-slot mean-zero face would make its positive
binary quadratic divide its cubic, contradicting THM-2824.

Consequently the factorial bad locus reduces exactly to

```text
F^o=0
or a shared cubic--quartic multipole line.              (40)
```

The reduction is global: it does not freeze one coordinate face or choose
one of the three conjugate pairs from THM-2843.  It also does not prove
that `(40)` is empty.

The `3 by 4` multipole incidence product is only a resultant/norm
aggregate.  It supplies neither THM-2655's connected `V_4` Kummer torsor
nor a graph-quartic realization for the Jacobian-conjecture resolvent
lane.  A shared discriminant or resultant is not that missing geometric
sidecar.

## 7. Strictly positive spherical multipole hostiles

Even strict positivity of a real density does not remove either branch of
`(40)` in an abstract moment model.

On the unit sphere with normalized area measure, use the density

```text
1+xyz.                                                   (41)
```

It is strictly positive because `|xyz|<=1/(3sqrt(3))`.  The law is
centered and isotropic.  For coefficient variables again denoted `x,y,z`,
its moment harmonics are

```text
C^o=(2/35)xyz,                    F^o=0.                (42)
```

Thus the vanishing-quartic-harmonic branch is genuinely compatible with
strict density positivity.

The second density

```text
1+xyz+xy(x^2-y^2)                                       (43)
```

is also strictly positive, since

```text
|xyz|+|xy(x^2-y^2)|
 <=1/(3sqrt(3))+1/2<1.
```

It has

```text
C^o=(2/35)xyz,
F^o=(8/315)xy(x^2-y^2).                                (44)
```

The two harmonics share the multipole lines `x=0` and `y=0`.  This is an
exact strictly positive hostile to excluding the shared-line branch by
Gram, Hankel, or density positivity alone.  These sphere laws are not
factorial-slot laws.

## 8. The factorial derivative witness and its missing access

Return to the factorial integral `(18)`.  Integration by parts gives, for
every polynomial `Z` and `r>=0`,

```text
L(Z^r Z')
 =[L(Z^(r+1))-Z(0)^(r+1)]/(r+1).                       (45)
```

At a putative first-window common point with

```text
L(Z)=L(Z^2)=L(Z^3)=L(Z^4)=0,
```

this becomes

```text
L(Z^r Z')=-Z(0)^(r+1)/(r+1),             0<=r<=3.      (46)
```

If `Z(0)!=0`, all four derivative pairings are nonzero and have a rigid
phase progression through the powers of one boundary value.  This is the
cleanest transverse factorial signal exposed by the whitening coordinate.

It is not yet a scalar-moment proof.  The original data provide
`L(Z^m)`, not the multiplier/lowering responses `L(Z^r Z')`.  THM-2815's
finite quotient has the same access boundary: a faithful carrier does not
make an unavailable multiplier observable.  If `Z(0)=0`, even `(46)`
vanishes.  A valid closure therefore needs either derivative access, a
lowering identity that stays inside the scalar moment ideal, or an
independent argument excluding both cases of `(40)`.

## 9. Exact companion and scope

The exact companion verifies:

1. the whitening expansion, cubic/quartic quotients, Pearson identity,
   conic coefficients, and both Bessel determinants;
2. the fifth-root Toeplitz rank, the full two-radius formulas, the six-atom
   equality, and eight exponential-circle Fourier controls;
3. the degree-three and degree-four harmonic projection constants;
4. the symbolic pair identity `(35)`, the full twelve-pair resultant
   constant `(36)`, transverse and coincident-line controls, and the
   separate `F^o=0` common-point branch;
5. the two strictly positive spherical hostile tensors; and
6. `(45)` for three unrelated exact polynomial controls and every
   `0<=r<=3`.

All truth-bearing checks use explicit exceptions, with no Python
`assert`, floating-point sign decision, or scratch dependency.  Reproduce
with

```text
python 04-computation/gmc_whitened_moving_plane_multipole_thm2848.py
python -O 04-computation/gmc_whitened_moving_plane_multipole_thm2848.py
```

Both modes byte-match the stored transcript.

An independent hostile audit rederived every whitening, quotient, Pearson,
Bessel, harmonic-projection, resultant, spherical-tensor, and
integration-by-parts constant; checked all zero-harmonic and strictness
branches; and independently replayed normal, optimized, and stored
evidence with the declared hashes.  It found no remaining defect.

This theorem proves neither unbounded four-slot SFC nor a new GMC case,
and it does not turn the derivative witness into a scalar-moment
observable.
