---
id: THM-2998
title: "Affine quartic reciprocal star-triangle wall, real sign chambers, and the cyclic X_0(14) quotient"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The derivative-square star/complement collision polynomial is an
  primitive irreducible matching-coordinate polynomial and weight-24 affine
  relative invariant, with exact norm formula and a strict
  real 0/2/4-root chamber classification.  Marking one colliding root gives
  two smooth ordered plane cubics.  Their free cyclic order-three quotient,
  not either ordered cubic itself, is the elliptic curve
  y^2=x(x^2+13x+128), hence the standard X_0(14) model.  This supplies no
  projective invariant, canonical quartic origin, Keller/Jelonek owner, or
  LRC carrier.
source: codex-quartic-cross-wall-2026-07-30
depends_on:
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
script: 04-computation/quartic_cross_wall_real_cyclic_quotient_thm2998.py
output: 05-knowledge/results/quartic_cross_wall_real_cyclic_quotient_thm2998.out
script_sha256: b469fce48f19e2bb8c7a92f9945b3ef8ce250efaae0f38ec3a7688b757d45749
output_sha256: a4f45a74b0c8ee27b0d84cab861877e338e3a73562b525a546297ead98216d81
hash_basis: LF-normalized bytes
---

# THM-2998 -- affine reciprocal wall and its cyclic modular quotient

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

[THM-2993](THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding.md)
constructs the four derivative-square star values and their four
discriminant-reciprocal triangle values.  It isolates the final cross wall
`H=0`, but does not determine its real geometry.  [THM-2864](THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product.md)
supplies the ordered edge/orientation viewpoint, while
[THM-2950](THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame.md)
warns that passing to the three-object quotient loses the affine `V_4`
position.

Let `k` have characteristic zero and let

```text
f(T)=T^4+pT^2+qT+r
```

be separable.  Its matching cubic is

```text
S(U)=U^3+aU^2+bU+c,
a=2p,                  b=p^2-4r,                  c=-q^2,
D=Disc(S)=Disc(f),      K=b^3-a^3c.                       (1)
```

Put

```text
H=a^6c^2-2a^4bc^2-2a^3b^3c-26a^3c^3+29a^2b^2c^2
  -2ab^4c-18abc^3+b^6-26b^3c^2+189c^4.                 (2)
```

The theorem has four parts.

1. `H` is primitive and irreducible in the matching coordinates
   `Q[a,b,c]`.  It is a weight-`24` affine relative invariant, but it is not
   a `PGL_2` invariant.  No irreducibility claim is made here for the
   pulled-back polynomial in the depressed-quartic coordinates `(p,q,r)`.
2. Over the reals, `H` is strictly positive on the zero-real-root chamber,
   nonnegative on the two-real-root chamber with four exact oriented affine
   equality representatives, and has an explicit signed four-gap formula on
   the four-real-root chamber.
3. Marking one colliding root produces the two smooth plane cubics

   ```text
   C_+ : xyz=(x-y)(x-z)(y-z),
   C_- : xyz=-(x-y)(x-z)(y-z).                           (3)
   ```

   They are exchanged by an odd flank permutation.  The even flank cycle
   `C_3` acts freely on either curve.
4. The quotient map has degree **exactly three**.  The quotient, not the
   ordered cubic, is

   ```text
   E_0 : y^2=x(x^2+13x+128),                             (4)
   ```

   which is isomorphic to the standard level-`14` modular elliptic curve

   ```text
   X_0(14): Y^2+XY+Y=X^3+4X-6.                           (5)
   ```

   The modular identification in `(5)` and the equivalent equation `(4)` are
   recorded explicitly by Choi--Li, *Quadratic twists of X_0(14)*,
   [J. Number Theory 224 (2021), 142--164, Sections 1--2](https://doi.org/10.1016/j.jnt.2021.01.011).
   The algebra below independently verifies the change of variables; the
   citation supplies the modular interpretation of the target equation.

No rational-point classification of `(3)--(5)` is asserted.

## 2. The exact norm wall and its covariance

Let `x_i` be the four roots of `f`.  THM-2993 identifies the star value at
`x_i` as `f'(x_i)^2` and its complementary triangle value as
`D/f'(x_i)^2`.  Therefore their collision is exactly

```text
f'(x_i)^4=D.
```

Taking the norm through the quartic algebra gives

```text
Res_T(f(T),f'(T)^4-D)
 =prod_i (f'(x_i)^4-D)
 =2^8 D^2 H.                                             (6)
```

This also recovers THM-2993's identity

```text
Res_Z(A_star,Z^2-D)=2^8D^2H.                             (7)
```

The polynomial `(2)` is primitive and irreducible in `Q[a,b,c]`.  A geometric
proof makes the irreducibility informative.  Mark a colliding root, translate
it to zero, and write the other three root differences as `x,y,z`.  Equation
`(6)` becomes one of `(3)`.  The smooth curve `C_+` is irreducible; `C_-` is
its image under an odd relabelling.  Relabelling the marked root does not
change the image in coefficient space.  Hence the closure of the marked
collision locus is irreducible.  The exact point

```text
(a,b,c)=(-14,49,-49)                                    (8)
```

lies on `H=0` and has nonzero `H`-gradient, so the defining polynomial occurs
with multiplicity one.  The companion independently factors `(2)` over `Q`.

If all root differences are multiplied by `lambda`, then

```text
(a,b,c) -> (lambda^2 a,lambda^4 b,lambda^6 c),
H -> lambda^24 H.                                       (9)
```

Translation disappears on depressing the quartic, so `H=0` is a real affine
configuration wall.  It is **not** projective.  Start from

```text
f=T^4-7T^2+7T,              (D,K,H)=(2401,-16807,0).    (10)
```

Apply the determinant-`29` fractional-linear change

```text
U=(T+14)/(1-2T).
```

Its pole `T=1/2` is not a root.  Substituting
`T=(U-14)/(2U+1)`, clearing the fourth power of the denominator, and
monic-normalizing gives the already-depressed integral quartic

```text
U^4-161U^2-581U+1274,                                   (11)
(D,K,H)=(1428170793721,-2238496245503,
         -641584271212871414216880).                     (12)
```

Thus neither a `PGL_2` invariant nor a binary-quartic hyperdeterminant may be
silently substituted for `H`.

## 3. The complete real sign classification

Assume now that `f` is real and separable.

### 3.1 No real roots: strict positivity

Write the roots, after a real translation, as

```text
u+-i beta,                 +-i delta,       beta,delta>0. (13)
```

For `z=u+i beta`, put

```text
A=(u+i(beta-delta))(u+i(beta+delta))
  =u^2-beta^2+delta^2+2i beta u,
f'(z)=2i beta A,                                       (14)
D=16 beta^2 delta^2 |A|^4>0.
```

If `f'(z)^4=D`, then `f'(z)` must be real or purely imaginary.  In the
purely-imaginary case `u=0`, and exact subtraction gives

```text
f'(z)^4-D=16 beta^2(beta-delta)^5(beta+delta)^5.         (15)
```

In the real case `u^2=beta^2-delta^2`, and reduction by that equation gives

```text
f'(z)^4-D=256 beta^6(beta-delta)^3(beta+delta)^3.        (16)
```

Either equality forces `beta=delta`, and then the quartic is inseparable.
Renaming the two conjugate pairs applies the same argument to a prospective
collision in either pair.  Thus no factor in `(6)` vanishes.  The four
factors form two positive norms, so

```text
H>0                         when f has no real root.     (17)
```

### 3.2 Two real roots: four oriented equality representatives

Here `D<0`.  The two real-root factors in `(6)` are positive, and the two
nonreal-root factors are complex conjugates.  Hence

```text
H>=0.                                                       (18)
```

Use the unique orientation-preserving real affine normalization sending the
nonreal pair to `+-i`.  Write the real roots as `y-h,y+h`, `h>0`.  At `i`,

```text
A=(i-y+h)(i-y-h),
f'(i)^4=16A^4,
D=-16h^2A^2 Abar^2.                                    (19)
```

Taking absolute values in `(19)` gives the exact identity

```text
|f'(i)^4|^2-D^2=2^8(1-h^4)|A|^8.                       (20a)
```

Since `A!=0` by separability and `h>0`, equality forces `h=1` before any
phase equation is used.  Then

```text
A=y^2-2-2iy,
A/Abar=+-i,
(y^2-2)^2=4y^2.                                        (20)
```

Consequently `H=0` has exactly the four oriented normalized representatives

```text
y in {-1-sqrt(3), -1+sqrt(3), 1-sqrt(3), 1+sqrt(3)}.    (21)
```

Reflection identifies `y` with `-y`, so these are two classes if the full
orientation-reversing affine group is also divided out.  One exact
coefficient control is

```text
(p,q,r)=(-2-sqrt(3), -2(1+sqrt(3)), (3+4sqrt(3))/4),
D=-4096(7+4sqrt(3)),      K=-512(12+7sqrt(3)),      H=0. (22)
```

### 3.3 Four real roots: four gap walls and the sign rule

Order the roots and let their three consecutive gaps be

```text
alpha,beta,gamma>0,                  L=alpha+beta+gamma.
```

Put `V=sqrt(D)>0` and

```text
rho_i=f'(x_i)^2/V.
```

Direct cancellation gives

```text
rho_1=alpha(alpha+beta)L/[beta gamma(beta+gamma)],
rho_2=alpha beta(beta+gamma)/[gamma(alpha+beta)L],
rho_3=beta gamma(alpha+beta)/[alpha(beta+gamma)L],
rho_4=gamma(beta+gamma)L/[alpha beta(alpha+beta)],
prod_i rho_i=1.                                          (23)
```

Equation `(6)` becomes

```text
H=D^2/2^8 product_i(rho_i^2-1).                          (24)
```

The four wall equations are therefore

```text
alpha(alpha+beta)L       = beta gamma(beta+gamma),
alpha beta(beta+gamma)   = gamma(alpha+beta)L,
beta gamma(alpha+beta)   = alpha(beta+gamma)L,
gamma(beta+gamma)L       = alpha beta(alpha+beta).       (25)
```

Off the walls, `(23)--(24)` give the complete sign test:

```text
H>0  iff exactly two rho_i exceed 1,
H<0  iff exactly one or three rho_i exceed 1.             (26)
```

Both signs and the wall occur among separable four-real quartics:

```text
T^4-7T^2+7T:       (D,K,H)=(2401,-16807,0),
T^4-4T^2+2T+1:     (D,K,H)=(592,-320,-154368),
T^4-10T^2+T:       (D,K,H)=(3973,992000,980121756189).   (27)
```

Thus `D>0`, `K!=0`, or four-real-rootedness alone does not determine the
reciprocal wall or its side.

## 4. The ordered wall is a smooth cubic

Mark one root `x_0` and put

```text
x=x_0-x_1,             y=x_0-x_2,             z=x_0-x_3.
```

Then

```text
f'(x_0)=xyz,
sqrt(D)=+-xyz(x-y)(x-z)(y-z).                           (28)
```

The collision equation `f'(x_0)^4=D` is exactly the union `(3)`.  Each sign
sheet is a smooth plane cubic, hence a genus-one curve.  Swapping two flank
roots exchanges `C_+` and `C_-`; the even three-cycle

```text
(x,y,z) -> (y,z,x)                                      (29)
```

preserves each.

The action in `(29)` is free on the projective curve.  A fixed point would
have the form

```text
[1:lambda:lambda^2],                 lambda^3=1,
```

and substitution in `(3)` has no common root with `lambda^3-1`.  Therefore
the quotient map has degree exactly `3` and is unramified.  After choosing an
origin it is a cyclic `3`-isogeny.  More concretely, the rational orbit

```text
[1:0:0] -> [0:0:1] -> [0:1:0] -> [1:0:0]
```

lies on `C_+`.  Taking `[1:0:0]` as origin, `(29)` is translation by its
rational `3`-torsion successor (the nonzero `j` in `(31)` rules out an
order-three origin-fixing automorphism).  Thus this is literally a rational
cyclic `3`-isogeny, not merely a geometric degree-three map.

Projection of `C_+` from `[1:0:0]` gives

```text
eta^2=tau^4-6tau^3+7tau^2-2tau+1.                       (30)
```

The binary-quartic invariants and discriminant are

```text
(I,J,disc)=(25,-506,-7168),
j(C_+)=-15625/28.                                       (31)
```

The same statements hold for `C_-` by its odd-permutation isomorphism.

## 5. The degree-three quotient is `X_0(14)`

Let `e_1,e_2,e_3` be the elementary symmetric functions of `x,y,z`.  On the
dense chart `e_1!=0`, set

```text
t=e_2/e_1^2,                    s=e_3/e_1^3.
```

Squaring `(3)` and using the cubic discriminant gives the quotient equation

```text
t^2-4t^3-4s+18ts-28s^2=0.                              (32)
```

These symmetric coordinates separate generic `C_3` orbits: on a fixed sign
sheet the three even orderings, and only those orderings, have one symmetric
image.  Thus `(32)` is the degree-three quotient, not a degree-six unordered
map.  Explicitly, over a generic point `(t,s)` of `(32)`, the roots of

```text
U^3-U^2+tU-s
```

recover the normalized triple.  Its discriminant is `s^2`; the six
orderings split into two parity orbits, and the Vandermonde sign in `(3)`
selects exactly the three even orderings on `C_+` (or exactly the three on
`C_-`).  Equivalently, the generic permutation group on the chosen sign
sheet is `A_3=C_3`.

Put

```text
v=28s+2-9t,
x_E=32-112t,                    y_E=-112v.               (33)
```

The inverse on this affine chart is explicit:

```text
t=(32-x_E)/112,             s=(64-9x_E-y_E)/3136.       (33a)
```

On `(32)` one obtains

```text
v^2=(2-7t)(16t^2-11t+2),
y_E^2=x_E(x_E^2+13x_E+128).                              (34)
```

The formulas extend from the dense chart to the smooth projective quotient.
Finally, the change of variables

```text
x_E=4X-4,                  y_E=8Y+4X+4                  (35)
```

with inverse

```text
X=(x_E+4)/4,               Y=(y_E-x_E-8)/8.             (35a)
```

transforms `(34)` into `(5)`; `(35a)` displays its inverse.  The exact
Weierstrass invariants are

```text
E_0:      (c4,c6,Delta)=(-3440,338624,-2^18*7^3),
X_0(14):  (c4,c6,Delta)=(-215,5291,-2^6*7^3).           (36)
```

They differ by the expected powers `2^4,2^6,2^12`.  In particular

```text
j(E_0)=9938375/21952 != -15625/28=j(C_+).               (37)
```

So the ordered cubic is not birational to `X_0(14)`.  The modular curve is
its degree-three cyclic quotient.  The standard modular interpretation of
equation `(5)` is used only to identify the target curve; no modular
parametrization of quartic roots or completeness statement for rational
points is inferred.

## 6. Connection contract and stopping boundaries

The precise bridge is

```text
source:     one separable depressed quartic on H=0, together with a marked
            colliding root and an ordered triple of its flank differences;
map:        projectivize the affine differences, then quotient their even
            cyclic permutation using elementary symmetric coordinates;
target:     the elliptic curve E_0, identified by (35) with X_0(14);
preserved:  the affine reciprocal collision, its sign sheet, and the marked
            star/complement equality;
destroyed:  cyclic flank phase, affine scale, odd ordering, a canonical
            choice of colliding root, PGL2 covariance, and integral owner;
sidecars:   a root/affine-origin section for the quartic and an independent
            graph-owner/Jelonek section for any Keller application;
cheapest hostiles:
            the inversion (10)--(12), the unequal j-invariants (37), and
            the two four-real sign controls in (27).                       (38)
```

Three non-implications are load-bearing.

```text
same elliptic target equation
    does not identify the ordered wall with its quotient;

affine H-covariance
    does not make H a PGL2 invariant or hyperdeterminant;

a marked algebraic star/complement collision
    does not construct a canonical quartic origin, graph sheet, Keller
    map, Jelonek owner, LRC current, or physical carry coordinate.          (39)
```

The appearance of the standard `X_0(14)` equation is therefore an exact
algebraic connection, not a proof of the speculative LRC modular functor.

## 7. Exact companion

Run

```text
python 04-computation/quartic_cross_wall_real_cyclic_quotient_thm2998.py
python -O 04-computation/quartic_cross_wall_real_cyclic_quotient_thm2998.py
```

Both modes LF-byte-match

```text
05-knowledge/results/quartic_cross_wall_real_cyclic_quotient_thm2998.out.
```

The companion checks the generic norm `(6)`, irreducibility and weight,
smoothness at `(8)`, the exact projective hostile, both zero-real obstruction
factors, all four two-real representatives, the four gap walls, three
four-real chamber controls, smoothness of both ordered cubics, freeness and
degree of the cyclic action, the ordered binary-quartic invariants, quotient
equations `(32)--(34)`, the changes and inverses `(33)--(35a)`, both
Weierstrass invariant triples,
and the unequal `j`-invariants.  Every truth gate uses explicit exceptions;
there is no truth-bearing `assert`, floating-point decision, or scratch-file
dependency.

Frozen LF-normalized audited hashes are

```text
script  b469fce48f19e2bb8c7a92f9945b3ef8ce250efaae0f38ec3a7688b757d45749
output  a4f45a74b0c8ee27b0d84cab861877e338e3a73562b525a546297ead98216d81
```

The independent hostile audit rederived the norm and real chambers, checked
smoothness and the free degree-three quotient, and required the explicit
inverse formulas `(33a)` and `(35a)` plus the matching-coordinate
irreducibility scope now recorded above.  QED.
