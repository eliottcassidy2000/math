---
id: THM-3933
title: "Centered finite-at-infinity degree-three root maps collapse to a non-unibranch octic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  centered trace-zero, finite-at-infinity degree-three repeated-root
  stratum of the one-place linear-color binary-cubic grammar, local trace and
  the finite
  Riemann--Hurwitz budget force either a shared pole address or one triple
  pole above a simple critical point. The shared-address branch is
  non-unibranch whenever it is genuine transposition ramification. The
  collision-free case has the normal form A=u^3+u^2 and
  t=(3 lambda u^2+3u-2)/u^3. Polynomiality of the color forces lambda=3.
  Its discriminant is a squared index line times an irreducible rational
  octic. Dividing the Delone--Faddeev basis element theta by the line gives
  the normal maximal order explicitly; it is globally monogenic, while the
  octic remains genuine ramification and has two exact two-address fibres.
  Hence the survivor violates THM-3801's required nonmonogenic-completion
  gate and the unibranch boundary gate. This closes the stated centered
  finite-at-infinity degree-three stratum, not the t(infinity)=infinity
  gauge, arbitrary root changes, higher root degree, or JC(2).
source: jc_zero_debt_lift / post-THM-3931 degree-three finite-root-pole stratum, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (two independent reads, 2026-08-23). The
  auditors reconstructed the completed-local trace table, finite
  Riemann--Hurwitz budget, finite-flat shared-address bridge, and the trace
  argument forcing every finite value t(infinity) to zero. They checked the
  primitive incidence row, the lambda=0 seam and
  unique lambda=3 color row, the line-versus-octic discriminant split, the
  octic normalization and both collision fibres. A separate hostile audit
  checked every Delone--Faddeev multiplication identity, the index ideal
  (27A-4), discriminant -H/16, and the R1+S2 maximal-order argument. The
  45-gate assertion-free companion LF-normalizes exactly to the frozen
  raw-LF output in normal and optimized mode.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
related:
  - THM-3927-unit-ideal-rational-sextic-affine-address-cap-two-place-boundary
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3930-two-pole-linear-color-aligned-one-place-branch-packet
  - THM-3931-degree-two-pole-cubic-principal-ramification-no-atlas
  - THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification
  - THM-3936-centered-degree-three-infinite-root-value-nonentry
script: 04-computation/jc2_centered_degree_three_root_map_octic_thm3933.py
output: 05-knowledge/results/jc2_centered_degree_three_root_map_octic_thm3933.out
script_sha256: 2a02f99465a8badcd32330f8e81ef46e7590793d85ce55957e17252f58e73604
output_sha256: 8753cb79c9214b4df581dcc1c79734be0b48bd1231e155b00c6b7ffe38f9d89b
semantic_sha256: 9c182abe7a8ac734b172dc6b64bcb4b4558cebf131d3be75a9e9a5ef750cd43a
hash_basis: raw LF bytes
---

# THM-3933 -- a centered finite-at-infinity degree-three pole divisor has nowhere left to go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. Consider the
linear-color binary cubic

```text
Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3.                       (1)
```

Assume that one irreducible one-place discriminant component has
normalization `A1_u`, that

```text
A=A(u), C=C(u) are polynomials,       deg A=3,              (2)
```

and that its repeated root `t=U/V in k(u)` satisfies

```text
deg(t:P1_u -> P1_t)=3,          t(infinity) in k.           (3)
```

We work in the centered root gauge: the primitive incidence equation is

```text
a(A)t^3-c(A)t-2d(A)=0,                                  (4)
```

so `Tr_(k(u)/k(A))(t)=0`. We also impose the coefficient-ideal gate

```text
(a,C,c,d)=k[A,C].                                         (5)
```

The normalization hypothesis means `k(A,C)=k(u)`. Thus the generic repeated
root generates the cubic extension over `k(A)` and `(4)` really is its
minimal polynomial after division by coefficient content.

The theorem closes precisely this centered, finite-at-infinity stratum. The
centering is the intrinsic missing-square term in `(4)`, not a root
translation. It does not assert that every degree-three homogeneous root map
can be moved into `(1)-(4)` while preserving the linear-color direction;
in particular, `t(infinity)=infinity` remains outside the theorem.

## 1. Local trace and the pole-partition trichotomy

Let `x` be a finite pole of `t`, let its pole order be `m`, and let `e` be
the ramification index of the cubic polynomial map

```text
pi_A:P1_u -> P1_A
```

at `x`. If no other pole lies over `A(x)`, trace zero forces the local trace
at `A(x)` to be regular. In completed tame coordinates one may write

```text
z=s^e.
```

For a Laurent series `sum b_j s^j`, the local trace retains exactly the
terms with `e|j`, multiplying them by `e`. In particular:

```text
m=1: an isolated pole requires e>=2;
m=2: an isolated pole requires e=3;
m=3: an isolated pole can only have e=2,
     and trace zero also removes its s^(-2) coefficient.    (6)
```

For example, the leading `s^(-2)` term makes `m=2,e=2`
impossible, while the leading `s^(-3)` term makes `m=3,e=3`
impossible. No cancellation with regular points over the same target value
can remove a negative trace term.

The pole divisor of `(3)` has total degree three, so its support partition is

```text
1+1+1,                   2+1,                   or 3.       (7)
```

A cubic polynomial map has total finite ramification budget two: the point
at infinity already contributes two of the four Riemann--Hurwitz units. If
all pole supports have distinct `A`-values, `(6)` makes the first partition
cost at least `1+1+1=3` units and the second at least `2+1=3` units. Both
are impossible. The only collision-free possibility is therefore

```text
(t)_infinity=3x,          e_x(pi_A)=2.                     (8)
```

Thus every other pole partition contains distinct pole supports `x!=y`
with

```text
A(x)=A(y).                                                   (9)
```

## 2. Why a shared pole address is a boundary obstruction

At a pole of `t`, the limiting repeated binary root is `[U:V]=[1:0]`.
Dehomogenizing `(1)` in `U=1` gives

```text
a+Cw+cw^2+dw^3,                    w=V/U.                  (10)
```

For `w=0` to be a multiple root, both its constant and linear coefficients
vanish. Hence every finite pole satisfies

```text
a(A(x))=0,                  C(x)=0.                         (11)
```

Consequently `(9)` gives two distinct normalization addresses over the same
target point `(A(x),0)`.

Here is the finite-completion bridge. The generic repeated root cannot be
triple. Indeed triple-root equations would give

```text
c=3at^2,                    C=-3at,
```

so the derivative `3aT^2-c` of the incidence polynomial `(4)` would vanish
at `T=t`. That contradicts separability of the minimal polynomial of `t`
over `k(A)`. Thus the generic root is exactly double. The universal
discriminant hypersurface is smooth there. Variation in the `C` direction
changes the cubic value at the repeated root by `t^2`; since generic `t` is
finite and nonzero, this direction is transverse. The pulled-back component
is therefore reduced, so its odd discriminant valuation cannot disappear
under an order-to-maximal-order index correction. It is genuine
transposition ramification in the maximal cubic completion.

Let `B` be the maximal normalization and `E` its residue-degree-one
ramification prime over the component. Since `B` is a normal finite module
over the regular target surface, it is Cohen--Macaulay of dimension two and
hence finite flat of rank three. Every closed point of `E` is non-etale, so
its local geometric fibre has length at least two. A length-three fibre
cannot contain two distinct such points. The two addresses in `(9)-(11)`
therefore land on the same point of `E`.

The normalization of `E` is finite birational to the normalization
`A1_u` of the branch component, hence isomorphic to it. The two addresses
are two branches of `E` at one point. THM-3920 says an irreducible boundary
curve of an affine-plane open in a normal affine surface is unibranch.
Therefore every shared-pole row is already incompatible with a Keller
`A2` open.

## 3. Normal form for the sole collision-free pole divisor

It remains to analyze `(8)`. Translate the pole and its `A`-value to zero.
Because the ramification is simple, the polynomial cubic has the form

```text
A=alpha u^3+beta u^2,             alpha beta !=0.
```

Scaling `u` and `A` gives the unique normal form

```text
A=u^3+u^2.                                                 (12)
```

A rational function finite at infinity and having only a triple pole at zero
is

```text
t=h+(b_2u^2+b_1u+b_0)/u^3,                 h in k.
```

Exact trace in the extension `u^3+u^2=A` is

```text
Tr(t)=3h+(3b_0+2b_1)/A.                                   (13)
```

Since `A` is nonconstant, trace zero separately forces

```text
h=0,                         3b_0+2b_1=0.
```

Thus the centered equation forces its own zero value at normalization
infinity; no root translation has been made. Triple-pole nondegeneracy gives
`b_0!=0`, so a root scaling puts every possibility in the one-parameter form

```text
t=(3lambda u^2+3u-2)/u^3.                                 (14)
```

The primitive resultant incidence is exactly `(4)` with

```text
a=-A^3,
c=-9A^2lambda^2-27A^2lambda+12Alambda+15A-4,
d=-(27A^2lambda^3-36Alambda^2+27Alambda
    +27A+12lambda-20)/2.                                  (15)
```

The triple `(a,c,d)` is primitive because `a=-A^3` and `c(0)=-4`. Any
other polynomial incidence defining the same minimal polynomial is
`h(A)(a,c,d)`. Its coefficient ideal together with `C` is `(C,h)` after
using primitivity, so `(5)` forces `h in k*`. We absorb that scalar.

## 4. Polynomial color singles out lambda=3

The derivative repeated-root equation is

```text
3at^2+2Ct+c=0.                                             (16)
```

Substituting `(12),(14),(15)` writes `C(u)=N_C/(2N)` where

```text
N=3lambda u^2+3u-2.                                       (17)
```

and `N_C` is the explicit polynomial frozen in the companion. Their exact
resultant is

```text
Res_u(N,N_C)=-124416 lambda^2(lambda-3)^3.                 (18)
```

If `lambda!=0`, polynomiality requires `N|N_C`, hence `(18)` forces
`lambda=3`. At `lambda=0`, the denominator is `3u-2` while

```text
N_C(2/3)=-512/243 !=0.                                    (19)
```

Conversely direct division at `lambda=3` gives the polynomial

```text
C(u)=u^3(3u+2)(9u^4+30u^3+24u^2-4)/2.                    (20)
```

Thus `(20)` is the unique polynomial-color row. Its binary cubic is

```text
a=-A^3,
c=-(6A-1)(27A-4),
d=-(27A-4)^2/2.                                           (21)
```

It still has unit coefficient ideal: `c(0)=-4`.

## 5. The squared line is index debt; the octic is genuine branch

Let `L=27A-4`. The exact discriminant is

```text
Disc(Phi)=-L^2 H/4,                                       (22)

H=19683A^8+87480A^7-60048A^6+5832A^5C+14688A^5
  -1836A^4C-1584A^4+144A^3C+64A^3
  -144A^2C^2+48AC^2-8C^3-4C^2.                           (23)
```

The cubic is generically irreducible. Indeed it is linear in `C`; any
factor independent of `C` would divide both `T^2` and
`aT^3+cT+d`, whose gcd is one because `d!=0`. Its nonsquare discriminant
then gives generic group `S3` once `H` is proved irreducible below.

The two discriminant factors play different roles. At the generic point of
`L=0`, `a` is a unit, `c=d=0`, and the monogenic root polynomial reduces to

```text
T^2(aT+C).                                                  (24)
```

Thus the displayed order has two distinct primes over `L`. Its discriminant
valuation is two, since

```text
H(4/27,C)=-4C^2(162C+1)/81 !=0                            (25)
```

generically. If the maximal order still had discriminant exponent two,
tame cubic ramification would be total of index three and would have only
one prime. Integral closure cannot merge the two primes in `(24)`. The
order/maximal-order discriminant formula changes valuation by twice the
local index, so the only remaining exponent is zero. Therefore `L^2` is
exactly index debt and `L` is not a ramification divisor of the maximal
completion.

At the generic point of `H=0`, the order discriminant has valuation one.
No nonzero even index correction is possible. Hence the order is already
maximal there and `H` is genuine tame `(2,1)` ramification with a unique
residue-degree-one ramification prime `E`.

## 6. The maximal completion is explicitly monogenic

The line debt can be removed globally, not just valuation by valuation. Let
`(1,omega,theta)` be the standard Delone--Faddeev basis of the cubic order
defined by `(21)`. Its multiplication laws are

```text
omega^2=-ac+C omega-a theta,
omega theta=-ad,
theta^2=-Cd+d omega-c theta.
```

Put

```text
e=theta/L,                         L=27A-4.
```

Substitution of `(21)` into the last two multiplication laws gives

```text
e^2=C/2-omega/2+(6A-1)e,          omega e=-A^3L/2.
```

In particular

```text
omega=C+2(6A-1)e-2e^2,            theta=Le,
```

so the original order `S=k[A,C,omega,theta]` lies in `T=k[A,C,e]`. The
standard theta equation, divided by `L^3`, is the monic equation

```text
P(e)=e^3-(6A-1)e^2-(C/2)e-A^3L/4=0.                     (25a)
```

Relative to the basis `(1,e,e^2)`, the columns of `(1,omega,theta)` form

```text
       [1       C       0]
M =    [0  2(6A-1)     L],             det(M)=2L.
       [0      -2       0]
```

Since `2` is a unit, the index ideal is exactly `(L)`. Direct discriminant
calculation gives

```text
Disc(P)=-H/16,
(2L)^2 Disc(P)=-L^2H/4=Disc(S).                           (25b)
```

The nonzero rank-three transition determinant also shows that `T` is a
domain with the same fraction field as `S`. This overorder is already the
full integral closure. Indeed `T` is a monic hypersurface, finite free of
rank three over the regular ring `R=k[A,C]`; hence it is Cohen--Macaulay and
satisfies `S2`. At every
height-one base prime away from `(H)`, `(25b)` is a unit and `T` is finite
etale. At `(H)`, its discriminant has valuation one. A proper further
normalization index of valuation `j>=1` would change that valuation to
`1-2j<0`, impossible for an integral order. Thus `T` is normal in
codimension one. By `R1+S2`, it is normal everywhere, and being a normal
finite `R`-order in the same cubic field, it is the maximal completion.

Therefore the maximal cubic completion is globally monogenic. THM-3801
excludes it from a polynomial Keller map independently of the boundary
argument below. The squared line was exactly order debt; removing it does
not remove the genuine octic branch.

## 7. The octic normalization and its two collision fibres

Substitution of `(12),(20)` gives `H=0`, and exact elimination gives

```text
Res_u(A-u^3-u^2, C-C(u))=-H/8.                             (26)
```

The parametrization is birational. We have `[k(u):k(A)]=3`. The repeated
root `(14)` at `lambda=3` cannot belong to `k(A)`, since its pole order three
at `u=0` is not divisible by the ramification index two of `A` there.
The derivative equation also gives

```text
C=-(3at^2+c)/(2t),
```

so the unique generic double root belongs to `k(A,C)` and conversely `C`
belongs to `k(A,t)`. Therefore

```text
k(A,C)=k(A,t)=k(u).                                        (27)
```

Therefore `H` is irreducible and `(12),(20)` is its normalization. Since
`u` is integral over `k[A]`, that normalization is the affine line; because
`deg A=3` and `deg C=8`, its projective closure has one normalization place
at infinity.

There are nevertheless two finite two-address fibres. For distinct `u,v`
put

```text
s=u+v,                  p=uv.
```

Dividing `A(u)-A(v)` by `u-v` gives

```text
p=s^2+s.                                                     (28)
```

After this substitution, dividing `C(u)-C(v)` by `u-v` gives

```text
s(3s^2+6s+2)^3/2.                                         (29)
```

The solution `s=0` is the diagonal cusp `u=v=0`. Each of the two roots of

```text
3s^2+6s+2=0                                                (30)
```

instead gives two distinct parameters: the pair discriminant is
`-3s^2-4s`, whose resultant with `(30)` is `-12`. These collision targets
do not lie on the index line; after `(28)` their `L`-value has resultant
`-54` with `(30)`. They are two distinct target fibres as well: modulo
`(30)`, `(s+1)^2=1/3` and hence their common `A`-value is `-s/3`, different
for the two distinct roots `s`.

Thus the normalization of the genuine ramification curve `E` has two
distinct points over each of two affine target points. As in Section 2, a
rank-three fibre contains at most one non-etale point, so each pair
coalesces on `E`. The curve `E` is not unibranch. THM-3920 excludes it as a
boundary curve of an affine-plane open, and the unique polynomial-color
survivor is not a planar Jacobian counterexample. This obstruction is
independent of Section 6's monogenic-maximal-order obstruction.

The conclusion is exactly the centered trace-zero, finite-at-infinity,
degree-three linear-color stratum `(1)-(5)`. The
`t(infinity)=infinity` gauge, arbitrary root changes, degree at least four,
other color directions, and JC(2) remain **OPEN**.

## Reproduction

```bash
python3 04-computation/jc2_centered_degree_three_root_map_octic_thm3933.py
python3 -O 04-computation/jc2_centered_degree_three_root_map_octic_thm3933.py
```

After platform newlines are normalized to LF, both runs must byte-match the
frozen raw-LF output
`05-knowledge/results/jc2_centered_degree_three_root_map_octic_thm3933.out`.
