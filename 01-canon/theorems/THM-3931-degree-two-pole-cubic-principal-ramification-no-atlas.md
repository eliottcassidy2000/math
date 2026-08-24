---
id: THM-3931
title: "Degree-two pole cubic has principal quintic ramification and no plane atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The
  THM-3930 degree-two finite-pole cubic has scalar units and class group Z.
  Its two visible quintic addresses u=+/-1 over (A,C)=(3kappa,0) land on
  the same unique length-two ramified point of the cubic fibre. Thus the
  quintic ramification prime is not unibranch and cannot be a boundary
  curve of an affine-plane open by THM-3920. Independently, its vertical
  ramification prime has nonprimitive class 2q, while its
  rational quintic ramification prime E is principal. Explicitly, if theta
  is the endpoint generator, its derivative factors as dF with
  F=3omega-12kappa theta-2C and div(F)=E. Every Keller plane open must
  delete both ramification primes; their classes (2q,0) cannot be a boundary
  basis. More directly, deleting E makes the nonconstant polynomial F a
  nowhere-zero function on the putative A2, contradicting scalar units.
  Thus the first degree-two finite-root-pole survivor is not a JC(2)
  counterexample. Higher pole maps and other coefficient grammars remain
  open.
source: jc_zero_debt_lift / THM-3930 class-and-different audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23). The audit reconstructed
  the factorial chart and its complete three-prime complement, the saturated
  Nagata relation lattice and scalar-unit kernel, the derivative/norm
  factorization div(F)=E, and the collision fibre. Independently,
  the two normalization tangent vectors at u=+/-1 have determinant
  256*kappa !=0, so the two addresses are genuinely distinct branches of E
  through its unique ramified fibre point. This verifies the THM-3920
  obstruction without relying on the class calculation; the principal E and
  nonprimitive P_0 classes then give a second THM-3922 obstruction. The
  assertion-free 53-gate companion LF-normalizes exactly to the frozen raw-LF
  output in normal and optimized mode; raw and semantic hashes and
  documentation checks pass. No repair was required.
depends_on:
  - THM-3930-two-pole-linear-color-aligned-one-place-branch-packet
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3927-unit-ideal-rational-sextic-affine-address-cap-two-place-boundary
script: 04-computation/jc2_degree_two_pole_cubic_class_different_thm3931.py
output: 05-knowledge/results/jc2_degree_two_pole_cubic_class_different_thm3931.out
script_sha256: cefec4c2a18fc3ac3302f34abcb17c54d1ce5f4521e9b99e85a9448f43c03a83
output_sha256: 1f70e1672ac2a11ecfd6718beb26b4586f911491b7945a9d7364c4d5c0672e81
semantic_sha256: 6ceb9cb4e9116d97fc68971a886206a7c3232280648cedf4a63160180fe8b298
hash_basis: raw LF bytes
---

# THM-3931 -- the first degree-two pole escape has a principal ramification debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. Fix `kappa in k`
with

```text
3 kappa^2+1=0.                                             (1)
```

Put `R=k[A,C]`. The THM-3930 binary cubic is

```text
Phi=a U^3+C U^2V+c UV^2+d V^3,                            (2)

a=(A-3kappa)^2,
c=6kappa A+2/3,
d=(A-kappa/3)/2.                                          (3)
```

Let `S` be its Delone--Faddeev algebra with basis `(1,omega,theta)`. Its
multiplication laws are

```text
omega^2=-ac+C omega-a theta,
omega theta=-ad,
theta^2=-Cd+d omega-c theta.                              (4)
```

THM-3930 proves that `S` is a normal finite-flat cubic domain with generic
group `S3`. Its squarefree binary discriminant is

```text
Delta=dG/6,                                                (5)
```

where `G` is an irreducible rational quintic and the two factors are the
vertical line and quintic branch. This theorem first resolves the two-address
collision left visible in THM-3930, then computes the independent divisor
class and different invoices that decide the plane-atlas question.

## 1. A two-boundary factorial chart

The monic equation of `theta` is

```text
P(T)=T^3+cT^2+CdT+ad^2.                                  (6)
```

On `d theta!=0`, equations `(4)` solve for the missing generators:

```text
omega=-ad/theta,
C=-(theta^3+c theta^2+ad^2)/(d theta).                    (7)
```

Conversely `(7)` satisfies `(4)`, so

```text
S_(d theta)=k[A,d^(-1),theta,theta^(-1)].                 (8)
```

This is a UFD, and its complete unit group is

```text
k* d^m theta^n,                          m,n in Z.          (9)
```

The chart complement has exactly three height-one primes:

```text
P_0=(d,omega,theta),
P_C=(d,theta,omega-C),
Q=(A-3kappa,theta,omega-C).                              (10)
```

Indeed, at `d=0`, relation `(1)` gives `a=-64/27`, `c=0`, and the monic
equation of `omega` becomes

```text
W^3-CW^2+acW-a^2d = W^2(W-C).                            (11)
```

Thus `P_0` is the double/ramified prime and `P_C` the simple companion.
The vertical discriminant factor in `(5)` is squarefree, so their
ramification indices are exactly two and one. Therefore

```text
div(d)=2P_0+P_C.                                          (12)
```

For `theta`, the mixed relation in `(4)` gives order one at both vertical
primes. At `P_C`, `omega=C` is a generic unit. At `P_0`, if `v` is its
normalized valuation, then `v(d)=2`; the first and mixed equations in `(4)`
force

```text
v(omega)=v(theta)=1.                                      (13)
```

At `Q`, the element `a=(A-3kappa)^2` has order two while `d` and `omega=C`
are generic units, so `theta` has order two. Hence

```text
div(theta)=P_0+P_C+2Q.                                   (14)
```

Equations `(6)` and `(8)` show that `(10)` lists every chart-boundary prime;
there are no omitted divisor components.

## 2. Class group and units

Nagata localization for `(8)` presents `Cl(S)` on `(P_0,P_C,Q)` with
relation rows

```text
        P_0 P_C Q
d        2   1  0
theta    1   1  2.                                       (15)
```

The minor on the first two columns is one, so the relation lattice is
saturated. Thus

```text
Cl(S)=Z q,
[Q]=q,                 [P_0]=2q,          [P_C]=-4q.       (16)
```

This also proves `S*=k*`. A chart unit `(9)` has divisor coefficients

```text
(2m+n, m+n, 2n).                                         (17)
```

They all vanish only for `m=n=0`.

## 3. The quintic ramification divisor is principal

The relation `(1)` turns the middle coefficient in `(3)` into

```text
c=12kappa d.                                              (18)
```

Differentiate `(6)` at `theta`:

```text
D_theta=P'(theta)=3theta^2+2c theta+Cd.                   (19)
```

Using the last relation in `(4)` and `(18)` gives the exact factorization
inside `S`:

```text
D_theta=dF,
F=3omega-12kappa theta-2C.                               (20)
```

The standard index-form identity and `(5)` give

```text
disc(P)=d^2 Delta,
Norm_(S/R)(D_theta)=-d^2 Delta,
Norm_(S/R)(F)=-Delta/d=-G/6.                              (21)
```

The minus sign is the cubic sign in `Norm(P'(theta))=-disc(P)`; nonzero
base scalars are irrelevant to divisors.

Let `E` be the ramification prime over the irreducible quintic `G=0`.
At the three chart-boundary primes, `(20)` reduces generically to

```text
F mod P_0 = -2C,
F mod P_C = C,
F mod Q   = C.                                            (22)
```

Thus `F` has no vertical or `Q` divisor. Away from `d=0`, `theta` is a
monogenic generator, and `P'(theta)` vanishes at the double root over `G`
but not at its simple companion. Since `G` is squarefree, `(21)` has total
valuation one there. Consequently

```text
div(F)=E,                         [E]=0 in Cl(S).           (23)
```

This identifies the two ramification primes without conflation:

```text
P_0 lies over the vertical line d=0;
E   lies over the rational quintic G=0.                    (24)
```

The simple prime `P_C` is the vertical companion, not ramification.

## 4. The two quintic addresses coalesce on the ramification curve

Write `u` for the THM-3930 normalization coordinate of `G=0`. Its exact
parametrization is

```text
A(u)=u^3+3kappa u^2-u,
C(u)=-(3/2)(u^2-1)(u+kappa)
             (u^2+8kappa u-11/3).                         (25)
```

The two distinct normalization points `u=1` and `u=-1` have the same affine
target address:

```text
A(1)=A(-1)=3kappa,               C(1)=C(-1)=0.             (26)
```

These are genuine finite-root-pole addresses: the repeated-root function
`t=(u+kappa)/(u^2-1)` has a pole at both points because
`(1+kappa)(1-kappa)=4/3`.

The cubic fibre over `(3kappa,0)` is completely explicit. There

```text
a=0,                  c=-16/3,                  d=4kappa/3,
P(T)=T^2(T-16/3),                                          (27)
```

and the Delone--Faddeev relations reduce to

```text
omega^2=omega theta=0,
theta^2=(4kappa/3)omega+(16/3)theta.                       (28)
```

Hence the fibre has exactly two closed points: the unique length-two
ramified point

```text
p_ram=(omega,theta)=(0,0)
```

and the simple etale point

```text
p_et=(omega,theta)=(0,16/3).
```

Indeed `P'(0)=0` and `P'(16/3)=256/9`. Equation `(20)` separates the two:

```text
F(p_ram)=0,                         F(p_et)=-64kappa !=0.   (29)
```

The ramification curve `E=V(F)` has residue degree one over the generic
point of `G=0`. Thus its normalization is finite birational, hence
isomorphic, to the normalization `A1_u` of the quintic. Both `u=1` and
`u=-1` consequently lie over `p_ram`, the only point of `E` in this fibre.
Therefore `E` has two normalization branches at `p_ram`: it is not
unibranch there.

THM-3920 proves that every irreducible boundary curve in a normal affine
surface containing an `A2` open is unibranch at every point. Since every
Keller source open must delete the finite-completion ramification curve
`E`, this already excludes an affine-plane Keller atlas.

Notice the distinction from the raw address cap. Two addresses do not exceed
the cubic bound of three; the obstruction is that both addresses coalesce
on the *same point of the ramification boundary curve*.

## 5. Independent class-and-unit obstruction

Suppose the same finite cubic completion came from a polynomial Keller map
`A2 -> A2_(A,C)`. Normalization-form Zariski Main identifies the source with
an affine-plane open in `X=Spec S`. Because the Keller map is etale, every
ramification divisor of the finite completion belongs to its boundary. In
particular, both `P_0` and `E` must be deleted.

Their classes are

```text
([E],[P_0])=(0,2q) in Cl(S)=Zq.                            (30)
```

They cannot be part of a free boundary basis: `E` is zero, `P_0` is
nonprimitive, and there are already two mandatory primes for a rank-one
class group. This violates THM-3922's necessary boundary-basis theorem.

There is an even shorter unit obstruction. Equation `(23)` says that the
nonconstant polynomial `F` vanishes precisely on `E`. Once `E` is deleted,
its restriction to the putative `A2` has no zero, hence is a unit. But

```text
k[x,y]*=k*,                                                (31)
```

while `F` is nonconstant because its norm is the quintic `-G/6`. This is a
contradiction independent of which additional boundary divisors are removed.

For the natural deletion of exactly the two ramification primes, the class
quotient is

```text
Cl(S)/<[E],[P_0]> = Z/2,                                  (32)
```

and the boundary relation already supplies the nonconstant unit `F`.

Therefore the THM-3930 degree-two pole packet is **not** a planar Jacobian
counterexample. The result does not close degree-two root-pole packets with
different coefficient dependence, higher-degree homogeneous root maps, or
JC(2).

## Reproduction

```bash
python3 04-computation/jc2_degree_two_pole_cubic_class_different_thm3931.py
python3 -O 04-computation/jc2_degree_two_pole_cubic_class_different_thm3931.py
```

Both runs must byte-match
`05-knowledge/results/jc2_degree_two_pole_cubic_class_different_thm3931.out`.
