---
id: THM-2846
title: "Arbitrary positive-cone moment-three transverse boundary"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; INDEPENDENT
  HOSTILE AUDIT REQUESTED.  Two positive adjacent-difference cones, both
  divisible by s, span an exact nonzero factorial moment-three null line.
  The common zero is unique and transverse in an explicit rational
  rectangle, while the fourth factorial moment has a fixed nonzero
  imaginary sign.  If promoted, this proves that Gaussian moments through
  six do not detect arbitrary positive cones, although moment eight detects
  the displayed witness.
source: root/audit-2809-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
script: 04-computation/gmc_positive_cone_moment3_transverse_boundary_thm2846.py
output: 05-knowledge/results/gmc_positive_cone_moment3_transverse_boundary_thm2846.out
script_sha256: c1e4d4432c8fe2e988bb0c877ca0d451ff1426d3e8c0ba0f103820dc21eb7526
output_sha256: b1193b34c28e1f49682c6de2e4b68c1b56ea14226df731afa37bcdb7f287b2ea
hash_basis: LF-normalized bytes
---

# THM-2846 -- arbitrary positive-cone moment-three transverse boundary

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; INDEPENDENT
HOSTILE AUDIT REQUESTED.**

This is the sharp failure boundary for extending THM-2830 from its
separated/transport-ordered positive cones to two arbitrary positive
adjacent-difference cones.

## 1. Exact algebraic parameters

Put

```text
L(s^n)=n!,                 f_n=s^n/n!,
d_i=f_(i+1)-f_i.
```

Let `y` be the unique positive root of

```text
p(Y)=
 933319220872000 Y^15
+5884528273360000 Y^14
+17194209549700000 Y^13
+30835666750660000 Y^12
+37890363245864000 Y^11
+33723003604779000 Y^10
+22402533775665500 Y^9
+11275117818448500 Y^8
+4315433784681000 Y^7
+1247681141032500 Y^6
+267329168334540 Y^5
+40861458466200 Y^4
+4132565626725 Y^3
+231008559525 Y^2
+2764978335 Y
-258946389.
```

All nonconstant coefficients are positive and the constant coefficient is
negative.  Descartes' rule and the intermediate value theorem therefore
give exactly one positive root.  Exact rational evaluation brackets it by

```text
2341/100000 < y < 2342/100000,

y=0.02341986411887604324463235332262655...
```

Define

```text
A(Y)=2(
 46665961043600 Y^10
+517246931163000 Y^9
+1734100192751000 Y^8
+2973460269848000 Y^7
+3076524208583550 Y^6
+2056143278933070 Y^5
+909053232722550 Y^4
+264045781471725 Y^3
+48371683260150 Y^2
+5052446879130 Y
+228165400071),

N(Y)=
 485253235285200 Y^10
+2203659617816000 Y^9
+4482965815747000 Y^8
+5358371582968500 Y^7
+4154724741540300 Y^6
+2179627177491390 Y^5
+783069490441350 Y^4
+190343609071200 Y^3
+29999951982225 Y^2
+2773709682480 Y
+114503538087,

x=N(y)/A(y)
 =0.26362138003485122740864992126146544....
```

Both `A` and `N` have positive coefficients, so `x>0`.

### A second, topological-degree proof

The large algebraic formulas above give explicit coordinates, but existence
has a much smaller robust certificate.  Let

```text
F=I1,                       G=6I2-7I1
```

and take the rational rectangle

```text
2636/10000 <= x <= 2637/10000,
23418/1000000 <= y <= 23421/1000000.                  (R)
```

Exact Bernstein-basis conversion on the four faces gives

```text
F(left,y)>0,           F(right,y)<0,
G(x,bottom)>0,         G(x,top)<0.                    (R1)
```

On the whole rectangle the same exact certificate gives

```text
F_x<0,       F_y>0,       G_x<0,       G_y<0,
J=F_x G_y-F_y G_x>0.                                  (R2)
```

For each fixed `y`, `(R1)--(R2)` make `F(x,y)=0` have one
root `x=x(y)`.  Along that curve,

```text
d/dy G(x(y),y)=J/F_x<0.                               (R3)
```

The bottom/top signs in `(R1)` therefore force one and only one common
zero of `F,G`, hence of `I1,I2`, inside `(R)`.  This is an independent
rectangle-degree proof of the hostile.  It also shows that the common zero
is transverse and persists under sufficiently small perturbations of the
Pascal tensors; it is not a numerical tangency.

## 2. The two positive cones

Take

```text
U=d_1+x d_3,                 V=d_2+y d_3.              (1)
```

These are nonzero, linearly independent, positive adjacent-difference
cones.  Both are divisible by `s`.

Write

```text
g11=L(U^2),       g12=L(UV),       g22=L(V^2),

t111=L(U^3),      t112=L(U^2V),
t122=L(UV^2),     t222=L(V^3).                         (2)
```

The THM-2824 division-free cubic remainders are

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,

I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.             (3)
```

Direct exact elimination gives

```text
I1(x,y)=I2(x,y)=0.                                     (4)
```

This assertion is not based on decimal approximation.  The exact companion
substitutes `x=N(y)/A(y)`, clears denominators, and obtains zero remainders
modulo `p(y)` for both polynomials.  It also reconstructs `p` as the unique
positive-root factor of the exact resultant.

## 3. Exact moment-three failure

For

```text
H=U+zV,

Q(z)=L(H^2)=g11+2g12 z+g22 z^2,

C(z)=L(H^3)=t111+3t112 z+3t122 z^2+t222 z^3,           (5)
```

the Gram form `Q` is positive definite over the reals because `U,V` are
independent real polynomials and `L(F^2)>0` for every nonzero real `F`.
Thus either root `z` of `Q` is nonreal.

Equations `(3)--(4)` are exactly the binary divisibility criterion

```text
C(z)=
 (t111/g11+(t222/g22)z) Q(z).                          (6)
```

Choose either root of `Q`.  Then `H` is nonzero and

```text
L(H)=L(H^2)=L(H^3)=0.                                  (7)
```

The first equality holds because every `d_i` has zero factorial mean.
Consequently arbitrary positive adjacent cones are **not** uniformly
detected by factorial moments one through three.

## 4. Gaussian consequence and exact scope

Since both directions in `(1)` are divisible by `s`, write `H=sh` and set

```text
P(Z,W)=W+Z h(ZW),                  W=conj(Z).            (8)
```

More explicitly, choose the upper-half-plane quadratic root

```text
zeta=(-g12+i sqrt(g11 g22-g12^2))/g22
    =-0.9033865066991233...+0.1804736111316870... i.
```

Then the finite hostile is the degree-seven polynomial

```text
P(Z,W)
 =W-Z
  +(1-zeta) Z^2 W/2
  +(-x+zeta(1-y)) Z^3 W^2/6
  +(x+zeta y) Z^4 W^3/24.                            (8a)
```

All parameters in `(8a)` are algebraic by the exact construction in
Section 1; the decimals only identify the chosen branches.

Charge balance gives

```text
E[P^(2k)]=binom(2k,k)L(H^k),       k=1,2,3,
E[P^(2k-1)]=0.                                            (9)
```

Therefore this explicit nonzero two-charge polynomial satisfies

```text
E[P^m]=0,                         1<=m<=6.              (10)
```

This is a lower-bound/no-go theorem for finite-moment approaches: the
degree-six conclusion of THM-2830 cannot be extended to arbitrary positive
adjacent cones.  It is **not** a counterexample to GMC2.

There is also an exact positive control at the next rung.  Reduce
`L((U+zV)^4)` modulo the quadratic `Q(z)`.  After multiplying by the
positive denominator `g22^3`, write its linear remainder as

```text
r0(x,y)+r1(x,y)z.
```

Exact Bernstein conversion proves `r1>0` throughout the entire rational
rectangle `(R)`.  At the upper-half-plane root of `Q`,

```text
Im L(H^4)=r1 Im(z)/g22^3>0.
```

The conjugate root gives the conjugate nonzero value.  Hence neither root
of `Q` kills the fourth factorial moment:

```text
L(H^4)!=0.                                             (11)
```

Thus the corresponding Gaussian eighth moment is nonzero.  The witness
separates the exact thresholds six and eight; it does not evade all
moments.

## 5. Equality geometry

The witness lies outside THM-2830's hypotheses.  Its coefficient profiles
are

```text
lambda=(0,1,0,x,...),
nu_i=mu_(i+1)=(0,1,y,...).
```

Their shifted `2 by 2` minors change sign, so neither monotone-likelihood
ratio order applies.  The two opposed transport contributions cancel
exactly on the algebraic hypersurface `(4)`.  This identifies the proper
next object: not a stronger one-sided sign estimate, but the higher-moment
transversality of this cancellation locus.

### 5.1. Concrete THM-2843 moving plane

Because `d_1,d_2,d_3` span the four normalized monomial slots
`f_1,f_2,f_3,f_4`, the plane

```text
E=span_R{U,V}
```

is an actual four-slot factorial plane in the sense of THM-2843.  Equations
`(4),(6)` say that `Q|_E` divides `C|_E`, while `(11)` says that
`Q|_E` does not divide the quartic restriction.  Thus this is a concrete
factorial realization of THM-2843's moving-plane geometry, not merely the
abstract non-factorial facewise hostile in that theorem.  The quartic
transverse sign is precisely the next complete-intersection exit.

## 6. Exact evidence and audit boundary

The exact companion independently:

1. reconstructs the quadratic and cubic adjacent-difference tensors;
2. verifies the resultant factor and the unique positive algebraic root;
3. substitutes `x=N(y)/A(y)` and checks both invariants and the cubic
   remainder modulo `p`;
4. certifies every face, derivative, and Jacobian sign in `(R1)--(R2)` by
   exact Bernstein coefficients; and
5. certifies positivity of the fourth-moment remainder coefficient `r1`
   on the same rectangle.

All truth-bearing gates are explicit and survive optimized execution.
Normal, optimized, and stored transcripts agree.  The theorem remains a
reserved provisional candidate until an independent hostile audit checks
the factorial/Gaussian typing, the Bernstein conversion, and the
fourth-remainder normalization.
