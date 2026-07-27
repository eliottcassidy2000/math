---
id: THM-2463
title: "Degree-twenty-two B-D plane square-lift closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the open
  first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane C=E=W=0 is empty. For B,D nonzero, the
  ratios p=B/y^2, v=u/y^2, and lambda=D/B^2 give an absolutely
  irreducible quartic moving resultant for every lambda!=0. The
  relation p=B/y^2 is then retained rather than discarded: adjoining
  Y^2=1/p gives a connected double cover. Its p=0 section is the fixed
  squarefree quintic L_5, so five smooth simple zeros force at least
  six branch points and genus at least two. Hence no rational
  trajectory survives. Together with the already closed B and D axes,
  this closes the third of ten support-two planes. It does not close
  the other seven planes, higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-26-degree-twenty-two-BD-square-lift
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
  - THM-2428-degree-twenty-two-B-axis-trigonal-ramification-closure
related:
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
script: 04-computation/jc2_degree22_bd_plane_square_lift_thm2463.py
output: 05-knowledge/results/jc2_degree22_bd_plane_square_lift_thm2463.out
script_sha256: 923454338d35f6b3191b33d2f3c4e20c88083a9325ab46de9d31bbbd919ce08a
output_sha256: ea3c9447c105c71457b033e0f5ae04596214b03d2993ed6023f995953e0b8488
independent_script: 04-computation/jc2_degree22_bd_plane_independent_referee_thm2463.py
independent_output: 05-knowledge/results/jc2_degree22_bd_plane_independent_referee_thm2463.out
independent_script_sha256: a2e4d351bab7c34193fc06b0346512cbd06a368133156929fcd6702fad740db5
independent_output_sha256: ec0b6583f51614f58398bd7cad430f11995f921a0bfb020d0df382ff9582c5ae
hash_basis: working-tree bytes (LF)
---

# THM-2463 -- the degree-twenty-two B-D plane is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The earlier `C,W` and `D,W` planes were closed by classifying every
exceptional member of a positive-genus family. The `B,D` plane initially
appears to require the same treatment: its quartic discriminant has a
moving degree-thirteen factor with a long exceptional divisor.

That is unnecessary. The quotient coordinate

```text
p=B/y^2
```

is not an arbitrary coordinate on the quartic. It remembers a square
class. Restoring that discarded coordinate produces a double cover with
five fixed simple branch points, uniformly in the weighted parameter.
The resulting genus obstruction is both shorter and stronger than a
case-by-case discriminant audit.

The exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
C=E=W=0
    => contradiction.                                           (1)
```

Thus the complete `B,D` support-two coefficient plane is empty.

## 1. Weighted quotient

Use the target-translated coordinates and weights of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `B=0` and `D=0` are respectively the already closed `D` and
`B` axes of THM-2425 and THM-2428. Hence suppose

```text
B,D!=0,

lambda=D/B^2 in C*.                                   (3)
```

First take `y!=0` in `C(x)` and put

```text
p=B/y^2,                v=u/y^2,                zeta=Z/y^3.
                                                               (4)
```

Then

```text
D/y^4=lambda p^2.                                    (5)
```

Dividing the first two fluxes `N_1=N_2=0` of THM-2411 by
`y^5,y^6` gives

```text
f_1
 =2044416 lambda p^2-2981440pv+819896p zeta+24640p
  +3689532v^2-1449459v zeta-101640v+83853zeta+252
 =0,                                                        (6)

f_2
 =-1978994688 lambda p^2v+16355328 lambda p^2
  +1443016960pv^2-71554560pv+65591680p zeta+98560p
  -1190488992v^3+147581280v^2-162339408v zeta-1219680v
  +15944049zeta^2+2236080zeta+672
 =0.                                                        (7)
```

The open first-flux chart is

```text
A=616p-1089v+63!=0.                                  (8)
```

Thus (6) reconstructs `zeta` uniquely from `(p,v)`.

## 2. Exact moving resultant

Define

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (9)
```

Exact elimination gives

```text
Res_zeta(f_1,f_2)=2^4 11^6 R_lambda(p,v),             (10)
```

where

```text
R_lambda
 =-7915978752 lambda(-297lambda+5929v+441)p^4

  +(165944641572864lambda v^2-1371443318784lambda v
    -392920399872lambda+34222590223360v^2
    +3959638538240v-44411530240)p^3

  +(-146683209961728lambda v^3+9698063468544lambda v^2
    +350653121280lambda v-8887076352lambda
    -149234938081152v^3-9500102156160v^2
    +695599766400v-6017413248)p^2

  +(206782580709936v^4+6246495741024v^3
    -1509756494400v^2+34466937120v-193496688)p

  -567L_5(v).                                          (11)
```

The fixed section has

```text
gcd(L_5,L_5')=1.                                      (12)
```

Moreover, the constant and linear coefficients in `p` are coprime as
polynomials in `v`. Hence `R_lambda` has no nonconstant vertical factor.

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `R_lambda` is absolutely
irreducible in `C[p,v]`.

### 3.1 Newton factor shapes

For

```text
lambda!=49/33,
```

the Newton polygon is

```text
Delta(4,1)=conv{(0,0),(4,0),(4,1),(0,5)}.             (13)
```

The exceptional coefficient is exactly

```text
[p^4v^0]R_lambda
 =71243808768 lambda(33lambda-49).                    (14)
```

At `lambda=49/33`, the polygon becomes

```text
conv{(0,0),(3,0),(4,1),(0,5)},                       (15)
```

and its `(3,0)` coefficient is `-627838790656`, so there is no further
drop.

For (13), distributing the primitive edge lengths in a Minkowski
decomposition gives summands `Delta(r,e)` with

```text
0<=r<=4,                    e in {0,1}.               (16)
```

For (15), its edge directions have lengths `(3,1,4,5)`. A summand has
lengths

```text
(a,b,a+b,a+2b),

0<=a<=3,                    b in {0,1}.               (17)
```

After excluding a vertical factor, both polygons therefore have the
same complete factor audit:

- a `1+3` split contains either

  ```text
  p+av+b                                             (18a)
  ```

  or

  ```text
  (v+a)p+(bv^2+cv+d);                               (18b)
  ```

- in a `2+2` split, exactly one factor has triangular polygon
  `Delta(2,0)` and can be normalized as

  ```text
  F=p^2+(av+b)p+cv^2+dv+e.                           (19)
  ```

Any factorization into more pieces contains one of these factors.

### 3.2 The two linear shapes are impossible

Substitute the root of (18a) into (11), collect the six coefficients
in `v`, and form their ideal in

```text
Q[a,b,lambda].
```

Its reduced Groebner basis is

```text
[1].                                                  (20)
```

For (18b), substitute

```text
p=-(bv^2+cv+d)/(v+a),
```

clear the fourth power of the denominator, and collect the ten
coefficients in `v`. Their ideal in

```text
Q[a,b,c,d,lambda]
```

also has reduced basis

```text
[1].                                                  (21)
```

These computations include `lambda=49/33`; no generic-polygon
specialization is being inferred there.

### 3.3 The quadratic shape is impossible

The total-degree-five part of (11) factors exactly as

```text
-38974342 v(56p-99v)^2
  (384lambda p^2-280pv+231v^2).                      (22)
```

Since (19) is monic in `p`, its homogeneous quadratic part must have
one of three forms:

```text
I.   (56p-99v)^2;

II.  384lambda p^2-280pv+231v^2;

III. (56p-99v)(p-hv),
     384lambda h^2-280h+231=0.                       (23)
```

For type I, substituting

```text
a=-99/28,                    c=9801/3136
```

in the eleven coefficient equations obtained by dividing (11) by
(19) gives the unit ideal in `Q[b,d,e,lambda]`.

For type II,

```text
a=-35/(48lambda),            c=77/(128lambda).        (24)
```

The top-next two equations are linear in `(b,d)` with determinant

```text
27102253467051586289664(891lambda-196)^2.             (25)
```

At `lambda=196/891`, the full coefficient ideal is already `[1]`.
Away from that value, the unique solution is

```text
b=-7(78682428lambda^2-34767117lambda+3841600)
      /(23232lambda(891lambda-196)^2),

d=147(2822688lambda^2-1246014lambda+137543)
      /(11264lambda(891lambda-196)^2).                (26)
```

After (26), the numerator ideal of the remaining equations in
`Q[e,lambda]` is `[1]`.

For type III, put

```text
lambda=(280h-231)/(384h^2),

a=-(99/56+h),                 c=(99/56)h.             (27)
```

Here `h!=0`. The two top-next equations depend on `b,d` through
`bh+d`. Their compatibility divisor, after removal of a nonzero
monomial in `h`, is

```text
(20h-33)(448h^2-1320h+1089)^2.                       (28)
```

At `h=33/20`, the full ideal in `Q[b,d,e]` is `[1]`. Adjoining

```text
448h^2-1320h+1089
```

to the full ideal in `Q[b,d,e,h]` again gives `[1]`. Thus all three
quadratic top types are impossible. Together with Sections 3.1--3.2,
this proves uniform absolute irreducibility.

## 4. Restore the square coordinate

Let `C_lambda` be the smooth projective normalization of

```text
R_lambda(p,v)=0.
```

Section 3 makes it connected. At `p=0`, equation (11) gives exactly
the five distinct points

```text
L_5(v)=0.                                             (29)
```

At each such point,

```text
partial_v R_lambda(0,v)=-567L_5'(v)!=0.               (30)
```

Hence the plane curve is smooth there, `p` is a local parameter, and
each point is a simple zero of the rational function `p` on
`C_lambda`.

Now form the quadratic function-field extension

```text
Y^2=1/p.                                              (31)
```

The five odd valuations in (29) show that `1/p` is not a square, so
(31) defines a connected double cover

```text
D_lambda -> C_lambda.                                 (32)
```

All five points (29) ramify. The branch divisor of a double cover has
even degree, so its total branch count `r` satisfies

```text
r>=6.                                                 (33)
```

Riemann--Hurwitz gives

```text
2g(D_lambda)-2
 =2(2g(C_lambda)-2)+r,

g(D_lambda)
 =2g(C_lambda)-1+r/2
 >=2.                                                 (34)
```

This bound is uniform in every `lambda!=0`. It needs no classification
of the degree-thirteen discriminant factor and no assumption that
`C_lambda` itself has positive genus.

## 5. Rational-trajectory and boundary closure

Choose `sqrt(B)` in the algebraically closed constant field. The actual
trajectory relation (4) supplies

```text
Y=y/sqrt(B),

Y^2=y^2/B=1/p.                                        (35)
```

Thus `(p,v,Y)` gives a rational map from `P^1` to `D_lambda`. If it
were nonconstant, it would contradict (34) by Riemann--Hurwitz.
Therefore `p,v,Y`, and hence `y,u`, are constant.

The open coefficient (8) is nonzero, so (6) reconstructs constant
`zeta`; consequently `Z=zeta y^3` is constant. THM-2411 has

```text
Z=T^2,                    T=q^2!=0.                  (36)
```

Thus `T` and `q` are algebraic over the algebraically closed constant
field and are constant. This contradicts the genuine deck, which fixes
the constant field but sends `q` to `-q`.

It remains to justify division by `y`. At `y=0`, the first flux in this
plane is

```text
N_1=1331(616B-1089u)Z.                               (37)
```

The open chart makes the parenthesis nonzero, so `Z=0`, contradicting
`T!=0`. This proves (1).

## 6. Structural boundary

Together with THM-2429 and THM-2437, this closes three of the ten
support-two planes in the inherited degree-twenty-two chart:

```text
(C,W),                  (D,W),                  (B,D).
```

The theorem does not address the other seven planes, higher mixed
strata, branches outside THM-2411, split/even short edges, integral
`2`-adic order raising, `JC(2)`, or `DC(2)`.

The reusable mechanism is a controlled-forgetting correction:

```text
quartic quotient R(p,v)=0
  loses p=B/y^2;

restore the square class Y^2=1/p
  -> five fixed simple branch points
  -> uniform positive genus.                         (38)
```

A difficult discriminant family can therefore become easy after
restoring the coordinate that made the quotient physical.

## 7. Independent hostile audit

The audit first reconstructed (6)--(7) directly from THM-2411's original
fluxes before eliminating `zeta` by the closed formula

```text
Res_z(a z+b,A z^2+B z+C)=A b^2-a B b+a^2 C.
```

This independently recovers the content and all `28` terms of (11). It then
enumerates every balanced allocation of primitive edge lengths, rather than
assuming the factor shapes from a polygon label. For the generic polygon the
allocations are exactly

```text
(a,b,a,a+b),                 0<=a<=4, b in {0,1},
```

and for the exceptional `lambda=49/33` polygon they are exactly

```text
(a,b,a+b,a+2b),              0<=a<=3, b in {0,1}.
```

After the one generic vertical split is removed by the coefficient gcd, each
proper decomposition in both lists either has one of the two linear factors
in Section 3.2 or is the unique `2+2` split containing the triangular factor
(19). Thus the exceptional polygon creates no untested factor shape.

The referee rebuilt all seven exact coefficient ideals from this independently
derived eliminant. A unit ideal is a polynomial identity over `Q`, so it
excludes solutions over every algebraically closed characteristic-zero
constant field, not merely rational parameter values. The denominator ledger
is complete:

- `lambda=0` is outside (3);
- the type-II generic solve has denominator
  `lambda(891lambda-196)^2`, and `lambda=196/891` has its own unit-ideal
  calculation;
- type III has only the denominator locus `h=0`, while
  `384lambda h^2-280h+231=0` has value `231` there; the separate value
  `h=231/280` gives `lambda=0` and is outside (3) (the exact unit-ideal
  computation in fact excludes the broader cleared system as well).

In particular `lambda=49/33` is nonzero, differs from `196/891`, and is
covered both by its exact Newton inventory and by the same factor ideals.

The geometric step was audited on the normalization rather than on the plane
model. At each of the five distinct roots of `L_5`, (30) is nonzero, so `p`
has valuation exactly one. Hence `1/p` has five odd valuations, is not a
square, and gives a connected double cover. Ramification occurs at all five;
the parity of a principal divisor makes the total number of odd valuations
even, hence at least six. Formula (34) then gives genus at least two without
assuming any genus floor for `C_lambda`.

Finally, the physical lift is exactly `Y=y/sqrt(B)`, so no artificial square
class was added. Constancy on the lift makes `y,u` constant; the open wall
reconstructs `zeta`, then `Z=T^2`, `T=q^2`, and `q` are constant, contradicting
the genuine deck. If `y` is the zero rational function, (37), the open wall,
and `T!=0` give the contradiction directly. This also distinguishes harmless
pointwise zeros of a nonzero rational function `y`, for which division in
`C(x)` is valid, from the separately audited identically-zero boundary.

Run the independent path with

```bash
python3 04-computation/jc2_degree22_bd_plane_independent_referee_thm2463.py
python3 -O 04-computation/jc2_degree22_bd_plane_independent_referee_thm2463.py
```

Its normal, optimized, and stored transcripts are byte-identical, and the
independent hashes are recorded in the frontmatter. No factor-shape,
exceptional-parameter, denominator, normalization, branch-parity, boundary,
constant-field, scope, or reproducibility defect remains.

## 8. Exact companion

Run

```text
python 04-computation/jc2_degree22_bd_plane_square_lift_thm2463.py
python -O 04-computation/jc2_degree22_bd_plane_square_lift_thm2463.py
```

The companion:

- reconstructs (6)--(11) by exact elimination;
- verifies the generic and exceptional Newton polygons;
- proves the absence of vertical factors;
- checks both linear factor ideals and all five quadratic subcase
  ideals are the unit ideal;
- verifies the type-II determinant and type-III compatibility divisor;
- checks the fixed squarefree quintic section and the square-lift genus
  floor;
- verifies the `y=0` first-flux boundary;
- uses explicit `require` checks which remain active under `-O`.

The normal, optimized, and stored transcripts are byte-identical.
LF-normalized hashes are recorded in the frontmatter.
