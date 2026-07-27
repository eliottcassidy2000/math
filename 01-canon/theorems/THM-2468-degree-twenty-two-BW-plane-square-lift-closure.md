---
id: THM-2468
title: "Degree-twenty-two B-W plane square-lift closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane C=D=E=0 is empty. For B,W nonzero, the
  ratios p=B/y^2, v=u/y^2, and lambda=W/B^3 give an absolutely
  irreducible plane quintic for every lambda!=0. Every factorization
  would contain a line or conic; the line coefficient ideal is unit,
  and the three possible conic top types are respectively unit or force
  the forbidden cubic-root boundary h=0. Restoring p=B/y^2 as the
  connected double cover Y^2=1/p is decisive. Its p=0 section is the
  fixed squarefree quintic L_5, so five smooth simple zeros force at
  least six branch places and genus at least two. No rational trajectory
  survives. Together with the axes, this closes the fourth of ten
  support-two planes. It does not close the other six planes, higher
  mixed strata, JC(2), or DC(2).
source: codex-2026-07-26-degree-twenty-two-BW-square-lift
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2428-degree-twenty-two-B-axis-trigonal-ramification-closure
related:
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
script: 04-computation/jc2_degree22_bw_plane_square_lift_thm2468.py
output: 05-knowledge/results/jc2_degree22_bw_plane_square_lift_thm2468.out
script_sha256: 00cf5ce18ffaca3adbe7352d9eb0fea5ab0b9a061b4ba84313e2a329badc45b6
output_sha256: cf08769351e18db72c72c16c1d7395e247bfe523175da2aba1fe97a30f4fc156
hash_basis: working-tree bytes (LF)
---

# THM-2468 -- the degree-twenty-two B-W plane is empty

**PROVED + VERIFIED-EXACT.**

THM-2463 closed the `B,D` plane by restoring a square class that the
weighted quotient had discarded. The same operation closes the `B,W` plane,
but the quotient curve is now a plane quintic rather than a quartic. Its
exceptional-fibre problem is replaced by a complete factor-shape audit.

The exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
C=D=E=0
    => contradiction.                                           (1)
```

Thus the complete `B,W` support-two coefficient plane is empty.

## 1. Weighted quotient and the quintic pencil

Use the target-translated coordinates and weights of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `B=0` and `W=0` are respectively the closed `W` and `B` axes of
THM-2423 and THM-2428. Hence suppose

```text
B,W!=0,                  lambda=W/B^3 in C*.          (3)
```

First take `y!=0` in `C(x)` and put

```text
p=B/y^2,                 v=u/y^2,                 zeta=Z/y^3.
                                                               (4)
```

Then `W/y^6=lambda p^3`. Dividing the first two fluxes `N_1=N_2=0`
of THM-2411 by `y^5,y^6` gives

```text
f_1
 =-2981440pv+819896p zeta+24640p
  +3689532v^2-1449459v zeta-101640v+83853zeta+252
 =0,                                                        (5)

f_2
 =-1319329792lambda p^3
  +1443016960pv^2-71554560pv+65591680p zeta+98560p
  -1190488992v^3+147581280v^2-162339408v zeta-1219680v
  +15944049zeta^2+2236080zeta+672
 =0.                                                        (6)
```

The open first-flux chart is

```text
A=616p-1089v+63!=0.                                  (7)
```

Thus (5) reconstructs `zeta` uniquely. Exact elimination gives

```text
Res_zeta(f_1,f_2)=2^4 11^6 R_lambda(p,v),             (8)
```

where

```text
R_lambda
 =-31289225347072lambda p^5

  +(110629761048576lambda v-6400068820992lambda)p^4

  +(-97788806641152lambda v^2+11314407379968lambda v
    -327276246528lambda+34222590223360v^2
    +3959638538240v-44411530240)p^3

  +(-149234938081152v^3-9500102156160v^2
    +695599766400v-6017413248)p^2

  +(206782580709936v^4+6246495741024v^3
    -1509756494400v^2+34466937120v-193496688)p

  -567L_5(v),                                           (9)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (10)
```

The parameter enters through one exact boundary-square pencil:

```text
R_lambda=R_0-82458112lambda p^3 A^2.                 (11)
```

This explains both the repeated wall direction at infinity and why the
fixed `p=0` section does not move with `lambda`.

## 2. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `R_lambda` is absolutely irreducible
in `C[p,v]`.

Its total degree is five, and its top homogeneous part factors as

```text
-38974342(56p-99v)^2
  (256lambda p^3-280pv^2+231v^3).                    (12)
```

At `v=0`, (12) is

```text
-31289225347072lambda p^5!=0.                         (13)
```

Consequently no factor can have top form divisible by `v` unless the other
top form cancels a factor, which is impossible in a polynomial product.

Any factorization of a degree-five polynomial has a factor of total degree
one or two. A line factor must therefore have the normalized form

```text
p-av-b.                                               (14)
```

Substitution in (9) and collection of the six coefficients of `v` gives an
ideal in `Q[a,b,lambda]` whose reduced Groebner basis is `[1]`. Thus there is
no line factor over any algebraic specialization of `lambda`.

A conic factor cannot have `p`-degree at most one: its degree-two top form
would be divisible by `v`, contrary to (13). It therefore has the normalized
form

```text
F=p^2+(av+b)p+cv^2+dv+e.                              (15)
```

Over the algebraic closure, the top quadratic of (15) selects two of the
five linear factors in (12), counted with multiplicity. There are exactly
three types.

### 2.1 The repeated wall pair

If the top of `F` is `(p-99v/56)^2`, then

```text
a=-99/28,                  c=9801/3136.               (16)
```

Dividing (9) by (15), collecting the nine remainder coefficients, and
substituting (16) gives the unit ideal in `Q[b,d,e,lambda]`.

### 2.2 One wall line and one cubic-root line

Let `h` be a root of the cubic factor in (12). Then

```text
lambda=(280h-231)/(256h^3),

a=-(99/56+h),             c=99h/56.                  (17)
```

Because `lambda!=0`, both `h` and `280h-231` are nonzero. After (17), the
reduced remainder ideal in `Q[b,d,e,h]` has Groebner basis

```text
7744b^2-1584b+81,   88be-9e,   e^2,   d,   h.        (18)
```

Thus divisibility forces `h=0`, outside the chart.

### 2.3 Two cubic-root lines

Omit one cubic root `h`. The product of the other two roots has

```text
a=h,                 c=-231h^2/(280h-231),           (19)
```

with the same expression for `lambda` as in (17). The reduced remainder
ideal now has basis

```text
b^2-e, b d, d^2, b e, d e, e^2,
2bh-d, d h, e h, h^2-11d.                            (20)
```

In particular `h^4=121d^2` lies in the ideal, so every common zero has
`h=0`, again impossible. The three top types exhaust every conic factor,
including repeated cubic roots and coincidences with the wall root. This
proves uniform absolute irreducibility.

Let `C_lambda` be the smooth projective normalization of (9). It is therefore
one connected curve for every `lambda in C*`.

## 3. The retained square class forces genus

The quotient coordinate `p` is not a free plane coordinate. Equation (4)
retains the physical relation

```text
y^2=B/p.                                               (21)
```

Choose `b_0 in C*` with `b_0^2=B` and put `Y=y/b_0`. Every physical
trajectory therefore lifts to

```text
Y^2=1/p                                                (22)
```

over `C_lambda`.

The section `p=0` is exactly

```text
R_lambda(0,v)=-567L_5(v).                              (23)
```

Exact Euclidean division gives

```text
gcd(L_5,L_5')=1.                                      (24)
```

At each of the five roots `alpha` of `L_5`,

```text
partial_v R_lambda(0,alpha)=-567L_5'(alpha)!=0.        (25)
```

Thus all five points are smooth on the plane curve and `p` is a local
uniformizer there. The top form at `p=0` is
`-88239118492602v^5`, so there is no additional `p=0` point at infinity.
Consequently `1/p` has odd valuation at five distinct places of
`C_lambda`. It is not a square, so (22) defines a connected double cover

```text
pi: D_lambda -> C_lambda.                              (26)
```

The branch divisor of a quadratic cover has even degree: equivalently, the
number of odd valuations in the principal divisor of `p` is even. The five
visible simple places therefore force

```text
deg Ram(pi)>=6.                                       (27)
```

Riemann--Hurwitz, without any knowledge of the genus of `C_lambda`, gives

```text
2g(D_lambda)-2
 =2(2g(C_lambda)-2)+deg Ram(pi)
 >=-4+6=2,

g(D_lambda)>=2.                                      (28)
```

This is the decisive gain from restoring the coefficient square class.

## 4. Excluding the rational trajectory

The rational functions `(v,p,Y)` supplied by (4) and (22) give a rational
map from `P^1` to `D_lambda`; it extends to the smooth projective
normalization. By (28), that map is constant. Hence `v,p,Y`, and therefore
`y,u`, are constants. The open first flux (5) reconstructs constant `zeta`,
so `Z,T,q` are constants. This contradicts the genuine deck: it fixes the
constant field but sends the nonzero `q` to `-q`.

It remains to justify division by `y`. At `y=0` with `C=D=E=0`, the original
first flux is

```text
N_1=1331(616B-1089u)Z.                                (29)
```

The open chart makes the parenthesized factor nonzero, so `Z=0`, contrary
to `T!=0`. Together with the two axis theorems, this proves (1).

## 5. Scope and structural lesson

THM-2429 and THM-2437 closed `C,W` and `D,W` by classifying positive-genus
families; THM-2463 and this theorem instead restore coefficient square
classes after proving the quotient curve connected. The fixed section

```text
p=0  ->  L_5(v)=0                                     (30)
```

is now a reusable boundary divisor: five transverse contacts plus branch
parity force genus on every connected quadratic lift. This closes the fourth
of ten support-two planes. The six planes `B,C`, `B,E`, `C,D`, `C,E`, `D,E`,
and `E,W`, higher mixed strata, branches outside the inherited reduction,
split/even short edges, and integral order raising remain open. Nothing here
proves `JC(2)` or `DC(2)`.

## 6. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_bw_plane_square_lift_thm2468.py
python3 -O 04-computation/jc2_degree22_bw_plane_square_lift_thm2468.py
```

The companion reconstructs (5)--(11), the 24-term quintic and its top form,
the unit line ideal, all three conic top charts and their exact Groebner
consequences, the fixed squarefree quintic section, the branch/genus floor,
and the `y=0` boundary. Normal and optimized transcripts byte-match the stored
output. All truth-bearing checks use explicit exceptions and remain active
under optimized Python. **QED.**
