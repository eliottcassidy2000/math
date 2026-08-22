---
id: THM-3639
title: "Russell-cylinder all-retained-cells universal second-stable debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the fixed minimal
  THM-3635 ordinary triple, every zero-stable/first-stable cell with retained
  determinant J_0=1 has the same second-stable retained cokernel -2072/81,
  whether its zero-stable tangent rank is one or two.  More generally the
  cokernel differs from this constant by an exact endpoint derivative of
  J_0.  The result is independent of J_1, all higher restriction lifts, and
  all target representatives.  It excludes only this fixed quadratic-fold
  compiler and makes no global Keller or JC(2) claim.
source: root / stable_debt_universal_formula all-cell continuation, 2026-08-21
audit: >
  PASS -- two independent hostile audits rederived the unrestricted retained
  residue, order-three representative invariance, exhaustive SL_2 rank split,
  division-free endpoint-curvature identity, and byte-identical normal,
  optimized, and stored 46-gate transcripts.
depends_on:
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
  - THM-3637-russell-cylinder-actual-rank-witness-second-stable-debt
related:
  - THM-3638-russell-cylinder-tangent-rank-one-universal-second-stable-debt
script: 04-computation/jc2_russell_cylinder_all_retained_cells_universal_second_stable_debt_thm3639.py
output: 05-knowledge/results/jc2_russell_cylinder_all_retained_cells_universal_second_stable_debt_thm3639.out
script_sha256: b2737e8654d088afd0c1928a0bb071851a33cc18c5928aa9149dd45226b5b6c8
output_sha256: 03e4a117ef2ab3bef9f30ea1c61e2092e12ec57a5aea65ba57cdc8047a84c680
hash_basis: raw LF bytes
---

# THM-3639 -- Russell-cylinder all-retained-cells universal second-stable debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem unifies the
rank-two cell of THM-3637 and the retained-tangent-rank-one boundary.  Its
load-bearing statement is the following stronger local identity:

```text
lambda(J_2)=-2072/81+(5/27)J_0'(-1)-(13/27)J_0'(1).     (1)
```

Here only the three retained **values** `J_0(-1)=J_0(0)=J_0(1)=1` are needed
to put the cell in normal form.  Thus if `J_0=1` as a polynomial, both
endpoint derivatives in `(1)` vanish and

```text
lambda(J_2)=-2072/81 != 0.                              (2)
```

In particular, any first-stable lift for which `J_1=0` fails at the next
stable coefficient.  In fact `(2)` does not require `J_1=0`; the quotient of
`J_2` is already forced by `J_0=1`.

All rings and closed points are over `C`.  The exact companion works over
`Q` and then base-changes.

## 0. Setup and exact scope

Use the exponent-two Danielewski surface and compiler

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,        c=xD(D+2),        e=q(D+3).       (3)
```

Fix the minimal THM-3635 Hermite polynomial

```text
Q=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4                    (4)
```

and write

```text
gamma(K)=K(b(x,Q(x)),c(x,Q(x)),e(x,Q(x))),
S=gamma(R),

delta(K)=partial_q K(b(x,q),c(x,q),e(x,q))|_(q=Q(x)).   (5)
```

Let

```text
F^#=sum_(i>=0) w^i F_i,       G^#=sum_(i>=0) w^i G_i,
F_i,G_i in R,                                             (6)
```

and pull back by the quadratic fold

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                              (7)
```

Put

```text
U=gamma(F_0),       V=gamma(G_0),
A=gamma(F_1),       B=gamma(G_1).                        (8)
```

The statement quantifies over every such `U,V,A,B in S`, with no degree
bound; every representative `F_0,G_0,F_1,G_1` of those restriction classes;
and every choice of `F_i,G_i` for `i>=2`.  Its hypothesis is

```text
J_0=U'B-AV'=1.                                           (9)
```

The retained derivative rank of `(U,V)` may be one or two.  Rank zero is
incompatible with `(9)`, so these are all cells.

## 1. Stable coefficients and the unrestricted retained quotient

Set

```text
C=delta(F_0)+gamma(F_2),      E=delta(G_0)+gamma(G_2),
D=delta(F_1)+gamma(F_3),      H=delta(G_1)+gamma(G_3).  (10)
```

The pulled-back functions have expansions

```text
f=U+tA+t^2C+t^3D+O(t^4),
g=V+tB+t^2E+t^3H+O(t^4),                                (11)
```

so

```text
J_0=U'B-AV',

J_1=2U'E+A'B-AB'-2CV',

J_2=3U'H+2A'E+C'B-AE'-2CB'-3DV'.                       (12)
```

At `x=-1,0,1`, the restricted curve has the common target point
`p=(b,c,e)=(0,0,-3)`.  Use regular local coordinates

```text
y=(c,z),                         z=e+3.                  (13)
```

The branch tangents, vertical compiler values, and one needed vertical
derivative are

```text
tau_c=(3,3,3),                 tau_z=(-9,4,9),

n_c=delta(c)=(2,0,-2),        n_z=delta(z)=(-2,4,-2),

n_c'=(-9,0,-9).                                         (14)
```

The pairwise tangent determinants are `(39,15,54)`.  Hence this is an
ordinary triple, its completed local image consists exactly of triples with
a common constant and derivative vector in

```text
T=span{tau_c,tau_z},                                    (15)
```

and its conductor in the three completed branches starts at branch order
two.  Fix the normalized cokernel row

```text
lambda(P)=(5P(-1)-18P(0)+13P(1))/18;                    (16)
```

it annihilates `T`.

Write `dK=(K'(-1),K'(0),K'(1))`, let `A_0,B_0` be the common retained
values of `A,B`, and put

```text
nu_K=(delta(K~)(-1),delta(K~)(0),delta(K~)(1)),
dot(nu_K)=((delta(K~))'(-1),
           (delta(K~))'(0),
           (delta(K~))'(1)),                            (17)
```

where `K~` is any target representative.  Applying `lambda` to `(12)` kills
every common value and every derivative in `T`.  Therefore all
`F_2,G_2,F_3,G_3` terms disappear and the exact unrestricted formula is

```text
lambda(J_2)=lambda(R_2),

R_2=3 dU*nu_B+2 dA*nu_V+B_0 dot(nu_U)-A_0 dot(nu_V)
       -2 nu_U*dB-3 nu_A*dV,                            (18)
```

where `*` means coordinatewise product.  Equivalently, for
`P=(U,V)` and `L=(A,B)`, branch by branch,

```text
R_2=d/dx det(delta(P),L)
       +3(det(P',delta(L))-det(delta(P),L')).            (19)
```

Thus before `(9)` is used, the quotient depends only on the common values of
`A,B`, the target two-jets of `F_0,G_0`, the target one-jets of `F_1,G_1`,
and the compiler data `(tau,n,n')`.  It does not depend on the zero-stable
values or on any coefficient from stable order two onward.

## 2. Representative invariance

Let `K in ker(gamma)`.  Its germ at `p` vanishes on all three retained
branches.  Its linear initial form vanishes on three distinct tangent lines
and is zero.  Its quadratic initial form also vanishes on those lines, while
a nonzero binary quadratic has at most two projective zeros.  Consequently

```text
ker(gamma) subset m_p^3.                                 (20)
```

It follows that the two-jets of `F_0,G_0` and the one-jets of `F_1,G_1` in
`(18)` are invariants of their restriction classes.  More explicitly,
differentiating a member of `m_p^3` once in a vertical target direction
leaves branch order at least two, so

```text
delta(K)(x_i)=0,            (delta(K))'(x_i)=0,
                            x_i=-1,0,1.                  (21)
```

This proves the full representative quantifier in the statement.  The
higher-lift quantifier was already discharged directly in `(18)`.

## 3. Exhaustive determinant-one normal form

Let `G` be the `2 x 2` target-gradient matrix whose rows are the gradients of
`F_0,G_0` in `(c,z)`, and let `s=(A_0,B_0)^T`.  For the three tangent columns
`t_i=(tau_c,i,tau_z,i)^T`, equation `(9)` at the retained points says

```text
det(G t_i,s)=1.                                          (22)
```

A constant output change `M in SL_2(C)`, applied to every coefficient pair
`(F_i,G_i)`, preserves the complete source Jacobian and acts by

```text
G |-> MG,                       s |-> Ms.                (23)
```

It also bijects all target representatives and higher lifts.

If `r=det(G) != 0`, take

```text
M=diag(1,r) G^(-1) in SL_2(C).                           (24)
```

Then `MG=diag(1,r)`.  Since the first tangent coordinate is always `3` and
the second takes the distinct values `-9,4,9`, equation `(22)` forces
`Ms=(0,1/3)^T`.

If `rank(G)=1`, factor `G=w alpha^T`.  Because
`det(w,s) alpha(t_i)=1` at three distinct second tangent coordinates,
`alpha` has zero `z`-component.  An element of `SL_2(C)` sends the remaining
nonzero output vector to `(1,0)^T`, giving `MG=diag(1,0)`.  Equation `(22)`
fixes the second component of `Ms` to `1/3`; the stabilizing shear
`[[1,h],[0,1]]` then kills its first component without changing `MG`.

Rank zero contradicts `(22)`.  Hence **every** cell admits, by a
Jacobian-preserving output change, the unified normal form

```text
G=diag(1,r),                   s=(0,1/3)^T,

r != 0  on tangent rank two,  r=0 on tangent rank one.  (25)
```

Output translations remove the common zero-stable values without affecting
any Jacobian coefficient.

## 4. Universal curvature identity

In the normal form `(25)`, write the only target jets that survive `(18)` as

```text
F_0=c+(1/2)u_cc c^2+u_cz cz+(1/2)u_zz z^2+O(m_p^3),
G_0=r z+O(m_p^2),

F_1=a_c c+a_z z+O(m_p^2),
G_1=1/3+b_c c+b_z z+O(m_p^2).                           (26)
```

No displayed remainder enters the retained quotient.  Let

```text
E_i=J_0'(x_i),                         x_i=-1,0,1.       (27)
```

Using

```text
c''|_triple=(157/2,0,-221/2)                            (28)
```

and `(14)`, direct two-jet expansion gives

```text
lambda(R_2)=
 (162r a_c+216r a_z-24b_c-162b_z
       -8u_cc-108u_cz-72u_zz-27)/9.                    (29)
```

The same variables in the derivative of `(9)` satisfy

```text
E_i=(1/3)(c_i''+t_i^T H_U t_i)
       +tau_c,i(b_c tau_c,i+b_z tau_z,i)
       -(a_c tau_c,i+a_z tau_z,i)r tau_z,i,             (30)

H_U=[[u_cc,u_cz],[u_cz,u_zz]].
```

Exact elimination has the polynomial identity, valid for every `r in C`,

```text
lambda(R_2)+2072/81=(5/27)E_- -(13/27)E_+.              (31)
```

This single identity covers `r!=0` and `r=0`; no division by the tangent
determinant occurs.  Because `SL_2` changes preserve `J_0,J_2` themselves,
`(31)` transports back to the original cell and proves `(1)`.  Under `(9)`,
`E_-=E_+=0`, proving `(2)`.

As active boundary controls, solving all three equations `E_i=0` gives, on
the rank-two stratum,

```text
a_c=(6b_z+4u_cz-7)/(6r),
a_z=(65u_zz-58)/(195r),
b_c=-u_cc/3-3658/1755,                                  (32)
```

while at rank one it gives

```text
u_cc=-3b_c-3658/585,
u_cz=7/4-(3/2)b_z,
u_zz=58/65.                                             (33)
```

Both substitutions return `-2072/81`.  These are checks on the two strata,
not extra hypotheses.

## 5. The role of `J_1` and stopping boundary

If the full pulled-back Jacobian were identically `1`, then its stable
coefficients would obey

```text
J_0=1,                         J_1=0,        J_2=0.       (34)
```

The theorem proves that the first equality alone forces the nonzero
quotient `(2)`, contradicting the third.  Therefore any cell which actually
reaches `J_1=0` has a genuine second-stable obstruction.  THM-3637 supplies
a globally exact rank-two `J_1` lift, so this implication is nonvacuous.
Local solvability of `J_1` on another cell is not asserted here.

The conclusion is confined to

```text
Q=Q_1, the compiler (3), the quadratic fold (7),
and coefficients in R.                                  (35)
```

It does **not** prove an obstruction for another collision polynomial, a
different fold, an arbitrary polynomial pair on `A^2`, the two-dimensional
Jacobian conjecture, or any Lonely Runner statement.  `JC(2)` remains OPEN.

## 6. Exact verification

The companion verifies the stable expansion, compiler and ordinary-triple
jets, the order-three representative lemma, higher-lift annihilation, both
`SL_2` normalizations, the unrestricted quotient `(18)`, the raw expression
`(29)`, the division-free identity `(31)`, both stratum eliminations, active
mutations, and an assertion-free AST.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_all_retained_cells_universal_second_stable_debt_thm3639.py
python3 -O 04-computation/jc2_russell_cylinder_all_retained_cells_universal_second_stable_debt_thm3639.py
```

Both streams must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_all_retained_cells_universal_second_stable_debt_thm3639.out
```

The frozen companion reports zero Python assertion statements.
