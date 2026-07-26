---
id: THM-2360
title: "Conditional degree-eighteen quadratic-ring cube descent"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Let the
  THM-2357 mixed residual satisfy
  (x^2-1)S_4^2=4(x-x_0)p_3^3+49q_5^2, with
  x_0^2!=1, Res(p_3,q_5)!=0, and q_5(x_0)!=0. In the Laurent UFD
  C[s,s^-1], s=x+sqrt(x^2-1), the conjugate norm factors are coprime.
  The polynomial B=s^5(7q_5+S_4 sqrt(x^2-1)) has exact degree ten,
  nonzero constant, and necessarily factors
  B=c(s-r_0)U_3(s)^3, where r_0 is one root above x_0 and U_3 has
  degree three. Repeated roots of p_3 and overlap p_3(x_0)=0 are
  allowed. A nonmaximal-generator hostile proves that THM-2357's branch
  signature alone does not supply Res(p_3,q_5)!=0. This is conditional
  structure, not an H_2 closure or a proof of JC(2).
source: codex-2026-07-25-degree-eighteen-quadratic-ring-cube
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2357-degree-eighteen-h2-moving-root-reduction
related:
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2347-degree-eighteen-double-zero-wall-saturation
script: 04-computation/jc2_degree18_quadratic_ring_cube_descent_thm2360.py
output: 05-knowledge/results/jc2_degree18_quadratic_ring_cube_descent_thm2360.out
script_sha256: f8756d29aa01b2951d6f471d50dbe19854c2c63a3ed2fd98d81879440da7071d
output_sha256: 6c1a05360651154c35b9f93d84af5005f1e0c9ea0e47e7b18e397778c378b87f
hash_basis: working-tree bytes (LF)
---

# THM-2360 -- the conditional quadratic-ring cube descent

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2357 turns the mixed degree-eighteen branch into

```text
H_2 S_4^2=4(x-x_0)p_3^3+49q_5^2,                 (1)

deg(H_2,p_3,q_5,S_4)=(2,3,5,4).
```

The squarefree quadratic `H_2` records the two transposition values and
`x_0` is the unique three-cycle value.  This theorem proves a strong
factorization of (1) under the two extra hypotheses

```text
Res_x(p_3,q_5)!=0,

q_5(x_0)!=0.                                       (2)
```

The first excludes nonmaximal common-root collisions.  The second is the
smooth polynomial-order chart; in THM-2357's moving-root coordinates it is
the condition `K!=0`.

Neither hypothesis is removed here.

## 1. Normalize the quadratic ring

Over `C`, an affine change of `x` sends the two distinct roots of `H_2` to
`-1,+1`.  Absorb the remaining nonzero scalar into `S_4`.  Equation (1)
becomes

```text
(x^2-1)S^2=4(x-x_0)p^3+49q^2,                    (3)

deg(p,q,S)=(3,5,4),

x_0^2!=1.                                          (4)
```

Put

```text
t^2=x^2-1,

s=x+t.                                             (5)
```

Since

```text
(x+t)(x-t)=1,
```

the inverse change is

```text
x=(s+s^(-1))/2,

t=(s-s^(-1))/2.                                    (6)
```

Consequently

```text
C[x,t]/(t^2-x^2+1)=C[s,s^(-1)],                   (7)
```

a Laurent UFD, and conjugation `t -> -t` is exactly

```text
s -> s^(-1).                                       (8)
```

## 2. The two norm factors are coprime

Define

```text
A=7q+St,

Abar=7q-St.                                        (9)
```

Equation (3) gives the exact norm identity

```text
A Abar
 =49q^2-(x^2-1)S^2
 =-4(x-x_0)p^3.                                   (10)
```

Suppose an irreducible Laurent prime `pi` divides both factors.  From
their sum and difference,

```text
pi divides q,

pi divides S t.                                    (11)
```

If `pi` divides `t`, its underlying `x`-value is `+1` or `-1`.  Equation
(10), together with `x_0^2!=1`, then makes `p` vanish there whenever `q`
does.  This contradicts `Res(p,q)!=0`.

If instead `pi` divides `S`, let `alpha` be its underlying `x`-value.
Since `q(alpha)=S(alpha)=0`, equation (10) gives

```text
(alpha-x_0)p(alpha)^3=0.                           (12)
```

The alternative `alpha=x_0` contradicts `q(x_0)!=0`; the other alternative
contradicts `Res(p,q)!=0`.  These cases exhaust (11), so

```text
gcd(A,Abar)=1 in C[s,s^(-1)].                     (13)
```

## 3. Clear the Laurent endpoints

Put

```text
B(s)=s^5 A(s).                                     (14)
```

Let `a,b,c` be the leading coefficients of `p,q,S`.  From (3),

```text
c^2=4a^3+49b^2.                                   (15)
```

Using (6), the two endpoint coefficients of `B` are

```text
lc(B)=(7b+c)/32,

B(0)=(7b-c)/32.                                    (16)
```

Their product is

```text
lc(B) B(0)
 =(49b^2-c^2)/1024
 =-a^3/256 !=0.                                   (17)
```

Thus no cancellation occurs at either Laurent endpoint:

```text
B belongs to C[s],

deg(B)=10,

B(0)!=0.                                           (18)
```

Write

```text
B^vee(s)=s^10 B(s^(-1))=s^5 Abar(s).              (19)
```

The factors in (13) remain coprime after multiplication by Laurent units,
so

```text
gcd(B,B^vee)=1.                                    (20)
```

## 4. Reciprocal-pair allocation

Define

```text
X_0(s)
 :=s(x-x_0)
  =(s^2-2x_0s+1)/2,

P_6(s)
 :=s^3p((s+s^(-1))/2).                             (21)
```

The polynomial `P_6` has exact degree six, nonzero constant, and is
reciprocal:

```text
s^6P_6(s^(-1))=P_6(s).                            (22)
```

Equations (10), (14), and (19) become

```text
B B^vee=-4X_0 P_6^3.                              (23)
```

The two roots of `X_0` are a reciprocal pair

```text
{r_0,r_0^(-1)}.                                    (24)
```

They are distinct because `x_0^2!=1`.  Every root of `p`, counted with
multiplicity, likewise gives one reciprocal pair of roots of `P_6`.
No root of `p` is `+1` or `-1`: equation (3) would then force `q` to
vanish at the same point, contrary to the resultant hypothesis.  Thus
none of these pairs collapses to a self-reciprocal root.

For a reciprocal pair `{rho,rho^(-1)}`, the exponent of `B^vee` at `rho`
is the exponent of `B` at `rho^(-1)`.  Coprimality (20) therefore forces
all the multiplicity in (23) onto exactly one member of each pair.

Choose the member of (24) occurring in `B` and call it `r_0`.  From each
of the three `p`-pairs, counted with multiplicity, choose the member
occurring in `B`, and let their product be `U_3`.  Then

```text
deg(U_3)=3,

B(s)=c_0(s-r_0)U_3(s)^3                           (25)
```

for a nonzero constant `c_0`.

This argument does not require `p` to be squarefree.  If a root of `p`
has multiplicity `m`, its selected Laurent root occurs in `U_3` with
multiplicity `m` and in `B` with multiplicity `3m`.

It also permits

```text
p(x_0)=0.                                          (26)
```

In that case the distinguished reciprocal pair is also a `p`-pair.
Coprimality forces the linear and cubic allocations to choose the same
member, so that member has multiplicity `1+3m` in (25).  The companion's
exact hostile-positive control realizes multiplicity four.

There is no branch-at-infinity exception: exact degrees of `p,q,S` and
the nonzero endpoint product (17) give both degree ten and a nonzero
constant before the root allocation.

## 5. Why THM-2357 does not supply the resultant

The field branch signature controls the normalization, not the chosen
polynomial order.  Start with

```text
w^3+p_0(x)w+q_0(x)=0
```

and choose `beta`.  The nonmaximal generator

```text
v=(x-beta)w
```

satisfies

```text
v^3+(x-beta)^2p_0(x)v+(x-beta)^3q_0(x)=0.          (27)
```

The function field and normalized cover are unchanged, but the new
depressed coefficients have a common root.  Moreover their polynomial
discriminant is

```text
(x-beta)^6[-4p_0^3-27q_0^2],                       (28)
```

so the defect is exactly an index square invisible to the field branch
signature.

The companion uses

```text
beta=2,             p_0=x+1,             q_0=x^2+1
```

and verifies

```text
gcd((x-2)^2p_0,(x-2)^3q_0)=(x-2)^2,

q(0)!=0.                                            (29)
```

Thus even the additional condition `q(x_0)!=0` does not recover
coprimality.  A degree-eighteen proof of

```text
Res(p_3,q_5)!=0                                    (30)
```

or an exact stripping theorem for the common index factor is a genuine
remaining obligation.

## 6. Frontier effect and scope

Under (2), the `H_2S_4^2` coefficient problem is no longer an arbitrary
six-equation factorization.  It has the rigid Laurent form

```text
s^5(7q_5+S_4t)=c_0(s-r_0)U_3^3.                  (31)
```

This replaces the square-class constraint by a linear-times-cube
coefficient comparison of degree ten.  It is a promising consumer of
THM-2357's moving-root and pivot coordinates.

However, THM-2357 does not currently prove (30), and its `K=0` boundary
does not satisfy the second hypothesis in (2).  Therefore this theorem
does not close either THM-2357 lane, the `H_2` stratum, the `H_4` stratum,
or any degree-eighteen Keller branch.  It proves neither `JC(2)` nor
`DC(2)`.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_quadratic_ring_cube_descent_thm2360.py
python3 -O 04-computation/jc2_degree18_quadratic_ring_cube_descent_thm2360.py
```

Both transcripts are byte-identical to

```text
05-knowledge/results/jc2_degree18_quadratic_ring_cube_descent_thm2360.out
```

The companion verifies the Laurent parametrization, conjugation, formal
norm, endpoint coefficients, leading relation, reciprocal sextic, generic
pair allocation, linear/cube overlap, repeated-root allocation, and the
nonmaximal-generator hostile discriminant multiplier.  No executable
check uses Python `assert`.

Independent audit is pending. QED.
