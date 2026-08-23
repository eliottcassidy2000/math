---
id: THM-3807
title: "Trinomial and full cubic R-repairs have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the c=1
  cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+r sum_(i=0)^3 b_i e^i whose profile has at least three nonzero
  coefficients has a critical point.  For b_3 outside {0,1}, the
  degree-seventeen logarithmic remainder has two top branches; their only
  surviving parameter seam has support two.  At b_3=1 the residual degree
  drops to fifteen or eleven, and separately recomputed divisions are still
  contradictory.  Thus no trinomial or full cubic r-repair has a regular
  Darboux mate.  Support at most two is outside this theorem but already
  closed by THM-3799 and THM-3806; together with THM-3805, the canon now
  closes every polynomial r-profile of degree at most three.  Degree at
  least four and mixed corrections remain open.
source: jc_sparse_direct_search / full-cubic logarithmic-resultant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn-boundary, 2026-08-23).  The
  universal P/Q resultant, all generic
  remainder factors, both b_3=1 degree drops, exact coprime resultants,
  support-two seams, homogeneous finite-root gate, excluded denominators,
  source reconstruction, and Casimir implication were rederived symbolically.
  An independently assembled fixed-degree Sylvester determinant reproduces
  the universal resultant.  The deterministic companion has 51 active gates;
  normal and optimized runs LF-normalize exactly to the frozen transcript
  and the declared raw-LF hashes match.  The audit independently checked the support
  partition, both branch exhaustions, all exceptional specializations,
  coprime resultants, multiplicity-safe boundary implication, and finite
  reconstruction.  The generic quotient/remainder common denominator is
  4(t-1)^3; the t=1 denominators are b_2 and 531441b_0^6 on the two cells,
  all invertible under their stated hypotheses.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points
related:
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3806-binomial-cubic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_trinomial_full_cubic_r_repair_thm3807.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_trinomial_full_cubic_r_repair_thm3807.out
script_sha256: 6a5c70d8bd728d29599315d4380d2a6e0c907cee94d08b829d354e06d87fc322
output_sha256: faf6f2a0933a54d85c5ee7968c6b31c9e5b00bc6ec43b156957f04977ba44da7
semantic_sha256: 692b1c1f9c260b8d9154b85d386afb82b6b4acff87093a41b533430c7dae92a7
hash_basis: raw LF bytes
---

# THM-3807 -- every cubic profile of support at least three remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

Let

```text
g(e)=b_0+b_1e+b_2e^2+b_3e^3,
A=e^2-z/3+r g(e).                                      (2)
```

If at least three of `b_0,b_1,b_2,b_3` are nonzero, then `A` has a critical
point on `Y`.  Consequently `A` has no regular Darboux mate.

When `b_3=0`, this is THM-3805.  We may therefore assume

```text
t=b_3!=0.                                               (3)
```

The proof below does not use the separately audited status of the support-two
rows in THM-3806.  Those rows appear only as excluded equality boundaries of
the present support-at-least-three statement.

## 1. Universal residual and the multiplicity-safe boundary test

Put

```text
u=re,                         K=1+2u,
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (4)
```

The Hamiltonian components are

```text
{A,r}=r^2-9z^2(2e+r g'),
{A,z}=3g r^2-3(1+2re)(2e+r g'),
{A,e}=9g z^2-(1+2re).                                 (5)
```

After `r=u/e`, the equation `P=0` is `{A,z}=0` up to `3/e^2`, while `Q=0`
is the compatibility equation for

```text
z^2=K/(9g),                       z^3=u(1+u)/e.         (6)
```

Keep independent symbols `G,D` in place of `g,g'` before eliminating `u`.
The universal resultant is

```text
Res_u(P,Q)=G^3e^4 H_univ(e,G,D).                       (7)
```

After `(G,D)=(g,g')`, write the residual polynomial as `H(e)`.

Suppose every root of a nonconstant `H` lay on the forbidden boundary
`V(e g)`.  Factoring with multiplicities over `k` gives

```text
H=mu product_alpha (e-alpha)^(n_alpha),
e g H'/H=sum_alpha n_alpha e g/(e-alpha) in k[e].      (8)
```

Thus boundary-only support necessarily implies

```text
H divides e g H'.                                      (9)
```

Only this necessary direction is used.  It remains valid when either `g` or
`H` has repeated roots.

## 2. The generic stratum `t notin {0,1}`

Assume first

```text
t(t-1)!=0.                                             (10)
```

Then

```text
deg H=17,                 LC(H)=8503056t^3(t-1)^2.    (11)
```

Divide `e g H'` by `H` over the localized coefficient ring in which
`t(t-1)` is invertible, and call the remainder `R`.

### 2.1 The cell `b_2=0`

Here support at least three forces `b_0b_1!=0`.  But

```text
[e^15]R=25509168 b_0b_1 t^3(t-1)!=0.                 (12)
```

So `(9)` is impossible.

### 2.2 The cell `b_2!=0`

The top remainder coefficient is

```text
[e^16]R=-1062882 b_2t^3 E_16/(t-1)^2,
E_16=16b_0(t-1)^2-4b_1b_2(t-1)^2+b_2^3(t-2).         (13)
```

If `(9)` held, `(13)` would give

```text
b_0=b_2 Phi/[16(t-1)^2],
Phi=4b_1(t-1)^2-b_2^2(t-2).                          (14)
```

Put also

```text
Psi=12b_1(t-1)^2-b_2^2(t-2).                         (15)
```

After `(14)`, the next coefficient factors exactly as

```text
[e^15]R=531441 b_2t^3 Phi Psi/[4(t-1)^3].            (16)
```

There are two branches.

If `Phi=0`, then

```text
b_0=0,                  b_1=b_2^2(t-2)/[4(t-1)^2],   (17)
[e^13]R=20412t^2(t-2)(9t^2-20t+20),
[e^12]R=2916b_2t(195t^4-838t^3+1380t^2-1080t+352)/(t-1).   (18)
```

The value `t=2` in `(17)` forces `b_0=b_1=0`, leaving exactly the two
monomials `b_2e^2+2e^3`; it is excluded by the support hypothesis.  Otherwise
vanishing in `(18)` would force the two polynomials

```text
9t^2-20t+20,
195t^4-838t^3+1380t^2-1080t+352                       (19)
```

to share a root.  Their exact resultant is

```text
210898944!=0.                                           (20)
```

If `Psi=0`, then

```text
b_0=-b_2^3(t-2)/[24(t-1)^2],
b_1= b_2^2(t-2)/[12(t-1)^2],                           (21)
[e^14]R=59049b_2^6t^3(t-2)^2(9t-10)/[8(t-1)^4].       (22)
```

Again `t=2` leaves support two.  Otherwise `(22)` forces `t=10/9`, but at
that value

```text
[e^13]R=-1792000/9!=0.                                 (23)
```

This closes every support-at-least-three profile under `(10)`.

## 3. The degree-drop stratum `t=1`

The generic leading coefficient and quotient used above vanish or acquire a
pole at `t=1`.  We therefore specialize `H` first and redo each Euclidean
division in its actual coefficient ring.

### 3.1 The cell `b_2!=0`

Now

```text
deg H=15,                    LC(H)=2125764b_2^2.       (24)
```

For the recomputed remainder `R_1`,

```text
[e^14]R_1=8503056b_0(b_0+b_1b_2).                    (25)
```

If `b_0=0`, the next coefficient is

```text
[e^13]R_1=-131220!=0.                                 (26)
```

On the other branch `b_0=-b_1b_2`, the next two coefficients are

```text
[e^13]R_1=26244(648b_1^2b_2^3-5),
[e^12]R_1=17496b_2(2916b_1^2b_2^3-5).                (27)
```

The first vanishing law in `(27)` makes the second coefficient

```text
306180b_2!=0.                                          (28)
```

### 3.2 The cell `b_2=0`

Support at least three now means `b_0b_1!=0`.  The residual degree drops
again:

```text
deg H=11,                    LC(H)=2125764b_0^2.       (29)
```

The top recomputed remainder coefficient is

```text
[e^10]R_1=
 4(3483891b_0^7-54675b_0^4b_1+1)/(81b_0^6).          (30)
```

If `(30)` vanished, then

```text
b_1=(3483891b_0^7+1)/(54675b_0^4).                    (31)
```

Set `X=b_0^7`.  After `(31)`, the next two coefficients are nonzero scalar
multiples of

```text
F(X)=1480774572985482X^2+867784104X-53,
G(X)= 579994041785274X^2- 28717497X+4.                (32)
```

More exactly, they are `-4F/(4100625b_0^8)` and
`-4G/(110716875b_0^10)`.  The exact resultant

```text
Res_X(F,G)=2408049135195443414554999563281250!=0       (33)
```

finishes the last degree-drop cell.  Equations `(24)--(33)` are specialized
divisions, not substitutions into the localized quotient of Section 2.

Sections 2 and 3 show that `(9)` is impossible whenever the cubic profile has
support at least three.  Hence `H` has a root `eta` with

```text
eta g(eta)!=0.                                         (34)
```

## 4. Every off-boundary root reconstructs a critical point

At `e=eta`, equations `(7)` and `(34)` make the fixed-degree projective
`u`-resultant vanish.  The quartic retains

```text
LC_u(Q)=-729g(eta)^3!=0,                               (35)
```

so the unique point at infinity cannot be a common root, even if the
quadratic leading coefficient of `P` vanishes.  Thus `P,Q` have a common
finite root `u_0`.  Direct substitution gives

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16.                           (36)
```

Therefore `u_0`, `1+u_0`, and `K_0=1+2u_0` are nonzero.  Define without a
root choice

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                      (37)
```

Before imposing `P=Q=0`, exact reduction gives

```text
z_0^2-K_0/(9g(eta))=-Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0=u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0=-Q/(eta^2K_0^2),
{A,z}|_0=3P/eta^2.                                    (38)
```

Thus `(r_0,z_0,eta)` lies on `Y` and the last two Hamiltonian components
vanish.  The surface Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                     (39)
```

and `K_0!=0` kill the remaining component.  Hence `(37)` is an actual
critical point.

The theorem itself closes precisely the cubic profiles with at least three
nonzero coefficients.  Combining THM-3792 (the zero profile), THM-3799
(monomials), THM-3805 (degree at most two), THM-3806 (cubic binomials), and
this theorem gives the useful synthesis

```text
deg g<=3  ==>  e^2-z/3+r g(e) has a critical point.   (40)
```

Degree at least four, mixed `z^2h(e)+r g(e)` corrections, other arm
profiles, and rational mates with poles remain open.  The exact companion
named in the metadata verifies `(4)--(39)`, both specialized divisions, the two coprime
resultants, support-drop seams, and reconstruction with 51 active gates.
Normal and optimized executions LF-normalize exactly to the frozen transcript.
**QED.**
