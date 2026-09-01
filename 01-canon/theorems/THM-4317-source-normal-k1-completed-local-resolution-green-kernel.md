---
id: THM-4317
title: "Source-normal k=1 completed-local resolution and Green kernel"
status: >
  PROVED RELATIVE TO THM-4312 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  On each of THM-4312's two k=1 equality rays, the completed-local first-
  carrier singularities have an exhaustive resolution by rational curves.
  The L1-nonzero carrier has two A55 boundary points and one 1/3(1,1)
  quotient point; the L1=0, L2-nonzero carrier has two A56 boundary points
  and two A1 quotient points.  Minimal resolution of the normalized strict
  transform adds respectively 111 and 114 rational curves.  Conditional on
  an actual Keller lift realizing the
  finite datum, every completed-local divisorial refinement is therefore
  constant.  Each A-chain is exactly the Dirichlet path Laplacian; its
  inverse is half the absorbing-walk occupation kernel, giving a precise
  valuation Poisson formula.  This is conditional local geometry, not an
  all-row lift, seam-entry, JC(2), or DC(2) theorem.
source: root / planar-Jacobian stochastic-process bridge session, 2026-09-01
audit: >
  PASS.  The primary SymPy certificate reconstructs the literal local
  equation, both weighted initials, finite boundary Hessians and critical
  orders, infinity covers and saturations, quotient actions, curve counts,
  and the path-Laplacian Green/Poisson identities for r=56,57.  A separate
  standard-library Fraction/sparse-polynomial implementation repeats these
  calculations, checks chart exhaustion, verifies the exact critical-value
  Euler identity, and audits additional path lengths.  Normal and optimized
  runs byte-match both frozen transcripts.
depends_on:
  - THM-4312-source-normal-cubic-corner-repeated-face-collapse
related:
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4316-source-normal-row-ten-cubic-corner-extinction
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
primary_script: 04-computation/jc2_source_normal_k1_completed_local_resolution_green_kernel_thm4317.py
primary_output: 05-knowledge/results/jc2_source_normal_k1_completed_local_resolution_green_kernel_thm4317.out
primary_script_sha256: 4cf9f3c25bdb80c518091a35cd0e2a9a8a7eff649a308dbc25b979c9ff34d5d4
primary_output_sha256: cd329b649d39ef81ac56cf5f4add01b3e10f8b728a6ee9d6e9b73a5129d3faf9
independent_audit_script: 04-computation/jc2_source_normal_k1_completed_local_resolution_thm4317_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_k1_completed_local_resolution_thm4317_independent_audit.out
independent_audit_script_sha256: 5b80c337b564c9773c35ec4adcc26f833a633265a4b126532eb0b291506aef69
independent_audit_output_sha256: 0aa731a69d3f5ddc2de0d1f691a37585ee71ca89667e249eea8b13f6f218d282
hash_basis: raw LF bytes
---

# THM-4317 -- Source-normal k=1 completed-local resolution and Green kernel

**PROVED RELATIVE TO THM-4312 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
CONDITIONAL ON AN ACTUAL KELLER LIFT REALIZING THE FINITE `k=1` DATUM, EVERY
COMPLETED-LOCAL DIVISOR ABOVE EITHER FIRST CARRIER IS CONSTANT.  THIS DOES
NOT ASSERT SUCH A LIFT, SEAM ENTRY, `JC(2)`, OR `DC(2)`.**

## 1. Literal local equation and the two rays

Retain THM-4312's exact cubic-corner `k=1` datum.  Write

```text
rho=-alpha_11/(2U),             c_2=U rho^2,
x=q-rho t,                      c_3=eta+zeta_3,
c_4=Delta+Theta.                                          (1)
```

In the local coordinates used there, direct expansion of the literal source
has the form

```text
F=-t^12/2+q(h(t)-z^12)+q^2 A(t)+q^3 Ccal(q,t),
h(t)=c_2t^2+c_3t^3+c_4t^4+...,
A(t)=O(t),                       Ccal(0,0)=U.             (2)
```

The first two critical-value coefficients are

```text
L_1=c_3+rho upsilon_5,

C=c_4-rho zeta_3
 =Delta+Theta+alpha_11 zeta_3/(2U),

L_2=C-upsilon_5^2/(4U)             when L_1=0.           (3)
```

THM-4312 proves that exactly one of `L_1!=0` or `L_1=0,L_2!=0` occurs and
that the corresponding first carrier is elliptic and Keller-constant for
any actual lift.  The task here is to resolve every point of those carriers,
including their weighted points at infinity.

The connection ledger is exact.  Its source is the intersection matrix and
valuation vector on the surface resolution; its target is a killed
nearest-neighbour chain.  The map is the discrete Laplacian in Section 6.
It preserves boundary orders and strict-transform contact, loses the cyclic
quotient points, and therefore carries those quotient curves as a separate
sidecar.  The cheapest hostile is the `1/3(1,1)` point: it is not part of an
unbiased `A`-chain.

## 2. The `L_1!=0` finite chart: two `A_55` points

For weights

```text
wt(x,t,z)=(6,4,1),
```

the exact lowest form of the literal equation is

```text
in_16(F)=rho t(Ux^2-z^12+L_1t^3).                       (4)
```

Set `t=wz^4`, `x=z^6Y`, and divide by `z^16`.  At `z=0` the strict
transform is

```text
rho w(UY^2-1+L_1w^3)=0.                                (5)
```

The elliptic component `UY^2=1-L_1w^3` meets the rational component `w=0`
at exactly two points

```text
w=0,                         Y=+/-U^(-1/2).              (6)
```

At either point the Hessian in `(w,Y)` has determinant

```text
-(2rho UY)^2=-4rho^2U !=0.                               (7)
```

The parameterized Morse lemma therefore leaves one transverse function of
`z`.  The exact critical equations in the original `(q,t)` coordinates give

```text
t=T_kappa(z)=kappa z^6+O(z^12),       c_2 kappa^2=1,
q=3kappa^12 z^60+O(z^66),                               (8)

F_crit=-T_kappa(z)^12/2+O(z^126)
      =-(kappa^12/2)z^72+O(z^78).                        (9)
```

The two roots of `c_2kappa^2=1` correspond to the two points in `(6)`.
After the `z^16` division, the transverse order is `72-16=56`.  Thus both
boundary points are

```text
uv=z^56,                         type A_55.               (10)
```

The distinction between the two remainders in `(9)` is load-bearing.  The
`O(z^126)` statement is for `F_crit+T_kappa^12/2`; after expanding
`T_kappa=kappa z^6+O(z^12)` in the literal `z`, the next possible term is
order 78.

## 3. The `L_1!=0` point at infinity

On the equality ray write `t=sigma z`; the actual weighted variables
`(sigma,z,x)` have weights `(3,1,6)`.  Use the cyclic cover

```text
sigma=lambda^3,        z=lambda u,        x=lambda^6X,
X=uV.                                                       (11)
```

Every literal monomial is divisible by the required exceptional power; after
the further common `u^3` saturation, the equation at `lambda=0` is

```text
UV^2+L_1u-u^10=0.                                      (12)
```

It is smooth.  The `mu_3` action keeping `(sigma,z,x)` fixed is

```text
(lambda,u,V) -> (omega lambda,omega^2u,omega V).         (13)
```

Its only fixed point on `(12)` is `(u,V)=(0,0)`.  Since `L_1!=0`, `u` is
solved there as a function of `V`, so local surface coordinates are
`(lambda,V)` with weights `(1,1)`.  Downstairs this is one cyclic quotient

```text
1/3(1,1),                                                (14)
```

whose minimal resolution is one rational curve of self-intersection `-3`.

## 4. The `L_1=0`, `L_2!=0` ray

Now use weights `(6,3,1)`.  Directly from the same literal `F`,

```text
in_15(F)=rho t(Ux^2+upsilon_5t^2x+Ct^4-z^12).           (15)
```

Set `t=wz^3`, `x=z^6Y`, and then

```text
Ybar=Y+upsilon_5w^2/(2U).
```

The exceptional equation becomes

```text
rho w(UYbar^2-1+L_2w^4)=0.                              (16)
```

There are again two points at `w=0`, with the same nondegenerate Hessian
`(7)`.  The critical value still has order 72, but the chart division is now
`z^15`; hence the transverse order is 57 and the two points have type

```text
uv=z^57,                         two A_56 points.         (17)
```

At infinity use

```text
 t=sigma z,             wt(sigma,z,x)=(2,1,6),
sigma=lambda^2,        z=lambda u,        x=lambda^6X,
X=u^2V.                                                     (18)
```

After the common `u^5` saturation, the exceptional equation is

```text
UV^2+upsilon_5V+C-u^8=0.                                (19)
```

At `u=0` its discriminant is

```text
upsilon_5^2-4UC=-4UL_2 !=0.                             (20)
```

Thus there are two smooth cover points.  The involution acts in local
coordinates as `(lambda,u)->(-lambda,-u)`, so each quotient is

```text
1/2(1,1)=A_1.                                           (21)
```

Each contributes one rational `-2` curve.

## 5. Exhaustion and completed-local consequence

The finite `z`-chart and the `sigma`-chart above are the two relevant charts
of the weighted carrier.  The cubic carrier has one point at infinity and
the quartic carrier has two; `(12)` and `(19)` account for them.  The third,
`x`-chart is empty: at `z=0` both homogeneous carrier equations reduce to
`Ux^2=0`, so `x=0` and the point lies in the `sigma`-chart.  Away from the
listed points, the strict transforms are smooth.

Consequently the minimal new exceptional inventories are

```text
L_1!=0:       2*(55 A-chain curves)+1 quotient curve=111,

L_1=0:        2*(56 A-chain curves)+2 quotient curves=114. (22)
```

Every curve in `(22)` is rational.  Any further divisorial valuation over a
smooth resolved surface is obtained after point blowups, which add only
rational curves.  THM-4312 already makes the elliptic first carrier constant;
every morphism from a rational curve to the good elliptic target is constant
in characteristic zero.  Hence, conditional on an actual Keller lift
realizing the finite datum, every completed-local divisorial refinement over
either `k=1` equality carrier is constant.

This completed-local conclusion is structural.  THM-4316 independently
shows that these finite corner data already fail the row-ten bracket gate.

## 6. Resolution chains as absorbing stochastic processes

For either boundary singularity, write

```text
A_(r-1): uv=z^r,                 r=56 or 57.             (23)
```

Orient its exceptional chain `C_1,...,C_(r-1)` so that

```text
nu_i(u)=r-i,          nu_i(v)=i,          nu_i(z)=1.     (24)
```

The negative intersection matrix is the Dirichlet path Laplacian

```text
L=2I-adjacency.                                            (25)
```

Its inverse is

```text
(L^(-1))_(ij)
 =min(i,j)(r-max(i,j))/r,               1<=i,j<r.        (26)
```

Let `X_n` be the unbiased nearest-neighbour walk on `{0,...,r}`, absorbed at
`tau=inf{n:X_n in {0,r}}`.  Then

```text
P_i(X_tau=r)=i/r,

E_i[number of visits to j before tau]
 =2(L^(-1))_(ij).                                        (27)
```

Thus every monomial valuation in `(24)` is affine in `i` and gives a stopped
martingale.  More generally, suppose divisor orders `a_i` and the effective
strict-transform contact vector `b_i` satisfy

```text
2a_i-a_(i-1)-a_(i+1)=b_i.                               (28)
```

Equations `(26)--(27)` give the exact Poisson representation

```text
a_i=E_i[a_(X_tau)]
    +(1/2)E_i[sum_(n<tau)b_(X_n)].                       (29)
```

Equivalently,

```text
a_(X_(n wedge tau))
+(1/2)sum_(k<n wedge tau)b_(X_k)                         (30)
```

is a martingale.  Effective divisor orders are superharmonic; monomials are
the harmonic case `b=0`.  This is a lawful stochastic model because the
geometry itself supplies the state graph, Laplacian, boundary values, and
contact sidecar.  It is not the automatic terminal-law realization of
THM-3163.

The cyclic quotient curve in `(14)` and the two `A_1` curves in `(21)` are
not part of the unbiased chains `(23)`.  A killed or sub-Markov extension
could encode the `-3` weight, but no such identification is used here.

## 7. Reproduction and scope

Run

```text
python3 -B 04-computation/jc2_source_normal_k1_completed_local_resolution_green_kernel_thm4317.py
python3 -B -O 04-computation/jc2_source_normal_k1_completed_local_resolution_green_kernel_thm4317.py
python3 -B 04-computation/jc2_source_normal_k1_completed_local_resolution_thm4317_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_k1_completed_local_resolution_thm4317_independent_audit.py
```

The normal and optimized streams must byte-match their frozen outputs.

This theorem resolves the completed-local surface geometry above THM-4312's
two finite first carriers.  It does not assert that either finite datum comes
from an all-row depth-module element or a polynomial Keller pair, prove seam
entry, cross the `U=0` or `Z=0` walls, or imply `JC(2)` or `DC(2)`.

**QED.**
