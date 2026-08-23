---
id: THM-3805
title: "Quadratic R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the c=1
  cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+r(b_0+b_1e+b_2e^2) has a critical point.  In the genuine
  quadratic case, the degree-fourteen residual resultant cannot have all
  roots on the forbidden boundary e g(e)=0: boundary-only support would
  force H to divide e g H', but the four leading remainder equations split
  into two explicitly contradictory branches.  Hence no quadratic r-repair
  has a regular Darboux mate.  Degree at least three and mixed corrections
  remain open.
source: jc_zero_debt_lift / quadratic-r boundary-divisibility lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The universal P/Q
  reduction, degree-fourteen residual resultant, repeated-root-safe
  divisibility necessity, exact quotient and full remainder, E13--E10
  two-branch contradiction, homogeneous finite-root seam, denominator
  exclusions, source reconstruction, and lower-degree inheritance were
  independently rederived.  Finite-field scans remain hostile controls only.
  The companion has 37 active gates; normal and optimized runs LF-normalize
  exactly to the frozen transcript and all declared raw-LF hashes match.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3803-affine-linear-r-repairs-of-nodal-carriers-have-critical-points
related:
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
  - THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points
script: 04-computation/jc2_cubic_pseudoplane_quadratic_r_repair_thm3805.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_quadratic_r_repair_thm3805.out
script_sha256: 39eabc9de1182aad0bcec183ba6dadcee9450197e3841773243ac08ee63e874c
output_sha256: ba8d59e6b7fd0e9b7c560a0fdc70d79ac3edde209ffee6d56f5fcdc11821fb39
semantic_sha256: 777f813d0d8f54910c09c04bd5b57a499b49a60a8e00d473225ab8951813b02e
hash_basis: raw LF bytes
---

# THM-3805 -- every quadratic r-repair remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary `b_0,b_1,b_2 in k`, the regular function

```text
A=e^2-z/3+r g(e),             g(e)=b_0+b_1e+b_2e^2   (2)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

If `b_2=0` and `(b_0,b_1)!=(0,0)`, this is THM-3803.  If all three
coefficients vanish, it is the pure first-normal row of THM-3792.  It
therefore remains only to prove the genuinely quadratic case

```text
b_2!=0.                                                (3)
```

## 1. The universal residual resultant

Put

```text
u=re,                         K=1+2u,                 (4)
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (5)
```

As in THM-3799 and THM-3803, `P/e^2={A,z}/3` after `r=u/e`, while `Q`
is the compatibility equation for

```text
z^2=K/(9g),                     z^3=u(1+u)/e.          (6)
```

Keep independent symbols `G,D` in place of `g,g'` for the elimination.
The universal `u`-resultant factors as

```text
Res_u(P,Q)=G^3e^4 H_univ(e,G,D).                       (7)
```

After substituting the quadratic `g` from `(2)`, write the residual factor
as `H(e)`.  Exact expansion gives

```text
deg H=14,                       LC(H)=8503056b_2^3.    (8)
```

Thus `(3)` makes `H` a nonzero polynomial of fixed degree fourteen.

## 2. Boundary-only roots force logarithmic divisibility

Suppose for contradiction that every root of `H` lies on the forbidden
boundary

```text
V(e g(e)).                                               (9)
```

This implication must include repeated roots, so a radical-only count is
not enough.  Factor over `k`:

```text
H=mu product_alpha (e-alpha)^(n_alpha).                 (10)
```

Every `alpha` in `(10)` is a zero of `e g`.  Hence

```text
e g H'/H=sum_alpha n_alpha e g/(e-alpha) in k[e].       (11)
```

Equivalently, boundary-only support necessarily implies

```text
H divides e g H'.                                       (12)
```

Only this necessary direction is used.  It is valid whether `g` is
squarefree or has a double root.

Divide `e g H'` by `H`.  The quotient is

```text
(16b_0+2b_1b_2+b_2^3
 +(22b_1+2b_2^2)e+28b_2e^2)/2.                         (13)
```

If `(12)` held, every coefficient of the remainder would vanish.  Its four
highest coefficients are nonzero scalar multiples of `b_2^4E_13`,
`b_2^3E_12`, `b_2^2E_11`, and `b_2E_10`, where

```text
E_13 = 8b_0-2b_1b_2-b_2^3,

E_12 = 72b_0b_1-4b_0b_2^2-12b_1^2b_2
       -4b_1b_2^3+b_2^5,

E_11 = 17496b_0^2b_2+29160b_0b_1^2-7776b_0b_1b_2^2
       -1458b_0b_2^4-2916b_1^3b_2+729b_1b_2^5+80,

E_10 = 58320b_0^2b_1b_2-3888b_0^2b_2^3+21384b_0b_1^3
       -14580b_0b_1^2b_2^2-2430b_0b_1b_2^4+243b_0b_2^6
       -972b_1^4b_2+972b_1^3b_2^3+729b_1^2b_2^5
       +64b_1-56b_2^2.                                 (14)
```

The exact scalar factors in descending degree are respectively
`-2125764`, `-1062882`, `-4374`, and `-4374`; they are nonzero in
characteristic zero.  Since `b_2!=0`, divisibility would force

```text
E_13=E_12=E_11=E_10=0.                                (15)
```

## 3. The two branches are both impossible

The first equation in `(15)` gives

```text
b_0=b_2(2b_1+b_2^2)/8.                                (16)
```

After this substitution,

```text
E_12=b_2(2b_1+b_2^2)(6b_1+b_2^2)/2.                  (17)
```

There are two branches.

* If `2b_1+b_2^2=0`, then `(16)` gives `b_0=0`, but direct substitution
  gives

  ```text
  E_11=80!=0.                                          (18)
  ```

* If `6b_1+b_2^2=0`, then `(16)` gives `b_0=b_2^3/12`, and

  ```text
  E_11=5(27b_2^7+32)/2,
  E_10=-5b_2^2(81b_2^7+80)/6.                         (19)
  ```

  The first equation would force `27b_2^7+32=0`; under that law the second
  becomes

  ```text
  E_10=40b_2^2/3!=0.                                  (20)
  ```

Both branches contradict `(15)`.  Therefore `(12)` is impossible, and `H`
has a root `eta` with

```text
eta g(eta)!=0.                                         (21)
```

## 4. The off-boundary root reconstructs a critical point

At `e=eta`, equations `(7)--(8)` make the `u`-resultant vanish.  The leading
coefficient of the quartic `Q` is

```text
LC_u(Q)=-729g(eta)^3!=0.                               (22)
```

Hence the homogenized pair has no common root at infinity, even if the
quadratic coefficient of `P` happens to vanish.  There is a common finite
root `u_0`.  Directly,

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16,                           (23)
```

so `u_0`, `1+u_0`, and `K_0=1+2u_0` are all nonzero.  Put

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                      (24)
```

No square or cube root is chosen.  Exact reduction gives

```text
z_0^2-K_0/(9g(eta)) = -Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0 = u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0 = -Q/(eta^2K_0^2),
{A,z}|_0 = 3P/eta^2.                                  (25)
```

Thus the point lies on `Y` and the last two Hamiltonian components vanish.
The surface Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                     (26)
```

and `K_0!=0` kill the remaining component.  This is an actual critical
point, so `(2)` has no regular Darboux mate.

The theorem closes every `g` of degree at most two.  Degree at least three,
mixed `z^2h(e)` plus `r g(e)` corrections, and different arm profiles remain
open.  The companion named in the metadata checks `(5)--(26)` and the full
remainder.  Exhaustive scans of the divisibility system over `F_5`, `F_7`,
and `F_11` find no survivors, but these are hostile controls only and play
no role in the characteristic-zero proof.  Normal and optimized executions
LF-normalize exactly to the frozen 37-gate transcript.  **QED.**
