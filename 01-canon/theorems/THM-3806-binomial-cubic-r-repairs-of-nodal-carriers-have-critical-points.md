---
id: THM-3806
title: "Binomial cubic R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On the
  c=1 cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+r(mu e^j+lambda e^3), j in {0,1,2}, mu*lambda!=0, has a
  critical point.  Boundary-only support for the degree-seventeen residual
  resultant would force H to divide e g H', but exact Euclidean remainders
  contradict this for all three lower exponents.  The genuine degree-drop
  seam lambda=1 and the surviving j=2 witness seam lambda=2 are recomputed
  in their specialized coefficient rings.  Hence no such binomial cubic
  r-repair has a regular Darboux mate.  Cubics with at least three monomials,
  higher degree, and mixed corrections remain open.
source: jc_zero_debt_lift / cubic-binomial boundary-divisibility lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn-boundary, 2026-08-23).  The exact
  companion has 46 active gates checking the
  universal P/Q resultant, all generic and exceptional degrees and leading
  coefficients, localized Euclidean divisions, decisive remainder
  coefficients, Hamiltonian signs, denominator exclusions, and explicit
  source reconstruction.  Normal and optimized runs byte-match the frozen
  transcript.  The audit re-derived the boundary-support divisibility with
  multiplicities, checked every localization and exceptional specialization,
  verified that the quartic Q excludes a root at infinity, replayed the
  source/critical reconstruction, ran normal and optimized companions against
  the frozen output, and matched the recorded raw-byte hashes.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points
related:
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
  - THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points
script: 04-computation/jc2_cubic_pseudoplane_binomial_cubic_r_repair_thm3806.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_binomial_cubic_r_repair_thm3806.out
script_sha256: bb98c4a4ccb4d3a50f4e716b6085201b2b9cc2601fb2415350c3270deef756e0
output_sha256: 8cde2e683a25152e12e9c65e2aeafaf1764c360939948e7a0e757c08acaab95e
semantic_sha256: 1127c20eb15e977d79c605186a7d17b848578963dba9ee6252219496cb4f77c7
hash_basis: raw LF bytes
---

# THM-3806 -- every binomial cubic r-repair remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),                         (1)
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.      (2)
```

For every

```text
j in {0,1,2},             mu,lambda in k*,            (3)
g(e)=mu e^j+lambda e^3,                                (4)
```

the regular function

```text
A=e^2-z/3+r g(e)                                      (5)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

The lower-degree and monomial faces are already closed by THM-3805 and
THM-3799.  The content here is that a nonzero cubic term cannot repair any
one nonzero lower monomial, including the degree-collision seams where the
generic residual resultant changes degree or its first witness vanishes.

## 1. The universal two-equation compression

Put

```text
u=re,                           K=1+2u,                (6)
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (7)
```

The Hamiltonian components of `(5)` are

```text
{A,r}=r^2-9z^2(2e+r g'),
{A,z}=3g r^2-3(1+2re)(2e+r g'),
{A,e}=9g z^2-(1+2re).                                (8)
```

After `r=u/e`, the equation `P=0` is exactly `{A,z}=0` up to the
nonzero factor `3/e^2`.  The equation `Q=0` is the compatibility of

```text
z^2=K/(9g),                    z^3=u(1+u)/e.           (9)
```

For an elimination that remains valid before specializing `g`, keep
independent symbols `G,D` in place of `g,g'`.  Exact elimination in `u`
gives

```text
Res_u(P,Q)=G^3e^4 H_univ(e,G,D).                     (10)
```

For each `j` in `(3)`, substitute `(G,D)=(g,g')` and call the resulting
residual polynomial `H_j(e)`.  Thus

```text
Res_u(P,Q)=g(e)^3 e^4 H_j(e).                        (11)
```

## 2. Boundary-only roots force logarithmic divisibility

Suppose every root of a nonconstant `H_j` lies on the forbidden boundary

```text
V(e g(e)).                                            (12)
```

This implication has to retain multiplicities.  Factor over `k` as

```text
H_j=a product_alpha (e-alpha)^(n_alpha).              (13)
```

Every `alpha` occurring in `(13)` is a zero of `e g`.  Therefore

```text
e g H_j'/H_j
  =sum_alpha n_alpha e g/(e-alpha) in k[e],           (14)
```

and boundary-only support necessarily implies

```text
H_j divides e g H_j'.                                 (15)
```

Only this necessary direction is used; `(14)` is valid with repeated roots
of either `H_j` or `g`.

## 3. The generic `lambda!=1` division

First assume `lambda!=1`.  Uniformly for `j=0,1,2`, exact expansion of
`(10)` gives

```text
deg H_j=17,
LC(H_j)=8503056 lambda^3(lambda-1)^2.                 (16)
```

Perform Euclidean division of `e g H_j'` by `H_j`.  The identities live in
the localized parameter ring

```text
Q[mu,lambda,mu^-1,lambda^-1,(lambda-1)^-1][e].       (17)
```

More sharply, common denominators for quotient and remainder may be taken
to be `lambda-1`, `1`, and `4(lambda-1)^3` for `j=0,1,2`, respectively.
Thus these are polynomial identities after every specialization allowed by
`mu lambda(lambda-1)!=0`; no generic parameter argument is being applied on
a pole of the quotient.

For `j=0` and `j=1`, the remainder has the same decisive coefficient:

```text
[e^6] rem_{H_j}(e g H_j')=70 lambda^2!=0.            (18)
```

Hence `(15)` is impossible in these two rows.  For `j=2`, the highest
decisive coefficient is

```text
[e^16] rem_{H_2}(e g H_2')
 =-1062882 mu^4 lambda^3(lambda-2)/(lambda-1)^2.      (19)
```

This is nonzero unless `lambda=2`.  At `lambda=2`, recompute the division
after specialization.  The specialized polynomial still has degree
seventeen and leading coefficient `68024448`, while

```text
[e^4] rem_{H_2}(e g H_2')=30mu^2!=0.                 (20)
```

So `(15)` fails for every `lambda!=1`.

## 4. The degree-drop seam `lambda=1`

Equation `(16)` loses its leading term at `lambda=1`, and the generic
quotients for `j=0,2` have poles there.  It would therefore be invalid to
substitute `lambda=1` into the generic division.  Instead substitute into
`H_j` first and repeat Euclidean division in `Q(mu)[e]`.  The three exact
rows are

```text
 j    deg H_j       LC(H_j)             nonzero remainder coefficient
 0       11         2125764mu^2          [e^2] = 361/(81mu)
 1       10         26244                [e^6] = 35
 2       15         2125764mu^2          [e^4] = 30mu^2.              (21)
```

For `j=0`, a common denominator for the recomputed quotient and remainder
is `531441mu^6`; for `j=1` it is `3`, and for `j=2` it is `1`.  All are
invertible under `(3)` in characteristic zero.  Thus every specialized
`H_j` is nonconstant and `(15)` again fails.  This separate computation is
essential: exceptional specializations are recomputed in their actual
coefficient rings, not inferred by substituting into a localized generic
quotient.

Combining Sections 3 and 4, `H_j` is always a nonconstant polynomial and
cannot have all its roots on `(12)`.  Algebraic closedness therefore gives
a root `eta` such that

```text
H_j(eta)=0,                      eta g(eta)!=0.        (22)
```

## 5. Every off-boundary residual root is an actual critical point

At `e=eta`, equations `(11)` and `(22)` make the `u`-resultant vanish.  The
leading coefficient of the quartic `Q` is

```text
LC_u(Q)=-729g(eta)^3!=0.                              (23)
```

Consequently the homogenized pair cannot acquire a common root only at
infinity, even if the quadratic leading coefficient of `P` vanishes.
There is a common finite root `u_0` of `P(eta,u)` and `Q(eta,u)`.  Directly,

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16,                         (24)
```

so `u_0`, `1+u_0`, and `K_0=1+2u_0` are all nonzero.  Define, without any
root choice,

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                    (25)
```

The following exact identities hold before setting `P=Q=0`:

```text
z_0^2-K_0/(9g(eta))=-Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0=u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0=-Q/(eta^2K_0^2),
{A,z}|_0=3P/eta^2.                                  (26)
```

Thus `(r_0,z_0,eta)` lies on `Y` and the last two Hamiltonian components
vanish.  The defining equation of `Y` is a Casimir, giving

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0.                  (27)
```

Since `K_0!=0`, `(27)` kills `{A,r}` as well.  This is an actual critical
point of `(5)`, so a regular function `B` with `{A,B}=1` cannot exist.

The theorem closes precisely the binomial cubics consisting of a nonzero
cubic term and one nonzero lower monomial.  Cubics with at least three
monomials, `deg g>=4`, simultaneous `z^2h(e)+r g(e)` corrections, other
arm profiles, and rational mates with poles remain open.  The exact
companion named in the metadata checks `(7)--(27)`, every generic and
exceptional division, the allowed denominator rings, and all reconstruction
identities.  Normal and optimized executions byte-match the frozen 46-gate
transcript.  **QED.**
