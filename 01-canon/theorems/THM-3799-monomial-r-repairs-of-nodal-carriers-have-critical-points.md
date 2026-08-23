---
id: THM-3799
title: "Monomial R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  On the c=1
  cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+lambda r e^m has a critical point for every integer m>=0 and
  lambda!=0.  Criticality reduces to two polynomials in e and u=re; their
  exact u-resultant has a residual factor J_m with nonzero constant and
  positive degree.  Every residual root survives all denominator and
  degree-drop boundaries and reconstructs an actual source critical point.
  Hence no monomial r-repair admits a Darboux mate.  Nonmonomial polynomial
  r g(e) and mixed corrections remain open.
source: jc_zero_debt_lift / nodal monomial-r resultant lane, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof rederives all Hamiltonian
  components, the P/Q compression, universal resultant, closed J_m formula,
  fixed leading coefficients, excluded u=0,-1,-1/2 strata, rational point
  reconstruction, and the remaining Hamiltonian via the surface Casimir.
  The m=3, lambda=1 collision is isolated explicitly.  The deterministic
  companion has 2,336 active gates; normal and optimized runs byte-match
  the frozen transcript.  Independent hostile audit remains due.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
script: 04-computation/jc2_cubic_pseudoplane_monomial_r_repair_thm3799.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_monomial_r_repair_thm3799.out
script_sha256: 7d08fa698e3379397807cb65e87d33a37aeaa0907fca1424917a06dd04c02abc
output_sha256: 3708cf945ad9f50711c82dce20d9b2848d0e157eb4df5db7d38449d78b876fdb
semantic_sha256: ca95b96e54ab327ec3b4c52837b5b074a0134df23e135e7489af23efaa672ebf
hash_basis: raw LF bytes
---

# THM-3799 -- every monomial r-repair remains critical

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),                         (1)
```

with Poisson packet

```text
{r,z}=3r^2,             {r,e}=9z^2,             {z,e}=3+6re.         (2)
```

For every integer `m>=0` and every `lambda in k*`, the regular function

```text
A=e^2-z/3+lambda r e^m                                (3)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

The correction in `(3)` is invisible on the arm `r=z=0`, and the
first-normal coefficient remains the forced nodal value `-1/3`.  Thus this
is the first genuinely `r`-coupled lane left open by THM-3795, not a change
of the boundary immersion.

## 1. Compress criticality to `(e,u=re)`

It is useful to carry a general polynomial `g(e)` for one step.  For

```text
A_g=e^2-z/3+r g(e),                  K=1+2re,          (4)
```

direct use of `(2)` gives

```text
{A_g,r}=r^2-9z^2(2e+r g'),
{A_g,z}=3g r^2-3K(2e+r g'),
{A_g,e}=9g z^2-K.                                      (5)
```

Put `u=re` and regard `g=g(e)`, `g'=g'(e)`.  Define

```text
P=g u^2-(1+2u)(2e^3+u e g'),                           (6)
Q=e^2(1+2u)^3-729g^3u^2(1+u)^2.                       (7)
```

The polynomial `P/e^2` is exactly `{A_g,z}/3` after `r=u/e`.  The other
equation expresses compatibility of

```text
z^2=(1+2u)/(9g),             z^3=u(1+u)/e,             (8)
```

the first identity coming from `{A_g,e}=0` and the second from the surface
equation.  Section 3 proves the converse constructively, after every
denominator has been audited.

## 2. Exact monomial resultant

Now set

```text
g=lambda e^m.                                          (9)
```

Exact elimination in `u` gives

```text
Res_u(P,Q)=lambda^4 e^(4m+4) J_m(e),                  (10)
```

where, for `m>=1`,

```text
J_m(e)=
 (1-2m)
 +23328(2m-1)e^7
 -17496(2m-1)(m-1)lambda e^(m+4)
 +2916(m-1)(3m^2-3m+1)lambda^2 e^(2m+1)
 +8503056lambda^2 e^(2m+8)
 -729m^2(m-1)^2lambda^3 e^(3m-2)
 -8503056(m-1)lambda^3 e^(3m+5)
 +2125764(m-1)^2lambda^4 e^(4m+2).                   (11)
```

For `m=0`, the same formula is used with the formally negative-exponent
term, whose coefficient contains `m^2`, omitted.  Thus `(11)` is a genuine
polynomial in every case.  One way to audit `(10)--(11)` without guessing
the pattern is to first keep independent symbols `G=g(e)` and `D=g'(e)`:
the universal resultant factors as `G^3e^4H(e,G,D)`, and substituting
`G=lambda e^m`, `D=m lambda e^(m-1)` gives

```text
H=lambda e^m J_m.                                     (12)
```

The residual factor always has a nonzero root.  First,

```text
J_m(0)=1-2m!=0.                                       (13)
```

For `m!=3`, no other exponent in `(11)` equals seven, so the coefficient

```text
[e^7]J_m=23328(2m-1)                                  (14)
```

is nonzero.  Hence `J_m` is nonconstant.  At the sole collision `m=3`,
the coefficients at degrees fourteen and seven are

```text
[e^14]J_3=8503056lambda^2(lambda-1)^2,                 (15)
```

and, when `lambda=1`,

```text
J_3(e)=26244e^7-5.                                    (16)
```

Thus `J_m` is nonconstant even on the exceptional seam.  Algebraic
closedness supplies a root `eta`, and `(13)` makes `eta!=0`.

## 3. Every residual root lifts to a source critical point

At `e=eta`, the two leading `u`-coefficients are

```text
LC_u(P)=lambda eta^m(1-2m),
LC_u(Q)=-729lambda^3eta^(3m).                         (17)
```

Both are nonzero.  Therefore the vanishing resultant in `(10)` is not a
degree-drop or common-root-at-infinity artifact: `P(eta,u)` and `Q(eta,u)`
have a common finite root `u_0 in k`.

None of the dangerous boundary values can occur.  Indeed

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729(lambda eta^m)^3/16,                  (18)
```

all nonzero.  Hence

```text
eta * g(eta) * u_0 * (1+u_0) * (1+2u_0) !=0.         (19)
```

Put `g_0=lambda eta^m`, `K_0=1+2u_0`, and define

```text
r_0=u_0/eta,
z_0=9g_0u_0(1+u_0)/(eta K_0).                        (20)
```

These are actual field elements, with no root choice left implicit.  Direct
reduction gives the identities

```text
z_0^2-K_0/(9g_0)=-Q/(9g_0 eta^2K_0^2),
r_0^2eta-z_0^3+r_0=u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0=-Q/(eta^2K_0^2),
{A,z}|_0=3P/eta^2.                                   (21)
```

Thus the reconstructed point lies on `Y` and the last two Hamiltonian
components vanish.  Finally the surface equation is a Casimir:

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0.                   (22)
```

Since `K_0!=0`, `(22)` kills the remaining component `{A,r}`.  The surface
is smooth and symplectic, so `(r_0,z_0,eta)` is a critical point of `(3)`.
At such a point `{A,C}=0` for every regular `C`; hence no Darboux mate can
satisfy `{A,C}=1`.

The theorem deliberately stops at monomial `g`.  A nonmonomial polynomial
can change both the residual resultant and the leading coefficient
`g-2eg'`; mixtures with `z^2h(e)` add another collision grammar.  Those are
the live repair coordinates, together with changing the arm profile.  The
deterministic companion named in the metadata verifies `(5)--(22)` and the
closed formula through `m=256` with 2,336 active assertion-free gates.
Normal and optimized executions byte-match the frozen transcript.
**QED, conditional only on independent hostile audit.**
