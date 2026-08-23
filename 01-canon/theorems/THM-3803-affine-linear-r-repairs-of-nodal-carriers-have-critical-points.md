---
id: THM-3803
title: "Affine-linear R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  On the c=1
  cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+r(b_0+b_1e), with (b_0,b_1) nonzero as a pair, has a critical
  point.  In the genuinely binomial case the degree-eleven residual
  resultant must have a root away from e(b_0+b_1e)=0: boundary-only support
  would force origin order eight, whereas its first three coefficients force
  origin order at most two.  Quadratic and higher nonmonomial repairs remain
  open.  This is a carrier obstruction, not a global JC(2) theorem.
source: root / affine-linear universal-resultant boundary lane, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The Hamiltonian reduction, full residual
  coefficient table, degree and boundary-support argument, homogeneous
  finite-root seam, excluded denominator values, explicit point recovery,
  and Casimir implication were rederived symbolically.  The deterministic
  companion has 279 active gates; normal and optimized runs byte-match the
  frozen transcript.  Independent hostile audit remains due.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
related:
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_affine_linear_r_repair_thm3803.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_affine_linear_r_repair_thm3803.out
script_sha256: 209d89b52a5a572812939cdd45fef27039e490d824b0416394a29c8a22ccba4c
output_sha256: ec0002d4d25846e09117f5fc42210ec0f6b51b56c1453c9f6df80c3f54bdf98c
semantic_sha256: be48a034fe61375ec1c1617aa6a2bf469d96f9cb3f8bdf5c3f7d3f41070456e6
hash_basis: raw LF bytes
---

# THM-3803 -- every affine-linear r-repair remains critical

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For every `(b_0,b_1) in k^2\{(0,0)}`, the regular function

```text
A=e^2-z/3+r g(e),                 g(e)=b_0+b_1e        (2)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.
The cases `b_1=0` and `b_0=0` are respectively the `m=0` and `m=1` rows of
THM-3799.  It remains to prove the genuinely binomial case
`b_0 b_1 !=0`.

## 1. Universal critical equations

Put `u=re` and `K=1+2u`.  Directly from `(1)`, criticality compresses to

```text
P=g u^2-K(2e^3+u e b_1),
Q=e^2K^3-729g^3u^2(1+u)^2.                              (3)
```

Indeed `P/e^2={A,z}/3` after `r=u/e`, while `Q` is the compatibility of

```text
z^2=K/(9g),                    z^3=u(1+u)/e.             (4)
```

Exact elimination in `u` gives

```text
Res_u(P,Q)=e^4 g(e)^3 H(e),              deg H=11.       (5)
```

Writing `H=sum_(j=0)^11 h_j e^j`, its complete coefficient table is

```text
h_0 = b_0(1-729b_0b_1^2)
h_1 = -2916b_0^3-b_1
h_2 = 2916b_0^2(729b_0^3+b_1)
h_3 = 2916b_0b_1(2187b_0^3-b_1)
h_4 = 8748b_0^2(729b_0b_1^2-2)
h_5 = 8748b_0(972b_0^3+243b_0b_1^3+2b_1)
h_6 = 25509168b_0^3b_1
h_7 = 11664b_0(2187b_0b_1^2-2)
h_8 = 11664(729b_0^3+729b_0b_1^3+2b_1)
h_9 = 25509168b_0^2b_1
h_10= 25509168b_0b_1^2
h_11= 8503056b_1^3.                                      (6)
```

In particular `h_11!=0`, and

```text
h_10/h_11=3b_0/b_1.                                    (7)
```

## 2. The residual resultant must leave the boundary

Suppose, for contradiction, that every root of `H` lies on the forbidden
boundary `e g(e)=0`.  Since `b_0b_1!=0`, the two boundary points are distinct.
Unique factorization over the algebraically closed field then gives

```text
H=mu e^a(b_0+b_1e)^(11-a),              0<=a<=11,      (8)
```

where `a=ord_(e=0) H`.  Comparing the top two coefficients of `(8)` with
`(7)` yields

```text
(11-a)b_0/b_1=3b_0/b_1,                  so a=8.        (9)
```

But `(6)` makes `ord_0(H)<=2`.  If `h_0!=0` or `h_1!=0` this is immediate.
If both vanish, `b_0!=0` and `h_1=0` give

```text
b_1=-2916b_0^3,
h_2=-2916*2187*b_0^5 !=0.                              (10)
```

This contradicts `(9)`.  Hence `H` has a root `eta` satisfying

```text
eta(b_0+b_1eta)!=0.                                    (11)
```

## 3. The off-boundary root is an actual critical point

At `e=eta`, equation `(5)` says that `P,Q` have a common projective
`u`-root.  It is finite: the homogenized value of `Q` at infinity is

```text
LC_u(Q)=-729g(eta)^3 !=0,                               (12)
```

even if the quadratic coefficient `b_0-b_1eta` of `P` happens to vanish.
Let `u_0` be a common finite root.  The three dangerous values are excluded
directly:

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16.                            (13)
```

Thus `eta`, `g(eta)`, `u_0`, `1+u_0`, and `K_0=1+2u_0` are all nonzero.
Define, without choosing a square or cube root,

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                        (14)
```

Straight reduction gives

```text
z_0^2-K_0/(9g(eta)) = -Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0 = u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0 = -Q/(eta^2K_0^2),
{A,z}|_0 = 3P/eta^2.                                   (15)
```

The common-root equations therefore put `(r_0,z_0,eta)` on `Y` and kill
the last two Hamiltonian components.  Finally the surface Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                      (16)
```

and `K_0!=0` kill `{A,r}`.  Since `(1)` is smooth and symplectic, this is an
actual critical point of `A`, excluding every Darboux mate.

The theorem closes all affine-linear `g`, not quadratic or higher
nonmonomial repairs, mixed `r`/`z^2` corrections, or a different nodal arm
profile.  The deterministic companion verifies `(3)--(16)`, the entire
coefficient table, and 81 exact specialized controls with 279 active gates.
Normal and optimized executions byte-match the frozen transcript.
**QED, conditional only on independent hostile audit.**
