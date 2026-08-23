---
id: THM-3814
title: "Nodal RZ profile repairs cannot be Darboux pairs"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  On the c=1 cubic pseudo-plane, no Darboux pair
  with the normalized nodal arm jet can have all higher canonical terms
  confined to r g(e), z^2 f(e), and rz p(e) in the two outputs.  The top
  Wronskian first forces a 5/4 Kummer tower.  A nonzero root of its parameter
  has incompatible orders 4m-1 and 5m-1; if every root is zero, the origin
  has incompatible orders 3+4d and 3+5d.  The constant tower instead leaves
  the exact residual 1/4.  This closes the exact rz-profile ansatz, not
  arbitrary higher canonical coefficients and not the planar JC.
source: jc_zero_debt_lift / nodal rz Kummer-tower lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The deterministic companion has 44 active
  gates checking the Poisson signs, unique monic z-normal reduction, the
  arm, Wronskian, and r^2z^2 buckets, both Kummer valuation patterns, the
  asymmetric boundary, the complete constant-tower divisions and residual
  1/4, and both all-degree local leading coefficients and order gaps.  As an
  independent hostile control only, it also closes the linear tower through
  two arm remainders, a uniform rank-five compatibility, and an exact unit
  ideal.  Normal/-O/frozen/hash/documentation replay is required before
  audit promotion.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_nodal_rz_kummer_gate_thm3814.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_nodal_rz_kummer_gate_thm3814.out
script_sha256: 2daac8be5a100771e12e0e99ffa70a59f378697df84ec63288a25ab63ff5220a
output_sha256: 9d47c17bb34a7655e9559c1ccdab3da0678fa1afa4d2918161f4ca733eabebe2
semantic_sha256: 93f4511b995aec57406202ed542a539f8c2d357489ee4c87fd35822f1161b619
hash_basis: raw LF bytes
---

# THM-3814 -- the first mixed nodal profile is still too small

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Let `k` be an algebraically closed field of
characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary `f,g,h,kappa,p,q in k[e]`, define

```text
A=e^2-z/3+r g(e)+z^2 f(e)+rz p(e),
C=e^3-e-ez/2+r h(e)+z^2 kappa(e)+rz q(e).             (2)
```

Then

```text
{A,C} != 1.                                            (3)
```

Thus the normalized nodal arm profile cannot be repaired by adding only
arm-dependent coefficients in the canonical `r`, `z^2`, and `rz` slots.

## 1. Three exact buckets

The monic presentation

```text
B=k[r,e][z]/(z^3-r^2e-r)                              (4)
```

makes `(1,z,z^2)` a free `k[r,e]`-basis.  Reduction to `z`-degree at most
two therefore loses no quotient information.  In the unique normal form of
`{A,C}-1`, three coefficients are

```text
[z] = (36e^2f-24e kappa-12f+1)/2,                     (5)

[r^3z] = 15e(pq'-qp'),                                (6)

[r^2z^2]
 = -3(-4efq'+4e kappa p'-5ep kappa'
      +5eqf'+2fq-2kappa p).                           (7)
```

If `(3)` were false, `(5)` would give the arm identity

```text
(3e^2-1)f-2e kappa=-1/12,              hence f(0)=1/12. (8)
```

The domain `k[e]` lets us cancel the factor `e` in `(6)`, giving

```text
pq'-qp'=0.                                             (9)
```

## 2. The one-sided mixed branch is empty

Suppose first that `p=0`.  If `q=0`, THM-3812 already contradicts the
putative bracket identity.  If `q!=0`, equation `(7)` becomes

```text
5e q f'-4e f q'+2fq=0.                                (10)
```

Here `f!=0` by `(8)`.  Dividing in `k(e)` and using that its differential
constant field is `k` gives

```text
(q^4/(e^2f^5))'=0,
q^4=gamma e^2f^5                         (gamma in k*). (11)
```

Unique factorization at `e` and at every other irreducible factor forces

```text
q=alpha e^3v^5,              f=beta e^2v^4             (12)
```

with `alpha beta!=0`.  This makes `f(0)=0`, contradicting `(8)`.  Hence any
hypothetical pair must have

```text
p!=0.                                                  (13)
```

## 3. Every remaining pair enters one Kummer tower

Equation `(9)` now says

```text
q=lambda p                                             (14)
```

for some `lambda in k`.  Put

```text
W=kappa-lambda f.                                      (15)
```

The polynomial `W` is nonzero: otherwise `(8)` would make the nonconstant
quadratic `3e^2-2lambda e-1` times a polynomial equal a nonzero scalar.
Substitution in `(7)` gives

```text
4eW p'-5epW'-2pW=0.                                   (16)
```

Equivalently in `k(e)`,

```text
(p^4/(e^2W^5))'=0,
p^4=gamma e^2W^5                         (gamma in k*). (17)
```

If `a=ord_pi(p)` and `b=ord_pi(W)`, unique factorization gives

```text
4a=5b                         for pi!=e,
4a=2+5b                       for pi=e.                (18)
```

The nonnegative solutions are `(a,b)=(5t,4t)` away from `e` and
`(a,b)=(3+5t,2+4t)` at `e`.  Consequently there are nonzero constants
`alpha,beta` and a nonzero polynomial `v` such that

```text
p=alpha e^3v^5,
W=beta e^2v^4.                                        (19)
```

It remains to show that neither constant nor nonconstant `v` survives.

## 4. The constant tower leaves an exact quarter

If `v` is constant, absorb it into `alpha,beta`.  Put

```text
D=3e^2-2lambda e-1,
p=alpha e^3,                     W=beta e^2.           (20)
```

The arm law `(8)` becomes

```text
D f=2beta e^3-1/12.                                   (21)
```

The remainder on division by `D` is

```text
4beta lambda/9-1/12
 + e(8beta lambda^2/9+2beta/3).                       (22)
```

It must vanish, so

```text
4lambda^2+3=0,              beta=-lambda/4,
f=1/12-lambda e/6,
kappa=lambda/12+e/8-lambda e^2/4.                     (23)
```

Write `H=h-lambda g`.  The pure `r^3` and then pure `z^2` buckets give

```text
eH'-2H=(3e+4lambda)/(120alpha),
H=delta e^2-e/(40alpha)-lambda/(60alpha),              (24)

Dg=(1-2lambda e+48eH)/24.                             (25)
```

Modulo `4lambda^2+3`, the remainder in `(25)` is

```text
(160alpha delta lambda+15alpha-6)/(360alpha)
 - e lambda(5alpha+4)/(60alpha).                      (26)
```

Since `lambda!=0`, this forces

```text
alpha=-4/5,                 delta=3lambda/16,
g=(3lambda e-1)/24,
h=(9lambda e^2-3e-lambda)/48.                         (27)
```

Substitute `(23),(27)` into the full bracket and reduce by `(4)`.  The
coefficient of the pure monomial `r e^0` is

```text
[r z^0 e^0]({A,C}-1)=1/4                              (28)
```

in `k[lambda]/(4lambda^2+3)`.  This is nonzero in characteristic zero, so
constant `v` is impossible.

## 5. Every nonconstant tower has a local order mismatch

For general `v`, the pure `r^3` bucket is, after division by `-3`,

```text
4e^2(Wf'-fW')+3eH p'-5epH'+Hp=0,       H=h-lambda g.  (29)
```

Suppose first that `v` has a nonzero root `rho` of multiplicity `m>=1`.
Equation `(8)`, written using `(15)`, is

```text
D f=2eW-1/12.                                         (30)
```

At `rho`, the right side is `-1/12`.  Thus `(30)` itself forces

```text
D(rho)!=0,                     f(rho)!=0.              (31)
```

This point matters: no root of `D` is a missing denominator case.

Take a local parameter `x=e-rho` and write
`v=x^m u` with `u(rho)!=0`.  Since `W=beta e^2v^4`, the first term in
`(29)` has exact order

```text
ord_rho(4e^2(Wf'-fW'))=4m-1,                          (32)
```

with leading coefficient

```text
-16 beta m rho^4 f(rho)u(rho)^4 !=0.                  (33)
```

But `p=alpha e^3v^5`, so every term involving `H,p` in `(29)` has order at
least `5m-1`.  The strict inequality

```text
4m-1 < 5m-1                                            (34)
```

makes cancellation impossible.  Hence `v` has no nonzero root.

Because `k` is algebraically closed, a nonconstant polynomial whose only
root is zero has the form `v=c e^d` with `d>=1`.  Absorb `c` into
`alpha,beta`.  Then

```text
W=beta e^(2+4d),                 p=alpha e^(3+5d).     (35)
```

Equation `(30)` gives `f(0)=1/12`.  Therefore the forcing in `(29)` has
exact origin order and leading coefficient

```text
ord_0(4e^2(Wf'-fW'))=3+4d,
[e^(3+4d)]=-beta(2+4d)/3 !=0.                         (36)
```

Every `H,p` term has order at least `3+5d`.  Since

```text
3+4d < 3+5d                         for d>=1,          (37)
```

equation `(29)` is again impossible.  This eliminates every nonconstant
`v`; Section 4 eliminated the remaining constant case.  The contradiction
proves `(3)`.

## 6. Exact scope and controls

This theorem is deliberately a canonical-profile result.  It says that a
live nodal Darboux construction must leave `(2)` by adding another
`r,z`-monomial layer or by changing the arm jet itself.  It does not permit
higher coefficients to be silently discarded.  Also, as THM-3812 records,
`r in (r,z)^3` and `rz in (r,z)^4`; `rz` is the first omitted *mixed
canonical coefficient*, not an `I^3` term.

The exact companion derives the arbitrary-function normal form and checks
`(5)--(7),(16),(24)--(29),(32)--(37)`.  As a logically redundant hostile
control, it also normalizes `v=e-t`, derives both arm remainders, proves the
remaining coefficient operator has uniform rank five (constant minor
`375000`), handles `t=0` separately from the generic left-null row, and
computes the exact unit ideal of the arm and compatibility packet.  The
all-degree proof is the local valuation argument in Section 5, not that
linear specialization.  No finite-field or bounded-degree inference is
used.  **QED.**
