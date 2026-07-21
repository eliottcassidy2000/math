---
id: THM-2023
title: "THE GMC(2) B-DOMINANT SHARP-BOUNDARY HYPER-BESSEL FUNCTIONS ARE LAGUERRE--POLYA: every zero of Phi_(p,q)(x)=sum x^k/((pk)!(qk)!) is strictly negative real"
status: >
  PROVED. Gauss multiplication identifies Phi_(p,q) exactly with a positive
  argument-rescaling of 0F_(p+q-1) whose denominator parameters are
  {1,1/p,...,(p-1)/p,1/q,...,(q-1)/q}. They are all positive, so the
  Baricz--Singh/Laguerre theorem for positive-parameter 0FQ functions gives
  only real-negative zeros; its order is below one, hence Phi is Laguerre--Polya
  type I. This confirms HYP-8775 and proves that THM-2017's rd-e=r boundary is
  NC2-clear whenever xi is off the negative real axis. It does not treat the
  opposite e-rd=r function Psi_r=sum y^j/(rj)!, nor the discrete negative-root
  parameters themselves.
source: codex-2026-07-21-gmc2-hyperbessel-lp
depends_on: [THM-2017]
related: [HYP-8775, HYP-8766, THM-2021]
script: 04-computation/gmc2_hyperbessel_laguerre_polya_thm2023.py
output: 05-knowledge/results/gmc2_hyperbessel_laguerre_polya_thm2023.out
external:
  - "Arpad Baricz and Sanjeev Singh, Zeros of some special entire functions, Proc. Amer. Math. Soc. 146 (2018), 2207--2216, doi:10.1090/proc/13927, arXiv:1702.00626"
---

# THM-2023 -- the sharp-boundary hyper-Bessel function is Laguerre--Polya

> **Global closure context.** THM-2022 now proves full NC2/GMC(2) by a
> finite-place lowest-face argument. THM-2023 remains a sharper analytic theorem
> locating the possible leading-limit boundary zeros; its negative-root
> transseries residual is no longer a logical NC2 gap.

## 1. Statement

For positive integers `p,q`, define

```text
Phi_(p,q)(x) = sum_(k>=0) x^k/((pk)!(qk)!).           (1)
```

Then every zero of `Phi_(p,q)` is real and strictly negative. More precisely,
`Phi_(p,q)` belongs to the type-I Laguerre--Polya class.

In THM-2017 the primitive charges are `p0,q0`; (1) is exactly the limiting
function at the `b`-dominant sharp degree boundary `rd-e=r`. Thus the boundary
limit can vanish only when its parameter `xi` is negative real.

## 2. Exact hypergeometric identification

For an integer `p>=1`, Gauss multiplication gives the elementary Pochhammer
identity

```text
(pk)! = p^(pk) k! product_(j=1)^(p-1) (j/p)_k.       (2)
```

Apply (2) to `p` and `q`. Let the positive parameter multiset be

```text
B_(p,q) = {1}
          union {j/p : 1<=j<p}
          union {j/q : 1<=j<q}.                      (3)
```

It has `p+q-1` entries. By the definition of the generalized hypergeometric
series,

```text
Phi_(p,q)(x)
 = 0F_(p+q-1)( - ; B_(p,q) ; x/(p^p q^q)).          (4)
```

Indeed the `k`th denominator on the right is

```text
k! (1)_k
 product_(j=1)^(p-1)(j/p)_k
 product_(j=1)^(q-1)(j/q)_k
 (p^p q^q)^k
 = (pk)!(qk)!.
```

No asymptotic or numerical zero matching enters this step.

## 3. Zero theorem

Baricz--Singh prove, by Laguerre's classical lemma (Theorem 2 and the paragraph
immediately following it), that

> if every denominator parameter is positive, then
> `0F_Q(-;b_1,...,b_Q;z)` has only real negative zeros.

Every entry of (3) is positive and the rescaling `p^p q^q` is positive.
Therefore (4) proves that every zero of (1) is real and negative. It is strictly
negative because `Phi_(p,q)(0)=1`.

Equivalently, in their normalized hyper-Bessel notation take
`alpha_i=B_i-1`. Positivity of `B_i` gives `alpha_i>-1`, exactly the hypothesis
of their Theorem 2; the sign change there is the sign appearing in the
positive-coefficient `0F_Q` used here.

For completeness, `0F_Q` has order `1/(Q+1)<1` here (`Q=p+q-1>=1`). Hence
its real-negative zero product has genus zero; together with positive value at
the origin this places (1) in the type-I Laguerre--Polya class. This is the
precise theorem numerically anticipated by HYP-8775.

## 4. GMC(2) consequence

At THM-2017's boundary `rd-e=r`,

```text
E[P^m]/L(b^m) -> Phi_(p0,q0)(xi),
xi = alpha/(beta^r d^r).                             (5)
```

If `xi` is not negative real, THM-2023 makes the limit in (5) nonzero. Therefore
`E[P^m]` is nonzero for all sufficiently large `m`, and NC2 holds on that part
of the boundary. In particular this proves the off-axis conclusion for every
complex primitive-charge pair, not just the `I_0` base `p0=q0=1`.

Positive `xi` was already immediate from the positive coefficients; the new
content is the whole complex plane minus the negative real axis and the exact
location of every possible exceptional zero.

## 5. Scope: there are two sharp boundaries

THM-2017's opposite boundary `e-rd=r` uses

```text
Psi_r(y)=sum_(j>=0)y^j/(rj)!,                         (6)
```

a Mittag--Leffler/root-of-unity-filter function. Formula (6) is **not** (1)
except in special low-order cases, and for `r>2` its zero geometry is not the
positive-parameter hyper-Bessel theorem above. Thus the slogan "the sharp
boundary is Laguerre--Polya" must be read as the `Phi`, or `b`-dominant,
boundary only. The negative-real zeros of `Phi` themselves also remain possible
leading-limit cancellations; THM-2017 removes them in its symmetric monomial
model by the next `1/m` derivative term, not universally.

## 6. Verification and methodology

The exact script checks (4), coefficient by coefficient through `k=30`, for the
eight parameter pairs used in HYP-8775:

```text
(1,1),(1,2),(1,3),(2,2),(2,3),(1,4),(3,4),(3,5).
```

It also verifies positivity and the exact count of the denominator parameters.
The zero-location theorem is cited, not replaced by finite numerics.

Tournament Analysis is not used here. Candidate vertices (zeros, denominator
parameters, the two boundary sides, and proof obligations) do not admit a
pairwise orientation preserving the unary Laguerre--Polya predicate. Imposing a
tie order would add structure and erase the exact hypergeometric identification.
This is a case where the honest carrier is the parameter multiset (3), not a
tournament quotient.

## 7. What changed

HYP-8775 treated general all-negative zeros as strong numerical evidence and
searched for a new Pólya--Schur multiplier-sequence proof. Equation (4) shows
that the exact statement is already a direct instance of the classical
positive-parameter `0F_Q` theorem. The numerical pattern is promoted to
THM-2023; only the exceptional negative-root transseries and the opposite
`Psi_r` boundary remain.
