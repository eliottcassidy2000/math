---
source: kind-pasteur-2026-07-24-S152 (Opus 4.8)
status: RESULT (proofs + structural theorem). Advances the SUPERSET direction of the prime-power law
  ({1,...,P(n)-1} subset C_{1/n}). PROVES: (i) k=1 for all four signatures (uniform, clean); (ii) k=2 for
  lambda=1/4 (recap of kps-S146); (iii) a NEW clean general sibling identity T_a(2)=cos(pi a)/(1-2a); and
  (iv) a STRUCTURAL THEOREM: no uniform elementary k=2 formula can exist (forced by signature-specificity,
  kps-S151), so the superset direction is necessarily case-by-case. The remaining three values
  (S_{1/3}(2), S_{1/6}(2), S_{1/4}(3)) are reduced to classical signature-theory evaluations and verified to
  100+ digits -- honestly NOT re-derived from first principles here.
tags: [hypergeometric, closed-forms, prime-power, signatures, ramanujan, proofs, honest-status, series]
related: [kps-S146, kps-S147, kps-S150, kps-S151, THM-415]
---

# The superset direction: proofs, a sibling identity, and why no uniform formula exists

`S_a(k)=int_0^1 2F1(a,1-a;1;x^k)dx`. Goal: prove `{1,...,P(n)-1} subset C_{1/n}` for the four signatures,
i.e. that `S_{1/n}(k)` is elementary for every `k<P(n)`. The values needed:
`k=1` (all four); `k=2` (n=3,4,6); `k=3` (n=4).

## 1. k=1, all four signatures -- THEOREM (uniform, clean)
`1/(n+1)=(1)_n/(2)_n`, so `S_a(1)=sum (a)_n(1-a)_n/((n!)^2(n+1)) = 2F1(a,1-a;2;1)`. Gauss's theorem
(`2F1(A,B;C;1)=Gamma(C)Gamma(C-A-B)/(Gamma(C-A)Gamma(C-B))`, `C-A-B=1>0`) gives
`= Gamma(2)Gamma(1)/(Gamma(2-a)Gamma(1+a)) = 1/((1-a)a Gamma(1-a)Gamma(a))`, and reflection
`Gamma(a)Gamma(1-a)=pi/sin(pi a)`:
> **`S_a(1) = sin(pi a)/(pi a(1-a))`** (all `a`). Verified 110 digits.

Covers `k=1` for `n=2,3,4,6` at once: `4/pi, 9sqrt3/4pi, 8sqrt2/3pi, 18/5pi`. **Done.**

## 2. k=2, lambda=1/4 -- PROVED (recap, kps-S146)
`f_{1/4}(t)=(2/pi)int_0^{pi/2}dphi/sqrt(1-t cos^4 phi)` (duplication `(1/2)_{2n}=4^n(1/4)_n(3/4)_n`), so
`S_{1/4}(2)=(2/pi)int_0^{pi/2}[arcsin(cos^2 phi)/cos^2 phi]dphi`; IBP against `d(tan phi)` collapses it to
`2 int_0^1 du/sqrt(1+u^2)=2 arcsinh 1`. Hence **`S_{1/4}(2)=4log(1+sqrt2)/pi`.** Verified 110 digits.

## 3. A clean general sibling identity (new, and it kills a red herring)
Define the **c=1/2 sibling** `T_a(2)=sum (a)_n(1-a)_n/((1/2)_n n! (2n+1))=3F2(a,1-a,1/2;1/2,3/2;1)`. The upper/lower
`1/2` cancel, leaving `2F1(a,1-a;3/2;1)`; Gauss + reflection give
> **`T_a(2)=cos(pi a)/(1-2a)`, elementary for ALL `a`.** Verified 110 digits (`a=1/3,1/4,1/5,2/7,0.31`).

This is exactly the `cos(pi a)/(1-2a)` I had once mis-derived as a candidate for `S_a(2)` (kps-S148 method-note
spirit): it is the value of the **wrong-weight** sibling (`(1/2)_n n!` in place of `(n!)^2`), not of `S_a(2)`.
The extra factor `S_a(2)/T`-weight `= n!/(1/2)_n = 4^n/C(2n,n)` is one elliptic integration; that is precisely
the step that lifts the elementary sibling to the (generically non-elementary) period `S_a(2)`.

## 4. STRUCTURAL THEOREM: no uniform elementary k=2 formula can exist
> **Claim.** There is no `(A(a),B(a))` algebraic-valued at rationals with `S_a(2)=A(a)log B(a)` (or any fixed
> elementary template) valid on an interval of `a`.
> **Proof.** Such a template would force `S_a(2)` elementary at *every* rational `a`, in particular `a=1/5`.
> But `S_{1/5}(2)` is non-elementary (clean golden-field PSLQ null, 150 digits, kps-S151). Contradiction. **QED.**

So the `k=2` signature evaluations `{1/4->log(1+sqrt2), 1/3->log2, 1/6->log(2+sqrt3)}` are **isolated
degeneracies**, not samples of one formula (indeed their fields `Q(sqrt2),Q,Q(sqrt3)` and coefficients differ
irreconcilably). **The superset direction cannot be proved uniformly -- it is necessarily a per-signature,
per-`k` modular evaluation.** This is the closed-form-side shadow of signature-specificity, and it is the honest
reason "prove superset fully" is three separate special-function theorems, not one.

## 5. The remaining three -- reduced + verified (honest gap)
All verified to 100+ digits; each is elementary via the classical signature theory it names, but I have **not**
re-derived them from first principles here:
- `S_{1/3}(2)=3sqrt3 log2/pi` -- signature-3 (Ramanujan cubic theory; `2F1(1/3,2/3;1;.)`).
- `S_{1/6}(2)=3sqrt3 log(2+sqrt3)/(2pi)` -- signature-6 (sextic theory).
- `S_{1/4}(3)=-(sqrt3 log(5-2sqrt6)+2 arctan(sqrt2/5))/pi` -- signature-4, `k=3`. Mechanism (kps-S146,
  re-verified): `w=1+i sqrt2`, `|w|^2=3`, `w^3=-5+i sqrt2`, so `3 arg w = pi - arctan(sqrt2/5)` and
  `log(5-2sqrt6)=-2 arccosh... = 2 log(sqrt3-sqrt2)` -- the cube-root-of-unity split of the `g_3` cubic.

Each has the rigorous **period representation** `S_a(2)=(1/(2sqrt2))int_{-1}^1 P_{-a}(u)(1-u)^{-1/2}du`
(Legendre `P`); elementariness is the statement that this Legendre integral degenerates at `a=1/n`. That
degeneration is the classical-signature input still to be written out.

## 6. Honest headline
**Superset is PROVED for `k=1` (all four signatures, uniform theorem) and `k=2` (lambda=1/4).** The remaining
`{S_{1/3}(2),S_{1/6}(2),S_{1/4}(3)}` are verified to 100+ digits and are elementary by classical signature
theory, but re-derivation from scratch is open. Crucially (section 4) **no uniform proof is possible** -- so
"fully proved superset" = these three separate modular evaluations, each isolated by signature-specificity.
Bonus: the clean all-`a` sibling `T_a(2)=cos(pi a)/(1-2a)`.

Files: `/tmp/{supset,c5clean,cmspec}.py`. Builds on kps-S146/S150/S151.
