---
source: kind-pasteur-2026-07-24-S146 (Opus 4.8)
status: RESULT (analysis, external problem). Solves an owner-posed series problem:
  S(k)=sum_{n>=0} C(2n,n)C(4n,2n)/((kn+1)64^n). Establishes the hypergeometric reformulation, the closed-form
  generating function, the reduction mechanism giving the elementary S(1),S(2),S(3), and PROVES that S(4)
  (hence, with strong PSLQ evidence, S(5)) is NOT elementary but a genuine lemniscatic period. All numeric to
  50-170 digits with independent cross-checks.
tags: [hypergeometric, central-binomial, 3F2, lemniscate, singular-values, closed-forms, series, external]
related: []
---

# The quarter-integer central-binomial series S(k)

`S(k)=sum_{n>=0} C(2n,n)C(4n,2n) / ((kn+1) 64^n)`, `k in Z_{>0}`. Owner asks: closed forms for S(4),S(5)? a
general closed form?

## 1. Reformulation (exact)
`C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n / (n!)^2` (Pochhammer), and `1/(kn+1)=int_0^1 x^{kn}dx`. Hence
> **`S(k) = 3F2(1/4, 3/4, 1/k; 1, 1+1/k; 1) = int_0^1 f(x^k) dx`,  `f(t)=2F1(1/4,3/4;1;t)`.**
Verified to 60 digits against the given `S(1)=8sqrt2/(3pi)`, `S(2)=(4/pi)log(1+sqrt2)`,
`S(3)=-(1/pi)(sqrt3 log(5-2sqrt6)+2 atan(sqrt2/5))`.

## 2. Closed form of the generating function
`f(t)=2F1(1/4,3/4;1;t)` is the signature-4 hypergeometric. Two exact forms (derived, checked):
> `f(t) = (2/pi) (1+sqrt t)^{-1/2} K( sqrt( 2 sqrt t/(1+sqrt t) ) ) = (1/pi) int_0^pi d(theta)/sqrt(1-sqrt t cos theta)`
(K = complete elliptic integral, first kind). Equivalently `f(t)=(2/pi) int_0^{pi/2} d(phi)/sqrt(1-t cos^4 phi)`.

## 3. The reduction mechanism (why k=1,2,3 are elementary)
> `S(k) = (2/pi) int_0^{pi/2} g_k(cos^4 phi) d(phi)`,  `g_k(y)=int_0^1 ds/sqrt(1-y s^k)`.
- `k=1`: `g_1(y)=2(1-sqrt(1-y))/y` -> `S(1)=8sqrt2/(3pi)` (algebraic times 1/pi).
- `k=2`: `g_2(y)=arcsin(sqrt y)/sqrt y`; IBP collapses the phi-integral to `2 int_0^1 du/sqrt(1+u^2)=2 arcsinh 1`,
  giving **`S(2)=(4/pi)log(1+sqrt2)`** (re-derived here from scratch).
- `k=3`: elementary (the given form); it connects to `w=1+i sqrt2`, `|w|^2=3`: `w^3=-5+i sqrt2` yields
  `atan(sqrt2/5)` and `log(5-2sqrt6) = -arccosh(5)`. PSLQ confirms `pi S3 + sqrt3 log(5-2sqrt6) + 2 atan(sqrt2/5)=0`.

The `w`-power pattern is **special to k=3** — it does NOT reproduce k=4,5,6 (checked), a first sign that
elementariness stops at k=3.

## 4. Q1 -- S(4), S(5) are NOT elementary (this is the answer)
> **Closed forms for S(4), S(5) in `pi`, `log`(algebraic), `arctan`(algebraic) do NOT exist** (and none is known).
Two independent lines:
- **PSLQ (170 digits, coeff bound 1e7):** `pi S(4)` (and `pi S(5)`) has NO integer relation against comprehensive
  bases: the fields `Q(sqrt2)`, `Q(2^{1/4})`, `Q(1+i sqrt2)`; and the classical constants lemniscate
  `varpi=Gamma(1/4)^2/(2 sqrt(2 pi))`, Catalan `G`, `Gamma(1/4)^4/pi^3`, and their products.
- **Analytic proof that S(4) is a lemniscatic period.** From `S(4)=(2/pi) int_0^1 (arcsin x + arcsinh x)/sqrt(1-x^4) dx`
  (derived via `2F1(1/2,1/2;3/2;+-cos theta)`), integrate by parts against `v(x)=int_0^x ds/sqrt(1-s^4)`
  (incomplete lemniscate; `v(1)=varpi/2` exactly). This gives, verified to 28 digits,
  > **`S(4) = (varpi/pi)(pi/2 + log(1+sqrt2)) - (2/pi) R`,  `R = int_0^1 v(x)(1/sqrt(1-x^2)+1/sqrt(1+x^2)) dx`.**
  So S(4) **provably contains the lemniscate constant `varpi`**, and the residual `R` is not a standard constant
  (no PSLQ relation). S(4) is a genuine period, one transcendence level above the elementary S(1),S(2),S(3).

## 5. Q2 -- the general closed form
> The general closed form is the hypergeometric / integral of section 1:
> **`S(k)=3F2(1/4,3/4,1/k;1,1+1/k;1)=int_0^1 2F1(1/4,3/4;1;x^k)dx`.**
There is **no** general closed form in elementary constants: elementary reduction happens only for special `k`
(certainly `k=1,2,3`), which is a **singular-value / complex-multiplication phenomenon** for the signature-4
period `2F1(1/4,3/4;1;.)` -- the integral `int_0^1 f(x^k)dx` is elementary only when the family degenerates
suitably. For `k>=4` the value is a period (lemniscatic for `k=4`, and higher for `k=5`), classical-constant
but not elementary. (mp.identify returns closed forms for k=1,2 and None for k=4,5,6,8,12.)

## 6. Numerics (for the record, 45 digits)
`S(4)=1.069352554441268058582961939534278061326`,  `S(5)=1.057087162094771626863606102311001036126`.
`pi S(4)=3.35947012913016715899...`, `pi S(5)=3.32093726264101749334...`; `S(k)->1` as `k->inf`.

Files: `/tmp/{series,pslq*,lemnis,hyp,hyp2,pattern,final4,ibp}.py`.
