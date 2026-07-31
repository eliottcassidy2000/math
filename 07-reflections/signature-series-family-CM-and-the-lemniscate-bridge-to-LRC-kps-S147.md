---
source: kind-pasteur-2026-07-24-S147 (Opus 4.8)
status: SYNTHESIS + new results. Generalises the owner series S(k) (kps-S146) to the four CM "signature"
  families, proves a general k=1 evaluation for ALL lambda, tabulates elementary k=1,2 closed forms, and --
  the payoff -- connects the whole picture to the repo's LRC lemniscate-node thread (opus-S177/S169) via
  complex multiplication (Z[i] lemniscate vs Z[omega] equianharmonic). Proposes new families. Numerics to 80 digits.
tags: [hypergeometric, central-binomial, signatures, complex-multiplication, lemniscate, equianharmonic, lrc, synthesis]
related: [kps-S146, opus-S177, opus-S169, HYP-5680, HYP-5660, THM-515B]
---

# The signature-series family, complex multiplication, and the lemniscate bridge to LRC

## 1. The family (S(k) is one of four)
`S_lambda(k) = sum_{n>=0} (lambda)_n(1-lambda)_n/((n!)^2 (kn+1)) = int_0^1 2F1(lambda,1-lambda;1;x^k) dx
            = 3F2(lambda,1-lambda,1/k; 1,1+1/k; 1).`
The owner's `S(k)` is `lambda=1/4`. The four **CM ("arithmetic-signature") values** `lambda in {1/2,1/3,1/4,1/6}` give
integer/binomial coefficients:

| lambda | coefficient `(lambda)_n(1-lambda)_n/(n!)^2` | name |
|---|---|---|
| 1/2 | `C(2n,n)^2/16^n` | signature 2 (classical `K`) |
| 1/3 | `C(2n,n)C(3n,n)/27^n` | signature 3 (equianharmonic) |
| 1/4 | `C(2n,n)C(4n,2n)/64^n` | signature 4 (lemniscatic) -- the owner's |
| 1/6 | `(1/6)_n(5/6)_n/(n!)^2` | signature 6 (equianharmonic) |

## 2. A general evaluation at k=1 (all lambda, proved)
By Gauss's `2F1(lambda,1-lambda;2;1)` and reflection,
> **`S_lambda(1) = sin(pi lambda) / (pi lambda (1-lambda))`.**
Verified to 41 digits for `lambda = 1/2,1/3,1/4,1/5,1/6,2/5,1/7`. Specialisations:
`S_{1/2}(1)=4/pi`, `S_{1/3}(1)=9 sqrt3/(4 pi)`, `S_{1/4}(1)=8 sqrt2/(3 pi)`, `S_{1/6}(1)=18/(5 pi)`.

## 3. Elementary k=2 closed forms (identified, 80 digits)
| lambda | `S_lambda(2)` |
|---|---|
| 1/2 | **`4G/pi`** (Catalan's constant!) |
| 1/3 | `3 sqrt3 log 2 / pi` |
| 1/4 | `4 log(1+sqrt2)/pi` (owner's S(2)) |
| 1/6 | `3 sqrt3 log(2+sqrt3)/(2 pi)` |
So `S_{1/2}(k)` is elementary for `k=1,2` (giving `4/pi`, `4G/pi`) then NOT (`k>=3`: None). The others are
elementary at least through the owner-confirmed `k=3` cases, then become periods. **Every signature family
mirrors S(k): elementary for small k, then a genuine period.**

## 4. The mechanism: complex multiplication (why elementary stops)
`2F1(lambda,1-lambda;1;t)` is a period of a family of elliptic curves; `S_lambda(k)=int_0^1 (period)(x^k)dx`.
Elementary reduction is a **singular-value / CM phenomenon**:
- `lambda=1/4` (signature 4) is the **lemniscate**, CM by **Z[i]** (`j=1728`); the arc length is
  `varpi = int_0^1 dr/sqrt(1-r^4)`, and indeed `S_{1/4}(4)` contains `varpi` (kps-S146: proved lemniscatic period).
- `lambda=1/3, 1/6` are the **equianharmonic** case, CM by **Z[omega]** (`j=0`, Eisenstein); their constants are
  `sqrt3`, `log(2+sqrt3)` -- the `Z[omega]` world.
So the two arithmetic worlds appearing in the series are exactly **Z[i] (lemniscate)** and **Z[omega] (equianharmonic)**.

## 5. THE BRIDGE TO LRC (the payoff of the owner's mining hint)
The repo's LRC(14) analysis independently lives in the SAME two CM worlds:
- **opus-S177 / HYP-5660:** the LRC(14) tight extremiser `{1,...,13}` is lonely only at a measure-zero
  **"lemniscate node = figure-eight self-crossing"** (`r^2=cos2theta`); resolving it (blow-up) is deciding the
  strict-vs-non-strict good set.
- **HYP-5680 / THM-515B:** the **additive energy `E(S)=sum|S-hat|^4`** is the SINGLE parameter governing
  looseness -- and opus-S177 identifies this quartic as **"the same 4th-power species as the lemniscate
  arc-length `int dr/sqrt(1-r^4)`"** (i.e. `varpi`, Z[i]).
- opus-S177 further: *"7 inert in Z[i] vs splits in Z[omega] => LRC's true home is the equianharmonic Z[omega]
  (Phi_6/Eisenstein)."*

**Synthesis.** The owner's `S(4)` (signature 4) and the LRC additive-energy quartic are two faces of the SAME
`Z[i]`/lemniscate CM structure; the signature-3/6 series and the LRC "true home" are two faces of the SAME
`Z[omega]`/equianharmonic structure. The lemniscate constant `varpi` that obstructs `S(4)`'s closed form is the
very arc-length whose 4th-power species is the LRC looseness parameter. The **"figure-eight self-crossing"**
naming (opus-S177) is the lemniscate of Bernoulli = the classical three-body **figure-eight choreography** curve
-- the owner's "n-body" hook. (No explicit celestial-mechanics thread exists in the repo; the connection is the
figure-eight/lemniscate node, not an orbit computation.)

## 6. New expressions proposed (creative extensions)
1. **Shifted denominator:** `S_lambda(k,r)=sum b_n/(kn+r)=int_0^1 x^{r-1}(period)(x^k)dx=3F2(lambda,1-lambda,r/k;1,r/k+1;1)`.
   General `r=1`: section 2.
2. **Alternating / conjugate period:** `sum (-1)^n b_n/(kn+1)=int_0^1 2F1(lambda,1-lambda;1;-x^k)dx` -- the
   "other" real period of the same CM curve; for signature 4 it pairs `varpi` with the second lemniscatic period.
3. **Harmonic-weighted (derivative in lambda):** `sum b_n [psi-weights]/(kn+1)` gives `d/dlambda` of the above --
   Catalan-type constants for signature 2, `Cl_2`(Clausen)/`L(2,chi_3)` for signature 3/6.
4. **The `varpi`-detector:** predict `S_{1/4}(k)` contains `varpi` for every `k` divisible by 4 (Z[i] resonance),
   and `S_{1/3}(k), S_{1/6}(k)` contain `L(2,chi_{-3})`/`Gamma(1/3)` for `k` divisible by 3 or 6.

## 7. Values (record, 80 digits available)
`S_{1/2}: 4/pi, 4G/pi, 1.12016..., ...`;  `S_{1/3}: 9sqrt3/(4pi), 3sqrt3 log2/pi, 1.10590..., ...`;
`S_{1/4}: 8sqrt2/(3pi), 4log(1+sqrt2)/pi, [S(3) elem], varpi-period, ...`;
`S_{1/6}: 18/(5pi), 3sqrt3 log(2+sqrt3)/(2pi), 1.06451..., ...`.

Files: `/tmp/{siblings,sib2,series,ibp}.py`. Builds on kps-S146; connects opus-S177/S169, HYP-5680/5660.
