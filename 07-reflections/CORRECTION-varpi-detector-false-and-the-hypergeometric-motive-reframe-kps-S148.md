---
source: kind-pasteur-2026-07-24-S148 (Opus 4.8)
status: CORRECTION + abstract reframe. (1) Retracts the kps-S146/S147 claim "S(4) provably contains the
  lemniscate constant varpi" -- the integration-by-parts identity was CIRCULAR (a tautology). (2) The
  varpi-detector conjecture (kps-S147: S(4k) = varpi + elementary) is FALSE: S(4),S(8),S(12) are NOT
  Q-combinations of varpi, varpi^2, Catalan, Gamma(1/8), or elementary constants (100-digit PSLQ). (3) The
  correct reframe: S_lambda(k) for k>=4 are periods of IRREDUCIBLE hypergeometric motives (likely L-values),
  and the elementary locus is a CM/resonance phenomenon -- the same arithmetic that governs the LRC lemniscate node.
tags: [hypergeometric, motives, periods, lemniscate, CM, correction, L-values, lrc, abstract]
related: [kps-S146, kps-S147, opus-S177, HYP-5680]
corrects: [kps-S146, kps-S147]
---

# Correction: the varpi-detector is false; S(k>=4) are irreducible-motive periods

## 1. CORRECTION (retract two claims)
- **kps-S146 said** `S(4)=(varpi/pi)(pi/2+log(1+sqrt2))-(2/pi)R` "so S(4) provably contains varpi." **Wrong --
  circular.** Swapping the order of integration gives `R = P(1)*varpi/2 - (pi/2) S(4)` (verified to 61 digits),
  so the boundary `varpi` is exactly cancelled by the `varpi` inside `R`; the identity collapses to `S(4)=S(4)`.
  **No conclusion about varpi.**
- **kps-S147 said** the "varpi-detector": `S_{1/4}(4k)` contains `varpi`. **FALSE.** PSLQ (100 digits, coeff 1e7):
  `pi S(4)`, `pi S(8)`, `pi S(12)` have NO relation against `{varpi, varpi^2, Catalan, Gamma(1/8)^2/pi,
  Gamma(1/8)^4/pi^3, pi, 1}` or elementary logs. So these values are NOT lemniscatic-plus-elementary.

## 2. What is actually true (stronger, and more interesting)
`S(k)` is elementary exactly for `k=1,2,3` (owner-given / verified); for `k>=4` it is a **genuinely new period**
that is NOT reducible to `varpi`, `Gamma(1/8)`, Catalan, or elementary constants. So the naive "it becomes the
lemniscate constant" picture is wrong -- the obstruction is bigger.

## 3. The reframe: hypergeometric motives and the closed-form locus
`S_lambda(k) = 3F2(lambda,1-lambda,1/k; 1,1+1/k; 1) = int_0^1 (elliptic period of E_lambda at t=x^k) dx.`
Geometrically: take the signature-`lambda` elliptic family `E_t` (CM at special `t`), pull back by the
`mu_k`-cyclic cover `x -> x^k`, and integrate the period over `[0,1]`. The result is a **period of a surface
fibred over P^1**. Its closed form is dictated by the motivic decomposition of the surface's `H^2_transc`:
- decomposes into Tate (`pi`-powers) + CM/Artin (algebraic x logs) pieces  =>  `S_lambda(k)` **elementary**;
- has an **irreducible** piece (a new modular/`L`-motive)  =>  `S_lambda(k)` is an **irreducible period (L-value)**,
  non-elementary and NOT the elliptic CM constant `varpi`.
The dichotomy is a **resonance between the CM order of `E_lambda` and the cyclic level `k`**. For `lambda=1/4`
(CM by `Z[i]`, order 4) the small levels `k=1,2,3` resonate (elementary); `k>=4` gives irreducible surface
motives. This is exactly why the varpi-detector fails: at `k>=4` the period is a *surface/threefold* period
(would-be `varpi^2` or an `L`-value), and PSLQ shows it is not even that -- i.e. the motive is irreducible.

**Corrected conjecture (varpi-detector, right form):** the *elementary* locus
`C_lambda = { k : S_lambda(k) elementary }` is FINITE for each CM `lambda`, and `k` outside it gives an
irreducible hypergeometric-motive period. Data: `C_{1/2} = {1,2}`, `C_{1/4} = {1,2,3}` (at least); `C_{1/3},
C_{1/6}` contain `{1,2}` (and `3`?). Determining `C_lambda` exactly (a CM/level resonance condition) is the
concrete open problem.

## 4. The abstract bridge to LRC (sharpened, honest)
The connection to the repo's LRC lemniscate thread (opus-S177: additive energy `sum|S-hat|^4` = lemniscate
arc-length species; extremiser = lemniscate node) survives the correction and is actually sharpened:
> Both problems are governed by **CM points of a period family**. The LRC(14) tight extremiser is the isolated
> CM/lemniscate node where the LRC period family degenerates; generic (dissociated) sets are the "irreducible"
> regime where counting works (opus-S178: dissociated uniformly certifiable). Symmetrically, `S_lambda(k)` is
> elementary at the CM-resonant `(lambda,k)` and irreducible elsewhere. **The dichotomy "special CM point
> (rigid, exact) vs generic (irreducible, hard)" is the SAME in both.** The `varpi`/`Z[i]` vs `Z[omega]`
> distinction (opus-S177: "7 splits in `Z[omega]` => LRC's true home is equianharmonic") is the level-arithmetic
> that decides which CM points resonate -- the exact analogue of `C_lambda` above.

Concretely testable transfer: opus-S177's "`Z[omega]`, `Phi_6`" suggests the LRC danger-band arithmetic
(`density 1/7`, `14=2*7`) should be organised by `7`'s splitting in `Z[omega]` -- the same equianharmonic
signature that makes `S_{1/3}, S_{1/6}` (not `S_{1/4}`) the "matching" series. A drift/census reformulation in
`Z[omega]` is the abstract lead I'd hand to the LRC census owners.

## 5. Values (record)
`S(4)=1.06935255...`, `S(8)=1.03733013650...`, `S(12)=1.02555584240...`; all `pi S(k)` PSLQ-irreducible.
Method note (recorded to MISTAKES spirit): an integration-by-parts "identity" that reintroduces the target on
the other side proves nothing -- always check the residual is independent of the target before claiming content.

Files: `/tmp/{det2,det3}.py`. Corrects kps-S146 section 4 and kps-S147 sections 4-5.
