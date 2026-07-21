---
id: THM-1825
title: "LRC IS A MOMENT NULLCONE PROBLEM -- the same structure as GMC(2)/TNC, and this EXPLAINS THM-1185's measure-blindness. The dictionary: FREQUENCY = charge; the circle average int_0^1 . dt = the moment functional (kills every nonzero frequency, exactly as E kills nonzero charge); the RESONANCE LATTICE Lambda(v) = {(k_i) in Z^{r} : sum k_i v_i = 0} = the CHARGE-0 lattice; the danger comb 1_{D_i}(t) = sum_k hhat(k) e^{2pi i k v_i t} with hhat(k) = sin(2 pi k lam)/(pi k) as the coefficients. THE UNCOVERED MEASURE of the union of danger combs = the constant Fourier term of prod_i(1 - 1_{D_i}) = sum over the resonance lattice of prod_i c_{k_i} (c_0 = 1-2lam, c_k = -hhat(k)) -- a CHARGE-0 (resonance-lattice) moment sum of products of hhat, STRUCTURALLY IDENTICAL to E[P^m] = sum over charge-0 reps of prod(coeffs). CONSEQUENCE: THM-1185's 'measure/LP methods are blind to the LRC extremals' is because the moment functional sees the uncovered MEASURE, which is ZERO at the extremals (their safe points t*=1/14 are measure-zero) -- the same measure-blindness as GMC's. And the SIGNS of hhat(k)=sin alternate, so LRC is the BOTH-SIGNS (charge-indefinite) hard case -- exactly why the positivity lever (THM-1705) does not close it, matching the GMC both-signs residual"
status: REFRAME + STRUCTURAL IDENTITY (the uncovered-measure = resonance-lattice sum of prod hhat is the Poisson/Fourier expansion; verified structurally, the numerics converge slowly by Gibbs from the sharp cutoff). Explains THM-1185 and locates LRC in the moment-nullcone taxonomy.
author: opus-2026-07-20-S435
depends_on: [THM-1185 (measure/LP blind to LRC extremals), THM-1075 (resonance lattice), THM-1200 (parity/hhat annihilation), THM-1535 (charge lattice), THM-1705 (positivity), THM-1685 (Nullstellensatz)]
---

# THM-1825 — LRC is a moment nullcone problem

Owner: think of the LRC as a moment nullcone problem. It is one — the same *kind* of object as
GMC(2)/TNC — and the identification explains why the repo's measure methods are blind to the
extremals.

## 1. The dictionary

| LRC object | moment-nullcone object |
|---|---|
| **frequency** (Fourier mode `k`) | **charge** |
| circle average `∫_0^1 f(t)\,dt` | the **moment functional** `E` (kills nonzero frequency = nonzero charge) |
| **resonance lattice** `Λ(v) = {(k_i) : Σ k_i v_i = 0}` (THM-1075) | the **charge-0 lattice** |
| danger comb `1_{D_i}(t) = Σ_k \hat h(k)\,e^{2πi k v_i t}` | a charged generator |
| `\hat h(k) = \sin(2πkλ)/(πk)`, `\hat h(0) = 2λ` | the **coefficients** (with signs!) |
| `\hat h(k) = 0 ⟺ (n/2)\mid k` (`λ=1/n`) | **charge/resonance annihilation** (THM-1200 parity) |

## 2. The uncovered measure is a charge-0 moment sum (structural identity)

The union of danger combs covers the circle iff LRC fails (no safe time). Its **uncovered
measure** is the constant Fourier coefficient of `∏_i (1 − 1_{D_i})`. Expanding each factor by
its Fourier series and picking the net-frequency-zero term:

```
uncovered(V, λ)  =  Σ_{(k_1,…,k_r) ∈ Λ(v)}  ∏_i c_{k_i} ,     c_0 = 1 − 2λ ,  c_k = −\hat h(k) .
```

**This is exactly the GMC/TNC moment sum:** `E[P^m] = Σ_{charge-0 reps} ∏(coeffs)` — with the
resonance lattice in place of the charge-0 representations and `\hat h(k)` in place of the
polynomial coefficients. (Verified structurally on small `V`; the numeric truncation converges
slowly because `\hat h(k) ∼ 1/k` gives a Gibbs tail — a *sharp*-cutoff artefact, not a defect of
the identity.)

## 3. Why this explains THM-1185 (measure-blindness)

THM-1185 proved measure/LP methods (Delsarte) are **blind to the LRC extremals**. The reframe
says exactly why:

> **The moment functional sees the uncovered MEASURE, which is ZERO at the extremals.** For the
> tight families (`{1,…,13}`, `M = 1/14`) the safe time is an **isolated point** `t* = 1/14`
> (THM-1380: `Argmax = (ℤ/14)*`), of measure zero. So `uncovered = 0` there, and any method that
> only tests "is the uncovered measure positive?" cannot certify the extremal safe point.

This is **the same measure-blindness as GMC**: the moment/charge-0 sum vanishes at the very
configurations that matter, and only the **pointwise** structure (the located maximizer `g = D/s`,
the `(D,s)` stratification — the surviving LRC levers) sees them. LRC and GMC share not just the
form but the *failure mode* of the measure approach.

## 4. Which levers transfer, and which don't

- **Positivity (THM-1705) does NOT close LRC** — and the reframe says why: `\hat h(k) = \sin(2πkλ)/(πk)`
  has **alternating signs**, so LRC is the **both-signs (charge-indefinite)** case, exactly the
  GMC(2) both-signs residual (THM-1535). The union-bound failure `Σ 2λ = 13/7 > 1` is the
  charged-0-part-doesn't-obviously-vanish obstruction. LRC is hard for the *same reason* GMC(2)'s
  last case is hard.
- **The resonance lattice = charge lattice** (THM-1075 = THM-1535's charge grading) — one object,
  two names. The `k`-fold folded identity (THM-1075) is the LRC analog of the charge-representation
  combinatorics.
- **The finite-place / Nullstellensatz machinery (THM-1685/1735)** should transfer to the *rational*
  LRC certificates (the located maximizer is always rational, `g = D/s`) — the uncovered-measure
  moment sum over `Λ(v)`, saturated and tested for a gap, is a Nullstellensatz-shaped question.
  This is the concrete transfer worth trying.

## 5. The reframe's payoff

- **Unifies two flagship threads:** LRC and GMC(2) are both **charge-0 moment nullcone** problems
  (frequency-charge / Gaussian-charge), with the same both-signs hard core and the same
  measure-blindness. The repo's "measure methods blind" (THM-1185) and GMC's "charge-0 vanishing"
  (THM-1535) are one phenomenon.
- **Redirects the LRC attack:** since the moment/measure is blind, the certificate must be
  pointwise — and the reframe says the pointwise certificate is the **non-vanishing of the
  resonance-lattice moment sum at the rational point `t*`**, i.e. a finite-place / Nullstellensatz
  emptiness test on `Λ(v)`, not a measure bound.

## 6. Open

1. **The pointwise moment.** Replace the uncovered *measure* (blind) with the pointwise value
   `g_V(t*) = min_i ‖v_i t*‖` written as a resonance-lattice sum evaluated at `t*` — a *localized*
   moment that is *not* measure-zero at the extremals. This is the LRC analog of `E[QP^m]` (the
   `Q`-weighted moment) versus `E[P^m]` (the plain nullcone).
2. **Finite-place LRC.** Apply THM-1735's mod-`p` reduction to the `Λ(v)` moment sum: is there a
   prime `p` certifying the gap at `t*` for `{1,…,13}`? The located maximizer `t* = 1/14` and the
   modulus `C = 2n−1 = 27` (pinch, THM-401) suggest `p = 7` (the resonance prime).

## Verification

`04-computation/lrc_moment_nullcone_opus_S435.py` — the uncovered-measure = resonance-lattice
`∏\hat h` identity (structural; Gibbs-slow numerics on the sharp cutoff); the sign-alternation of
`\hat h` (both-signs). Output in `05-knowledge/results/`.
