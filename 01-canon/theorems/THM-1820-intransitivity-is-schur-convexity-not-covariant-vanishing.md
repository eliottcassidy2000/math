---
id: THM-1820
title: "MAX INTRANSITIVITY IS SCHUR-CONVEXITY, NOT COVARIANT VANISHING -- and it is a DIFFERENT extremal problem from H-maximisation. CORRECTS THM-1810/HYP-8600's conflation. (1) By Kendall-Babington-Smith, c_3(T) = C(n,3) - sum_i C(s_i,2) (s_i the scores), so MAXIMISING the 3-cycle count (intransitivity) = MINIMISING sum_i C(s_i,2), which is SCHUR-CONVEX in the score sequence -> minimised at the REGULAR tournament (all scores (n-1)/2), doubly-regular/Paley when n=3 mod 4. Verified n=3,4,5: max c_3 = 1,2,5 at the regular/balanced scores, min 0 at transitive. This is a Schur-convexity/variance statement, NOT a covariant vanishing. (2) H-MAXIMAL IS A DISTINCT PROBLEM: the repo already knows the doubly-regular Paley is provably BEATEN for large n as an H-maximiser (H is Schur-CONCAVE, a different extremal). So HYP-8600 conflated two extremal problems; they COINCIDE at small n but DIVERGE for large n. (3) THE CORRECTED SL(2) STATEMENT: the score-variance sum C(s_i,2) is the QUADRATIC (catalecticant/apolar) invariant of the score-generating form; max intransitivity = the APOLAR/HARMONIC stratum where it is minimised = roots maximally spread (regular polygon), of which the n=4 EQUIANHARMONIC j=0 is the harmonic 4-point config and the Gauss-sum character tournament is the arithmetic realisation. So 'max intransitivity = SL(2)-special covariant vanishing' is corrected to 'max intransitivity = the apolar/harmonic locus = minimal quadratic invariant (Schur-convex minimum)'"
status: PROVED (the c_3 = C(n,3) - sum C(s_i,2) formula is classical Schur-convex; regular = minimiser, verified n=3,4,5). Corrects THM-1810 Q1's 'covariant vanishing' framing and separates the two extremal problems in HYP-8600.
author: opus-2026-07-20-S434
corrects: THM-1810 (Q1 covariant-vanishing framing), HYP-8600 (H-extremal vs 3-cycle-maximal conflation)
depends_on: [THM-1800/1810 (dictionary + sub-questions), the repo's Schur-concavity-of-H and 'Paley beaten for large n' threads]
---

# THM-1820 — Intransitivity is Schur-convexity, not covariant vanishing

> **PARTIALLY CORRECTED by THM-1865 (opus-S438).** §1 (`c_3` Schur-convex, regular = intransitivity
> maximiser) stands. But §2's "`H` is Schur-CONCAVE" is **ill-posed**: `H` (Ham-path count) is *not*
> a function of the score sequence at all (one `n=6` score carries six `H` values). The "Paley
> beaten for large `n`" claim is nonetheless **correct** — via the circulant census (LEM-004:
> rotation beats Paley at `n≥13`), *not* via Schur-concavity. Open Q1 ("what maximises `H`?") is
> answered by existing canon (LEM-004/THM-128/THM-212). See THM-1865.

A correction and sharpening of THM-1810 Q1 / HYP-8600, prompted by mining the repo's
H-extremal threads.

## 1. The correction: max intransitivity = regular tournament (Schur-convex)

By the **Kendall–Babington–Smith identity**, for a tournament with score sequence
`(s_1,…,s_n)`,

```
c_3(T)  =  \binom{n}{3} − Σ_i \binom{s_i}{2} .
```

So **maximising the 3-cycle count (intransitivity) = minimising `Σ_i \binom{s_i}{2}`**. That
sum is **Schur-convex** in the score sequence, hence minimised at the **most balanced** point —
the **regular tournament** (all `s_i = (n−1)/2`), which exists for odd `n` and is
doubly-regular (Paley) when `n ≡ 3 mod 4`. Verified `n = 3,4,5`: `max c_3 = 1, 2, 5` at the
regular/near-regular scores, `min = 0` at the transitive (score `(0,…,n−1)`).

> **This is a Schur-convexity / score-variance statement — NOT a covariant vanishing.**
> THM-1810 Q1 said "intransitivity = a covariant vanishing"; that framing is corrected here.
> The right object is the *convexity* of the score-imbalance functional.

## 2. Max intransitivity ≠ H-maximal (two different problems)

The repo's census already establishes that the **doubly-regular Paley tournament is provably
beaten for large `n` as an `H`-maximiser** (`H` = the determinant invariant, which is
**Schur-concave** — a *different* extremal). Therefore:

| extremal problem | maximiser | monotonicity |
|---|---|---|
| **3-cycles (intransitivity)** | regular / doubly-regular (Paley) | `Σ\binom{s_i}{2}` **Schur-convex** (minimise) |
| **`H` (determinant)** | *not* Paley for large `n` | `H` **Schur-concave** (a distinct optimum) |

> **HYP-8600 conflated these.** They coincide at small `n` (where regular = Paley = also
> H-good) but **diverge for large `n`**. The clean statement is only about **intransitivity**:
> that is the regular tournament, period. The `H`-extremiser is a separate, subtler object.

## 3. The corrected SL(2) / binary-form statement: apolar = harmonic

The score sequence `s_i` is the vector of **first moments** (out-degrees) of the configuration.
Its imbalance `Σ\binom{s_i}{2}` is (up to affine normalisation) the **quadratic /
catalecticant invariant** of the score-generating form — the "variance." Hence:

> **Max intransitivity = the APOLAR / HARMONIC stratum**, where the quadratic invariant of the
> score data is minimised = the roots are **maximally spread** (the regular-polygon / harmonic
> configuration). At `n = 4` this is the **equianharmonic `j = 0`** (harmonic 4-point) config;
> the **Gauss-sum character (Paley) tournament** is its arithmetic realisation at `n ≡ 3 mod 4`.

So the slogan "max intransitivity = `SL(2)`-special covariant vanishing" (THM-1810) is corrected
to:

> **Max intransitivity = the apolar/harmonic locus = the minimum of the quadratic (catalecticant)
> invariant** — a Schur-convex minimum, achieved at the maximal-symmetry (regular-polygon)
> configuration, realised arithmetically by the character/Gauss-sum tournament.

The equianharmonic `j = 0` at `n = 4` *is* an `SL(2)`-special point, so the earlier intuition was
directionally right (special stratum = max symmetry = max intransitivity); the *mechanism* is
apolarity/Schur-convexity, not a covariant zero.

## 4. What survives of THM-1810, corrected

| sub-question | corrected answer |
|---|---|
| Q1 (covariant ↔ cycle statistic) | intransitivity = **Schur-convex minimum** of the score-variance = the **apolar/harmonic** stratum (`j=0` at `n=4`); *not* a covariant vanishing |
| Q2 (Paley `d` as Gauss sum) | **stands**: `d(Paley) = ((p+1)/4)^{(p−1)/2}` from the `±i√p` spectrum |
| Q3 (Rédei = discriminant mod 2) | **stands**: Rédei parity = the mod-2 shadow (THM-1425) |
| HYP-8600 (`H`-extremal = character?) | **REFUTED as stated**: `H`-extremal ≠ 3-cycle-extremal; Paley is *not* `H`-max for large `n`. Only the *intransitivity* extremal is the regular/character tournament |

## 5. Open

1. **The `H`-extremiser's own SL(2) character.** Since `H` is Schur-concave and its maximiser
   drifts off Paley for large `n`, what configuration *does* maximise `H`? The repo's
   "interval/most-concentrated-spectrum maximises `H`" suggests the `H`-extremiser is the
   *spectrally concentrated* (not the balanced) tournament — a different binary-form stratum
   (coalescing roots?), worth identifying.
2. **Apolarity made exact.** State `Σ\binom{s_i}{2}` precisely as the apolar/catalecticant
   invariant of an explicit form attached to the tournament, so "harmonic = max intransitivity"
   is a theorem, not an analogy.

## Verification

`04-computation/intransitivity_schur_apolar_opus_S434.py` — the `c_3 = \binom{n}{3} − Σ\binom{s_i}{2}`
maximisation for `n=3,4,5` (regular = maximiser); the separation of the two extremal problems.
Output in `05-knowledge/results/`.
