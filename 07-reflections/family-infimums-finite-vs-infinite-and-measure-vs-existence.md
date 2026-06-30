# Family infimums: finite vs infinite, measure vs existence — why inf R' = 0 but the floor still holds

*klein-2026-06-29-S16. Asked to understand family infimums better — the caveat behind the finalization (inf R' = 0.344 was a SCAN-inf, not the family-inf). Probing the covering family corrected the picture and sharpened why the descent is essential.*

## The correction: inf R' over the full family is 0, not 0.344

I had recorded `inf R' = 114382/332563 = 0.344` (HYP-3593) as "the floor constant." That is the minimum
over my *scan* (`R subset {1..13}`, size-valid `Q`). The covering family is INFINITE, and along it `R'`
goes to `0`. Holding `Q = {1,2}` and adding high speeds to the binding `R = {1..13}\{7}`:
```
R = {1..13}\{7}            m_R=0.0839   R'=0.344
   + {15}                  m_R=0.0756   R'=0.253
   + {15,16}               m_R=0.0686   R'=0.173
   + {15,16,17}            m_R=0.0627   R'=0.107
   + {15,16,17,18}         ...          R'=0   (m_S = 0: the danger combs cover; lonely measure vanishes)
```
So **`inf R' over the unconstrained covering family = 0`.** The scan-inf `0.344` was a slice; the true
family infimum is `0`, approached as `S` grows into a trivial over-covering.

## Why `inf R' = 0` does NOT break LRC (measure vs existence)

`R' = m_S/(m_R m_Q)` is built from the *measure* `m_S = meas(lonely(S))`. `R' -> 0` means the lonely set's
MEASURE vanishes — it shrinks to isolated points. But LRC is an EXISTENCE statement (`lonely(S) != empty`),
not a positive-measure statement. A measure-zero lonely set is still nonempty (a finite set of lonely
times), so LRC can hold while `R' = 0`. **The floor `R' > 0` (positive measure) is SUFFICIENT for existence,
but not necessary; the measure can legitimately vanish.** Moreover, the `R' = 0` configurations are
over-large `S` (17+ speeds) that trivially cover — they are OUTSIDE the relevant LRC(14) family (the kps
HYP-3415 reduction uses specific size-~13–14 `S`, not arbitrary coverings). So `inf R' = 0` over the naive
infinite family is an artifact of including non-relevant over-coverings.

## The real lesson: a family infimum is provable only over a FINITE (or finitized) family

Two infimums, two families:

| quantity | family | size | infimum | status |
|---|---|---|---|---|
| `rho_j` (per-level decorrelation) | `Z_7`-cores `O` | **FINITE** (`2^7`) | `4cos^2(3pi/7) = 0.198` | **PROVED** (THM-590), a true MIN, attained at the doublet |
| `R'` (floor measure ratio) | all coverings `R u 14Q` | **INFINITE** | `0` (cusp limit, measure vanishes) | NOT the LRC object (measure != existence) |

- Over an **infinite** family, a scan only UPPER-bounds the infimum; the true inf is the limit (here `0`),
  and the measure-zero limit is fine for existence but useless as a positive floor.
- Over a **finite** family, the infimum IS the minimum — computable exactly, attained, provable. The `Z_7`-
  core family is finite, so `inf rho_j = 4cos^2(3pi/7)` is a theorem (THM-590), not a scan.

**The 2-adic descent is exactly the FINITIZATION** that converts the infinite covering family into the
finite `Z_7`-core family. That is *why* the descent is the right move and *why* it makes the bound
rigorous: it replaces an infinite-family infimum (which is `0`, and about measure) with a finite-family
minimum (which is `4cos^2(3pi/7) > 0`, and about the discrete/cyclotomic structure). The rigor of the floor
lives entirely in the finitization.

## The measure/existence split is the σ-even/σ-odd split, at the cusp

At the cusp (`m_S -> 0`), the MEASURE (`R'`, the σ-even floor) vanishes — and EXISTENCE (a nonempty,
measure-zero lonely set) is carried by the DISCRETE side: the count of lonely times, the Redei-odd witness
(σ-odd), the finite `Z_7`-core cyclotomic structure. So the floor-measure approach is sufficient in the
bulk but *cannot* be uniform at the cusp — there the measure is genuinely `0`, and existence comes from the
discrete/witness/finitized side. This is the precise reason "the proof lives at the cusps" (mac-mini S29)
and why the measure-floor framing alone is insufficient there: at the cusp you must count, not measure.

## What this corrects and sharpens

- **Correction:** `inf R'` over coverings is `0` (measure vanishes), not `0.344`; `0.344` is the scan-inf
  over `R subset {1..13}`. The HYP-3593/3596 "floor constant `0.344`" should read "scan-inf over the
  size-valid `R subset {1..13}` slice" (a proxy), NOT a family infimum. (The PROVED bound is THM-590's
  finite-family `4cos^2(3pi/7)`, untouched.)
- **Sharpening:** the binding single-drop is the apex `7` (`R'=0.344`); even-speed drops give `R'=1.27`
  (the 2-adic signature); the floor-measure is the σ-even bulk and the cusp needs σ-odd existence.
- **The clean statement:** the LRC floor's rigor is the FINITE-family minimum `4cos^2(3pi/7) > 0` (the
  finitized, discrete, existence-carrying object), not the infinite-family measure infimum (which is `0`).

So "understanding family infimums" resolves to: **only finite families have provable (attained) infimums;
the descent finitizes the covering family; the measure-infimum over the infinite family is `0` and is the
wrong object (measure != existence); the right, provable object is the finite-family minimum, which is the
positive cyclotomic atom.**

See HYP-3597 (this), THM-590 (the finite-family minimum), HYP-3596 (the ledger, corrected),
HYP-3593 (the scan-inf, now re-scoped), HYP-3554 (`CV` unbounded — the same infinite-family pathology),
HYP-3580 (the cusps), [[the-right-frame-audit-when-the-proof-becomes-finite]] (finitization = the right frame).
