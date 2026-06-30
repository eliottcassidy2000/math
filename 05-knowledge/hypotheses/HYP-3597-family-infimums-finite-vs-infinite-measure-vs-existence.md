---
id: HYP-3597
title: Understanding family infimums -- inf R' over the INFINITE covering family is 0 (the lonely MEASURE vanishes at the cusp; verified 0.344->0.253->0.173->0.107->0 as high speeds are added), NOT 0.344 (that was the scan-inf over R subset {1..13}); but inf R'=0 does NOT break LRC because LRC is EXISTENCE not positive measure; the provable object is the FINITE-family minimum inf rho_j = 4cos^2(3pi/7) (THM-590, the 2^7 Z_7-cores), and the 2-adic descent is the FINITIZATION that makes the infimum a true (attained, positive) minimum
status: VERIFIED (covering-family probe) + the conceptual correction. Corrects the scope of HYP-3593/3596 (0.344 = scan-inf, not family-inf). THM-590 (the finite-family minimum) is unaffected.
source: klein-2026-06-29-S16
depends_on:
  - THM-590    # the finite Z_7-core family: inf = min = 4cos^2(3pi/7) (PROVED)
  - HYP-3593   # the scan-inf 0.344 (now re-scoped)
related:
  - HYP-3596   # the finalization ledger (caveat corrected)
  - HYP-3554   # CV(N_R)^2 unbounded -- the same infinite-family pathology
  - HYP-3580   # the cusps (where the measure vanishes / existence takes over)
  - HYP-3415   # the kps reduction (which constrains the relevant family)
results:
  - 04-computation/family_infimum_structure_klein.py
  - 05-knowledge/results/family_infimum_structure_klein.out
---

# HYP-3597 — family infimums: finite vs infinite, measure vs existence

## The correction (verified)

`inf R'` over the SCAN (`R subset {1..13}`, size-valid `Q`) is `114382/332563 = 0.344` (HYP-3593). But the
covering family is INFINITE, and `R' -> 0` along it: with `Q={1,2}` fixed, adding high speeds to the binding
`R={1..13}\{7}` gives `R' = 0.344, 0.253, 0.173, 0.107, ... -> 0` (at `{1..13}\{7}+{15,16,17,18}`, `m_S=0`).
So **`inf R' over the full covering family = 0`**, not `0.344`. The `0.344` is the scan-inf over a slice.

## Why `inf R'=0` does NOT break LRC

`R'=m_S/(m_R m_Q)` uses the MEASURE `m_S=meas(lonely(S))`. `R'->0` means the lonely set's MEASURE vanishes
(it shrinks to isolated points). LRC is EXISTENCE (`lonely(S) != empty`), not positive measure -- a
measure-zero lonely set is still nonempty. So the floor `R'>0` (positive measure) is SUFFICIENT for
existence but NOT necessary; the measure can legitimately vanish. The `R'=0` cases are over-large `S`
(17+ speeds) that trivially cover -- OUTSIDE the relevant LRC(14) family (kps HYP-3415 uses size-~13-14 `S`).

## The lesson: provable infimum <=> finite (or finitized) family

| quantity | family | size | infimum | status |
|---|---|---|---|---|
| `rho_j` (per-level) | `Z_7`-cores | FINITE (`2^7`) | `4cos^2(3pi/7)=0.198` | PROVED min (THM-590), attained (doublet) |
| `R'` (floor measure) | all coverings | INFINITE | `0` (cusp limit) | not the LRC object (measure != existence) |

Over an infinite family a scan only UPPER-bounds the inf; the true inf is the cusp limit (`0`). Over a
finite family the inf IS the min -- exact, attained, provable. **The 2-adic descent is the FINITIZATION**
that converts the infinite covering family into the finite `Z_7`-core family; that is why it is the right
move and why the bound is rigorous only there. The rigor lives in the finitization.

## Measure/existence = sigma-even/sigma-odd at the cusp

At the cusp (`m_S->0`) the MEASURE (`R'`, sigma-even floor) vanishes; EXISTENCE (a nonempty measure-zero
lonely set) is carried by the DISCRETE side -- the lonely-time count, the Redei-odd witness (sigma-odd),
the finite `Z_7`-core cyclotomic structure. So the measure-floor is uniform in the bulk but NOT at the
cusp; there one must COUNT, not measure. This is the precise content of "the proof lives at the cusps."

## Sharpenings (verified)

- Binding single-drop = the APEX 7 (`R'=0.344`); even-speed drops give `R'=1.27` (the 2-adic signature);
  odd non-7 drops `0.70-1.0`.
- The same infinite-family pathology is `CV(N_R)^2 -> infinity` (HYP-3554) and `R' -> 0` (here): both say
  the per-set measure object is the wrong (infinite-family) coordinate; the finite cyclotomic core is right.

## Net

Only finite families have provable, attained infimums. The covering family is infinite and its measure-
infimum is `0` (the wrong object). The descent finitizes it to the `Z_7`-cores, where the infimum is the
positive cyclotomic atom `4cos^2(3pi/7)` (THM-590). "The floor constant" is that finite-family minimum, NOT
the infinite-family measure infimum.
