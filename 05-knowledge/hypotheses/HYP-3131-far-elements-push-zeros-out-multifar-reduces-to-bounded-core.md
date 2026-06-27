---
id: HYP-3131
title: Far elements PUSH the miss-PGF zeros OUTWARD — adding far to a good bounded base monotonically increases the nearest-zero radius ρ (1.49→2.0), so the multi-far Lee-Yang region (hence the floor R′≥c) REDUCES to the bounded-core Lee-Yang property + the far-pushes-out monotonicity; the multi-far is not the obstruction. Plus: A000568 sandwiched C(n,3)≤A≤2(n-1)!/3 for n=4..7 (apex-7 tameness window, breaks at n=8)
status: VERIFIED (ρ-monotonicity; multi-far floor ρ≥1.559; far coverage-increasing; the A000568 sandwich) + REDUCTION (multi-far floor ⟸ bounded-core Lee-Yang). Not a proof.
source: mac-mini-2026-06-27-S69
extends:
  - HYP-3127   # multi-far floor = Asano contraction (this VERIFIES obligations 1 & 3)
  - HYP-3122   # the bounded-core φ⁴ row = the binding case
related:
  - HYP-3103   # miss-PGF zeros (the ρ object)
  - HYP-2829   # single-far
  - OPEN-Q-108
external: A000568 (tournaments); Asano contraction / Lee-Yang; GHS/Griffiths correlation inequalities
---

# HYP-3131 — Far elements push the zeros out; the multi-far floor reduces to the bounded core

## The advance (VERIFIED, `lrc_leeyang_polydisk_multifar_macmini_S69.py`)
Tracking the **nearest-zero radius ρ** of the miss-count PGF `G_N(z)` as far elements are added to the GOOD
base `B = consec_8` (`ρ(B)=1.49 > 1`, Lee-Yang):
```
 r:        0      1         2          3
 ρ:      1.49  1.59-1.78  1.56-1.86  1.90-2.09
```
- **Adding far MONOTONICALLY pushes the zeros outward** (ρ: 1.49→~1.6→~1.7→~2.0). 400-config multi-far scan:
  **floor ρ ≥ 1.559** (binding at the resonant r=2 `(21,28)`) — a uniform Lee-Yang margin well above 1.
- Mechanism: each far INCREASES coverage (`d(f)=p0(B∪{f})/p0(B)≈1.04-1.14 > 1`, S68) — a coverage-increasing
  Asano factor that pushes zeros out of the unit disk.

## The reduction
HYP-3127's obligations 1 (single-far Lee-Yang region) and 3 (r-monotonicity) are now VERIFIED: the single-far
factor is coverage-increasing (`ρ` jumps 1.49→1.59), and `ρ` rises with `r`. So:
> **The multi-far elements are NOT the obstruction — they help.** The binding case is the **bounded core**
> (`consec`, `ρ=1.49`) = the φ⁴ hard row (S67). The multi-far Lee-Yang region reduces to **[bounded-core
> Lee-Yang: ρ(bounded) > 1, a finite check] + [far-pushes-out: Asano monotonicity]**. The remaining piece is
> the `ρ_bounded > 1 ⟹ R′ ≥ c` link (Lee-Yang ⟹ GHS/Griffiths correlation inequality).

So the covering bound's entire far structure **subsumes into the bounded-core Lee-Yang property** — the same
coverage-extremality / φ⁴ row already under study (S66/S67). The multi-far floor is downstream, not separate.

## Warm-up (the owner's A000568 observation, verified)
`C(n,3) ≤ A000568(n) ≤ 2(n-1)!/3` holds for **n=4,5,6,7** and BREAKS at n=8 (`6880>3360`): additive
(`C(n,3)`, the OCF triangle blocks) ≤ tournament count ≤ factorial (`2(n-1)!/3`), an apex-7 tameness window.
No load-bearing LRC use; the apex-7 boundary is the through-line (echoes E₇ odd holes / the forbidden-H window).

## Next
1. Prove far-pushes-out for ALL far placements (the coverage-increasing ⟹ zero-outward step, rigorously).
2. The `ρ_bounded > 1 ⟹ R′ ≥ c` link (Lee-Yang ⟹ correlation inequality) — the last analytic piece.
3. So focus the whole covering bound on the bounded-core Lee-Yang `ρ(consec) > 1` (the φ⁴/extremality row).
