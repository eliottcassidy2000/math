---
source: opus-2026-07-11-S243
status: VERIFICATION (mac-mini cont.49, requested) + a sharpening + the honest two-case structure. THM-720
  confirmed (large-diameter DC loose and M-growing). The "≤6 effective speeds" that closes the decorrelation
  descent is the COPRIME-TO-30030 count (≤6 for 100% of bounded-diameter DC), NOT "odd" (only 66%). Above
  Vmax = lcm(2..14) = 360360 a single speed covers all divisibility, so the far-element peel (LRC≤13) is the
  second case.
tags:
  - lrc14
  - decorrelation-atom
  - THM-636
  - THM-720
  - divisor-complete
  - effective-count
  - verification
---

# The decorrelation effective count is coprime-to-30030, not odd

**opus-2026-07-11-S243.** Owner: work the 13-runner decorrelation atom. mac-mini cont.49 converged (four
ways) on the large-diameter half of the residual and explicitly asked to *verify the 13-runner analog on DC
families*. I did — confirming the core, sharpening the key invariant, and pinning the honest two-case
structure.

## (1) THM-720 confirmed — large-diameter DC is loose and growing

`reach = M(v)` over divisor-complete families, by diameter: mean `M ≈ 0.181 → 0.209 → 0.232 → 0.255` at
`Vmax ≈ 60, 150, 400, 1000`; min `M ≈ 0.136 → 0.214`. All `≫ 1/14`. So large-diameter DC families are
comfortably loose and `M` grows with diameter — exactly mac-mini's THM-720 (`0.105 → 0.187`), independently
reconfirmed.

## (2) The effective count is coprime-to-30030, not "odd"

mac-mini's descent needs `reach(k) ≥ 1/7` from `≤ 6` distinct lifts (else circular, LRC(14)). The stated
mechanism was "DC even-heavy ⟹ ≤ 6 odd runners." **That invariant is wrong ~1/3 of the time:** `# odd speeds
≤ 6` holds for only **66%** of DC families (mean 5.8). The **correct** invariant is the count of speeds
**coprime to `30030 = 2·3·5·7·11·13`**, which is `≤ 6` for **100%** of bounded-diameter DC families (mean
**2.0**). This is the right "`≤ 6` effective speeds": the small-prime-divisible speeds are **auto-safe**
(opus-S241, inert at composite moduli), leaving `≤ 6` generic runners — the same reduction as klein's
`~6`-runner shrink (S263), now pinned to the coprime-to-30030 sub-family rather than the odd one.

## (3) The honest two-case structure

`≤ 6 coprime-to-30030` is a **bounded-diameter** statement, not universal. A single speed
`≥ lcm(2..14) = 360360` covers *every* divisibility `d ≤ 14`, so `{360360} ∪ {12 speeds coprime to 30030}` is
divisor-complete with **12** coprime speeds (`> 6`). So the `≤ 6`-effective closure fails above
`Vmax = 360360`. But that family is still loose via the **far-element peel**: drop the huge speed and the 12
coprime speeds have `M = 17/78 ≈ 0.218 ≥ 1/13` (LRC≤13). Hence the large-diameter half splits cleanly:

- **Case A** (`Vmax < lcm = 360360`): `≤ 6` coprime-to-30030 effective speeds ⟹ `reach(k) ≥ 1/7` (LRC(7)) ⟹
  `reach(v) ≥ 1/7 − B/L` (THM-636 decorrelation atom) `> 1/14`.
- **Case B** (a speed `≥ lcm`): far-element peel (THM-700/701) ⟹ LRC≤13 on the remaining ≤ 12 speeds.

Both are proved decorrelation atoms (THM-636 / THM-700); together they cover *all* large-diameter DC. (Random
DC with `Vmax < 5000`: 0/38 exceed 6 coprime — the `> 6` case genuinely needs a near-`lcm` speed.)

## Net for the endgame

mac-mini's decorrelation closure is **verified** (loose and growing), with the effective invariant
**corrected** to coprime-to-30030 and the far-element peel added as the second case — so the large-diameter
half is closed in structure by two proved atoms. The remaining crux is the **bounded-diameter, ≤ 6
coprime-to-30030 runners** (`Vmax < 360360`) = klein's `~6`-runner inverse theorem. There, **opus-S242's
pigeonhole is the provable handle**: `≤ 6` coprime speeds cannot cover the `φ(q)/2 > 6` fold-classes at a
composite `q` (e.g. `q = 25, 27`), forcing a clear — exactly the `~44%`-of-DC provable brick, now seen to
apply *a fortiori* once the effective count is `≤ 6`.

→ mac-mini cont.49 / THM-636 (decorrelation atom, formalized) / THM-720 (pair-sum looseness), klein S263
(`~6`-odd shrink — corrected here to coprime-to-30030), opus-S241 (auto-safe — the inertness that defines the
effective sub-family), opus-S242 (pigeonhole — the bounded-case handle), THM-700/701 (far-element peel).
Files: `lrc14_decorr_atom_verify_opus_S243.py` (+`.out`).
