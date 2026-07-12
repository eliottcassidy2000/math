---
source: opus-2026-07-11-S247
status: A CORRECTION to opus-S246. The "Farey-window rigidity" [M < 2/27 ⟹ dilated interval {1..13}] I
  proposed as the LRC(14) crux is FALSE at k=13: the family {1,…,11,13,36} has M = 3/41 = 0.0732 ∈
  (1/14, 2/27). HYP-4306's window-empty was verified only for k=12 (mod-13 PRIME); it does NOT extend to
  k=13 (mod-14 COMPOSITE). The correct crux is "AP minimizes M" (no M < 1/14), not the window-empty. S246's
  E3/divisor-complete-loose/all-levers-unify core STANDS; only the specific window-empty form was wrong.
tags:
  - lrc14
  - correction
  - farey-window
  - k13-composite
  - AP-extremality
  - crux
---

# Correction: the k=13 Farey-window is NOT empty — the crux is "AP minimizes M"

**opus-2026-07-11-S247.** Working S246's proposed closing target (extend HYP-4151's Farey-window rigidity to
k=13) refutes the target itself. Honest correction.

## The refutation

S246 claimed: the M-spectrum window `(1/14, 2/27)` is empty, so `M < 2/27 ⟹ dilated interval {1..13}`, which
would prove LRC(14). **This is false at k=13.** The family

`{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 36}` has `M = 3/41 = 0.07317 ∈ (1/14, 2/27)`

(attained at `(r,Q) = (3, 41)`). It is `{1..13} ∖ {12}` plus `36` — a near-AP with a far element, **lonely**
(`M > 1/14`), not a counterexample. So the window is **non-empty**, and the rigidity `M < 2/27 ⟹ interval`
does **not** hold at k=13.

## Why S246 was wrong: k=12 prime vs k=13 composite

HYP-4306 verified the empty window `(1/(k+1), 2/(2k+1))` for **k = 3..7 and k = 12** — I over-extrapolated to
k=13. The obstruction is exactly the recurring **`14 = 2·7` composite vs `13` prime** difficulty:

- **k=12** (window governed by mod-13, a **field**): HYP-4151 proves `M < 2/25 ⟹` residues form the AP mod 13
  ⟹ `M = 1/13`. The window is empty.
- **k=13** (mod-14, **composite**, zero-divisor at 7): the spectrum near `1/14` is **denser** — `3/41` (and
  its kin) fill the interval `(1/14, 2/27)`. No clean second-rung gap.

So the equioscillation/three-distance rigidity that closes k=12 does **not** transfer to k=13; the composite
modulus admits the extra rung `3/41`. (The `3/41` value is itself a repo landmark — the `K_{3,3}` incidence
wall / near-miss `12→36`, HYP-2934 — so this is the known composite-side structure, not a new anomaly.)

## The correct crux (and what of S246 survives)

LRC(14) is `M ≥ 1/14` for all families = **no `M < 1/14`** = **the AP `{1..13}` minimizes `M`** (global min
`1/14`, verified: min over near-AP perturbations is exactly `1/14`, at the AP). This is the **extremality /
inverse theorem** — *not* the window-empty. The window-fillers (`M = 3/41`, etc.) are lonely and irrelevant
to LRC(14); they only show the spectrum is denser than the clean ladder.

**What survives from S246 (still correct):**
- **E₃ is the right additive invariant** (translation-sensitive; LEM-015-proved max at the interval); R2
  fails. ✓
- **Divisor-complete families are loose** (`M ≥ 2/27`, verified) — the clean-ruler "residual" is a
  certificate issue, not loneliness. ✓
- **All levers unify onto the AP extremality** ("the AP is the unique `M`-minimizer / best coverer") — this
  is the genuine crux. ✓ Only its *specific form* was mis-stated: it is "AP minimizes M", **not** the
  (false) "window `(1/14,2/27)` empty".

## Net

Honest correction: the clean Farey-window rigidity is **false at k=13** (`3/41` realized); the k=12→k=13
extension fails on the composite-14 apex, so HYP-4151 does **not** transfer as I claimed. The crux remains the
**AP `M`-extremality** ("no family beats `1/14`"), which is what the three-gap, pigeonhole, E₃, and energy
levers all reduce to — the correct, still-open unification. The lesson: the composite `14` genuinely
separates k=13 from the proved k=12 case; any clean-ladder/window shortcut that works at k=12 must be checked
at k=13 before use.

→ opus-S246 (corrected here), HYP-4306 (window-empty — k≤12 and k=12 only), HYP-4151 (k=12 rigidity — does NOT
transfer to k=13), HYP-2934 (3/41 = the K33 composite wall), LEM-015 (E₃ — still the right invariant),
opus-S239 (AP extremality — the correct crux). Files: (verification inline; families `{1..11,13,36}` → 3/41).
