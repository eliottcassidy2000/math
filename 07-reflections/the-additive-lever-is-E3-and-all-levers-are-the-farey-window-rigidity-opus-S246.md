---
source: opus-2026-07-11-S246
status: A reframing synthesis. R2 fails (translation-invariant); E3 (Schur, LEM-015-proved) is the right
  additive lever, cleanly separating divisor-complete (E3 ≤ 30, M ≥ 2/27) from the tight interval (E3 = 36,
  M = 1/14). KEY REFRAME: divisor-complete families are LOOSE (M ≥ 2/27 ≫ 1/14) — the clean-ruler "residual"
  was a CERTIFICATE issue, not loneliness. LRC(14) reduces to the 13-runner FAREY-WINDOW RIGIDITY
  [M < 2/27 ⟹ dilated interval {1..13}] (HYP-4151 at k=13, proved r=1), where ALL levers meet and the hard
  case is near-AP, not divisor-complete.
tags:
  - lrc14
  - E3-schur
  - farey-window
  - rigidity
  - reframe
  - unification
  - HYP-4151
---

# The additive lever is E₃, and all levers are the Farey-window rigidity

**opus-2026-07-11-S246.** Owner: keep working the remaining open math. Pursuing the additive-energy lever
(S245's Lever 3) resolves *which* additive invariant works and, via the Farey ladder, reframes the whole arc.

## R2 fails; E₃ is the right lever

The reduced additive energy `R2` is **translation-invariant**: `R2({1..13}) = R2({2..14}) = 1300`, yet
`M({1..13}) = 1/14` (tight) and `M({2..14}) = 1/8` (loose). So `R2` cannot discriminate (opus-S181's
"necessary-not-sufficient"). The **Schur-triple count `E₃ = #{i<j : vᵢ+vⱼ ∈ S}`** is **translation-sensitive**
and is the right invariant: **LEM-015 (PROVED)** — the interval uniquely maximizes `E₃`, with `E₃({1..13}) =
36`; and `E₃({2..14}) = 30` (the shift drops it). Divisor-complete families have `E₃ ≤ 30 ≪ 36` (mean 3.7),
`corr(E₃, M) = −0.34…−0.45` (vs R2's +0.15). A clean E₃-deficit separation.

## The reframe: the clean-ruler residual is LOOSE

Verified: **every divisor-complete family has `M ≥ 2/27 = 0.074`** (min `M = 28/191 = 0.147`, *far* above
`1/14 = 0.0714`). So the divisor-complete "residual" of the clean-ruler route is **loose** — LRC(14) holds
for these families *trivially by margin*. Their difficulty (my S230–S245 arc) was the **bounded-modulus
clean-ruler certificate `hB5`** (a Lean-formalization obligation), **not loneliness**. The genuine
loneliness-hard families are the **near-AP** ones (`M` near `1/14`), which are *not* divisor-complete.

## The unification: all levers are the Farey-window rigidity

The achievable `M`-spectrum for 13 speeds has an **empty window** `(1/14, 2/27)` (HYP-4306, verified;
`2/27 = mediant(1/14, 1/13)`). So the single statement

> **`M < 2/27` ⟹ the family is a dilated interval `{1..13}` (with `M = 1/14`).**

**proves LRC(14)**: a counterexample (`M < 1/14 < 2/27`) would be forced to be the dilated interval
(`M = 1/14`), a contradiction — so no `M < 1/14`. This is **HYP-4151's k=12 rigidity** (`M < 2/25 ⟹ dilated
AP`, PROVED for `r=1` via equioscillation/three-distance, "residues form the AP") **at k = 13**. And **every
lever the fleet has developed reduces to it**:

- three-gap / AP-is-best-coverer (S239) — the coverage face of the rigidity;
- pigeonhole / fold-class clustering (S242/S245) — the discrete face;
- E₃ / additive (this session) — the additive face (`E₃ = 36` ⟺ interval ⟺ tight);
- Farey-ladder (HYP-4306) — the Diophantine face.

One rigidity, four faces. The hard case is **near-AP** (`M` near `1/14`), and `14 = 2·7` composite (vs `13`
prime for the proved k=12 case) is the extra difficulty (the apex zero-divisor).

## Net

- **The additive lever is E₃, not R2** (translation-sensitive; LEM-015-proved max at the interval).
  Divisor-complete cleanly separates (`E₃ ≤ 30`, `M ≥ 2/27`).
- **Reframe:** the clean-ruler residual (divisor-complete) is **loose** — its hardness was the certificate,
  not loneliness. This is honest context for the whole S230–S245 arc: it built the certificate side; the
  loneliness side is the Farey rigidity.
- **The sharpest closing target** is the **13-runner Farey-window rigidity** `M < 2/27 ⟹ {1..13}` — a single
  Diophantine equioscillation statement, **proved for k=12/r=1** (HYP-4151), where all levers meet. Extending
  HYP-4151's residues-form-the-AP equioscillation argument from k=12 (mod-13, a field) to k=13 (mod-14,
  composite) is the concrete open step.

→ LEM-015 (E₃ extremal, proved), opus-S181/S182 (energy necessary-not-sufficient; E₃ not R2), HYP-4306 (Farey
ladder / empty window), HYP-4151 (the k=12 rigidity, proved r=1 — the closing target at k=13), opus-S239
(three-gap face), opus-S242/S245 (pigeonhole face). Files: `lrc14_E3_farey_window_synthesis_opus_S246.py`
(+`.out`).
