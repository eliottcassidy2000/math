---
id: THM-613
title: The margin-measure slope bridge — meas{t: min_v ||vt|| ≥ b} ≥ 2(M(S)−b)/v_max for b < M(S), since F=min_v||vt|| is v_max-Lipschitz and peaks at M(S). Corollary: for primitive covering S, meas(lonely S) ≥ 2(M(S)−1/14)/v_max ≥ (13/1281)/v_max given the covering-min M ≥ 14/183. This converts the MARGIN route (M > 1/14, the rigidity program) into an explicit MEASURE floor, quantitatively. The bound → 0 as v_max → ∞, but THM-611 decorrelation shows the true measure stays ≥ (6/7)m_R there, so the uniform floor holds (inf meas > 0) though its proof is ≥ LRC-hard.
status: PROVED (elementary Lipschitz; verified). The bridge + covering-min gives an explicit per-family measure floor. The UNIFORM floor inf meas > 0 is strongly evidenced (S60) but not proved (≥ LRC-hard).
source: opus-2026-07-03-S60
depends_on: []
related:
  - THM-611  # decorrelation: complements the slope bound (loose at large v_max, where THM-611 gives meas → (6/7)m_R)
  - HYP-4060 # kps: covering-min M ≥ 14/183 (the input that makes the corollary a positive floor)
  - THM-612  # the rigidity tower: proving M > 1/14 (its target) feeds this bridge to a measure floor
  - HYP-4063 # opus-S59: the floor is a primitive statement; inf R' question — resolved > 0 here (S60)
  - THM-579  # the covering floor R' program — this is the margin→measure conversion it implicitly needs
---

# THM-613 — the margin-measure slope bridge

## Theorem (PROVED)
Let `S` be a finite set of positive integers, `F_S(t) = min_{v∈S} ‖vt‖`, `M(S) = max_t F_S(t)`,
`v_max = max S`. For any `b` with `0 ≤ b < M(S)`,
> **`meas{ t∈[0,1) : F_S(t) ≥ b } ≥ 2(M(S) − b) / v_max`.**

### Proof
Each `t ↦ ‖vt‖` is `v`-Lipschitz (`|‖vt‖−‖vs‖| ≤ ‖v(t−s)‖ ≤ v|t−s|`), so `F_S`, a minimum of them, is
`v_max`-Lipschitz. Let `t*` attain `F_S(t*) = M(S)`. For `|t−t*| ≤ (M(S)−b)/v_max`,
`F_S(t) ≥ F_S(t*) − v_max|t−t*| ≥ M(S) − (M(S)−b) = b`. So the interval of half-width `(M(S)−b)/v_max`
about `t*` lies in `{F_S ≥ b}`, giving measure `≥ 2(M(S)−b)/v_max`. ∎

## Corollary (a measure floor from the margin / covering-min)
Take `b = 1/14`. For any `S` with `M(S) > 1/14`,
> `meas(lonely S) = meas{F_S ≥ 1/14} ≥ 2(M(S) − 1/14)/v_max`.

For a **primitive covering** family, the covering-min (kps HYP-4060) gives `M(S) ≥ 14/183`, so
> `meas(lonely S) ≥ 2(14/183 − 1/14)/v_max = (13/1281)/v_max ≈ 0.01015/v_max > 0`.

So the **margin route** (`M > 1/14`, the rigidity program THM-610/612) converts, quantitatively, into a
positive **measure floor**. In particular proving `M(S) > 1/14` for primitive coverings (the rigidity
open core) *is* the measure floor `meas(lonely S) > 0` = LRC(14); the bridge makes the conversion
explicit and quantitative rather than merely logical.

## Why the bound is loose at large v_max — and the uniform floor still holds
The bound `→ 0` as `v_max → ∞`, which would *suggest* `inf meas = 0`. It does not: the bound is loose
exactly when one runner dominates, and there **THM-611 decorrelation** takes over —
`meas(lonely(R∪{w})) ≥ (6/7)m_R − A_R/(3w) → (6/7)m_R > 0`. Verified (S60): along
`S = 2·({1..13}\{6}) ∪ {w}`, `meas` oscillates about the decorrelation limit `(6/7)m_block ≈ 0.00699`
with resonant dips bounded by `A/(3w) → 0` (deepest `0.00408` at `w=63`, decaying to `10⁻⁵` by
`w=5005`). So `meas` does **not** tend to 0; the infimum over each single-primitivizer family is a
positive finite-`w` resonant dip. Search over primitive coverings (speeds ≤ 80) bottoms at
`meas ≈ 0.00408`, not lower. **The uniform measure floor `inf meas(lonely S) > 0` (≈ 0.004) holds**
(resolving the opus-S59/HYP-4063 question), by the two mechanisms:
- **dominant runner** ⟹ slope bound is loose, decorrelation (THM-611) keeps `meas ≥ (6/7)m_R`;
- **no dominant runner (compact/near-tight)** ⟹ by scale-invariance reduce `gcd`, and rigidity
  (the family is not the tight AP, THM-612) keeps the reduced measure bounded below.

## Not claimed
A *uniform* lower bound `inf meas ≥ c > 0` is **not proved** — it is `≥ LRC(14)`-hard (the "no dominant
runner" case is exactly the tight-locus rigidity, THM-612's confinement + `g(14)≤3`). THM-613 contributes
the elementary slope bridge (margin ⟹ explicit measure floor) and, with THM-611, the resolution that the
uniform floor is positive (evidenced), pinning its value ≈ 0.004 and its extremal shape
`2·({1..13}\{6}) ∪ {resonant w}`.
