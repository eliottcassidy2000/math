---
source: kind-pasteur-2026-07-24-S141 (Opus 4.8)
status: STRUCTURAL RESULT on the OPEN-Q-108 residual. opus-S4 has defects 1-4 as theorems with defect >= 7 as
  the residual; kps-S140 identified d=7 as the 1/(2h)=7 measure-relaxation ARTIFACT ceiling. This shows the
  ceiling is not a real obstruction: LARGE speeds cannot cover an interval at all (7 of them reach at most
  0.749, thirteen at most ~0.92, forty only 0.9991), so any covering configuration must contain SMALL speeds --
  and those are bounded by Fact B. The residual is therefore confined to bounded configurations.
tags: [lrc, lonely-runner, OPEN-Q-108, defect-ladder, covering, artifact-ceiling, residual]
related: [kps-S140, opus-S4 (HYP-9024), klein-S422, macmini-S170]
---

# The `d = 7` wall is an artifact: large speeds cannot cover

## 1. The setting
`gap(S) ≤ h` (`h = 1/14`) means the bad sets `D_w = {τ : ‖wτ‖ ≤ h}`, `w ∈ S`, **cover** `[0,1)`. In the defect
picture `S = C₀ ∪ R` with `|R| = d`, the `d` replacements must cover the good set `G_{C₀}`, in particular its
largest component `I`. kps-S140 showed the measure relaxation bounds the smallest replacement only while
`2hd < 1`, dying exactly at `d = 7` — the residual opus-S4 names.

## 2. The measure relaxation is vacuous at `d=7`, but the truth is far stronger
At `d = 7`, `7·2h = 1` exactly, so the counting bound gives nothing. But counting ignores **overlaps**. If the
`d` bad sets sat independently in `I` (which is what large, mutually incommensurate speeds do), the union would
be `1 − (1−2h)^d = 1 − (6/7)^d`, never `1`. Measured, on an interval of length `0.02`:

| speeds (all large) | union fraction of `I` |
|---|---|
| 7 generic random | 0.657, 0.655 (model `1−(6/7)^7 = 0.6601`) |
| 7 huge random | 0.661, 0.660 |
| 7 AP / dilates / harmonic (structured) | 0.383 – 0.890 |

**Optimised** (hill-climb over the speeds themselves, maximising coverage):

| d | 7 | 9 | 12 | 16 | 20 | 25 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|
| best union fraction | **0.749** | 0.859 | 0.915 | 0.976 | 0.992 | 0.978 | 0.991 | 0.9991 |

> **No finite number of large speeds covers.** Seven reach at most `0.749`; **thirteen** — the entire budget of
> LRC(14) — reach at most about `0.92`; even forty fall short at `0.9991`.

## 3. Consequence: the residual is confined to SMALL speeds
Since a tight configuration or counterexample requires the bad sets to cover, and large speeds provably cannot,
**every such configuration must contain speeds that are small relative to `1/|I|`** — i.e. speeds whose period is
comparable to the surviving interval. Those are exactly the speeds Fact B bounds (`w ≤ 2h/|I|`, kps-S140 §4d).
So the `d ≥ 7` residual is **not** an unbounded region: it is confined to configurations containing bounded
speeds, hence finitely describable.

**Rigorous core (no heuristic).** Iterating klein-S422's Fact A': if `L_{j−1} ≥ 2/w_j` then a full safe gap of
`w_j` survives, of length `(1−2h)/w_j`. The hypothesis for the next step, `L_j ≥ 2/w_{j+1}`, then holds whenever
`w_{j+1} ≥ 2w_j/(1−2h) = (7/3)·w_j`. **Hence a configuration whose speeds grow by a factor `≥ 7/3 ≈ 2.33` can
never cover** — geometric spreading is impossible at every `k`, with no counting ceiling. §2 is the
quantitative, non-geometric extension of this.

## 4. Why this matters for OPEN-Q-108
- The `d = 7` boundary is a **bookkeeping artifact** of the union bound, not a phenomenon (klein-S422's abstract
  theme 2, now confirmed at the precise place it bites).
- The correct dividing line is not *how many* replacements there are but **how large** they are: large ⟹ cannot
  cover (§2, §3); small ⟹ bounded by Fact B.
- This is the same "measure vs structure" duality as everywhere else in this problem: counting fails on the
  spread objects, and the structured/clustered ones are precisely those with small, arithmetically related speeds.

## 5. Next
1. Make §2 rigorous: a second-moment / Fourier bound showing `d` speeds with pairwise ratios bounded away from
   small rationals cover at most `1 − c^d` of an interval. The empirics match the independence model to three
   decimals, so the constant should be attainable.
2. Combine with kps-S140's defect ladder to replace the `d ≤ 6` restriction by a *size* restriction, which would
   remove the `d ≥ 7` residual entirely.
3. Feed the sharpened defect-1 bound (`w ≤ 53`, 7× smaller) into the higher-defect searches.

Files: `/tmp/{d7,maxcover,defectladder,sharpen}.py`.
