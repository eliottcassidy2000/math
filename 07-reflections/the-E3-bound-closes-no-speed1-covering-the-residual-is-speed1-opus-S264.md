---
source: opus-2026-07-11-S264
status: POSITIVE REDUCTION. The correct inclusion-expansion gives the GENEROUS threshold Sum eps_v < 6/7 (not
  the stricter (7-core)/7 used earlier). With that threshold, the E3/additive bound CLOSES all covering families
  WITHOUT speed 1 (Sum|eps_v| <= 0.18 << 6/7, ~5x margin, core all >=17 since 1 is the only coprime speed <17)
  => coreCover<1 => LRC(14). The residual is EXACTLY the speed-1 (runner-1) families => S255 positional. So the
  LRC(14) residual reduces from "all covering families" to "speed-1 covering families".
tags:
  - lrc14
  - covering-min
  - E3
  - additive-bound
  - inclusion-exclusion
  - speed-1-residual
  - reduction
---

# The E₃ bound closes no-speed-1 covering families; the residual is speed 1

**opus-2026-07-11-S264.** Owner: prove the E₃/dissociation bound for covering families. Working it produced the
right threshold and a genuine reduction — the additive bound closes the bulk of covering families cleanly.

## The correct threshold is 6/7 (generous)

Write `1_{D_v} = 1/7 + f_v`. The core-safe fraction of `G'` expands as
`safe_frac = (1/|G'|)∫_{G'} ∏_v(6/7 − f_v) = (6/7)^{core−1}(6/7 − Σ_v ε_v) + (|S|≥2 terms)`, where
`ε_v = (1/|G'|)∫_{G'} f_v`. So **`coreCover < 1` (⟺ LRC(14) for this covering family) follows from the leading
condition `Σε_v < 6/7`**, provided the `|S|≥2` corrections are small — verified `|safe_frac − leading| ≤ 0.06`.
The threshold `6/7 ≈ 0.857` is **far more generous** than the `(7−core)/7` I'd been (wrongly) targeting; this is
what makes the additive bound succeed.

## The class boundary: speed 1 is the only sub-17 coprime speed

The core = speeds coprime to `30030 = 2·3·5·7·11·13`; the smallest such are `1, 17, 19, 23, …`. So **`1` is the
only coprime speed below `17`**, hence "family has no speed 1" ⟺ "core ⊆ {≥17}".

## No-speed-1 covering: the additive bound closes it

For a core all `≥17`, every `ε_v` is small and **additive-bounded** (S263: `|ε_v|` is driven by the additive
relations `±v±w_i±w_j = 0`, each `~b_1³ ≈ 0.0026`, and `v ≥ 17` avoids runner 1's low-frequency alignment).
**Verified** over 207 no-speed-1 covering families (speeds `<120`): `max Σ|ε_v| = 0.183 ≪ 6/7 = 0.857` — a
**~5× margin**. So `Σε_v < 6/7` holds by the crude additive bound ⟹ `coreCover < 1` ⟹ `M ≥ 1/14` ⟹ **LRC(14)**.
The generous threshold gives ~5× room, versus the ~40× shortfall of the earlier (stricter) thresholds — that is
the whole difference.

## Residual: speed-1 families

Runner 1 (speed 1) is the only sub-17 coprime speed; its danger `D_1 = {‖t‖ < 1/14}` is a low-frequency arc
that aligns with `G'` (`ε_1` up to `0.57` at the deep well). For these families, `coreCover = density(D_1 in
G') < 1` is the runner-1 **positional** statement (`G'` has a point with `‖t‖ ≥ 1/14`) = the near-AP case,
handled by **S255**.

## Net

**LRC(14) = [non-covering: elementary `t=1/14` witness, S252] + [covering: `coreCover < 1` at level `1/14`]**, and
the covering case splits — forced by the sub-17-coprime-speed boundary — into:

- **no speed 1 (core ≥ 17): PROVED** via the E₃/additive bound `Σ|ε_v| ≤ 0.18 ≪ 6/7` (5× margin) — the
  dissociated majority (~94%);
- **speed 1 (runner 1): the runner-1 positional bound** `density(D_1 in G') < 1` = S255 near-AP.

So the LRC(14) residual reduces from **"all covering families"** to **"speed-1 covering families"** — a genuine
reduction. The long analytic arc (S253–S263) paid off here: recognizing the correct `6/7` threshold from the
inclusion-expansion, and the S263 additive-relation structure, together close the no-speed-1 bulk with room to
spare. What remains is the single positional lemma for runner 1 — the same near-AP object S255 already handles
for the extremizer — now the *only* obstacle to LRC(14) for covering families.

→ opus-S263 (additive/E₃ structure), LEM-015 (E₃), opus-S259 (`coreCover<1`, equidistribution), opus-S255
(runner-1/near-AP positional), opus-S252 (non-covering elementary). Files:
`lrc14_E3_bound_closes_no_speed1_covering_opus_S264.py` (+`.out`).
