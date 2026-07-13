---
source: opus-2026-07-11-S259
status: ADVANCE + correction of S258. Attacking the <=6-core anti-concentration gives a WORKING favorable
  mechanism: the coprime <=6-core (coprime to 30030) EQUIDISTRIBUTES in the good set G' (Weyl), so the
  within-G' union bound Sum_core density(D_v in G') < 1 holds for ALL covering families (coreCover <= 0.65 < 1,
  robust) -- giving a safe point => LRC(14). This SUCCEEDS where the single even-fold (S558o) failed (its
  non-coprime odds concentrate), and CORRECTS S258 (which missed the within-G' union bound, only ruling out the
  global one + second moment). Rigor reduces to an Erdos-Turan discrepancy estimate + a runner-1 positional
  bound (the near-AP alignment = S255 territory).
tags:
  - lrc14
  - covering-min
  - anti-concentration
  - six-core
  - equidistribution
  - union-bound
  - correction
  - advance
---

# The ≤6-core equidistributes in the good set: a working union bound (correcting S258)

**opus-2026-07-11-S259.** Owner: attack the ≤6-core anti-concentration against the good set. Doing so produces
a **working favorable mechanism** and corrects my own S258 pessimism.

## The setup (iterated fold)

A covering family `v` splits, via auto-safe (S241/S243), into its **core** = speeds coprime to
`30030 = 2·3·5·7·11·13` (`≤6` of them, S258) and the **non-core** = the rest (each divisible by a prime `≤13`).
The **good set** `G' = {t : ‖w t‖ ≥ 1/14 ∀ non-core w}` has positive measure (iterated LRC(≤13)). Then

> LRC(14) for `v` ⟺ the core does **not** cover `G'`: `coreCover := |G' ∩ ⋃_core D_v| / |G'| < 1`
> (`D_v = {‖v t‖ < 1/14}`, `|D_v| = 1/7`).

## The mechanism (the attack)

Each core `v` is coprime to `30030`, hence **coprime to the small-prime structure** of the non-core / `G'`. So
`D_v` **equidistributes** in `G'` (Weyl): `density(D_v in G') = |D_v ∩ G'|/|G'| ≈ |D_v| = 1/7`. With `≤6` core
arcs, the **within-`G'` union bound** gives

> `coreCover ≤ Σ_core density(D_v in G') ≈ ≤6/7 < 1` ⟹ safe point ⟹ LRC(14).

This **succeeds where the single even-fold (S558o) failed.** The single even-fold's leftover odd runners
include `3, 5, 7, …` — *not* coprime to the even structure — so they **concentrate** in the even-good set
(density → 1 at the wall). The **iterated** fold's core is coprime to **all** primes `≤ 13`, so it
**equidistributes** (density `1/7`). The move from one fold to the iterated fold is exactly what buys the
decorrelation.

## Verified (all covering families, speeds < 150)

| `|core|` | mean coreCover | max coreCover | mean Σ-density | union-bound `< 1` |
|---|---|---|---|---|
| 1 | 0.152 | 0.211 | 0.152 | all |
| 2 | 0.280 | 0.357 | 0.301 | all |
| 3 | 0.390 | 0.471 | 0.450 | all |
| 4 | 0.497 | 0.568 | 0.606 | all |
| 5 | 0.554 | 0.597 | 0.732 | all |
| 6 | — | 0.649 | 0.926 | all |

- **`coreCover < 1` for every covering family** (max `0.65`) — the ≤6-core anti-concentration **holds
  robustly**.
- **`Σ density < 1` for every covering family** (max `0.926`, at `|core|=6`) — the equidistribution union bound
  **works**. (Empirically more core ⟹ lower per-arc density, so the discrepancies do not add adversarially.)
- Per-arc density `≈ 1/7`: large core runners (>50) mean `0.154`, mid (17–50) mean `0.152` — Weyl confirmed.
  Runner 1 (`v=1`, no comb) is `0.15–0.21` for generic families and `→1` **only at the AP**, which is
  **non-covering** (no multiple of 14), hence outside the target.

## This corrects S258

S258 ruled out the **global** union bound (`|G| > o/7`, S558o) and the **second moment** (Paley-Zygmund is
wrong-direction) and concluded "loose is hard." That **missed** the **within-`G'` union bound with the coprime
core's equidistribution** — which is a working favorable mechanism. **The loose stratum is favorable after
all**, via equidistribution of the ≤6-core in `G'`.

## Residual (empirical → rigorous)

1. **The Weyl discrepancy.** `density(D_v in G') = 1/7 + ε(v)`; the rigorous bound is **Erdős–Turán**
   (`ε` bounded by the Fourier coefficients of `1_{G'}` and `1/v`) — clean for the large core runners. The
   observed `Σ density < 1` (max `0.926`) says the discrepancies do not accumulate, but a rigorous
   `Σ density < 1` needs the Erdős–Turán estimate summed over the core.
2. **Runner 1** (`v=1`): no equidistribution (single arc, coarsest scale). Its density is low (`0.15–0.21`) for
   generic covering families and `→1` only near the AP — the **same near-AP alignment that S255 handles** for
   the tight stratum. It needs a positional bound, not equidistribution.

## Net

The ≤6-core anti-concentration **holds** (`coreCover < 1`, robust) via the **coprime core's equidistribution
in `G'`** — a working favorable mechanism (the within-`G'` union bound), which **corrects S258**. Rigor reduces
to (i) an **Erdős–Turán discrepancy** estimate giving `Σ_core density < 1`, and (ii) a **positional bound on
runner 1** (the near-AP alignment, = S255 territory). This turns the loose stratum from S258's "hard, no clean
tool" into a **concrete favorable route**: equidistribution + Erdős–Turán, with the single non-equidistributing
runner (`v=1`) absorbed into the already-proved tight/near-AP case. The crux is now sharply and *favorably*
located: an effective equidistribution (discrepancy) estimate for the coprime ≤6-core against the structured
good set.

→ opus-S258 (corrected), s558o (even-fold: single fold fails, iterated fold succeeds via coprimality),
opus-S255 (deep-well/tight = the runner-1 alignment), opus-S241/S243 (auto-safe, ≤6 core), LRC(≤13). Files:
`lrc14_six_core_equidistribution_union_bound_opus_S259.py` (+`.out`).
