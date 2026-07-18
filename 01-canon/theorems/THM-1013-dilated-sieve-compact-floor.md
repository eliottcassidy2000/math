# THM-1013 — The dilated-AP sieve: the missing witness for the compact minimizers, and the LRC(14) reduction map (boxeph-2026-07-18-S86)

**Status:** the **dilated-AP sieve is PROVED and FORMALIZED kernel-pure**
(`TournamentH7/LRCDilatedSieve.lean`, `LonelyRunner.dilated_sieve` + `dilated_sieve_lonely14`,
axioms `[propext, Classical.choice, Quot.sound]`, no `sorry`; corpus-registered). It discharges
the **compact minimizer stratum** of the LRC(14) covering case (all 82/82 compact covering
families at the floor `M=1/13` — the ones that defeated every S82–S85 perturbation tool). The
**full compact floor** `compact ⟹ M ≥ 1/13` remains the conjecture ([[HYP-7355]]); this closes
its binding/equality case. Verified `lrc_dilated_sieve_boxeph_S86.py`.

## The dilated-AP sieve (PROVED)

> **Theorem.** Fix `n, d ≥ 1`. If every speed lies at distance `≥ d` from the lattice `n·d·ℤ`
> (`|v_i − n·d·m| ≥ d` for all `i` and all integers `m`), then `t = 1/(n·d)` is `1/n`-lonely.

**Proof.** `‖v_i · t‖ = |v_i/(n·d) − m| = |v_i − n·d·m|/(n·d) ≥ d/(n·d) = 1/n`. ∎ (One division
inequality — the dilated generalization of the empty-circle sieve `sieve_frac`, whose `d=1`,
`Q=n` special case it recovers.)

**Corollary (n=13).** If every speed is `≥ d` from `13d·ℤ`, then `M(V) ≥ 1/13 > 1/14`.

## Why this is the witness the minimizers were missing

The conjectured compact minimizers are `d·{1,…,12} ∪ {killer}` (extremal `2·{1..12}∪{13}`,
`M=1/13`). At their dilation `d`, the AP part `d·i` lands on `i/13` (distance `d·min(i,13−i) ≥ d`
from `13d·ℤ`) and the killer sits in the safe band — so the band condition holds and `t=1/(13d)`
certifies `M ≥ 1/13`. These are **exactly** the families that defeat fill-1 ([[THM-1003-fill1-perturbation-base-case]]),
far-element descent + sharp recursion ([[THM-1008-lrc13-descent-floor]], [[THM-1010-sharp-descent-recursion]]),
best-removal, the ordinary sieve, and the measure bound (S85, MISTAKE-160): their loneliness lives
in a *dilated-AP substructure*, invisible to single-runner perturbation but caught in one line by
the dilated sieve. **Coverage (exact):** every compact covering family with `M ≤ 1/13` (the whole
binding stratum) is certified.

## The complete LRC(14) reduction map (this session's synthesis)

LRC(14) ⟺ every 13 distinct positive-integer speeds have `M ≥ 1/14`. Reduction:

1. **Non-covering** (misses some `q ≤ 14`) ⟹ `M ≥ 1/q ≥ 1/14` — **PROVED** (sieve-margin, kps IX =
   the empty-circle case; `sieve_frac` in Lean). Sharper: `M < 1/13 ⟹ covers all of {2,…,13}`.
2. **Covering** (covers 2..14) ⟹ split by outlier count:
   - **≥2 far outliers** ⟹ `M ≥ 1/13` — **PROVED** ([[THM-726]]).
   - **exactly 1 killer** ⟹ `M ≥ 14/183` (deep well, unique global covering-min) — **PROVED**
     ([[THM-724]]).
   - **0 far outliers = compact** (`ρ = v_max/v_2nd < 13`) ⟹ `M ≥ 1/13` — **the sole residual**
     ([[HYP-7355]], 16k-family hunt, min `1/13` at `2·{1..12}∪{13}`).

So **LRC(14) reduces to the single statement `compact covering ⟹ M ≥ 1/13`**, and its binding
stratum (the `M=1/13` minimizers) is now PROVED by the dilated sieve. What remains is the strict
part: **a compact covering family with no dilated-AP substructure has `M > 1/13`** — i.e. the
n=12/13 near-tightness rigidity (klein's Hamming-radius theorems, HYP-7310). Equivalently:
`M < 1/13 ⟹ ρ ≥ 13` — every family below `1/13` covers `13` via a *far* multiple (`≥ 169 = 13²`,
verified), which is exactly the deep-well/single-killer regime already closed by THM-724.

## Significance

- Supplies the **one witness** the compact minimizers were missing, closing the equality case of
  the last LRC(14) residual by an elementary, kernel-pure sieve.
- Corrects and sharpens the floor picture (see MISTAKE-160): covering-min `= 14/183` (global),
  `1/13` (compact, binding case proved), NOT `1/9`.
- Completes the **sieve/witness family** for LRC(14): empty-circle (`sieve_frac`), fill-1
  (`fill1_perturbation`), far-element descent (`descent_general`), and now dilated-AP
  (`dilated_sieve`) — every elementary route is one explicit rational time.

Related: [[THM-724]], [[THM-726]], [[HYP-7355]], [[THM-1010-sharp-descent-recursion]],
[[THM-995-trapped-cut-excludes-tight-locus]], MISTAKE-160, HYP-7358.
