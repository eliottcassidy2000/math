# The official Finset target of LRC(14) now reduces to `INVcov` in the kernel

*boxeph-2026-07-18-S109; corrected codex-2026-07-18 (MISTAKE-166). The bridge now proves
`LRC14_finset_of_INVcov : LRC(≤13)[cited] + INVcov[open] ⟹ LRC14.LRC14`,
kernel-pure. `INVcov` explicitly assumes `Covering(2..14)`; the original unqualified `INV`
was false, although the transport from the working shape to the Finset target was sound.*

## The gap bridged

`LRC14_of_INVcov` (corrected S108) concludes in the working shape:
`∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t : ℝ, Lonely 14 v t`.
The ledger's target is the human-facing shape:
```
def LRC14.LRC14 : Prop :=
  ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
    ∃ t ∈ Icc (0:ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, 1/14 ≤ |(w:ℝ)*t - a|
```
Two elementary transfers close the gap:

1. **Enumerate the Finset.** `W.equivFinOfCardEq hWcard : W ≃ Fin 13` gives `v i := ((e i : ℕ) : ℤ)`, a
   `Fin 13 → ℤ` whose range is exactly `W` (positivity transfers from `hWpos`). Feed it to
   `LRC14_of_INVcov`.

2. **Reduce the time to `[0,1)`.** Loneliness is invariant under integer shifts, so the new lemma
   **`lonely_fract`** replaces `t` by `Int.fract t ∈ [0,1)`:
   `v_i·fract t − m = v_i·t − (m + v_i·⌊t⌋)`, and `m + v_i⌊t⌋ ∈ ℤ` is absorbed by the universal `m`.
   `Int.fract_nonneg` / `Int.fract_lt_one` place it in `Icc 0 1`.

Then `∀ w ∈ W` is discharged pointwise: `w = e (e.symm ⟨w,hw⟩)`, so `(v i : ℝ) = (w : ℝ)` and the
lonely bound at that index is exactly the target's.

## What is now in the kernel

> **`LRC14_finset_of_INVcov` (corrected S109):
> `LRC(≤13)`[cited] + `INVcov`[open] ⟹ `LRC14.LRC14`** (the ledger target).

Chaining corrected S108's `LRC14_of_INVcov`, the *official* LRC(14) statement — the one written as the ledger's target
Prop, over `Finset ℕ` with the lonely time normalized to `[0,1]` — is now a kernel-checked consequence of
LRC(≤13) and the covering inverse theorem `INVcov`. There is no shape-mismatch caveat between the working
form and the target; the reduction lands on the actual goal Prop.

## The formalization program (S105–S109), summarized

Every rung kernel-pure (`[propext, Classical.choice, Quot.sound]`), all in the root aggregator:

| rung | theorem | content |
|---|---|---|
| S105 | `ap_core_bridge` | ρ≥13 + LRC(≤13) ⟹ Lonely 14 (compact/AP-core) |
| S106 | `sieve_dispatch` | ¬Covering ⟹ Lonely 14 (non-covering sieve) |
| S107 | `density_far_extension` | d≥91V + frame lonely ⟹ Lonely 14 (density Φ>0) |
| S108 | `LRC14_of_INVcov` | LRC(≤13) + INVcov ⟹ ∀v, Lonely 14 (covering M-split) |
| S109 | `LRC14_finset_of_INVcov` | LRC(≤13) + INVcov ⟹ `LRC14.LRC14` (the ledger target) |

The elementary formalization of LRC(14) is complete end-to-end: the ledger's own target Prop reduces, in
the kernel, to LRC(≤13) plus the single open `INVcov`. Connecting this sufficient dominance target to
global n=12 AP uniqueness requires a separate structural bridge.

## Net

- **Delivered:** `LRCFinsetBridge.lean`, 2 kernel-pure theorems (`lonely_fract`, `LRC14_finset_of_INVcov`),
  built and registered. The Lonely-form of LRC(14) now connects to the ledger's Finset target with no
  remaining shape gap.
- **Milestone closed:** the official LRC(14) target Prop is a kernel-checked corollary of
  LRC(≤13) + `INVcov`.
- **Honest:** `INVcov` remains open. The corpus does not prove LRC(14); it proves the
  sufficient reduction to `INVcov` at the level of the official goal statement.

Cross-links:
[[the-M-split-and-the-complete-kernel-checked-reduction-of-lrc14-to-INV-boxeph-S108]],
[[the-density-route-discharge-is-now-kernel-pure-lean-boxeph-S107]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]].
