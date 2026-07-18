# The official Finset target of LRC(14) now reduces to INV in the kernel

*boxeph-2026-07-18-S109. Owner: connect `LRC14_of_INV`'s Lonely-form to the ledger's Finset target.
Done — `LRCFinsetBridge.lean` proves `LRC14_finset_of_INV : LRC(≤13)[cited] + INV[open] ⟹ LRC14.LRC14`,
kernel-pure (`[propext, Classical.choice, Quot.sound]`, no sorry, built first try). So the **official**
LRC(14) statement of `LRC14Ledger` — `Finset ℕ`, card 13, lonely time in `[0,1]` — now reduces, in the
kernel, to LRC(≤13) plus the single open inverse theorem.*

## The gap bridged

`LRC14_of_INV` (S108) concludes in the working shape:
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
   `LRC14_of_INV`.

2. **Reduce the time to `[0,1)`.** Loneliness is invariant under integer shifts, so the new lemma
   **`lonely_fract`** replaces `t` by `Int.fract t ∈ [0,1)`:
   `v_i·fract t − m = v_i·t − (m + v_i·⌊t⌋)`, and `m + v_i⌊t⌋ ∈ ℤ` is absorbed by the universal `m`.
   `Int.fract_nonneg` / `Int.fract_lt_one` place it in `Icc 0 1`.

Then `∀ w ∈ W` is discharged pointwise: `w = e (e.symm ⟨w,hw⟩)`, so `(v i : ℝ) = (w : ℝ)` and the
lonely bound at that index is exactly the target's.

## What is now in the kernel

> **`LRC14_finset_of_INV`  (S109):  `LRC(≤13)`[cited] + `INV`[open]  ⟹  `LRC14.LRC14`** (the ledger target).

Chaining S108's `LRC14_of_INV`, the *official* LRC(14) statement — the one written as the ledger's target
Prop, over `Finset ℕ` with the lonely time normalized to `[0,1]` — is now a kernel-checked consequence of
LRC(≤13) and the single inverse theorem `INV`. There is no shape-mismatch caveat left between the working
form and the target; the reduction lands on the actual goal Prop.

## The formalization program (S105–S109), summarized

Every rung kernel-pure (`[propext, Classical.choice, Quot.sound]`), all in the root aggregator:

| rung | theorem | content |
|---|---|---|
| S105 | `ap_core_bridge` | ρ≥13 + LRC(≤13) ⟹ Lonely 14 (compact/AP-core) |
| S106 | `sieve_dispatch` | ¬Covering ⟹ Lonely 14 (non-covering sieve) |
| S107 | `density_far_extension` | d≥91V + frame lonely ⟹ Lonely 14 (density Φ>0) |
| S108 | `LRC14_of_INV` | LRC(≤13) + INV ⟹ ∀v, Lonely 14 (M-split assembly) |
| S109 | `LRC14_finset_of_INV` | LRC(≤13) + INV ⟹ `LRC14.LRC14` (the ledger target) |

The elementary formalization of LRC(14) is complete end-to-end: the ledger's own target Prop reduces, in
the kernel, to LRC(≤13) plus the single open `INV` (= the additive inverse theorem = Tao n=12, which
S101–S104 showed is beyond the elementary and additive toolkits).

## Net

- **Delivered:** `LRCFinsetBridge.lean`, 2 kernel-pure theorems (`lonely_fract`, `LRC14_finset_of_INV`),
  built and registered. The Lonely-form of LRC(14) now connects to the ledger's Finset target with no
  remaining shape gap.
- **Milestone closed:** the official LRC(14) target Prop is a kernel-checked corollary of LRC(≤13) + INV.
  The whole S105–S109 program certifies, at the target level, exactly what LRC(14) rests on.
- **Honest:** `INV` remains open (the research crux). The corpus does not prove LRC(14); it proves the
  reduction to `INV` at the level of the official goal statement.

Cross-links:
[[the-M-split-and-the-complete-kernel-checked-reduction-of-lrc14-to-INV-boxeph-S108]],
[[the-density-route-discharge-is-now-kernel-pure-lean-boxeph-S107]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]].
