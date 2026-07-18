# The M-split, and the complete kernel-checked reduction of LRC(14) to a single inverse theorem

*boxeph-2026-07-18-S108. Owner: formalize the M-split separating M≥1/14 from M<1/13. Done — and it
completes the assembly. `LRCMSplit.lean` adds four kernel-pure theorems (`[propext, Classical.choice,
Quot.sound]`, no sorry), culminating in `LRC14_of_INV`: **LRC(14) ⟸ LRC(≤13)[cited] + INV[open], fully
kernel-checked.** The maximizer `M` never needed to be defined as a supremum. Built into the corpus
(8480 jobs).*

## The idea: M-thresholds are loneliness propositions

The obstacle to formalizing an "M-split" was that `M(V) = sup_t min_i ‖v_i·t‖` is a supremum (needs
compactness to be achieved). The realization that removes it entirely:

> **`M(V) < 1/n  ⟺  ¬ ∃ t, Lonely n V t`.**

No sup, no max, no achievement lemma — the threshold `M < 1/n` *is* the proposition "no `n`-lonely time".
So the split "`M ≥ 1/14` (immediate) vs `M < 1/13` (crux)" is nothing but `Classical` case analysis on
`∃ t, Lonely 14 V t`, glued by the band-shrink monotonicity `Lonely 13 ⟹ Lonely 14`.

## What was formalized

`LRCMSplit.lean`, `namespace LonelyRunner`, all kernel-pure:

**`M_split`** (PROVED): `(¬∃t, Lonely 13 v t) → ∃t, Lonely 14 v t) → ∃t, Lonely 14 v t`. To prove a family
is `1/14`-lonely it suffices to handle the `M<1/13` sub-case. `by_cases` on `∃t Lonely 14`: if yes, done;
if no, then no `1/13`-lonely time either (a `1/13`-lonely time would be `1/14`-lonely, since `1/14<1/13`),
so the crux applies. Four lines.

**`crux_of_dominance`** (PROVED): the crux `¬∃Lonely13 → ∃Lonely14` follows from the inverse theorem in
dominance form (`M<1/13 ⟹ ρ≥13`) plus `ap_core_bridge` (S105).

**`lonely14_of_dominance`** / **`LRC14_of_INV`** (PROVED — the capstone):
```
theorem LRC14_of_INV (cite : LRCUpTo13)
    (INV : ∀ v, (∀ i, 0 < v i) → (¬ ∃ t, Lonely 13 v t) →
             ∃ vstar, ∀ i ≠ vstar, 13 * v i ≤ v vstar) :
    ∀ v, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t
```
**Every 13-family of positive speeds is `1/14`-lonely, given LRC(≤13) and the single inverse theorem
`INV`.** The M-split reduces to the crux; the crux to `INV + ap_core_bridge`. Non-covering families are
absorbed for free: `¬∃Lonely13` already entails covering (a non-covering family has a sieve witness that
is `1/13`-lonely), so no separate sieve branch is needed in this assembly.

`INV` here is the honest crux in its cleanest form: `M<1/13 ⟹ ρ≥13` (dominance). It is the LRC(14)
covering crux (≡ Tao n=12, S94/S104), entered as a named hypothesis, never a `sorry`.

## The formalization program, complete (S105–S108)

The Lean corpus now reduces LRC(14) to one open theorem, entirely in the kernel:

> **`LRC14_of_INV`  (S108):  LRC(14)  ⟸  LRC(≤13)[cited]  +  INV[open].**

with every supporting bridge kernel-pure:
- **`ap_core_bridge`** (S105): `ρ≥13 + LRC(≤13) ⟹ Lonely 14` — the descent discharge of the compact case.
- **`sieve_dispatch`** (S106): `¬Covering ⟹ Lonely 14` — the non-covering sieve.
- **`density_far_extension`** (S107): `d≥91V + frame lonely ⟹ Lonely 14` — the density Φ>0 discharge.
- **`M_split`** (S108): reduces the whole thing to the `M<1/13` crux, no `M` sup needed.

Three independent proofs of the far-element/covering discharges (descent, sieve, density), and the M-split
assembling them, all certify: **the only thing between the current Lean corpus and LRC(14) is `INV`** — the
additive inverse theorem that S101–S104 proved is beyond the elementary and additive toolkits (= Tao n=12,
open).

## Net

- **Delivered:** `LRCMSplit.lean`, 4 kernel-pure theorems. The M-split via `M<1/n ⟺ ¬∃Lonely n` (no sup),
  and `LRC14_of_INV` — the complete, kernel-checked reduction of LRC(14) to LRC(≤13) + the single open
  inverse theorem.
- **Milestone:** the elementary formalization of LRC(14) is now **assembled end-to-end**. What LRC(14)
  rests on is no longer a prose claim but a single Lean theorem with one named open hypothesis.
- **Honest:** `INV` is open. The corpus does not close LRC(14); it certifies, in the kernel, that LRC(14)
  is *exactly* LRC(≤13) plus the inverse theorem — the constructive terminus of the S96–S107 program.

Cross-links:
[[the-density-route-discharge-is-now-kernel-pure-lean-boxeph-S107]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]].
