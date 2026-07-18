# The elementary half of LRC(14) is now kernel-pure Lean: the AP-core bridge, and exactly what LRC(14) rests on

*boxeph-2026-07-18-S105. Owner: Lean-formalize the completed elementary half. Done — `LRCAPCoreBridge.lean`
adds three kernel-pure theorems (`[propext, Classical.choice, Quot.sound]`, no sorry, no custom axioms):
`ap_core_bridge` (THM-1017(II): ρ≥13 + LRC(≤13) ⟹ Lonely 14), `ap_core_bridge_of_shape` (the explicit
dilated-AP + `lcm(13,14)` far-element mechanism), and `lonely14_of_INV` (records
`LRC(14) ⟸ LRC(≤13)[cited] + INV[open] + {sieve,descent}[proved]`). Built into the corpus (8478 jobs).*

## What was formalized

`LRCAPCoreBridge.lean`, three theorems in `namespace LonelyRunner`, all kernel-pure:

**1. `ap_core_bridge` — the elementary half's payoff (THM-1017(II)).**
```
theorem ap_core_bridge (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (vstar : Fin 13)
    (hdom : ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t : ℝ, Lonely 14 v t
```
If the largest speed `vstar` dominates 13-fold (ρ = v_max/v_2nd ≥ 13), the family is 1/14-lonely. Proof:
reindex the 12 non-max speeds through `vstar.succAbove : Fin 12 ↪ Fin 13`; the **LRC(≤13) citation**
(`LRCUpTo13`, klein's node) makes them 1/13-lonely at some `t₀`; the **descent floor** `descent_dominant`
(THM-1008, already in the corpus) lifts that to a 1/14-lonely time. The proof is ~10 lines: the only
subtleties are the `Fin.exists_succAbove_eq` reindexing and the `(1:ℝ)/(12+1) = 1/13` cast.

**2. `ap_core_bridge_of_shape` — the THM-1017 mechanism made explicit.**
```
theorem ap_core_bridge_of_shape (... ) (d : ℤ)
    (hcore : ∀ i, i ≠ vstar → v i ≤ 12 * d) (hfar : 156 * d ≤ v vstar) :
    ∃ t : ℝ, Lonely 14 v t
```
A dilated-AP core `d·{1,…,12}` has `v_i ≤ 12d`, and the inverse theorem forces the far element to be a
multiple of `lcm(13,14)·d = 182d ≥ 156d = 13·(12d)`, so `13·v_i ≤ 156d ≤ v_max` — the dominance holds.
This records the deep-well arithmetic (`182 = 13·14`, boxeph-S103) in Lean.

**3. `lonely14_of_INV` — the reduction, recorded.**
```
def INV (Compact : (Fin 13 → ℤ) → Prop) : Prop :=
  ∀ v, (∀ i, 0 < v i) → Compact v → ∃ vstar, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar
theorem lonely14_of_INV (cite : LRCUpTo13) (inv : INV Compact)
    (v ...) (hc : Compact v) : ∃ t, Lonely 14 v t
```
`INV` is the inverse theorem in **dominance (ρ≥13) form** — the ONE open piece (≡ LRC(14) covering crux ≡
Tao n=12, boxeph-S94/S104) — entered as a **named hypothesis, never a `sorry`** (exactly the CLAUDE.md
policy for LRC(≤13)). The composition is one line: `INV` supplies the dominance, `ap_core_bridge` supplies
the loneliness.

## What LRC(14) now formally rests on

The corpus now makes the dependency explicit and machine-checked:

> **`LRC(14)` ⟸ `LRC(≤13)` [cited] + `INV` [open] + { `sieve_frac`, `fill1_perturbation`,
> `descent_general/dominant`, `dilated_sieve` } [all kernel-pure Lean].**

Every arrow except `INV` is discharged: the elementary witness/descent family is kernel-pure (fill-1,
descent, dilated sieve, and now the AP-core bridge assembling them with the citation), and the density
route is discharged for separated far elements (S96–S100, analytic — not yet in Lean). `INV` is the single
remaining hypothesis, and it is precisely the additive inverse theorem that S101–S104 showed is beyond the
elementary toolkit (maximality / sieve / CF descent) and the additive toolkit (BSG/PFR needs the energy
input `M<1/13` doesn't supply). The Lean file names it and records that nothing else is missing on the
elementary side.

## Net

- **Delivered:** `LRCAPCoreBridge.lean`, 3 kernel-pure theorems, built into the corpus (`lake build
  TournamentH7`, 8478 jobs; `#print axioms` = the standard trio on all three). Registered in the root
  aggregator.
- **Records:** the exact reduction `LRC(14) ⟸ LRC(≤13) + INV + {elementary, all in Lean}`, with `INV`
  a named open hypothesis. The elementary half of LRC(14) is now formally complete and machine-checked.
- **Honest:** `INV` is open (= LRC(14) covering crux). The Lean file does not close it; it certifies that
  the *only* thing between the current corpus and LRC(14) is that single additive inverse theorem.

This is the constructive counterpart to S101–S104: having shown the crux is the irreducible open core, the
formalization now pins, in the kernel, exactly where LRC(14) stands and what it waits on.

Cross-links:
[[bsg-pfr-attack-the-wrong-half-the-crux-is-the-diophantine-to-energy-bridge-boxeph-S104]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-1017 (AP-core bridge), THM-1008/1010 (descent floor), LRC13Citation (`LRCUpTo13`).
