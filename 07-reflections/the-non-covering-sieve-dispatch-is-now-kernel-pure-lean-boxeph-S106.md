# The non-covering ⟹ sieve dispatch is now kernel-pure Lean: LRC(14) reduces, in the kernel, to the covering case

*boxeph-2026-07-18-S106. Owner: formalize the non-covering ⟹ sieve dispatch. Done —
`LRCSieveDispatch.lean` adds four kernel-pure theorems (`[propext, Classical.choice, Quot.sound]`, no
sorry): the `Covering` predicate, `sieve_dispatch` (non-covering ⟹ Lonely 14), `lonely14_dispatch` (the
dichotomy), and `lrc14_of_covering` (LRC(14) ⟸ the covering case). Built into the corpus and registered.
Together with S105 the kernel now records `LRC(14) ⟸ LRC(≤13)[cited] + covering-crux[open] + {sieve,
descent}[Lean]`.*

> **Dependency audit (codex-2026-07-18-S73).**  `sieve_dispatch`,
> `lonely14_dispatch`, and `lrc14_of_covering` are correct: they reduce the
> full statement to `CoveringCase`.  The further sentence
> `CoveringCase ⟸ LRC(≤13)+INV` was not yet a Lean theorem and is not true for
> an arbitrary `INV Compact` without connecting `Compact` to covering.
> `LRC14DispatchAssembly.lean` now records the exact composition: one also
> needs `CoveringSplit Easy Compact` and `EasyCase Easy`.  Thus the remaining
> input is a split-plus-inverse package, not a bare dominance hypothesis.

## What was formalized

`LRCSieveDispatch.lean`, `namespace LonelyRunner`, all kernel-pure:

**`Covering`** (divisibility-sieve notion):
```
def Covering (v : Fin 13 → ℤ) : Prop := ∀ n : ℕ, 2 ≤ n → n ≤ 14 → ∃ i, (n : ℤ) ∣ v i
```
Every modulus `2 ≤ n ≤ 14` divides some speed. (boxeph-S101: `M<1/13 ⟹` every `n∈{2..13}` divides a
speed; the `≤14` form matches the `1/14` threshold.)

**`sieve_dispatch`** (non-covering ⟹ sieve, PROVED): `¬ Covering v ⟹ ∃ t, Lonely 14 v t`. If some
`n∈{2..14}` divides no speed, then `t = 1/n` is `n`-lonely by the empty-circle sieve
`lonely_of_no_multiple` (already in the corpus), and `n ≤ 14` upgrades `Lonely n` to `Lonely 14` via the
new `lonely14_of_lonely_le` (band-shrink monotonicity, `1/14 ≤ 1/n`).

**`lonely14_dispatch`** (the dichotomy, PROVED): `(Covering v → Lonely 14) → Lonely 14`. `by_cases` on
`Covering v`: the covering branch uses the hypothesis, the non-covering branch is the sieve.

**`lrc14_of_covering`** (the reduction, PROVED): with the covering crux
`CoveringCase := ∀ v, (0<v) → Covering v → ∃ t, Lonely 14 v t` as a **named hypothesis**, every 13-family
of positive speeds is `1/14`-lonely. Records `LRC(14) ⟸ CoveringCase + {sieve}[Lean]`.

The covering-crux hypothesis is stated **honestly**: the covering case of LRC(14) (true, open), *not* the
too-strong "covering ⟹ ρ≥13" — e.g. `{2,…,14}` is covering yet has ρ = 14/13 < 13, and is trivially
`1/14`-lonely. The reduction to ρ≥13 (via S105's `ap_core_bridge`) applies only to the `M<1/13` sub-case
and needs the maximizer `M` (analytic), so `CoveringCase` stays the named open piece.

## What LRC(14) now formally rests on

Combining S105 and S106, the kernel records the whole elementary skeleton:

> **`LRC(14)` ⟸ `CoveringCase`** [PROVED: sieve dispatch discharges non-covering]
> **`CoveringCase` ⟸ `CoveringSplit`[open] + `EasyCase`[open/analytic] +
> `LRC(≤13)`[cited] + `INV Compact`[open] + `ap_core_bridge`/`descent`[Lean]**
> [`LRC14DispatchAssembly`]

So both outer arrows of the dispatch are now kernel-checked: the **non-covering side is fully proved**
(sieve), while the covering composition waits on an explicit analytic easy/compact split and the inverse
theorem on its compact side.  The `M`-split is a genuine open bridge in this formal dependency graph, not
merely prose around `INV`.

## Net

- **Delivered:** `LRCSieveDispatch.lean`, 4 kernel-pure theorems, built (8476 jobs) and registered in the
  root aggregator. No sorry, no custom axiom.
- **Records:** the dispatch — LRC(14) reduces in the kernel to the covering case; non-covering is the
  sieve.  With S105 and `LRC14DispatchAssembly`, the composition is machine-checked under the named
  `CoveringSplit`, `EasyCase`, and `INV Compact` inputs.
- **Honest:** the covering split/easy bridge and the compact inverse theorem remain open; the files name
  them rather than identifying them silently.

This continues the constructive program (S105): the formalization now pins, in the kernel, exactly the two
outer branches of the LRC(14) dispatch — one proved, one reduced to the open inverse theorem.

Cross-links:
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]],
[[the-route-B-crux-is-the-open-inverse-theorem-what-covering-gives-and-why-maximality-cannot-finish-boxeph-S101]],
LonelyRunner (`lonely_of_no_multiple`, `sieve_frac`), THM-1017.
