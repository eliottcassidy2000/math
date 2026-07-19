# Corrected M-split: residual target exposed; historical global-INV claim withdrawn

> **CORRECTION (codex-S67, 2026-07-18).**  The universal `INV` described below
> is known false, not the open Tao `n=12` inverse.  The row `{1,...,13}` has
> exact maximum `M=1/14`, hence no `Lonely 13` time, but has no 13-dominant
> speed.  The original Lean theorem was a kernel-valid implication from this
> over-strong premise.  `LRCMSplit.lean` now exposes `ResidualINV`: positivity,
> `Covering`, and the actual outer `not exists Lonely 14` branch imply
> dominance.  Its consumer is `LRC14_of_residual_INV`.  This is an honest exact
> residual target and, with the AP bridge, is equivalent to LRC(14); it is not
> a reduction to Tao `n=12` without a separate theorem.  The historical account
> below is retained to show precisely where the quantifier error entered.

*boxeph-2026-07-18-S108; corrected codex-2026-07-18 (MISTAKE-166). `LRCMSplit.lean`
now culminates in `LRC14_of_INVcov`:
**LRC(14) ⟸ LRC(≤13)[cited] + INVcov[open]**, kernel-checked. Here `INVcov` explicitly
assumes divisibility coverage through 14. The original `INV` omitted that premise and was false;
Lean had proved a valid implication from an overstrong premise, not the intended reduction target.*

## The idea: M-thresholds are loneliness propositions

The obstacle to formalizing an "M-split" was that `M(V) = sup_t min_i ‖v_i·t‖` is a supremum (needs
compactness to be achieved). The realization that removes it entirely:

> **`M(V) < 1/n  ⟺  ¬ ∃ t, Lonely n V t`.**

No sup, no max, no achievement lemma — the threshold `M < 1/n` *is* the proposition "no `n`-lonely time".
So the honest split is "`M ≥ 1/14` (immediate) vs
`Covering(2..14) ∧ M < 1/13` (crux)". It is `Classical` case analysis on
`∃ t, Lonely 14 V t`: in the negative branch, `counterexample_needs_all_divisors 14`
supplies coverage through 14 and band-shrink monotonicity supplies no `Lonely 13` time.

## What was formalized

`LRCMSplit.lean`, `namespace LonelyRunner`, all kernel-pure:

**`M_split`** (PROVED):
`(Covering v → ¬∃ Lonely13 → ∃ Lonely14) → ∃ Lonely14`. In the negative branch it derives
`Covering v` from the 14-level divisor sieve and derives `¬∃ Lonely13` by band monotonicity.

**`crux_of_dominance`** (PROVED): the covering crux follows from
`Covering(2..14) ∧ M<1/13 ⟹ ρ≥13` plus `ap_core_bridge` (S105).

**`lonely14_of_dominance`** / **`LRC14_of_INVcov`** (PROVED — corrected capstone):
```
theorem LRC14_of_INVcov (cite : LRCUpTo13)
    (inv : ∀ v, (∀ i, 0 < v i) → Covering v → (¬ ∃ t, Lonely 13 v t) →
             ∃ vstar, ∀ i ≠ vstar, 13 * v i ≤ v vstar) :
    ∀ v, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t
```
**Every 13-family of positive speeds is `1/14`-lonely, given LRC(≤13) and the single inverse theorem
`INVcov`.** The M-split reduces to the covering crux; the crux to `INVcov + ap_core_bridge`.

The covering premise is essential. `¬∃ Lonely13` yields only coverage through 13. The family
`{1,2,3,5,7,8,9,10,11,12,17,19,104}` has `M=8/105<1/13`, misses 14, and has
`ρ=104/19<13`, so it refutes the original unqualified `INV`.

## The formalization program, complete (S105–S108)

The Lean corpus now reduces LRC(14) to one open theorem, entirely in the kernel:

> **`LRC14_of_INVcov`  (corrected S108):
> LRC(14) ⟸ LRC(≤13)[cited] + INVcov[open].**

with every supporting bridge kernel-pure:
- **`ap_core_bridge`** (S105): `ρ≥13 + LRC(≤13) ⟹ Lonely 14` — the descent discharge of the compact case.
- **`sieve_dispatch`** (S106): `¬Covering ⟹ Lonely 14` — the non-covering sieve.
- **`density_far_extension`** (S107): `d≥91V + frame lonely ⟹ Lonely 14` — the density Φ>0 discharge.
- **`M_split`** (S108): derives `Covering(2..14)` in the no-`Lonely14` branch and
  reduces the whole thing to the covering `M<1/13` crux, no `M` sup needed.

Three independent proofs of the far-element/covering discharges (descent, sieve, density), and the M-split
assembling them, certify the sufficient reduction to **`INVcov`**. Relating `INVcov` to the
global n=12 AP-uniqueness/sporadic-emptiness program still requires the appropriate structural bridge;
the two open statements should not be identified without that bridge.

## Net

- **Delivered:** `LRCMSplit.lean`, the covering M-split via `M<1/n ⟺ ¬∃Lonely n`
  (no sup), and `LRC14_of_INVcov` — the kernel-checked reduction of LRC(14) to
  LRC(≤13) + the single open covering inverse theorem.
- **Milestone:** the elementary formalization of LRC(14) is now **assembled end-to-end**. What LRC(14)
  rests on is no longer a prose claim but a single Lean theorem with one named open hypothesis.
- **Honest:** `INVcov` is open and is a sufficient structural target, not an asserted equivalence.
  The corpus does not close LRC(14).

Cross-links:
[[the-density-route-discharge-is-now-kernel-pure-lean-boxeph-S107]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]].
