/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S199)
-/
import TournamentH7.LRCResidualFromLedger

/-!
# The LRC(14) completion audit — the machine-checked "how complete is it" certificate (opus-S199)

This file is the up-to-date axiom-footprint audit of the LRC(14) finish line. It supersedes the coverage
of `LRC14AxiomAudit.lean` (klein-S113, which audited the older `lrc14_of_covering_far_22` surface): the
current tightest top theorem is `LonelyRunner.LRC14Grand.lrc14_from_B5`.

## The certified state (owner directive: "get the entire LRC(14) proof completely formalized")

`lrc14_from_B5 : LRCUpTo13 → hB5 → LRC14.LRC14Statement`, where

* `cite : LRCUpTo13` — the owner-sanctioned LRC(≤13) citation node (a HYPOTHESIS, discharged by citation
  per the 2026-07-02 owner directive — never by `sorry`);
* `hB5` — the SINGLE remaining obligation: every residual-class family (covering, scale-gapped,
  compressed, distinct-|·|, some `|vᵢ| ≥ 23`, divisor-closed, not coarse-≤12-reducible) has a pair-sum
  ruler `q` with positive depth-5 discrete Bonferroni bound `0 < B5 v q`.

The `#print axioms` below is the machine's verdict. It reads:

    'lrc14_from_B5' depends on axioms:
      [propext, Classical.choice, Quot.sound,                       -- Lean foundations (always allowed)
       WindowData.winData22_complete._native.native_decide.ax_1_1,  -- the window-≤22 census (native_decide)
       WindowData.winData22_ok._native.native_decide.ax_1_1]

with **NO `sorryAx`** and the two named obligations taken as HYPOTHESES, not axioms. So:

> **LRC(14) is fully formalized and kernel-sound, sorry-free, modulo (i) the owner-sanctioned LRC(≤13)
> citation and (ii) the single analytic obligation `hB5`.** Everything else — the q-witness sieve, the
> five grand-assembly branches (non-covering ∣ ratio-≤13 ∣ dominant-peel ∣ window-≤22 ∣ repeated-|speeds|),
> `M(AP) = 1/14` exact, the consumer chain `B5 > 0 ⟹ Mreach ≥ 1/14 ⟹ ∃ lonely` — is discharged.

`hB5` is the fleet's open ANALYTIC frontier (the residual-family liveness / signed diagonal-suppression
gate — confirmed empirically over 849 covering sets but not yet proved a-priori; klein/monad/kps active).
It is a genuine number-theoretic ingredient, NOT a formalization gap: no amount of Lean wiring closes it.

The only trusted extra-foundational axioms are the two `native_decide` window-census facts; a purity
refinement (replacing them by kernel `decide`) is possible but orthogonal to closing `hB5`.
-/

open LonelyRunner LonelyRunner.LRC14Grand LonelyRunner.LRC14

/-- The tightest LRC(14) finish-line theorem restated as the audit's subject: `LRC14Statement` from the
LRC(≤13) citation and the single `B5 > 0` residual obligation (definitionally `lrc14_from_B5`). -/
theorem lrc14_reduced_to_the_B5_obligation
    (cite : LRCUpTo13)
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    LRC14.LRC14Statement :=
  lrc14_from_B5 cite hB5

-- The machine-checked footprint: Lean foundations + the two window-22 native_decide facts, NO sorryAx.
#print axioms lrc14_reduced_to_the_B5_obligation
#print axioms lrc14_from_B5
