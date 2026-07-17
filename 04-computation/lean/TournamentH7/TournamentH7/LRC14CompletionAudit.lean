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

The `#print axioms` below is the machine's verdict. As of **opus-S202 it reads**:

    'lrc14_from_B5' depends on axioms: [propext, Classical.choice, Quot.sound]

— **Lean's foundational axioms ONLY**, with **NO `sorryAx`** and **NO `native_decide`**, and the two named
obligations taken as HYPOTHESES, not axioms. So:

> **LRC(14) is fully formalized, kernel-sound, sorry-free, and FOUNDATIONAL-AXIOMS-ONLY, modulo (i) the
> owner-sanctioned LRC(≤13) citation and (ii) the single analytic obligation `hB5`.** Everything else — the
> q-witness sieve, the grand-assembly branches (non-covering ∣ ratio-≤13 ∣ dominant-peel ∣ window-≤22 ∣
> repeated-|speeds| ∣ detuned ∣ coarse ∣ common-residue), `M(AP) = 1/14` exact, and the consumer chain
> `B5 > 0 ⟹ Mreach ≥ 1/14 ⟹ ∃ lonely` — is discharged.

`hB5` is the fleet's open ANALYTIC frontier (the residual-family liveness / signed diagonal-suppression
gate — confirmed empirically over 849 covering sets but not yet proved a-priori; klein/monad/kps active).
It is a genuine number-theoretic ingredient, NOT a formalization gap: no amount of Lean wiring closes it.

## How the last two axioms were removed (opus-S200 → S202)

The window-≤22 branch used to rest on two `native_decide` census facts (`winData22_ok`,
`winData22_complete`). Kernel `decide` on them is **infeasible** (MISTAKE-135: `winData22_complete` ranges
over `C(22,13) = 497 420` candidates; measured >13 h just to generate the sublists, plus OOM). And THM-665
does **not** shrink that census — its own Consequence 2 proves the a-priori/density-floor route never fires
on covering clusters, which *is* the census case (opus-S201).

The honest route was a **mathematical proof of the census** — delivered as **LEM-024** and formalized in
`LRCWindow22Census.lean` (kernel-pure): after `spread13` peels `ratio ≤ 13`, the census's necessary domain
is the `min = 1` covering ≤22 tuples, and *every* such tuple is far at one of **six** fixed witnesses
`{12/25, 9/26, 7/27, 11/28, 4/23, 11/26}` — whose danger sets over `[1,22]` are the tiny
`{2},{3},{4},{5},{6,17},{7,19}`. Failing all six forces `{1,2,3,4,5}` + one of `{6,17}` + one of `{7,19}`;
covering forces `{12,13,14}` + one each of `{8,16},{9,18},{10,20},{11,22}` — **14 pairwise-distinct
elements in a 13-set**, a contradiction. So one witness always works and
`KernelGate.lonely_of_kernelWitness` finishes. No enumeration of tuples ever occurs; the only `decide`s are
the 22 × 6 tiny per-speed `speedOK` facts.

The 1-line assembly swap `hwindow22_closed cite` → `hwindowW_closed 22 cite hdistinct22_kernel`
(`LRC14GrandAssembly.lean`, window branch) retired both axioms. `lrc14_grand_assembly_pure` (which bought
foundational purity by *enlarging* the residual obligation) is now redundant as a purity route: the main
`lrc14_grand_assembly` is foundational-only **and** keeps the smaller residual (`Vmax ≥ 23`).
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

-- The machine-checked footprint: Lean foundations only, with NO native_decide and NO sorryAx.
#print axioms lrc14_reduced_to_the_B5_obligation
#print axioms lrc14_from_B5
