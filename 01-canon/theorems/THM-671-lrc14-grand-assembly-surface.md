---
id: THM-671
title: The LRC(14) grand assembly surface (v3) — LRC14Statement follows in Lean from the LRC(≤13) citation plus ONE residual obligation (covering ∧ scale-gapped ∧ compressed ∧ distinct-|speeds| ∧ max ≥ 23 ∧ NO detuned-harmonic structure ∧ NO multi-scale decomposition), with SEVEN branches discharged by the sorry-free corpus (incl. THM-668 formalized in LRCDetunedDispatch.lean, kernel-pure); plus the machine-checked completeness of the 966-witness [1,18] base case
status: PROVED/BUILT (Lean: TournamentH7/LRC14GrandAssembly.lean, build 8506 jobs). TWO variants: lrc14_grand_assembly (sharpest surface; axioms = [propext, Classical.choice, Quot.sound] + the two winData22 native_decide certificate axioms, inherited from the window-22 branch) and lrc14_grand_assembly_pure (no window branch; KERNEL-PURE [propext, Classical.choice, Quot.sound]; residual additionally contains max ≤ 22 families). covering18_complete carries its own native_decide axiom (8568-case sweep). The residual obligations are the honest open surface, NOT asserted.
source: monad-explorer-2026-07-09-S6 (HYP-5757).
depends_on:
  - LRCSpread13 (kps-S28), LRCEndgameAssembly/LRCDominantPeel (opus), LRCWindowData22 (window ≤ 22, cite-conditional), LRC13Citation, LRCCovering966 (kps-S115), LonelyRunner sieve
related:
  - THM-663   # the covering-case prose closure this Lean surface mirrors
  - THM-665/667/668/669/670  # the analytic program aimed at the residual class
---

# THM-671 — the LRC(14) grand assembly surface

## The theorem (Lean, built, kernel-pure)

```
theorem lrc14_grand_assembly (cite : LRCUpTo13) (hresidual : ResidualObligation) :
    LRC14.LRC14Statement
```

where `ResidualObligation` is the single remaining obligation:

> every 13-family that is **covering** (all q ∈ {2,…,14} divide some speed), **scale-gapped**
> (some pair of |speeds| beyond factor 13), **compressed** (no runner exceeds 13× all others),
> with **pairwise-distinct |speeds|**, and **max |speed| ≥ 23**, is lonely.

Discharged internally: (1) non-covering → `sieve_one_div` (unconditional); (2) no scale
gap → `spread13_lonely` (unconditional); (3) dominant → `hdom_discharged` (cite);
(4) window ≤ 22 → `hwindow22_closed` (cite; sharp variant only); (5) repeated |speed| →
`lonely_of_abs` + `lonely14_of_repeat` (cite); (6) DETUNED HARMONIC (one g ≥ 2 divides
all but one speed) → `DetunedDispatch.lonely14_of_detuned` — THM-668 FORMALIZED here
(kernel-pure: the quarter-window + Bezout branch shift + triangle pigeonhole); (7)
MULTI-SCALE (a coarse decomposition v = a + L·k with ≤ 12 distinct coarse values and
budget A/L ≤ 1/182) → `CoarseReduction.lonely14_of_coarse_le12` (cite).

## Why this is the sharpest surface in the corpus

Prior best: `lrc14_of_compressed` (residual = covering ∧ compressed) and `lrc14_endgame`
(opaque `witnessG2` hypotheses). This surface adds the scale-gap, distinct-|speeds|, and
max ≥ 23 carve-outs — each removing an infinite family class from the obligation — and
consumes only the project-sanctioned citation. The five branches use both unconditional
witnesses (sieve, spread13) and cite-conditional machinery (window22, dominant peel,
repeat), matching the owner's LRC(≤13) policy.

## The [1,18] base case, pinned end-to-end

`covering18_complete` (one native_decide over the C(18,13) = 8568 subsets): every
primitive covering 13-subset of [1,18] appears in kps-S115's `coveringWitnesses` list;
with `coveringWitnesses_lonely` (all 966 have Mreach ≥ 1/14) and `lonely_comp_perm`
(permutation invariance, proved here), the [1,18] slice is machine-closed. (The
window-22 branch subsumes it in the assembly; the completeness theorem additionally
certifies kps's list itself.)

## What remains (the honest surface)

`ResidualObligation` — the class the analytic program (THM-665/667/668/669/670 +
klein-S205's drift embed + the certificates C0–C3) attacks. Its known-closed slices
(detuned harmonics via THM-668; the composition-certified batteries) are not yet
wired into Lean quantified form; they shrink the *mathematical* residue, not yet the
*Lean* residue. The remaining formalization work is exactly: formalize the analytic
program's uniform statements over this one named class.
