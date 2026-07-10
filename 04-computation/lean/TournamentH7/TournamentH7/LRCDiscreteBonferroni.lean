/-
  TournamentH7.LRCDiscreteBonferroni — THM-671's HISTOGRAM BOUND in Lean
  (death-star-2026-07-09-S5, HYP-5780; klein-S210's named handoff "Lean next:
  LRCDiscreteBonferroni.lean (decide-shaped, no analysis)").

  THE OBJECT (klein-S210, THM-671).  Fix a modulus `q > 0` and count, for every nonzero
  multiplier `p ∈ (0, q)`, how many runners' residues `(v i · p) mod q` FALL OUTSIDE the safe
  band `[q/14, 13q/14]` — the coverage count `C(p) = bandCount`.  A LIVE multiplier (`C(p) = 0`)
  is an explicit loneliness witness `t = p/q` (kps-S114 `mreach_ge_of_pairsum_band`).  The
  quintic Bonferroni functional

      B5(v, q) = Σ_{d ≤ 5} (−1)^d · S_d,     S_d = Σ_{p ≠ 0} C(C(p), d)

  is computable in ONE pass over the coverage histogram, and `B5 ≤ liveCount` — so `B5 > 0` is a
  DECIDABLE loneliness certificate.  The depth ladder explains the certificate history: at the
  iid heuristic, `B1(13) = −6/7`, `B3(13) = −34/343` (why every union-bound/depth-2 ledger
  broke), `B5(13) = 2052/7⁵ = +0.1221` — depth 5 is the first truncation that clears 13 runners
  at threshold `1/7`.

  THE PROOF SHAPE.  Pointwise engine: mac-mini-S101's `BonferroniTruncation`
  (`odd_truncation_le_uncovered`: for odd `D`, `Σ_{d≤D} (−1)^d C(c,d) ≤ 1{c=0}` — THM-599's
  core, NOT re-proved here).  Sum it over the multipliers, swap the two finite sums, and the
  moment form `B5` appears; positivity then EXHIBITS a live `p`, and kps-S114's pair-sum
  dispatch turns it into `Mreach ≥ 1/14`.

  END-TO-END DEMO.  The AP extremal `{1,…,13}` at its unique live ruler `q = 14` (mac-mini
  THM-668: the AP kills every other ruler): coverage histogram `{0: 6, 1: 6, 6: 1}` — the six
  units of `(ℤ/14)ˣ` are live, the six even multipliers each cover only `v = 7`, and `p = 7`
  covers the six even speeds — so `B5 = 6·1 + 6·0 + 1·(−C(5,5)) = 5 > 0`, and the machine
  certifies `Mreach ≥ 1/14` for the equality extremal by `decide` + the dispatch.

  Kernel-pure target (plain `decide` for the demo; no `native_decide`, no `sorry`).
-/
import Mathlib
import TournamentH7.BonferroniTruncation
import TournamentH7.LRCPairSumDispatch

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Runner `i` passes the safe band at multiplier `p`, modulus `q`:
`(v i · p) mod q ∈ [q/14, 13q/14]`, as exact integer inequalities (the exact
shape `mreach_ge_of_pairsum_band` consumes). -/
def inBand (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13) : Prop :=
  (q : ℤ) ≤ 14 * ((v i * (p : ℤ)) % (q : ℤ)) ∧
    14 * ((v i * (p : ℤ)) % (q : ℤ)) ≤ 13 * (q : ℤ)

instance (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13) : Decidable (inBand v q p i) :=
  inferInstanceAs (Decidable (_ ∧ _))

/-- **The coverage count `C(p)`** (klein-S210): how many runners fail the safe band at
multiplier `p`.  `C(p) = 0` ⟺ `p` is a live multiplier ⟺ `t = p/q` is a loneliness witness. -/
def bandCount (v : Fin 13 → ℤ) (q p : ℕ) : ℕ :=
  (univ.filter fun i : Fin 13 => ¬ inBand v q p i).card

/-- The live-multiplier count `LM(q)`: nonzero multipliers whose coverage vanishes. -/
def liveCount (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 0).card

/-- The `d`-th factorial moment `S_d` of the coverage histogram (one pass). -/
def momentS (v : Fin 13 → ℤ) (q : ℕ) (d : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d

/-- **THM-671's quintic Bonferroni functional** `B5 = S₀ − S₁ + S₂ − S₃ + S₄ − S₅`. -/
def B5 (v : Fin 13 → ℤ) (q : ℕ) : ℤ :=
  ∑ d ∈ range 6, (-1 : ℤ) ^ d * (momentS v q d : ℤ)

/-- **The histogram bound (THM-671, klein-S210): `B5 ≤ LM(q)`.**  Sum mac-mini-S101's
pointwise odd-depth truncation over the multipliers and swap the finite sums. -/
theorem B5_le_liveCount (v : Fin 13 → ℤ) (q : ℕ) : B5 v q ≤ (liveCount v q : ℤ) := by
  have hswap : B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold B5 momentS
    calc ∑ d ∈ range 6, (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range 6, ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = ∑ p ∈ Finset.Ioo 0 q,
            ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) :=
          Finset.sum_comm
  rw [hswap]
  have hpoint : ∀ p ∈ Finset.Ioo 0 q,
      ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
        ≤ (if bandCount v q p = 0 then 1 else 0) := fun p _ =>
    TournamentH7.BonferroniTruncation.odd_truncation_le_uncovered (bandCount v q p) 5 rfl
  calc ∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
      ≤ ∑ p ∈ Finset.Ioo 0 q, (if bandCount v q p = 0 then (1 : ℤ) else 0) :=
        Finset.sum_le_sum hpoint
    _ = (liveCount v q : ℤ) := by
        unfold liveCount
        rw [Finset.sum_boole]

/-- **Positivity exhibits a live multiplier.** -/
theorem exists_live_of_B5_pos (v : Fin 13 → ℤ) (q : ℕ) (h : 0 < B5 v q) :
    ∃ p ∈ Finset.Ioo 0 q, bandCount v q p = 0 := by
  by_contra hno
  push_neg at hno
  have hzero : liveCount v q = 0 := by
    unfold liveCount
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    exact hno
  have hle := B5_le_liveCount v q
  rw [hzero] at hle
  omega

/-- **The decidable loneliness certificate (THM-671 → THM-668 dispatch): `B5(v, q) > 0` ⟹
`Mreach v ≥ 1/14`.**  A positive quintic Bonferroni functional at ANY modulus `q > 0` exhibits a
live multiplier `p`, and the rational `t = p/q` clears every runner (kps-S114).  Everything to
the left of `Mreach` is integer arithmetic — `decide`-shaped, no analysis. -/
theorem mreach_ge_of_B5_pos (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q) (h : 0 < B5 v q) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  obtain ⟨p, _, hlive⟩ := exists_live_of_B5_pos v q h
  apply mreach_ge_of_pairsum_band v (p : ℤ) (q : ℤ) (by exact_mod_cast hq)
  intro i
  have hempty : (univ.filter fun i : Fin 13 => ¬ inBand v q p i) = ∅ :=
    Finset.card_eq_zero.mp hlive
  have hi : ¬ ¬ inBand v q p i := by
    intro hni
    have hmem : i ∈ (univ.filter fun i : Fin 13 => ¬ inBand v q p i) :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ i, hni⟩
    rw [hempty] at hmem
    exact absurd hmem (Finset.notMem_empty i)
  exact not_not.mp hi

/-! ## End-to-end demo: the AP extremal at its unique live ruler -/

/-- The LRC(14) equality extremal `{1, …, 13}`. -/
def apSpeeds : Fin 13 → ℤ := fun i => (i : ℕ) + 1

/-- `B5(AP, 14) = 5 > 0`: the coverage histogram at the AP's unique live ruler `q = 14` is
`{0: 6, 1: 6, 6: 1}` (units live; even multipliers cover only `v = 7`; `p = 7` covers the six
even speeds), so `B5 = 6 − C(5,5) = 5`.  Kernel `decide` — 13 multipliers × 13 runners. -/
theorem b5_ap_fourteen_pos : 0 < B5 apSpeeds 14 := by decide

/-- **The histogram machine certifies the equality extremal**: `Mreach(AP) ≥ 1/14`, recovering
kps-S110's `mreach_AP_ge` through the THM-671 → THM-668 pipe with a `decide` certificate. -/
theorem mreach_AP_via_histogram : (1 : ℝ) / 14 ≤ Mreach apSpeeds :=
  mreach_ge_of_B5_pos apSpeeds 14 (by norm_num) b5_ap_fourteen_pos

/-! ## Axiom audit -/
#print axioms B5_le_liveCount
#print axioms exists_live_of_B5_pos
#print axioms mreach_ge_of_B5_pos
#print axioms b5_ap_fourteen_pos
#print axioms mreach_AP_via_histogram

end LRC14Concrete
end LonelyRunner
