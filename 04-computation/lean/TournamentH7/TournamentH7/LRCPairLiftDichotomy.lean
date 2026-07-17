/-
  TournamentH7.LRCPairLiftDichotomy — THE PAIR LIFT of the chain-split dichotomy
  (death-star-2026-07-16-S34, HYP-7171; lifts THM-937 through kps-S20's pair dodge).

  THM-937's dense core is a family whose LAST sub-3 ratio pair blocks every singles
  chain.  kps-S20 (`LRCPairBlock`) proved the PAIR DODGE: a near-equal pair
  `w ≤ w' < 3w` can serve as ONE chain level (entry `w·L ≥ 13/7`, output window
  `1/(14·w')`).  This file inserts exactly that pair AT the dense position and re-runs
  the split search through `cite_blockchain_lonely`:

  * `bChainOK_of_chainOK` — a singles `ChainOK` ledger is a `BChainOK` ledger on
    mapped single levels (the S19 arithmetic reused verbatim).
  * `lonely_of_pair_split` — the S20 engine consumed at a sorted split whose first
    level is the dense pair and whose tail is the mapped singles drop.
  * `lonely_or_tripleCore` — THE LIFTED DICHOTOMY: every sorted positive 13-family is
    lonely (cite + a singles split, OR cite + the pair-at-the-dense-position split)
    or lies in `TripleDenseCore`: the THM-937 certificate PLUS
      (the compressed triple  `w(j+2) < 21·w(j+1)`,  or
       the pair-entry failure `(13−j)·w(j) < 13·(j+1)·w(j−1)` at `j ≥ 1`).
    Dense cores at `j = 11` (pair on top, nothing above) close OUTRIGHT when the
    pair entry fee holds; at `j = 0` the pair entry is FREE (`w₀ ≥ 1`).
  * `residualObligation_of_tripleCore` / `lrc14_of_tripleCore` — the wire into the
    grand assembly: strictly sharper than THM-937's `DenseCoreObligation`.

  Transition fees behind the normal form: single→single `3`, pair→single `21`
  (the pair's output window is 7× tighter), pair entry `13(k+1)/(13−k)` over the
  cited base (vs `21(k+1)/(2(13−k))` for a single entry — the pair enters CHEAPER
  by the factor 42/26 = 21/13, but exits 7× more expensively: the dodge trades
  downstream room for upstream admission).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCPairBlock
import TournamentH7.LRCChainDichotomy

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- A singles `ChainOK` ledger is a `BChainOK` ledger on mapped single levels. -/
theorem bChainOK_of_chainOK : ∀ (l : List ℤ) (L : ℝ), (∀ x ∈ l, 0 < x) →
    LRC14.ChainOK l L → LRC14.BChainOK (l.map LRC14.BLevel.single) L := by
  intro l
  induction l with
  | nil => intro L _ _; trivial
  | cons w ws ih =>
    intro L hpos hchain
    obtain ⟨hentry, hrest⟩ := hchain
    have hw : 0 < w := hpos w (List.mem_cons_self ..)
    refine ⟨hw, hentry, ?_⟩
    have hout : LRC14.BLevel.out (LRC14.BLevel.single w) = 1 / (2 * (w : ℝ)) := rfl
    rw [hout]
    exact ih (1 / (2 * (w : ℝ)))
      (fun x hx => hpos x (List.mem_cons_of_mem _ hx)) hrest

/-- **The pair-at-the-dense-position split.**  Cite the `j` smallest speeds over
`B`, dodge the dense pair `(w(j), w(j+1))` as one pair level, chain the rest as
singles.  The `j = 11` instance has an empty singles tail. -/
theorem lonely_of_pair_split (cite : LRCUpTo13) (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) (j : Fin 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < (j : ℕ) → w i ≤ B)
    (hpair : w j.succ < 3 * w j.castSucc)
    (hentry : (13 : ℝ) / 7 ≤ (w j.castSucc : ℝ)
      * (2 * (((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1) * (B : ℝ)))))
    (htail : LRC14.BChainOK (((List.ofFn w).drop ((j : ℕ) + 2)).map LRC14.BLevel.single)
      (1 / (14 * (w j.succ : ℝ)))) :
    ∃ t : ℝ, Lonely 14 w t := by
  set ℓs : List LRC14.BLevel :=
    LRC14.BLevel.pair (w j.castSucc) (w j.succ)
      :: ((List.ofFn w).drop ((j : ℕ) + 2)).map LRC14.BLevel.single with hℓs
  have hj13 : (j : ℕ) < 13 := by have := j.isLt; omega
  apply LRC14.cite_blockchain_lonely cite w (fun i => (hpos i).ne')
    (j : ℕ) (by have := j.isLt; omega) B hB
    (fun i hi => by rw [abs_of_pos (hpos i)]; exact hcited i hi)
    ℓs
  · -- the split: cited, the pair, or a mapped single
    intro i
    rcases lt_trichotomy ((i : ℕ)) ((j : ℕ)) with hlt | heq | hgt
    · exact Or.inl hlt
    · right
      refine ⟨LRC14.BLevel.pair (w j.castSucc) (w j.succ), List.mem_cons_self .., ?_⟩
      have hieq : i = j.castSucc := Fin.ext heq
      rw [hieq, abs_of_pos (hpos j.castSucc)]
      show w j.castSucc ∈ [w j.castSucc, w j.succ]
      exact List.mem_cons_self ..
    · rcases Nat.eq_or_lt_of_le hgt with heq1 | hgt1
      · right
        refine ⟨LRC14.BLevel.pair (w j.castSucc) (w j.succ), List.mem_cons_self .., ?_⟩
        have hieq : i = j.succ := Fin.ext (by
          show (i : ℕ) = (j : ℕ) + 1
          omega)
        rw [hieq, abs_of_pos (hpos j.succ)]
        show w j.succ ∈ [w j.castSucc, w j.succ]
        exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
      · -- i ≥ j+2: the mapped single
        right
        have hki : (j : ℕ) + 2 ≤ (i : ℕ) := by omega
        have hlen : (i : ℕ) - ((j : ℕ) + 2)
            < ((List.ofFn w).drop ((j : ℕ) + 2)).length := by
          simp only [List.length_drop, List.length_ofFn]
          have := i.isLt
          omega
        have hval : ((List.ofFn w).drop ((j : ℕ) + 2))[(i : ℕ) - ((j : ℕ) + 2)] = w i := by
          rw [List.getElem_drop]
          simp only [List.getElem_ofFn]
          congr 1
          exact Fin.ext (by
            show (j : ℕ) + 2 + ((i : ℕ) - ((j : ℕ) + 2)) = (i : ℕ)
            omega)
        refine ⟨LRC14.BLevel.single (w i), ?_, ?_⟩
        · rw [hℓs]
          apply List.mem_cons_of_mem
          apply List.mem_map.mpr
          exact ⟨w i, by rw [← hval]; exact List.getElem_mem hlen, rfl⟩
        · rw [abs_of_pos (hpos i)]
          show w i ∈ [w i]
          exact List.mem_cons_self ..
  · -- the ledger: pair OK + pair entry + singles tail
    refine ⟨⟨hpos j.castSucc, hmono Fin.castSucc_lt_succ.le, hpair⟩, hentry, ?_⟩
    have hout : LRC14.BLevel.out (LRC14.BLevel.pair (w j.castSucc) (w j.succ))
        = 1 / (14 * (w j.succ : ℝ)) := rfl
    rw [hout]
    exact htail

/-- **The lifted dense core**: the THM-937 certificate plus the pair-split
obstruction — either the compressed triple (`w(j+2) < 21·w(j+1)`) or the pair-entry
fee failure at `j ≥ 1`. -/
def TripleDenseCore (w : Fin 13 → ℤ) : Prop :=
  ∃ j : Fin 12,
    w j.succ < 3 * w j.castSucc ∧
    (∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ) ∧
    (∀ k : Fin 12, j ≤ k →
      2 * (12 - ((k : ℕ) : ℤ)) * w k.succ < 21 * (((k : ℕ) : ℤ) + 2) * w k.castSucc) ∧
    ((∃ h2 : (j : ℕ) + 2 < 13, w ⟨(j : ℕ) + 2, h2⟩ < 21 * w j.succ) ∨
     (∃ _h1 : 1 ≤ (j : ℕ),
       (13 - ((j : ℕ) : ℤ)) * w j.castSucc
         < 13 * (((j : ℕ) : ℤ) + 1) * w ⟨(j : ℕ) - 1, by omega⟩))

/-- **THE LIFTED DICHOTOMY.**  Every sorted positive 13-family is lonely — via a
singles split (THM-937) or via the pair dodge at the dense position — or lies in the
explicit `TripleDenseCore` normal form. -/
theorem lonely_or_tripleCore (cite : LRCUpTo13) (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) :
    (∃ t : ℝ, Lonely 14 w t) ∨ TripleDenseCore w := by
  rcases lonely_or_denseCore cite w hpos hmono with hlonely | hcore
  · exact Or.inl hlonely
  obtain ⟨j, hjbad, hladder, hfee⟩ := hcore
  -- the two pair-split conditions
  by_cases htail : ∀ h2 : (j : ℕ) + 2 < 13, 21 * w j.succ ≤ w ⟨(j : ℕ) + 2, h2⟩
  · by_cases hentry : (1 ≤ (j : ℕ)) →
        13 * (((j : ℕ) : ℤ) + 1) * w ⟨(j : ℕ) - 1, by omega⟩
          ≤ (13 - ((j : ℕ) : ℤ)) * w j.castSucc
    · -- both hold: the pair split closes the family
      left
      have hj13 : (j : ℕ) < 12 := j.isLt
      -- the singles tail ledger (empty at j = 11)
      have htailOK : LRC14.BChainOK
          (((List.ofFn w).drop ((j : ℕ) + 2)).map LRC14.BLevel.single)
          (1 / (14 * (w j.succ : ℝ))) := by
        by_cases hj11 : (j : ℕ) + 2 < 13
        · apply bChainOK_of_chainOK
          · intro x hx
            obtain ⟨i, rfl⟩ := List.mem_ofFn.mp (List.mem_of_mem_drop hx)
            exact hpos i
          · apply chainOK_drop_aux w hpos ((j : ℕ) + 1)
              (fun k hk => hladder k (by
                have : (j : ℕ) < (k : ℕ) := by omega
                exact this))
              (12 - ((j : ℕ) + 2)) ((j : ℕ) + 2) hj11 (by omega) (by omega)
            -- entry: 3/2 ≤ w(j+2) · (1/(14·w(j+1)))
            have h21 := htail hj11
            have hwj1 : (0 : ℝ) < (w j.succ : ℝ) := by exact_mod_cast hpos j.succ
            have h21R : (21 : ℝ) * (w j.succ : ℝ) ≤ (w ⟨(j : ℕ) + 2, hj11⟩ : ℝ) := by
              exact_mod_cast h21
            rw [mul_one_div, le_div_iff₀ (by positivity)]
            linarith
        · -- j = 11: drop 13 = [], empty ledger
          have h13 : (j : ℕ) + 2 = 13 := by have := j.isLt; omega
          have hnil : (List.ofFn w).drop ((j : ℕ) + 2) = [] := by
            rw [h13]
            apply List.drop_eq_nil_of_le
            simp only [List.length_ofFn]
            omega
          rw [hnil]
          trivial
      by_cases hj0 : (j : ℕ) = 0
      · -- j = 0: free entry over B = 1
        apply lonely_of_pair_split cite w hpos hmono j 1 one_pos
          (fun i hi => absurd hi (by omega))
          hjbad ?_ htailOK
        have hw0 : (1 : ℝ) ≤ (w j.castSucc : ℝ) := by
          exact_mod_cast hpos j.castSucc
        rw [hj0]
        push_cast
        norm_num
        nlinarith [hw0]
      · -- j ≥ 1: entry from the fee over B = w(j−1)
        have hj1 : 1 ≤ (j : ℕ) := by omega
        set B : ℤ := w ⟨(j : ℕ) - 1, by omega⟩ with hBdef
        have hB : 0 < B := hpos _
        have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
        apply lonely_of_pair_split cite w hpos hmono j B hB
          (fun i hi => hmono (by
            show (i : ℕ) ≤ (j : ℕ) - 1
            omega))
          hjbad ?_ htailOK
        have hfeeZ := hentry hj1
        have hfeeR : (13 : ℝ) * (((j : ℕ) : ℝ) + 1) * (B : ℝ)
            ≤ (13 - ((j : ℕ) : ℝ)) * (w j.castSucc : ℝ) := by
          exact_mod_cast hfeeZ
        have hwjR : (0 : ℝ) < (w j.castSucc : ℝ) := by exact_mod_cast hpos j.castSucc
        have hden : (0 : ℝ) < 14 * (((j : ℕ) : ℝ) + 1) * (B : ℝ) := by positivity
        rw [show (w j.castSucc : ℝ)
              * (2 * (((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1) * (B : ℝ))))
            = (2 * ((13 : ℝ) - (j : ℕ)) * (w j.castSucc : ℝ))
              / (14 * (((j : ℕ) : ℝ) + 1) * (B : ℝ)) from by ring,
          le_div_iff₀ hden]
        nlinarith [hfeeR]
    · -- entry fee fails at j ≥ 1
      right
      push Not at hentry
      obtain ⟨h1, hlt⟩ := hentry
      exact ⟨j, hjbad, hladder, hfee, Or.inr ⟨h1, by linarith⟩⟩
  · -- the compressed triple
    right
    push Not at htail
    obtain ⟨h2, hlt⟩ := htail
    exact ⟨j, hjbad, hladder, hfee, Or.inl ⟨h2, by linarith⟩⟩

/-- **The triple-core obligation**: the residual clauses verbatim, plus the
`TripleDenseCore` certificate on the sorted absolute speeds. -/
def TripleCoreObligation : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    LRC14.CoveringFamily v →
    GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      TripleDenseCore (fun i => |v (σ i)|)) →
    ∃ t : ℝ, Lonely 14 v t

/-- **The residual obligation from the triple core** — strictly sharper than
THM-937's `DenseCoreObligation`: the pair dodge additionally closes every dense
core whose successor jumps by `21×` (or which sits at the top). -/
theorem residualObligation_of_tripleCore (cite : LRCUpTo13)
    (htriple : TripleCoreObligation) : ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
  set va : Fin 13 → ℤ := fun i => |v i| with hva
  have hva_pos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
  set σ : Equiv.Perm (Fin 13) := Tuple.sort va with hσ
  set w : Fin 13 → ℤ := va ∘ σ with hw
  have hw_mono : Monotone w := Tuple.monotone_sort va
  have hw_pos : ∀ i, 0 < w i := fun i => hva_pos (σ i)
  rcases lonely_or_tripleCore cite w hw_pos hw_mono with ⟨t, ht⟩ | hcore
  · refine ⟨t, (LRC14.lonely_abs_iff 14 v t).mp ?_⟩
    exact (LRC14.lonely_comp_equiv 14 va t σ).mp ht
  · exact htriple v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
      ⟨σ, hw_mono, hcore⟩

/-- **LRC(14) from the citation node and the triple-core obligation.** -/
theorem lrc14_of_tripleCore (cite : LRCUpTo13) (htriple : TripleCoreObligation) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_tripleCore cite htriple)

/-! ## Axiom audit -/
#print axioms bChainOK_of_chainOK
#print axioms lonely_of_pair_split
#print axioms lonely_or_tripleCore
#print axioms lrc14_of_tripleCore

end LRC14Grand
end LonelyRunner
