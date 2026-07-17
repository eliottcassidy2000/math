/-
  TournamentH7.LRCBlockSplitLift — THE GLEVEL LIFT of the chain-split dichotomy
  (death-star-2026-07-16-S35, HYP-7172; lifts THM-938 through kps-S22's fat-block
  window step).

  THM-938's residual (`TripleDenseCore`) is a dense pair whose successor stays within
  `21×` — three clustered runners no singles/pair chain passes.  kps-S22
  (`LRCFatBlockChain.block_window_step`) certifies a common good SUB-WINDOW of width
  `ε` for ANY block of runners whose fat-tooth mass fits:

      Σ_{w ∈ ws} [ L/7 + 3/(7w) + ε·(w·L + 3) ]  <  L − ε.

  This file composes that step with the S19 citation transport and singles tail:

  * `lonely_of_block_split` — THE GENERIC BLOCK SPLIT: cite `k ≤ 12` runners ≤ B at
    gap `1/(k+1)`; ONE block (arbitrary `ws`, mass fee at the transported window
    `L = 2δ`) opens a good `ε`-window; the S19 singles ledger (`ChainOK`) runs INSIDE
    that window; every runner is cited, in the block, or in the tail.  The engines
    compose because `block_window_step` certifies the WHOLE `[x, x+ε]` window and
    `good_chain` finds its point inside it.
  * `lonely_or_quadCore` — THE LIFTED DICHOTOMY: on a `TripleDenseCore` family, put
    the dense TRIPLE `{w(j), w(j+1), w(j+2)}` in the block with tail width
    `ε = 3/(2·w(j+3))` (top-block `ε` free when `j ≥ 10`): either the triple's mass
    fee passes — lonely — or the family lies in `QuadDenseCore`: TripleDenseCore
    PLUS the explicit triple-fee failure.  The 7-wall is now the only wall left
    standing between the fee ledger and the dense core.
  * `residualObligation_of_quadCore` / `lrc14_of_quadCore` — the wire, strictly
    sharper than THM-938's `TripleCoreObligation`.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCFatBlockChain
import TournamentH7.LRCPairLiftDichotomy

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- **The generic block split.**  Cite `k ≤ 12` runners bounded by `B`; one fat block
`ws` whose mass fee fits the transported window `2δ` at width `ε`; a singles `ChainOK`
tail inside the `ε`-window; every runner cited, in the block, or in the tail. -/
theorem lonely_of_block_split (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ws : List ℤ) (hwpos : ∀ w ∈ ws, 0 < w)
    (tail : List ℤ) (htpos : ∀ w ∈ tail, 0 < w)
    (ε : ℝ) (hε : 0 ≤ ε)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ |v i| ∈ ws ∨ |v i| ∈ tail)
    (hmass : ((ws.map fun (w : ℤ) =>
        (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) / 7
          + 3 / (7 * (w : ℝ))
          + ε * ((w : ℝ) * (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) + 3))).sum
      < 2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))) - ε)
    (htail : LRC14.ChainOK tail ε) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hk13 : k ≤ 13 := by omega
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  have hkR : (k : ℝ) ≤ 12 := by exact_mod_cast hk
  set δ : ℝ := ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by
    rw [hδ]
    apply div_pos (by linarith) (by positivity)
  -- cite the first k
  set wc : Fin k → ℤ := fun j => v (Fin.castLE hk13 j) with hwc
  have hwcne : ∀ j, wc j ≠ 0 := fun j => hv _
  obtain ⟨t₀, hcite⟩ := cite k hk wc hwcne
  -- the block opens an ε-window inside [t₀ − δ, t₀ + δ]
  obtain ⟨x, hx1, hx2, hblock⟩ := LRC14.block_window_step ws hwpos ε
    (t₀ - δ) (t₀ + δ) hε (by linarith)
    (by
      have h2δ : (t₀ + δ) - (t₀ - δ) = 2 * δ := by ring
      rw [h2δ]
      exact hmass)
  -- the singles tail digs its point inside the ε-window
  obtain ⟨τ, hτ1, hτ2, hτtail⟩ := LRC14.good_chain tail x ε hε htpos htail
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | hmem | hmem
  · -- cited: Lipschitz transport of the 1/(k+1) margin
    have hidx : Fin.castLE hk13 ⟨(i : ℕ), hlt⟩ = i := by
      apply Fin.ext
      rfl
    have h0 : (1 : ℝ) / (k + 1 : ℕ) ≤ |(v i : ℝ) * t₀ - m| := by
      simpa [hwc, hidx] using hcite ⟨(i : ℕ), hlt⟩ m
    have hvB : |(v i : ℝ)| ≤ (B : ℝ) := by
      rw [← Int.cast_abs]
      exact_mod_cast hcited i hlt
    have htri : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by
      calc |(v i : ℝ) * t₀ - m|
          = |((v i : ℝ) * τ - m) + (v i : ℝ) * (t₀ - τ)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ) * (t₀ - τ)| := abs_add_le _ _
        _ = |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by rw [abs_mul]
    have hwin : |t₀ - τ| ≤ δ := by
      rw [abs_le]
      constructor <;> linarith
    have hlip : |(v i : ℝ)| * |t₀ - τ| ≤ (B : ℝ) * δ := by
      apply mul_le_mul hvB hwin (abs_nonneg _) hBR.le
    have hBδ : (B : ℝ) * δ = ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1)) := by
      rw [hδ]
      field_simp
    have hmargin : (1 : ℝ) / (k + 1 : ℕ) - ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1))
        = 1 / 14 := by
      have hcast : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
      rw [hcast]
      field_simp
      ring
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ = 1 / (k + 1 : ℕ) - ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1)) := hmargin.symm
      _ ≤ |(v i : ℝ) * τ - m| := by
          rw [hBδ] at hlip
          linarith
  · -- in the block: the whole ε-window is good, and τ sits inside it
    have hgood := hblock τ hτ1 hτ2 |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood m
      rw [hcast] at this
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ this
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood (-m)
      rw [hcast] at this
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at this
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ this
  · -- in the tail: the chain's point serves it
    have hgood := hτtail |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood m
      rw [hcast] at this
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ this
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood (-m)
      rw [hcast] at this
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at this
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ this


/-- **The lifted dense core**: the THM-938 certificate plus the triple-block
obstruction — the explicit mass-fee failure at tail width `ε = 3/(2·w(j+3))`
(interior positions `j ≤ 9`), or the deferred top-corner `j ≥ 10`. -/
def QuadDenseCore (w : Fin 13 → ℤ) : Prop :=
  ∃ j : Fin 12,
    w j.succ < 3 * w j.castSucc ∧
    (∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ) ∧
    (∀ k : Fin 12, j ≤ k →
      2 * (12 - ((k : ℕ) : ℤ)) * w k.succ < 21 * (((k : ℕ) : ℤ) + 2) * w k.castSucc) ∧
    (((∃ h2 : (j : ℕ) + 2 < 13, w ⟨(j : ℕ) + 2, h2⟩ < 21 * w j.succ) ∨
      (∃ _h1 : 1 ≤ (j : ℕ),
        (13 - ((j : ℕ) : ℤ)) * w j.castSucc
          < 13 * (((j : ℕ) : ℤ) + 1) * w ⟨(j : ℕ) - 1, by omega⟩))) ∧
    (10 ≤ (j : ℕ) ∨
      ∃ h3 : (j : ℕ) + 3 < 13,
        ¬ ((([w j.castSucc, w j.succ, w ⟨(j : ℕ) + 2, by omega⟩].map fun (u : ℤ) =>
            (2 * (((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1)
              * ((if (j : ℕ) = 0 then 1 else w ⟨(j : ℕ) - 1, by omega⟩) : ℤ) : ℝ))) / 7
              + 3 / (7 * (u : ℝ))
              + (3 / (2 * (w ⟨(j : ℕ) + 3, h3⟩ : ℝ)))
                * ((u : ℝ) * (2 * (((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1)
                  * ((if (j : ℕ) = 0 then 1 else w ⟨(j : ℕ) - 1, by omega⟩) : ℤ) : ℝ))) + 3))).sum
          < 2 * (((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1)
              * ((if (j : ℕ) = 0 then 1 else w ⟨(j : ℕ) - 1, by omega⟩) : ℤ) : ℝ))
            - 3 / (2 * (w ⟨(j : ℕ) + 3, h3⟩ : ℝ))))

/-- **THE GLEVEL-LIFTED DICHOTOMY.**  Every sorted positive 13-family is lonely —
via the THM-938 routes or via the dense-TRIPLE block split — or lies in the
explicit `QuadDenseCore` normal form. -/
theorem lonely_or_quadCore (cite : LRCUpTo13) (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) :
    (∃ t : ℝ, Lonely 14 w t) ∨ QuadDenseCore w := by
  rcases lonely_or_tripleCore cite w hpos hmono with hlonely | hcore
  · exact Or.inl hlonely
  obtain ⟨j, hjbad, hladder, hfee, hobstruct⟩ := hcore
  by_cases hj10 : 10 ≤ (j : ℕ)
  · exact Or.inr ⟨j, hjbad, hladder, hfee, hobstruct, Or.inl hj10⟩
  have h3 : (j : ℕ) + 3 < 13 := by omega
  have h2 : (j : ℕ) + 2 < 13 := by omega
  set B : ℤ := (if (j : ℕ) = 0 then 1 else w ⟨(j : ℕ) - 1, by omega⟩) with hBdef
  have hB : 0 < B := by
    rw [hBdef]
    split
    · norm_num
    · exact hpos _
  set δ : ℝ := ((13 : ℝ) - (j : ℕ)) / (14 * (((j : ℕ) : ℝ) + 1) * (B : ℝ)) with hδdef
  set ε : ℝ := 3 / (2 * (w ⟨(j : ℕ) + 3, h3⟩ : ℝ)) with hεdef
  have hεpos : 0 < ε := by
    rw [hεdef]
    have : (0 : ℝ) < (w ⟨(j : ℕ) + 3, h3⟩ : ℝ) := by exact_mod_cast hpos _
    positivity
  by_cases hmass : (([w j.castSucc, w j.succ, w ⟨(j : ℕ) + 2, by omega⟩].map fun (u : ℤ) =>
      (2 * δ) / 7 + 3 / (7 * (u : ℝ)) + ε * ((u : ℝ) * (2 * δ) + 3))).sum
      < 2 * δ - ε
  · -- the triple block passes: lonely
    left
    apply lonely_of_block_split cite w (fun i => (hpos i).ne') (j : ℕ)
      (by have := j.isLt; omega) B hB
      (fun i hi => by
        rw [abs_of_pos (hpos i)]
        rw [hBdef]
        split
        · omega
        · exact hmono (by show (i : ℕ) ≤ (j : ℕ) - 1; omega))
      [w j.castSucc, w j.succ, w ⟨(j : ℕ) + 2, by omega⟩]
      (by
        intro u hu
        rcases List.mem_cons.mp hu with rfl | hu
        · exact hpos _
        rcases List.mem_cons.mp hu with rfl | hu
        · exact hpos _
        rcases List.mem_singleton.mp hu with rfl
        · exact hpos _)
      ((List.ofFn w).drop ((j : ℕ) + 3))
      (fun x hx => by
        obtain ⟨i, rfl⟩ := List.mem_ofFn.mp (List.mem_of_mem_drop hx)
        exact hpos i)
      ε hεpos.le
      ?_ hmass ?_
    · -- the split: cited, block, or tail
      intro i
      rcases lt_trichotomy ((i : ℕ)) ((j : ℕ)) with hlt | heq | hgt
      · exact Or.inl hlt
      · refine Or.inr (Or.inl ?_)
        have hieq : i = j.castSucc := Fin.ext heq
        rw [hieq, abs_of_pos (hpos j.castSucc)]
        exact List.mem_cons_self ..
      · rcases Nat.eq_or_lt_of_le hgt with heq1 | hgt1
        · refine Or.inr (Or.inl ?_)
          have hieq : i = j.succ := Fin.ext (by show (i : ℕ) = (j : ℕ) + 1; omega)
          rw [hieq, abs_of_pos (hpos j.succ)]
          exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
        · rcases Nat.eq_or_lt_of_le hgt1 with heq2 | hgt2
          · refine Or.inr (Or.inl ?_)
            have hieq : i = ⟨(j : ℕ) + 2, h2⟩ := Fin.ext (by
              show (i : ℕ) = (j : ℕ) + 2
              omega)
            rw [hieq, abs_of_pos (hpos ⟨(j : ℕ) + 2, h2⟩)]
            exact List.mem_cons_of_mem _
              (List.mem_cons_of_mem _ (List.mem_singleton.mpr rfl))
          · -- i ≥ j+3: the tail
            refine Or.inr (Or.inr ?_)
            rw [abs_of_pos (hpos i)]
            have hki : (j : ℕ) + 3 ≤ (i : ℕ) := by omega
            have hlen : (i : ℕ) - ((j : ℕ) + 3)
                < ((List.ofFn w).drop ((j : ℕ) + 3)).length := by
              simp only [List.length_drop, List.length_ofFn]
              have := i.isLt
              omega
            have hval : ((List.ofFn w).drop ((j : ℕ) + 3))[(i : ℕ) - ((j : ℕ) + 3)]
                = w i := by
              rw [List.getElem_drop]
              simp only [List.getElem_ofFn]
              congr 1
              exact Fin.ext (by
                show (j : ℕ) + 3 + ((i : ℕ) - ((j : ℕ) + 3)) = (i : ℕ)
                omega)
            rw [← hval]
            exact List.getElem_mem hlen
    · -- the tail ledger at ε = 3/(2·w(j+3)): entry is exact, ratios from the ladder
      refine chainOK_drop_aux w hpos ((j : ℕ) + 1)
        (fun k hk => hladder k (by
          have : (j : ℕ) < (k : ℕ) := by omega
          exact this))
        (12 - ((j : ℕ) + 3)) ((j : ℕ) + 3) h3 (by omega) (by omega) _ ?_
      have hwR : (0 : ℝ) < (w ⟨(j : ℕ) + 3, h3⟩ : ℝ) := by exact_mod_cast hpos _
      rw [hεdef]
      rw [show (w ⟨(j : ℕ) + 3, h3⟩ : ℝ) * (3 / (2 * (w ⟨(j : ℕ) + 3, h3⟩ : ℝ)))
          = 3 / 2 from by field_simp]
  · -- the fee fails: the quad core
    right
    exact ⟨j, hjbad, hladder, hfee, hobstruct, Or.inr ⟨h3, hmass⟩⟩

/-- **The quad-core obligation**: the residual clauses verbatim, plus the
`QuadDenseCore` certificate on the sorted absolute speeds. -/
def QuadCoreObligation : Prop :=
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
      QuadDenseCore (fun i => |v (σ i)|)) →
    ∃ t : ℝ, Lonely 14 v t

/-- **The residual obligation from the quad core** — strictly sharper than THM-938's
`TripleCoreObligation`. -/
theorem residualObligation_of_quadCore (cite : LRCUpTo13)
    (hquad : QuadCoreObligation) : ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
  set va : Fin 13 → ℤ := fun i => |v i| with hva
  have hva_pos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
  set σ : Equiv.Perm (Fin 13) := Tuple.sort va with hσ
  set w : Fin 13 → ℤ := va ∘ σ with hw
  have hw_mono : Monotone w := Tuple.monotone_sort va
  have hw_pos : ∀ i, 0 < w i := fun i => hva_pos (σ i)
  rcases lonely_or_quadCore cite w hw_pos hw_mono with ⟨t, ht⟩ | hcore
  · refine ⟨t, (LRC14.lonely_abs_iff 14 v t).mp ?_⟩
    exact (LRC14.lonely_comp_equiv 14 va t σ).mp ht
  · exact hquad v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
      ⟨σ, hw_mono, hcore⟩

/-- **LRC(14) from the citation node and the quad-core obligation.** -/
theorem lrc14_of_quadCore (cite : LRCUpTo13) (hquad : QuadCoreObligation) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_quadCore cite hquad)

/-! ## Axiom audit -/
#print axioms lonely_of_block_split
#print axioms lonely_or_quadCore
#print axioms lrc14_of_quadCore

end LRC14Grand
end LonelyRunner
