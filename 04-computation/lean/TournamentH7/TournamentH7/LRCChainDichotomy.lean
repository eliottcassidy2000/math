/-
  TournamentH7.LRCChainDichotomy — THE SPLIT-SEARCH DICHOTOMY: wiring the S19 chain
  engine into the grand assembly's residual (death-star-2026-07-16-S33, HYP-7161).

  kind-pasteur's S19 `cite_chain_lonely` (LRCChainPeel) closes any family that splits as
  [cited base of k ≤ 12 speeds ≤ B] ++ [ratio-3 chain with a valid entry fee], but the
  engine was never wired into `ResidualObligation`: no theorem SEARCHED for the split.
  This file performs the search, kernel-pure, and shrinks the residual surface:

  * `chainOK_drop` — a sorted tail with consecutive ratios ≥ 3 and a valid entry is
    `ChainOK` (the recursive fee ledger holds down the whole tail).
  * `lonely_of_singles_split` — the S19 engine consumed at a sorted split, membership
    and positivity bookkeeping discharged.
  * `lonely_or_denseCore` — THE DICHOTOMY: every sorted positive 13-family is either
    lonely outright (via cite + some split k ∈ {0,…,12}) or `ChainDenseCore`: there is
    a witness index j with ratio < 3 at j, ratios ≥ 3 strictly above j, and the entry
    fee FAILING at every split above j — the explicit integer inequalities
        2·(12 − k)·w(k+1) < 21·(k + 2)·w(k)   for all k ≥ j.
    The middle band is now EXPLICIT in Lean: above the last dense pair, every
    consecutive ratio lies in [3, 21(k+2)/(2(12−k))).
  * `residualObligation_of_denseCore` / `lrc14_of_denseCore` — the residual obligation
    (hence `LRC14Statement`) from the citation node plus loneliness of the DENSE-CORE
    families only: residual families whose sorted absolute speeds satisfy
    `ChainDenseCore`.  Strictly sharper than `ResidualObligation`: every family with
    any citable ratio-3 tail is now machine-closed.

  Split-fee arithmetic (split at position m = k+1 over base bound B = w(k), δ =
  (13−m)/(14(m+1)B)): entry `w(m)·2δ ≥ 3/2` ⟺ `21(m+1)B ≤ 2(13−m)w(m)` ⟺ (in the
  Fin-12 pair index k = m−1) `21(k+2)w(k) ≤ 2(12−k)w(k+1)`.  At k = 0 with B = 1 the
  entry is free (`26·w₀ ≥ 26 > 21`), so an all-ratios-≥-3 family closes with NO cited
  runners — the pure lacunary case.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCChainPeel
import TournamentH7.LRC14GrandAssembly

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- **The dense core**: some pair index `j` has ratio `< 3`, all pairs strictly above
`j` have ratio `≥ 3`, and the chain entry fee fails at every split at or above `j` —
the explicit middle-band normal form.  (Pair index `k : Fin 12` compares `w k.castSucc`
= w(k) against `w k.succ` = w(k+1); the split at position k+1 uses base bound w(k).) -/
def ChainDenseCore (w : Fin 13 → ℤ) : Prop :=
  ∃ j : Fin 12,
    w j.succ < 3 * w j.castSucc ∧
    (∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ) ∧
    (∀ k : Fin 12, j ≤ k →
      2 * (12 - ((k : ℕ) : ℤ)) * w k.succ < 21 * (((k : ℕ) : ℤ) + 2) * w k.castSucc)

/-- `ChainOK` for the dropped sorted tail, by downward induction: a valid entry at `m`
plus consecutive ratios `≥ 3` at every pair index at or above `m₀ ≤ m` closes the
recursive fee ledger down the whole tail.  `n` is the fuel `12 − m`. -/
theorem chainOK_drop_aux (w : Fin 13 → ℤ) (hpos : ∀ i, 0 < w i) (m₀ : ℕ)
    (hratio : ∀ j : Fin 12, m₀ ≤ (j : ℕ) → 3 * w j.castSucc ≤ w j.succ) :
    ∀ (n m : ℕ) (hm : m < 13), m₀ ≤ m → m + n = 12 → ∀ (L : ℝ),
      (3 : ℝ) / 2 ≤ (w ⟨m, hm⟩ : ℝ) * L →
      LRC14.ChainOK ((List.ofFn w).drop m) L := by
  intro n
  induction n with
  | zero =>
    intro m hm _ hsum L hentry
    have hm12 : m = 12 := by omega
    subst hm12
    have hdrop : (List.ofFn w).drop 12 = [w ⟨12, hm⟩] := by
      have hlen : 12 < (List.ofFn w).length := by
        simp only [List.length_ofFn]
        omega
      rw [List.drop_eq_getElem_cons hlen]
      have h13 : (List.ofFn w).drop 13 = [] := by
        apply List.drop_eq_nil_of_le
        simp only [List.length_ofFn]
        omega
      rw [h13]
      simp
    rw [hdrop]
    exact ⟨hentry, trivial⟩
  | succ n ih =>
    intro m hm hm₀ hsum L hentry
    have hm' : m + 1 < 13 := by omega
    have hmlt : m < 12 := by omega
    have hdrop : (List.ofFn w).drop m = w ⟨m, hm⟩ :: (List.ofFn w).drop (m + 1) := by
      have hlen : m < (List.ofFn w).length := by simp; omega
      rw [List.drop_eq_getElem_cons hlen]
      simp
    rw [hdrop]
    refine ⟨hentry, ?_⟩
    have hwm : (0 : ℝ) < (w ⟨m, hm⟩ : ℝ) := by exact_mod_cast hpos ⟨m, hm⟩
    have hr : 3 * w ⟨m, hmlt⟩.castSucc ≤ w ⟨m, hmlt⟩.succ :=
      hratio ⟨m, hmlt⟩ hm₀
    have hrR : 3 * (w ⟨m, hm⟩ : ℝ) ≤ (w ⟨m + 1, hm'⟩ : ℝ) := by
      have h1 : (3 : ℝ) * ((w ⟨m, hmlt⟩.castSucc : ℤ) : ℝ) ≤ ((w ⟨m, hmlt⟩.succ : ℤ) : ℝ) := by
        exact_mod_cast hr
      have hcs : (⟨m, hmlt⟩ : Fin 12).castSucc = (⟨m, hm⟩ : Fin 13) := rfl
      have hsc : (⟨m, hmlt⟩ : Fin 12).succ = (⟨m + 1, hm'⟩ : Fin 13) := rfl
      rwa [hcs, hsc] at h1
    have hentry' : (3 : ℝ) / 2 ≤ (w ⟨m + 1, hm'⟩ : ℝ) * (1 / (2 * (w ⟨m, hm⟩ : ℝ))) := by
      rw [mul_one_div, le_div_iff₀ (by positivity)]
      linarith
    exact ih (m + 1) hm' (by omega) (by omega) (1 / (2 * (w ⟨m, hm⟩ : ℝ))) hentry'

/-- **The S19 chain engine at a sorted split.**  A positive family whose speeds at
indices `≥ k` form a `ChainOK` tail over the cited base bound `B` is lonely.  All the
membership/positivity bookkeeping of `cite_chain_lonely` is discharged here. -/
theorem lonely_of_singles_split (cite : LRCUpTo13) (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → w i ≤ B)
    (hchain : LRC14.ChainOK ((List.ofFn w).drop k)
      (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))))) :
    ∃ t : ℝ, Lonely 14 w t := by
  apply LRC14.cite_chain_lonely cite w (fun i => (hpos i).ne') k hk B hB
    (fun i hi => by rw [abs_of_pos (hpos i)]; exact hcited i hi)
    ((List.ofFn w).drop k)
    (fun x hx => by
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp (List.mem_of_mem_drop hx)
      exact hpos i)
    (fun i => by
      by_cases hik : (i : ℕ) < k
      · exact Or.inl hik
      · right
        rw [abs_of_pos (hpos i)]
        have hki : k ≤ (i : ℕ) := le_of_not_gt hik
        have hlen : (i : ℕ) - k < ((List.ofFn w).drop k).length := by
          simp only [List.length_drop, List.length_ofFn]
          omega
        have hval : ((List.ofFn w).drop k)[(i : ℕ) - k] = w i := by
          rw [List.getElem_drop]
          simp only [List.getElem_ofFn]
          congr 1
          exact Fin.ext (by show k + ((i : ℕ) - k) = (i : ℕ); omega)
        rw [← hval]
        exact List.getElem_mem hlen)
    hchain

/-- **THE SPLIT-SEARCH DICHOTOMY.**  Every sorted positive 13-family is either lonely
(the citation node plus some chain split `k ∈ {0,…,12}` — searched here) or lies in
the explicit `ChainDenseCore` middle-band normal form. -/
theorem lonely_or_denseCore (cite : LRCUpTo13) (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) :
    (∃ t : ℝ, Lonely 14 w t) ∨ ChainDenseCore w := by
  classical
  by_cases hbad : ∃ j : Fin 12, w j.succ < 3 * w j.castSucc
  · -- some ratio < 3: locate the LAST one
    obtain ⟨j₀, hj₀⟩ := hbad
    set bad : Finset (Fin 12) := Finset.univ.filter (fun j => w j.succ < 3 * w j.castSucc)
      with hbadset
    have hne : bad.Nonempty := ⟨j₀, by simp [hbadset, hj₀]⟩
    set js : Fin 12 := bad.max' hne with hjs
    have hjs_mem : js ∈ bad := bad.max'_mem hne
    have hjs_bad : w js.succ < 3 * w js.castSucc :=
      (Finset.mem_filter.mp hjs_mem).2
    have habove : ∀ k : Fin 12, js < k → 3 * w k.castSucc ≤ w k.succ := by
      intro k hk
      by_contra hlt
      push_neg at hlt
      have hmem : k ∈ bad := by simp [hbadset, hlt]
      exact absurd (Finset.le_max' bad k hmem) (not_le.mpr hk)
    -- search the entry fee at splits ≥ js
    by_cases hgood : ∃ k : Fin 12, js ≤ k ∧
        21 * (((k : ℕ) : ℤ) + 2) * w k.castSucc ≤ 2 * (12 - ((k : ℕ) : ℤ)) * w k.succ
    · -- a split works: chain at m = k+1 over B = w k.castSucc
      left
      obtain ⟨k, hjk, hentryZ⟩ := hgood
      have hk12 : (k : ℕ) < 12 := k.isLt
      set m : ℕ := (k : ℕ) + 1 with hmdef
      have hm13 : m < 13 := by omega
      set B : ℤ := w k.castSucc with hBdef
      have hB : 0 < B := hpos _
      have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
      have hjkv : (js : ℕ) ≤ (k : ℕ) := hjk
      apply lonely_of_singles_split cite w hpos m (by omega) B hB
      · -- cited bound: sorted below the split
        intro i hi
        have hik : i ≤ k.castSucc := by
          have hv : (i : ℕ) ≤ (k : ℕ) := by omega
          exact hv
        exact hmono hik
      · -- the ledger: entry from the fee inequality, ratios from maximality of js
        refine chainOK_drop_aux w hpos m
          (fun j hj => habove j (by
            have hjv : (js : ℕ) < (j : ℕ) := by omega
            exact hjv))
          (12 - m) m hm13 (le_refl m) (by omega) _ ?_
        -- entry: 3/2 ≤ w(m) · 2δ with δ = (13−m)/(14(m+1)B)
        have hwm : (⟨m, hm13⟩ : Fin 13) = k.succ := rfl
        rw [hwm]
        have hwkR : (0 : ℝ) < (w k.succ : ℝ) := by exact_mod_cast hpos k.succ
        have hZR : 21 * (((k : ℕ) : ℤ) + 2) * B ≤ 2 * (12 - ((k : ℕ) : ℤ)) * w k.succ :=
          hentryZ
        have hZRr : (21 : ℝ) * (((k : ℕ) : ℝ) + 2) * (B : ℝ)
            ≤ 2 * (12 - ((k : ℕ) : ℝ)) * (w k.succ : ℝ) := by
          exact_mod_cast hZR
        have h13m : ((13 : ℝ) - m) = 12 - ((k : ℕ) : ℝ) := by
          rw [hmdef]; push_cast; ring
        have hm1 : ((m : ℝ) + 1) = ((k : ℕ) : ℝ) + 2 := by
          rw [hmdef]; push_cast; ring
        have hden : (0 : ℝ) < 14 * ((m : ℝ) + 1) * (B : ℝ) := by
          rw [hm1]; positivity
        rw [show (w k.succ : ℝ) * (2 * (((13 : ℝ) - m) / (14 * ((m : ℝ) + 1) * (B : ℝ))))
              = (2 * ((13 : ℝ) - m) * (w k.succ : ℝ)) / (14 * ((m : ℝ) + 1) * (B : ℝ))
            from by ring,
          le_div_iff₀ hden]
        rw [h13m, hm1]
        nlinarith [hZRr]
    · -- no split works above js: the dense core
      right
      push_neg at hgood
      exact ⟨js, hjs_bad, habove, fun k hk => hgood k hk⟩
  · -- every ratio ≥ 3: the pure lacunary split at k = 0 over B = 1
    left
    push_neg at hbad
    have hZ : (0 : ℕ) < 13 := by omega
    apply lonely_of_singles_split cite w hpos 0 (by omega) 1 one_pos
      (fun i hi => absurd hi (by omega))
    refine chainOK_drop_aux w hpos 0 (fun j _ => hbad j) 12 0 hZ (le_refl 0) rfl _ ?_
    -- entry at m = 0, B = 1: w₀ ≥ 1 gives w₀·(2·13/14) = w₀·13/7 ≥ 13/7 > 3/2
    have hw0 : (1 : ℝ) ≤ (w ⟨0, hZ⟩ : ℝ) := by exact_mod_cast hpos ⟨0, hZ⟩
    push_cast
    nlinarith [hw0]

/-- **The dense-core obligation**: the residual clauses verbatim, plus the explicit
`ChainDenseCore` certificate on the sorted absolute speeds. -/
def DenseCoreObligation : Prop :=
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
      ChainDenseCore (fun i => |v (σ i)|)) →
    ∃ t : ℝ, Lonely 14 v t

/-- **The residual obligation from the dense core**: the split search closes every
residual family whose sorted speeds admit a citable chain; only the dense core
remains.  Strictly sharper than `ResidualObligation`. -/
theorem residualObligation_of_denseCore (cite : LRCUpTo13)
    (hdense : DenseCoreObligation) : ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
  set va : Fin 13 → ℤ := fun i => |v i| with hva
  have hva_pos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
  set σ : Equiv.Perm (Fin 13) := Tuple.sort va with hσ
  set w : Fin 13 → ℤ := va ∘ σ with hw
  have hw_mono : Monotone w := Tuple.monotone_sort va
  have hw_pos : ∀ i, 0 < w i := fun i => hva_pos (σ i)
  rcases lonely_or_denseCore cite w hw_pos hw_mono with ⟨t, ht⟩ | hcore
  · -- the chain split closed it: transport back through sort and signs
    refine ⟨t, (lonely_abs_iff 14 v t).mp ?_⟩
    exact (lonely_comp_equiv 14 va t σ).mp ht
  · -- the dense core: hand the certificate to the obligation
    exact hdense v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
      ⟨σ, hw_mono, hcore⟩

/-- **LRC(14) from the citation node and the dense-core obligation** — the sharpest
composed surface on this route: every residual family with any citable ratio-3 tail
is machine-closed; what remains is the explicit middle-band dense core. -/
theorem lrc14_of_denseCore (cite : LRCUpTo13) (hdense : DenseCoreObligation) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_denseCore cite hdense)

/-! ## Axiom audit -/
#print axioms lonely_or_denseCore
#print axioms residualObligation_of_denseCore
#print axioms lrc14_of_denseCore

end LRC14Grand
end LonelyRunner
