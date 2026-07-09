/-
  TournamentH7.LRCSchurRigidity — the equality case of the Schur-triple maximiser
  (kind-pasteur-2026-07-09-S114).

  opus-S183 (`LRCSchurTriples.schurTriple_card_le`) proved the upper bound `E₃(S) ≤ C(k,2)` via the
  injection `(a,b) ↦ {a, a+b}` of Schur pairs into the 2-subsets of `S`.  This file characterises
  **equality**: `E₃(S) = C(k,2)` iff `S` is a **dilated interval** `{d, 2d, …, kd}`.

  The bridge is closure under positive differences.  The injection is a *bijection* onto the 2-subsets
  exactly when every 2-subset `{x,y}` (`x<y`) is hit, i.e. `(x, y−x)` is a Schur pair, i.e. `y−x ∈ S` —
  so `E₃ = C(k,2) ⟺ ClosedUnderDiff S`.  And a finite difference-closed set of positive naturals is
  rigid: normalising by `d = min S` gives a difference-closed set containing `1`, which is forced to be
  `{1,…,k}`, so `S = {d,2d,…,kd}` (`dilated_of_closedUnderDiff`).  The converse is immediate.

  This is the discrete rigidity companion of `LRCAPTight.mreach_AP_eq` (`M(AP)=1/14`): the interval is
  the unique maximiser of Schur triples, as the AP is the unique LRC(14) equality extremal.
  Self-contained (imports Mathlib; the E₃ bound is opus's, cited).
-/
import Mathlib

open Finset

namespace LRCSchurRigidity

/-- `S` is **closed under positive differences**: `x < y` in `S` implies `y − x ∈ S`. -/
def ClosedUnderDiff (S : Finset ℕ) : Prop := ∀ x ∈ S, ∀ y ∈ S, x < y → y - x ∈ S

/-- `S` is a **dilated interval** `{d, 2d, …, kd}` (`k = |S|`, `d > 0`). -/
def DilatedInterval (S : Finset ℕ) : Prop :=
  ∃ d : ℕ, 0 < d ∧ S = (Finset.Icc 1 S.card).image (· * d)

/-- **Every element of a difference-closed set is a multiple of its minimum.** -/
theorem dvd_of_closedUnderDiff {S : Finset ℕ} (hne : S.Nonempty) (h0 : 0 ∉ S)
    (hcl : ClosedUnderDiff S) : ∀ s ∈ S, S.min' hne ∣ s := by
  set d := S.min' hne with hd
  have hd0 : 0 < d := Nat.pos_of_ne_zero (fun h => h0 (h ▸ S.min'_mem hne))
  intro s
  induction s using Nat.strong_induction_on with
  | _ s ih =>
    intro hs
    rcases eq_or_lt_of_le (S.min'_le s hs) with heq | hlt
    · exact heq ▸ dvd_refl d
    · have hsd : s - d ∈ S := hcl d (S.min'_mem hne) s hs hlt
      have hdvd : d ∣ (s - d) := ih (s - d) (by omega) hsd
      have : d ∣ (s - d) + d := Nat.dvd_add hdvd (dvd_refl d)
      rwa [Nat.sub_add_cancel (le_of_lt hlt)] at this

/-- **Rigidity: a difference-closed set of positive naturals is a dilated interval.** -/
theorem dilated_of_closedUnderDiff {S : Finset ℕ} (hne : S.Nonempty) (h0 : 0 ∉ S)
    (hcl : ClosedUnderDiff S) : DilatedInterval S := by
  set d := S.min' hne with hd
  have hd0 : 0 < d := Nat.pos_of_ne_zero (fun h => h0 (h ▸ S.min'_mem hne))
  have hdvd := dvd_of_closedUnderDiff hne h0 hcl
  have hdvd_sub : ∀ a b : ℕ, d ∣ a → d ∣ b → (a - b) / d = a / d - b / d := by
    rintro a b ⟨p, rfl⟩ ⟨q, rfl⟩
    rw [Nat.mul_div_cancel_left _ hd0, Nat.mul_div_cancel_left _ hd0,
        mul_comm d p, mul_comm d q, ← Nat.sub_mul, Nat.mul_div_cancel _ hd0]
  -- normalise: D = S / d
  set D : Finset ℕ := S.image (· / d) with hDdef
  -- the scaling `· * d` is a bijection `D ↔ S`
  have hDS : D.image (· * d) = S := by
    ext s
    simp only [hDdef, Finset.mem_image]
    constructor
    · rintro ⟨t, ⟨s', hs', rfl⟩, rfl⟩
      rwa [Nat.div_mul_cancel (hdvd s' hs')]
    · intro hs
      exact ⟨s / d, ⟨s, hs, rfl⟩, Nat.div_mul_cancel (hdvd s hs)⟩
  have hcardD : D.card = S.card := by
    rw [← hDS]
    exact (Finset.card_image_of_injOn (fun a _ b _ hab => Nat.eq_of_mul_eq_mul_right hd0 hab)).symm
  -- `D` is difference-closed and contains `1 = min`
  have h1D : 1 ∈ D := by
    simp only [hDdef, Finset.mem_image]
    exact ⟨d, S.min'_mem hne, by rw [Nat.div_self hd0]⟩
  have hstep : ∀ a ∈ D, 2 ≤ a → a - 1 ∈ D := by
    intro a ha ha2
    simp only [hDdef, Finset.mem_image] at ha ⊢
    obtain ⟨s, hs, rfl⟩ := ha
    -- s = (s/d)*d ∈ S with s/d = a ≥ 2, so s - d ∈ S and (s-d)/d = s/d - 1
    have hsdvd := hdvd s hs
    have hsval : s / d * d = s := Nat.div_mul_cancel hsdvd
    have hdlt : d < s := by
      have h2d : 2 * d ≤ s := by rw [← hsval]; exact Nat.mul_le_mul_right d ha2
      omega
    have hsd : s - d ∈ S := hcl d (S.min'_mem hne) s hs hdlt
    refine ⟨s - d, hsd, ?_⟩
    rw [hdvd_sub s d hsdvd (dvd_refl d), Nat.div_self hd0]
  -- `D` closed downward from `1` ⟹ `D = Icc 1 (max D)`
  have hDne : D.Nonempty := ⟨1, h1D⟩
  have hmem : ∀ i, 1 ≤ i → i ≤ D.max' hDne → i ∈ D := by
    have key : ∀ t, D.max' hDne - t ∈ D ∨ D.max' hDne - t < 1 := by
      intro t
      induction t with
      | zero => left; simpa using D.max'_mem hDne
      | succ n ih =>
        rcases ih with hin | hlt
        · rcases Nat.lt_or_ge (D.max' hDne - n) 2 with hlt2 | hge2
          · right; omega
          · left
            have := hstep _ hin hge2
            have he : D.max' hDne - n - 1 = D.max' hDne - (n + 1) := by omega
            rwa [he] at this
        · right; omega
    intro i hi1 hiM
    have := key (D.max' hDne - i)
    have he : D.max' hDne - (D.max' hDne - i) = i := by omega
    rw [he] at this
    rcases this with h | h
    · exact h
    · omega
  have hDeq : D = Finset.Icc 1 (D.max' hDne) := by
    apply Finset.Subset.antisymm
    · intro a ha
      rw [Finset.mem_Icc]
      refine ⟨?_, D.le_max' a ha⟩
      by_contra hlt
      push_neg at hlt
      interval_cases a
      exact h0 (by rw [← hDS]; simp only [Finset.mem_image]; exact ⟨0, ha, by simp⟩)
    · intro a ha
      rw [Finset.mem_Icc] at ha
      exact hmem a ha.1 ha.2
  -- card pins `max D = k`, so `D = Icc 1 k` and `S = (Icc 1 k).image (·*d)`
  have hmaxk : D.max' hDne = S.card := by
    have := congrArg Finset.card hDeq
    rw [Nat.card_Icc, hcardD] at this
    omega
  refine ⟨d, hd0, ?_⟩
  conv_lhs => rw [← hDS]
  rw [hDeq, hmaxk]

/-- **Converse: a dilated interval is closed under positive differences.** -/
theorem closedUnderDiff_of_dilated {S : Finset ℕ} (h : DilatedInterval S) : ClosedUnderDiff S := by
  obtain ⟨d, hd0, hSeq⟩ := h
  intro x hx y hy hxy
  rw [hSeq] at hx hy ⊢
  simp only [Finset.mem_image, Finset.mem_Icc] at hx hy ⊢
  obtain ⟨i, ⟨hi1, hik⟩, rfl⟩ := hx
  obtain ⟨j, ⟨hj1, hjk⟩, rfl⟩ := hy
  have hij : i < j := by
    by_contra hcon; push_neg at hcon
    exact absurd hxy (not_lt.mpr (Nat.mul_le_mul_right d hcon))
  exact ⟨j - i, ⟨by omega, by omega⟩, by rw [Nat.sub_mul]⟩

end LRCSchurRigidity
