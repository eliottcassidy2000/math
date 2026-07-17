/- Thm866Order.lean — klein-2026-07-16-S317.
   THM-866 rung two, part II: the order theorem and the exchange step.
   * `arc_iff_lt_of_injective` — distinct scores force the arc relation to BE the
     score order (proved by strong induction on the upper-set cardinality; the
     score bijection comes from pigeonhole, no subtournament machinery needed).
   * `transitive_of_injective` — distinct scores ⟹ transitive (THM-866 step 4).
   * `exists_tie_of_not_transitive` — the contrapositive: off the transitive
     order a tied pair exists.
   * `exists_plus_eight_flip` — THM-866 steps (3)+(4) combined: every
     non-transitive tournament admits a tie-splitting flip raising xLevel by
     exactly 8 (the F3 exchange step).  No sorries. -/
import TournamentH7.Thm866Flip

open Finset

namespace Tournament

variable {n : ℕ}

/-- Scores live strictly below `n`. -/
theorem outDegree_lt_n (T : Tournament n) (v : Fin n) : T.outDegree v < n := by
  have hsub : Finset.univ.filter (fun w => T.arc v w = true) ⊆ Finset.univ.erase v := by
    intro j hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
    refine Finset.mem_erase.mpr ⟨?_, Finset.mem_univ j⟩
    rintro rfl
    rw [T.irrefl j] at hj
    exact absurd hj (by simp)
  have h1 : T.outDegree v ≤ (Finset.univ.erase v).card := Finset.card_le_card hsub
  rw [Finset.card_erase_of_mem (Finset.mem_univ v), Finset.card_univ,
      Fintype.card_fin] at h1
  have := v.pos
  omega

/-- The strict upper set of a vertex in score order. -/
def upperSet (T : Tournament n) (u : Fin n) : Finset (Fin n) :=
  Finset.univ.filter (fun x => T.outDegree u < T.outDegree x)

theorem mem_upperSet {T : Tournament n} {u x : Fin n} :
    x ∈ T.upperSet u ↔ T.outDegree u < T.outDegree x := by
  simp [upperSet]

/-- With injective scores, the upper set has exactly `n − 1 − s(u)` members
    (pigeonhole: the score map is a bijection of `Fin n`). -/
theorem card_upperSet (T : Tournament n) (hinj : Function.Injective T.outDegree)
    (u : Fin n) : (T.upperSet u).card = n - 1 - T.outDegree u := by
  set sFin : Fin n → Fin n := fun v => ⟨T.outDegree v, T.outDegree_lt_n v⟩ with hsFin
  have hinj' : Function.Injective sFin := by
    intro a b hab
    exact hinj (by simpa [hsFin, Fin.ext_iff] using hab)
  have hbij : Function.Bijective sFin := (Finite.injective_iff_bijective).mp hinj'
  have hcard : (T.upperSet u).card = (Finset.Ioi (sFin u)).card := by
    apply Finset.card_bij (fun x _ => sFin x)
    · intro x hx
      rw [mem_upperSet] at hx
      exact Finset.mem_Ioi.mpr (by simpa [hsFin, Fin.lt_def] using hx)
    · intro a _ b _ hab
      exact hinj' hab
    · intro y hy
      obtain ⟨x, hx⟩ := hbij.2 y
      refine ⟨x, ?_, hx⟩
      rw [mem_upperSet]
      have := Finset.mem_Ioi.mp hy
      rw [← hx] at this
      simpa [hsFin, Fin.lt_def] using this
  rw [hcard, Fin.card_Ioi]

/-- **The order theorem**: injective scores force `arc = score order`. -/
theorem arc_iff_lt_of_injective (T : Tournament n)
    (hinj : Function.Injective T.outDegree) :
    ∀ u v : Fin n, v ≠ u → (T.arc u v = true ↔ T.outDegree v < T.outDegree u) := by
  -- the one-step engine: if every strictly-higher vertex already obeys the
  -- order theorem, then u does too
  have step : ∀ u : Fin n,
      (∀ x, x ∈ T.upperSet u → ∀ v, v ≠ x →
        (T.arc x v = true ↔ T.outDegree v < T.outDegree x)) →
      ∀ v, v ≠ u → (T.arc u v = true ↔ T.outDegree v < T.outDegree u) := by
    intro u IHloc v hvu
    have hbeats : ∀ x, x ∈ T.upperSet u → T.arc x u = true := by
      intro x hx
      have hux : u ≠ x := by
        rintro rfl
        exact absurd (mem_upperSet.mp hx) (lt_irrefl _)
      exact (IHloc x hx u hux).mpr (mem_upperSet.mp hx)
    -- the out-set of u avoids u and its upper set
    have hsub : Finset.univ.filter (fun j => T.arc u j = true)
        ⊆ (Finset.univ.erase u) \ (T.upperSet u) := by
      intro j hj
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
      refine Finset.mem_sdiff.mpr ⟨Finset.mem_erase.mpr ⟨?_, Finset.mem_univ j⟩, ?_⟩
      · rintro rfl
        rw [T.irrefl j] at hj
        exact absurd hj (by simp)
      · intro hjU
        exact T.asym u j ⟨hj, hbeats j hjU⟩
    have hUsub : T.upperSet u ⊆ Finset.univ.erase u := by
      intro x hx
      refine Finset.mem_erase.mpr ⟨?_, Finset.mem_univ x⟩
      rintro rfl
      exact absurd (mem_upperSet.mp hx) (lt_irrefl _)
    have hslt : T.outDegree u < n := T.outDegree_lt_n u
    have htarget : ((Finset.univ.erase u) \ (T.upperSet u)).card = T.outDegree u := by
      rw [Finset.card_sdiff hUsub, Finset.card_erase_of_mem (Finset.mem_univ u),
          Finset.card_univ, Fintype.card_fin, card_upperSet T hinj u]
      omega
    have heq : Finset.univ.filter (fun j => T.arc u j = true)
        = (Finset.univ.erase u) \ (T.upperSet u) := by
      apply Finset.eq_of_subset_of_card_le hsub
      rw [htarget]
      exact le_of_eq rfl
    constructor
    · intro h
      have hv : v ∈ (Finset.univ.erase u) \ (T.upperSet u) := by
        rw [← heq]
        simp [h]
      have hnotU : ¬ T.outDegree u < T.outDegree v :=
        fun hc => (Finset.mem_sdiff.mp hv).2 (mem_upperSet.mpr hc)
      have hne : T.outDegree v ≠ T.outDegree u := fun hc => hvu (hinj hc)
      omega
    · intro h
      have hv : v ∈ (Finset.univ.erase u) \ (T.upperSet u) := by
        refine Finset.mem_sdiff.mpr ⟨Finset.mem_erase.mpr ⟨hvu, Finset.mem_univ v⟩, ?_⟩
        intro hc
        exact absurd (mem_upperSet.mp hc) (by omega)
      rw [← heq] at hv
      simpa using hv
  -- strong induction on the upper-set cardinality
  have main : ∀ m : ℕ, ∀ u : Fin n, (T.upperSet u).card ≤ m →
      ∀ v, v ≠ u → (T.arc u v = true ↔ T.outDegree v < T.outDegree u) := by
    intro m
    induction m with
    | zero =>
        intro u hm
        apply step u
        intro x hx
        have : (T.upperSet u).card = 0 := Nat.le_zero.mp hm
        rw [Finset.card_eq_zero] at this
        rw [this] at hx
        exact absurd hx (Finset.notMem_empty x)
    | succ m IH =>
        intro u hm
        apply step u
        intro x hx
        apply IH x
        have hss : T.upperSet x ⊂ T.upperSet u := by
          constructor
          · intro y hy
            rw [mem_upperSet] at hy ⊢
            exact lt_trans (mem_upperSet.mp hx) hy
          · intro hcon
            exact absurd (mem_upperSet.mp (hcon hx)) (lt_irrefl _)
        have := Finset.card_lt_card hss
        omega
  intro u v hvu
  exact main (T.upperSet u).card u (le_refl _) v hvu

/-- **THM-866 step (4)**: distinct scores ⟹ transitive. -/
theorem transitive_of_injective (T : Tournament n)
    (hinj : Function.Injective T.outDegree) :
    ∀ u v x : Fin n, T.arc u v = true → T.arc v x = true → T.arc u x = true := by
  intro u v x huv hvx
  have hvu : v ≠ u := by
    rintro rfl
    rw [T.irrefl v] at huv
    exact absurd huv (by simp)
  have hxv : x ≠ v := by
    rintro rfl
    rw [T.irrefl x] at hvx
    exact absurd hvx (by simp)
  have hxu : x ≠ u := by
    rintro rfl
    exact T.asym u v ⟨huv, hvx⟩
  have h1 := (T.arc_iff_lt_of_injective hinj u v hvu).mp huv
  have h2 := (T.arc_iff_lt_of_injective hinj v x hxv).mp hvx
  exact (T.arc_iff_lt_of_injective hinj u x hxu).mpr (lt_trans h2 h1)

/-- The contrapositive: a non-transitive tournament has a tied pair. -/
theorem exists_tie_of_not_transitive (T : Tournament n)
    (hnt : ¬ (∀ u v x : Fin n, T.arc u v = true → T.arc v x = true →
      T.arc u x = true)) :
    ∃ u v : Fin n, u ≠ v ∧ T.outDegree u = T.outDegree v := by
  by_contra hno
  push_neg at hno
  apply hnt
  apply transitive_of_injective
  intro a b hab
  by_contra hne
  exact (hno a b hne) hab

/-- **THM-866, steps (3)+(4) combined — the F3 exchange step**: every
    non-transitive tournament admits a tie-splitting arc flip that raises the
    axis level by exactly 8. -/
theorem exists_plus_eight_flip (T : Tournament n)
    (hnt : ¬ (∀ u v x : Fin n, T.arc u v = true → T.arc v x = true →
      T.arc u x = true)) :
    ∃ (w l : Fin n) (hne : w ≠ l), T.arc w l = true ∧
      (T.flip w l hne).xLevel = T.xLevel + 8 := by
  obtain ⟨u, v, hne, htie⟩ := T.exists_tie_of_not_transitive hnt
  have hdev : T.dev u = T.dev v := (T.dev_eq_iff_outDegree_eq).mpr htie
  rcases T.total u v hne with h | h
  · exact ⟨u, v, hne, h, T.xLevel_flip_tie hne h hdev⟩
  · exact ⟨v, u, hne.symm, h, T.xLevel_flip_tie hne.symm h hdev.symm⟩

end Tournament
