/-
  TournamentH7.LRC14Dispatch — the dispatch layer (kind-pasteur-2026-07-02-S10, HYP-3967).

  Four pieces, closing the named remaining wiring of the module map:

  1. THE MREACH INSTANTIATION (module 6 → skeleton): `Mreach_ge_of_lonely` — the reverse
     of the skeleton's compactness handoff.  With `lonely_of_Mreach_ge` this makes
     `1/14 ≤ Mreach v  ⟺  ∃ t, Lonely 14 v t` (`Mreach_ge_iff_lonely`), so every certified
     family FEEDS the skeleton's reach vocabulary: the hpartA-style flow is derivable for
     the whole corpus, with no opaque intermediary.

  2. THE FUEL-CHECKER GATE (module-6 soundness instantiation): `LadderPack` + the decidable
     `ladderOK` + `ladder_speeds_lonely`/`lonely_of_ladder_mem` — ONE decidable check per
     multi-cluster pack; soundness = kps `cert_ladder` behind the norm→`Lonely` bridge.
     The recursion's mathematical content is `SepChain` (already structural); the gate is
     its packaging as a single `decide`-shaped predicate over pack data.

  3. THE CASE-SPLIT WIRING (shapes → censuses): executable row matchers `matchesRow13/12`
     and the page dispatcher `dispatch` (decidable per family), with soundness
     `lonely_of_dispatch` / `Mreach_ge_of_dispatch`: any 13-family that syntactically
     matches a census row at its admissible reference speed is lonely.

  4. THE THM-602 COMPLETENESS SURFACE, stated concretely: `DispatchComplete W` — every
     covering family, after the normalization moves (signs, permutation, integer scale),
     either dispatches to the corpus or lies in the bounded window `≤ W`.  The reduction
     `lrc14_of_dispatchComplete` : DispatchComplete W → (bounded-window check) →
     LRC14Statement is sorry-free glue consuming the S9 normalization trio.  HONESTY NOTE:
     over the CURRENT pages `DispatchComplete` is a target, not a fact — wide clusters
     (diameter > the census window) and multi-cluster shapes need the pack rows and the
     THM-602 case tree; the Prop is stated so the pages can grow under it.  This is the
     last parameter, in final vocabulary.
-/
import TournamentH7.LRC14CertRoute
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14

/-! ## 1. The Mreach instantiation -/

/-- The reverse compactness bridge: a lonely time forces the skeleton's concrete
max-min reach to clear `1/14`.  (Periodicity moves the witness into `[0,1]`.) -/
theorem Mreach_ge_of_lonely {v : Fin 13 → ℤ} (h : ∃ t : ℝ, Lonely 14 v t) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v := by
  obtain ⟨t, ht⟩ := h
  have hmr : (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v t :=
    (LRC14Concrete.lonely_iff_minReach_ge v t).mp ht
  have h1 : Int.fract t + (⌊t⌋ : ℝ) = t := by
    linarith [Int.floor_add_fract t]
  have hper : LRC14Concrete.minReach v (Int.fract t) = LRC14Concrete.minReach v t := by
    calc LRC14Concrete.minReach v (Int.fract t)
        = LRC14Concrete.minReach v (Int.fract t + (⌊t⌋ : ℝ)) :=
          (LRC14Concrete.minReach_int_periodic v (Int.fract t) ⌊t⌋).symm
      _ = LRC14Concrete.minReach v t := by rw [h1]
  have hbdd : BddAbove (LRC14Concrete.minReach v '' Set.Icc (0 : ℝ) 1) := by
    refine ⟨1 / 2, ?_⟩
    rintro y ⟨s, _, rfl⟩
    exact LRC14Concrete.minReach_le_half v s
  have hmem : LRC14Concrete.minReach v (Int.fract t)
      ∈ LRC14Concrete.minReach v '' Set.Icc (0 : ℝ) 1 :=
    ⟨Int.fract t, ⟨Int.fract_nonneg t, (Int.fract_lt_one t).le⟩, rfl⟩
  calc (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v (Int.fract t) := by rw [hper]; exact hmr
    _ ≤ LRC14Concrete.Mreach v := le_csSup hbdd hmem

/-- **The reach equivalence**: the skeleton's `Mreach ≥ 1/14` and the certificate
corpus's `∃ lonely time` are the same statement.  This is the module-6 → skeleton
instantiation: every certified family now clears the skeleton's reach node. -/
theorem Mreach_ge_iff_lonely {v : Fin 13 → ℤ} (hv : ∀ i, v i ≠ 0) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v ↔ ∃ t : ℝ, Lonely 14 v t :=
  ⟨LRC14Concrete.lonely_of_Mreach_ge v hv, Mreach_ge_of_lonely⟩

/-! ## 2. The fuel-checker gate (module-6 soundness as one decidable predicate) -/

/-- A multi-cluster certificate pack: the working band, the budget, the window,
the bounded speeds, and the ladder levels. -/
structure LadderPack where
  h      : ℚ
  μ      : ℚ
  lo     : ℚ
  hi     : ℚ
  P      : List ℤ
  levels : List WitnessCert.Level

/-- The decidable pack gate: EVERY hypothesis of `cert_ladder`, bundled.  This is
the fuel-checker: one `decide`/`native_decide` per pack certifies the whole
multi-cluster family. -/
def ladderOK (D : LadderPack) : Prop :=
  ((1 : ℚ) / 14 ≤ D.h) ∧ (0 ≤ D.h) ∧ (0 < D.μ) ∧ (D.lo < D.hi) ∧
  (∀ s ∈ D.P, 0 ≤ s) ∧ (∀ s ∈ D.P, WitnessCert.arcSafe D.h s 0 D.lo D.hi) ∧
  WitnessCert.SepChain D.μ (D.hi - D.lo) D.levels ∧
  (∀ L ∈ D.levels, ∀ o ∈ L.offs, 0 ≤ o ∧ WitnessCert.arcSafe (D.h + D.μ) o L.c D.lo D.hi)

instance (D : LadderPack) : Decidable (ladderOK D) := by
  unfold ladderOK
  infer_instance

/-- The full speed list certified by a pack: the bounded part plus every
`L.V − o` across the levels. -/
def LadderPack.speeds (D : LadderPack) : List ℤ :=
  D.P ++ D.levels.flatMap fun L => L.offs.map fun o => L.V - o

/-- **Fuel-checker soundness**: a passing pack certifies its entire speed list at
a common time, at the `1/14` band. -/
theorem ladder_speeds_lonely (D : LadderPack) (hOK : ladderOK D) :
    ∃ τ : ℝ, ∀ s ∈ D.speeds,
      (1 : ℝ) / 14 ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  obtain ⟨h14, h0, hμ, hlohi, hPpos, hPsafe, hsep, hLsafe⟩ := hOK
  obtain ⟨τ, hP, hL⟩ := WitnessCert.cert_ladder h0 hμ D.P D.levels hPpos hPsafe hlohi hsep hLsafe
  have h14R : (1 : ℝ) / 14 ≤ ((D.h : ℚ) : ℝ) := by
    have hq := (Rat.cast_le (K := ℝ)).mpr h14
    norm_num at hq
    exact hq
  refine ⟨τ, fun s hs => ?_⟩
  unfold LadderPack.speeds at hs
  rcases List.mem_append.mp hs with hsP | hsL
  · exact le_trans h14R (hP s hsP)
  · simp only [List.mem_flatMap, List.mem_map] at hsL
    obtain ⟨L, hLmem, o, homem, rfl⟩ := hsL
    exact le_trans h14R (hL L hLmem o homem)

/-- **The pack dispatcher**: any 13-family drawn from a passing pack's speed list is
lonely — and hence clears the skeleton's reach node. -/
theorem lonely_of_ladder_mem (D : LadderPack) (hOK : ladderOK D)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ∈ D.speeds) : ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨τ, hτ⟩ := ladder_speeds_lonely D hOK
  exact ⟨τ, lonely_of_norm_forall fun i => hτ _ (hv i)⟩

/-! ## 3. The case-split wiring: executable census matchers -/

/-- The reference speed a family would have if it matched a 13-offset row: read it
off the first runner. -/
def rowRef13 (v : Fin 13 → ℤ) (K : WitnessCert.CertRow) : ℤ := v 0 + K.offs.getD 0 0

/-- Decidable syntactic match against a 13-offset census row at the admissible
reference speed. -/
def matchesRow13 (v : Fin 13 → ℤ) (K : WitnessCert.CertRow) : Prop :=
  K.Vs ≤ rowRef13 v K ∧ ∀ i : Fin 13, v i = rowRef13 v K - K.offs.getD (i : ℕ) 0

instance (v : Fin 13 → ℤ) (K : WitnessCert.CertRow) : Decidable (matchesRow13 v K) := by
  unfold matchesRow13
  infer_instance

/-- Decidable syntactic match against a 1 + 12 census row. -/
def matchesRow12 (v : Fin 13 → ℤ) (K : WitnessCert.CertRow) : Prop :=
  K.Vs ≤ rowRef13 v K ∧ ∀ i : Fin 13,
    v i = if (i : ℕ) < 12 then rowRef13 v K - K.offs.getD (i : ℕ) 0 else K.P.getD 0 1

instance (v : Fin 13 → ℤ) (K : WitnessCert.CertRow) : Decidable (matchesRow12 v K) := by
  unfold matchesRow12
  infer_instance

theorem lonely_of_matchesRow13 (v : Fin 13 → ℤ) (K : WitnessCert.CertRow)
    (hK : K ∈ WitnessCert.page13) (hm : matchesRow13 v K) : ∃ t : ℝ, Lonely 14 v t := by
  have hfun : v = rowFamily13 K (rowRef13 v K) := funext fun i => hm.2 i
  rw [hfun]
  exact rowFamily13_lonely K hK _ hm.1

theorem lonely_of_matchesRow12 (v : Fin 13 → ℤ) (K : WitnessCert.CertRow)
    (hK : K ∈ WitnessCert.page12) (hm : matchesRow12 v K) : ∃ t : ℝ, Lonely 14 v t := by
  have hfun : v = rowFamily12 K (rowRef13 v K) := funext fun i => hm.2 i
  rw [hfun]
  exact rowFamily12_lonely K hK _ hm.1

/-- **The page dispatcher** (decidable per family): the family syntactically matches
some census row at its admissible reference speed. -/
def dispatch (v : Fin 13 → ℤ) : Prop :=
  (∃ K ∈ WitnessCert.page13, matchesRow13 v K) ∨
  (∃ K ∈ WitnessCert.page12, matchesRow12 v K)

instance (v : Fin 13 → ℤ) : Decidable (dispatch v) := by
  unfold dispatch
  infer_instance

/-- Dispatch soundness: matched families are lonely. -/
theorem lonely_of_dispatch (v : Fin 13 → ℤ) (hd : dispatch v) : ∃ t : ℝ, Lonely 14 v t := by
  rcases hd with ⟨K, hK, hm⟩ | ⟨K, hK, hm⟩
  · exact lonely_of_matchesRow13 v K hK hm
  · exact lonely_of_matchesRow12 v K hK hm

/-- Matched families clear the skeleton's reach node. -/
theorem Mreach_ge_of_dispatch (v : Fin 13 → ℤ) (hd : dispatch v) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v :=
  Mreach_ge_of_lonely (lonely_of_dispatch v hd)

/-! ## 4. The THM-602 completeness surface (the last parameter, concretely) -/

/-- **The completeness target** (THM-602's Lean face): every covering family, after
the normalization moves — signs (`|·|`), relabeling (`σ`), integer scale (`c`) —
either dispatches to the certified corpus or lies inside the bounded window `≤ W`.
Over the current pages this is a TARGET, not a fact: the pages must grow (wide
clusters, multi-cluster pack rows) before it can be proved.  It is stated so they
can grow under it. -/
def DispatchComplete (W : ℤ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∃ (c : ℤ) (w : Fin 13 → ℤ) (σ : Equiv.Perm (Fin 13)),
      c ≠ 0 ∧ (∀ i, v i = c * w i) ∧ dispatch ((fun i => |w i|) ∘ σ)) ∨
    (∀ i, |v i| ≤ W)

/-- **The final reduction (sorry-free glue).**  Completeness at window `W`, plus the
finite bounded-window check, gives LRC(14).  Consumes the S9 normalization trio:
signs, relabeling, scale. -/
theorem lrc14_of_dispatchComplete (W : ℤ)
    (hcomp : DispatchComplete W)
    (hwindow : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ W) →
      ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hc
  rcases hcomp v hv hc with ⟨c, w, σ, hc0, hvw, hd⟩ | hbnd
  · -- dispatched after normalization: unwind σ, then |·|, then the scale
    obtain ⟨t, ht⟩ := lonely_of_dispatch _ hd
    have h1 : Lonely 14 (fun i => |w i|) t :=
      (lonely_comp_equiv 14 (fun i => |w i|) t σ).mp ht
    have h2 : ∃ s : ℝ, Lonely 14 w s := ⟨t, (lonely_abs_iff 14 w t).mp h1⟩
    have h3 : ∃ s : ℝ, Lonely 14 (fun i => c * w i) s :=
      lonely_exists_of_scale 14 w c hc0 h2
    obtain ⟨s, hs⟩ := h3
    refine ⟨s, ?_⟩
    have hfun : v = fun i => c * w i := funext hvw
    rw [hfun]
    exact hs
  · exact hwindow v hv hbnd

end LRC14
end LonelyRunner
