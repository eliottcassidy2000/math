/-
  TournamentH7.LRC14CertRoute — THE CONCRETE ENDGAME (kind-pasteur-2026-07-02-S9).

  The assembly surface (`LRC14Assembly.lrc14_endgame`) states LRC(14) from two node
  hypotheses phrased over the skeleton's OPAQUE measure vocabulary (`witnessG2`, `shapeOf`).
  Those parameters cannot be discharged by bookkeeping: `opaque` functions have no
  provable properties, and concretizing them means Lebesgue measure on the critical
  path — against the playbook's T1 (all-ℚ, no MeasureTheory).

  This file states the SAME frontier in fully CONCRETE vocabulary and proves the
  reduction sorry-free:

    `LRC14Statement  ⟺  every covering family is lonely`

  where `CoveringFamily v` is the concrete arithmetic predicate (for every
  `q ∈ {2,…,14}` some speed is divisible by `q`).  The forward glue is the
  denominator sieve (`sieve_one_div`, sorry-free): a non-covering family is lonely
  at `t = 1/q` for a missed `q`.  Nothing here is conditional, abstract, or opaque —
  the remaining content of LRC(14) is EXACTLY the covering branch, and it is now a
  single named `Prop` in the language of integers and `Lonely`.

  SECOND HALF: the first concrete covering discharges.  The witness-certificate
  corpus (LRCWitnessCert / LRCCertTable / ladder packs) produces loneliness in the
  `UnitAddCircle`-norm form; `lonely_of_norm_forall` bridges it to `Lonely`.
  Applied to `certAP_tail` this machine-checks: **every consecutive block
  `{V−12, …, V}` (V ≥ 15) is 14-lonely** — in particular the COVERING blocks
  (`V % 14 ≠ 13`), the first infinite family of covering LRC(14) instances proved
  lonely end-to-end in Lean.
-/
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCWitnessCert

namespace LonelyRunner
namespace LRC14

/-! ## The concrete covering predicate and the frontier equivalence -/

/-- A 13-speed family is **covering** if every modulus `q ∈ {2,…,14}` divides some
speed.  By the sieve, non-covering families are lonely; covering families are the
entire remaining content of LRC(14). -/
def CoveringFamily (v : Fin 13 → ℤ) : Prop :=
  ∀ q : ℕ, 2 ≤ q → q ≤ 14 → ∃ i, (q : ℤ) ∣ v i

/-- **The concrete endgame reduction (sorry-free).**  If every covering family is
lonely, LRC(14) holds: a non-covering family misses some `q ∈ {2,…,14}` and is
lonely at `t = 1/q` by the denominator sieve. -/
theorem lrc14_of_covering_lonely
    (hcov : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  intro v hv
  by_cases hc : CoveringFamily v
  · exact hcov v hv hc
  · unfold CoveringFamily at hc
    push_neg at hc
    obtain ⟨q, hq2, hq14, hdiv⟩ := hc
    exact ⟨(1 : ℝ) / q, sieve_one_div 14 q v hq14 (by omega) hdiv⟩

/-- **The frontier, as an equivalence.**  LRC(14) is EXACTLY the statement that
every covering family is lonely — in concrete integer/`Lonely` vocabulary, with
no opaque functions.  This is the honest remaining surface of the formalization. -/
theorem lrc14_statement_iff_covering_lonely :
    LRC14Statement ↔
      (∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
        ∃ t : ℝ, Lonely 14 v t) :=
  ⟨fun h v hv _ => h v hv, lrc14_of_covering_lonely⟩

/-! ## The norm → `Lonely` bridge

The certificate corpus proves loneliness in the `UnitAddCircle` norm form
`1/14 ≤ ‖v·τ‖`.  `Lonely` asks `1/14 ≤ |v·τ − m|` for every integer `m`.  The
bridge: the round-point is the nearest integer, and any other integer is at
distance `≥ 1/2 ≥ 1/14`. -/

theorem lonely_of_norm_forall {v : Fin 13 → ℤ} {t : ℝ}
    (h : ∀ i, (1 : ℝ) / 14 ≤ ‖(((v i : ℝ) * t : ℝ) : UnitAddCircle)‖) :
    Lonely 14 v t := by
  intro i m
  have hn := h i
  rw [UnitAddCircle.norm_eq] at hn
  set x : ℝ := (v i : ℝ) * t with hx
  rcases eq_or_ne m (round x) with rfl | hne
  · calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ ≤ |x - round x| := hn
  · have hr : |x - (round x : ℝ)| ≤ 1 / 2 := abs_sub_round x
    have hmr : (1 : ℝ) ≤ |((m : ℝ)) - (round x : ℝ)| := by
      have hne' : m - round x ≠ 0 := sub_ne_zero.mpr hne
      have h1 : (1 : ℤ) ≤ |m - round x| := Int.one_le_abs hne'
      have h1' : ((1 : ℤ) : ℝ) ≤ ((|m - round x| : ℤ) : ℝ) := by exact_mod_cast h1
      rw [Int.cast_abs] at h1'
      push_cast at h1'
      linarith
    have htri : |((m : ℝ)) - (round x : ℝ)| ≤ |(m : ℝ) - x| + |x - (round x : ℝ)| :=
      abs_sub_le _ _ _
    have habs : |x - (m : ℝ)| = |(m : ℝ) - x| := abs_sub_comm _ _
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ ≤ |x - (m : ℝ)| := by
          rw [habs]
          linarith

/-! ## The first covering discharges: consecutive blocks -/

/-- The consecutive block `{V, V−1, …, V−12}` as a 13-speed family. -/
def blockFamily (V : ℤ) : Fin 13 → ℤ := fun i => V - ((i : ℕ) : ℤ)

theorem blockFamily_nonzero (V : ℤ) (hV : 15 ≤ V) : ∀ i, blockFamily V i ≠ 0 := by
  intro i
  unfold blockFamily
  have : ((i : ℕ) : ℤ) ≤ 12 := by
    have := i.isLt
    omega
  omega

/-- Thirteen consecutive integers cover every modulus `q ≤ 13` automatically, and
cover `q = 14` exactly when `V % 14 ≠ 13`. -/
theorem blockFamily_covering (V : ℤ) (hV : 15 ≤ V) (hmod : V % 14 ≠ 13) :
    CoveringFamily (blockFamily V) := by
  intro q hq2 hq14
  have hqz : (0 : ℤ) < (q : ℤ) := by exact_mod_cast (by omega : 0 < q)
  have hr0 : 0 ≤ V % (q : ℤ) := Int.emod_nonneg V (ne_of_gt hqz)
  have hrq : V % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos V hqz
  -- the residue is an admissible offset: ≤ 12 always for q ≤ 13, and for q = 14 by hmod
  have hr12 : V % (q : ℤ) ≤ 12 := by
    rcases eq_or_ne q 14 with rfl | hq13
    · have : V % (14 : ℤ) ≠ 13 := hmod
      omega
    · have : (q : ℤ) ≤ 13 := by
        have : q ≤ 13 := by omega
        exact_mod_cast this
      omega
  refine ⟨⟨(V % (q : ℤ)).toNat, by omega⟩, ?_⟩
  unfold blockFamily
  have htn : (((V % (q : ℤ)).toNat : ℕ) : ℤ) = V % (q : ℤ) := Int.toNat_of_nonneg hr0
  rw [htn]
  exact ⟨V / (q : ℤ), by linarith [Int.emod_add_ediv V (q : ℤ)]⟩

/-- **Every consecutive block `{V−12,…,V}` (V ≥ 15) is 14-lonely** — bridged from
the witness-arc certificate `certAP_tail` (the V-independent tail row). -/
theorem blockFamily_lonely (V : ℤ) (hV : 15 ≤ V) :
    ∃ t : ℝ, Lonely 14 (blockFamily V) t := by
  obtain ⟨τ, hτ⟩ := WitnessCert.certAP_tail V hV
  refine ⟨τ, lonely_of_norm_forall (fun i => ?_)⟩
  have hi : ((i : ℕ) : ℤ) ∈ WitnessCert.offsAP := by
    have := i.isLt
    unfold WitnessCert.offsAP
    interval_cases h : (i : ℕ) <;> simp
  have hbound := hτ ((i : ℕ) : ℤ) hi
  have hval : ((WitnessCert.h14 : ℚ) : ℝ) = 1 / 14 := by
    norm_num [WitnessCert.h14]
  rw [hval] at hbound
  exact hbound

/-- **The first infinite family of COVERING LRC(14) instances machine-checked
lonely end-to-end**: consecutive blocks with `V % 14 ≠ 13` are covering AND lonely.
(The `V % 14 = 13` blocks are also lonely — they are simply handled by the sieve
instead, being non-covering.) -/
theorem covering_block_lonely (V : ℤ) (hV : 15 ≤ V) (hmod : V % 14 ≠ 13) :
    CoveringFamily (blockFamily V) ∧ ∃ t : ℝ, Lonely 14 (blockFamily V) t :=
  ⟨blockFamily_covering V hV hmod, blockFamily_lonely V hV⟩

end LRC14
end LonelyRunner
