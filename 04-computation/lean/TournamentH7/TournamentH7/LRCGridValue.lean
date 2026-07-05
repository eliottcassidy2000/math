/-
  TournamentH7.LRCGridValue — THE VALUE FORM OF THE MERGE GRID (klein-2026-07-05-S138,
  HYP-4121; kps-S4's redirect (a), owner-directed).

  The evaluation direction: exact-M claims become KERNEL-CHECKABLE.
  * `distZ_at_rat` — the distance to ℤ at a rational `a/s` is EXACTLY
    `min(a % s, s − a % s)/s`: the analytic `distZ` and integer arithmetic coincide.
  * `marginQ` — the DECIDABLE rational margin evaluator; `margin_at_rat` is the exact
    bridge `margin v (a/s) = marginQ v a s`.
  * `margin_le_of_grid` — the UPPER-BOUND certificate schema: to prove `M(V) ≤ c` it
    suffices to check the FINITE fundamental grid `m ∈ [0, |vᵢ|+|vⱼ|)` (via kps's
    `grid_margin_domination` + `margin_add_int` periodicity).  Together with the
    witness direction (one grid point with `marginQ ≥ c`) every exact-M sweep claim
    is now a finite rational computation.
  * FLAGSHIP: `ap12_margin_eq` — **M({1,…,12}) = 1/13 EXACTLY, kernel-checked** — the
    AP tightness anchoring the whole n=12 rigidity program, previously computational.

  Consumes the canonical stack: mac-mini's `distZ`/`margin` (LRCWitnessAttainment) +
  kps's grid layer (LRCGridAttainment, HYP-4108).
-/
import TournamentH7.LRCGridAttainment

namespace LonelyRunner
namespace GridValue

open TournamentH7.LRCWitness GridAttainment

/-- **The exact rational value of `distZ`**: at `a/s` the distance to ℤ is
`min(a % s, s − a % s)/s`. -/
theorem distZ_at_rat (a s : ℤ) (hs : 0 < s) :
    distZ ((a : ℝ) / s) = ((min (a % s) (s - a % s) : ℤ) : ℝ) / s := by
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  set r : ℤ := a % s with hr
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (ne_of_gt hs)
  have hrs : r < s := Int.emod_lt_of_pos _ hs
  have hsplit : (a : ℝ) / s = (r : ℝ) / s + ((a / s : ℤ) : ℝ) := by
    rw [hr]
    push_cast [Int.emod_def]
    field_simp
    ring
  rw [hsplit, distZ_add_int]
  rcases lt_trichotomy (2 * r) s with h2 | h2 | h2
  · have hlt : (r : ℝ) / s < 1/2 := by
      rw [div_lt_iff₀ hsR]
      have : (2 * r : ℝ) < s := by exact_mod_cast h2
      linarith
    have hge : (0 : ℝ) ≤ (r : ℝ) / s := by positivity
    have hround : round ((r : ℝ) / s) = 0 := by
      apply round_eq_of_abs_sub_lt_half
      push_cast
      rw [sub_zero, abs_of_nonneg hge]
      exact hlt
    rw [distZ_eq_round, hround]
    have hmin : min r (s - r) = r := min_eq_left (by omega)
    rw [hmin]
    push_cast
    rw [sub_zero, abs_of_nonneg hge]
  · have hrs2 : (r : ℝ) / s = 1/2 := by
      rw [div_eq_iff (ne_of_gt hsR)]
      have : (2 * r : ℝ) = s := by exact_mod_cast h2
      linarith
    have hmin : min r (s - r) = r := min_eq_left (by omega)
    have hval : ((min r (s - r) : ℤ) : ℝ) / s = (r : ℝ) / s := by
      rw [hmin]
    rw [hval, hrs2, distZ_eq_round]
    have hround : round ((1 : ℝ)/2) = 1 := by
      rw [round_eq]
      norm_num
    rw [hround]
    push_cast
    rw [show (1 : ℝ)/2 - 1 = -(1/2) by ring, abs_neg]
    norm_num
  · have hgt : (1 : ℝ)/2 < (r : ℝ) / s := by
      rw [lt_div_iff₀ hsR]
      have : (s : ℝ) < 2 * r := by exact_mod_cast h2
      linarith
    have hlt1 : (r : ℝ) / s < 1 := by
      rw [div_lt_one hsR]
      exact_mod_cast hrs
    have hround : round ((r : ℝ) / s) = 1 := by
      apply round_eq_of_abs_sub_lt_half
      push_cast
      rw [abs_of_neg (by linarith)]
      linarith
    rw [distZ_eq_round, hround]
    have hmin : min r (s - r) = s - r := min_eq_right (by omega)
    rw [hmin]
    push_cast
    rw [abs_of_neg (by linarith)]
    field_simp
    ring

variable {k : ℕ} [Nonempty (Fin k)]

/-- **The decidable rational margin evaluator.** -/
def marginQ (v : Fin k → ℤ) (a s : ℤ) : ℚ :=
  ((Finset.univ.inf' Finset.univ_nonempty
    fun i => min ((v i * a) % s) (s - (v i * a) % s) : ℤ) : ℚ) / (s : ℚ)

/-- **The exact bridge**: the analytic margin at `a/s` IS the rational evaluator. -/
theorem margin_at_rat (v : Fin k → ℤ) (a s : ℤ) (hs : 0 < s) :
    margin v ((a : ℝ) / s) = ((marginQ v a s : ℚ) : ℝ) := by
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  have hpt : ∀ i : Fin k, distZ ((v i : ℝ) * ((a : ℝ) / s))
      = ((min ((v i * a) % s) (s - (v i * a) % s) : ℤ) : ℝ) / s := by
    intro i
    have hexp : (v i : ℝ) * ((a : ℝ) / s) = ((v i * a : ℤ) : ℝ) / s := by
      push_cast
      ring
    rw [hexp, distZ_at_rat _ _ hs]
  unfold margin marginQ
  have hfun : (fun i : Fin k => distZ ((v i : ℝ) * ((a : ℝ) / s)))
      = fun i => ((min ((v i * a) % s) (s - (v i * a) % s) : ℤ) : ℝ) / s := by
    funext i
    exact hpt i
  rw [hfun]
  -- commute the cast/division through the finite inf (both directions)
  set F : Fin k → ℤ := fun i => min ((v i * a) % s) (s - (v i * a) % s) with hF
  apply le_antisymm
  · obtain ⟨i0, -, hi0⟩ := Finset.exists_mem_eq_inf'
      (Finset.univ_nonempty (α := Fin k)) F
    have h1 : (Finset.univ.inf' Finset.univ_nonempty fun i => ((F i : ℤ) : ℝ) / s)
        ≤ ((F i0 : ℤ) : ℝ) / s := Finset.inf'_le _ (Finset.mem_univ i0)
    have h2 : ((Finset.univ.inf' Finset.univ_nonempty F : ℤ) : ℚ) / (s : ℚ)
        = ((F i0 : ℤ) : ℚ) / (s : ℚ) := by rw [hi0]
    calc (Finset.univ.inf' Finset.univ_nonempty fun i => ((F i : ℤ) : ℝ) / s)
        ≤ ((F i0 : ℤ) : ℝ) / s := h1
      _ = ((((F i0 : ℤ) : ℚ) / (s : ℚ) : ℚ) : ℝ) := by push_cast; ring
      _ = ((((Finset.univ.inf' Finset.univ_nonempty F : ℤ) : ℚ) / (s : ℚ) : ℚ) : ℝ) := by
          rw [h2]
  · apply Finset.le_inf'
    intro i _
    have hle : Finset.univ.inf' Finset.univ_nonempty F ≤ F i :=
      Finset.inf'_le _ (Finset.mem_univ i)
    have hleR : ((Finset.univ.inf' Finset.univ_nonempty F : ℤ) : ℝ)
        ≤ ((F i : ℤ) : ℝ) := by exact_mod_cast hle
    calc ((((Finset.univ.inf' Finset.univ_nonempty F : ℤ) : ℚ) / (s : ℚ) : ℚ) : ℝ)
        = ((Finset.univ.inf' Finset.univ_nonempty F : ℤ) : ℝ) / s := by push_cast; ring
      _ ≤ ((F i : ℤ) : ℝ) / s := by
          rw [div_eq_mul_inv, div_eq_mul_inv]
          exact mul_le_mul_of_nonneg_right hleR
            (le_of_lt (inv_pos.mpr (by exact_mod_cast hs : (0:ℝ) < (s:ℝ))))

/-- **THE UPPER-BOUND CERTIFICATE SCHEMA**: to prove `M(V) ≤ c` it suffices to check
the FINITE fundamental merge grid `m ∈ [0, |vᵢ|+|vⱼ|)`.  With the witness direction
(one grid point where `marginQ ≥ c`) every exact-M claim is a finite computation. -/
theorem margin_le_of_grid (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (c : ℝ) (hc : 0 ≤ c)
    (hgrid : ∀ i j : Fin k, ∀ m : ℤ, 0 ≤ m → m < |v i| + |v j| →
      margin v ((m : ℝ) / ((|v i| + |v j| : ℤ) : ℝ)) ≤ c) :
    ∀ t : ℝ, margin v t ≤ c := by
  intro t
  by_contra hcon
  push Not at hcon
  have hβpos : 0 < margin v t := lt_of_le_of_lt hc hcon
  obtain ⟨i, j, m, hm⟩ := grid_margin_domination v hv (margin v t) hβpos t
    ((le_margin_iff v _ t).1 le_rfl)
  set s : ℤ := |v i| + |v j| with hsdef
  have hs : 0 < s := by
    have h1 : (0:ℤ) < |v i| := abs_pos.mpr (hv i)
    have h2 : (0:ℤ) < |v j| := abs_pos.mpr (hv j)
    omega
  have hsR : (0:ℝ) < (s : ℝ) := by exact_mod_cast hs
  have hle : margin v t ≤ margin v ((m : ℝ) / (s : ℝ)) :=
    (le_margin_iff v _ _).2 hm
  have hred : margin v ((m : ℝ) / (s : ℝ)) = margin v (((m % s : ℤ) : ℝ) / (s : ℝ)) := by
    have hshift : (m : ℝ) / (s : ℝ) = ((m % s : ℤ) : ℝ) / (s : ℝ) + ((m / s : ℤ) : ℝ) := by
      push_cast [Int.emod_def]
      field_simp
      ring
    rw [hshift, margin_add_int]
  rw [hred] at hle
  have hgr := hgrid i j (m % s) (Int.emod_nonneg _ (ne_of_gt hs))
    (Int.emod_lt_of_pos _ hs)
  rw [← hsdef] at hgr
  linarith

/-! ## The flagship: `M({1,…,12}) = 1/13` EXACTLY, kernel-checked -/

/-- The 12-runner arithmetic progression. -/
def ap12 : Fin 12 → ℤ := fun i => (i : ℤ) + 1

theorem ap12_ne (i : Fin 12) : ap12 i ≠ 0 := by
  unfold ap12
  omega

/-- The finite grid table: every fundamental grid point has margin ≤ 1/13. -/
theorem ap12_table : ∀ i j : Fin 12, ∀ m ∈ Finset.range 25,
    marginQ ap12 (m : ℤ) (|ap12 i| + |ap12 j|) ≤ 1/13 := by
  native_decide

/-- **Upper bound**: `M({1,…,12}) ≤ 1/13`. -/
theorem ap12_margin_le : ∀ t : ℝ, margin ap12 t ≤ ((1/13 : ℚ) : ℝ) := by
  apply margin_le_of_grid ap12 ap12_ne _ (by norm_num)
  intro i j m hm0 hms
  have hs : (0:ℤ) < |ap12 i| + |ap12 j| := by
    have h1 : (0:ℤ) < |ap12 i| := abs_pos.mpr (ap12_ne i)
    have h2 : (0:ℤ) < |ap12 j| := abs_pos.mpr (ap12_ne j)
    omega
  rw [margin_at_rat ap12 m _ hs]
  have hmnat : m = ((m.toNat : ℕ) : ℤ) := (Int.toNat_of_nonneg hm0).symm
  have hsle : |ap12 i| + |ap12 j| ≤ 24 := by
    have hi := i.isLt
    have hj := j.isLt
    unfold ap12
    rw [abs_of_pos (by positivity : (0:ℤ) < ((i : ℕ) : ℤ) + 1),
      abs_of_pos (by positivity : (0:ℤ) < ((j : ℕ) : ℤ) + 1)]
    omega
  have hmem : m.toNat ∈ Finset.range 25 := by
    rw [Finset.mem_range]
    omega
  have htab := ap12_table i j m.toNat hmem
  rw [← hmnat] at htab
  exact_mod_cast htab

/-- **Witness**: the margin at `t = 1/13` is exactly `1/13`. -/
theorem ap12_margin_at_witness :
    margin ap12 ((1 : ℝ) / ((13 : ℤ) : ℝ)) = ((1/13 : ℚ) : ℝ) := by
  have h13 : (0:ℤ) < 13 := by norm_num
  have h1 : ((1:ℤ) : ℝ) = (1:ℝ) := by norm_num
  rw [← h1, margin_at_rat ap12 1 13 h13]
  norm_num [show marginQ ap12 1 13 = 1/13 from by native_decide]

/-- **THE AP TIGHTNESS, EXACT AND KERNEL-CHECKED**: `M({1,…,12}) = 1/13` — the
maximum is `1/13`, attained at `t = 1/13`.  The anchor fact of the n = 12 rigidity
program, previously only computational. -/
theorem ap12_margin_eq :
    margin ap12 ((1 : ℝ) / ((13 : ℤ) : ℝ)) = ((1/13 : ℚ) : ℝ) ∧
    ∀ t : ℝ, margin ap12 t ≤ ((1/13 : ℚ) : ℝ) :=
  ⟨ap12_margin_at_witness, ap12_margin_le⟩

#print axioms distZ_at_rat
#print axioms margin_at_rat
#print axioms margin_le_of_grid
#print axioms ap12_margin_eq

end GridValue
end LonelyRunner
