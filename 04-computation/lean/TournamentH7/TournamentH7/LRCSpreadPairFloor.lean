/-
  TournamentH7.LRCSpreadPairFloor — THE SPREAD PAIR-FLOOR, Stage 1 (klein-2026-07-02-S118,
  HYP-4023).  The first Hunter pair credit, lower-bounded on the PRISTINE window for a
  spread consecutive pair — the quantitative input kps-S24's clustered/spread repair
  routes to the pair lane.

  THE MECHANISM (Python-pinned, lrc14_spread_pair_floor_pin_klein_S118.py: identity
  PASS 200/200 random cases, adversarial floor ≥ 0.92·(L/49) at D·L ≥ 1):
   * TWO-TOOTH TRAPEZOID IDENTITY (this stage): the overlap of `tooth w₁ n` and
     `tooth w₂ m` is EXACTLY `max 0 (min (2h/w₂) ((S − |r|)/(w₁w₂)))` with
     `r = n·w₂ − m·w₁` the lattice residue and `S = h(w₁+w₂)` the support half-width.
     The proof is case-free: min-sub-max turns the clip length into a four-way min.
   * THE RESIDUE WALK (stage 3): `r` advances by `−D (mod w₂)`, `D = w₂ − w₁`, per
     w₂-tooth; `D·L ≥ 2` forces wraps, each landing in the half-support and staying
     `~S/(2D)` consecutive teeth, each worth `≥ h/w₂` — the floor `≥ L/392`.
  The clustered complement (`D·L < 2`) is kps-S24's C0-combo lane; the wide complement
  (`14D > w₁`) is mac-mini's dense-events regime.
-/
import TournamentH7.LRCRealRegions

namespace LonelyRunner
namespace SpreadPairFloor

open LRC14 RealRegion

/-! ## Stage 1: the two-tooth trapezoid identity -/

/-- Four-way min: `min a b − max c d` distributes into the four differences. -/
theorem min_sub_max (a b c d : ℝ) :
    min a b - max c d = min (min (a - c) (a - d)) (min (b - c) (b - d)) := by
  simp only [min_def, max_def]
  split_ifs <;> linarith

/-- Balanced min: `min (S + r) (S − r) = S − |r|`. -/
theorem min_add_sub_abs (S r : ℝ) : min (S + r) (S - r) = S - |r| := by
  rcases le_total r 0 with hr | hr
  · rw [abs_of_nonpos hr, min_eq_left (by linarith)]; ring
  · rw [abs_of_nonneg hr, min_eq_right (by linarith)]

/-- The trapezoid: the overlap profile of two teeth at lattice residue `r`.
Plateau `2·(1/14)/w₂` (the narrower tooth), support `|r| < (1/14)(w₁+w₂)`. -/
noncomputable def trap (w₁ w₂ r : ℝ) : ℝ :=
  max 0 (min (2 * (1/14) / w₂) (((1/14) * (w₁ + w₂) - |r|) / (w₁ * w₂)))

theorem trap_nonneg (w₁ w₂ r : ℝ) : 0 ≤ trap w₁ w₂ r := le_max_left _ _

/-- **The two-tooth trapezoid identity**: the clip length of `tooth w₁ n` against the
window given by `tooth w₂ m` is exactly the trapezoid at the residue `n·w₂ − m·w₁`. -/
theorem clipLen_tooth_tooth (w₁ w₂ : ℤ) (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) (n m : ℤ) :
    clipLen (tooth w₁ n) ((tooth w₂ m).1) ((tooth w₂ m).2)
      = trap (w₁ : ℝ) (w₂ : ℝ) (((n * w₂ - m * w₁ : ℤ) : ℝ)) := by
  have hw₁R : (0 : ℝ) < (w₁ : ℝ) := by exact_mod_cast hw₁
  have hw₂R : (0 : ℝ) < (w₂ : ℝ) := lt_of_lt_of_le hw₁R (by exact_mod_cast hw₁₂)
  have hw₁₂R : (w₁ : ℝ) ≤ (w₂ : ℝ) := by exact_mod_cast hw₁₂
  have hmul : (0 : ℝ) < (w₁ : ℝ) * w₂ := mul_pos hw₁R hw₂R
  unfold clipLen tooth trap
  simp only
  rw [min_sub_max]
  push_cast
  congr 1
  have e₁ : ((n : ℝ) + 1/14) / w₁ - ((n : ℝ) - 1/14) / w₁ = 2 * (1/14) / (w₁ : ℝ) := by
    field_simp; ring
  have e₂ : ((m : ℝ) + 1/14) / w₂ - ((m : ℝ) - 1/14) / w₂ = 2 * (1/14) / (w₂ : ℝ) := by
    field_simp; ring
  have e₃ : ((n : ℝ) + 1/14) / w₁ - ((m : ℝ) - 1/14) / w₂
      = ((1/14) * ((w₁ : ℝ) + w₂) + ((n : ℝ) * w₂ - (m : ℝ) * w₁)) / ((w₁ : ℝ) * w₂) := by
    field_simp; ring
  have e₄ : ((m : ℝ) + 1/14) / w₂ - ((n : ℝ) - 1/14) / w₁
      = ((1/14) * ((w₁ : ℝ) + w₂) - ((n : ℝ) * w₂ - (m : ℝ) * w₁)) / ((w₁ : ℝ) * w₂) := by
    field_simp; ring
  rw [e₁, e₂, e₃, e₄]
  set P : ℝ := (w₁ : ℝ) * w₂ with hP
  set S : ℝ := (1/14) * ((w₁ : ℝ) + w₂) with hS
  set r : ℝ := (n : ℝ) * w₂ - (m : ℝ) * w₁ with hr
  -- RHS's (S − |r|)/P is the min of the two leg values
  have habs : (S - |r|) / P = min ((S + r) / P) ((S - r) / P) := by
    rw [← min_add_sub_abs S r]
    rcases le_total (S + r) (S - r) with hxy | hxy
    · rw [min_eq_left hxy, min_eq_left (by apply div_le_div_of_nonneg_right hxy hmul)]
    · rw [min_eq_right hxy, min_eq_right (by apply div_le_div_of_nonneg_right hxy hmul)]
  have hplateau : 2 * (1/14) / (w₂ : ℝ) ≤ 2 * (1/14) / (w₁ : ℝ) :=
    div_le_div_of_nonneg_left (by norm_num) hw₁R hw₁₂R
  rw [habs]
  -- min (min (2h/w₁) ((S+r)/P)) (min ((S−r)/P) (2h/w₂)) = min (2h/w₂) (min ((S+r)/P) ((S−r)/P))
  apply le_antisymm
  · apply le_min
    · exact (min_le_right _ _).trans (min_le_right _ _)
    · apply le_min
      · exact (min_le_left _ _).trans (min_le_right _ _)
      · exact (min_le_right _ _).trans (min_le_left _ _)
  · apply le_min
    · apply le_min
      · exact (min_le_left _ _).trans hplateau
      · exact (min_le_right _ _).trans (min_le_left _ _)
    · apply le_min
      · exact (min_le_right _ _).trans (min_le_right _ _)
      · exact min_le_left _ _

/-! ## Stage 2: membership + the per-tooth lower bound -/

/-- Membership: `tooth w m ∈ teeth w a b` for `m` in the enumeration range. -/
theorem tooth_mem_teeth (w : ℤ) (a b : ℝ) (m : ℤ)
    (hlo : ⌈(w : ℝ) * a⌉ - 1 ≤ m) (hhi : m ≤ ⌊(w : ℝ) * b⌋ + 1) :
    tooth w m ∈ teeth w a b := by
  unfold teeth
  rw [List.mem_map]
  refine ⟨(m - (⌈(w : ℝ) * a⌉ - 1)).toNat, ?_, ?_⟩
  · rw [List.mem_range]
    omega
  · congr 1
    omega

/-- The window-clip of an interior interval is itself. -/
theorem rclip_window_interior (a b : ℝ) (q : ℝ × ℝ) (h₁ : a ≤ q.1) (h₂ : q.2 ≤ b) :
    rclip (a, b) q = q := by
  unfold rclip
  simp only
  rw [max_eq_right h₁, min_eq_right h₂]

/-- **Stage-2 plumbing**: the first-credit term of one interior `w₂`-tooth dominates
the trapezoid at any partner `n` inside the `w₁` enumeration range. -/
theorem per_tooth_ge_trap (w₁ w₂ : ℤ) (hw₁ : 0 < w₁) (h₁₂ : w₁ ≤ w₂)
    (a b : ℝ) (m n : ℤ)
    (hint₁ : a ≤ (tooth w₂ m).1) (hint₂ : (tooth w₂ m).2 ≤ b)
    (hnlo : ⌈(w₁ : ℝ) * a⌉ - 1 ≤ n) (hnhi : n ≤ ⌊(w₁ : ℝ) * b⌋ + 1) :
    trap (w₁ : ℝ) (w₂ : ℝ) (((n * w₂ - m * w₁ : ℤ) : ℝ))
      ≤ rlength (rinter (rinter [(a, b)] [tooth w₂ m]) (teeth w₁ a b)) := by
  have hstep₁ : rinter [(a, b)] [tooth w₂ m] = [tooth w₂ m] := by
    unfold rinter
    simp only [List.flatMap_cons, List.map_cons, List.map_nil, List.flatMap_nil,
      List.append_nil]
    rw [rclip_window_interior a b _ hint₁ hint₂]
  rw [hstep₁]
  have hstep₂ : rinter [tooth w₂ m] (teeth w₁ a b)
      = (teeth w₁ a b).map (fun q => rclip (tooth w₂ m) q) := by
    unfold rinter
    simp
  rw [hstep₂]
  unfold rlength
  rw [List.map_map]
  have hmem : tooth w₁ n ∈ teeth w₁ a b := tooth_mem_teeth w₁ a b n hnlo hnhi
  have hmem' : ((fun p : ℝ × ℝ => max 0 (p.2 - p.1)) ∘ fun q => rclip (tooth w₂ m) q)
        (tooth w₁ n)
      ∈ (teeth w₁ a b).map
        ((fun p : ℝ × ℝ => max 0 (p.2 - p.1)) ∘ fun q => rclip (tooth w₂ m) q) :=
    List.mem_map_of_mem hmem
  have hval : ((fun p : ℝ × ℝ => max 0 (p.2 - p.1)) ∘ fun q => rclip (tooth w₂ m) q)
        (tooth w₁ n)
      = trap (w₁ : ℝ) (w₂ : ℝ) (((n * w₂ - m * w₁ : ℤ) : ℝ)) := by
    rw [← clipLen_tooth_tooth w₁ w₂ hw₁ h₁₂ n m]
    unfold clipLen rclip
    simp only [Function.comp]
    rw [min_comm, max_comm]
  rw [← hval]
  apply List.single_le_sum _ _ hmem'
  intro x hx
  rw [List.mem_map] at hx
  obtain ⟨q, -, rfl⟩ := hx
  exact le_max_left _ _

/-! ## Stage 3: the integer residue walk -/

/-- The walk formula: the remainder of `(m₀+j)·w₁` mod `w₂` is the start remainder
shifted down by `j·D`, mod `w₂`.  Pure `emod` algebra, no induction. -/
theorem walk_formula (w₁ w₂ D : ℤ) (hD : w₂ = w₁ + D) (m₀ j : ℤ) :
    ((m₀ + j) * w₁) % w₂ = ((m₀ * w₁) % w₂ - j * D) % w₂ := by
  have h1 : (m₀ + j) * w₁ = m₀ * w₁ - j * D + j * w₂ := by rw [hD]; ring
  rw [h1, Int.add_mul_emod_self]
  conv_rhs => rw [Int.sub_emod, Int.emod_emod_of_dvd _ dvd_rfl, ← Int.sub_emod]

/-- **One wrap of the walk**: from any start `m₀` the descending remainder enters
`[0, D)` within `w₂/D + 1` teeth, and the following `k`-run (while `k·D ≤ T`) gives
consecutive teeth whose lattice residue against SOME partner `n` is `≤ T` in size. -/
theorem walk_one_wrap (w₁ w₂ D T : ℤ) (hw₁ : 0 < w₁) (hD : w₂ = w₁ + D)
    (hDpos : 0 < D) (hDT : D ≤ T) (hTw : 2 * T ≤ w₁) (m₀ : ℤ) :
    ∃ mstar : ℤ, m₀ ≤ mstar ∧ mstar ≤ m₀ + w₂ / D ∧
      ∀ k : ℤ, 0 ≤ k → k * D ≤ T →
        ∃ n : ℤ, |n * w₂ - (mstar + k) * w₁| ≤ T := by
  have hw₂ : 0 < w₂ := by omega
  set v₀ : ℤ := (m₀ * w₁) % w₂ with hv₀
  have hv₀nn : 0 ≤ v₀ := Int.emod_nonneg _ (by omega)
  have hv₀lt : v₀ < w₂ := Int.emod_lt_of_pos _ hw₂
  set j₀ : ℤ := v₀ / D with hj₀
  have hj₀nn : 0 ≤ j₀ := Int.ediv_nonneg hv₀nn (le_of_lt hDpos)
  have hj₀le : j₀ ≤ w₂ / D := Int.ediv_le_ediv hDpos (le_of_lt hv₀lt)
  refine ⟨m₀ + j₀, by omega, by omega, ?_⟩
  -- the remainder at mstar is v₀ % D ∈ [0, D)
  set vs : ℤ := v₀ % D with hvs
  have hvsnn : 0 ≤ vs := Int.emod_nonneg _ (by omega)
  have hvslt : vs < D := Int.emod_lt_of_pos _ hDpos
  have hstar : ((m₀ + j₀) * w₁) % w₂ = vs := by
    rw [walk_formula w₁ w₂ D hD m₀ j₀, ← hv₀]
    have : v₀ - j₀ * D = vs := by
      rw [hvs, hj₀, Int.emod_def]; ring
    rw [this]
    exact Int.emod_eq_of_lt hvsnn (by omega)
  intro k hk hkT
  rcases eq_or_lt_of_le hk with hk0 | hkpos
  · -- k = 0: bottom hit, r = -vs, partner n = the Euclidean quotient
    refine ⟨((m₀ + j₀) * w₁) / w₂, ?_⟩
    have hdiv : ((m₀ + j₀) * w₁) / w₂ * w₂ = (m₀ + j₀) * w₁ - vs := by
      rw [← hstar]
      have := Int.emod_def ((m₀ + j₀) * w₁) w₂
      linarith [this]
    rw [← hk0]
    have : ((m₀ + j₀) * w₁) / w₂ * w₂ - (m₀ + j₀ + 0) * w₁ = -vs := by
      rw [hdiv]; ring
    rw [this, abs_neg, abs_of_nonneg hvsnn]
    omega
  · -- k ≥ 1: top hit, r = k·D - vs ∈ (0, T], partner = quotient + 1
    have hwk : ((m₀ + j₀ + k) * w₁) % w₂ = vs - k * D + w₂ := by
      have h1 : ((m₀ + j₀ + k) * w₁) % w₂ = (vs - k * D) % w₂ := by
        have := walk_formula w₁ w₂ D hD (m₀ + j₀) k
        rw [hstar] at this
        convert this using 3
        ring
      rw [h1]
      have h2 : (vs - k * D) % w₂ = (vs - k * D + w₂) % w₂ := (Int.add_emod_self).symm
      rw [h2]
      apply Int.emod_eq_of_lt
      · omega
      · omega
    refine ⟨((m₀ + j₀ + k) * w₁) / w₂ + 1, ?_⟩
    have hdiv : ((m₀ + j₀ + k) * w₁) / w₂ * w₂
        = (m₀ + j₀ + k) * w₁ - (vs - k * D + w₂) := by
      rw [← hwk]
      have := Int.emod_def ((m₀ + j₀ + k) * w₁) w₂
      linarith [this]
    have hr : (((m₀ + j₀ + k) * w₁) / w₂ + 1) * w₂ - (m₀ + j₀ + k) * w₁
        = k * D - vs := by
      have expand : (((m₀ + j₀ + k) * w₁) / w₂ + 1) * w₂
          = ((m₀ + j₀ + k) * w₁) / w₂ * w₂ + w₂ := by ring
      rw [expand, hdiv]; ring
    have hkD : D ≤ k * D := by
      have hk1 : (1 : ℤ) ≤ k := by omega
      calc D = 1 * D := (one_mul D).symm
        _ ≤ k * D := mul_le_mul_of_nonneg_right hk1 hDpos.le
    rw [hr, abs_of_nonneg (by omega)]
    omega

end SpreadPairFloor
end LonelyRunner
