/-
  TournamentH7.LRCScaleSeparation  (opus-2026-07-03-S50, THM-608)

  Lean formalization of mac-mini-S23's SCALE-SEPARATION / CLUSTER-ABSORPTION LEMMA — the
  rigorous single-step core of the deep-cluster renormalization (HYP-3901) and the hinge of
  the renormalization-depth architecture (HYP-4041), the LARGE-magnitude side of kps's
  two-sided `lrc14_of_magnitude_split`.

  If a base `R` is lonely at `t₀` with slack `δ` (‖r t₀‖ ≥ 1/14 + δ, V = max|R|), and a fast
  near-equal cluster `C ⊆ [N, N+D]` satisfies (i) `2δN ≥ V` and (ii) `D·(t₀+δ/V) < 6/7`, then
  the whole family `R ∪ C` is lonely — the cluster's magnitude `N` enters ONLY through (i),
  and larger `N` makes it EASIER (the fast phase). This is why the aligned band-blockers that
  defeat the bounded-denominator census are nonetheless lonely (they are `R ∪ C` with fast `C`).

  ELEMENTARY, IVT-FREE PROOF (the linear sweep phase is solved explicitly):
   * window `W = [t₀−η, t₀+η]`, `η = δ/V`; `‖·‖` is 1-Lipschitz ⟹ `R` safe on all of `W`;
   * the phase `N·t` is linear, so the EXACT witness `t* = (k + 1/14)/N`, `k = ⌈N(t₀−η) − 1/14⌉`,
     lands `t* ∈ W` with `fract(N t*) = 1/14` (the interval has length `N·2η ≥ 1` by (i), so it
     contains such a `k`);
   * every cluster phase `c t* = N t* + (c−N)t*` stays in `[1/14, 13/14]` by (ii) (no wrap).
  Kernel-pure (no native_decide). Purely real-analytic; builds independently of the census packs.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace ScaleSeparation

/-- **THM-608 (Scale-Separation / Cluster-Absorption).** Base `R` lonely at `t₀` with slack `δ`
(`V = max|R|`), fast near-equal cluster `C ⊆ [N, N+D]`, `N > 0`. Conditions (i) `2δN ≥ V`,
(ii) `D·(t₀+δ/V) < 6/7`, and the window sits in `[0,∞)` (`δ/V ≤ t₀`). Then there is a common
safe time `t` for every speed in `R ∪ C`: the family is `1/14`-lonely. -/
theorem scale_separation
    (t0 δ V : ℝ) (hδ : 0 < δ) (hV : 0 < V)
    (R : List ℤ)
    (hRV : ∀ r ∈ R, |(r : ℝ)| ≤ V)
    (hRsafe : ∀ r ∈ R, ∀ m : ℤ, (1 : ℝ) / 14 + δ ≤ |(r : ℝ) * t0 - m|)
    (N : ℤ) (hN : 0 < N) (D : ℝ) (hD : 0 ≤ D)
    (C : List ℤ)
    (hClo : ∀ c ∈ C, (N : ℝ) ≤ (c : ℝ))
    (hChi : ∀ c ∈ C, (c : ℝ) ≤ (N : ℝ) + D)
    (hi : V ≤ 2 * δ * (N : ℝ))
    (hii : D * (t0 + δ / V) < 6 / 7)
    (hwin : δ / V ≤ t0) :
    ∃ t : ℝ,
      (∀ r ∈ R, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(r : ℝ) * t - m|) ∧
      (∀ c ∈ C, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(c : ℝ) * t - m|) := by
  set η : ℝ := δ / V with hηdef
  have hη0 : 0 < η := div_pos hδ hV
  have hVη : V * η = δ := by rw [hηdef]; field_simp
  set A : ℝ := (N : ℝ) with hAdef
  have hA0 : 0 < A := by rw [hAdef]; exact_mod_cast hN
  have hs0 : (0 : ℝ) ≤ t0 - η := by linarith [hwin]
  -- (i) ⟹ the sweep window covers a full period
  have hlen : (1 : ℝ) ≤ A * (2 * η) := by
    have hstep : A * (2 * η) = 2 * δ * A / V := by rw [hηdef]; ring
    rw [hstep, le_div_iff₀ hV]
    have : 2 * δ * A = 2 * δ * (N : ℝ) := by rw [hAdef]
    linarith [hi, this]
  -- the integer offset of the witness
  set k : ℤ := ⌈A * (t0 - η) - 1 / 14⌉ with hkdef
  have hk_lo : A * (t0 - η) - 1 / 14 ≤ (k : ℝ) := Int.le_ceil _
  have hk_hi : (k : ℝ) < A * (t0 - η) - 1 / 14 + 1 := Int.ceil_lt_add_one _
  -- the explicit witness
  set t : ℝ := ((k : ℝ) + 1 / 14) / A with htdef
  have hAt : A * t = (k : ℝ) + 1 / 14 := by
    rw [htdef, mul_div_cancel₀ _ (ne_of_gt hA0)]
  have ht_ge : t0 - η ≤ t := by
    rw [htdef, le_div_iff₀ hA0]; nlinarith [hk_lo]
  have ht_le : t ≤ t0 + η := by
    rw [htdef, div_le_iff₀ hA0]
    have hexp : (t0 + η) * A = A * (t0 - η) + A * (2 * η) := by ring
    nlinarith [hk_hi, hlen, hexp]
  have ht0 : |t - t0| ≤ η := by rw [abs_le]; constructor <;> linarith [ht_ge, ht_le]
  have ht_nonneg : 0 ≤ t := by linarith [ht_ge, hs0]
  refine ⟨t, ?_, ?_⟩
  · -- R stays safe on the whole window (1-Lipschitz)
    intro r hr m
    have hrV := hRV r hr
    have hsafe := hRsafe r hr m
    have htri : |(r : ℝ) * t0 - m| - |(r : ℝ) * t0 - (r : ℝ) * t| ≤ |(r : ℝ) * t - m| := by
      have h := abs_sub_abs_le_abs_sub ((r : ℝ) * t0 - (m : ℝ)) ((r : ℝ) * t - (m : ℝ))
      have he : |((r : ℝ) * t0 - (m : ℝ)) - ((r : ℝ) * t - (m : ℝ))|
          = |(r : ℝ) * t0 - (r : ℝ) * t| := by congr 1; ring
      rw [he] at h; linarith [h]
    have hdist : |(r : ℝ) * t0 - (r : ℝ) * t| ≤ δ := by
      rw [show (r : ℝ) * t0 - (r : ℝ) * t = (r : ℝ) * (t0 - t) by ring, abs_mul]
      have h1 : |(r : ℝ)| * |t0 - t| ≤ V * η := by
        apply mul_le_mul hrV _ (abs_nonneg _) hV.le
        rw [abs_sub_comm]; exact ht0
      rw [hVη] at h1; exact h1
    linarith [htri, hsafe, hdist]
  · -- the cluster fits the safe band (fast phase + (ii), no wrap)
    intro c hc
    have hclo := hClo c hc
    have hchi := hChi c hc
    have hcN : (0 : ℝ) ≤ (c : ℝ) - A := by rw [hAdef]; linarith [hclo]
    have hcND : (c : ℝ) - A ≤ D := by rw [hAdef]; linarith [hchi]
    have hoff_nonneg : 0 ≤ ((c : ℝ) - A) * t := mul_nonneg hcN ht_nonneg
    have hoff_ub : ((c : ℝ) - A) * t ≤ D * (t0 + η) := by
      calc ((c : ℝ) - A) * t ≤ D * t := mul_le_mul_of_nonneg_right hcND ht_nonneg
        _ ≤ D * (t0 + η) := mul_le_mul_of_nonneg_left ht_le hD
    have hDtmax : D * (t0 + η) < 6 / 7 := by rw [hηdef]; exact hii
    have hct : (c : ℝ) * t = (1 / 14 + ((c : ℝ) - A) * t) + (k : ℝ) := by
      have hsplit : (c : ℝ) * t = A * t + ((c : ℝ) - A) * t := by ring
      rw [hsplit, hAt]; ring
    have hfr : Int.fract ((c : ℝ) * t) = 1 / 14 + ((c : ℝ) - A) * t := by
      rw [hct, Int.fract_add_intCast]
      apply Int.fract_eq_self.mpr
      constructor
      · linarith [hoff_nonneg]
      · linarith [hoff_ub, hDtmax]
    rw [far_iff_fract, hfr]
    constructor
    · linarith [hoff_nonneg]
    · linarith [hoff_ub, hDtmax]

/-- **THM-608, family form** — directly pluggable into the endgame `Lonely 14` vocabulary.
A 13-speed family `v` whose every speed lies in the lonely base `R` or the fast near-equal
cluster `C` (under the same (i),(ii),window hypotheses) is `Lonely 14`. This is the renormalization
step as a loneliness-transport: base lonely (with slack) + fast near-equal cluster ⟹ whole family
lonely. -/
theorem lonely_of_scale_separation
    (t0 δ V : ℝ) (hδ : 0 < δ) (hV : 0 < V)
    (R : List ℤ)
    (hRV : ∀ r ∈ R, |(r : ℝ)| ≤ V)
    (hRsafe : ∀ r ∈ R, ∀ m : ℤ, (1 : ℝ) / 14 + δ ≤ |(r : ℝ) * t0 - m|)
    (N : ℤ) (hN : 0 < N) (D : ℝ) (hD : 0 ≤ D)
    (C : List ℤ)
    (hClo : ∀ c ∈ C, (N : ℝ) ≤ (c : ℝ))
    (hChi : ∀ c ∈ C, (c : ℝ) ≤ (N : ℝ) + D)
    (hi : V ≤ 2 * δ * (N : ℝ))
    (hii : D * (t0 + δ / V) < 6 / 7)
    (hwin : δ / V ≤ t0)
    {ι : Type*} (v : ι → ℤ) (hcover : ∀ i, v i ∈ R ∨ v i ∈ C) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨t, hRs, hCs⟩ :=
    scale_separation t0 δ V hδ hV R hRV hRsafe N hN D hD C hClo hChi hi hii hwin
  refine ⟨t, fun i m => ?_⟩
  have key : (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - (m : ℝ)| := by
    rcases hcover i with hR | hC
    · exact hRs (v i) hR m
    · exact hCs (v i) hC m
  simpa using key

/-- **Slack from LRC(13)** — how the LRC(≤13) citation supplies scale_separation's base slack `δ`.
A family lonely at the `1/13` band has slack `1/182` at the `1/14` band, since `1/13 = 1/14 + 1/182`.
A `≤ 12`-speed base is `Lonely 13` by the citation node (`LRCUpTo13` with `k+1 = 13`), hence lonely at
`1/14` with slack `δ = 1/182 > 0` — exactly the `hRsafe` hypothesis of `scale_separation`. So the
renormalization step composes: (LRC(≤13) base) + (fast near-equal cluster) ⟹ `Lonely 14`. -/
theorem slack_of_lonely13 {ι : Type*} (w : ι → ℤ) (t0 : ℝ) (h : Lonely 13 w t0) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / 14 + 1 / 182 ≤ |(w i : ℝ) * t0 - (m : ℝ)| := by
  intro i m
  have hi := h i m
  have hcast : (1 : ℝ) / 14 + 1 / 182 = (1 : ℝ) / (13 : ℕ) := by norm_num
  rw [hcast]; exact hi

/-- Loneliness is invariant under an integer shift of `t` (integer speeds). -/
theorem lonely_add_int {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) (kk : ℤ) :
    Lonely n v (t + (kk : ℝ)) ↔ Lonely n v t := by
  unfold Lonely
  constructor
  · intro h i m
    have hm := h i (m + v i * kk)
    have he : |(v i : ℝ) * (t + (kk : ℝ)) - ((m + v i * kk : ℤ) : ℝ)|
        = |(v i : ℝ) * t - (m : ℝ)| := by push_cast; ring_nf
    rw [he] at hm; exact hm
  · intro h i m
    have hm := h i (m - v i * kk)
    have he : |(v i : ℝ) * (t + (kk : ℝ)) - (m : ℝ)|
        = |(v i : ℝ) * t - ((m - v i * kk : ℤ) : ℝ)| := by push_cast; ring_nf
    rw [he]; exact hm

/-- **(b) THE RENORMALIZATION STEP, COMPOSED WITH THE LRC(≤13) CITATION** — discharges the
fast-near-equal-cluster slice of `lrc14_of_magnitude_split`'s `hlarge`. A family `v` splitting into a
`≤12`-speed base (closed wholesale by `LRCUpTo13`, giving slack `δ = 1/(m+1) − 1/14 > 0` for free) and a
fast (i) near-equal (ii) cluster `C ⊆ [N, N+D]` is `Lonely 14`. The base point is shifted into `[1,2)` by
periodicity so the sweep window sits in `[0,∞)`; hence (ii) is stated with the a-priori bound `t₀' < 2`. -/
theorem lonely14_of_citation_cluster
    (cite : LRCUpTo13)
    {ι : Type*} (v : ι → ℤ)
    (m : ℕ) (hm : m ≤ 12) (base : Fin m → ℤ) (hbase0 : ∀ i, base i ≠ 0)
    (V : ℝ) (hV : 0 < V) (hbaseV : ∀ i, |(base i : ℝ)| ≤ V)
    (hη1 : (1 / ((m : ℝ) + 1) - 1 / 14) / V ≤ 1)
    (N : ℤ) (hN : 0 < N) (D : ℝ) (hD : 0 ≤ D) (cluster : List ℤ)
    (hClo : ∀ c ∈ cluster, (N : ℝ) ≤ (c : ℝ))
    (hChi : ∀ c ∈ cluster, (c : ℝ) ≤ (N : ℝ) + D)
    (hcover : ∀ i, (∃ j, base j = v i) ∨ v i ∈ cluster)
    (hi : V ≤ 2 * (1 / ((m : ℝ) + 1) - 1 / 14) * (N : ℝ))
    (hii : D * (2 + (1 / ((m : ℝ) + 1) - 1 / 14) / V) < 6 / 7) :
    ∃ t : ℝ, Lonely 14 v t := by
  set δ : ℝ := 1 / ((m : ℝ) + 1) - 1 / 14 with hδdef
  have hm14 : ((m : ℝ) + 1) < 14 := by
    have : (m : ℝ) ≤ 12 := by exact_mod_cast hm
    linarith
  have hmpos : (0 : ℝ) < (m : ℝ) + 1 := by positivity
  have hδ0 : 0 < δ := by
    rw [hδdef]
    have h1 : (1 : ℝ) / 14 < 1 / ((m : ℝ) + 1) := one_div_lt_one_div_of_lt hmpos hm14
    linarith
  -- base lonely point from the citation, shifted into [1,2)
  obtain ⟨t0, ht0⟩ := cite m hm base hbase0
  set t0' : ℝ := t0 + ((1 - ⌊t0⌋ : ℤ) : ℝ) with ht0'def
  have ht0'lonely : Lonely (m + 1) base t0' :=
    (lonely_add_int (m + 1) base t0 (1 - ⌊t0⌋)).mpr ht0
  have hfloor_le : (⌊t0⌋ : ℝ) ≤ t0 := Int.floor_le t0
  have hlt_floor : t0 < (⌊t0⌋ : ℝ) + 1 := Int.lt_floor_add_one t0
  have ht0'ge1 : (1 : ℝ) ≤ t0' := by rw [ht0'def]; push_cast; linarith
  have ht0'lt2 : t0' < 2 := by rw [ht0'def]; push_cast; linarith
  -- slack: Lonely (m+1) ⟹ 1/14 + δ ≤ |·| (since 1/(m+1) = 1/14 + δ)
  have hslack : ∀ i, ∀ mm : ℤ, (1 : ℝ) / 14 + δ ≤ |(base i : ℝ) * t0' - (mm : ℝ)| := by
    intro i mm
    have hi0 := ht0'lonely i mm
    have hcast : (1 : ℝ) / 14 + δ = (1 : ℝ) / ((m + 1 : ℕ)) := by
      rw [hδdef]; push_cast; ring
    rw [hcast]; exact hi0
  -- window and (ii) at t0'
  set η : ℝ := δ / V with hηdef
  have hwin : δ / V ≤ t0' := le_trans (by rw [← hηdef]; exact hη1) ht0'ge1
  have hii' : D * (t0' + δ / V) < 6 / 7 := by
    have hle : D * (t0' + δ / V) ≤ D * (2 + δ / V) :=
      mul_le_mul_of_nonneg_left (by linarith [ht0'lt2]) hD
    linarith [hii, hle]
  -- assemble via the family form (base as a list, plus the cluster)
  refine lonely_of_scale_separation t0' δ V hδ0 hV (List.ofFn base) ?_ ?_
    N hN D hD cluster hClo hChi hi hii' hwin v ?_
  · intro r hr
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hr
    exact hbaseV i
  · intro r hr mm
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hr
    exact hslack i mm
  · intro i
    rcases hcover i with ⟨j, hj⟩ | hc
    · exact Or.inl (List.mem_ofFn.mpr ⟨j, hj⟩)
    · exact Or.inr hc

end ScaleSeparation
end LonelyRunner

