/-
  TournamentH7.LRCWitnessCert — the witness-arc certificate layer.
  kind-pasteur-2026-07-02-S1 (HYP-3958; the (R)-node of HYP-3953/THM-599 made formal).

  THE REASONING IMPROVEMENT.  The witness route's formal interface needs NO integral, NO average,
  NO measure of unions: for a shape (P, offsets) a certificate is ONE rational arc [lo, hi] and a
  ruler phase c such that every small speed is safe on the arc (target 0) and every offset is safe
  on the arc (target c) — each a DECIDABLE rational-arithmetic check.  Then for EVERY reference
  speed V with V·(hi − lo) > 1 the ruler time τ = (j + c)/V lands in the arc for some integer j,
  and the EXACT c-ruler identity ‖(V − o)τ‖ = ‖o·τ − c‖ converts offset-safety into loneliness of
  the fast runners.  One certificate covers the whole V-tail.  The A(U)/Fubini layer (THM-599)
  remains the mathematics that FINDS certificates; it is no longer an obligation of the proof.

  Contents (all sorry-free):
   * `coe_add_int`, `cruler_identity` — the exact identity (pure algebra on `UnitAddCircle`);
   * `arcSafe` — the decidable per-speed safety check; `norm_ge_of_arcSafe` — its semantics;
   * `exists_int_in_long_interval` — the V-tail pigeonhole;
   * `cert_lonely_tail` — the glue: a checked certificate + V·(hi − lo) > 1 gives a time at which
     every listed runner is ≥ h from the origin.
-/
import Mathlib.Analysis.Normed.Group.AddCircle

namespace LonelyRunner
namespace WitnessCert

noncomputable section

open scoped Real

/-- Adding an integer does not move a point of the unit circle. -/
theorem coe_add_int (x : ℝ) (j : ℤ) : ((x + j : ℝ) : UnitAddCircle) = (x : UnitAddCircle) := by
  have hj : ((j : ℝ) : UnitAddCircle) = 0 := by
    have : ((j : ℝ)) ∈ AddSubgroup.zmultiples (1 : ℝ) := by
      refine ⟨j, ?_⟩
      simp
    simpa [AddCircle] using (QuotientAddGroup.eq_zero_iff _).mpr this
  calc ((x + j : ℝ) : UnitAddCircle) = ((x : ℝ) : UnitAddCircle) + ((j : ℝ) : UnitAddCircle) := rfl
    _ = (x : UnitAddCircle) := by rw [hj, add_zero]

/-- **The exact c-ruler identity**: at the ruler time `τ = (j + c)/V`, the runner at speed
`V − o` sits exactly where `c − o·τ` sits on the circle. -/
theorem cruler_identity (V o j : ℤ) (c : ℝ) (hV : (V : ℝ) ≠ 0) :
    ((((V - o : ℤ) : ℝ) * (((j : ℝ) + c)/V) : ℝ) : UnitAddCircle)
      = ((c - (o : ℝ) * (((j : ℝ) + c)/V) : ℝ) : UnitAddCircle) := by
  have key : ((V - o : ℤ) : ℝ) * (((j : ℝ) + c)/V)
      = (c - (o : ℝ) * (((j : ℝ) + c)/V)) + j := by
    push_cast
    field_simp
    ring
  rw [key, coe_add_int]

/-- Norms agree at the two presentations (corollary form used by the glue). -/
theorem cruler_norm (V o j : ℤ) (c : ℝ) (hV : (V : ℝ) ≠ 0) :
    ‖((((V - o : ℤ) : ℝ) * (((j : ℝ) + c)/V) : ℝ) : UnitAddCircle)‖
      = ‖((c - (o : ℝ) * (((j : ℝ) + c)/V) : ℝ) : UnitAddCircle)‖ := by
  rw [cruler_identity V o j c hV]

end

/-! ## The decidable safety check -/

/-- `arcSafe h s cTar lo hi`: for every `x ∈ [lo, hi]`, the point `s·x − cTar` stays at distance
`≥ h` from every integer — checked by fitting `[s·lo − cTar, s·hi − cTar]` inside
`[k + h, k + 1 − h]` for `k = ⌊s·lo − cTar⌋`.  (With `s = 0` this checks `‖cTar‖ ≥ h`, which is
exactly the reference runner's condition.) -/
def arcSafe (h : ℚ) (s : ℤ) (cTar lo hi : ℚ) : Prop :=
  (s : ℚ) * lo - cTar ≤ (s : ℚ) * hi - cTar ∧
  h ≤ ((s : ℚ) * lo - cTar) - ⌊(s : ℚ) * lo - cTar⌋ ∧
  (s : ℚ) * hi - cTar ≤ (⌊(s : ℚ) * lo - cTar⌋ : ℚ) + 1 - h

instance (h : ℚ) (s : ℤ) (cTar lo hi : ℚ) : Decidable (arcSafe h s cTar lo hi) := by
  unfold arcSafe; infer_instance

/-- Distance to the integers is at least `h` on `[k + h, k + 1 − h]`. -/
theorem norm_ge_of_band {y : ℝ} {k : ℤ} {h : ℝ} (h0 : 0 ≤ h)
    (hlo : (k : ℝ) + h ≤ y) (hhi : y ≤ (k : ℝ) + 1 - h) :
    h ≤ ‖((y : ℝ) : UnitAddCircle)‖ := by
  rw [UnitAddCircle.norm_eq]
  rcases le_or_gt (round y) k with hr | hr
  · have : (round y : ℝ) ≤ (k : ℝ) := by exact_mod_cast hr
    have : h ≤ y - round y := by linarith
    calc h ≤ y - round y := this
      _ ≤ |y - round y| := le_abs_self _
  · have hk1 : (k : ℤ) + 1 ≤ round y := hr
    have : ((k : ℤ) + 1 : ℝ) ≤ (round y : ℝ) := by exact_mod_cast hk1
    have hneg : y - round y ≤ -h := by push_cast at this; linarith
    calc h ≤ -(y - round y) := by linarith
      _ ≤ |y - round y| := neg_le_abs _

/-- **Semantics of `arcSafe`**: a checked speed keeps distance `≥ h` from the integers at every
real time in the arc (for `s ≥ 0`). -/
theorem norm_ge_of_arcSafe {h : ℚ} {s : ℤ} {cTar lo hi : ℚ}
    (hs : 0 ≤ s) (h0 : 0 ≤ h) (hsafe : arcSafe h s cTar lo hi)
    {x : ℝ} (hxlo : (lo : ℝ) ≤ x) (hxhi : x ≤ (hi : ℝ)) :
    (h : ℝ) ≤ ‖(((s : ℝ) * x - (cTar : ℝ) : ℝ) : UnitAddCircle)‖ := by
  obtain ⟨-, hA, hB⟩ := hsafe
  set k : ℤ := ⌊(s : ℚ) * lo - cTar⌋ with hkdef
  -- the two rational band facts
  have q1 : (k : ℚ) + h ≤ (s : ℚ) * lo - cTar := by linarith
  have q2 : (s : ℚ) * hi - cTar ≤ (k : ℚ) + 1 - h := hB
  -- cast to ℝ (component form)
  have q1R : (k : ℝ) + (h : ℝ) ≤ (s : ℝ) * (lo : ℝ) - (cTar : ℝ) := by exact_mod_cast q1
  have q2R : (s : ℝ) * (hi : ℝ) - (cTar : ℝ) ≤ (k : ℝ) + 1 - (h : ℝ) := by exact_mod_cast q2
  have hsR : (0 : ℝ) ≤ (s : ℝ) := by exact_mod_cast hs
  have hmul1 : (s : ℝ) * (lo : ℝ) ≤ (s : ℝ) * x := mul_le_mul_of_nonneg_left hxlo hsR
  have hmul2 : (s : ℝ) * x ≤ (s : ℝ) * (hi : ℝ) := mul_le_mul_of_nonneg_left hxhi hsR
  exact norm_ge_of_band (k := k) (by exact_mod_cast h0) (by linarith) (by linarith)

/-! ## The V-tail pigeonhole and the glue -/

/-- An interval of length `> 1` contains an integer. -/
theorem exists_int_in_long_interval {a b : ℝ} (hab : a + 1 ≤ b) :
    ∃ j : ℤ, a ≤ (j : ℝ) ∧ (j : ℝ) ≤ b := by
  refine ⟨⌈a⌉, Int.le_ceil a, ?_⟩
  have h1 : (⌈a⌉ : ℝ) < a + 1 := Int.ceil_lt_add_one a
  linarith

/-- **The witness-arc certificate glue.**  If every small speed `s ∈ P` is `arcSafe` at target 0
and every offset `o ∈ offs` is `arcSafe` at target `c` on the SAME rational arc `[lo, hi]`, then
for EVERY integer reference speed `V` with `(V : ℚ) * (hi − lo) > 1` there is a real time `τ` at
which every small runner and every fast runner `V − o` is at circle-distance `≥ h` from the
origin.  One decidable certificate covers the whole `V`-tail. -/
theorem cert_lonely_tail {h : ℚ} (h0 : 0 ≤ h) {cTar lo hi : ℚ}
    (P offs : List ℤ)
    (hPpos : ∀ s ∈ P, 0 ≤ s) (hOpos : ∀ o ∈ offs, 0 ≤ o)
    (hPsafe : ∀ s ∈ P, arcSafe h s 0 lo hi)
    (hOsafe : ∀ o ∈ offs, arcSafe h o cTar lo hi)
    (V : ℤ) (hVlen : 1 < (V : ℚ) * (hi - lo)) (hVpos : 0 < V) :
    ∃ τ : ℝ,
      (∀ s ∈ P, (h : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs, (h : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hVpos
  have hVQ : (V : ℝ) ≠ 0 := ne_of_gt hVR
  -- find the ruler index: an integer in [V·lo − c, V·hi − c]
  obtain ⟨j, hj1, hj2⟩ := exists_int_in_long_interval
    (a := (V : ℝ) * (lo : ℝ) - (cTar : ℝ)) (b := (V : ℝ) * (hi : ℝ) - (cTar : ℝ))
    (by
      have : (1 : ℝ) < (V : ℝ) * ((hi : ℝ) - (lo : ℝ)) := by exact_mod_cast hVlen
      nlinarith)
  set τ : ℝ := ((j : ℝ) + (cTar : ℝ)) / (V : ℝ) with hτ
  have hτlo : (lo : ℝ) ≤ τ := by
    rw [hτ, le_div_iff₀ hVR]
    linarith
  have hτhi : τ ≤ (hi : ℝ) := by
    rw [hτ, div_le_iff₀ hVR]
    linarith
  refine ⟨τ, ?_, ?_⟩
  · intro s hs
    have := norm_ge_of_arcSafe (hPpos s hs) h0 (hPsafe s hs) hτlo hτhi
    simpa using this
  · intro o ho
    have hsafe := norm_ge_of_arcSafe (hOpos o ho) h0 (hOsafe o ho) hτlo hτhi
    have hid := cruler_norm V o j (cTar : ℝ) hVQ
    rw [hτ]
    rw [hid]
    have hswap : ‖(((cTar : ℝ) - (o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) : ℝ) : UnitAddCircle)‖
        = ‖(((o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) - (cTar : ℝ) : ℝ) : UnitAddCircle)‖ := by
      have : ((cTar : ℝ) - (o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) : ℝ)
          = -(((o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) - (cTar : ℝ))) := by ring
      rw [this]
      have hneg : ((-(((o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) - (cTar : ℝ))) : ℝ) : UnitAddCircle)
          = -((((o : ℝ) * (((j : ℝ) + (cTar : ℝ))/(V : ℝ)) - (cTar : ℝ) : ℝ) : UnitAddCircle)) := rfl
      rw [hneg, norm_neg]
    rw [hswap]
    have := hsafe
    rw [hτ] at this
    simpa using this

/-! ## Certificate rows (generated by the exact Python engine; checked here by `native_decide`)

Each row: a shape `(P, offs)`, a ruler phase `c`, and ONE safe arc `[lo, hi]`; the tail theorem
then covers EVERY reference speed `V ≥ V*`.  Rows 3 and 4 are genuine 13-runner families
(|P| + |offs| = 13): infinite families of LRC(14) instances, machine-checked. -/

def h14 : ℚ := 1/14

/-- Row 2 — the consecutive (AP) cluster `{V, V−1, …, V−12}`. -/
def offsAP : List ℤ := [0,1,2,3,4,5,6,7,8,9,10,11,12]

theorem certAP_safe : ∀ o ∈ offsAP, arcSafe h14 o (59/800) (62413/67200) 1 := by
  native_decide

/-- **For every `V ≥ 15`, the thirteen runners `V, V−1, …, V−12` share a `1/14`-lonely time.**
(The dilated-AP top-cluster family, formally covered by one certificate.) -/
theorem certAP_tail (V : ℤ) (hV : 15 ≤ V) :
    ∃ τ : ℝ, ∀ o ∈ offsAP, ((h14 : ℚ) : ℝ) ≤
      ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  have hVQ : (15 : ℚ) ≤ (V : ℚ) := by exact_mod_cast hV
  obtain ⟨τ, -, hO⟩ := cert_lonely_tail (h := h14) (by norm_num [h14])
    ([] : List ℤ) offsAP
    (by intro s hs; simp at hs) (by native_decide)
    (by intro s hs; simp at hs) certAP_safe
    V (by norm_num; nlinarith) (by omega)
  exact ⟨τ, hO⟩

/-- Row 3 — the admissible 13-runner family `{1,2,3} ∪ {V − o : o ∈ offs3}` (3 + 10 = 13). -/
def P3 : List ℤ := [1,2,3]
def offs3 : List ℤ := [0,1,2,3,5,8,11,14,17,20]

theorem cert3_safeP : ∀ s ∈ P3, arcSafe h14 s 0 (29/42) (79729/112000) := by native_decide
theorem cert3_safeO : ∀ o ∈ offs3, arcSafe h14 o (247/800) (29/42) (79729/112000) := by
  native_decide

/-- **An infinite 13-runner LRC(14) family, machine-checked**: for every `V ≥ 47`, the runners
`{1, 2, 3, V−20, V−17, V−14, V−11, V−8, V−5, V−3, V−2, V−1, V}` share a `1/14`-lonely time. -/
theorem cert3_tail (V : ℤ) (hV : 47 ≤ V) :
    ∃ τ : ℝ,
      (∀ s ∈ P3, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs3, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  have hVQ : (47 : ℚ) ≤ (V : ℚ) := by exact_mod_cast hV
  exact cert_lonely_tail (h := h14) (by norm_num [h14]) P3 offs3
    (by native_decide) (by native_decide)
    cert3_safeP cert3_safeO
    V (by norm_num; nlinarith) (by omega)

/-- Row 4 — the admissible 13-runner family `{2,3,4,5,6} ∪ {V − o : o even ≤ 14}` (5 + 8 = 13). -/
def P4 : List ℤ := [2,3,4,5,6]
def offs4 : List ℤ := [0,2,4,6,8,10,12,14]

theorem cert4_safeP : ∀ s ∈ P4, arcSafe h14 s 0 (34413/78400) (13/28) := by native_decide
theorem cert4_safeO : ∀ o ∈ offs4, arcSafe h14 o (59/800) (34413/78400) (13/28) := by
  native_decide

/-- **A second infinite 13-runner LRC(14) family, machine-checked**: for every `V ≥ 40`, the
runners `{2,3,4,5,6} ∪ {V, V−2, V−4, …, V−14}` share a `1/14`-lonely time. -/
theorem cert4_tail (V : ℤ) (hV : 40 ≤ V) :
    ∃ τ : ℝ,
      (∀ s ∈ P4, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs4, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  have hVQ : (40 : ℚ) ≤ (V : ℚ) := by exact_mod_cast hV
  exact cert_lonely_tail (h := h14) (by norm_num [h14]) P4 offs4
    (by native_decide) (by native_decide)
    cert4_safeP cert4_safeO
    V (by norm_num; nlinarith) (by omega)

/-! ## The two-level nested certificate (HYP-3959 — the bridging lemma's first instance)

Two scales compose: level 1 (reference `V₁`, offsets `offs₁`) is certified at the INFLATED band
`h + μ`; level 2 (reference `V₂ > V₁/μ`, offsets `offs₂`) samples the window `[τ₀, τ₀ + μ/V₁]`
around a level-1 ruler point.  At the final time τ (a level-2 ruler point):
 * `P` and `offs₂` are exact (τ ∈ [lo, hi]; the c-ruler identity at `V₂`);
 * only the `V₁`-PHASE drifts, by at most `V₁ · (μ/V₁) = μ` — absorbed by the inflated band
   (the circle norm is 1-Lipschitz).
Iterating this step is the `n`-uniform ladder; no other analytic input appears. -/

/-- The circle norm is 1-Lipschitz against real perturbation. -/
theorem norm_ge_norm_sub_abs (y x : ℝ) :
    ‖((y : ℝ) : UnitAddCircle)‖ - |x - y| ≤ ‖((x : ℝ) : UnitAddCircle)‖ := by
  have hxy : ((y : ℝ) : UnitAddCircle)
      = ((x : ℝ) : UnitAddCircle) + ((y - x : ℝ) : UnitAddCircle) := by
    have hsum : (x : ℝ) + (y - x) = y := by ring
    calc ((y : ℝ) : UnitAddCircle) = (((x + (y - x) : ℝ)) : UnitAddCircle) := by rw [hsum]
      _ = ((x : ℝ) : UnitAddCircle) + ((y - x : ℝ) : UnitAddCircle) := rfl
  have hle : ‖((y : ℝ) : UnitAddCircle)‖
      ≤ ‖((x : ℝ) : UnitAddCircle)‖ + ‖((y - x : ℝ) : UnitAddCircle)‖ := by
    rw [hxy]; exact norm_add_le _ _
  have habs : ‖((y - x : ℝ) : UnitAddCircle)‖ ≤ |y - x| := by
    rw [UnitAddCircle.norm_eq]
    simpa using round_le (y - x) 0
  have habs' : |y - x| = |x - y| := abs_sub_comm y x
  linarith

/-- **The two-level nested certificate.** -/
theorem cert_two_level {h μ : ℚ} (h0 : 0 ≤ h) (hμpos : 0 < μ) {cTar₁ cTar₂ lo hi : ℚ}
    (P offs₁ offs₂ : List ℤ)
    (hPpos : ∀ s ∈ P, 0 ≤ s) (hO1pos : ∀ o ∈ offs₁, 0 ≤ o) (hO2pos : ∀ o ∈ offs₂, 0 ≤ o)
    (hPsafe : ∀ s ∈ P, arcSafe h s 0 lo hi)
    (hO1safe : ∀ o ∈ offs₁, arcSafe (h + μ) o cTar₁ lo hi)
    (hO2safe : ∀ o ∈ offs₂, arcSafe h o cTar₂ lo hi)
    (V₁ V₂ : ℤ) (hV1pos : 0 < V₁) (hV2pos : 0 < V₂)
    (hV1len : 1 < (V₁ : ℚ) * (hi - lo - μ / V₁))
    (hV2len : 1 < (V₂ : ℚ) * (μ / V₁)) :
    ∃ τ : ℝ,
      (∀ s ∈ P, (h : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs₁, (h : ℝ) ≤ ‖((((V₁ - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs₂, (h : ℝ) ≤ ‖((((V₂ - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  have hV1R : (0 : ℝ) < (V₁ : ℝ) := by exact_mod_cast hV1pos
  have hV2R : (0 : ℝ) < (V₂ : ℝ) := by exact_mod_cast hV2pos
  have hV2Q : (V₂ : ℝ) ≠ 0 := ne_of_gt hV2R
  set δ : ℝ := (μ : ℝ) / (V₁ : ℝ) with hδdef
  have hδpos : 0 < δ := by
    apply div_pos (by exact_mod_cast hμpos) hV1R
  have hδcast : ((μ / V₁ : ℚ) : ℝ) = δ := by rw [hδdef]; push_cast; ring
  -- level-1 ruler point τ₀ ∈ [lo, hi − δ]
  obtain ⟨j₁, hj1, hj2⟩ := exists_int_in_long_interval
    (a := (V₁ : ℝ) * (lo : ℝ) - (cTar₁ : ℝ))
    (b := (V₁ : ℝ) * ((hi : ℝ) - δ) - (cTar₁ : ℝ))
    (by
      have h1 : (1 : ℝ) < (V₁ : ℝ) * ((hi : ℝ) - (lo : ℝ) - ((μ / V₁ : ℚ) : ℝ)) := by
        exact_mod_cast hV1len
      rw [hδcast] at h1
      nlinarith)
  set τ₀ : ℝ := ((j₁ : ℝ) + (cTar₁ : ℝ)) / (V₁ : ℝ) with hτ₀
  have hτ₀lo : (lo : ℝ) ≤ τ₀ := by rw [hτ₀, le_div_iff₀ hV1R]; linarith
  have hτ₀hi : τ₀ ≤ (hi : ℝ) - δ := by rw [hτ₀, div_le_iff₀ hV1R]; linarith
  have hV₁τ₀ : (V₁ : ℝ) * τ₀ = (j₁ : ℝ) + (cTar₁ : ℝ) := by
    rw [hτ₀]; field_simp
  -- level-2 ruler point τ ∈ [τ₀, τ₀ + δ]
  obtain ⟨j₂, hk1, hk2⟩ := exists_int_in_long_interval
    (a := (V₂ : ℝ) * τ₀ - (cTar₂ : ℝ)) (b := (V₂ : ℝ) * (τ₀ + δ) - (cTar₂ : ℝ))
    (by
      have h2 : (1 : ℝ) < (V₂ : ℝ) * ((μ / V₁ : ℚ) : ℝ) := by exact_mod_cast hV2len
      rw [hδcast] at h2
      nlinarith)
  set τ : ℝ := ((j₂ : ℝ) + (cTar₂ : ℝ)) / (V₂ : ℝ) with hτdef
  have hτ1 : τ₀ ≤ τ := by rw [hτdef, le_div_iff₀ hV2R]; linarith
  have hτ2 : τ ≤ τ₀ + δ := by rw [hτdef, div_le_iff₀ hV2R]; linarith
  have hτlo : (lo : ℝ) ≤ τ := le_trans hτ₀lo hτ1
  have hτhi : τ ≤ (hi : ℝ) := by linarith
  refine ⟨τ, ?_, ?_, ?_⟩
  · -- small speeds: direct
    intro s hs
    have := norm_ge_of_arcSafe (hPpos s hs) h0 (hPsafe s hs) hτlo hτhi
    simpa using this
  · -- level-1 cluster: phase-drift absorption
    intro o ho
    -- (V₁ − o)·τ = [(c₁ − o·τ) + V₁·(τ − τ₀)] + j₁
    have hkey : ((V₁ - o : ℤ) : ℝ) * τ
        = (((cTar₁ : ℝ) - (o : ℝ) * τ) + (V₁ : ℝ) * (τ - τ₀)) + (j₁ : ℝ) := by
      push_cast
      nlinarith [hV₁τ₀]
    rw [hkey, coe_add_int]
    -- the inflated-band safety of o at τ (τ ∈ [lo, hi])
    have hsafe := norm_ge_of_arcSafe (hO1pos o ho) (by linarith [hμpos.le] : (0:ℚ) ≤ h + μ)
      (hO1safe o ho) hτlo hτhi
    have hcast : (((h + μ : ℚ)) : ℝ) = (h : ℝ) + (μ : ℝ) := by push_cast; ring
    rw [hcast] at hsafe
    -- flip the sign inside the norm: ‖o·τ − c₁‖ = ‖c₁ − o·τ‖
    have hflip : ‖(((o : ℝ) * τ - (cTar₁ : ℝ) : ℝ) : UnitAddCircle)‖
        = ‖(((cTar₁ : ℝ) - (o : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
      have hne : ((cTar₁ : ℝ) - (o : ℝ) * τ : ℝ) = -(((o : ℝ) * τ - (cTar₁ : ℝ))) := by ring
      rw [hne]
      have : ((-(((o : ℝ) * τ - (cTar₁ : ℝ))) : ℝ) : UnitAddCircle)
          = -((((o : ℝ) * τ - (cTar₁ : ℝ) : ℝ) : UnitAddCircle)) := rfl
      rw [this, norm_neg]
    -- Lipschitz absorption of the V₁-phase drift (0 ≤ V₁(τ − τ₀) ≤ μ)
    have hdrift : |(((cTar₁ : ℝ) - (o : ℝ) * τ) + (V₁ : ℝ) * (τ - τ₀))
        - ((cTar₁ : ℝ) - (o : ℝ) * τ)| ≤ (μ : ℝ) := by
      have h1 : (((cTar₁ : ℝ) - (o : ℝ) * τ) + (V₁ : ℝ) * (τ - τ₀))
          - ((cTar₁ : ℝ) - (o : ℝ) * τ) = (V₁ : ℝ) * (τ - τ₀) := by ring
      rw [h1, abs_of_nonneg (by nlinarith)]
      have hδV : (V₁ : ℝ) * δ = (μ : ℝ) := by
        rw [hδdef]; field_simp
      nlinarith
    have hlip := norm_ge_norm_sub_abs ((cTar₁ : ℝ) - (o : ℝ) * τ)
      ((((cTar₁ : ℝ) - (o : ℝ) * τ) + (V₁ : ℝ) * (τ - τ₀)))
    rw [hflip] at hsafe
    linarith
  · -- level-2 cluster: exact c-ruler identity
    intro o ho
    have hid := cruler_norm V₂ o j₂ (cTar₂ : ℝ) hV2Q
    rw [hτdef, hid]
    have hsafe := norm_ge_of_arcSafe (hO2pos o ho) h0 (hO2safe o ho) hτlo hτhi
    rw [hτdef] at hsafe
    have hswap : ‖(((cTar₂ : ℝ) - (o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) : ℝ) : UnitAddCircle)‖
        = ‖(((o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) - (cTar₂ : ℝ) : ℝ) : UnitAddCircle)‖ := by
      have hne : ((cTar₂ : ℝ) - (o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) : ℝ)
          = -(((o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) - (cTar₂ : ℝ))) := by ring
      rw [hne]
      have : ((-(((o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) - (cTar₂ : ℝ))) : ℝ) : UnitAddCircle)
          = -((((o : ℝ) * (((j₂ : ℝ) + (cTar₂ : ℝ))/(V₂ : ℝ)) - (cTar₂ : ℝ) : ℝ) : UnitAddCircle)) := rfl
      rw [this, norm_neg]
    rw [hswap]
    simpa using hsafe

end WitnessCert
end LonelyRunner
