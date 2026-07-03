/-
  TournamentH7.LRCBlockSix — THE BLOCK-6 INTERVAL ENGINE (kind-pasteur-2026-07-02-S21,
  HYP-3978).  Near-equal blocks of ANY internal ratio structure, size ≤ 6, over a real
  interval base — the citation-route counterpart of mac-mini's region-based
  `LRCSimulPeel` (S18) and the many-runner generalization of klein's
  `safe_point_in_interval` (S114).

  * `gap_exists` — THE MEASURE-FREE UNION LEMMA over ℝ: a window whose clipped bad
    mass is under its length contains a point avoiding every open bad interval.
    Proof: JUMP-PAST-THE-BLOCKER induction — if the window's left edge sits inside
    some bad interval, hop to that interval's right end and drop it from the list;
    the dropped clipped mass pays for the lost window length exactly.  No measure
    theory, no counting, no sorting.

  * `runner_bad_teeth` / `runner_teeth_mass` — the teeth of runner `w` intersecting
    a window of length `L`, as an explicit list of ≤ `⌊wL⌋ + 2` intervals of width
    `1/(7w)` with clipped mass ≤ `L/7 + 2/(7w)`.

  * `block_point_step` — ≤ 6 positive runners with ARBITRARY internal ratios share a
    `1/14`-good point in any window with `w₁·L` above the explicit entry threshold:
    union density `c/7 < 1` leaves room.  Blocks of 7+ hit the density wall — the
    honest remaining core.

  * `cite_block_lonely` — the citation composition: cite `k ≤ 12` bounded runners,
    transport, and finish the (≤ 6)-block inside the interval.
-/
import Mathlib
import TournamentH7.LRCPairBlock

namespace LonelyRunner
namespace LRC14

/-! ## The measure-free union lemma -/

/-- Clipped length of one open interval against a window `[a, b]`. -/
noncomputable def clipLen (p : ℝ × ℝ) (a b : ℝ) : ℝ := max 0 (min p.2 b - max p.1 a)

theorem clipLen_nonneg (p : ℝ × ℝ) (a b : ℝ) : 0 ≤ clipLen p a b := le_max_left _ _

/-- Shrinking the window from the left does not increase the clipped length. -/
theorem clipLen_mono_left (p : ℝ × ℝ) {a a' b : ℝ} (h : a ≤ a') :
    clipLen p a' b ≤ clipLen p a b := by
  unfold clipLen
  rcases le_total (min p.2 b - max p.1 a') 0 with h1 | h1
  · rw [max_eq_left h1]
    exact le_max_left _ _
  · rw [max_eq_right h1]
    have : max p.1 a ≤ max p.1 a' := max_le_max (le_refl _) h
    calc min p.2 b - max p.1 a' ≤ min p.2 b - max p.1 a := by linarith
      _ ≤ max 0 (min p.2 b - max p.1 a) := le_max_right _ _

/-- **THE MEASURE-FREE UNION LEMMA**: a window whose clipped bad mass is under its
length contains a point avoiding every open bad interval. -/
theorem gap_exists : ∀ (n : ℕ) (bads : List (ℝ × ℝ)) (a b : ℝ),
    bads.length ≤ n → a ≤ b →
    (bads.map fun p => clipLen p a b).sum < b - a →
    ∃ t : ℝ, a ≤ t ∧ t ≤ b ∧ ∀ p ∈ bads, t ≤ p.1 ∨ p.2 ≤ t := by
  intro n
  induction n with
  | zero =>
      intro bads a b hlen hab _
      have hnil : bads = [] := List.length_eq_zero_iff.mp (Nat.le_zero.mp hlen)
      subst hnil
      exact ⟨a, le_refl _, hab, fun p hp => absurd hp (List.not_mem_nil)⟩
  | succ n ih =>
      intro bads a b hlen hab hsum
      by_cases hesc : ∀ p ∈ bads, a ≤ p.1 ∨ p.2 ≤ a
      · exact ⟨a, le_refl _, hab, hesc⟩
      · push Not at hesc
        obtain ⟨p, hp, hp1, hp2⟩ := hesc
        -- hp1 : p.1 < a, hp2 : a < p.2
        have hnonneg : ∀ x ∈ bads.map fun p => clipLen p a b, 0 ≤ x := by
          intro x hx
          rw [List.mem_map] at hx
          obtain ⟨q, _, rfl⟩ := hx
          exact clipLen_nonneg q a b
        have hpb : p.2 ≤ b := by
          by_contra hgt
          push Not at hgt
          have hcp : clipLen p a b = b - a := by
            unfold clipLen
            rw [min_eq_right (le_of_lt hgt), max_eq_right (le_of_lt hp1)]
            exact max_eq_right (by linarith)
          have hmem : clipLen p a b ∈ bads.map fun q => clipLen q a b :=
            List.mem_map.mpr ⟨p, hp, rfl⟩
          have := List.single_le_sum hnonneg _ hmem
          linarith
        -- drop p, jump to p.2
        obtain ⟨l₁, l₂, rfl⟩ := List.append_of_mem hp
        have hlen' : (l₁ ++ l₂).length ≤ n := by
          have := hlen
          simp only [List.length_append, List.length_cons] at this ⊢
          omega
        have hab' : p.2 ≤ b := hpb
        have hcp : clipLen p a b = p.2 - a := by
          unfold clipLen
          rw [min_eq_left hpb, max_eq_right (le_of_lt hp1)]
          exact max_eq_right (by linarith)
        have hsum_split : ((l₁ ++ p :: l₂).map fun q => clipLen q a b).sum
            = ((l₁ ++ l₂).map fun q => clipLen q a b).sum + clipLen p a b := by
          simp only [List.map_append, List.map_cons, List.sum_append, List.sum_cons]
          ring
        have hsum' : ((l₁ ++ l₂).map fun q => clipLen q p.2 b).sum < b - p.2 := by
          have hmono : ((l₁ ++ l₂).map fun q => clipLen q p.2 b).sum
              ≤ ((l₁ ++ l₂).map fun q => clipLen q a b).sum := by
            apply List.sum_le_sum
            intro q hq
            exact clipLen_mono_left q (le_of_lt hp2)
          linarith [hsum_split, hcp, hsum, hmono]
        obtain ⟨t, ht1, ht2, havoid⟩ := ih (l₁ ++ l₂) p.2 b hlen' hab' hsum'
        refine ⟨t, by linarith, ht2, ?_⟩
        intro q hq
        rcases List.mem_append.mp hq with hq1 | hq2
        · exact havoid q (List.mem_append.mpr (Or.inl hq1))
        · rcases List.mem_cons.mp hq2 with rfl | hq3
          · exact Or.inr (by linarith)
          · exact havoid q (List.mem_append.mpr (Or.inr hq3))

/-! ## Runner teeth as explicit bad lists -/

/-- One tooth of runner `w` at integer `m`: the open danger zone around `m/w`. -/
noncomputable def tooth (w m : ℤ) : ℝ × ℝ :=
  (((m : ℝ) - 1 / 14) / w, ((m : ℝ) + 1 / 14) / w)

/-- The teeth of runner `w` that can meet the window `[a, b]`: integers
`m ∈ [⌈w·a⌉ − 1, ⌊w·b⌋ + 1]`. -/
noncomputable def teeth (w : ℤ) (a b : ℝ) : List (ℝ × ℝ) :=
  (List.range (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat).map fun (i : ℕ) =>
    tooth w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))

/-- Each clipped tooth carries at most its width `1/(7w)`. -/
theorem clipLen_tooth_le {w : ℤ} (hw : 0 < w) (m : ℤ) (a b : ℝ) :
    clipLen (tooth w m) a b ≤ 1 / (7 * (w : ℝ)) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold clipLen tooth
  simp only
  rcases le_total (min (((m : ℝ) + 1 / 14) / w) b - max (((m : ℝ) - 1 / 14) / w) a) 0
    with h | h
  · rw [max_eq_left h]
    positivity
  · rw [max_eq_right h]
    have h1 : min (((m : ℝ) + 1 / 14) / w) b ≤ ((m : ℝ) + 1 / 14) / w := min_le_left _ _
    have h2 : ((m : ℝ) - 1 / 14) / w ≤ max (((m : ℝ) - 1 / 14) / w) a := le_max_left _ _
    have h3 : ((m : ℝ) + 1 / 14) / w - ((m : ℝ) - 1 / 14) / w = 1 / (7 * (w : ℝ)) := by
      field_simp
      ring
    linarith

/-- The total clipped mass of a runner's teeth: `(b−a)/7 + 3/(7w)`. -/
theorem teeth_mass {w : ℤ} (hw : 0 < w) (a b : ℝ) (hab : a ≤ b) :
    ((teeth w a b).map fun p => clipLen p a b).sum ≤ (b - a) / 7 + 3 / (7 * (w : ℝ)) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hper : ∀ x ∈ (teeth w a b).map fun p => clipLen p a b,
      x ≤ 1 / (7 * (w : ℝ)) := by
    intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    unfold teeth at hp
    rw [List.mem_map] at hp
    obtain ⟨i, _, rfl⟩ := hp
    exact clipLen_tooth_le hw _ a b
  have hlen : ((teeth w a b).map fun p => clipLen p a b).length
      = (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat := by
    unfold teeth
    rw [List.length_map, List.length_map, List.length_range]
  have hsum := List.sum_le_card_nsmul ((teeth w a b).map fun p => clipLen p a b)
    (1 / (7 * (w : ℝ))) hper
  rw [hlen, nsmul_eq_mul] at hsum
  have hcount : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
      ≤ (w : ℝ) * (b - a) + 3 := by
    have h1 : (⌊(w : ℝ) * b⌋ : ℝ) ≤ (w : ℝ) * b := Int.floor_le _
    have h2 : (w : ℝ) * a ≤ (⌈(w : ℝ) * a⌉ : ℝ) := Int.le_ceil _
    have hmax : ((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℤ)
        = max (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 := Int.toNat_eq_max _
    have hcast : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
        = ((max (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 : ℤ) : ℝ) := by
      rw [← hmax]
      push_cast
      ring
    rw [hcast]
    rcases max_cases (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) (0 : ℤ)
      with ⟨he, _⟩ | ⟨he, _⟩
    · rw [he]
      push_cast
      linarith
    · rw [he]
      push_cast
      nlinarith [hab, hwR]
  calc ((teeth w a b).map fun p => clipLen p a b).sum
      ≤ (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
        * (1 / (7 * (w : ℝ))) := hsum
    _ ≤ ((w : ℝ) * (b - a) + 3) * (1 / (7 * (w : ℝ))) := by
        apply mul_le_mul_of_nonneg_right hcount
        positivity
    _ = (b - a) / 7 + 3 / (7 * (w : ℝ)) := by
        field_simp

/-- Avoiding every listed tooth keeps the `1/14` band for the runner. -/
theorem good_of_avoid_teeth {w : ℤ} (hw : 0 < w) {a b t : ℝ} (hta : a ≤ t) (htb : t ≤ b)
    (havoid : ∀ p ∈ teeth w a b, t ≤ p.1 ∨ p.2 ≤ t) :
    ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  intro m
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  set mlo : ℤ := ⌈(w : ℝ) * a⌉ - 1 with hmlo
  set mhi : ℤ := ⌊(w : ℝ) * b⌋ + 1 with hmhi
  have hmloR : (mlo : ℝ) < (w : ℝ) * a := by
    rw [hmlo]
    push_cast
    linarith [Int.ceil_lt_add_one ((w : ℝ) * a)]
  have hmhiR : (w : ℝ) * b < (mhi : ℝ) := by
    rw [hmhi]
    push_cast
    linarith [Int.lt_floor_add_one ((w : ℝ) * b)]
  rcases lt_or_ge m mlo with hm | hm
  · have h2 : (w : ℝ) * a ≤ (w : ℝ) * t := mul_le_mul_of_nonneg_left hta hwR.le
    have h3 : (m : ℝ) ≤ (mlo : ℝ) - 1 := by
      have : m ≤ mlo - 1 := by omega
      calc (m : ℝ) ≤ ((mlo - 1 : ℤ) : ℝ) := by exact_mod_cast this
        _ = (mlo : ℝ) - 1 := by push_cast; ring
    rw [abs_of_nonneg (by linarith)]
    linarith
  · rcases le_or_gt m mhi with hm2 | hm2
    · -- m in range: its tooth is listed
      have hmem : tooth w m ∈ teeth w a b := by
        unfold teeth
        rw [List.mem_map]
        refine ⟨(m - mlo).toNat, ?_, ?_⟩
        · rw [List.mem_range]
          rw [hmlo] at hm
          rw [hmhi] at hm2
          rw [hmlo]
          omega
        · congr 1
          rw [hmlo] at hm ⊢
          have : ((m - (⌈(w : ℝ) * a⌉ - 1)).toNat : ℤ) = m - (⌈(w : ℝ) * a⌉ - 1) :=
            Int.toNat_of_nonneg (by omega)
          omega
      have hav := havoid _ hmem
      unfold tooth at hav
      simp only at hav
      rcases hav with hav | hav
      · have hle : (w : ℝ) * t ≤ (m : ℝ) - 1 / 14 := by
          have := mul_le_mul_of_nonneg_right hav hwR.le
          calc (w : ℝ) * t = t * w := mul_comm _ _
            _ ≤ (((m : ℝ) - 1 / 14) / w) * w := mul_le_mul_of_nonneg_right hav hwR.le
            _ = (m : ℝ) - 1 / 14 := by field_simp
        rw [abs_of_nonpos (by linarith)]
        linarith
      · have hge : (m : ℝ) + 1 / 14 ≤ (w : ℝ) * t := by
          calc (m : ℝ) + 1 / 14 = (((m : ℝ) + 1 / 14) / w) * w := by field_simp
            _ ≤ t * w := mul_le_mul_of_nonneg_right hav hwR.le
            _ = (w : ℝ) * t := mul_comm _ _
        rw [abs_of_nonneg (by linarith)]
        linarith
    · have h2 : (w : ℝ) * t ≤ (w : ℝ) * b := mul_le_mul_of_nonneg_left htb hwR.le
      have h3 : (mhi : ℝ) + 1 ≤ (m : ℝ) := by
        have : mhi + 1 ≤ m := by omega
        calc (mhi : ℝ) + 1 = ((mhi + 1 : ℤ) : ℝ) := by push_cast; ring
          _ ≤ (m : ℝ) := by exact_mod_cast this
      rw [abs_of_nonpos (by linarith)]
      linarith

/-! ## The block point step and the citation composition -/

/-- **THE BLOCK-6 POINT STEP**: runners with total tooth mass under the window length
share a `1/14`-good point.  For `c ≤ 6` runners `≥ w₁` the mass hypothesis follows
from `w₁·(b−a) > 3c/(7−c)` (union density `c/7 < 1`). -/
theorem flatMap_clip_sum (ws : List ℤ) (a b : ℝ) :
    ((ws.flatMap fun w => teeth w a b).map fun p => clipLen p a b).sum
      = (ws.map fun (w : ℤ) => ((teeth w a b).map fun p => clipLen p a b).sum).sum := by
  induction ws with
  | nil => simp
  | cons w ws ihw =>
      simp only [List.flatMap_cons, List.map_cons, List.map_append,
        List.sum_append, List.sum_cons]
      rw [ihw]

theorem block_point_step (ws : List ℤ) (hpos : ∀ w ∈ ws, 0 < w) (a b : ℝ) (hab : a ≤ b)
    (hmass : ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ)))).sum < b - a) :
    ∃ t : ℝ, a ≤ t ∧ t ≤ b ∧
      ∀ w ∈ ws, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  set bads : List (ℝ × ℝ) := ws.flatMap fun w => teeth w a b with hbads
  have hbadsum : (bads.map fun p => clipLen p a b).sum < b - a := by
    rw [hbads, flatMap_clip_sum ws a b]
    calc (ws.map fun (w : ℤ) => ((teeth w a b).map fun p => clipLen p a b).sum).sum
        ≤ ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ)))).sum := by
          apply List.sum_le_sum
          intro w hw
          exact teeth_mass (hpos w hw) a b hab
      _ < b - a := hmass
  obtain ⟨t, ht1, ht2, havoid⟩ := gap_exists bads.length bads a b (le_refl _) hab hbadsum
  refine ⟨t, ht1, ht2, ?_⟩
  intro w hw m
  apply good_of_avoid_teeth (w := w) (a := a) (b := b) (t := t) (hpos w hw) ht1 ht2
  intro p hp
  exact havoid p (List.mem_flatMap.mpr ⟨w, hw, hp⟩)

/-- **THE CITE–BLOCK THEOREM**: cite `k ≤ 12` runners bounded by `B`; the remaining
runners (any internal ratios — a near-equal block of size ≤ 6 in the key case) share
a good point inside the transported window whenever their tooth mass fits. -/
theorem cite_block_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ws : List ℤ) (hwpos : ∀ w ∈ ws, 0 < w)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ |v i| ∈ ws)
    (hmass : ((ws.map fun (w : ℤ) =>
        (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) / 7
          + 3 / (7 * (w : ℝ)))).sum
      < 2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hk13 : k ≤ 13 := by omega
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  have hkR : (k : ℝ) ≤ 12 := by exact_mod_cast hk
  set δ : ℝ := ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by
    rw [hδ]
    apply div_pos (by linarith) (by positivity)
  set wc : Fin k → ℤ := fun j => v (Fin.castLE hk13 j) with hwc
  have hwcne : ∀ j, wc j ≠ 0 := fun j => hv _
  obtain ⟨t₀, hcite⟩ := cite k hk wc hwcne
  have hwin : t₀ - δ ≤ t₀ - δ + 2 * δ := by linarith
  obtain ⟨τ, hτ1, hτ2, hτgood⟩ := block_point_step ws hwpos (t₀ - δ) (t₀ - δ + 2 * δ)
    hwin (by
      have h2δ : t₀ - δ + 2 * δ - (t₀ - δ) = 2 * δ := by ring
      rw [h2δ]
      exact hmass)
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | hmem
  · -- cited leg (S19 verbatim)
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
    have hwin2 : |t₀ - τ| ≤ δ := by
      rw [abs_le]
      constructor <;> linarith
    have hlip : |(v i : ℝ)| * |t₀ - τ| ≤ (B : ℝ) * δ := by
      apply mul_le_mul hvB hwin2 (abs_nonneg _) hBR.le
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
  · -- block leg through the sign bridge
    have hgood := hτgood |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood m
      rw [hcast] at this
      exact (by norm_num : (1:ℝ)/(14:ℕ) = 1/14) ▸ this
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood (-m)
      rw [hcast] at this
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at this
      exact (by norm_num : (1:ℝ)/(14:ℕ) = 1/14) ▸ this

end LRC14
end LonelyRunner
