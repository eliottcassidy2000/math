/- FragmentationLemma.lean -- mac-mini-2026-07-16-S123 (statements),
   death-star-2026-07-16-S30 (periodicity lemma + independent catch of the draft plan's
   flaw), klein-2026-07-16-S316 (complete proofs; concurrent-edit merge).

   THM-883's core counting inequality, the first Lean target of the covering program.
   Statement: an arc-grid of modulus w >= 1 and radius lam meets an interval I
   of length L in measure at most (L*w + 1) * (2*lam/w).

   NOTE (both sessions independently): the draft's arc-COUNTING route is subtly flawed --
   the number of arcs meeting I can reach floor(L*w) + 2 (two partial edge arcs), which
   overshoots the claimed constant.  The correct route is windowing:

   WINDOW LEMMA -- badArcs meets any half-open window of length 1/w in measure at most
   2*lam/w (split the window at the unique grid-cell boundary inside it; each piece sees
   only one arc, and the two clipped arcs are complementary translates of a single arc,
   so their lengths sum to at most one full arc).  Tile [x, x+L] by floor(L*w) + 1
   windows.  The case lam > 1/2 is trivial (the bound exceeds L).  No sorries. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

open MeasureTheory

namespace LRC14

/-- The bad arc set of modulus `w` at radius `lam`, within the unit circle lifted to ℝ. -/
def badArcs (w : ℕ) (lam : ℝ) : Set ℝ :=
  ⋃ a : ℤ, Set.Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)

/-- `badArcs` is `(1/w)`-periodic: translation by `1/w` permutes the arcs.
    (death-star-S30; API lemma, not needed by the proofs below.) -/
theorem badArcs_periodic (w : ℕ) (hw : 1 ≤ w) (lam : ℝ) :
    (fun y => y + 1 / (w : ℝ)) ⁻¹' badArcs w lam = badArcs w lam := by
  ext y
  simp only [badArcs, Set.mem_preimage, Set.mem_iUnion, Set.mem_Ioo]
  constructor
  · rintro ⟨a, h1, h2⟩
    refine ⟨a - 1, ?_, ?_⟩
    · have : ((a - 1 : ℤ) : ℝ) = (a : ℝ) - 1 := by push_cast; ring
      rw [this]
      rw [sub_div]
      linarith [h1]
    · have : ((a - 1 : ℤ) : ℝ) = (a : ℝ) - 1 := by push_cast; ring
      rw [this]
      rw [sub_div]
      linarith [h2]
  · rintro ⟨a, h1, h2⟩
    refine ⟨a + 1, ?_, ?_⟩
    · have : ((a + 1 : ℤ) : ℝ) = (a : ℝ) + 1 := by push_cast; ring
      rw [this, add_div]
      linarith [h1]
    · have : ((a + 1 : ℤ) : ℝ) = (a : ℝ) + 1 := by push_cast; ring
      rw [this, add_div]
      linarith [h2]

/-- The window lemma: the arc-grid meets any half-open window of length `1/w`
    in measure at most `2*lam/w` (for `lam ≤ 1/2`). -/
lemma window_bound (w : ℕ) (hw : 1 ≤ w) (lam : ℝ) (hlam : 0 < lam)
    (hlam2 : lam ≤ 1/2) (c : ℝ) :
    volume (badArcs w lam ∩ Set.Ico c (c + 1 / w))
      ≤ ENNReal.ofReal (2 * lam / w) := by
  have hw0 : (0:ℝ) < w := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hw
  -- the unique cell around c and the boundary point inside the window
  set a0 : ℤ := ⌊c * w + 1/2⌋ with ha0
  have h1 : (a0 : ℝ) ≤ c * w + 1/2 := Int.floor_le _
  have h2 : c * w + 1/2 < (a0 : ℝ) + 1 := Int.lt_floor_add_one _
  set bd : ℝ := ((a0 : ℝ) + 1/2) / w with hbd
  have hcbd : c ≤ bd := by
    rw [hbd, le_div_iff₀ hw0]; linarith
  have hbdc : bd ≤ c + 1 / w := by
    rw [hbd, div_le_iff₀ hw0]
    have hne : (w:ℝ) ≠ 0 := ne_of_gt hw0
    have : (c + 1 / w) * w = c * w + 1 := by field_simp
    rw [this]; linarith
  -- names for the central arc endpoints
  set bot : ℝ := (a0 : ℝ) / w - lam / w with hbot
  set top : ℝ := (a0 : ℝ) / w + lam / w with htop
  -- piece 1: on [c, bd) only arc a0 appears
  have P1 : badArcs w lam ∩ Set.Ico c bd ⊆ Set.Icc (max bot c) top := by
    rintro t ⟨htb, htc, htbd⟩
    obtain ⟨a, hta⟩ := Set.mem_iUnion.mp htb
    obtain ⟨hta1, hta2⟩ := hta
    have m1 : (a : ℝ) - lam < t * w := by
      have h' : ((a : ℝ) - lam) / w < t := by rw [sub_div]; exact hta1
      exact (div_lt_iff₀ hw0).mp h'
    have m2 : t * w < (a : ℝ) + lam := by
      have h' : t < ((a : ℝ) + lam) / w := by rw [add_div]; exact hta2
      exact (lt_div_iff₀ hw0).mp h'
    have mc : c * w ≤ t * w := mul_le_mul_of_nonneg_right htc (le_of_lt hw0)
    have mb : t * w < (a0 : ℝ) + 1/2 := by
      have := htbd
      rw [hbd] at this
      exact (lt_div_iff₀ hw0).mp this
    have haa0 : a = a0 := by
      have e1 : (a : ℝ) < (a0 : ℝ) + 1 := by linarith
      have e2 : (a0 : ℝ) < (a : ℝ) + 1 := by linarith
      have e1' : a < a0 + 1 := by exact_mod_cast e1
      have e2' : a0 < a + 1 := by exact_mod_cast e2
      omega
    subst haa0
    refine ⟨max_le (le_of_lt hta1) htc, le_of_lt hta2⟩
  -- piece 2: on [bd, c + 1/w) only arc a0+1 appears
  have P2 : badArcs w lam ∩ Set.Ico bd (c + 1 / w) ⊆
      Set.Icc (((a0 : ℝ) + 1) / w - lam / w)
              (min (((a0 : ℝ) + 1) / w + lam / w) (c + 1 / w)) := by
    rintro t ⟨htb, htbd, htc⟩
    obtain ⟨a, hta⟩ := Set.mem_iUnion.mp htb
    obtain ⟨hta1, hta2⟩ := hta
    have m1 : (a : ℝ) - lam < t * w := by
      have h' : ((a : ℝ) - lam) / w < t := by rw [sub_div]; exact hta1
      exact (div_lt_iff₀ hw0).mp h'
    have m2 : t * w < (a : ℝ) + lam := by
      have h' : t < ((a : ℝ) + lam) / w := by rw [add_div]; exact hta2
      exact (lt_div_iff₀ hw0).mp h'
    have mb : (a0 : ℝ) + 1/2 ≤ t * w := by
      have := htbd
      rw [hbd] at this
      exact (div_le_iff₀ hw0).mp this
    have mc : t * w < c * w + 1 := by
      have h' : t < (c * w + 1) / w := by
        have hne : (w:ℝ) ≠ 0 := ne_of_gt hw0
        have e : (c * w + 1) / w = c + 1 / w := by field_simp
        rw [e]; exact htc
      exact (lt_div_iff₀ hw0).mp h'
    have haa1 : a = a0 + 1 := by
      have e1 : (a : ℝ) < (a0 : ℝ) + 2 := by linarith
      have e2 : (a0 : ℝ) < (a : ℝ) := by linarith
      have e1' : a < a0 + 2 := by exact_mod_cast e1
      have e2' : a0 < a := by exact_mod_cast e2
      omega
    subst haa1
    have hcast : ((a0 + 1 : ℤ) : ℝ) = (a0 : ℝ) + 1 := by push_cast; ring
    rw [hcast] at hta1 hta2
    exact ⟨le_of_lt hta1, le_min (le_of_lt hta2) (le_of_lt htc)⟩
  -- split the window
  have hsplit : Set.Ico c (c + 1 / w) = Set.Ico c bd ∪ Set.Ico bd (c + 1 / w) :=
    (Set.Ico_union_Ico_eq_Ico hcbd hbdc).symm
  have hsub : badArcs w lam ∩ Set.Ico c (c + 1 / w) =
      (badArcs w lam ∩ Set.Ico c bd) ∪ (badArcs w lam ∩ Set.Ico bd (c + 1 / w)) := by
    rw [hsplit, Set.inter_union_distrib_left]
  -- measures of the two pieces
  set B1 : ℝ := top - max bot c with hB1def
  set B2 : ℝ := min (((a0 : ℝ) + 1) / w + lam / w) (c + 1 / w)
                  - (((a0 : ℝ) + 1) / w - lam / w) with hB2def
  have M1 : volume (badArcs w lam ∩ Set.Ico c bd) ≤ ENNReal.ofReal B1 := by
    calc volume (badArcs w lam ∩ Set.Ico c bd)
        ≤ volume (Set.Icc (max bot c) top) := measure_mono P1
      _ = ENNReal.ofReal B1 := by rw [Real.volume_Icc]
  have M2 : volume (badArcs w lam ∩ Set.Ico bd (c + 1 / w)) ≤ ENNReal.ofReal B2 := by
    calc volume (badArcs w lam ∩ Set.Ico bd (c + 1 / w))
        ≤ volume (Set.Icc (((a0 : ℝ) + 1) / w - lam / w)
            (min (((a0 : ℝ) + 1) / w + lam / w) (c + 1 / w))) := measure_mono P2
      _ = ENNReal.ofReal B2 := by rw [Real.volume_Icc]
  -- the real inequalities
  have hgap : top - bot = 2 * lam / w := by rw [htop, hbot]; ring
  have hB1le : B1 ≤ 2 * lam / w := by
    have := le_max_left bot c
    rw [hB1def]; linarith
  have hB2le : B2 ≤ 2 * lam / w := by
    have hm : min (((a0 : ℝ) + 1) / w + lam / w) (c + 1 / w)
        ≤ ((a0 : ℝ) + 1) / w + lam / w := min_le_left _ _
    have e : (((a0 : ℝ) + 1) / w + lam / w) - (((a0 : ℝ) + 1) / w - lam / w)
        = 2 * lam / w := by ring
    rw [hB2def]; linarith
  have hsum : B1 + B2 ≤ 2 * lam / w := by
    -- B2 = min top c - bot  after shifting by 1/w
    have e2 : B2 = min top c - bot := by
      rw [hB2def, htop, hbot]
      have ea : ((a0 : ℝ) + 1) / w + lam / w = ((a0 : ℝ) / w + lam / w) + 1 / w := by
        push_cast; ring
      have eb : ((a0 : ℝ) + 1) / w - lam / w = ((a0 : ℝ) / w - lam / w) + 1 / w := by
        push_cast; ring
      rw [ea, eb]
      rw [min_add_add_right]
      ring
    have hminmax : min top c ≤ max bot c := by
      rcases le_total bot c with h | h
      · have : max bot c = c := max_eq_right h
        rw [this]; exact min_le_right _ _
      · have : max bot c = bot := max_eq_left h
        rw [this]
        calc min top c ≤ c := min_le_right _ _
          _ ≤ bot := h
    rw [hB1def, e2]; linarith
  -- assemble in ENNReal
  calc volume (badArcs w lam ∩ Set.Ico c (c + 1 / w))
      = volume ((badArcs w lam ∩ Set.Ico c bd) ∪
                (badArcs w lam ∩ Set.Ico bd (c + 1 / w))) := by rw [hsub]
    _ ≤ volume (badArcs w lam ∩ Set.Ico c bd)
        + volume (badArcs w lam ∩ Set.Ico bd (c + 1 / w)) := measure_union_le _ _
    _ ≤ ENNReal.ofReal B1 + ENNReal.ofReal B2 := add_le_add M1 M2
    _ ≤ ENNReal.ofReal (2 * lam / w) := by
        rcases le_total B1 0 with hb1 | hb1
        · rw [ENNReal.ofReal_of_nonpos hb1, zero_add]
          exact ENNReal.ofReal_le_ofReal hB2le
        · rcases le_total B2 0 with hb2 | hb2
          · rw [ENNReal.ofReal_of_nonpos hb2, add_zero]
            exact ENNReal.ofReal_le_ofReal hB1le
          · rw [← ENNReal.ofReal_add hb1 hb2]
            exact ENNReal.ofReal_le_ofReal hsum

/-- THM-883 Lemma 1 (fragmentation): the arc-grid meets an interval `I = [x, x+L]`
    in measure at most `(L * w + 1) * (2 * lam / w)`. -/
theorem fragmentation (w : ℕ) (hw : 1 ≤ w) (lam L x : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    volume (badArcs w lam ∩ Set.Icc x (x + L))
      ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := by
  have hw0 : (0:ℝ) < w := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hw
  have hy0 : (0:ℝ) ≤ 2 * lam / w := by positivity
  rcases le_total lam (1/2) with hsmall | hbig
  · -- tile [x, x+L] by N = ⌊L*w⌋ + 1 windows of length 1/w
    set N : ℕ := Nat.floor (L * w) + 1 with hN
    have hLw0 : 0 ≤ L * w := mul_nonneg hL (le_of_lt hw0)
    have hcover : Set.Icc x (x + L) ⊆
        ⋃ k ∈ Finset.range N, Set.Ico (x + k / w) (x + k / w + 1 / w) := by
      intro t ht
      obtain ⟨ht1, ht2⟩ := ht
      have h0 : 0 ≤ (t - x) * w := mul_nonneg (by linarith) (le_of_lt hw0)
      set k : ℕ := Nat.floor ((t - x) * w) with hk
      have hkN : k < N := by
        have hle : (t - x) * w ≤ L * w :=
          mul_le_mul_of_nonneg_right (by linarith) (le_of_lt hw0)
        have := Nat.floor_le_floor hle
        omega
      refine Set.mem_biUnion (Finset.mem_range.mpr hkN) ⟨?_, ?_⟩
      · -- x + k/w ≤ t
        have hkle : (k : ℝ) ≤ (t - x) * w := Nat.floor_le h0
        have : (k : ℝ) / w ≤ t - x := by
          rw [div_le_iff₀ hw0]; exact hkle
        linarith
      · -- t < x + k/w + 1/w
        have hklt : (t - x) * w < (k : ℝ) + 1 := Nat.lt_floor_add_one _
        have : t - x < ((k : ℝ) + 1) / w := by
          rw [lt_div_iff₀ hw0]; exact hklt
        have e : ((k : ℝ) + 1) / w = (k : ℝ) / w + 1 / w := by ring
        rw [e] at this
        linarith
    calc volume (badArcs w lam ∩ Set.Icc x (x + L))
        ≤ volume (⋃ k ∈ Finset.range N,
            badArcs w lam ∩ Set.Ico (x + k / w) (x + k / w + 1 / w)) := by
          apply measure_mono
          rintro t ⟨htb, hti⟩
          obtain ⟨s, hs, hts⟩ := Set.mem_iUnion₂.mp (hcover hti)
          exact Set.mem_biUnion hs ⟨htb, hts⟩
      _ ≤ ∑ k ∈ Finset.range N,
            volume (badArcs w lam ∩ Set.Ico (x + k / w) (x + k / w + 1 / w)) :=
          measure_biUnion_finset_le _ _
      _ ≤ ∑ _k ∈ Finset.range N, ENNReal.ofReal (2 * lam / w) :=
          Finset.sum_le_sum fun k _ => window_bound w hw lam hlam hsmall (x + k / w)
      _ = (N : ENNReal) * ENNReal.ofReal (2 * lam / w) := by
          rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
      _ ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := by
          rw [← ENNReal.ofReal_natCast N, ← ENNReal.ofReal_mul (by positivity)]
          apply ENNReal.ofReal_le_ofReal
          apply mul_le_mul_of_nonneg_right _ hy0
          have hfl : (Nat.floor (L * w) : ℝ) ≤ L * w := Nat.floor_le hLw0
          have : (N : ℝ) = (Nat.floor (L * w) : ℝ) + 1 := by
            rw [hN]; push_cast; ring
          linarith
  · -- lam ≥ 1/2: the bound already exceeds the whole interval length
    calc volume (badArcs w lam ∩ Set.Icc x (x + L))
        ≤ volume (Set.Icc x (x + L)) := measure_mono Set.inter_subset_right
      _ = ENNReal.ofReal L := by rw [Real.volume_Icc]; try ring_nf
      _ ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := by
          apply ENNReal.ofReal_le_ofReal
          have hne : (w:ℝ) ≠ 0 := ne_of_gt hw0
          have e : (L * w + 1) * (2 * lam / w) = 2 * lam * L + 2 * lam / w := by
            field_simp
            try ring
          rw [e]
          have h1 : L ≤ 2 * lam * L := by nlinarith
          have h2 : (0:ℝ) < 2 * lam / w := by positivity
          linarith

/-- Corollary (killer budget, THM-883 Lemma 2 shape): if a family of moduli
    `w₁ ≤ … ≤ w_j` covers `[x, x+L]` by their arc-grids, then
    `L * (1 - 2*j*lam) ≤ 2*lam * Σ 1/wᵢ`. -/
theorem killer_budget (j : ℕ) (ws : Fin j → ℕ) (hws : ∀ i, 1 ≤ ws i)
    (lam L x : ℝ) (hlam : 0 < lam) (hL : 0 ≤ L)
    (hcover : Set.Icc x (x + L) ⊆ ⋃ i, badArcs (ws i) lam) :
    L * (1 - 2 * j * lam) ≤ 2 * lam * ∑ i, (1 : ℝ) / ws i := by
  have hwpos : ∀ i : Fin j, (0:ℝ) < ws i := fun i => by
    exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one (hws i)
  -- measure chain: L = |I| ≤ Σ |bad_i ∩ I| ≤ Σ (L wᵢ + 1)(2 lam / wᵢ)
  have hchain : ENNReal.ofReal L
      ≤ ∑ i : Fin j, ENNReal.ofReal ((L * ws i + 1) * (2 * lam / ws i)) := by
    calc ENNReal.ofReal L
        = volume (Set.Icc x (x + L)) := by rw [Real.volume_Icc]; try ring_nf
      _ ≤ volume (⋃ i : Fin j, badArcs (ws i) lam ∩ Set.Icc x (x + L)) := by
          apply measure_mono
          intro t ht
          have := hcover ht
          rw [Set.mem_iUnion] at this
          obtain ⟨i, hi⟩ := this
          exact Set.mem_iUnion.mpr ⟨i, hi, ht⟩
      _ ≤ ∑' i : Fin j, volume (badArcs (ws i) lam ∩ Set.Icc x (x + L)) :=
          measure_iUnion_le _
      _ = ∑ i : Fin j, volume (badArcs (ws i) lam ∩ Set.Icc x (x + L)) :=
          tsum_fintype _
      _ ≤ ∑ i : Fin j, ENNReal.ofReal ((L * ws i + 1) * (2 * lam / ws i)) :=
          Finset.sum_le_sum fun i _ =>
            fragmentation (ws i) (hws i) lam L x hlam hL
  -- convert to a real inequality
  have hterm : ∀ i : Fin j,
      (L * ws i + 1) * (2 * lam / ws i) = 2 * lam * L + 2 * lam / ws i := by
    intro i
    have hne : ((ws i : ℝ)) ≠ 0 := ne_of_gt (hwpos i)
    field_simp
    try ring
  have hnonneg : ∀ i ∈ Finset.univ, (0:ℝ) ≤ (L * ws i + 1) * (2 * lam / ws i) := by
    intro i _
    rw [hterm i]
    have h1 : (0:ℝ) ≤ 2 * lam * L := by positivity
    have h2 : (0:ℝ) < 2 * lam / ws i := by
      have := hwpos i; positivity
    linarith
  have hsum_eq : ∑ i : Fin j, ENNReal.ofReal ((L * ws i + 1) * (2 * lam / ws i))
      = ENNReal.ofReal (∑ i : Fin j, (L * ws i + 1) * (2 * lam / ws i)) :=
    (ENNReal.ofReal_sum_of_nonneg hnonneg).symm
  rw [hsum_eq] at hchain
  have hsumnn : (0:ℝ) ≤ ∑ i : Fin j, (L * ws i + 1) * (2 * lam / ws i) :=
    Finset.sum_nonneg hnonneg
  have hreal : L ≤ ∑ i : Fin j, (L * ws i + 1) * (2 * lam / ws i) :=
    (ENNReal.ofReal_le_ofReal_iff hsumnn).mp hchain
  -- expand the sum
  have hexp : ∑ i : Fin j, (L * ws i + 1) * (2 * lam / ws i)
      = (j : ℝ) * (2 * lam * L) + 2 * lam * ∑ i : Fin j, (1:ℝ) / ws i := by
    calc ∑ i : Fin j, (L * ws i + 1) * (2 * lam / ws i)
        = ∑ i : Fin j, (2 * lam * L + 2 * lam / ws i) :=
          Finset.sum_congr rfl fun i _ => hterm i
      _ = ∑ _i : Fin j, 2 * lam * L + ∑ i : Fin j, 2 * lam / ws i :=
          Finset.sum_add_distrib
      _ = (j : ℝ) * (2 * lam * L) + 2 * lam * ∑ i : Fin j, (1:ℝ) / ws i := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
          congr 1
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun i _ => by
            rw [mul_one_div]
  rw [hexp] at hreal
  nlinarith [hreal]

end LRC14
