/-
  TournamentH7.LRCOneWindowPeel — THE ONE-WINDOW PEEL LEMMA + THE PRIMITIVITY SPLIT
  (klein-2026-07-05-S132, HYP-4095).

  The uniform threshold side of hcomp's loose-base case (mac-mini HYP-4094 alignment
  bands are the complementary per-base side):

  * `good_point_in_long_interval` — any interval longer than one bad arc `1/(7u)`
    contains a point at circle-distance ≥ 1/14 from every integer multiple; the
    witness is explicit: the left endpoint, or `(m₀ + 1/14)/u` past the blocking arc.
  * `lonely_of_window_margin` — if the 12 base runners have margin β > 1/14 at one
    point t* and the killer satisfies `B < 14(β − 1/14)·|v*|` (B = base speed bound),
    the family is lonely: the Lipschitz window `[t* ± (β−1/14)/B]` keeps the base
    ≥ 1/14 and outlasts one killer bad arc.  At β = 1/13 the threshold is EXACTLY
    the dominant 13× line; at β = 2/25 it is (25/3)·B; at β ≥ 1/7 every killer passes.
    Verified constructively on 400 random instances (lrc14_one_window_peel_klein_S132).
  * `hcomp_of_primitive` — klein-S131's split: the compressed leaf reduces to
    PRIMITIVE families; a non-primitive family scales down by its gcd, and the
    quotient is either covering (primitive case) or sieved at 1/14.

  Kernel-pure; no native_decide.
-/
import TournamentH7.LRC14CertRoute
import TournamentH7.LRC14WindowWiring

namespace LonelyRunner
namespace OneWindowPeel

open LRC14

/-- **A long interval contains a 1/14-good point**: if `ℓ > 1/(7u)` (one bad arc),
some `t ∈ [c, c+ℓ]` has `|u·t − m| ≥ 1/14` for every integer `m`.  Constructive:
either the left endpoint is good, or the blocking arc's right edge `(m₀ + 1/14)/u`
is inside the interval. -/
theorem good_point_in_long_interval (u c ℓ : ℝ) (hu : 0 < u) (hℓ : 1 / (7 * u) < ℓ) :
    ∃ t : ℝ, c ≤ t ∧ t ≤ c + ℓ ∧ ∀ m : ℤ, 1/14 ≤ |u * t - m| := by
  have hℓpos : 0 < ℓ := lt_trans (by positivity) hℓ
  by_cases h0 : ∀ m : ℤ, 1/14 ≤ |u * c - m|
  · exact ⟨c, le_refl c, by linarith, h0⟩
  · push Not at h0
    obtain ⟨m₀, hm₀⟩ := h0
    have h1 : u * c - m₀ < 1/14 := (abs_lt.mp hm₀).2
    have h2 : -(1/14) < u * c - m₀ := (abs_lt.mp hm₀).1
    have h3 : 1 / 7 < ℓ * u := by
      rw [div_lt_iff₀ (by positivity : (0:ℝ) < 7 * u)] at hℓ
      linarith
    refine ⟨((m₀ : ℝ) + 1/14) / u, ?_, ?_, ?_⟩
    · rw [le_div_iff₀ hu]
      nlinarith
    · rw [div_le_iff₀ hu]
      nlinarith
    · intro m
      have hval : u * (((m₀ : ℝ) + 1/14) / u) = (m₀ : ℝ) + 1/14 := by
        field_simp
      rw [hval]
      rcases eq_or_ne m m₀ with rfl | hne
      · rw [show (m : ℝ) + 1/14 - m = 1/14 by ring]
        norm_num
      · have h4 : (1 : ℝ) ≤ |(m₀ : ℝ) - m| := by
          have hZ : (1 : ℤ) ≤ |m₀ - m| :=
            Int.one_le_abs (sub_ne_zero.mpr (fun h => hne h.symm))
          calc (1 : ℝ) = ((1 : ℤ) : ℝ) := by norm_num
            _ ≤ ((|m₀ - m| : ℤ) : ℝ) := by exact_mod_cast hZ
            _ = |((m₀ - m : ℤ) : ℝ)| := by rw [Int.cast_abs]
            _ = |(m₀ : ℝ) - m| := by push_cast; ring_nf
        have h5 : |(m₀ : ℝ) - m| ≤ |(m₀ : ℝ) + 1/14 - m| + 1/14 := by
          calc |(m₀ : ℝ) - m| = |((m₀ : ℝ) + 1/14 - m) + (-(1/14))| := by ring_nf
            _ ≤ |(m₀ : ℝ) + 1/14 - m| + |(-(1/14) : ℝ)| := abs_add_le _ _
            _ = |(m₀ : ℝ) + 1/14 - m| + 1/14 := by norm_num
        linarith

/-- **THE ONE-WINDOW PEEL LEMMA**: base margin β at a single point beats any killer
above the threshold `B/(14(β − 1/14))`.  At β = 1/13 (the LRC(12) citation floor)
the threshold is EXACTLY 13·B — the dominant/compressed line; every improvement of
the base margin buys down the killer threshold continuously. -/
theorem lonely_of_window_margin (v : Fin 13 → ℤ) (istar : Fin 13) (tstar β B : ℝ)
    (hβ : 1/14 < β)
    (hbase : ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)
    (hB : ∀ i, i ≠ istar → |(v i : ℝ)| ≤ B) (hBpos : 0 < B)
    (hkne : v istar ≠ 0)
    (hkill : B < 14 * (β - 1/14) * |(v istar : ℝ)|) :
    ∃ t : ℝ, Lonely 14 v t := by
  set u : ℝ := |(v istar : ℝ)| with hu
  have hupos : 0 < u := abs_pos.mpr (by exact_mod_cast hkne)
  set δ : ℝ := (β - 1/14) / B with hδ
  have hδpos : 0 < δ := div_pos (by linarith) hBpos
  have hlen : 1 / (7 * u) < 2 * δ := by
    rw [div_lt_iff₀ (by positivity : (0:ℝ) < 7 * u), hδ]
    rw [show 2 * ((β - 1/14) / B) * (7 * u) = 14 * (β - 1/14) * u / B by ring]
    rw [lt_div_iff₀ hBpos]
    linarith
  obtain ⟨t, ht1, ht2, htk⟩ :=
    good_point_in_long_interval u (tstar - δ) (2 * δ) hupos hlen
  have htw : |t - tstar| ≤ δ := by
    rw [abs_le]
    constructor <;> linarith
  refine ⟨t, fun i m => ?_⟩
  show (1 : ℝ) / (14 : ℕ) ≤ |(v i : ℝ) * t - m|
  push_cast
  by_cases hi : i = istar
  · subst hi
    rcases abs_cases ((v i : ℝ)) with ⟨habs, -⟩ | ⟨habs, -⟩
    · have := htk m
      rw [hu, habs] at this
      linarith
    · have := htk (-m)
      rw [hu, habs] at this
      have heq : |-(v i : ℝ) * t - (-m : ℤ)| = |(v i : ℝ) * t - m| := by
        push_cast
        rw [show -(v i : ℝ) * t - (-(m : ℝ)) = -((v i : ℝ) * t - m) by ring, abs_neg]
      rw [heq] at this
      linarith
  · have hbi := hbase i hi m
    have hBi := hB i hi
    have key : |(v i : ℝ) * tstar - m| ≤ |(v i : ℝ) * t - m| + |(v i : ℝ)| * |tstar - t| := by
      calc |(v i : ℝ) * tstar - m|
          = |((v i : ℝ) * t - m) + (v i : ℝ) * (tstar - t)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * t - m| + |(v i : ℝ) * (tstar - t)| := abs_add_le _ _
        _ = |(v i : ℝ) * t - m| + |(v i : ℝ)| * |tstar - t| := by rw [abs_mul]
    have habs' : |(v i : ℝ)| * |tstar - t| ≤ B * δ := by
      apply mul_le_mul hBi ?_ (abs_nonneg _) (le_of_lt hBpos)
      rwa [abs_sub_comm]
    have hBδ : B * δ = β - 1/14 := by
      rw [hδ]
      field_simp
    linarith

/-- **THE PRIMITIVITY SPLIT** (klein-S131): the compressed leaf `hcomp` reduces to
PRIMITIVE compressed covering families.  A non-primitive family is `g · w` with `w`
primitive; `w` is still compressed; if `w` is covering it is handled by the primitive
hypothesis, and if not it is lonely by the denominator sieve — in both cases loneliness
transports back through the scale `g`. -/
theorem hcomp_of_primitive
    (hprim : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∃ t : ℝ, Lonely 14 v t) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t := by
  intro v hv hcov hcomp
  set g : ℕ := tupleGcd v with hg
  have hgpos : 0 < g := tupleGcd_pos v hv
  have hgZ : (0 : ℤ) < (g : ℤ) := by exact_mod_cast hgpos
  set w : Fin 13 → ℤ := primPart v with hw
  have hvw : ∀ i, v i = (g : ℤ) * w i := fun i => primPart_mul v i
  have hwne : ∀ i, w i ≠ 0 := by
    intro i hzero
    exact hv i (by rw [hvw i, hzero, mul_zero])
  have hvfun : v = fun i => (g : ℤ) * w i := funext hvw
  have htrans : (∃ t : ℝ, Lonely 14 w t) → ∃ t : ℝ, Lonely 14 v t := by
    intro h
    rw [hvfun]
    exact lonely_exists_of_scale 14 w (g : ℤ) (ne_of_gt hgZ) h
  have hwcomp : ∀ i, ∃ j, j ≠ i ∧ |w i| ≤ 13 * |w j| := by
    intro i
    obtain ⟨j, hj, hij⟩ := hcomp i
    refine ⟨j, hj, ?_⟩
    have h1 : |v i| = (g : ℤ) * |w i| := by
      rw [hvw i, abs_mul, abs_of_nonneg (le_of_lt hgZ)]
    have h2 : |v j| = (g : ℤ) * |w j| := by
      rw [hvw j, abs_mul, abs_of_nonneg (le_of_lt hgZ)]
    rw [h1, h2] at hij
    have h3 : (g : ℤ) * |w i| ≤ (g : ℤ) * (13 * |w j|) := by linarith
    exact le_of_mul_le_mul_left h3 hgZ
  by_cases hcw : CoveringFamily w
  · exact htrans (hprim w hwne hcw hwcomp (primPart_gcd_eq_one v hv))
  · unfold CoveringFamily at hcw
    push Not at hcw
    obtain ⟨q, hq2, hq14, hdiv⟩ := hcw
    exact htrans ⟨(1 : ℝ) / q, sieve_one_div 14 q w hq14 (by omega) hdiv⟩

#print axioms lonely_of_window_margin
#print axioms hcomp_of_primitive

end OneWindowPeel
end LonelyRunner
