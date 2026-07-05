/-
  TournamentH7.LRCHarmonicGate — THE RATIONAL-POINT MARGIN CERTIFICATE and the
  RATIO-GATE REDUCTION of TightLooseDichotomy (kind-pasteur-2026-07-05-S2, HYP-4102).

  THE ATOM (`rational_point_margin`): for integers `s > 0`, `k`, `μ` and any speed `v`
  with `μ ≤ (v·k) % s ≤ s − μ`, the point `t = k/s` has `|v·t − m| ≥ μ/s` for EVERY
  integer `m`.  Pure integer arithmetic: `|v·k − m·s| ≥ min(r, s−r) ≥ μ` where
  `r = (v·k) % s`.  The hypothesis is DECIDABLE per family — every loose-branch
  witness of the dichotomy is an instance (M is attained at rationals), so this is
  the universal certificate language for the loose side, the exact complement of
  mac-mini's killer-side band certificates (THM-619/620, IntervalEscape).

  Instances formalized here:
  * `k = 1`, `s = a+b`, `μ = a` — the reciprocal band point (kps-S1
    `band_margin_reciprocal` re-derived through the atom, hypothesis-decidable form);
  * `k = 2`, `s = a+b`, `μ = 2a` — the SECOND HARMONIC with the middle hole
    `2|w| ∉ (b−a, 3a+b)`: covers the rigidity spectrum's second value
    `{1,…,11,24}` at `t = 2/25`, margin exactly `2/25`;
  * `μ = 1` — the sieve shape (`s ∤ v` ⟹ margin `1/s` at `1/s`) is the degenerate
    instance, unifying `sieve_one_div` with the band lemmas.

  THE REDUCTION (`tightLooseDichotomyAt_of_spread`): the dichotomy's loose branch
  holds OUTRIGHT at margin `β` whenever the base band `[a,b]` is narrow
  (`β·(a+b) ≤ a`), via the `k=1` instance.  Hence `TightLooseDichotomyAt β` only
  needs proving on SPREAD bases (`a < β·(a+b)`, i.e. ratio `> (1−β)/β`):
  at `β = 2/25` that is ratio `> 23/2`; the AP `c·{1,…,12}` (ratio 12) stays on
  the tight side, `{1,…,11,24}·c` (ratio 24) on the second-harmonic side.
  At `β = 1/13` the tight branch itself forces ratio `≤ 12` ⟹ the gate SUBSUMES
  the tight branch (`dichotomyAt_13_of_spread_loose`): the 1/13-dichotomy is
  purely a spread-base margin statement.

  VERIFIED (lrc_harmonic_gate_kps_S2.py): atom 20k random instances 0 violations;
  gate never lies on all 1820 primitive 12-subsets of [1,16]; every one of those
  1820 subsets except {1,…,12} itself is atom-certified loose at margin 2/25
  (hard survivors: 0).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHcompSurface

namespace LonelyRunner
namespace HarmonicGate

open LRC14 HcompSurface

/-! ## The emod shift helper -/

/-- Compute `x % s` via an explicit shift: if `x − c·s ∈ [0, s)` then `x % s = x − c·s`. -/
theorem emod_eq_sub_of_range (x s c : ℤ) (h0 : 0 ≤ x - c * s) (h1 : x - c * s < s) :
    x % s = x - c * s := by
  conv_lhs => rw [show x = (x - c * s) + s * c by ring]
  rw [Int.add_mul_emod_self_left, Int.emod_eq_of_lt h0 h1]

/-! ## The atom -/

/-- **THE RATIONAL-POINT MARGIN CERTIFICATE.**  If `μ ≤ (v i · k) % s ≤ s − μ` for
every runner (a pure integer, decidable condition), then at `t = k/s` every runner is
at distance `≥ μ/s` from every integer.  Integer core: writing `v·k = q·s + r`,
`|v·k − m·s| = |(q−m)·s + r| ≥ min(r, s−r) ≥ μ`. -/
theorem rational_point_margin {ι : Type*} (v : ι → ℤ) (k s μ : ℤ) (hs : 0 < s)
    (hcond : ∀ i, μ ≤ (v i * k) % s ∧ (v i * k) % s ≤ s - μ) :
    ∀ i, ∀ m : ℤ, (μ : ℝ) / s ≤ |(v i : ℝ) * ((k : ℝ) / s) - m| := by
  intro i m
  obtain ⟨h1, h2⟩ := hcond i
  -- the integer core: μ ≤ |v i * k − m * s|
  have hint : μ ≤ |v i * k - m * s| := by
    set r : ℤ := (v i * k) % s with hr
    set q : ℤ := (v i * k) / s with hq
    have hde : s * q + r = v i * k := Int.mul_ediv_add_emod (v i * k) s
    have key : v i * k - m * s = (q - m) * s + r := by rw [← hde]; ring
    rcases le_or_gt m q with hmq | hmq
    · -- m ≤ q: the value is ≥ r ≥ μ
      have h0 : 0 ≤ (q - m) * s := mul_nonneg (by omega) hs.le
      have hval : μ ≤ v i * k - m * s := by rw [key]; linarith
      exact le_trans hval (le_abs_self _)
    · -- m ≥ q+1: the value is ≤ −s + r ≤ −μ
      have hle : (q - m) * s ≤ (-1) * s :=
        mul_le_mul_of_nonneg_right (by omega) hs.le
      have hval : v i * k - m * s ≤ -μ := by rw [key]; linarith
      have habs : μ ≤ -(v i * k - m * s) := by linarith
      exact le_trans habs (neg_le_abs _)
  -- cast to ℝ and divide by s
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  have hne : (s : ℝ) ≠ 0 := ne_of_gt hsR
  have heq : (v i : ℝ) * ((k : ℝ) / s) - m = ((v i * k - m * s : ℤ) : ℝ) / s := by
    push_cast
    field_simp
  have hXR : (μ : ℝ) ≤ |((v i * k - m * s : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hint
  rw [heq, abs_div, abs_of_pos hsR, div_eq_mul_inv, div_eq_mul_inv]
  exact mul_le_mul_of_nonneg_right hXR (inv_nonneg.mpr hsR.le)

/-- **The dichotomy-shaped consumable**: a certificate `(k, s, μ)` for the base
(all runners except `istar`) with `β ≤ μ/s` produces exactly the loose branch of
`TightLooseDichotomyAt β`. -/
theorem loose_branch_of_cert (v : Fin 13 → ℤ) (istar : Fin 13) (β : ℝ) (k s μ : ℤ)
    (hs : 0 < s) (hβ : β ≤ (μ : ℝ) / s)
    (hcond : ∀ i, i ≠ istar → μ ≤ (v i * k) % s ∧ (v i * k) % s ≤ s - μ) :
    ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m| := by
  refine ⟨(k : ℝ) / s, fun i hi m => ?_⟩
  calc β ≤ (μ : ℝ) / s := hβ
    _ ≤ |(v i : ℝ) * ((k : ℝ) / s) - m| :=
        rational_point_margin (fun j : {j : Fin 13 // j ≠ istar} => v j) k s μ hs
          (fun j => hcond j j.2) ⟨i, hi⟩ m

/-! ## The k = 1 instance: band residues and the ratio gate -/

/-- Band residues, first harmonic: `a ≤ |w| ≤ b` with `0 < a` puts `w % (a+b)`
inside `[a, b] = [μ, s−μ]` for `μ = a`. -/
theorem emod_band (w a b : ℤ) (ha : 0 < a) (hlo : a ≤ |w|) (hhi : |w| ≤ b) :
    a ≤ w % (a + b) ∧ w % (a + b) ≤ (a + b) - a := by
  rcases le_or_gt 0 w with hw0 | hw0
  · have habs : |w| = w := abs_of_nonneg hw0
    rw [habs] at hlo hhi
    rw [emod_eq_sub_of_range w (a + b) 0 (by omega) (by omega)]
    omega
  · have habs : |w| = -w := abs_of_neg hw0
    rw [habs] at hlo hhi
    rw [emod_eq_sub_of_range w (a + b) (-1) (by omega) (by omega)]
    omega

/-- **THE RATIO GATE (k = 1).**  A narrow base band — `β·(a+b) ≤ a` — yields the
loose branch of the dichotomy outright, at the reciprocal point `t = 1/(a+b)`.
At `β = 2/25` the gate condition is `2b ≤ 23a` (ratio `≤ 23/2`). -/
theorem loose_branch_of_narrow_band (v : Fin 13 → ℤ) (istar : Fin 13) (β : ℝ)
    (a b : ℤ) (ha : 0 < a)
    (hlo : ∀ i, i ≠ istar → a ≤ |v i|) (hhi : ∀ i, i ≠ istar → |v i| ≤ b)
    (hab : a ≤ b) (hgate : β * ((a : ℝ) + b) ≤ a) :
    ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m| := by
  have hs : (0 : ℤ) < a + b := by omega
  have hsR : (0 : ℝ) < ((a + b : ℤ) : ℝ) := by exact_mod_cast hs
  refine loose_branch_of_cert v istar β 1 (a + b) a hs ?_ ?_
  · rw [le_div_iff₀ hsR]
    push_cast
    linarith
  · intro i hi
    have := emod_band (v i) a b ha (hlo i hi) (hhi i hi)
    rw [mul_one]
    exact this

/-! ## The k = 2 instance: the second harmonic with the middle hole -/

/-- Hole residues, second harmonic: `a ≤ |w| ≤ b` avoiding the middle window —
`2|w| ≤ b − a` or `2|w| ≥ 3a + b` — puts `(w·2) % (a+b)` inside
`[2a, (a+b) − 2a]`. -/
theorem emod_hole (w a b : ℤ) (ha : 0 < a) (hlo : a ≤ |w|) (hhi : |w| ≤ b)
    (hhole : 2 * |w| ≤ b - a ∨ 3 * a + b ≤ 2 * |w|) :
    2 * a ≤ (w * 2) % (a + b) ∧ (w * 2) % (a + b) ≤ (a + b) - 2 * a := by
  rcases le_or_gt 0 w with hw0 | hw0
  · have habs : |w| = w := abs_of_nonneg hw0
    rw [habs] at hlo hhi hhole
    rcases hhole with hlow | hhigh
    · rw [emod_eq_sub_of_range (w * 2) (a + b) 0 (by omega) (by omega)]
      omega
    · rw [emod_eq_sub_of_range (w * 2) (a + b) 1 (by omega) (by omega)]
      omega
  · have habs : |w| = -w := abs_of_neg hw0
    rw [habs] at hlo hhi hhole
    rcases hhole with hlow | hhigh
    · rw [emod_eq_sub_of_range (w * 2) (a + b) (-1) (by omega) (by omega)]
      omega
    · rw [emod_eq_sub_of_range (w * 2) (a + b) (-2) (by omega) (by omega)]
      omega

/-- **THE SECOND-HARMONIC GATE (k = 2).**  A base band `[a,b]` avoiding the middle
window `(b−a, 3a+b)` in doubled absolute value yields the loose branch at
`t = 2/(a+b)` whenever `β·(a+b) ≤ 2a` (ratio `≤ 24` at `β = 2/25`).  This closes
the rigidity spectrum's second value: `{1,…,11,24}` at `t = 2/25`, margin `2/25`. -/
theorem loose_branch_of_hole (v : Fin 13 → ℤ) (istar : Fin 13) (β : ℝ)
    (a b : ℤ) (ha : 0 < a)
    (hlo : ∀ i, i ≠ istar → a ≤ |v i|) (hhi : ∀ i, i ≠ istar → |v i| ≤ b)
    (hhole : ∀ i, i ≠ istar → 2 * |v i| ≤ b - a ∨ 3 * a + b ≤ 2 * |v i|)
    (hab : a ≤ b) (hgate : β * ((a : ℝ) + b) ≤ 2 * a) :
    ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m| := by
  have hs : (0 : ℤ) < a + b := by omega
  have hsR : (0 : ℝ) < ((a + b : ℤ) : ℝ) := by exact_mod_cast hs
  refine loose_branch_of_cert v istar β 2 (a + b) (2 * a) hs ?_ ?_
  · rw [le_div_iff₀ hsR]
    push_cast
    linarith
  · intro i hi
    exact emod_hole (v i) a b ha (hlo i hi) (hhi i hi) (hhole i hi)

/-! ## The reduction: the dichotomy only needs proving on spread bases -/

/-- **THE RATIO-GATE REDUCTION.**  `TightLooseDichotomyAt β` holds outright if the
dichotomy disjunction holds for every family whose base is SPREAD: base min `a`
and max `b` with `a < β·(a+b)` (ratio `> (1−β)/β`).  Narrow bases are loose by the
`k = 1` gate.  At `β = 2/25`: only bases with ratio `> 23/2` remain. -/
theorem tightLooseDichotomyAt_of_spread (β : ℝ)
    (hspread : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      ∀ a b : ℤ, 0 < a → a ≤ b →
        (∀ i, i ≠ istar → a ≤ |v i|) → (∀ i, i ≠ istar → |v i| ≤ b) →
        (a : ℝ) < β * ((a : ℝ) + b) →
        (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
        (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)) :
    TightLooseDichotomyAt β := by
  intro v hv hcov hcomp hprim istar hmax
  -- base min and max over the 12 non-killer runners
  have hne : (Finset.univ.erase istar).Nonempty := by
    rw [← Finset.card_pos, Finset.card_erase_of_mem (Finset.mem_univ istar)]
    simp
  obtain ⟨imin, himmem, himmin⟩ :=
    Finset.exists_min_image (Finset.univ.erase istar) (fun i => |v i|) hne
  obtain ⟨imax, hixmem, hixmax⟩ :=
    Finset.exists_max_image (Finset.univ.erase istar) (fun i => |v i|) hne
  set a : ℤ := |v imin| with hadef
  set b : ℤ := |v imax| with hbdef
  have ha : 0 < a := abs_pos.mpr (hv imin)
  have hlo : ∀ i, i ≠ istar → a ≤ |v i| := fun i hi =>
    himmin i (Finset.mem_erase.mpr ⟨hi, Finset.mem_univ i⟩)
  have hhi : ∀ i, i ≠ istar → |v i| ≤ b := fun i hi =>
    hixmax i (Finset.mem_erase.mpr ⟨hi, Finset.mem_univ i⟩)
  have hab : a ≤ b := hlo imax (Finset.ne_of_mem_erase hixmem)
  by_cases hgate : β * ((a : ℝ) + b) ≤ a
  · exact Or.inr (loose_branch_of_narrow_band v istar β a b ha hlo hhi hab hgate)
  · exact hspread v hv hcov hcomp hprim istar hmax a b ha hab hlo hhi
      (not_le.mp hgate)

/-- The `β = 2/25` bridge to the surface's non-parametric predicate. -/
theorem tightLooseDichotomy_of_at (h : TightLooseDichotomyAt (2 / 25)) :
    TightLooseDichotomy :=
  h

/-- **The concrete 2/25 reduction**: the surface's `TightLooseDichotomy` from the
spread case only — bases with `23a < 2b` (ratio `> 23/2`). -/
theorem tightLooseDichotomy_of_spread
    (hspread : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      ∀ a b : ℤ, 0 < a → a ≤ b →
        (∀ i, i ≠ istar → a ≤ |v i|) → (∀ i, i ≠ istar → |v i| ≤ b) →
        23 * a < 2 * b →
        (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
        (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2 / 25 ≤ |(v i : ℝ) * tstar - m|)) :
    TightLooseDichotomy := by
  apply tightLooseDichotomy_of_at
  apply tightLooseDichotomyAt_of_spread
  intro v hv hcov hcomp hprim istar hmax a b ha hab hlo hhi hgap
  refine hspread v hv hcov hcomp hprim istar hmax a b ha hab hlo hhi ?_
  -- a < (2/25)(a+b)  ⟹  25a < 2a + 2b  ⟹  23a < 2b
  have h25 : (25 : ℝ) * a < 2 * ((a : ℝ) + b) := by linarith
  have : (23 : ℝ) * a < 2 * b := by linarith
  exact_mod_cast this

/-- **The `β = 1/13` subsumption**: the tight branch (values in `c·{1,…,12}`)
forces base ratio `≤ 12`, which passes the 1/13-gate — so `TightLooseDichotomyAt
(1/13)` needs NO tight branch at all: it follows from spread bases (`12a < b`)
having a `1/13`-margin point.  The 1/13-dichotomy is purely a spread-base margin
statement. -/
theorem dichotomyAt_13_of_spread_loose
    (hspread : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      ∀ a b : ℤ, 0 < a → a ≤ b →
        (∀ i, i ≠ istar → a ≤ |v i|) → (∀ i, i ≠ istar → |v i| ≤ b) →
        12 * a < b →
        ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 1 / 13 ≤ |(v i : ℝ) * tstar - m|) :
    TightLooseDichotomyAt (1 / 13) := by
  apply tightLooseDichotomyAt_of_spread
  intro v hv hcov hcomp hprim istar hmax a b ha hab hlo hhi hgap
  refine Or.inr (hspread v hv hcov hcomp hprim istar hmax a b ha hab hlo hhi ?_)
  -- a < (1/13)(a+b)  ⟹  13a < a + b  ⟹  12a < b
  have h13 : (13 : ℝ) * a < (a : ℝ) + b := by linarith
  have : (12 : ℝ) * a < b := by linarith
  exact_mod_cast this

/-! ## Propagation to the pinned surface -/

/-- **THE SHARPENED SURFACE**: LRC(14) ≤ citation + `CornerLonely` + the dichotomy
on SPREAD bases only (`23a < 2b`, ratio `> 23/2`).  The ratio gate discharges every
narrow-based compressed family outright; the AP tight locus (ratio `12`) and the
second value (ratio `24`) both live on the spread side, where the remaining open
mathematics now provably is. -/
theorem lrc14_of_spread_dichotomy_and_corner (cite : LonelyRunner.LRCUpTo13)
    (hspread : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      ∀ a b : ℤ, 0 < a → a ≤ b →
        (∀ i, i ≠ istar → a ≤ |v i|) → (∀ i, i ≠ istar → |v i| ≤ b) →
        23 * a < 2 * b →
        (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
        (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2 / 25 ≤ |(v i : ℝ) * tstar - m|))
    (hcorner : CornerLonely) :
    LRC14Statement :=
  lrc14_of_dichotomy_and_corner cite (tightLooseDichotomy_of_spread hspread) hcorner

/-- **The β-parametric sharpened surface**: the same propagation at any loose margin
`β > 1/14` — spread threshold `a < β·(a+b)`, corner threshold `B/(14(β−1/14))`. -/
theorem lrc14_of_spread_dichotomy_and_corner_at (cite : LonelyRunner.LRCUpTo13)
    (β : ℝ) (hβ : 1 / 14 < β)
    (hspread : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
      ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      ∀ a b : ℤ, 0 < a → a ≤ b →
        (∀ i, i ≠ istar → a ≤ |v i|) → (∀ i, i ≠ istar → |v i| ≤ b) →
        (a : ℝ) < β * ((a : ℝ) + b) →
        (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
        (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|))
    (hcorner : CornerLonelyAt β) :
    LRC14Statement :=
  lrc14_of_dichotomy_and_corner_at cite β hβ
    (tightLooseDichotomyAt_of_spread β hspread) hcorner

/-! ## The second value, formal -/

/-- **The rigidity spectrum's second value as a theorem**: any base shaped like
`±c·{1,…,11,24}` — more generally `|v i| ∈ [a, 24a]` avoiding the middle window
`(11a, 14a)` — is loose at margin `2/25`, at `t = 2/(25a)`.  The gate inequality
holds with EQUALITY (`(2/25)·25a = 2a`): the second value is the second harmonic's
extremal instance, exactly as the AP is the first harmonic's. -/
theorem second_value_loose (v : Fin 13 → ℤ) (istar : Fin 13) (a : ℤ) (ha : 0 < a)
    (hlo : ∀ i, i ≠ istar → a ≤ |v i|) (hhi : ∀ i, i ≠ istar → |v i| ≤ 24 * a)
    (hhole : ∀ i, i ≠ istar → |v i| ≤ 11 * a ∨ 14 * a ≤ |v i|) :
    ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2 / 25 ≤ |(v i : ℝ) * tstar - m| := by
  apply loose_branch_of_hole v istar (2 / 25) a (24 * a) ha hlo hhi ?_ (by omega) ?_
  · intro i hi
    rcases hhole i hi with h | h
    · left; omega
    · right; omega
  · push_cast
    linarith

/-! ## The sieve as the degenerate instance -/

/-- The sieve shape recovered from the atom (`μ = 1`, `k = 1`): if `s ∤ v i` for
all runners, margin `1/s` at `t = 1/s`.  Unifies `sieve_one_div` with the band
lemmas: sieve, reciprocal band, and harmonics are ONE primitive. -/
theorem sieve_margin {ι : Type*} (v : ι → ℤ) (s : ℤ) (hs : 1 < s)
    (hnd : ∀ i, ¬ (s ∣ v i)) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / s ≤ |(v i : ℝ) * ((1 : ℝ) / s) - m| := by
  have hs0 : (0 : ℤ) < s := by omega
  have h := rational_point_margin v 1 s 1 hs0 (fun i => by
    have hr0 : (v i * 1) % s ≠ 0 := by
      rw [mul_one]
      exact fun h => hnd i (Int.dvd_of_emod_eq_zero h)
    have hnn : 0 ≤ (v i * 1) % s := Int.emod_nonneg _ (by omega)
    have hlt : (v i * 1) % s < s := Int.emod_lt_of_pos _ hs0
    omega)
  intro i m
  have := h i m
  simpa using this

#print axioms rational_point_margin
#print axioms tightLooseDichotomyAt_of_spread
#print axioms tightLooseDichotomy_of_spread
#print axioms loose_branch_of_hole
#print axioms lrc14_of_spread_dichotomy_and_corner
#print axioms lrc14_of_spread_dichotomy_and_corner_at
#print axioms second_value_loose

end HarmonicGate
end LonelyRunner
