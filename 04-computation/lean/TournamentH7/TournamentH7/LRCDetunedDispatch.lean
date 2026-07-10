/-
  TournamentH7.LRCDetunedDispatch — THM-668 in Lean (monad-explorer-2026-07-09-S7, HYP-5777).

  THE DETUNED-HARMONIC DISPATCH (d = 1): a 13-family in which one modulus `g ≥ 2` divides
  every speed EXCEPT one is lonely, unconditionally — from the LRC(≤13) citation.

  MECHANISM (THM-668, monad-S3; formalized via the TRIANGLE SHORTCUT):
   · the 12 harmonic speeds are `g·w`; the citation gives `u` with `Lonely 13 w u`; at every
     branch time `τ_c = (u + c)/g` (c : ℤ) the harmonic clearance is UNCHANGED:
     `(g·w)·τ_c = w·u + (w·c)`, an integer shift — clearance stays ≥ 1/13.
   · the detuned speed `δ = v i₀` (`g ∤ δ`): with `d = gcd(δ, g)`, `q = g/d ≥ 2`, the window
     multiple `m = d·(q/2)` satisfies `g ≤ 4m ∧ 2m ≤ g`, so `m/g ∈ [1/4, 1/2]` and every
     integer is ≥ 1/4 from `m/g` (`quarter_window`).  Bezout gives `c` with
     `c·δ = m − g·e`, so the branch phases `δ·τ_0`, `δ·τ_c` differ by `m/g` mod 1, and by
     the triangle inequality one branch has detuned clearance ≥ 1/8 (`branch_pigeonhole`).
   · at that branch every runner clears ≥ min(1/13, 1/8) = 1/13 > 1/14.

  Sorry-free; kernel-pure target (no native_decide).
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace DetunedDispatch

open LonelyRunner

/-- The `m/g ∈ [1/4, 1/2]` window keeps every integer at distance ≥ 1/4. -/
theorem quarter_window (g m : ℤ) (hg : 0 < g) (h4 : g ≤ 4 * m) (h2 : 2 * m ≤ g) :
    ∀ n : ℤ, (1 : ℝ) / 4 ≤ |(m : ℝ) / g - n| := by
  intro n
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg
  have hlo : (1 : ℝ) / 4 ≤ (m : ℝ) / g := by
    rw [le_div_iff₀ hgR]
    have h4R : (g : ℝ) ≤ 4 * m := by exact_mod_cast h4
    linarith
  have hhi : (m : ℝ) / g ≤ 1 / 2 := by
    rw [div_le_iff₀ hgR]
    have h2R : (2 : ℝ) * m ≤ g := by exact_mod_cast h2
    linarith
  by_cases hn : (n : ℝ) ≤ 0
  · rw [abs_of_nonneg (by linarith)]
    linarith
  · have hn1 : (1 : ℝ) ≤ n := by
      have : (0 : ℤ) < n := by exact_mod_cast not_le.mp hn
      exact_mod_cast this
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- **The branch-pair pigeonhole (triangle shortcut).**  If two reals differ by
`m/g + integer` with `m/g` in the quarter window, one of them is ≥ 1/8 from every
integer. -/
theorem branch_pigeonhole (x y : ℝ) (g m e : ℤ) (hg : 0 < g)
    (h4 : g ≤ 4 * m) (h2 : 2 * m ≤ g)
    (hdiff : y = x + (m : ℝ) / g + (e : ℝ)) :
    (∀ n : ℤ, (1 : ℝ) / 8 ≤ |x - n|) ∨ (∀ n : ℤ, (1 : ℝ) / 8 ≤ |y - n|) := by
  by_contra hcon
  push_neg at hcon
  obtain ⟨⟨n₁, hn₁⟩, ⟨n₂, hn₂⟩⟩ := hcon
  have key := quarter_window g m hg h4 h2 (n₂ - e - n₁)
  have hexp : (m : ℝ) / g - ((n₂ - e - n₁ : ℤ) : ℝ) = (y - n₂) + -(x - n₁) := by
    rw [hdiff]; push_cast; ring
  have htri : |(m : ℝ) / g - ((n₂ - e - n₁ : ℤ) : ℝ)| ≤ |y - n₂| + |x - n₁| := by
    rw [hexp]
    calc |(y - (n₂:ℝ)) + -(x - n₁)| ≤ |y - (n₂:ℝ)| + |-(x - (n₁:ℝ))| := abs_add_le _ _
    _ = |y - (n₂:ℝ)| + |x - (n₁:ℝ)| := by rw [abs_neg]
  linarith

/-- **THM-668 (Lean): the detuned-harmonic dispatch.**  If `g ≥ 2` divides every speed
except the one at `i₀` (which it does not divide), the family is lonely — with an
explicit witness `τ ∈ {u/g, (u+c)/g}` built from the LRC(13) time `u` of the quotient
family and a Bezout branch shift `c`. -/
theorem lonely14_of_detuned (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g) (i₀ : Fin 13)
    (hH : ∀ j, j ≠ i₀ → g ∣ v j) (hδ : ¬ g ∣ v i₀) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg0 : (0 : ℤ) < g := by omega
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  have hgR' : (g : ℝ) ≠ 0 := ne_of_gt hgR
  -- the quotient family on the 12 harmonics
  set w : Fin 12 → ℤ := fun k => v (i₀.succAbove k) / g with hw
  have hwv : ∀ k, v (i₀.succAbove k) = g * w k := by
    intro k
    exact (Int.mul_ediv_cancel' (hH _ (Fin.succAbove_ne i₀ k))).symm
  have hwnz : ∀ k, w k ≠ 0 := by
    intro k hk
    exact hv (i₀.succAbove k) (by rw [hwv k, hk, mul_zero])
  -- the citation: an LRC(13) time for the quotient family
  obtain ⟨u, hu⟩ := cite 12 (by norm_num) w hwnz
  -- the detuned data
  set δ : ℤ := v i₀ with hδdef
  have hδnz : δ ≠ 0 := hv i₀
  set d : ℤ := (Int.gcd δ g : ℤ) with hd
  have hdpos : (0 : ℤ) < d := by
    rw [hd]
    have hne : Int.gcd δ g ≠ 0 := fun h => hδnz (Int.gcd_eq_zero_iff.mp h).1
    exact_mod_cast Nat.pos_of_ne_zero hne
  have hddvdg : d ∣ g := by rw [hd]; exact Int.gcd_dvd_right δ g
  have hddvdδ : d ∣ δ := by rw [hd]; exact Int.gcd_dvd_left δ g
  set q : ℤ := g / d with hq
  have hgdq : g = d * q := (Int.mul_ediv_cancel' hddvdg).symm
  have hdq : (0 : ℤ) < d * q := by rw [← hgdq]; exact hg0
  have hqpos : (0 : ℤ) < q := by
    by_contra h
    push_neg at h
    nlinarith
  have hq2 : (2 : ℤ) ≤ q := by
    by_contra h
    push_neg at h
    have hq1 : q = 1 := by omega
    rw [hq1, mul_one] at hgdq
    exact hδ (hgdq ▸ hddvdδ)
  -- the window multiple m = d * (q / 2)
  set s : ℤ := q / 2 with hs
  have hsq : 2 * s + q % 2 = q := by
    have h := Int.ediv_add_emod q 2
    omega
  have hmod : q % 2 = 0 ∨ q % 2 = 1 := Int.emod_two_eq_zero_or_one q
  have hs1 : (1 : ℤ) ≤ s := by rcases hmod with h | h <;> omega
  set m : ℤ := d * s with hm
  have h2m : 2 * m ≤ g := by
    rw [hm, hgdq]
    rcases hmod with h | h <;> nlinarith
  have h4m : g ≤ 4 * m := by
    rw [hm, hgdq]
    rcases hmod with h | h <;> nlinarith
  -- Bezout: c with c·δ = m − g·(s·B)
  set A : ℤ := Int.gcdA δ g with hA
  set B : ℤ := Int.gcdB δ g with hB
  have hbez : d = δ * A + g * B := Int.gcd_eq_gcd_ab δ g
  set c : ℤ := s * A with hc
  have hkey : c * δ = m - g * (s * B) := by
    have : c * δ = s * (δ * A) := by rw [hc]; ring
    rw [this, show δ * A = d - g * B by linarith [hbez], hm]
    ring
  -- the two branch times
  set τ₀ : ℝ := u / g with hτ₀
  set τc : ℝ := (u + (c : ℝ)) / g with hτc
  -- detuned phases differ by m/g + integer
  have hph : (δ : ℝ) * τc = (δ : ℝ) * τ₀ + (m : ℝ) / g + ((-(s * B) : ℤ) : ℝ) := by
    have hcδR : ((c : ℝ) * δ) = (m : ℝ) - (g : ℝ) * ((s : ℝ) * B) := by
      have := hkey
      have hcast : ((c * δ : ℤ) : ℝ) = ((m - g * (s * B) : ℤ) : ℝ) := by
        exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) this
      push_cast at hcast
      linarith
    rw [hτc, hτ₀]
    push_cast
    field_simp
    linarith [hcδR]
  -- one branch clears the detuned runner by 1/8
  rcases branch_pigeonhole ((δ : ℝ) * τ₀) ((δ : ℝ) * τc) g m (-(s * B)) hg0 h4m h2m hph
    with hbr | hbr
  -- τ₀ fires
  · refine ⟨τ₀, fun i n => ?_⟩
    rcases eq_or_ne i i₀ with heq | hne
    · have h8 := hbr n
      rw [heq]
      calc (1 : ℝ) / 14 ≤ 1 / 8 := by norm_num
      _ ≤ |(δ : ℝ) * τ₀ - n| := h8
    · obtain ⟨k, hk⟩ := Fin.exists_succAbove_eq hne
      have hval : (v i : ℝ) * τ₀ - n = (w k : ℝ) * u - n := by
        rw [← hk, hwv k, hτ₀]
        push_cast
        field_simp
      rw [hval]
      calc (1 : ℝ) / 14 ≤ 1 / 13 := by norm_num
      _ ≤ |(w k : ℝ) * u - n| := hu k n
  -- τc fires
  · refine ⟨τc, fun i n => ?_⟩
    rcases eq_or_ne i i₀ with heq | hne
    · have h8 := hbr n
      rw [heq]
      calc (1 : ℝ) / 14 ≤ 1 / 8 := by norm_num
      _ ≤ |(δ : ℝ) * τc - n| := h8
    · obtain ⟨k, hk⟩ := Fin.exists_succAbove_eq hne
      have hval : (v i : ℝ) * τc - n = (w k : ℝ) * u - ((n - w k * c : ℤ) : ℝ) := by
        rw [← hk, hwv k, hτc]
        push_cast
        field_simp
        ring
      rw [hval]
      calc (1 : ℝ) / 14 ≤ 1 / 13 := by norm_num
      _ ≤ |(w k : ℝ) * u - ((n - w k * c : ℤ) : ℝ)| := hu k (n - w k * c)

/-! ## Axiom audit -/
#print axioms quarter_window
#print axioms branch_pigeonhole
#print axioms lonely14_of_detuned

end DetunedDispatch
end LonelyRunner
