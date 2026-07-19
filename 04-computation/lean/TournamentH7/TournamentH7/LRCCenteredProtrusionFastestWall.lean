import Mathlib.Tactic

/-!
# Centered-protrusion / fastest-wall arithmetic (THM-1273)

This is the arithmetic and propositional consumer for the paper bridge from
THM-1267's endpoint protrusion to a phase-located fastest wall.  The paper
providers are THM-1198's strict one-comb load bounds, containment of the
five-comb survivor in the protrusion, selection of regular points away from
the finite wall set, continuity of the fastest circle-depth function, and
the tooth endpoint calculation.

The module checks the exact `2/9`, `1/8`, and `7/72` mass ledger, the
bare-wall-versus-low-crosser logic, the inherited lcm quantum, and all rational
comparisons used by the `c=140` positive control.  It introduces no covering
axiom and contains no `sorry` or `native_decide`.
-/

namespace LRC14.CenteredProtrusionFastestWall

def oneCombCap : ℚ := 7 / 36
def fourPrefixFloor : ℚ := 2 / 9
def smallTailCap : ℚ := 1 / 8
def insideFloor : ℚ := 7 / 72
def fivePrefixFloor : ℚ := 11 / 360

theorem bridge_constant_ledger :
    4 * oneCombCap = 7 / 9 ∧
      1 - 4 * oneCombCap = fourPrefixFloor ∧
      (3 / 4 : ℚ) * (1 / 6) = smallTailCap ∧
      fourPrefixFloor - smallTailCap = insideFloor ∧
      1 - (4 * oneCombCap + 23 / 120) = fivePrefixFloor := by
  norm_num [oneCombCap, fourPrefixFloor, smallTailCap, insideFloor,
    fivePrefixFloor]

/-- Four strict one-comb load bounds leave strictly more than `2/9` of the
six-bin density for the four-prefix survivor. -/
theorem four_loads_leave_two_ninths
    (a b c d : ℚ)
    (ha : a < oneCombCap) (hb : b < oneCombCap)
    (hc : c < oneCombCap) (hd : d < oneCombCap) :
    fourPrefixFloor < 1 - (a + b + c + d) := by
  dsimp [oneCombCap, fourPrefixFloor] at *
  linarith

/-- In the small-tail branch the endpoint sixth has constant density `3/4`,
so the complete protrusion has mass below `1/8`. -/
theorem small_endpoint_tail_mass
    (ell tailMass : ℚ) (hell : ell < 1 / 6)
    (hmass : tailMass = (3 / 4) * ell) :
    tailMass < smallTailCap := by
  dsimp [smallTailCap]
  rw [hmass]
  linarith

/-- If the total four-prefix survivor has mass above `2/9` and at most the
tail mass can lie outside `K`, then more than `7/72` remains inside `K`. -/
theorem four_prefix_inside_K
    (survivorMass tailMass insideMass : ℚ)
    (hsurvivor : fourPrefixFloor < survivorMass)
    (htail : tailMass < smallTailCap)
    (hpartition : survivorMass ≤ insideMass + tailMass) :
    insideFloor < insideMass := by
  dsimp [fourPrefixFloor, smallTailCap, insideFloor] at *
  linarith

/-- The separated fifth load leaves THM-1267's strict outside obligation. -/
theorem five_prefix_survivor_positive
    (load : ℚ) (hload : load < 349 / 360) :
    fivePrefixFloor < 1 - load := by
  dsimp [fivePrefixFloor]
  linarith

/-- The density maximum converts the inside mass obligation to normalized
Lebesgue length greater than `1/12`. -/
theorem inside_obligation_length
    (mass length : ℚ) (hmass : insideFloor < mass)
    (hdensity : mass ≤ (7 / 6) * length) :
    1 / 12 < length := by
  dsimp [insideFloor] at hmass
  linarith

/-- On the endpoint sixth the density is exactly `3/4`, so the outside
five-prefix survivor has normalized length greater than `11/270`. -/
theorem outside_obligation_length
    (mass length : ℚ) (hmass : fivePrefixFloor < mass)
    (hdensity : mass = (3 / 4) * length) :
    11 / 270 < length := by
  dsimp [fivePrefixFloor] at hmass
  rw [hdensity] at hmass
  linarith

theorem needle_separation_constant :
    (1 : ℚ) / 12 + 11 / 270 = 67 / 540 := by
  norm_num

/-- Choosing one obligation deeper than its measure from each side of the
interface yields the exact `67/540` oriented needle separation. -/
theorem obligation_points_separated
    (insideDistance outsideDistance : ℚ)
    (hinside : 1 / 12 < insideDistance)
    (houtside : 11 / 270 < outsideDistance) :
    67 / 540 < insideDistance + outsideDistance := by
  linarith

/-- If `W` fastest walls split the oriented needle into `W+1` wall-free
cells, each of normalized length at most `d1/h`, the `67/540` separation
forces the displayed exact integer wall-count invoice. -/
theorem wall_count_invoice
    (W : ℕ) (d₁ h needleLength : ℚ)
    (hh : 0 < h)
    (hlength : 67 / 540 < needleLength)
    (hcells : needleLength ≤ ((W : ℚ) + 1) * d₁ / h) :
    67 * h < 540 * ((W : ℚ) + 1) * d₁ := by
  have hmul : needleLength * h ≤ ((W : ℚ) + 1) * d₁ :=
    (le_div_iff₀ hh).mp hcells
  nlinarith

/-- At ratio `h/d1 >= 1080/67`, the needle cannot contain only one fastest
wall; this is the first useful multi-event threshold. -/
theorem two_walls_at_large_fastest_ratio
    (W : ℕ) (d₁ h : ℚ)
    (hd₁ : 0 < d₁)
    (hinvoice : 67 * h < 540 * ((W : ℚ) + 1) * d₁)
    (hratio : 1080 * d₁ ≤ 67 * h) :
    2 ≤ W := by
  by_contra hnot
  have hW : W ≤ 1 := by omega
  have hcast : (W : ℚ) ≤ 1 := by exact_mod_cast hW
  nlinarith

/-- At a fastest wall the slowest component owner and the fastest owner are
both strict-safe.  If no one of the four lower owners crosses the wall, the
cover predicate excludes the wall from `K`; along the oriented bridge it must
therefore be in the protrusion. -/
theorem bare_wall_or_low_crosser
    (inK inTail d₁Danger d₂Danger d₃Danger d₄Danger d₅Danger hDanger : Prop)
    (hlocation : inK ∨ inTail)
    (hcover : inK →
      d₁Danger ∨ d₂Danger ∨ d₃Danger ∨ d₄Danger ∨ d₅Danger ∨ hDanger)
    (hd₁safe : ¬d₁Danger) (hwallSafe : ¬hDanger) :
    (inTail ∧ ¬d₂Danger ∧ ¬d₃Danger ∧ ¬d₄Danger ∧ ¬d₅Danger) ∨
      d₂Danger ∨ d₃Danger ∨ d₄Danger ∨ d₅Danger := by
  by_cases hd₂ : d₂Danger
  · exact Or.inr (Or.inl hd₂)
  by_cases hd₃ : d₃Danger
  · exact Or.inr (Or.inr (Or.inl hd₃))
  by_cases hd₄ : d₄Danger
  · exact Or.inr (Or.inr (Or.inr (Or.inl hd₄)))
  by_cases hd₅ : d₅Danger
  · exact Or.inr (Or.inr (Or.inr (Or.inr hd₅)))
  left
  refine ⟨?_, hd₂, hd₃, hd₄, hd₅⟩
  rcases hlocation with hinK | hinTail
  · exfalso
    rcases hcover hinK with hd₁ | hd₂' | hd₃' | hd₄' | hd₅' | hh
    · exact hd₁safe hd₁
    · exact hd₂ hd₂'
    · exact hd₃ hd₃'
    · exact hd₄ hd₄'
    · exact hd₅ hd₅'
    · exact hwallSafe hh
  · exact hinTail

/-- Re-export the exact endpoint arithmetic consumer: a positive compatible
overlap numerator divisible by the owner gcd pays `1/(14*lcm)`. -/
theorem crossed_wall_pays_lcm_quantum
    {h j numerator : ℕ} (hh : 0 < h) (hj : 0 < j)
    (hpos : 0 < numerator) (hdiv : Nat.gcd h j ∣ numerator) :
    (1 : ℚ) / (14 * Nat.lcm h j) ≤ numerator / (14 * h * j) := by
  have hl : 0 < Nat.lcm h j := Nat.lcm_pos hh hj
  have hleNat : Nat.gcd h j ≤ numerator := Nat.le_of_dvd hpos hdiv
  have hle : (Nat.gcd h j : ℚ) ≤ numerator := by exact_mod_cast hleNat
  have hidentity :
      (Nat.gcd h j : ℚ) / (14 * h * j) =
        1 / (14 * Nat.lcm h j) := by
    have hh0 : (h : ℚ) ≠ 0 := by positivity
    have hj0 : (j : ℚ) ≠ 0 := by positivity
    have hl0 : (Nat.lcm h j : ℚ) ≠ 0 := by positivity
    field_simp
    exact_mod_cast Nat.gcd_mul_lcm h j
  rw [← hidentity]
  have hden : (0 : ℚ) < 14 * h * j := by positivity
  exact (div_le_div_iff_of_pos_right hden).2 hle

theorem sharp_control_tail :
    (33 : ℚ) / 280 < 1 / 6 ∧
      (11 : ℚ) / 270 < 33 / 280 := by
  norm_num

theorem sharp_control_oriented_points :
    (2045 : ℚ) / 3556 < 7476011 / 12938240 ∧
      (7476011 : ℚ) / 12938240 < 1133 / 1960 ∧
      (1133 : ℚ) / 1960 < 7425603 / 12837160 ∧
      (7425603 : ℚ) / 12837160 < 2057 / 3556 := by
  norm_num

theorem sharp_control_wall_order :
    (7476011 : ℚ) / 12938240 < 14603 / 25270 ∧
      (14603 : ℚ) / 25270 < 2923 / 5054 ∧
      (2923 : ℚ) / 5054 < 14617 / 25270 ∧
      (14617 : ℚ) / 25270 < 7425603 / 12837160 := by
  norm_num

/-- The first two fastest walls in the sharp row lie strictly inside the
`256@148` tooth; the last wall lies after that tooth. -/
theorem sharp_control_low_wall_incidence :
    (2071 : ℚ) / 3584 < 14603 / 25270 ∧
      (14603 : ℚ) / 25270 < 2073 / 3584 ∧
      (2071 : ℚ) / 3584 < 2923 / 5054 ∧
      (2923 : ℚ) / 5054 < 2073 / 3584 ∧
      (2073 : ℚ) / 3584 < 14617 / 25270 := by
  norm_num

/-- At the last sharp-row wall, `h=1805` is exactly safe and `d1,d2,...,d5`
are all strictly beyond the `1/14` danger radius. -/
theorem sharp_control_bare_wall_depths :
    (1805 : ℚ) * (14617 / 25270) - 1044 = 1 / 14 ∧
      |(254 : ℚ) * (14617 / 25270) - 147| = 986 / 12635 ∧
      |(255 : ℚ) * (14617 / 25270) - 148| = 2525 / 5054 ∧
      |(256 : ℚ) * (14617 / 25270) - 148| = 996 / 12635 ∧
      |(257 : ℚ) * (14617 / 25270) - 149| = 8661 / 25270 ∧
      |(292 : ℚ) * (14617 / 25270) - 169| = 1233 / 12635 ∧
      1 / 14 < (986 : ℚ) / 12635 ∧
      1 / 14 < (2525 : ℚ) / 5054 ∧
      1 / 14 < (996 : ℚ) / 12635 ∧
      1 / 14 < (8661 : ℚ) / 25270 ∧
      1 / 14 < (1233 : ℚ) / 12635 := by
  norm_num [abs_of_nonneg, abs_of_neg]

theorem sharp_control_bare_wall_in_tail :
    (1133 : ℚ) / 1960 < 14617 / 25270 ∧
      (14617 : ℚ) / 25270 < 2057 / 3556 := by
  norm_num

/-- The two crossed walls in the sharp row pay respectively 213 and 325
copies of the primitive `lcm(1805,256)=1805*256` seam quantum. -/
theorem sharp_control_seam_quanta :
    (213 : ℚ) / 6469120 = 213 * (1 / 6469120) ∧
      (65 : ℚ) / 1293824 = 325 * (1 / 6469120) ∧
      1 / 6469120 = 1 / (14 * (1805 * 256)) := by
  norm_num

#print axioms bridge_constant_ledger
#print axioms four_loads_leave_two_ninths
#print axioms small_endpoint_tail_mass
#print axioms four_prefix_inside_K
#print axioms five_prefix_survivor_positive
#print axioms inside_obligation_length
#print axioms outside_obligation_length
#print axioms needle_separation_constant
#print axioms obligation_points_separated
#print axioms wall_count_invoice
#print axioms two_walls_at_large_fastest_ratio
#print axioms bare_wall_or_low_crosser
#print axioms crossed_wall_pays_lcm_quantum
#print axioms sharp_control_tail
#print axioms sharp_control_oriented_points
#print axioms sharp_control_wall_order
#print axioms sharp_control_low_wall_incidence
#print axioms sharp_control_bare_wall_depths
#print axioms sharp_control_bare_wall_in_tail
#print axioms sharp_control_seam_quanta

end LRC14.CenteredProtrusionFastestWall
