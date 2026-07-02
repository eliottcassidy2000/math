/-
  TournamentH7.PolygonMirskyNewman  (mac-mini-2026-07-01-S95)

  THE REGULAR-POLYGON COLORING THEOREM (discrete Mirsky–Newman / Davenport–Rado):

    if the vertices of a regular n-gon are colored so that every color class
    is itself the vertex set of a regular polygon, and there are at least two
    classes, then two classes form CONGRUENT polygons.

  Formulation: vertices = `Finset.range n`; a class with difference `d ∣ n` and
  offset `a < d` is `{x < n : x % d = a}` (a regular (n/d)-gon).  A coloring
  `c` assigns each vertex its class; the theorem says the difference function
  cannot be injective.

  Proof (the root-of-unity pole argument, the discrete twin of THM-594(C)):
  let `D` be the largest difference and `ζ` a primitive `D`-th root of unity.
  Each class sums `ζ^x` to `ζ^a · Σ_{t < n/d} (ζ^d)^t`, which vanishes unless
  `D ∣ d` (forcing `d = D` by maximality), where it equals `ζ^a · (n/D) ≠ 0`.
  The total `Σ_{x<n} ζ^x = 0`, so the `d = D` classes must cancel among
  themselves — impossible for a single class.  Hence two classes share `D`.

  Continuous twin: THM-594 Part C (no finite distinct-speed danger-arc system
  tiles the circle; a divisor-minimal frequency carries `sin(2πr)/π ≠ 0`).
-/
import Mathlib.Tactic
import Mathlib.RingTheory.RootsOfUnity.Complex

namespace TournamentH7.PolygonMirskyNewman

open Finset Complex

/-- The color class with difference `d` and offset `a`, inside `range n`. -/
def cls (n d a : ℕ) : Finset ℕ := (range n).filter (fun x => x % d = a)

/-- Re-indexing: the class `{x < n : x % d = a}` is the image of
`t ↦ a + t * d` on `range (n / d)`, when `d ∣ n` and `a < d`. -/
lemma cls_eq_image (n d a : ℕ) (hd : d ∣ n) (ha : a < d) :
    cls n d a = (range (n / d)).image (fun t => a + t * d) := by
  have hd0 : 0 < d := lt_of_le_of_lt (Nat.zero_le a) ha
  ext x
  simp only [cls, mem_filter, mem_range, mem_image]
  constructor
  · rintro ⟨hxn, hxa⟩
    have hx := Nat.div_add_mod x d
    rw [hxa] at hx
    refine ⟨x / d, Nat.div_lt_div_of_lt_of_dvd hd hxn, ?_⟩
    rw [Nat.mul_comm]
    omega
  · rintro ⟨t, ht, rfl⟩
    have hnd : n / d * d = n := Nat.div_mul_cancel hd
    refine ⟨?_, by simp [Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt ha]⟩
    have htd : t + 1 ≤ n / d := ht
    calc a + t * d < d + t * d := by omega
      _ = (t + 1) * d := by ring
      _ ≤ n / d * d := Nat.mul_le_mul_right d htd
      _ = n := hnd

/-- The class character sum: `Σ_{x ∈ cls} ζ^x = ζ^a · Σ_{t < n/d} (ζ^d)^t`. -/
lemma cls_charSum (n d a : ℕ) (hd : d ∣ n) (ha : a < d) (ζ : ℂ) :
    ∑ x ∈ cls n d a, ζ ^ x = ζ ^ a * ∑ t ∈ range (n / d), (ζ ^ d) ^ t := by
  have hd0 : 0 < d := lt_of_le_of_lt (Nat.zero_le a) ha
  rw [cls_eq_image n d a hd ha, Finset.sum_image ?hinj, Finset.mul_sum]
  · exact Finset.sum_congr rfl fun t _ => by rw [pow_add, mul_comm t d, pow_mul]
  case hinj =>
    intro s _ t _ h
    simp only at h
    exact Nat.eq_of_mul_eq_mul_right hd0 (by omega)

/-- **The regular-polygon coloring theorem** (discrete Mirsky–Newman).
Let `n > 0`, and let a `k`-coloring of `{0,…,n−1}` be given by a class map
`c : ℕ → Fin k` whose fiber over each `i` is the regular-(n / d i)-gon class
`cls n (d i) (a i)`.  If `k ≥ 2`, two classes have the same difference —
congruent polygons. -/
theorem two_congruent_classes (n k : ℕ) (d a : Fin k → ℕ) (c : ℕ → Fin k)
    (hn : 0 < n) (hk : 2 ≤ k)
    (hdvd : ∀ i, d i ∣ n) (hoff : ∀ i, a i < d i)
    (hfiber : ∀ i, (range n).filter (fun x => c x = i) = cls n (d i) (a i)) :
    ∃ i j : Fin k, i ≠ j ∧ d i = d j := by
  by_contra hcon
  push_neg at hcon
  have hinj : ∀ i j : Fin k, d i = d j → i = j := by
    intro i j hij
    by_contra hne
    exact hcon i j hne hij
  have hd0 : ∀ i, 0 < d i := fun i => lt_of_le_of_lt (Nat.zero_le (a i)) (hoff i)
  -- the largest difference D at index i₀
  have hne : (Finset.univ : Finset (Fin k)).Nonempty := ⟨⟨0, by omega⟩, mem_univ _⟩
  obtain ⟨i₀, -, hmax⟩ := Finset.exists_max_image Finset.univ d hne
  set D := d i₀ with hDdef
  -- a second index, to rule out D = 1
  have hj : ∃ j : Fin k, j ≠ i₀ := by
    by_cases h0 : i₀ = ⟨0, by omega⟩
    · exact ⟨⟨1, by omega⟩, by simp [h0, Fin.ext_iff]⟩
    · exact ⟨⟨0, by omega⟩, fun h => h0 h.symm⟩
  have hD2 : 2 ≤ D := by
    by_contra hD1
    have hDone : D = 1 := by have := hd0 i₀; omega
    obtain ⟨j, hjne⟩ := hj
    have h1 : d j ≤ D := hmax j (mem_univ j)
    have h2 : d j = 1 := by have := hd0 j; omega
    exact hjne (hinj j i₀ (by rw [h2, ← hDone]))
  have hDn : D ∣ n := hdvd i₀
  have hD0 : D ≠ 0 := by omega
  -- the primitive D-th root of unity
  have hζ : IsPrimitiveRoot (Complex.exp (2 * Real.pi * Complex.I / D)) D :=
    Complex.isPrimitiveRoot_exp D hD0
  set ζ := Complex.exp (2 * Real.pi * Complex.I / D) with hζdef
  have hζne1 : ζ ≠ 1 := hζ.ne_one (by omega)
  have hζn : ζ ^ n = 1 := (hζ.pow_eq_one_iff_dvd n).mpr hDn
  -- total sum over all vertices is zero
  have htotal : ∑ x ∈ range n, ζ ^ x = 0 := by
    rw [geom_sum_eq hζne1, hζn]
    simp
  -- fiber decomposition of the total sum
  have hsplit : ∑ x ∈ range n, ζ ^ x
      = ∑ i : Fin k, ∑ x ∈ (range n).filter (fun x => c x = i), ζ ^ x :=
    (Finset.sum_fiberwise_of_maps_to (fun x _ => mem_univ (c x)) _).symm
  -- each fiber's value via the dichotomy
  have hfval : ∀ i : Fin k,
      ∑ x ∈ (range n).filter (fun x => c x = i), ζ ^ x
        = if D ∣ d i then ζ ^ (a i) * (n / d i : ℕ) else 0 := by
    intro i
    rw [hfiber i, cls_charSum n (d i) (a i) (hdvd i) (hoff i) ζ]
    by_cases hDd : D ∣ d i
    · rw [if_pos hDd, (hζ.pow_eq_one_iff_dvd (d i)).mpr hDd]
      simp
    · rw [if_neg hDd]
      have hζd : ζ ^ (d i) ≠ 1 := fun h => hDd (hζ.dvd_of_pow_eq_one (d i) h)
      rw [geom_sum_eq hζd, ← pow_mul, Nat.mul_div_cancel' (hdvd i), hζn]
      simp
  -- D ∣ d i forces i = i₀ (maximality + injectivity)
  have hDdi : ∀ i : Fin k, D ∣ d i → i = i₀ := by
    intro i hDd
    have h1 : d i ≤ D := hmax i (mem_univ i)
    have h2 : D ≤ d i := Nat.le_of_dvd (hd0 i) hDd
    exact hinj i i₀ (by omega)
  -- the total collapses to the single i₀ term
  have hsingle : ∑ i : Fin k, ∑ x ∈ (range n).filter (fun x => c x = i), ζ ^ x
      = ζ ^ (a i₀) * (n / D : ℕ) := by
    rw [Finset.sum_eq_single i₀]
    · rw [hfval i₀, if_pos (dvd_refl D)]
    · intro j _ hjne
      rw [hfval j]
      rw [if_neg (fun hDd => hjne (hDdi j hDd))]
    · intro h
      exact absurd (mem_univ i₀) h
  -- contradiction: a nonzero single term equals zero
  have hzero : ζ ^ (a i₀) * (n / D : ℕ) = 0 := by
    rw [← hsingle, ← hsplit, htotal]
  have hζa : ζ ^ (a i₀) ≠ 0 := pow_ne_zero _ (hζ.ne_zero hD0)
  have hnD : (0 : ℕ) < n / D := Nat.div_pos (Nat.le_of_dvd hn hDn) (by omega)
  have : ((n / D : ℕ) : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  exact this (by
    rcases mul_eq_zero.mp hzero with h | h
    · exact absurd h hζa
    · exact h)

end TournamentH7.PolygonMirskyNewman
