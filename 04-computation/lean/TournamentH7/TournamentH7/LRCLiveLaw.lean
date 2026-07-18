/-
  TournamentH7.LRCLiveLaw — THE LIVE LAW OF THE CANONICAL FAMILY
  (death-star-2026-07-17-S55, HYP-7280).

  **liveCount((1,…,13), q) = 6 if 14 ∣ q, else 0** — with THM-987's deep
  count this closes the ENTIRE canonical census in closed form.

  * `resonant_live` — at `q = 14m`, `p = o·m` with `o` a unit mod 14 is live:
    the mod-scaling `(x·m) % (14m) = (x % 14)·m` makes every band check the
    single fact `14 ∤ (i+1)·o` (coprimality + `i+1 ≤ 13`).
  * `live_residues_distinct` / `live_gap` — a live multiplier makes the 14
    residues `{i·p mod q : i ∈ [0,13]}` pairwise ≥ `⌈q/14⌉` apart — the
    canonical family is DIFFERENCE-CLOSED (every difference of indices in
    [0,13] is a family speed), which is exactly what the rigidity consumes.
  * `live_forces_dvd` — the block injection: 14 elements pairwise ≥ c apart
    inject into `(max−min)/c + 1` blocks, giving `max − min ≥ 13c`; the wrap
    gap (`min + q − max ≥ c`, again a family speed) closes `q ≥ 14·⌈q/14⌉`,
    which forces `14 ∣ q`.  NO sorting machinery — the block map
    `x ↦ (x − min)/c` is injective outright.
  * `liveCount_eq_zero_of_not_dvd` / `liveCount_ge_six_of_dvd` — the census-
    critical directions of the law.  (The exact-6 upper bound at resonance —
    the m-net extraction with the δ-chain — is the named next piece; recon
    confirms = 6 with the exact set {o·m : o a unit}.)

  Recon (`livelaw_recon_deathstar_S55.out`): count AND the exact live set
  verified for every q ∈ [2, 3000).

  Kernel-pure: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCTwoCircleII

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The ℕ-residue bridge for `inBand` on the canonical family. -/
theorem inBand_canonical_iff (q p : ℕ) (i : Fin 13) :
    inBand canonical q p i ↔
      (q ≤ 14 * ((((i : ℕ) + 1) * p) % q) ∧
        14 * ((((i : ℕ) + 1) * p) % q) ≤ 13 * q) := by
  unfold inBand canonical
  have hcast : ((i : ℤ) + 1) * (p : ℤ) % (q : ℤ)
      = ((((i : ℕ) + 1) * p % q : ℕ) : ℤ) := by
    push_cast
    ring_nf
  rw [hcast]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩
  · rintro ⟨h1, h2⟩
    exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩

/-- **The resonant direction**: at `q = 14m`, every unit multiple `o·m` is
live. -/
theorem resonant_live (m o : ℕ) (hm : 0 < m) (ho : o < 14)
    (hcop : Nat.gcd o 14 = 1) :
    bandCount canonical (14 * m) (o * m) = 0 := by
  unfold bandCount
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro i _
  rw [not_not, inBand_canonical_iff]
  -- the residue scales: ((i+1)·o·m) % (14m) = (((i+1)·o) % 14)·m
  have hscale : (((i : ℕ) + 1) * (o * m)) % (14 * m)
      = ((((i : ℕ) + 1) * o) % 14) * m := by
    have h1 : ((i : ℕ) + 1) * (o * m) = (((i : ℕ) + 1) * o) * m := by ring
    rw [h1, Nat.mul_mod_mul_right]
  rw [hscale]
  -- the scaled residue is in [m, 13m] because (i+1)·o is never 0 mod 14
  have hne : (((i : ℕ) + 1) * o) % 14 ≠ 0 := by
    intro h0
    have hdvd : 14 ∣ ((i : ℕ) + 1) * o := Nat.dvd_of_mod_eq_zero h0
    have hcop' : Nat.Coprime 14 o := by
      rw [Nat.coprime_comm]
      exact hcop
    have h14 : 14 ∣ (i : ℕ) + 1 := hcop'.dvd_of_dvd_mul_right hdvd
    have := Nat.le_of_dvd (by omega) h14
    omega
  have hlt : (((i : ℕ) + 1) * o) % 14 < 14 := Nat.mod_lt _ (by omega)
  constructor
  · nlinarith [Nat.one_le_iff_ne_zero.mpr hne, hm]
  · nlinarith [hlt, hm]

/-- A live multiplier's safety, in ℕ-residue form, for every family speed. -/
theorem live_safe (q p : ℕ) (hlive : bandCount canonical q p = 0)
    (v : ℕ) (hv1 : 1 ≤ v) (hv13 : v ≤ 13) :
    q ≤ 14 * ((v * p) % q) ∧ 14 * ((v * p) % q) ≤ 13 * q := by
  have hib : inBand canonical q p ⟨v - 1, by omega⟩ := by
    by_contra hcon
    unfold bandCount at hlive
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff] at hlive
    exact hlive (Finset.mem_univ (⟨v - 1, by omega⟩ : Fin 13)) hcon
  have h := (inBand_canonical_iff q p ⟨v - 1, by omega⟩).mp hib
  have hval : ((⟨v - 1, by omega⟩ : Fin 13) : ℕ) + 1 = v := by
    show (v - 1) + 1 = v
    omega
  rwa [hval] at h

/-- The symmetric pair-gap law: distinct indices `k, l ≤ 13` of a live
multiplier have residues at distance ≥ `⌈q/14⌉` AND co-distance ≥ `⌈q/14⌉`
— both edges of the safe band, through the difference speed `l − k`. -/
theorem live_gap (q p : ℕ) (hq : 0 < q)
    (hlive : bandCount canonical q p = 0)
    (k l : ℕ) (hkl : k < l) (hl13 : l ≤ 13) :
    ((q + 13) / 14 ≤ (l * p) % q - (k * p) % q ∨
     (q + 13) / 14 ≤ (k * p) % q - (l * p) % q) ∧
    (l * p) % q ≠ (k * p) % q ∧
    (q + 13) / 14 ≤ q - ((l * p) % q - (k * p) % q) ∧
    (q + 13) / 14 ≤ q - ((k * p) % q - (l * p) % q) := by
  have hd1 : 1 ≤ l - k := by omega
  have hd13 : l - k ≤ 13 := by omega
  obtain ⟨hs1, hs2⟩ := live_safe q p hlive (l - k) hd1 hd13
  set r : ℕ := ((l - k) * p) % q with hr_def
  have hr_lt : r < q := Nat.mod_lt _ hq
  have hfk_lt : (k * p) % q < q := Nat.mod_lt _ hq
  have hfl_lt : (l * p) % q < q := Nat.mod_lt _ hq
  -- the split (l−k)p + kp = lp gives f l = (r + f k) % q
  have hsplit : (l - k) * p + k * p = l * p := by
    rw [Nat.sub_mul]
    have : k * p ≤ l * p := Nat.mul_le_mul_right p (le_of_lt hkl)
    omega
  have heq : (l * p) % q = (r + (k * p) % q) % q := by
    conv_lhs => rw [← hsplit]
    rw [Nat.add_mod, ← hr_def]
  rcases Nat.lt_or_ge (r + (k * p) % q) q with hcase | hcase
  · -- no wraparound: f l = r + f k
    have hfl : (l * p) % q = r + (k * p) % q := by
      rw [heq]
      exact Nat.mod_eq_of_lt hcase
    refine ⟨Or.inl (by omega), by omega, by omega, by omega⟩
  · -- wraparound: f l = r + f k − q
    have hlt2 : r + (k * p) % q - q < q := by omega
    have hfl : (l * p) % q = r + (k * p) % q - q := by
      rw [heq, Nat.mod_eq_sub_mod hcase]
      exact Nat.mod_eq_of_lt hlt2
    refine ⟨Or.inr (by omega), by omega, by omega, by omega⟩

open Classical in
/-- **The rigidity**: a live multiplier forces `14 ∣ q` — the block injection
(pairwise gaps ≥ c pack the 14 residues into `M/c + 1` blocks, so the spread
is ≥ 13c) plus the wrap gap close `q ≥ 14·⌈q/14⌉`. -/
theorem live_forces_dvd (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hlive : bandCount canonical q p = 0) :
    14 ∣ q := by
  set c : ℕ := (q + 13) / 14 with hc_def
  have hc1 : 1 ≤ c := by omega
  set f : ℕ → ℕ := fun k => (k * p) % q with hf_def
  have hf0 : f 0 = 0 := by simp [hf_def]
  obtain ⟨kM, hkM_mem, hkM⟩ := Finset.exists_max_image (Finset.range 14) f
    ⟨0, Finset.mem_range.mpr (by omega)⟩
  set M : ℕ := f kM with hM_def
  -- distinct indices land in distinct c-blocks
  have hgap_sep : ∀ k l : ℕ, k < l → l ≤ 13 → f k / c ≠ f l / c := by
    intro k l hkl hl13 hEq
    obtain ⟨hor, -, -, -⟩ := live_gap q p hq hlive k l hkl hl13
    have hdk := Nat.div_add_mod (f k) c
    have hdl := Nat.div_add_mod (f l) c
    have hrk : f k % c < c := Nat.mod_lt _ (by omega)
    have hrl : f l % c < c := Nat.mod_lt _ (by omega)
    rw [← hEq] at hdl
    -- the common block term cancels; the residues are within c of each other
    have hcancel : f k + f l % c = f l + f k % c := by
      linarith [hdk, hdl]
    have hfk_eq : f k = (k * p) % q := rfl
    have hfl_eq : f l = (l * p) % q := rfl
    rcases hor with h | h <;> omega
  have hinj : Set.InjOn (fun k => f k / c) (Finset.range 14) := by
    intro k hk l hl hEq
    rw [Finset.coe_range, Set.mem_Iio] at hk hl
    rcases lt_trichotomy k l with hkl | rfl | hlk
    · exact absurd hEq (hgap_sep k l hkl (by omega))
    · rfl
    · exact absurd hEq.symm (hgap_sep l k hlk (by omega))
  -- the blocks live inside range (M/c + 1)
  have hmaps : ∀ k ∈ Finset.range 14, f k / c ∈ Finset.range (M / c + 1) := by
    intro k hk
    rw [Finset.mem_range]
    have := Nat.div_le_div_right (c := c) (hkM k hk)
    omega
  have hcard := Finset.card_le_card_of_injOn (fun k => f k / c) hmaps
    (by
      intro a ha b hb hab
      exact hinj (by simpa using ha) (by simpa using hb) hab)
  rw [Finset.card_range, Finset.card_range] at hcard
  -- spread ≥ 13c
  have hM13 : 13 * c ≤ M := by
    have h13 : 13 ≤ M / c := by omega
    exact (Nat.le_div_iff_mul_le (by omega)).mp h13
  -- the wrap gap: q − M ≥ c
  have hkM13 : kM ≤ 13 := by
    rw [Finset.mem_range] at hkM_mem
    omega
  have hkM0 : 0 < kM := by
    by_contra h0
    have hkm0 : kM = 0 := by omega
    have : M = 0 := by rw [hM_def, hkm0, hf0]
    omega
  obtain ⟨-, -, hco, -⟩ := live_gap q p hq hlive 0 kM hkM0 hkM13
  -- hco : c ≤ q − (f kM − f 0) = q − M
  have hco' : c ≤ q - M := by
    have h0 : (0 * p) % q = 0 := by simp
    rw [h0, Nat.sub_zero] at hco
    exact hco
  -- q ≥ M + c ≥ 14c = 14·⌈q/14⌉ forces 14 ∣ q
  omega

open Classical in
/-- **The live law, zero direction**: off resonance there are NO live
multipliers. -/
theorem liveCount_eq_zero_of_not_dvd (q : ℕ) (hq : 0 < q) (hndvd : ¬ 14 ∣ q) :
    liveCount canonical q = 0 := by
  unfold liveCount
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro p hp
  rw [Finset.mem_Ioo] at hp
  intro hlive
  exact hndvd (live_forces_dvd q p hq hp.1 hp.2 hlive)

open Classical in
/-- **The live law, resonant direction**: at `q = 14m` the six unit multiples
are live — `liveCount ≥ 6`. -/
theorem liveCount_ge_six_of_dvd (m : ℕ) (hm : 0 < m) :
    6 ≤ liveCount canonical (14 * m) := by
  unfold liveCount
  have hsub : (({1, 3, 5, 9, 11, 13} : Finset ℕ).image (fun o => o * m))
      ⊆ (Finset.Ioo 0 (14 * m)).filter
          fun p => bandCount canonical (14 * m) p = 0 := by
    intro x hx
    rw [Finset.mem_image] at hx
    obtain ⟨o, ho_mem, rfl⟩ := hx
    rw [Finset.mem_filter, Finset.mem_Ioo]
    fin_cases ho_mem <;>
      exact ⟨⟨by omega, by omega⟩,
        resonant_live m _ hm (by norm_num) (by decide)⟩
  have hcard : (({1, 3, 5, 9, 11, 13} : Finset ℕ).image (fun o => o * m)).card
      = 6 := by
    rw [Finset.card_image_of_injOn]
    · decide
    · intro a _ b _ hab
      exact Nat.eq_of_mul_eq_mul_right hm hab
  calc 6 = _ := hcard.symm
    _ ≤ _ := Finset.card_le_card hsub

/-! ## Axiom audit -/
#print axioms inBand_canonical_iff
#print axioms resonant_live
#print axioms live_safe
#print axioms live_gap
#print axioms live_forces_dvd
#print axioms liveCount_eq_zero_of_not_dvd
#print axioms liveCount_ge_six_of_dvd

end LRC14Concrete
end LonelyRunner
