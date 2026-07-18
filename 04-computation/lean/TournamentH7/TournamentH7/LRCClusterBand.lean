/- LRCClusterBand.lean -- kind-pasteur-2026-07-18-S128 (cont.57), THM-1032 / THM-1018(II).

   THE BAND CERTIFICATE, and the explicit modulus that closes its existence step.

   THM-1018(II) (the band lemma): fix a modulus `q` and let each speed reduce to
   `±eᵢ (mod q)`.  If every `eᵢ · a` lands in the band `[q/14, 13q/14]` then
   `t = a/q` is 14-lonely -- each effective speed sits in the safe band with no
   wraparound, and `‖·‖` is even so the sign of the residue is irrelevant.

   THM-1032 (the existence step): for a clustered family `V = P ∪ {v₁ … v_r}` with
   core aspect ratio `max P ≤ 12 · min P` and killer spread `≤ max P − min P`, the
   modulus `q = v₁ + max P` puts EVERY killer residue back inside `[min P, max P]`
   (because `v_i < q` and `q − v_i = max P − δᵢ`), so the band hypothesis is
   inherited from the core rather than assumed.  No divisor count is needed.

   This module proves the band certificate in the reusable `±e (mod q)` form, plus
   the two speed shapes the construction produces (`e` itself and `q − e`), and
   instantiates it on a critical covering family.

   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

section Band

/-- **The band inequality.**  If `e · a` lies in `[q/14, 13q/14]` then `e · (a/q)`
is at distance `≥ 1/14` from every integer.  This is the arithmetic heart of
THM-1018(II): the point never wraps, so both the "too small" and "too large"
failure modes are excluded by one interval condition. -/
theorem band_dist (q a e : ℕ) (hq : 0 < q)
    (hlo : q ≤ 14 * (e * a)) (hhi : 14 * (e * a) ≤ 13 * q) (m : ℤ) :
    (1 : ℝ) / 14 ≤ |(e : ℝ) * ((a : ℝ) / q) - m| := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hlo' : (q : ℝ) ≤ 14 * ((e : ℝ) * (a : ℝ)) := by exact_mod_cast hlo
  have hhi' : 14 * ((e : ℝ) * (a : ℝ)) ≤ 13 * (q : ℝ) := by exact_mod_cast hhi
  have hxeq : (e : ℝ) * ((a : ℝ) / q) = ((e : ℝ) * (a : ℝ)) / q := by ring
  rw [hxeq]
  set y : ℝ := (e : ℝ) * (a : ℝ) with hy
  have h1 : (1 : ℝ) / 14 ≤ y / q := by
    rw [le_div_iff₀ hqR]; linarith
  have h2 : y / q ≤ 13 / 14 := by
    rw [div_le_iff₀ hqR]; linarith
  by_cases hm : m ≤ 0
  · have hmR : (m : ℝ) ≤ 0 := by exact_mod_cast hm
    have : (1 : ℝ) / 14 ≤ y / q - m := by linarith
    exact this.trans (le_abs_self _)
  · have hm1 : (1 : ℤ) ≤ m := by omega
    have hmR : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm1
    have : (1 : ℝ) / 14 ≤ -(y / q - m) := by linarith
    exact this.trans (neg_le_abs _)

/-- **THE BAND CERTIFICATE (THM-1018 II).**  If every speed is `≡ ±eᵢ (mod q)` and
every `eᵢ · a` lies in the band `[q/14, 13q/14]`, then `a/q` is a 14-lonely time. -/
theorem lonely_of_band {ι : Type*} (q a : ℕ) (hq : 0 < q) (v : ι → ℤ) (e : ι → ℕ)
    (hlo : ∀ i, q ≤ 14 * (e i * a)) (hhi : ∀ i, 14 * (e i * a) ≤ 13 * q)
    (hres : ∀ i, ∃ k : ℤ, v i = k * q + (e i : ℤ) ∨ v i = k * q - (e i : ℤ)) :
    Lonely 14 v ((a : ℝ) / q) := by
  intro i m
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hqne : (q : ℝ) ≠ 0 := ne_of_gt hqR
  obtain ⟨k, hk⟩ := hres i
  have hnum : ((14 : ℕ) : ℝ) = 14 := by norm_num
  simp only [hnum]
  rcases hk with hk | hk
  · -- v i = k q + e i :  v i · (a/q) = k a + e i · (a/q)
    have hval : ((v i : ℤ) : ℝ) * ((a : ℝ) / q)
        = (e i : ℝ) * ((a : ℝ) / q) + (k : ℝ) * (a : ℝ) := by
      rw [hk]; push_cast; field_simp; try ring
    rw [hval]
    have : (e : ι → ℕ) i * ((a : ℝ) / q) + (k : ℝ) * (a : ℝ) - (m : ℝ)
        = (e i : ℝ) * ((a : ℝ) / q) - ((m - k * a : ℤ) : ℝ) := by push_cast; ring
    rw [this]
    exact band_dist q a (e i) hq (hlo i) (hhi i) _
  · -- v i = k q - e i :  v i · (a/q) = k a - e i · (a/q)
    have hval : ((v i : ℤ) : ℝ) * ((a : ℝ) / q)
        = (k : ℝ) * (a : ℝ) - (e i : ℝ) * ((a : ℝ) / q) := by
      rw [hk]; push_cast; field_simp; try ring
    rw [hval]
    have : (k : ℝ) * (a : ℝ) - (e i : ℝ) * ((a : ℝ) / q) - (m : ℝ)
        = -((e i : ℝ) * ((a : ℝ) / q) - ((k * a - m : ℤ) : ℝ)) := by push_cast; ring
    rw [this, abs_neg]
    exact band_dist q a (e i) hq (hlo i) (hhi i) _

end Band

section Construction

/-- A core speed is its own residue: `p = 0·q + p`. -/
theorem res_self (q : ℕ) (p : ℕ) : ∃ k : ℤ, (p : ℤ) = k * q + (p : ℤ) :=
  ⟨0, by ring⟩

/-- **The killer residue identity (THM-1032 b).**  With `q = v₁ + M` and a killer
`v = v₁ + δ` where `δ ≤ M`, we have `v = 1·q − (M − δ)`: the killer, however large,
reduces to the small residue `M − δ`, which lies in the core's own window whenever
`δ ≤ M − μ`.  This is what makes the band hypothesis inherited rather than assumed. -/
theorem killer_res (v₁ M δ : ℕ) (hδ : δ ≤ M) :
    ((v₁ + δ : ℕ) : ℤ) = 1 * ((v₁ + M : ℕ) : ℤ) - ((M - δ : ℕ) : ℤ) := by
  have : (((M - δ : ℕ)) : ℤ) = (M : ℤ) - (δ : ℤ) := by
    exact_mod_cast Int.ofNat_sub hδ
  push_cast [this]
  ring

end Construction

section CGeneralized

/-! ### THM-1041: letting the offset `c` float

`THM-1032` used `q = v₁ + M`.  For a general `q = v₁ + c` the killer residues are
`|δᵢ − c|`: below `c` the killer sits under `q` and `q − vᵢ = c − δᵢ`, above `c` it sits
over `q` and `vᵢ − q = δᵢ − c`.  Both shapes are recorded here, and then the midpoint
choice `c = ⌈D/2⌉` is shown to land every offset of a two-killer cluster inside `[μ, M]`
exactly when `2μ ≤ D ≤ 2M` — which is what extends the explicit reach to `D ≤ 2M`. -/

/-- Killer below the offset: `v₁ + δ = 1·(v₁ + c) − (c − δ)` for `δ ≤ c`. -/
theorem killer_res_below (v₁ c δ : ℕ) (hδ : δ ≤ c) :
    ((v₁ + δ : ℕ) : ℤ) = 1 * ((v₁ + c : ℕ) : ℤ) - ((c - δ : ℕ) : ℤ) := by
  have h : (((c - δ : ℕ)) : ℤ) = (c : ℤ) - (δ : ℤ) := by exact_mod_cast Int.ofNat_sub hδ
  push_cast [h]; ring

/-- Killer above the offset: `v₁ + δ = 1·(v₁ + c) + (δ − c)` for `c ≤ δ`. -/
theorem killer_res_above (v₁ c δ : ℕ) (hδ : c ≤ δ) :
    ((v₁ + δ : ℕ) : ℤ) = 1 * ((v₁ + c : ℕ) : ℤ) + ((δ - c : ℕ) : ℤ) := by
  have h : (((δ - c : ℕ)) : ℤ) = (δ : ℤ) - (c : ℤ) := by exact_mod_cast Int.ofNat_sub hδ
  push_cast [h]; ring

/-- **THE MIDPOINT CHOICE (THM-1041 II).**  For a two-killer cluster of spread `D` with
`2μ ≤ D ≤ 2M`, the offset `c = ⌈D/2⌉ = (D+1)/2` puts BOTH killer residues — `c` itself and
`D − c` — inside the core window `[μ, M]`.  This is what makes the reach `D ≤ 2M`. -/
theorem midpoint_in_window (μ M D : ℕ) (hlo : 2 * μ ≤ D) (hhi : D ≤ 2 * M) :
    μ ≤ (D + 1) / 2 ∧ (D + 1) / 2 ≤ M ∧ μ ≤ D - (D + 1) / 2 ∧ D - (D + 1) / 2 ≤ M := by
  omega

/-- **Sharpness (THM-1041 III).**  At `D = 2M + 1` no offset works: any `c ≤ M` leaves the
far killer at distance `D − c ≥ M + 1`, outside the window. -/
theorem no_offset_at_two_M_succ (μ M c : ℕ) (hμ : 1 ≤ μ) (hc : μ ≤ c) (hcM : c ≤ M) :
    M + 1 ≤ (2 * M + 1) - c := by
  omega

end CGeneralized

section Demo

/-- A worked instance of THM-1032 on a critical covering family of THM-1011(VII):
core `{2,…,12}`, killers `v₁ = 168`, `v₂ = 169`, so `M = 12`, `μ = 2`,
`q = v₁ + M = 180`, `a = ⌈q/(14μ)⌉ = 7`.  The killers reduce to `e₁ = 12`,
`e₂ = 11`, and every effective speed `e ∈ {2,…,12}` satisfies the band
`180 ≤ 14·(7e)` and `14·(7e) ≤ 13·180`.  Checked here for the extreme
residues, which are the binding ones. -/
example : (180 : ℕ) ≤ 14 * (2 * 7) ∧ 14 * (12 * 7) ≤ 13 * 180 := by norm_num

/-- The band certificate applied to the single speed `168` at `t = 7/180`:
`168 = 1·180 − 12`, and `12·7 = 84` lies in `[180/14, 13·180/14]`. -/
example (m : ℤ) : (1 : ℝ) / 14 ≤ |(12 : ℝ) * ((7 : ℝ) / 180) - m| :=
  band_dist 180 7 12 (by norm_num) (by norm_num) (by norm_num) m

end Demo

end LonelyRunner

/-! ### Axiom audit -- kernel purity check (expect the standard trio only). -/
#print axioms LonelyRunner.band_dist
#print axioms LonelyRunner.lonely_of_band
#print axioms LonelyRunner.killer_res
#print axioms LonelyRunner.res_self
#print axioms LonelyRunner.killer_res_below
#print axioms LonelyRunner.killer_res_above
#print axioms LonelyRunner.midpoint_in_window
#print axioms LonelyRunner.no_offset_at_two_M_succ
