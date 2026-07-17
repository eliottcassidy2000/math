/-
  TournamentH7.LRCRungLock — THE RUNG LOCK AND THE INTEGER COMPONENT COLLAPSE
  (death-star-2026-07-17-S46, HYP-7225; the dense-core discrete side).

  A runner FAILS the safe band at multiplier `p`, modulus `q` iff some integer
  witness `w` has `14·|v·p − w·q| < q` (the `¬inBand` normal form, THM-947).
  Three theorems sharpen the adaptive-q component count (THM-956) on exact-
  ratio strata — the dense core's single ≥7-comparable blocks:

  * `rung_lock` — two failing runners with EXACT integer speed ratio
    `v_j = m·v_i`, `1 ≤ m ≤ 13`, have EXACTLY locked witnesses: `w_j = m·w_i`.
    Zero choices — the strict chain `14|w_j − m·w_i|·q < q + m·q = (m+1)·q ≤
    14·q` survives THROUGH m = 13 because both band bounds are strict.  The
    threshold is EXACTLY the lonely-runner denominator: lock ⟺ m ≤ 13 = n;
    `rung_lock_sharp_at_14` exhibits the break at m = 14 in-kernel
    (v = 1 → 14, q = 15, p = 1: witnesses 0 vs 1).  On any block with
    consecutive integer ratios ≤ 13 the bottom witness determines ALL
    witnesses — sharper than the conjectured two-choice tower bound.
    HONEST SCOPE (recon S46): the companion divisor-collapse heuristic is
    REFUTED on compressed block towers (deep instants are COMMON there, K
    huge, no scale divisibility); the salvage is the CARRIER-COUNT pigeonhole
    (worst 134 carriers vs 299 window moduli on 60 block families — deep-free
    q exists in every recon window).
  * `rung_lock_chain3` — composition: three runners with ratios `m₁, m₂`
    lock multiplicatively, `w₃ = m₂·m₁·w₁`.
  * `ray_fails_above` — an exact ray (`q ∣ v_i·p`) forces every integer
    multiple `v_j = k·v_i` to FAIL the band: upward band propagation on
    towers (deep sets over rays come for free).
  * `carrier_instant_eq_of_shared_witness` — the ALL-INTEGER Farey collapse:
    two moduli `q₁, q₂ ≤ D` whose multipliers share one witness `w` for one
    runner `v_i` with `D² ≤ 7·v_i` have THE SAME instant
    (`p₁·q₂ = p₂·q₁`) — the exact identity
    `v_i·(p₁q₂ − p₂q₁) = (v_i p₁ − w q₁)q₂ − (v_i p₂ − w q₂)q₁`
    plus `1 ≤ |Δ|` gives `14 v_i ≤ 2 q₁ q₂ ≤ 2 D²`, contradiction.
    Composed with `rung_lock`, the component map of THM-956's pigeonhole is
    CONCRETE on exact-ratio strata: component = the bottom failing witness.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDeepCount

namespace LonelyRunner
namespace LRC14Concrete

/-- **THE RUNG LOCK**: failing runners with exact integer speed ratio
`m ≤ 12` lock witnesses exactly. -/
theorem rung_lock (vi vj wi wj m : ℤ) (q p : ℕ)
    (hm1 : 1 ≤ m) (hm13 : m ≤ 13) (hratio : vj = m * vi)
    (hi : 14 * |vi * (p : ℤ) - wi * q| < q)
    (hj : 14 * |vj * (p : ℤ) - wj * q| < q) :
    wj = m * wi := by
  have hmpos : (0 : ℤ) < m := by linarith
  -- scale the low-rung bound by m
  have hj' : 14 * |vj * (p : ℤ) - m * wi * q| < m * q := by
    have hrw : vj * (p : ℤ) - m * wi * q = m * (vi * p - wi * q) := by
      rw [hratio]; ring
    rw [hrw, abs_mul, abs_of_pos hmpos]
    calc 14 * (m * |vi * (p : ℤ) - wi * q|)
        = m * (14 * |vi * (p : ℤ) - wi * q|) := by ring
      _ < m * q := by
          exact mul_lt_mul_of_pos_left hi hmpos
  -- triangle through vj·p
  have htri : |wj * (q : ℤ) - m * wi * q|
      ≤ |vj * (p : ℤ) - wj * q| + |vj * (p : ℤ) - m * wi * q| := by
    have h1 : wj * (q : ℤ) - m * wi * q
        = (vj * (p : ℤ) - m * wi * q) - (vj * (p : ℤ) - wj * q) := by ring
    rw [h1]
    calc |(vj * (p : ℤ) - m * wi * q) - (vj * (p : ℤ) - wj * q)|
        ≤ |vj * (p : ℤ) - m * wi * q| + |vj * (p : ℤ) - wj * q| := abs_sub _ _
      _ = |vj * (p : ℤ) - wj * q| + |vj * (p : ℤ) - m * wi * q| := by ring
  -- 14·|wj − m·wi|·q < (m+1)·q ≤ 13·q < 14·q  ⟹  |wj − m·wi| = 0
  have hfact : |wj * (q : ℤ) - m * wi * q| = |wj - m * wi| * q := by
    have h2 : wj * (q : ℤ) - m * wi * q = (wj - m * wi) * q := by ring
    rw [h2, abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))]
  have hlt : 14 * (|wj - m * wi| * (q : ℤ)) < 14 * q := by
    calc 14 * (|wj - m * wi| * (q : ℤ))
        = 14 * |wj * (q : ℤ) - m * wi * q| := by rw [hfact]
      _ ≤ 14 * |vj * (p : ℤ) - wj * q| + 14 * |vj * (p : ℤ) - m * wi * q| := by
          linarith [htri]
      _ < q + m * q := by linarith [hj, hj']
      _ ≤ 14 * q := by nlinarith [hm13, Int.natCast_nonneg q]
  set x : ℤ := wj - m * wi with hx
  have hxabs : 14 * (|x| * (q : ℤ)) < 14 * q := hlt
  have hx0 : x = 0 := by
    by_contra hne
    have h1 : 1 ≤ |x| := Int.one_le_abs hne
    have hq1 : 1 ≤ (q : ℤ) := by
      rcases Nat.eq_zero_or_pos q with rfl | hq
      · simp at hxabs
      · exact_mod_cast hq
    nlinarith [h1, hq1]
  linarith [hx]

/-- **Chain composition**: three runners with consecutive exact ratios lock
multiplicatively — the bottom witness determines the whole tower. -/
theorem rung_lock_chain3 (v₁ v₂ v₃ w₁ w₂ w₃ m₁ m₂ : ℤ) (q p : ℕ)
    (hm₁ : 1 ≤ m₁) (hm₁13 : m₁ ≤ 13) (hm₂ : 1 ≤ m₂) (hm₂13 : m₂ ≤ 13)
    (hr₁ : v₂ = m₁ * v₁) (hr₂ : v₃ = m₂ * v₂)
    (h₁ : 14 * |v₁ * (p : ℤ) - w₁ * q| < q)
    (h₂ : 14 * |v₂ * (p : ℤ) - w₂ * q| < q)
    (h₃ : 14 * |v₃ * (p : ℤ) - w₃ * q| < q) :
    w₃ = m₂ * m₁ * w₁ := by
  have hlock₁ : w₂ = m₁ * w₁ := rung_lock v₁ v₂ w₁ w₂ m₁ q p hm₁ hm₁13 hr₁ h₁ h₂
  have hlock₂ : w₃ = m₂ * w₂ := rung_lock v₂ v₃ w₂ w₃ m₂ q p hm₂ hm₂13 hr₂ h₂ h₃
  rw [hlock₂, hlock₁]; ring

/-- **Upward band propagation**: an exact ray (`q ∣ v_i·p`) forces every
integer multiple `v_j = k·v_i` to fail the safe band. -/
theorem ray_fails_above (v : Fin 13 → ℤ) (q p : ℕ) (hq : 0 < q)
    (i j : Fin 13) (k : ℤ) (hratio : v j = k * v i)
    (hray : (v i * (p : ℤ)) % (q : ℤ) = 0) :
    ¬ inBand v q p j := by
  have hdvd : (q : ℤ) ∣ v i * p := Int.dvd_of_emod_eq_zero hray
  have hdvd2 : (q : ℤ) ∣ v j * p := by
    rw [hratio]
    have : k * v i * (p : ℤ) = k * (v i * p) := by ring
    rw [this]
    exact Dvd.dvd.mul_left hdvd k
  have hmod : (v j * (p : ℤ)) % (q : ℤ) = 0 := Int.emod_eq_zero_of_dvd hdvd2
  intro hband
  have h1 := hband.1
  rw [hmod] at h1
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  linarith

/-- **THE ALL-INTEGER FAREY COLLAPSE**: two window moduli sharing one witness
for one runner with `D² ≤ 7·v_i` carry THE SAME instant. -/
theorem carrier_instant_eq_of_shared_witness (vi w : ℤ) (q₁ p₁ q₂ p₂ D : ℕ)
    (hvi : 0 < vi) (hq₁ : 0 < q₁) (hq₂ : 0 < q₂)
    (hD₁ : q₁ ≤ D) (hD₂ : q₂ ≤ D)
    (hsuper : (D : ℤ) * D ≤ 7 * vi)
    (h₁ : 14 * |vi * (p₁ : ℤ) - w * q₁| < q₁)
    (h₂ : 14 * |vi * (p₂ : ℤ) - w * q₂| < q₂) :
    (p₁ : ℤ) * q₂ = (p₂ : ℤ) * q₁ := by
  by_contra hne
  have hq₁Z : (0 : ℤ) < (q₁ : ℤ) := by exact_mod_cast hq₁
  have hq₂Z : (0 : ℤ) < (q₂ : ℤ) := by exact_mod_cast hq₂
  -- the exact cross-multiplication identity
  have hkey : vi * ((p₁ : ℤ) * q₂ - p₂ * q₁)
      = (vi * p₁ - w * q₁) * q₂ - (vi * p₂ - w * q₂) * q₁ := by ring
  have habs : |vi * ((p₁ : ℤ) * q₂ - p₂ * q₁)|
      ≤ |vi * (p₁ : ℤ) - w * q₁| * q₂ + |vi * (p₂ : ℤ) - w * q₂| * q₁ := by
    rw [hkey]
    calc |(vi * (p₁ : ℤ) - w * q₁) * q₂ - (vi * (p₂ : ℤ) - w * q₂) * q₁|
        ≤ |(vi * (p₁ : ℤ) - w * q₁) * q₂| + |(vi * (p₂ : ℤ) - w * q₂) * q₁| :=
          abs_sub _ _
      _ = |vi * (p₁ : ℤ) - w * q₁| * q₂ + |vi * (p₂ : ℤ) - w * q₂| * q₁ := by
          rw [abs_mul, abs_mul,
            abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q₂ : ℤ)),
            abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q₁ : ℤ))]
  -- scale the band bounds by the opposite modulus
  have hm1 : 14 * (|vi * (p₁ : ℤ) - w * q₁| * q₂) < q₁ * q₂ := by
    have := mul_lt_mul_of_pos_right h₁ hq₂Z
    linarith [this]
  have hm2 : 14 * (|vi * (p₂ : ℤ) - w * q₂| * q₁) < q₂ * q₁ := by
    have := mul_lt_mul_of_pos_right h₂ hq₁Z
    linarith [this]
  -- 1 ≤ |Δ| pushes 14·vi below 2·D²
  have hΔ : 1 ≤ |(p₁ : ℤ) * q₂ - p₂ * q₁| := Int.one_le_abs (sub_ne_zero.mpr hne)
  have hvineq : vi ≤ |vi * ((p₁ : ℤ) * q₂ - p₂ * q₁)| := by
    rw [abs_mul, abs_of_pos hvi]
    calc vi = vi * 1 := by ring
      _ ≤ vi * |(p₁ : ℤ) * q₂ - p₂ * q₁| :=
          mul_le_mul_of_nonneg_left hΔ hvi.le
  have hDD : (q₁ : ℤ) * q₂ ≤ (D : ℤ) * D := by
    have hd₁ : (q₁ : ℤ) ≤ (D : ℤ) := by exact_mod_cast hD₁
    have hd₂ : (q₂ : ℤ) ≤ (D : ℤ) := by exact_mod_cast hD₂
    exact mul_le_mul hd₁ hd₂ hq₂Z.le (by positivity)
  -- chain: 14·vi ≤ 14·|vi·Δ| ≤ 14(|X₁|q₂ + |X₂|q₁) < 2·q₁q₂ ≤ 2·D² ≤ 14·vi
  linarith [habs, hm1, hm2, hvineq, hDD, hsuper]

/-- **Sharpness at `m = 14`**: the lock breaks one step past 13 — the
concrete witness `v = 1 → 14, q = 15, p = 1` has witnesses `0` vs `1`
(kernel-checked). -/
theorem rung_lock_sharp_at_14 :
    ∃ (vi wi wj : ℤ) (q p : ℕ),
      14 * |vi * (p : ℤ) - wi * q| < q ∧
      14 * |(14 * vi) * (p : ℤ) - wj * q| < q ∧
      wj ≠ 14 * wi := by
  refine ⟨1, 0, 1, 15, 1, ?_, ?_, ?_⟩ <;> decide

/-! ## Axiom audit -/
#print axioms rung_lock
#print axioms rung_lock_sharp_at_14
#print axioms rung_lock_chain3
#print axioms ray_fails_above
#print axioms carrier_instant_eq_of_shared_witness

end LRC14Concrete
end LonelyRunner
