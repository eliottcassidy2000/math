/-
  TournamentH7.LRCCertCompleteness — CERTIFICATE COMPLETENESS (the atom's converse,
  with an EXPLICIT MODULUS BOUND) and THE GAP FILTERS
  (kind-pasteur-2026-07-05-S3, HYP-4105).

  COMPLETENESS (the S2 "next brick"): the rational-point margin certificate
  (`rational_point_margin`, HYP-4102) is COMPLETE — every real margin witness
  converts to an integer certificate with bounded modulus:

    `cert_of_margin`: margin `β'` at ANY real `t*`, speeds bounded by `B`, and any
    modulus `s` with `B ≤ 2s(β'−β)` ⟹ the atom conditions hold at
    `(k, s, μ) = (round(s·t*), s, ⌈β·s⌉)`.

  Proof is pure Lipschitz transfer + round approximation — no equioscillation /
  breakpoint structure theory (THM-592) needed.  Instantiated at the verified
  rigidity slack (second value `14/169` vs margin `2/25`, gap `12/4225`):
  any `s ≥ (4225/24)·B ≈ 176·B` suffices (`loose_branch_cert_exists` packages the
  existential with the bound `s ≤ B/(2(β'−β)) + 1`).  Together with HYP-4102 the
  loose branch of `TightLooseDichotomy` is a BOUNDED integer search: the margin
  language and the integer-certificate language are one and the same.

  THE GAP FILTERS: what the FAILURE of the loose branch forces.  If no margin-β
  point exists (`∀ t, ∃ i m, |vᵢt − m| < β`), then at every rational `a/q`:
    * `not_loose_dvd`  (`β·q ≤ 1`; q ≤ 12 at β = 2/25):  `q ∣ vᵢ·a` for some `i`;
    * `not_loose_near_unit` (`β·q ≤ 2`; q = 13, 14 at 2/25):  `vᵢ·a ≡ 0, ±1 (mod q)`
      for some `i` — residue pinning INTO THE GAP.  opus's HYP-4098 pinning is the
      exact-tight (`M = 1/13`) version; these run at any `M < 2/25`, i.e. on the
      whole hypothetical spectral gap `(1/13, 2/25)`.
    * `pair_pinning_13`: with no 13-multiples, every unit ±pair class mod 13 is hit.
  These are the dichotomy prover's working hypotheses: a proof of
  `TightLooseDichotomy` argues from exactly this congruence package plus the
  spread gate (HYP-4102).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace CertCompleteness

open HarmonicGate

/-! ## Lipschitz margin transfer -/

/-- Margin transfers between nearby points: margin `β'` at `t`, speeds `≤ B`, and
`B·|t' − t| ≤ β' − β` give margin `β` at `t'`. -/
theorem margin_transfer {ι : Type*} (v : ι → ℤ) (t t' β β' B : ℝ)
    (hB : ∀ i, |(v i : ℝ)| ≤ B)
    (hmargin : ∀ i, ∀ m : ℤ, β' ≤ |(v i : ℝ) * t - m|)
    (hclose : B * |t' - t| ≤ β' - β) :
    ∀ i, ∀ m : ℤ, β ≤ |(v i : ℝ) * t' - m| := by
  intro i m
  have h1 := hmargin i m
  have htri : |(v i : ℝ) * t - m| ≤ |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t' - t| := by
    calc |(v i : ℝ) * t - m|
        = |((v i : ℝ) * t' - m) + (v i : ℝ) * (t - t')| := by
          rw [show (v i : ℝ) * t - m = ((v i : ℝ) * t' - m) + (v i : ℝ) * (t - t') by ring]
      _ ≤ |(v i : ℝ) * t' - m| + |(v i : ℝ) * (t - t')| := abs_add_le _ _
      _ = |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t - t'| := by rw [abs_mul]
      _ = |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t' - t| := by rw [abs_sub_comm t t']
  have hBt : |(v i : ℝ)| * |t' - t| ≤ B * |t' - t| :=
    mul_le_mul_of_nonneg_right (hB i) (abs_nonneg _)
  linarith

/-! ## The completeness theorem -/

/-- **CERTIFICATE COMPLETENESS.**  A real margin witness converts to an integer
certificate at ANY sufficiently large modulus: margin `β'` at `t*`, speeds `≤ B`,
`B ≤ 2s(β'−β)` ⟹ the atom's mod conditions hold at `k = round(s·t*)` with
`μ = ⌈β·s⌉`.  Converse of `rational_point_margin`: the margin language and the
integer-certificate language coincide, with explicit modulus bounds. -/
theorem cert_of_margin {ι : Type*} (v : ι → ℤ) (tstar : ℝ) (β β' : ℝ) (B s : ℤ)
    (hs : 0 < s) (hB0 : 0 ≤ B)
    (hB : ∀ i, |v i| ≤ B)
    (hmargin : ∀ i, ∀ m : ℤ, β' ≤ |(v i : ℝ) * tstar - m|)
    (hsbound : (B : ℝ) ≤ 2 * s * (β' - β)) :
    ∀ i, ⌈β * s⌉ ≤ (v i * round ((s : ℝ) * tstar)) % s ∧
         (v i * round ((s : ℝ) * tstar)) % s ≤ s - ⌈β * s⌉ := by
  set k : ℤ := round ((s : ℝ) * tstar) with hk
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  have h2s : (0 : ℝ) < 2 * (s : ℝ) := by positivity
  -- |k/s − t*| ≤ 1/(2s)
  have hround : |(k : ℝ) / s - tstar| ≤ 1 / (2 * s) := by
    have h1 : |(s : ℝ) * tstar - k| ≤ 1 / 2 := abs_sub_round ((s : ℝ) * tstar)
    have h2 : (k : ℝ) / s - tstar = ((k : ℝ) - (s : ℝ) * tstar) / s := by
      field_simp
    rw [h2, abs_div, abs_of_pos hsR, abs_sub_comm, div_le_div_iff₀ hsR h2s]
    calc |(s : ℝ) * tstar - k| * (2 * s) ≤ (1 / 2) * (2 * s) :=
          mul_le_mul_of_nonneg_right h1 (le_of_lt h2s)
      _ = 1 * (s : ℝ) := by ring
  -- the Lipschitz budget: B·|k/s − t*| ≤ β' − β
  have hclose : ((B : ℤ) : ℝ) * |(k : ℝ) / s - tstar| ≤ β' - β := by
    have hB0R : (0 : ℝ) ≤ ((B : ℤ) : ℝ) := by exact_mod_cast hB0
    have hstep : ((B : ℤ) : ℝ) * |(k : ℝ) / s - tstar| ≤ (B : ℝ) * (1 / (2 * s)) :=
      mul_le_mul_of_nonneg_left hround hB0R
    have hdiv : (B : ℝ) * (1 / (2 * s)) ≤ β' - β := by
      rw [mul_one_div, div_le_iff₀ h2s]
      calc (B : ℝ) ≤ 2 * s * (β' - β) := hsbound
        _ = (β' - β) * (2 * s) := by ring
    linarith
  -- margin β at k/s (real form)
  have hm : ∀ i, ∀ m : ℤ, β ≤ |(v i : ℝ) * ((k : ℝ) / s) - m| :=
    margin_transfer v tstar ((k : ℝ) / s) β β' ((B : ℤ) : ℝ)
      (fun i => by rw [← Int.cast_abs]; exact_mod_cast hB i) hmargin hclose
  -- margin ⟹ mod conditions
  intro i
  have hint : ∀ m : ℤ, β * s ≤ ((|v i * k - m * s| : ℤ) : ℝ) := by
    intro m
    have h := hm i m
    have heq : (v i : ℝ) * ((k : ℝ) / s) - m = ((v i * k - m * s : ℤ) : ℝ) / s := by
      push_cast
      field_simp
    rw [heq, abs_div, abs_of_pos hsR, le_div_iff₀ hsR, ← Int.cast_abs] at h
    exact h
  set r : ℤ := (v i * k) % s with hrdef
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (ne_of_gt hs)
  have hrs : r < s := Int.emod_lt_of_pos _ hs
  have hX0 : v i * k - ((v i * k) / s) * s = r := by
    rw [hrdef, Int.emod_def]; ring
  have hX1 : v i * k - ((v i * k) / s + 1) * s = r - s := by
    rw [hrdef, Int.emod_def]; ring
  constructor
  · -- ⌈βs⌉ ≤ r
    have h0 := hint ((v i * k) / s)
    rw [hX0, abs_of_nonneg hr0] at h0
    exact Int.ceil_le.mpr h0
  · -- r ≤ s − ⌈βs⌉  ⟸  ⌈βs⌉ ≤ s − r
    have h1 := hint ((v i * k) / s + 1)
    rw [hX1, abs_of_nonpos (by omega), neg_sub] at h1
    have : ⌈β * s⌉ ≤ s - r := Int.ceil_le.mpr (by exact_mod_cast h1)
    omega

/-- **The packaged loose-branch certificate with the explicit modulus bound**:
a real margin-`β'` witness for the base yields an integer certificate
`(k, s, ⌈β·s⌉)` with `s ≤ B/(2(β'−β)) + 1`.  The loose branch of the dichotomy
is a bounded integer search. -/
theorem loose_branch_cert_exists (v : Fin 13 → ℤ) (istar : Fin 13) (β β' : ℝ) (B : ℤ)
    (hB0 : 0 < B) (hββ' : β < β')
    (hBb : ∀ i, i ≠ istar → |v i| ≤ B)
    (hloose : ∃ t : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β' ≤ |(v i : ℝ) * t - m|) :
    ∃ k s : ℤ, 0 < s ∧ (s : ℝ) ≤ (B : ℝ) / (2 * (β' - β)) + 1 ∧
      ∀ i, i ≠ istar →
        ⌈β * s⌉ ≤ (v i * k) % s ∧ (v i * k) % s ≤ s - ⌈β * s⌉ := by
  obtain ⟨tstar, hm⟩ := hloose
  have hδ : (0 : ℝ) < β' - β := by linarith
  have hδ2 : (0 : ℝ) < 2 * (β' - β) := by linarith
  set s : ℤ := max ⌈(B : ℝ) / (2 * (β' - β))⌉ 1 with hsdef
  have hs : 0 < s := lt_of_lt_of_le one_pos (le_max_right _ _)
  refine ⟨round ((s : ℝ) * tstar), s, hs, ?_, ?_⟩
  · -- s ≤ B/(2(β'−β)) + 1
    have hBpos : (0 : ℝ) < (B : ℝ) / (2 * (β' - β)) :=
      div_pos (by exact_mod_cast hB0) hδ2
    have hceil1 : 1 ≤ ⌈(B : ℝ) / (2 * (β' - β))⌉ := Int.ceil_pos.mpr hBpos
    have hsval : s = ⌈(B : ℝ) / (2 * (β' - β))⌉ := by
      rw [hsdef]; omega
    rw [hsval]
    exact le_of_lt (Int.ceil_lt_add_one _)
  · -- the certificate conditions via cert_of_margin over the base subtype
    have hsbound : (B : ℝ) ≤ 2 * s * (β' - β) := by
      have h1 : (B : ℝ) / (2 * (β' - β)) ≤ ((⌈(B : ℝ) / (2 * (β' - β))⌉ : ℤ) : ℝ) :=
        Int.le_ceil _
      have h2 : ((⌈(B : ℝ) / (2 * (β' - β))⌉ : ℤ) : ℝ) ≤ (s : ℝ) := by
        exact_mod_cast le_max_left _ _
      have h3 : (B : ℝ) / (2 * (β' - β)) ≤ (s : ℝ) := le_trans h1 h2
      rw [div_le_iff₀ hδ2] at h3
      calc (B : ℝ) ≤ s * (2 * (β' - β)) := h3
        _ = 2 * s * (β' - β) := by ring
    have := cert_of_margin (fun j : {j : Fin 13 // j ≠ istar} => v j) tstar β β' B s
      hs (le_of_lt hB0) (fun j => hBb j j.2) (fun j m => hm j j.2 m) hsbound
    intro i hi
    exact this ⟨i, hi⟩

/-! ## The gap filters -/

/-- **The evaluation lemma**: if NO margin-`β` point exists, then at every rational
`a/q` some runner is `β`-close in cleared-denominator form: `|vᵢ·a − q·m| < β·q`. -/
theorem not_loose_eval {ι : Type*} (v : ι → ℤ) (β : ℝ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < β)
    (q a : ℤ) (hq : 0 < q) :
    ∃ i, ∃ m : ℤ, ((|v i * a - q * m| : ℤ) : ℝ) < β * q := by
  obtain ⟨i, m, h⟩ := hnl ((a : ℝ) / q)
  refine ⟨i, m, ?_⟩
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have heq : (v i : ℝ) * ((a : ℝ) / q) - m = ((v i * a - q * m : ℤ) : ℝ) / q := by
    push_cast
    field_simp
  rw [heq, abs_div, abs_of_pos hqR, div_lt_iff₀ hqR, ← Int.cast_abs] at h
  exact h

/-- **Covering filter** (`β·q ≤ 1`): failure of the loose branch forces a multiple of
`q` in every unit direction `a` — at `β = 2/25` this is every `q ≤ 12`. -/
theorem not_loose_dvd {ι : Type*} (v : ι → ℤ) (β : ℝ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < β)
    (q a : ℤ) (hq : 0 < q) (hβq : β * q ≤ 1) :
    ∃ i, (q : ℤ) ∣ (v i * a) := by
  obtain ⟨i, m, h⟩ := not_loose_eval v β hnl q a hq
  have h1 : ((|v i * a - q * m| : ℤ) : ℝ) < 1 := lt_of_lt_of_le h hβq
  have h2 : |v i * a - q * m| < 1 := by exact_mod_cast h1
  have h4 := abs_lt.mp h2
  exact ⟨i, m, by omega⟩

/-- **Near-unit filter** (`β·q ≤ 2`): failure of the loose branch forces
`vᵢ·a ≡ 0, ±1 (mod q)` — at `β = 2/25` this runs at `q = 13` and `q = 14`,
extending the residue pinning into the whole gap `M < 2/25`. -/
theorem not_loose_near_unit {ι : Type*} (v : ι → ℤ) (β : ℝ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < β)
    (q a : ℤ) (hq : 0 < q) (hβq : β * q ≤ 2) :
    ∃ i, ∃ m : ℤ, |v i * a - q * m| ≤ 1 := by
  obtain ⟨i, m, h⟩ := not_loose_eval v β hnl q a hq
  have h1 : ((|v i * a - q * m| : ℤ) : ℝ) < 2 := lt_of_lt_of_le h hβq
  have h2 : |v i * a - q * m| < 2 := by exact_mod_cast h1
  exact ⟨i, m, by omega⟩

/-- **The mod-13 gap pinning**: no margin-`2/25` point ⟹ for every direction `a`,
some runner has `vᵢ·a ≡ 0, 1, or 12 (mod 13)`.  (opus's HYP-4098 pinning at the
exact-tight threshold `1/13`, extended to the whole gap `M < 2/25`.) -/
theorem not_loose_pinning_13 {ι : Type*} (v : ι → ℤ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) (a : ℤ) :
    ∃ i, (v i * a) % 13 = 0 ∨ (v i * a) % 13 = 1 ∨ (v i * a) % 13 = 12 := by
  obtain ⟨i, m, h⟩ := not_loose_near_unit v (2 / 25) hnl 13 a (by norm_num) (by norm_num)
  have h' := abs_le.mp h
  exact ⟨i, by omega⟩

/-- **The mod-14 gap pinning**: same at the composite modulus 14. -/
theorem not_loose_pinning_14 {ι : Type*} (v : ι → ℤ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) (a : ℤ) :
    ∃ i, (v i * a) % 14 = 0 ∨ (v i * a) % 14 = 1 ∨ (v i * a) % 14 = 13 := by
  obtain ⟨i, m, h⟩ := not_loose_near_unit v (2 / 25) hnl 14 a (by norm_num) (by norm_num)
  have h' := abs_le.mp h
  exact ⟨i, by omega⟩

/-- **The ±pair covering**: no margin-`2/25` point and no 13-multiples ⟹ every unit
±pair class mod 13 is hit by some runner: given `u·a ≡ 1 (mod 13)`, some `vᵢ ≡ ±u`. -/
theorem pair_pinning_13 {ι : Type*} (v : ι → ℤ)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (hnd : ∀ i, ¬ ((13 : ℤ) ∣ v i))
    (u a : ℤ) (hua : (u * a) % 13 = 1) :
    ∃ i, (13 : ℤ) ∣ (v i - u) ∨ (13 : ℤ) ∣ (v i + u) := by
  obtain ⟨i, hcase⟩ := not_loose_pinning_13 v hnl a
  have hp13 : Prime (13 : ℤ) := by
    rw [Int.prime_iff_natAbs_prime]; norm_num
  have hnda : ¬ ((13 : ℤ) ∣ a) := by
    intro hdvd
    have h13 : (u * a) % 13 = 0 := by
      obtain ⟨c, rfl⟩ := hdvd
      rw [show u * (13 * c) = 13 * (u * c) by ring, Int.mul_emod_right]
    omega
  refine ⟨i, ?_⟩
  rcases hcase with h0 | h1 | h12
  · -- vᵢ·a ≡ 0 is impossible: 13 prime, 13 ∤ vᵢ, 13 ∤ a
    exfalso
    have hdvd : (13 : ℤ) ∣ v i * a := Int.dvd_of_emod_eq_zero h0
    rcases hp13.dvd_mul.mp hdvd with hv | ha
    · exact hnd i hv
    · exact hnda ha
  · -- vᵢ·a ≡ 1 ≡ u·a ⟹ 13 ∣ a(vᵢ − u) ⟹ 13 ∣ vᵢ − u
    left
    have hdvd : (13 : ℤ) ∣ (v i - u) * a := by
      have : (v i * a - u * a) % 13 = 0 := by omega
      have h' : (13 : ℤ) ∣ v i * a - u * a := Int.dvd_of_emod_eq_zero this
      have heq : v i * a - u * a = (v i - u) * a := by ring
      rwa [heq] at h'
    rcases hp13.dvd_mul.mp hdvd with h | ha
    · exact h
    · exact absurd ha hnda
  · -- vᵢ·a ≡ 12 ≡ −1 ≡ −u·a ⟹ 13 ∣ a(vᵢ + u) ⟹ 13 ∣ vᵢ + u
    right
    have hdvd : (13 : ℤ) ∣ (v i + u) * a := by
      have : (v i * a + u * a) % 13 = 0 := by omega
      have h' : (13 : ℤ) ∣ v i * a + u * a := Int.dvd_of_emod_eq_zero this
      have heq : v i * a + u * a = (v i + u) * a := by ring
      rwa [heq] at h'
    rcases hp13.dvd_mul.mp hdvd with h | ha
    · exact h
    · exact absurd ha hnda

#print axioms margin_transfer
#print axioms cert_of_margin
#print axioms loose_branch_cert_exists
#print axioms not_loose_dvd
#print axioms not_loose_pinning_13
#print axioms pair_pinning_13

end CertCompleteness
end LonelyRunner
