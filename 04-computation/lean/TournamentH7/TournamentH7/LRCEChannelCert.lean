/-
  TournamentH7.LRCEChannelCert — the e-channel certificate format
  (death-star-2026-07-19-S59i; HYP-7970, the THM-1271 lead (xx)).

  A finite, DECIDABLE certificate whose validity pins the exact loneliness
  supremum of an integer family.  Format:

    `EChannelCert v D Q`: for every merge-grid modulus `S = |vᵢ| + |vⱼ|` and
    every multiplier `kk ∈ [0, S)`, some element `l` has its residue
    `(v l · kk) mod S` within `D·S/Q` of `0` or `S` — the KILL.

  Soundness (`margin_le_of_cert`): a valid certificate forces
  `margin v t ≤ D/Q` for EVERY real `t`.  The proof consumes the repo's
  formalized Pinch chain: a putative `t` with `margin > D/Q` yields a global
  maximizer (`exists_global_max_margin`), which lives on the merge grid
  (`maximizer_on_grid` via `grid_margin_domination`), where the margin
  converts to residue bands (`mod_of_margin_at_rat`) — and the certificate's
  killing element sits in the band AND below the kill bound: contradiction.

  Together with `LRCRungFloor.rung_floor_witness` (the L1 floor) this yields
  the EXACT-VALUE package `margin_sSup_eq_of_cert`, and the concrete
  kernel-checked value

    `M({1,2,3,4,5,7,18}) = sSup (margin v '' Icc 0 1) = 3/23`

  — the N=7 canonical mediant, the first first-gap member exact by kernel
  (`member_3_23_exact`), its certificate discharged by `decide` through a
  ℕ/List reflection checker.  (The 4/127 member's checker is ~250× heavier
  and exceeds this box's memory under kernel `decide`; the per-modulus
  splitting of the checker is the named scaling path.)
-/
import TournamentH7.LRCGridAttainment
import TournamentH7.LRCRungFloor

namespace LonelyRunner
namespace EChannelCert

open GridAttainment TournamentH7.LRCWitness

variable {k : ℕ} [Nonempty (Fin k)]

/-- **The e-channel certificate.**  For every merge-grid modulus
`S = |vᵢ| + |vⱼ|` and every multiplier `kk ∈ [0, S)`, some element's residue
is within `D·S/Q` of `0` or `S` (in the integer form
`min(r, S−r)·Q ≤ D·S`).  Decidable for concrete data. -/
def Cert (v : Fin k → ℤ) (D Q : ℤ) : Prop :=
  ∀ i j : Fin k, ∀ kk ∈ Finset.Ico (0:ℤ) (v i + v j),
    ∃ l : Fin k,
      (v l * kk) % (v i + v j) * Q ≤ D * (v i + v j) ∨
      ((v i + v j) - (v l * kk) % (v i + v j)) * Q ≤ D * (v i + v j)

noncomputable instance (v : Fin k → ℤ) (D Q : ℤ) : Decidable (Cert v D Q) := by
  unfold Cert; infer_instance

-- kernel-reduction probe: does `decide` work despite the noncomputable marking?
example : Cert (fun i : Fin 2 => if i = 0 then 3 else 5) 1 2 := by decide

/-- Residues only depend on the multiplier mod `S`. -/
private lemma res_reduce (x kk S : ℤ) :
    (x * kk) % S = (x * (kk % S)) % S := by
  rw [Int.mul_emod, Int.mul_emod x (kk % S), Int.emod_emod_of_dvd _ dvd_rfl]

/-- **Soundness.**  A valid e-channel certificate bounds the margin at EVERY
real time by `D/Q`. -/
theorem margin_le_of_cert (v : Fin k → ℤ) (hv : ∀ i, 0 < v i)
    (D Q : ℤ) (hD : 0 < D) (hQ : 0 < Q)
    (hcert : Cert v D Q) :
    ∀ t : ℝ, margin v t ≤ (D : ℝ) / (Q : ℝ) := by
  intro t
  by_contra hcon
  push Not at hcon
  set β : ℝ := margin v t with hβdef
  have hDQ : (0:ℝ) < (D : ℝ) / (Q : ℝ) := by
    apply div_pos <;> exact_mod_cast ‹_›
  have hβpos : 0 < β := lt_trans hDQ hcon
  -- margin ≥ β at t, pointwise form
  have hm : ∀ i, ∀ m : ℤ, β ≤ |(v i : ℝ) * t - m| :=
    (le_margin_iff v β t).1 (le_of_eq hβdef)
  -- domination: β at a merge-grid point
  have hvne : ∀ i, v i ≠ 0 := fun i => ne_of_gt (hv i)
  obtain ⟨i, j, kk, hgridm⟩ := grid_margin_domination v hvne β hβpos t hm
  -- positivity turns |vᵢ|+|vⱼ| into vᵢ+vⱼ
  rw [abs_of_pos (hv i), abs_of_pos (hv j)] at hgridm
  set S : ℤ := v i + v j with hSdef
  have hS0 : (0:ℤ) < S := by
    have h1 := hv i; have h2 := hv j; omega
  -- residue band at the grid point
  have hband := mod_of_margin_at_rat v β kk S hS0 hgridm
  -- the certificate's kill at the reduced multiplier
  have hkkmem : kk % S ∈ Finset.Ico (0:ℤ) S := by
    simp only [Finset.mem_Ico]
    exact ⟨Int.emod_nonneg _ (ne_of_gt hS0), Int.emod_lt_of_pos _ hS0⟩
  obtain ⟨l, hl⟩ := hcert i j (kk % S) hkkmem
  rw [← res_reduce (v l) kk S] at hl
  -- the killing element's residue r obeys both the band and the kill
  set r : ℤ := (v l * kk) % S with hrdef
  obtain ⟨hb1, hb2⟩ := hband l
  -- kill (either side) + band ⟹ ⌈β·S⌉·Q ≤ D·S
  have hkill : (⌈β * (S:ℝ)⌉ : ℤ) * Q ≤ D * S := by
    rcases hl with h | h
    · calc (⌈β * (S:ℝ)⌉ : ℤ) * Q ≤ r * Q :=
            mul_le_mul_of_nonneg_right (by omega) (le_of_lt hQ)
        _ ≤ D * S := h
    · calc (⌈β * (S:ℝ)⌉ : ℤ) * Q ≤ (S - r) * Q :=
            mul_le_mul_of_nonneg_right (by omega) (le_of_lt hQ)
        _ ≤ D * S := h
  -- but β > D/Q ⟹ β·S > D·S/Q ⟹ ⌈β·S⌉·Q > D·S: contradiction
  have hceil : (β * (S:ℝ)) ≤ ((⌈β * (S:ℝ)⌉ : ℤ) : ℝ) := Int.le_ceil _
  have hSR : (0:ℝ) < (S : ℝ) := by exact_mod_cast hS0
  have hQR : (0:ℝ) < (Q : ℝ) := by exact_mod_cast hQ
  have hchain : (D : ℝ) * S < ((⌈β * (S:ℝ)⌉ : ℤ) : ℝ) * Q := by
    have h1 : (D : ℝ) / Q * S < β * S := by
      exact mul_lt_mul_of_pos_right hcon hSR
    have h2 : (D : ℝ) * S = ((D : ℝ) / Q * S) * Q := by field_simp
    rw [h2]
    have h3 : ((D : ℝ) / Q * S) * Q < (β * S) * Q :=
      mul_lt_mul_of_pos_right h1 hQR
    calc ((D : ℝ) / Q * S) * Q < (β * S) * Q := h3
      _ ≤ ((⌈β * (S:ℝ)⌉ : ℤ) : ℝ) * Q :=
          mul_le_mul_of_nonneg_right hceil (le_of_lt hQR)
  -- integer side: ⌈β·S⌉·Q ≤ min·Q ≤ D·S
  have hint : ((⌈β * (S:ℝ)⌉ : ℤ) : ℝ) * Q ≤ (D : ℝ) * S := by
    exact_mod_cast hkill
  linarith

/-! ## The exact-value package -/

/-- **Exact loneliness value from a floor witness plus an e-channel
certificate.**  If some element list realizes the band floor at `a/Q` (via
`RungFloor.rung_floor_witness`-style data: every element at distance `≥ D/Q`)
and the certificate holds, then the margin supremum over the period is
EXACTLY `D/Q`. -/
theorem margin_sSup_eq_of_cert (v : Fin k → ℤ) (hv : ∀ i, 0 < v i)
    (D Q : ℤ) (hD : 0 < D) (hQ : 0 < Q)
    (t₀ : ℝ) (ht₀ : t₀ ∈ Set.Icc (0:ℝ) 1)
    (hfloor : ∀ i, ∀ m : ℤ, (D : ℝ) / (Q : ℝ) ≤ |(v i : ℝ) * t₀ - m|)
    (hcert : Cert v D Q) :
    sSup (margin v '' Set.Icc (0:ℝ) 1) = (D : ℝ) / (Q : ℝ) := by
  have hub : ∀ x ∈ margin v '' Set.Icc (0:ℝ) 1, x ≤ (D : ℝ) / (Q : ℝ) := by
    rintro x ⟨t, -, rfl⟩
    exact margin_le_of_cert v hv D Q hD hQ hcert t
  have hne : (margin v '' Set.Icc (0:ℝ) 1).Nonempty :=
    ⟨margin v 0, ⟨0, by constructor <;> norm_num, rfl⟩⟩
  have hbdd : BddAbove (margin v '' Set.Icc (0:ℝ) 1) :=
    ⟨(D : ℝ) / (Q : ℝ), fun x hx => hub x hx⟩
  apply le_antisymm
  · exact Real.sSup_le (fun x hx => hub x hx) (by positivity)
  · have hmem : margin v t₀ ∈ margin v '' Set.Icc (0:ℝ) 1 := ⟨t₀, ht₀, rfl⟩
    have hfl : (D : ℝ) / (Q : ℝ) ≤ margin v t₀ :=
      (le_margin_iff v _ t₀).2 hfloor
    exact le_trans hfl (le_csSup hbdd hmem)

/-! ## The concrete kernel-exact tower member: M({1..29, 31, 120}) = 4/127 -/

/-- The computable ℕ/List form of the certificate check (kernel-fast:
Nat arithmetic is GMP-backed in the kernel). -/
def certCheck (l : List ℕ) (D Q : ℕ) : Bool :=
  l.all fun vi => l.all fun vj =>
    (List.range (vi + vj)).all fun kk =>
      l.any fun w =>
        decide ((w * kk) % (vi + vj) * Q ≤ D * (vi + vj)) ||
        decide (((vi + vj) - (w * kk) % (vi + vj)) * Q ≤ D * (vi + vj))

/-- **Reflection**: a passing ℕ-checker yields the `Cert` proposition for any
positive family whose values all appear in the list and vice versa. -/
theorem cert_of_check (l : List ℕ) (D Q : ℕ)
    (hchk : certCheck l D Q = true)
    (v : Fin k → ℤ) (hv : ∀ i, 0 < v i)
    (hmem : ∀ i, (v i).toNat ∈ l)
    (hall : ∀ x ∈ l, ∃ i, v i = (x : ℤ)) :
    Cert v (D : ℤ) (Q : ℤ) := by
  intro i j kk hkk
  simp only [Finset.mem_Ico] at hkk
  obtain ⟨hkk0, hkkS⟩ := hkk
  have hvi := hv i
  have hvj := hv j
  set Sn : ℕ := (v i).toNat + (v j).toNat with hSn
  have hScast : (Sn : ℤ) = v i + v j := by
    simp only [hSn]
    push_cast [Int.toNat_of_nonneg (le_of_lt hvi), Int.toNat_of_nonneg (le_of_lt hvj)]
    ring
  have hkkn : (kk.toNat : ℤ) = kk := Int.toNat_of_nonneg hkk0
  have hkklt : kk.toNat < Sn := by
    have hlt : kk < (Sn : ℤ) := by rw [hScast]; exact hkkS
    omega
  simp only [certCheck, List.all_eq_true, List.any_eq_true] at hchk
  obtain ⟨w, hwmem, hwkill⟩ :=
    hchk _ (hmem i) _ (hmem j) kk.toNat (List.mem_range.mpr hkklt)
  obtain ⟨lidx, hlidx⟩ := hall w hwmem
  refine ⟨lidx, ?_⟩
  rw [Bool.or_eq_true, decide_eq_true_eq, decide_eq_true_eq] at hwkill
  have hSn0 : 0 < Sn := by
    have := hv i
    omega
  rcases hwkill with h | h
  · left
    have hz : ((w * kk.toNat % Sn * Q : ℕ) : ℤ) ≤ ((D * Sn : ℕ) : ℤ) := by
      exact_mod_cast h
    push_cast at hz
    rw [hkkn, hScast, ← hlidx] at hz
    exact hz
  · right
    have hle : w * kk.toNat % Sn ≤ Sn := le_of_lt (Nat.mod_lt _ hSn0)
    have hz : (((Sn - w * kk.toNat % Sn) * Q : ℕ) : ℤ) ≤ ((D * Sn : ℕ) : ℤ) := by
      exact_mod_cast h
    push_cast [hle] at hz
    rw [hkkn, hScast, ← hlidx] at hz
    exact hz

/-! ### The flagship kernel-exact member: M({1,2,3,4,5,7,18}) = 3/23
(the N=7 canonical mediant -- the original first-gap member).  Its checker is
small enough for a direct kernel `decide`; the floor comes element-by-element
through `rung_floor_single` at the closed-form witness t0 = 19/23. -/

/-- The N=7 mediant member as a `Fin 7` family: `{1,2,3,4,5,7,18}`. -/
def v7 : Fin 7 → ℤ := ![1, 2, 3, 4, 5, 7, 18]

/-- The value list of the 3/23 member. -/
def l7 : List ℕ := [1, 2, 3, 4, 5, 7, 18]

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8000 in
/-- The ℕ-checker passes on the 3/23 member (kernel computation). -/
theorem check_3_23 : certCheck l7 3 23 = true := by decide

/-- The e-channel certificate of the 3/23 member. -/
theorem cert_3_23 : Cert v7 3 23 := by
  have h := cert_of_check l7 3 23 check_3_23 v7
    (by decide) (by decide) (by decide)
  exact_mod_cast h

/-- **The first kernel-exact first-gap member.**
`sSup (margin v7 '' [0,1]) = 3/23`: the exact loneliness supremum of
`{1,2,3,4,5,7,18}` -- the N=7 canonical mediant, machine-checked end to end.
Floor: the residue band of the closed-form witness `t0 = 19/23`
(`a = 3·5⁻¹ mod 23`), element-by-element through `rung_floor_single`.
Ceiling: the e-channel certificate through the formal Pinch chain. -/
theorem member_3_23_exact :
    sSup (margin v7 '' Set.Icc (0:ℝ) 1) = (3 : ℝ) / (23 : ℝ) := by
  have h3 : ((3:ℤ) : ℝ) = (3:ℝ) := by norm_num
  have h23 : ((23:ℤ) : ℝ) = (23:ℝ) := by norm_num
  rw [← h3, ← h23]
  apply margin_sSup_eq_of_cert v7 (by decide) 3 23 (by norm_num) (by norm_num)
    ((19 : ℝ) / 23) (by constructor <;> norm_num)
  · -- the floor: each element's residue at a = 19 lies in the band [3, 20]
    intro i m
    have hband : 3 ≤ (v7 i * 19) % 23 ∧ (v7 i * 19) % 23 ≤ 23 - 3 := by
      fin_cases i <;> decide
    have h := RungFloor.rung_floor_single 23 3 (v7 i) 19 (by norm_num) hband m
    have hc : ((23 : ℤ) : ℝ) = (23 : ℝ) := by norm_num
    rw [hc] at h
    convert h using 2 <;> norm_num
  · exact cert_3_23

end EChannelCert
end LonelyRunner

#print axioms LonelyRunner.EChannelCert.margin_le_of_cert
#print axioms LonelyRunner.EChannelCert.margin_sSup_eq_of_cert
#print axioms LonelyRunner.EChannelCert.cert_of_check
#print axioms LonelyRunner.EChannelCert.check_3_23
#print axioms LonelyRunner.EChannelCert.cert_3_23
#print axioms LonelyRunner.EChannelCert.member_3_23_exact
