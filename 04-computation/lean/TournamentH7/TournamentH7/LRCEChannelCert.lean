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

/-- Per-modulus checker (the computable ℕ/List form; kernel-GMP arithmetic):
all multipliers at ONE modulus `S` are killed. -/
def certCheckS (l : List ℕ) (D Q S : ℕ) : Bool :=
  (List.range S).all fun kk =>
    l.any fun w =>
      decide ((w * kk) % S * Q ≤ D * S) ||
      decide ((S - (w * kk) % S) * Q ≤ D * S)

def certCheck (l : List ℕ) (D Q : ℕ) : Bool :=
  l.all fun vi => l.all fun vj => certCheckS l D Q (vi + vj)

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
  simp only [certCheck, certCheckS, List.all_eq_true, List.any_eq_true] at hchk
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


/-! ### The kernel-exact 4/127 TOWER member, via the per-modulus split.
Each modulus is its own small `decide` (the monolithic checker exceeded a
7 GB box); an ordered coverage lemma reassembles `certCheck`. -/

/-- The first D=4 tower member as a `Fin 31` family: `{1,…,29, 31, 120}`. -/
def v31 : Fin 31 → ℤ := fun i =>
  if (i : ℕ) < 29 then (i : ℕ) + 1 else if (i : ℕ) = 29 then 31 else 120

/-- The value list of the 4/127 member. -/
def l31 : List ℕ :=
  [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,31,120]

/-- The 91 distinct pair-sum moduli of `l31`. -/
def moduli31 : List ℕ := [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 151, 240]

section PerModulusChecks
set_option maxRecDepth 4000
theorem chk127_2 : certCheckS l31 4 127 2 = true := by decide
theorem chk127_3 : certCheckS l31 4 127 3 = true := by decide
theorem chk127_4 : certCheckS l31 4 127 4 = true := by decide
theorem chk127_5 : certCheckS l31 4 127 5 = true := by decide
theorem chk127_6 : certCheckS l31 4 127 6 = true := by decide
theorem chk127_7 : certCheckS l31 4 127 7 = true := by decide
theorem chk127_8 : certCheckS l31 4 127 8 = true := by decide
theorem chk127_9 : certCheckS l31 4 127 9 = true := by decide
theorem chk127_10 : certCheckS l31 4 127 10 = true := by decide
theorem chk127_11 : certCheckS l31 4 127 11 = true := by decide
theorem chk127_12 : certCheckS l31 4 127 12 = true := by decide
theorem chk127_13 : certCheckS l31 4 127 13 = true := by decide
theorem chk127_14 : certCheckS l31 4 127 14 = true := by decide
theorem chk127_15 : certCheckS l31 4 127 15 = true := by decide
theorem chk127_16 : certCheckS l31 4 127 16 = true := by decide
theorem chk127_17 : certCheckS l31 4 127 17 = true := by decide
theorem chk127_18 : certCheckS l31 4 127 18 = true := by decide
theorem chk127_19 : certCheckS l31 4 127 19 = true := by decide
theorem chk127_20 : certCheckS l31 4 127 20 = true := by decide
theorem chk127_21 : certCheckS l31 4 127 21 = true := by decide
theorem chk127_22 : certCheckS l31 4 127 22 = true := by decide
theorem chk127_23 : certCheckS l31 4 127 23 = true := by decide
theorem chk127_24 : certCheckS l31 4 127 24 = true := by decide
theorem chk127_25 : certCheckS l31 4 127 25 = true := by decide
theorem chk127_26 : certCheckS l31 4 127 26 = true := by decide
theorem chk127_27 : certCheckS l31 4 127 27 = true := by decide
theorem chk127_28 : certCheckS l31 4 127 28 = true := by decide
theorem chk127_29 : certCheckS l31 4 127 29 = true := by decide
theorem chk127_30 : certCheckS l31 4 127 30 = true := by decide
theorem chk127_31 : certCheckS l31 4 127 31 = true := by decide
theorem chk127_32 : certCheckS l31 4 127 32 = true := by decide
theorem chk127_33 : certCheckS l31 4 127 33 = true := by decide
theorem chk127_34 : certCheckS l31 4 127 34 = true := by decide
theorem chk127_35 : certCheckS l31 4 127 35 = true := by decide
theorem chk127_36 : certCheckS l31 4 127 36 = true := by decide
theorem chk127_37 : certCheckS l31 4 127 37 = true := by decide
theorem chk127_38 : certCheckS l31 4 127 38 = true := by decide
theorem chk127_39 : certCheckS l31 4 127 39 = true := by decide
theorem chk127_40 : certCheckS l31 4 127 40 = true := by decide
theorem chk127_41 : certCheckS l31 4 127 41 = true := by decide
theorem chk127_42 : certCheckS l31 4 127 42 = true := by decide
theorem chk127_43 : certCheckS l31 4 127 43 = true := by decide
theorem chk127_44 : certCheckS l31 4 127 44 = true := by decide
theorem chk127_45 : certCheckS l31 4 127 45 = true := by decide
theorem chk127_46 : certCheckS l31 4 127 46 = true := by decide
theorem chk127_47 : certCheckS l31 4 127 47 = true := by decide
theorem chk127_48 : certCheckS l31 4 127 48 = true := by decide
theorem chk127_49 : certCheckS l31 4 127 49 = true := by decide
theorem chk127_50 : certCheckS l31 4 127 50 = true := by decide
theorem chk127_51 : certCheckS l31 4 127 51 = true := by decide
theorem chk127_52 : certCheckS l31 4 127 52 = true := by decide
theorem chk127_53 : certCheckS l31 4 127 53 = true := by decide
theorem chk127_54 : certCheckS l31 4 127 54 = true := by decide
theorem chk127_55 : certCheckS l31 4 127 55 = true := by decide
theorem chk127_56 : certCheckS l31 4 127 56 = true := by decide
theorem chk127_57 : certCheckS l31 4 127 57 = true := by decide
theorem chk127_58 : certCheckS l31 4 127 58 = true := by decide
theorem chk127_59 : certCheckS l31 4 127 59 = true := by decide
theorem chk127_60 : certCheckS l31 4 127 60 = true := by decide
theorem chk127_62 : certCheckS l31 4 127 62 = true := by decide
theorem chk127_121 : certCheckS l31 4 127 121 = true := by decide
theorem chk127_122 : certCheckS l31 4 127 122 = true := by decide
theorem chk127_123 : certCheckS l31 4 127 123 = true := by decide
theorem chk127_124 : certCheckS l31 4 127 124 = true := by decide
theorem chk127_125 : certCheckS l31 4 127 125 = true := by decide
theorem chk127_126 : certCheckS l31 4 127 126 = true := by decide
theorem chk127_127 : certCheckS l31 4 127 127 = true := by decide
theorem chk127_128 : certCheckS l31 4 127 128 = true := by decide
theorem chk127_129 : certCheckS l31 4 127 129 = true := by decide
theorem chk127_130 : certCheckS l31 4 127 130 = true := by decide
theorem chk127_131 : certCheckS l31 4 127 131 = true := by decide
theorem chk127_132 : certCheckS l31 4 127 132 = true := by decide
theorem chk127_133 : certCheckS l31 4 127 133 = true := by decide
theorem chk127_134 : certCheckS l31 4 127 134 = true := by decide
theorem chk127_135 : certCheckS l31 4 127 135 = true := by decide
theorem chk127_136 : certCheckS l31 4 127 136 = true := by decide
theorem chk127_137 : certCheckS l31 4 127 137 = true := by decide
theorem chk127_138 : certCheckS l31 4 127 138 = true := by decide
theorem chk127_139 : certCheckS l31 4 127 139 = true := by decide
theorem chk127_140 : certCheckS l31 4 127 140 = true := by decide
theorem chk127_141 : certCheckS l31 4 127 141 = true := by decide
theorem chk127_142 : certCheckS l31 4 127 142 = true := by decide
theorem chk127_143 : certCheckS l31 4 127 143 = true := by decide
theorem chk127_144 : certCheckS l31 4 127 144 = true := by decide
theorem chk127_145 : certCheckS l31 4 127 145 = true := by decide
theorem chk127_146 : certCheckS l31 4 127 146 = true := by decide
theorem chk127_147 : certCheckS l31 4 127 147 = true := by decide
theorem chk127_148 : certCheckS l31 4 127 148 = true := by decide
theorem chk127_149 : certCheckS l31 4 127 149 = true := by decide
theorem chk127_151 : certCheckS l31 4 127 151 = true := by decide
theorem chk127_240 : certCheckS l31 4 127 240 = true := by decide
end PerModulusChecks

/-- All 91 listed moduli pass (assembled from the per-modulus facts —
no case dispatch, just the Bool conjunction). -/
theorem moduli_ok : moduli31.all (certCheckS l31 4 127) = true := by
  simp only [moduli31, List.all_cons, List.all_nil, Bool.and_eq_true]
  exact ⟨chk127_2, chk127_3, chk127_4, chk127_5, chk127_6, chk127_7, chk127_8, chk127_9, chk127_10, chk127_11, chk127_12, chk127_13, chk127_14, chk127_15, chk127_16, chk127_17, chk127_18, chk127_19, chk127_20, chk127_21, chk127_22, chk127_23, chk127_24, chk127_25, chk127_26, chk127_27, chk127_28, chk127_29, chk127_30, chk127_31, chk127_32, chk127_33, chk127_34, chk127_35, chk127_36, chk127_37, chk127_38, chk127_39, chk127_40, chk127_41, chk127_42, chk127_43, chk127_44, chk127_45, chk127_46, chk127_47, chk127_48, chk127_49, chk127_50, chk127_51, chk127_52, chk127_53, chk127_54, chk127_55, chk127_56, chk127_57, chk127_58, chk127_59, chk127_60, chk127_62, chk127_121, chk127_122, chk127_123, chk127_124, chk127_125, chk127_126, chk127_127, chk127_128, chk127_129, chk127_130, chk127_131, chk127_132, chk127_133, chk127_134, chk127_135, chk127_136, chk127_137, chk127_138, chk127_139, chk127_140, chk127_141, chk127_142, chk127_143, chk127_144, chk127_145, chk127_146, chk127_147, chk127_148, chk127_149, chk127_151, chk127_240, trivial⟩

theorem chk127_mem : ∀ S ∈ moduli31, certCheckS l31 4 127 S = true := by
  have h := moduli_ok
  rw [List.all_eq_true] at h
  exact h

set_option maxHeartbeats 2000000 in
/-- Every pair sum of `l31` is a listed modulus (Bool `contains` sweep). -/
theorem sums_covered31 : ∀ vi ∈ l31, ∀ vj ∈ l31, (vi + vj) ∈ moduli31 := by
  have h : (l31.all fun vi => l31.all fun vj => moduli31.contains (vi+vj)) = true := by
    decide
  simp only [List.all_eq_true] at h
  intro vi hvi vj hvj
  exact List.contains_iff_mem.mp (h vi hvi vj hvj)

/-- The full checker passes on the 4/127 member — assembled, not monolithic. -/
theorem check_4_127 : certCheck l31 4 127 = true := by
  simp only [certCheck, List.all_eq_true]
  intro vi hvi vj hvj
  exact chk127_mem _ (sums_covered31 vi hvi vj hvj)

/-- The e-channel certificate of the 4/127 member. -/
theorem cert_4_127 : Cert v31 4 127 := by
  have h := cert_of_check l31 4 127 check_4_127 v31
    (by decide) (by decide) (by decide)
  exact_mod_cast h

/-- **The first kernel-exact TOWER member.**
`sSup (margin v31 '' [0,1]) = 4/127`: the exact loneliness supremum of
`{1,…,29, 31, 120}` — THM-1285's D=4 gap member — machine-checked end to
end.  Floor: the closed-form witness `t₀ = 55/127` element-by-element through
`rung_floor_single`; ceiling: the per-modulus e-channel certificate through
the formal Pinch chain. -/
theorem member_4_127_exact :
    sSup (margin v31 '' Set.Icc (0:ℝ) 1) = (4 : ℝ) / (127 : ℝ) := by
  have h4 : ((4:ℤ) : ℝ) = (4:ℝ) := by norm_num
  have h127 : ((127:ℤ) : ℝ) = (127:ℝ) := by norm_num
  rw [← h4, ← h127]
  apply margin_sSup_eq_of_cert v31 (by decide) 4 127 (by norm_num) (by norm_num)
    ((55 : ℝ) / 127) (by constructor <;> norm_num)
  · intro i m
    have hband : 4 ≤ (v31 i * 55) % 127 ∧ (v31 i * 55) % 127 ≤ 127 - 4 := by
      fin_cases i <;> decide
    have h := RungFloor.rung_floor_single 127 4 (v31 i) 55 (by norm_num) hband m
    have hc : ((127 : ℤ) : ℝ) = (127 : ℝ) := by norm_num
    rw [hc] at h
    convert h using 2 <;> norm_num
  · exact cert_4_127

end EChannelCert
end LonelyRunner

#print axioms LonelyRunner.EChannelCert.margin_le_of_cert
#print axioms LonelyRunner.EChannelCert.margin_sSup_eq_of_cert
#print axioms LonelyRunner.EChannelCert.cert_of_check
#print axioms LonelyRunner.EChannelCert.check_3_23
#print axioms LonelyRunner.EChannelCert.cert_3_23
#print axioms LonelyRunner.EChannelCert.member_3_23_exact
#print axioms LonelyRunner.EChannelCert.check_4_127
#print axioms LonelyRunner.EChannelCert.cert_4_127
#print axioms LonelyRunner.EChannelCert.member_4_127_exact
