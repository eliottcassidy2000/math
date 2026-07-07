/-
  TournamentH7.LRCDecorrelation  (mac-mini-2026-07-06-S38)

  THE DECORRELATION ATOM for the escape families — a height-descent handle on the
  crux (C)'s obstruction (the compressed `≡ AP mod L` families, mac-mini-S36).

  Owner's inspiration: for `V` with `vᵢ = bᵢ + L·kᵢ` (a bounded base `bᵢ` plus an
  `L`-scaled lift `kᵢ`), the coarse lift dominates:

        M(V) ≥ M(K) − B/L,     K = the lift family (i ↦ kᵢ),  B = maxᵢ |bᵢ|,

  by the reverse triangle inequality at the witness `t = t_K / L` (`distZ` is
  1-Lipschitz). Two consequences for an escape family (`bᵢ ∈ {1,…,12}`, `B = 12`,
  `kᵢ ≥ 1`):
   • **r = #distinct k ≤ 11** (a repeated lift = a close pair): `K` has ≤ 11 speeds,
     so LRC(≤12) gives `M(K) ≥ 1/12`, hence `M(V) ≥ 1/12 − 12/L > 2/25` for
     `L > 3600` — LOOSE, no recursion. (`1/12 − 12/3600 = 2/25` exactly.)
   • **r = 12** (all distinct): `K` is a smaller 12-family (height ÷ L) — descend.

  This file formalizes the atom (`reach_decorr`) and the `r ≤ 11` looseness
  (`escape_loose_of_lift_floor`), kernel-pure. `M(K) ≥ 1/12` enters as a cited
  hypothesis (LRC(≤12), settled) — exactly like LRC(≤13).
-/
import TournamentH7.LRCLadderD1

open TournamentH7.LRCWitness

namespace TournamentH7.LRCDecorr

/-- `distZ` is 1-Lipschitz: `distZ y − |x − y| ≤ distZ x`. -/
lemma distZ_lipschitz (x y : ℝ) : distZ y - |x - y| ≤ distZ x := by
  have h : distZ y ≤ distZ x + dist y x :=
    Metric.infDist_le_infDist_add_dist
  rw [Real.dist_eq, abs_sub_comm y x] at h
  linarith

/-- **Per-time decorrelation.**  If `Vᵢ = bᵢ + L·Kᵢ` with `|bᵢ| ≤ B`, `0 < L`, and
`|t| ≤ 1`, then `margin K t − B/L ≤ margin V (t/L)`. -/
theorem margin_decorr (V K b : Fin 12 → ℤ) (L B : ℤ) (hL : 0 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B) (t : ℝ) (ht : |t| ≤ 1) :
    margin K t - (B : ℝ) / L ≤ margin V (t / L) := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  rw [le_margin_iff]
  intro i m
  -- distZ (V i · t/L) ≥ distZ (K i · t) − B/L ≥ margin K t − B/L
  have hcast : ((V i : ℝ)) * (t / L) = (K i : ℝ) * t + (b i : ℝ) * t / L := by
    have hc : ((V i : ℤ) : ℝ) = (b i : ℝ) + (L : ℝ) * (K i : ℝ) := by
      rw [hV i]; push_cast; ring
    rw [hc]; field_simp; ring
  have hdiff : |((V i : ℝ)) * (t / L) - (K i : ℝ) * t| ≤ (B : ℝ) / L := by
    rw [hcast, show (K i : ℝ) * t + (b i : ℝ) * t / L - (K i : ℝ) * t = (b i : ℝ) * t / L from by ring,
      abs_div, abs_mul, abs_of_pos hLR]
    have hbB : |(b i : ℝ)| ≤ (B : ℝ) := by rw [← Int.cast_abs]; exact_mod_cast hb i
    have h1 : |(b i : ℝ)| * |t| ≤ (B : ℝ) :=
      le_trans (mul_le_of_le_one_right (abs_nonneg _) ht) hbB
    have hLinv : (0 : ℝ) ≤ (L : ℝ)⁻¹ := by positivity
    rw [div_eq_mul_inv, div_eq_mul_inv]
    exact mul_le_mul_of_nonneg_right h1 hLinv
  have hlip : distZ ((K i : ℝ) * t) - |((V i : ℝ)) * (t / L) - (K i : ℝ) * t|
      ≤ distZ ((V i : ℝ) * (t / L)) := distZ_lipschitz _ _
  have hmarginK : margin K t ≤ distZ ((K i : ℝ) * t) :=
    Finset.inf'_le _ (Finset.mem_univ i)
  have key_i : margin K t - (B : ℝ) / L ≤ distZ ((V i : ℝ) * (t / L)) := by linarith
  exact key_i.trans ((le_distZ_iff (distZ ((V i : ℝ) * (t / L))) _).1 le_rfl m)

/-- **The decorrelation atom (reach form).**  `Vᵢ = bᵢ + L·Kᵢ`, `|bᵢ| ≤ B`, `0 < L`
⟹ `reach(K) − B/L ≤ reach(V)`, where `reach X = sSup (margin X '' [0,1])`. -/
theorem reach_decorr (V K b : Fin 12 → ℤ) (L B : ℤ) (hL : 0 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B) :
    sSup (margin K '' Set.Icc (0 : ℝ) 1) - (B : ℝ) / L
      ≤ sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  obtain ⟨t0, ht0mem, ht0max⟩ := exists_max_margin K
  -- reach K = margin K t0 (max attained)
  have hbddK : BddAbove (margin K '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨s, _, rfl⟩; exact LRCLadderD1.margin_le_half K s⟩
  have hsupK : sSup (margin K '' Set.Icc (0 : ℝ) 1) = margin K t0 := by
    apply le_antisymm
    · refine csSup_le ((show (Set.Icc (0 : ℝ) 1).Nonempty from ⟨t0, ht0mem⟩).image (margin K)) ?_
      rintro _ ⟨s, hs, rfl⟩; exact ht0max s hs
    · exact le_csSup hbddK ⟨t0, ht0mem, rfl⟩
  -- t0/L ∈ [0,1] and |t0| ≤ 1
  have ht0_01 : t0 ∈ Set.Icc (0 : ℝ) 1 := ht0mem
  have ht0abs : |t0| ≤ 1 := by rw [abs_le]; exact ⟨by linarith [ht0_01.1], ht0_01.2⟩
  have hL1 : (1 : ℝ) ≤ (L : ℝ) := by exact_mod_cast hL
  have htL_01 : t0 / L ∈ Set.Icc (0 : ℝ) 1 := by
    refine Set.mem_Icc.mpr ⟨div_nonneg (by linarith [ht0_01.1]) (le_of_lt hLR), ?_⟩
    rw [div_le_one hLR]; linarith [ht0_01.2]
  -- reach V ≥ margin V (t0/L) ≥ margin K t0 − B/L = reach K − B/L
  have hbddV : BddAbove (margin V '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨s, _, rfl⟩; exact LRCLadderD1.margin_le_half V s⟩
  have hmemV : margin V (t0 / L) ∈ margin V '' Set.Icc (0 : ℝ) 1 := ⟨t0 / L, htL_01, rfl⟩
  have hle : margin V (t0 / L) ≤ sSup (margin V '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbddV hmemV
  have hdec := margin_decorr V K b L B hL hV hb t0 ht0abs
  rw [hsupK]; linarith

/-- **Escape-family looseness, `r ≤ 11` (close-pair case).**  If `Vᵢ = bᵢ + L·Kᵢ`
with `|bᵢ| ≤ 12`, `0 < L`, and the lift family `K` clears the LRC(≤12) floor
`1/12 ≤ reach K` (settled: `K` has ≤ 11 distinct nonzero speeds), then for `L >
3600`, `reach V > 2/25` — the family is loose, not in the gap. -/
theorem escape_loose_of_lift_floor (V K b : Fin 12 → ℤ) (L : ℤ) (hL : 3600 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ 12)
    (hKfloor : (1 : ℝ) / 12 ≤ sSup (margin K '' Set.Icc (0 : ℝ) 1)) :
    (2 : ℝ) / 25 < sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  have hL0 : 0 < L := by omega
  have hdec := reach_decorr V K b L 12 hL0 hV hb
  have hLR : (3600 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  have hbnd : (12 : ℝ) / L < (1 : ℝ) / 12 - 2 / 25 := by
    rw [div_lt_iff₀ (by linarith : (0:ℝ) < (L:ℝ))]
    nlinarith
  linarith

#print axioms distZ_lipschitz
#print axioms margin_decorr
#print axioms reach_decorr
#print axioms escape_loose_of_lift_floor

end TournamentH7.LRCDecorr
