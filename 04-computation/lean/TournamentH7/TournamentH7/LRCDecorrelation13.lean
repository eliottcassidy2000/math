/-
  TournamentH7.LRCDecorrelation13  (kind-pasteur-2026-07-11-S127 cont.48)

  THE 13-RUNNER DECORRELATION ATOM — the direct-LRC(14) form of THM-636 (mac-mini-S38,
  12-runner / 2/25-gap). For a 13-speed family `V` with `vᵢ = bᵢ + L·kᵢ` (bounded base
  `|bᵢ| ≤ B` + `L`-scaled lift `kᵢ`):

        reach(V) ≥ reach(K) − B/L,       K = the lift family (i ↦ kᵢ),

  by the reverse triangle inequality at the witness `t = t_K / L` (`distZ` is 1-Lipschitz).
  Consequence (the escape / large-diameter looseness, mac-mini cont.49): a large-diameter
  divisor-complete family collapses at its scale `L` to a lift family `K` with FEW distinct
  speeds, so `LRC(≤13)` (settled) gives a floor `reach(K) ≥ f`, and

        reach(V) ≥ f − B/L > 1/14   once   B/L < f − 1/14.

  With `f = 1/13` (≤ 12 distinct lifts ⟹ LRC(≤13)) and `B = 13`, this needs `L > 2366`;
  with `f = 1/7` (≤ 6 distinct lifts, DC even-heaviness ⟹ LRC(7)) only `L > 182`. So every
  large-diameter DC family is LOOSE (`reach > 1/14`) by decorrelation descent — the large-
  diameter half of the endgame dichotomy; the small-`L` base is the bounded-diameter finite
  check (kps cont.47). The lift floor enters as a cited `LRC(≤13)` hypothesis, exactly as in
  THM-636.
-/
import TournamentH7.LRCDecorrelation

open TournamentH7.LRCWitness TournamentH7.LRCDecorr

namespace TournamentH7.LRCDecorr13

/-- `margin v t ≤ 1/2` for a 13-speed family (generic bound; the `inf'` is `≤` its first
term, and `distZ ≤ 1/2`). -/
lemma margin_le_half13 (v : Fin 13 → ℤ) (t : ℝ) : margin v t ≤ 1 / 2 :=
  le_trans (Finset.inf'_le _ (Finset.mem_univ (0 : Fin 13))) (LRCLadderD1.distZ_le_half _)

/-- **Per-time decorrelation (13 runners).**  `Vᵢ = bᵢ + L·Kᵢ`, `|bᵢ| ≤ B`, `0 < L`, `|t| ≤ 1`
⟹ `margin K t − B/L ≤ margin V (t/L)`. -/
theorem margin_decorr13 (V K b : Fin 13 → ℤ) (L B : ℤ) (hL : 0 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B) (t : ℝ) (ht : |t| ≤ 1) :
    margin K t - (B : ℝ) / L ≤ margin V (t / L) := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  rw [le_margin_iff]
  intro i m
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

/-- **The 13-runner decorrelation atom (reach form).**  `Vᵢ = bᵢ + L·Kᵢ`, `|bᵢ| ≤ B`, `0 < L`
⟹ `reach(K) − B/L ≤ reach(V)`. -/
theorem reach_decorr13 (V K b : Fin 13 → ℤ) (L B : ℤ) (hL : 0 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B) :
    sSup (margin K '' Set.Icc (0 : ℝ) 1) - (B : ℝ) / L
      ≤ sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  obtain ⟨t0, ht0mem, ht0max⟩ := exists_max_margin K
  have hbddK : BddAbove (margin K '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨s, _, rfl⟩; exact margin_le_half13 K s⟩
  have hsupK : sSup (margin K '' Set.Icc (0 : ℝ) 1) = margin K t0 := by
    apply le_antisymm
    · refine csSup_le ((show (Set.Icc (0 : ℝ) 1).Nonempty from ⟨t0, ht0mem⟩).image (margin K)) ?_
      rintro _ ⟨s, hs, rfl⟩; exact ht0max s hs
    · exact le_csSup hbddK ⟨t0, ht0mem, rfl⟩
  have ht0_01 : t0 ∈ Set.Icc (0 : ℝ) 1 := ht0mem
  have ht0abs : |t0| ≤ 1 := by rw [abs_le]; exact ⟨by linarith [ht0_01.1], ht0_01.2⟩
  have htL_01 : t0 / L ∈ Set.Icc (0 : ℝ) 1 := by
    refine Set.mem_Icc.mpr ⟨div_nonneg (by linarith [ht0_01.1]) (le_of_lt hLR), ?_⟩
    rw [div_le_one hLR]; nlinarith [ht0_01.2, (by exact_mod_cast hL : (1:ℝ) ≤ (L:ℝ))]
  have hbddV : BddAbove (margin V '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨s, _, rfl⟩; exact margin_le_half13 V s⟩
  have hmemV : margin V (t0 / L) ∈ margin V '' Set.Icc (0 : ℝ) 1 := ⟨t0 / L, htL_01, rfl⟩
  have hle : margin V (t0 / L) ≤ sSup (margin V '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbddV hmemV
  have hdec := margin_decorr13 V K b L B hL hV hb t0 ht0abs
  rw [hsupK]; linarith

/-- **The escape looseness, general floor.**  If `Vᵢ = bᵢ + L·Kᵢ` with `|bᵢ| ≤ 13`, `0 < L`,
the lift family clears a floor `f ≤ reach K` (from `LRC(≤13)` on `K`'s distinct-speed count),
and `13/L < f − 1/14`, then `reach V > 1/14` — the family is LOOSE. -/
theorem escape_loose13_of_floor (V K b : Fin 13 → ℤ) (L : ℤ) (f : ℝ) (hL : 0 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ 13)
    (hKfloor : f ≤ sSup (margin K '' Set.Icc (0 : ℝ) 1))
    (hmargin : (13 : ℝ) / L < f - 1 / 14) :
    (1 : ℝ) / 14 < sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  have hdec := reach_decorr13 V K b L 13 hL hV hb
  linarith

/-- **Escape looseness, ≤ 12 distinct lifts (`f = 1/13`, `L > 2366`).**  `LRC(≤13)` gives
`reach K ≥ 1/13`; then any `L > 2366` makes `V` loose. -/
theorem escape_loose13_le12 (V K b : Fin 13 → ℤ) (L : ℤ) (hL : 2366 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ 13)
    (hKfloor : (1 : ℝ) / 13 ≤ sSup (margin K '' Set.Icc (0 : ℝ) 1)) :
    (1 : ℝ) / 14 < sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  refine escape_loose13_of_floor V K b L (1 / 13) (by omega) hV hb hKfloor ?_
  have hLR : (2366 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  rw [div_lt_iff₀ (by linarith : (0 : ℝ) < (L : ℝ))]
  nlinarith

/-- **Escape looseness, ≤ 6 distinct lifts (`f = 1/7`, `L > 182`).**  DC even-heaviness
collapses the lifts to ≤ 6 distinct speeds, so `LRC(7)` gives `reach K ≥ 1/7`; then any
`L > 182` makes `V` loose (the sharp DC threshold, mac-mini cont.49). -/
theorem escape_loose13_le6 (V K b : Fin 13 → ℤ) (L : ℤ) (hL : 182 < L)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ 13)
    (hKfloor : (1 : ℝ) / 7 ≤ sSup (margin K '' Set.Icc (0 : ℝ) 1)) :
    (1 : ℝ) / 14 < sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  refine escape_loose13_of_floor V K b L (1 / 7) (by omega) hV hb hKfloor ?_
  have hLR : (182 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  rw [div_lt_iff₀ (by linarith : (0 : ℝ) < (L : ℝ))]
  nlinarith

#print axioms reach_decorr13
#print axioms escape_loose13_le12
#print axioms escape_loose13_le6

end TournamentH7.LRCDecorr13
