/-
  TournamentH7.LRCFarPeelGood — THE FAR-PEEL POSITIVITY on the ACTUAL good region
  (kind-pasteur-2026-07-03-S30).

  opus-S49 (LRCFarPeelCore) proved far-peel positivity in the single-`comb`
  formulation; mac-mini/kps carry the good region as `goodRegion2` (the wrapped
  danger, `dangerPair`).  This file lands the positivity DIRECTLY on `goodRegion2`,
  via the already-proved `damped_peel` (the wrapped-`dangerPair` rate) — so no
  comb/dangerPair bridge is needed.  This is opus's steps 2+3 fused in the right
  formulation:

      length (goodRegion2 (E ++ [w]) h) > 0

  once the base has positive length (step 1 — mac-mini's floor) and the far speed `w`
  clears the peel threshold `(#pieces)·4h < (1−2h)·length·w`.  With
  `exists_lonely_of_goodRegion2_pos` (step 4, kps-S14) this is the far-element peel of
  `CoveringFarLonely 22` for large `w`, minus the finite small-`w` window (step 5).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import TournamentH7.LRCPeelAssembly

namespace LonelyRunner
namespace RatIntervals

/-- **FAR-PEEL POSITIVITY (good-region form).**  If the base good region has positive
length and the far speed `w` clears the peel threshold
`(#pieces)·4h < (1−2h)·length·w`, the peeled good region `goodRegion2 (E ++ [w])`
still has positive length.  Direct consequence of `damped_peel` (the wrapped
`dangerPair` rate). -/
theorem goodRegion2_peel_pos (E : List ℤ) (w : ℤ) (hw : 0 < w) {h : ℚ}
    (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hbig : ((goodRegion2 E h).length : ℚ) * (4 * h)
      < (1 - 2 * h) * length (goodRegion2 E h) * (w.toNat : ℚ)) :
    0 < length (goodRegion2 (E ++ [w]) h) := by
  have hd := damped_peel E w hw hh0 hh1
  have hwnat : 0 < w.toNat := by omega
  set G : ℚ := length (goodRegion2 E h) with hG
  set c : ℚ := ((goodRegion2 E h).length : ℚ) with hc
  set W : ℚ := (w.toNat : ℚ) with hW
  have hWpos : (0 : ℚ) < W := by rw [hW]; exact_mod_cast hwnat
  -- hd : (1 − 2h)·G − c·(4·(h/W)) ≤ length (good (E ++ [w]))
  -- hbig : c·(4h) < (1 − 2h)·G·W
  have hthr : c * (4 * h / W) < (1 - 2 * h) * G := by
    have he : c * (4 * h / W) = c * (4 * h) / W := by ring
    rw [he, div_lt_iff₀ hWpos]
    linarith [hbig]
  linarith [hd, hthr]

/-- **FAR-PEEL POSITIVITY, division-cleared integer threshold.**  Same, with the
threshold as `(#pieces)·4h < (1−2h)·length·w` — the integer-friendly `w`-bound for
certifying covering-far families. -/
theorem goodRegion2_peel_pos_of_gt (E : List ℤ) (w : ℤ) (hw : 0 < w) {h : ℚ}
    (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hlen : 0 < length (goodRegion2 E h))
    (hwbig : ((goodRegion2 E h).length : ℚ) * (4 * h)
      < (1 - 2 * h) * length (goodRegion2 E h) * (w.toNat : ℚ)) :
    0 < length (goodRegion2 (E ++ [w]) h) :=
  goodRegion2_peel_pos E w hw hh0 hh1 hwbig

/-- **THE FAR-ELEMENT PEEL CLOSER** (opus steps 2+3+4 composed): a positive-speed
`13`-family whose last runner `w = v (Fin.last 12)` clears the peel threshold against
the base of the first twelve is lonely.  Reduces `CoveringFarLonely` for a large far
runner to the ONE remaining lemma — the base good-region floor (step 1, from
`LRC(≤13)`), which makes the threshold `hbig` satisfiable.  (No separate positivity of
the base is needed: `hbig` forces it, since a zero-length base makes the strict
inequality `(#pieces)·4h < 0` impossible.) -/
theorem far_peel_lonely (v : Fin 13 → ℤ) (hv : ∀ i, 0 < v i)
    (hbig : ((goodRegion2 (List.ofFn (Fin.init v)) (1 / 14)).length : ℚ) * (4 * (1 / 14))
      < (1 - 2 * (1 / 14)) * length (goodRegion2 (List.ofFn (Fin.init v)) (1 / 14))
          * ((v (Fin.last 12)).toNat : ℚ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hw : 0 < v (Fin.last 12) := hv _
  have hpeel := goodRegion2_peel_pos (List.ofFn (Fin.init v)) (v (Fin.last 12)) hw
    (by norm_num) (by norm_num) hbig
  have heq : List.ofFn v = List.ofFn (Fin.init v) ++ [v (Fin.last 12)] := by
    rw [List.ofFn_succ', List.concat_eq_append]
    rfl
  rw [← heq] at hpeel
  exact exists_lonely_of_goodRegion2_pos v hv hpeel

/-- **A member forces positive length** (reverse of `exists_mem_of_length_pos`): if a
rational `x` lies in a region `L`, then `L` has positive length — because the piece
containing `x` has `p.1 ≤ x < p.2`, hence is non-degenerate.  This reduces the base
good-region floor (step 1) to producing ONE rational good point. -/
theorem length_pos_of_mem {L : Region} {x : ℚ} (hx : mem x L) : 0 < length L := by
  obtain ⟨p, hpL, hp1, hp2⟩ := hx
  have hpos : (0 : ℚ) < p.2 - p.1 := by linarith
  unfold length
  have hnn : ∀ y ∈ L.map fun q => max 0 (q.2 - q.1), (0 : ℚ) ≤ y := by
    intro y hy
    rw [List.mem_map] at hy
    obtain ⟨q, _, rfl⟩ := hy
    exact le_max_left _ _
  have hmem : (max 0 (p.2 - p.1)) ∈ L.map fun q => max 0 (q.2 - q.1) :=
    List.mem_map_of_mem hpL
  have hle := List.single_le_sum hnn _ hmem
  have heq : max 0 (p.2 - p.1) = p.2 - p.1 := max_eq_right hpos.le
  rw [heq] at hle
  linarith

end RatIntervals
end LonelyRunner
