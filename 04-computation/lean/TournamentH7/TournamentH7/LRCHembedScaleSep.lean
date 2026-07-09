/-
  TournamentH7.LRCHembedScaleSep — hembed via the cluster-absorption embedding, with the `e = Vmax − v`
  binding (kind-pasteur-2026-07-09-S106).

  `ScaleSeparation.scale_separation_phase` (THM-608, opus-S50) is a sorry-free ruler embedding: given a
  base `R` safe at `t0` and a cluster `C` near `N` with bounded phase spread `Δφ` and drift `Dd`, it
  produces a real time `t` at which every runner is `1/14`-safe.  This file INSTANTIATES it for the
  LRC(14) reach conclusion, supplying the missing binding: the 13 speeds `v` are the cluster `C`, the
  ruler is `N = Vmax`, and the co-offsets `e_i = Vmax − v_i` enter as the cluster's drift/phase data
  `(v_i − Vmax) = −e_i`.  The base is empty (a single-scale cluster).  The conclusion `∀ c∈C, ∀ m,
  1/14 ≤ |c·t − m|` is exactly `minReach v t ≥ 1/14`, closing to `1/14 ≤ Mreach v` (kps-S99b
  `Mreach_ge_of_lonely_instant`).

  This discharges `hembed` in the DRIFT-CONTROLLED regime `Δφ + Dd·(δ/Vmax) ≤ 3/7` — the large-ruler /
  small-relative-spread case (the cluster-absorption / near-AP branch), where the slow-fast drift
  `spread/Vmax` is small.  (The good-period *window* `Vmax ≈ 7·spread/6` saturates the `3/7` budget and
  needs the density/equidistribution route instead — kps-S105 reflection.)  Self-contained on
  `LRCScaleSeparation` + `LRCReachWitness`.
-/
import Mathlib
import TournamentH7.LRCScaleSeparation
import TournamentH7.LRCReachWitness

namespace LonelyRunner
namespace LRC14Concrete

/-- **`∀ m, c ≤ |x − m|` ⟹ `c ≤ nearInt x`.**  If every integer is `≥ c` from `x`, the nearest one is,
so the nearest-integer distance is `≥ c`.  (The two nearest integers are `⌊x⌋` and `⌊x⌋+1`.) -/
theorem le_nearInt_of_forall_int (c x : ℝ) (h : ∀ m : ℤ, c ≤ |x - (m : ℝ)|) :
    c ≤ nearInt x := by
  have hfr : x - (⌊x⌋ : ℝ) = Int.fract x := by linarith [Int.floor_add_fract x]
  have h1 : c ≤ Int.fract x := by
    have hm := h ⌊x⌋
    rwa [hfr, abs_of_nonneg (Int.fract_nonneg x)] at hm
  have h2 : c ≤ 1 - Int.fract x := by
    have hm := h (⌊x⌋ + 1)
    have he : x - (((⌊x⌋ + 1 : ℤ)) : ℝ) = Int.fract x - 1 := by push_cast; linarith [hfr]
    rw [he, abs_of_nonpos (by linarith [Int.fract_lt_one x])] at hm
    linarith
  exact le_min h1 h2

/-- **hembed via cluster absorption (with the `e = Vmax − v` binding).**  For the 13 speeds `v`, ruler
`Vmax > 0`, a slow time `t0`, slack `δ`, phase spread `Δφ` and drift `Dd`, IF the co-offsets are
bounded (`|v_i − Vmax| ≤ Dd`, i.e. `e_i ≤ Dd`), their phases cluster at `t0` (`|(v_i − Vmax)·t0 − k| ≤
Δφ` — the teeth clearance), and the size conditions `Vmax ≤ 2δ·Vmax`, `Δφ + Dd·(δ/Vmax) ≤ 3/7` hold,
THEN `Mreach v ≥ 1/14`: the family is lonely.  This is `hembed` in the drift-controlled regime. -/
theorem mreach_ge_via_scale_separation
    (v : Fin 13 → ℤ) (Vmax : ℤ) (hVmax : 0 < Vmax)
    (t0 δ Δφ Dd : ℝ) (hδ : 0 < δ) (hΔφ : 0 ≤ Δφ) (hDd : 0 ≤ Dd)
    (hbind : ∀ i, |((v i : ℝ)) - (Vmax : ℝ)| ≤ Dd)
    (hclust : ∀ i, ∃ k : ℤ, |(((v i : ℝ)) - (Vmax : ℝ)) * t0 - (k : ℝ)| ≤ Δφ)
    (hi : (Vmax : ℝ) ≤ 2 * δ * (Vmax : ℝ))
    (hii : Δφ + Dd * (δ / (Vmax : ℝ)) ≤ 3 / 7) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hVmax
  obtain ⟨t, -, hC⟩ :=
    ScaleSeparation.scale_separation_phase t0 δ (Vmax : ℝ) hδ hVR
      ([] : List ℤ) (by simp) (by simp)
      Vmax hVmax Δφ Dd hΔφ hDd
      (List.ofFn v)
      (fun c hc => by
        obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hc
        exact hclust i)
      (fun c hc => by
        obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hc
        exact hbind i)
      hi hii
  apply Mreach_ge_of_lonely_instant
  refine ⟨t, fun i => ?_⟩
  have hmem : (v i) ∈ List.ofFn v := List.mem_ofFn.mpr ⟨i, rfl⟩
  exact le_nearInt_of_forall_int _ _ (fun m => hC (v i) hmem m)

end LRC14Concrete
end LonelyRunner
