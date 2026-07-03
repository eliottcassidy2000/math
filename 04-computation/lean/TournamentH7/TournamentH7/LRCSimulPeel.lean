/-
  TournamentH7.LRCSimulPeel  (mac-mini-2026-07-02-S18, HYP-3875)

  THE SIMULTANEOUS UNION-BOUND FAR-PEEL — opus-S32's simultaneous-peel lemma (HYP-3900),
  the tractable middle-band closer.  Where `damped_peel` (LRCPeelAssembly) peels ONE far
  runner and re-fragments the region (the `#pieces` fee compounds), this peels ALL far
  runners from the FIXED window good region in ONE `length_diffF_ge` — a union bound, no
  iteration, no fragmentation, no blowup:

      length (good (B ++ far))  ≥  (1 − 2h·j)·length (good B)  −  c_B·4h·Σ_{w∈far} 1/w ,

  where `j = |far|`, `c_B = #pieces of good B`.  The floor `(1 − 2h·j) = (1 − j/7)` at
  `h = 1/14` is POSITIVE for `j ≤ 6`, and the fee `→ 0` as the far speeds `→ ∞`, so this
  closes the intermediate band `22 < N < N*` for `≤ 6` far elements by a FINITE sweep —
  WITHOUT the sharp joint rate_core telescoping.  Kernel-pure, built on landed lemmas.
-/
import TournamentH7.LRCPeelAssembly

namespace LonelyRunner
namespace RatIntervals

/-- **Peeling a whole far-list at once.**  `good (B ++ far) = good B  with all far dangers
subtracted` (the list generalization of `goodRegion2_append`). -/
theorem goodRegion2_append_list (B far : List ℤ) (h : ℚ) :
    goodRegion2 (B ++ far) h
      = diffF (goodRegion2 B h) (far.flatMap fun s => dangerPair s.toNat h) := by
  unfold goodRegion2 diffF
  rw [List.flatMap_append, List.foldl_append]

/-- **`length ∘ inter` distributes over `flatMap`** (danger regions concatenate; `length`
adds over `++`, so the intersection length is the sum of per-element intersection lengths). -/
theorem length_inter_flatMap (L : Region) (l : List ℤ) (f : ℤ → Region) :
    length (inter L (l.flatMap f)) = (l.map fun s => length (inter L (f s))).sum := by
  induction l with
  | nil => simp only [List.flatMap_nil, List.map_nil, List.sum_nil, length_inter_nil]
  | cons a t ih =>
      rw [List.flatMap_cons, length_inter_append_right, ih, List.map_cons, List.sum_cons]

/-- **The per-runner danger bound on a fixed region** (extracted from `damped_peel`'s core):
one speed's danger pair clips at most density `2h` of `G` plus `c_G·4h/w`. -/
theorem length_inter_dangerPair_le (G : Region)
    (hG : ∀ p ∈ G, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1)
    (w : ℕ) (hw : 0 < w) {h : ℚ} (hh0 : 0 < h) (hh1 : h ≤ 1 / 2) :
    length (inter G (dangerPair w h)) ≤ 2 * h * length G + (G.length : ℚ) * (4 * h / w) := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  have hsplit : length (inter G (dangerPair w h))
      = length (inter G (comb w (h / 2) (h / 2)))
        + length (inter G (comb w (h / 2) (1 - h / 2))) := by
    unfold dangerPair; exact length_inter_append_right G _ _
  have hrate1 := length_inter_comb_near_region (w := w) hw (r := h / 2) (φ := h / 2)
    (by linarith) (le_refl _) (by linarith) G hG
  have hrate2 := length_inter_comb_near_region (w := w) hw (r := h / 2) (φ := 1 - h / 2)
    (by linarith) (by linarith) (by linarith) G hG
  have habs1 := (abs_le.mp hrate1).2
  have habs2 := (abs_le.mp hrate2).2
  have e1 : 2 * (h / 2) * length G + (G.length : ℚ) * (4 * (h / 2) / (w : ℚ))
      = h * length G + (G.length : ℚ) * (2 * h / (w : ℚ)) := by ring
  have hb1 : length (inter G (comb w (h / 2) (h / 2)))
      ≤ h * length G + (G.length : ℚ) * (2 * h / (w : ℚ)) := by rw [← e1]; linarith
  have hb2 : length (inter G (comb w (h / 2) (1 - h / 2)))
      ≤ h * length G + (G.length : ℚ) * (2 * h / (w : ℚ)) := by rw [← e1]; linarith
  have hcomb : (G.length : ℚ) * (2 * h / (w : ℚ)) * 2 = (G.length : ℚ) * (4 * h / (w : ℚ)) := by
    ring
  rw [hsplit]; linarith

/-- Every piece of a speed's danger pair is a live (non-empty-width) interval. -/
theorem dangerPair_live (w : ℕ) (hw : 0 < w) {h : ℚ} (hh0 : 0 ≤ h) :
    ∀ q ∈ dangerPair w h, q.1 ≤ q.2 := by
  intro q hq
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  unfold dangerPair at hq
  rcases List.mem_append.mp hq with hq' | hq' <;>
  · unfold comb at hq'
    simp only [List.map_map, List.mem_map, List.mem_range, Function.comp_apply] at hq'
    obtain ⟨k, _, rfl⟩ := hq'
    simp only
    apply div_le_div_of_nonneg_right ?_ hwQ.le
    linarith

/-- Pointwise `≤` on a list gives `≤` on the mapped sums. -/
private theorem map_sum_le' (l : List ℤ) (p q : ℤ → ℚ) (h : ∀ x ∈ l, p x ≤ q x) :
    (l.map p).sum ≤ (l.map q).sum := by
  induction l with
  | nil => simp
  | cons a t ih =>
      simp only [List.map_cons, List.sum_cons]
      have ha := h a (List.mem_cons_self ..)
      have ih' := ih (fun x hx => h x (List.mem_cons_of_mem _ hx))
      linarith

/-- The far-danger fee sum evaluates: `Σ (A + D·4h/wᵢ) = |far|·A + D·4h·Σ(1/wᵢ)`. -/
theorem sum_map_far (far : List ℤ) (A D hh : ℚ) :
    (far.map fun s => A + D * (4 * hh / (s.toNat : ℚ))).sum
      = (far.length : ℚ) * A + D * (4 * hh) * (far.map fun w => 1 / (w.toNat : ℚ)).sum := by
  induction far with
  | nil => simp
  | cons a t ih =>
      simp only [List.map_cons, List.sum_cons, ih, List.length_cons]
      push_cast; ring

/-- **THE SIMULTANEOUS UNION-BOUND FAR-PEEL** (opus-S32's simultaneous-peel lemma, HYP-3900).
Peeling a whole far-list from the fixed window good region loses at most the union of the
per-runner fees — NO iteration, NO fragmentation:

    length (good (B ++ far)) ≥ (1 − 2h·|far|)·length (good B) − c_B·4h·Σ_{w∈far} 1/w.

The floor `(1 − 2h·|far|)` is positive for `|far| < 1/(2h) = 7`, and the fee `→ 0` as the far
speeds grow: this closes the intermediate band for `≤ 6` far elements by a finite sweep. -/
theorem goodRegion2_simul_peel (B far : List ℤ) {h : ℚ} (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hpos : ∀ w ∈ far, 0 < w) :
    (1 - 2 * h * (far.length : ℚ)) * length (goodRegion2 B h)
      - ((goodRegion2 B h).length : ℚ) * (4 * h)
        * (far.map fun w => 1 / (w.toNat : ℚ)).sum
    ≤ length (goodRegion2 (B ++ far) h) := by
  set L := goodRegion2 B h with hL
  have hLbounds : ∀ p ∈ L, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1 := by
    rw [hL]; unfold goodRegion2
    apply diffF_window_bounds
    intro p hp; rcases List.mem_singleton.mp hp with rfl; norm_num
  have hlive : ∀ q ∈ (far.flatMap fun s => dangerPair s.toNat h), q.1 ≤ q.2 := by
    intro q hq
    rw [List.mem_flatMap] at hq
    obtain ⟨s, hs, hqs⟩ := hq
    exact dangerPair_live s.toNat (by have := hpos s hs; omega) hh0.le q hqs
  rw [goodRegion2_append_list]
  have hge := length_diffF_ge L (far.flatMap fun s => dangerPair s.toNat h) hlive
  rw [length_inter_flatMap] at hge
  have hterm : ∀ s ∈ far, length (inter L (dangerPair s.toNat h))
      ≤ 2 * h * length L + (L.length : ℚ) * (4 * h / (s.toNat : ℚ)) := by
    intro s hs
    exact length_inter_dangerPair_le L hLbounds s.toNat
      (by have := hpos s hs; omega) hh0 hh1
  have hmapbound : (far.map fun s => length (inter L (dangerPair s.toNat h))).sum
      ≤ (far.length : ℚ) * (2 * h * length L)
        + (L.length : ℚ) * (4 * h) * (far.map fun w => 1 / (w.toNat : ℚ)).sum := by
    calc (far.map fun s => length (inter L (dangerPair s.toNat h))).sum
        ≤ (far.map fun s => 2 * h * length L + (L.length : ℚ) * (4 * h / (s.toNat : ℚ))).sum :=
          map_sum_le' far _ _ hterm
      _ = _ := sum_map_far far (2 * h * length L) (L.length : ℚ) h
  have hgoal_eq : (1 - 2 * h * (far.length : ℚ)) * length L
        - (L.length : ℚ) * (4 * h) * (far.map fun w => 1 / (w.toNat : ℚ)).sum
      = length L - ((far.length : ℚ) * (2 * h * length L)
        + (L.length : ℚ) * (4 * h) * (far.map fun w => 1 / (w.toNat : ℚ)).sum) := by ring
  rw [hgoal_eq]
  linarith [hge, hmapbound]

end RatIntervals
end LonelyRunner
