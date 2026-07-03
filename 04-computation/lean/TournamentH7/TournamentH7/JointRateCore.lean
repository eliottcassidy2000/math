/-
  TournamentH7.JointRateCore  (mac-mini-2026-07-02-S17, HYP-3874)

  THE j-ARC JOINT rate_core ENGINE — the sharp-F3 middle-band write-out (the program's
  last Tier-0 item, klein-S106 HYP-4011).  The single-comb `rate_core` (kps,
  LRCFarElementRate.lean) bounds ONE interval against ONE ρ-comb with error `2ρ` (two
  boundary cells).  F3-sharp needs the JOINT version: `j` far runners `N + c_i`, error
  `2 c_B (j+1)/N` — crucially `Δ`-FREE (no dependence on the offset spread `Δ = max c_i`).

  THE MECHANISM (mac-mini-S17, sharpened): the `Δ`-free bound is NOT a bounded-variation /
  Riemann-sum bound (that gives the weaker `O(Δ·j/N)`, since `TV(D_c) = O(Δ·j)`).  It is a
  SIGNED TELESCOPING: the per-fast-cell drift residual has the form `h (k+1) − h k`, so the
  interior contributions CANCEL and only `O(j+1)` boundary terms survive — independent of `Δ`.

  This file formalizes the ENGINE of that telescoping, SORRY-FREE, and reduces the joint
  `rate_core` to a single crisp per-cell obligation (`joint_rate_core_reduced`), exactly the
  way the single-comb `rate_core` is `interior-exact + two boundary cells`.
-/
import Mathlib.Tactic

namespace TournamentH7.JointRateCore

/-! ## List helpers (explicit inductions — no fragile API guesses) -/

/-- `Σ (p − q) = Σ p − Σ q` over a mapped list. -/
private theorem sum_map_sub (l : List ℕ) (p q : ℕ → ℚ) :
    (l.map fun k => p k - q k).sum = (l.map p).sum - (l.map q).sum := by
  induction l with
  | nil => simp
  | cons a t ih => simp only [List.map_cons, List.sum_cons, ih]; ring

/-- `|Σ| ≤ Σ|·|` for a rational list. -/
private theorem abs_sum_le (l : List ℚ) : |l.sum| ≤ (l.map abs).sum := by
  induction l with
  | nil => simp
  | cons a t ih =>
      simp only [List.sum_cons, List.map_cons]
      exact (abs_add_le a t.sum).trans (by linarith)

/-- Monotonicity of mapped sums: pointwise `≤` on the index list gives `≤` on the sums. -/
private theorem map_sum_le (l : List ℕ) (p q : ℕ → ℚ)
    (h : ∀ x ∈ l, p x ≤ q x) : (l.map p).sum ≤ (l.map q).sum := by
  induction l with
  | nil => simp
  | cons a t ih =>
      simp only [List.map_cons, List.sum_cons]
      have ha : p a ≤ q a := h a (List.mem_cons_self ..)
      have ht : (t.map p).sum ≤ (t.map q).sum :=
        ih (fun x hx => h x (List.mem_cons_of_mem a hx))
      linarith

/-! ## The telescoping engine (the source of the `Δ`-free bound) -/

/-- **Telescoping sum over a range.**  `Σ_{k<w} (h (k+1) − h k) = h w − h 0`.  This is the
identity whose SIGNED cancellation makes the joint rate bound independent of the offset spread
`Δ` — the interior of every fast cell cancels against its neighbour. -/
theorem sum_range_telescope (h : ℕ → ℚ) (w : ℕ) :
    ((List.range w).map fun k => h (k + 1) - h k).sum = h w - h 0 := by
  induction w with
  | zero => simp
  | succ n ih =>
      rw [List.range_succ, List.map_append, List.sum_append, ih]
      simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
      ring

/-- **The joint-rate telescoping bound (the ENGINE).**  If each cell term `f k` equals the
telescoping increment `h (k+1) − h k` up to a per-cell error `e k`, then the total `Σ f`
equals the boundary value `h w − h 0` up to `Σ e`.  The `h`-increments cancel EXACTLY, so only
the accumulated per-cell errors `e` survive — in the joint rate_core these are supported on the
`O(j)` boundary-crossing cells, NOT the total variation.  This is why the bound is `Δ`-free. -/
theorem telescope_error_bound (f h e : ℕ → ℚ) (w : ℕ)
    (hcell : ∀ k ∈ List.range w, |f k - (h (k + 1) - h k)| ≤ e k) :
    |((List.range w).map f).sum - (h w - h 0)| ≤ ((List.range w).map e).sum := by
  rw [← sum_range_telescope h w, ← sum_map_sub]
  calc |((List.range w).map fun k => f k - (h (k + 1) - h k)).sum|
      ≤ ((List.range w).map fun k => |f k - (h (k + 1) - h k)|).sum := by
        have := abs_sum_le ((List.range w).map fun k => f k - (h (k + 1) - h k))
        rwa [List.map_map] at this
    _ ≤ ((List.range w).map e).sum :=
        map_sum_le _ _ _ (fun x hx => hcell x hx)

/-! ## The joint rate_core, reduced to its per-cell obligation -/

/-- **The joint `rate_core`, reduced.**  Package of the engine: given the far-safe fast-cell
measures `f k`, the drift/coverage function `h` (whose boundary values `h w − h 0` carry the
`D_c`-integral), and per-cell errors `e k` supported on the `O(j)` boundary-crossing cells with
`Σ e ≤ 2 (J+1) / N`, the joint measure `Σ f` equals the boundary value up to `2 (J+1)/N`.

The two hypotheses are EXACTLY the sharp-F3 mechanism:
 * `hcell` — the per-fast-cell residual is a telescoping increment up to a boundary error (the
   generalization of the single-comb `rate_core`'s `interior-exact` cells);
 * `herr`  — those boundary errors are supported on `O(j)` cells and sum to `2(J+1)/N` (the
   `2(j+1)` boundary curves: `2j` arc-edges + the `2 c_B` component ends), `Δ`-free.
This is the crisp remaining Lean obligation for the intermediate band; the engine discharges
everything downstream of it. -/
theorem joint_rate_core_reduced (f h e : ℕ → ℚ) (w : ℕ) (J N : ℚ)
    (hcell : ∀ k ∈ List.range w, |f k - (h (k + 1) - h k)| ≤ e k)
    (herr : ((List.range w).map e).sum ≤ 2 * (J + 1) / N) :
    |((List.range w).map f).sum - (h w - h 0)| ≤ 2 * (J + 1) / N :=
  (telescope_error_bound f h e w hcell).trans herr

end TournamentH7.JointRateCore
