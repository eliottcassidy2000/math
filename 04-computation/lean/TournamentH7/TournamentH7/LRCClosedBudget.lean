/-
  TournamentH7.LRCClosedBudget — the ALGEBRAIC SPINE of THM-750 (was THM-748;
  opus-2026-07-14-S283/S284): the closed budget of the tail lane,
  `W (L − Area) = PHI + QPOT + KAP`.

  Every algebraic atom of the identity, sorry-free, organized by the PERSPECTIVE FRAME
  that produced it (the owner's `(n−1)² = spokes ⊕ pairs` decomposition):

    • ORIGIN BAND — `strong_no_wrap`: exposed heights are fenced `≥ 1/14`, and a
      full-inside crossing sits `≥ 1/W` before its segment's death, so `h ≥ 1/14 + j/W`:
      the height march never wraps (S280).
    • SPOKES — `per_crossing_exact`: the strand crossing `σ* = hW/(W+j)` vs the strip
      average `h − j/(2W)`: the difference is `(j/W)ψ(h) + h j²/(W(W+j))`, NO remainder.
    • TELESCOPING — `F_step_nowrap` / `F_step_wrap`: `αψ(h) = −ΔF − α²/2 (+ wrap term)`,
      `F(x) = x(1−x)/2` (THM-745's atom); `F_mirror`, `psi_mirror` (time reversal).
    • PAIR EVENTS — the κ closed forms (death / birth / swap / cut), the `φ = 0`
      boundary degeneracy (`kappaD_at_zero`), and the THM-743-style magnitude bound.
    • MIRROR / GRID — `floor_W_sub`, `ceil_W_sub`, `grid_count_match` (S279's Lemma C)
      and `fract_int_shift` (integer time makes the two lines of an arc exchange).
    • THE BUDGET SPEC — computable ℚ definitions of `PHI`, `QPOT`, `KAP` over explicit
      segment/event data, and the master identity `ClosedBudgetSpec` as a DECIDABLE
      Prop.  Its semantic witness is the exact-rational referee battery
      (lrc14_closed_budget_thm748_opus_S283.py: Fraction equality 6/6, both shapes);
      any Lean-native instance re-verification is a `decide` away.

  No measure theory anywhere; ℝ for the analytic atoms, ℚ for the computable spec.
-/
import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14
namespace ClosedBudget

/-! ### The telescoping potential and the sawtooth -/

/-- `F(x) = x(1−x)/2` — the quadratic telescoping potential (THM-744/745). -/
noncomputable def Fpot (x : ℝ) : ℝ := x * (1 - x) / 2

/-- `ψ(x) = 1/2 − x` — the (linearized) sawtooth. -/
noncomputable def psi (x : ℝ) : ℝ := 1 / 2 - x

/-- Mirror symmetry of the potential: `F(1−x) = F(x)` (time reversal fixes `F`). -/
theorem F_mirror (x : ℝ) : Fpot (1 - x) = Fpot x := by
  unfold Fpot; ring

/-- Mirror antisymmetry of the sawtooth: `ψ(1−x) = −ψ(x)`. -/
theorem psi_mirror (x : ℝ) : psi (1 - x) = -psi x := by
  unfold psi; ring

/-- **THM-745's atom (no-wrap step).**  `α ψ(h) = −(F(h−α) − F(h)) − α²/2` — exact,
because `F'' = −1` and the expansion has no tail. -/
theorem F_step_nowrap (α h : ℝ) :
    α * psi h = -(Fpot (h - α) - Fpot h) - α ^ 2 / 2 := by
  unfold Fpot psi; ring

/-- **THM-745's atom (wrap step).**  When the march wraps (`h < α`, next height
`h − α + 1`), the correction is exactly `α − h`. -/
theorem F_step_wrap (α h : ℝ) :
    α * psi h = -(Fpot (h - α + 1) - Fpot h) - α ^ 2 / 2 + (α - h) := by
  unfold Fpot psi; ring

/-! ### The origin band: the strong no-wrap lemma (S280) -/

/-- **Strong no-wrap (THM-745 §2c).**  If the segment's death height is `≥ 1/14`
(the origin band fences all exposed heights) and a crossing sits at `u`-distance
`≥ 1/W` before the death, its height is `≥ 1/14 + j/W > j/W`: the wrap indicator can
never fire at a full-inside crossing — at ANY `W`. -/
theorem strong_no_wrap {j W h2 d : ℝ} (hj : 0 ≤ j) (_hW : 0 < W)
    (hfence : 1 / 14 ≤ h2) (hdist : 1 / W ≤ d) :
    1 / 14 + j / W ≤ h2 + j * d := by
  have hjd : j / W ≤ j * d := by
    have := mul_le_mul_of_nonneg_left hdist hj
    calc j / W = j * (1 / W) := by ring
    _ ≤ j * d := this
  linarith

/-! ### The spoke sector: the exact per-crossing wedge -/

/-- **The per-crossing formula is EXACT (THM-750 §PHI+QPOT).**  The strand crosses a
slope-`−j` boundary entering at height `h` at `σ* = hW/(W+j)`; the strip average of the
boundary is `h − j/(2W)`; their difference is `(j/W)ψ(h) + h j²/(W(W+j))` with no
remainder. -/
theorem per_crossing_exact {W j : ℝ} (h : ℝ) (hW : W ≠ 0) (hWj : W + j ≠ 0) :
    h * W / (W + j) - (h - j / (2 * W))
      = (j / W) * psi h + h * j ^ 2 / (W * (W + j)) := by
  unfold psi
  field_simp
  ring

/-- The crossing is the fixed point of the drift: `σ* = h − j σ*/W` (self-referential
drift resolved exactly — the S275 "exact fiber" reformulation's kernel). -/
theorem sigma_fixed_point {W j : ℝ} (h : ℝ) (hW : W ≠ 0) (hWj : W + j ≠ 0) :
    h * W / (W + j) = h - j * (h * W / (W + j)) / W := by
  field_simp
  ring

/-! ### The pair sector: the event corrections κ -/

/-- Death correction: the run closing at height `x` with relative speed
`δ = j_up − j_lo`, event phase `φ` in its strip. -/
noncomputable def kappaD (W jup jlo x φ : ℝ) : ℝ :=
  W * (jup - jlo) * max (φ - x) 0 / ((W + jup) * (W + jlo))
    - (jup - jlo) * φ ^ 2 / (2 * W)

/-- Birth correction (the mirror of death). -/
noncomputable def kappaB (W jup jlo x φ : ℝ) : ℝ :=
  W * (jlo - jup) * max (x - φ) 0 / ((W + jup) * (W + jlo))
    - (jlo - jup) * (1 - φ) ^ 2 / (2 * W)

/-- The swap crossing: the kinked path's strand crossing on the branch of slope `−j`. -/
noncomputable def sigmaCross (W j x φ : ℝ) : ℝ := (x * W + j * φ) / (W + j)

/-- **The `φ = 0` boundary degeneracy (the S283 referee's second bug).**  When the event
sits exactly on a strip boundary, the death correction vanishes — the neighbouring
full-inside sums already carry everything. -/
theorem kappaD_at_zero {W jup jlo x : ℝ} (hx : 0 ≤ x) :
    kappaD W jup jlo x 0 = 0 := by
  unfold kappaD
  have hmax : max (0 - x) 0 = 0 := by
    apply max_eq_right; linarith
  rw [hmax]
  ring

/-- **THM-743-style magnitude bound.**  For `0 ≤ j_lo ≤ j_up`, `0 ≤ x`, `0 ≤ φ ≤ 1`:
`|κ_D| ≤ 3 δ/(2W)` with `δ = j_up − j_lo` — the per-event cost is the PAIR DIFFERENCE,
not the worst line. -/
theorem kappaD_bound {W jup jlo x φ : ℝ} (hW : 0 < W) (hlo : 0 ≤ jlo)
    (hj : jlo ≤ jup) (hx : 0 ≤ x) (hφ0 : 0 ≤ φ) (hφ1 : φ ≤ 1) :
    |kappaD W jup jlo x φ| ≤ 3 * (jup - jlo) / (2 * W) := by
  set δ := jup - jlo with hδ
  have hδ0 : 0 ≤ δ := by simp [hδ]; linarith
  have hWup : 0 < W + jup := by linarith
  have hWlo : 0 < W + jlo := by linarith
  have hden : 0 < (W + jup) * (W + jlo) := mul_pos hWup hWlo
  have hmax0 : (0 : ℝ) ≤ max (φ - x) 0 := le_max_right _ _
  have hmax1 : max (φ - x) 0 ≤ 1 := by
    apply max_le _ (by norm_num)
    linarith
  have hfirst0 : 0 ≤ W * δ * max (φ - x) 0 / ((W + jup) * (W + jlo)) := by
    apply div_nonneg _ (le_of_lt hden)
    have := mul_nonneg (mul_nonneg (le_of_lt hW) hδ0) hmax0
    linarith
  have hfirst1 : W * δ * max (φ - x) 0 / ((W + jup) * (W + jlo)) ≤ δ / W := by
    rw [div_le_div_iff₀ hden hW]
    have hWW : W * W ≤ (W + jup) * (W + jlo) := by nlinarith
    calc W * δ * max (φ - x) 0 * W ≤ W * δ * 1 * W := by
          have h1 : (0:ℝ) ≤ W * δ := mul_nonneg (le_of_lt hW) hδ0
          have := mul_le_mul_of_nonneg_left hmax1 h1
          nlinarith
      _ = δ * (W * W) := by ring
      _ ≤ δ * ((W + jup) * (W + jlo)) := by
          exact mul_le_mul_of_nonneg_left hWW hδ0
  have hsecond0 : 0 ≤ δ * φ ^ 2 / (2 * W) := by
    apply div_nonneg _ (by linarith)
    positivity
  have hsecond1 : δ * φ ^ 2 / (2 * W) ≤ δ / (2 * W) := by
    have h2W : (0 : ℝ) < 2 * W := by linarith
    rw [div_le_div_iff₀ h2W h2W]
    have hφsq : φ ^ 2 ≤ 1 := by nlinarith
    have hkey := mul_le_mul_of_nonneg_left hφsq (mul_nonneg hδ0 (le_of_lt h2W))
    nlinarith [hkey]
  unfold kappaD
  rw [← hδ]
  rw [abs_sub_comm]
  calc |δ * φ ^ 2 / (2 * W) - W * δ * max (φ - x) 0 / ((W + jup) * (W + jlo))|
      ≤ |δ * φ ^ 2 / (2 * W)| + |W * δ * max (φ - x) 0 / ((W + jup) * (W + jlo))| :=
        abs_sub _ _
    _ = δ * φ ^ 2 / (2 * W) + W * δ * max (φ - x) 0 / ((W + jup) * (W + jlo)) := by
        rw [abs_of_nonneg hsecond0, abs_of_nonneg hfirst0]
    _ ≤ δ / (2 * W) + δ / W := by linarith
    _ = 3 * δ / (2 * W) := by ring

/-- The swap crossing solves the drift fixed-point on its branch. -/
theorem sigmaCross_fixed {W j : ℝ} (x φ : ℝ) (hW : W ≠ 0) (hWj : W + j ≠ 0) :
    sigmaCross W j x φ = x + j * (φ - sigmaCross W j x φ) / W := by
  unfold sigmaCross
  field_simp
  ring

/-- Branch validity: the swap crossing lands before the vertex iff `x ≤ φ`. -/
theorem sigmaCross_le_iff {W j x φ : ℝ} (hW : 0 < W) (hWj : 0 < W + j) :
    sigmaCross W j x φ ≤ φ ↔ x ≤ φ := by
  unfold sigmaCross
  rw [div_le_iff₀ hWj]
  constructor
  · intro hle
    nlinarith
  · intro hle
    nlinarith

/-! ### The mirror / grid sector (S279's Lemma C; integer time) -/

/-- `⌊W − y⌋ = W − ⌈y⌉` for integer `W`. -/
theorem floor_W_sub (W : ℤ) (y : ℝ) : ⌊(W : ℝ) - y⌋ = W - ⌈y⌉ := by
  rw [show (W : ℝ) - y = -y + (W : ℤ) by ring, Int.floor_add_intCast, Int.floor_neg]
  ring

/-- `⌈W − y⌉ = W − ⌊y⌋` for integer `W`. -/
theorem ceil_W_sub (W : ℤ) (y : ℝ) : ⌈(W : ℝ) - y⌉ = W - ⌊y⌋ := by
  rw [show (W : ℝ) - y = -y + (W : ℤ) by ring, Int.ceil_add_intCast, Int.ceil_neg]
  ring

/-- **Grid-count matching (S279 Lemma C).**  Mirror-paired segments have equal
full-inside crossing counts: `⌊u₂W⌋ − ⌈u₁W⌉ = ⌊(1−u₁)W⌋ − ⌈(1−u₂)W⌉`. -/
theorem grid_count_match (W : ℤ) (u1 u2 : ℝ) :
    ⌊u2 * W⌋ - ⌈u1 * W⌉ = ⌊(1 - u1) * W⌋ - ⌈(1 - u2) * W⌉ := by
  have h1 : (1 - u1) * (W : ℝ) = (W : ℝ) - u1 * W := by ring
  have h2 : (1 - u2) * (W : ℝ) = (W : ℝ) - u2 * W := by ring
  rw [h1, h2, floor_W_sub, ceil_W_sub]
  ring

/-- **Integer time exchanges the two lines of one arc** (`j − 1/14 ≡ −1/14 (mod 1)`):
`Int.fract ((j : ℝ) + x) = Int.fract x` — the S278 mirror-pairing's arithmetic kernel. -/
theorem fract_int_shift (j : ℤ) (x : ℝ) : Int.fract ((j : ℝ) + x) = Int.fract x := by
  rw [add_comm, Int.fract_add_intCast]

/-! ### The computable budget spec over ℚ -/

/-- Segment data: line index `j`, orientation `±1`, endpoints `u1 < u2`, and the
line's height intercept `c` (so the height at `u` is the fractional part of `c − j u`). -/
structure Seg where
  j : ℤ
  orient : ℚ
  u1 : ℚ
  u2 : ℚ
  c : ℚ
deriving Repr

/-- Event data (the pair sector): closed-form κ inputs. -/
inductive Ev where
  | death (jup jlo : ℤ) (x u : ℚ)
  | birth (jup jlo : ℤ) (x u : ℚ)
  | swap (o : ℚ) (jA jB : ℤ) (x u : ℚ)
  | cutR (runs : List (ℚ × ℤ × ℚ × ℤ)) (uG : ℚ)
  | cutL (runs : List (ℚ × ℤ × ℚ × ℤ)) (uG : ℚ)
deriving Repr

/-- Fractional part over ℚ (computable). -/
def fracQ (q : ℚ) : ℚ := q - ⌊q⌋

/-- Per-segment `PHI + QPOT` at integer time `W`: the wrap-free arithmetic series over
the full-inside crossings (empty when `⌊u₂W⌋ − 1 < ⌈u₁W⌉`). -/
def segPQ (W : ℤ) (s : Seg) : ℚ :=
  let m1 : ℤ := ⌈s.u1 * W⌉
  let m2 : ℤ := ⌊s.u2 * W⌋ - 1
  if m2 < m1 then 0 else
    let K1 : ℚ := (m2 - m1 + 1 : ℤ)
    let α : ℚ := (s.j : ℚ) / (W : ℚ)
    let hf : ℚ := fracQ (s.c - (s.j : ℚ) * (m1 : ℚ) / (W : ℚ))
    let sumh : ℚ := K1 * hf - α * (K1 - 1) * K1 / 2
    let sumpsi : ℚ := K1 * (1/2 : ℚ) - sumh
    s.orient * (α * sumpsi + ((s.j : ℚ)^2 / ((W : ℚ) * ((W : ℚ) + (s.j : ℚ)))) * sumh)

/-- Per-event κ at integer time `W` (the `φ = 0` degeneracy included). -/
def evK (W : ℤ) : Ev → ℚ
  | .death jup jlo x u =>
      let φ := fracQ (u * W)
      if φ = 0 then 0 else
        (W : ℚ) * ((jup : ℚ) - jlo) * max (φ - x) 0
            / (((W : ℚ) + jup) * ((W : ℚ) + jlo))
          - ((jup : ℚ) - jlo) * φ^2 / (2 * W)
  | .birth jup jlo x u =>
      let φ := fracQ (u * W)
      if φ = 0 then 0 else
        (W : ℚ) * ((jlo : ℚ) - jup) * max (x - φ) 0
            / (((W : ℚ) + jup) * ((W : ℚ) + jlo))
          - ((jlo : ℚ) - jup) * (1 - φ)^2 / (2 * W)
  | .swap o jA jB x u =>
      let φ := fracQ (u * W)
      if φ = 0 then 0 else
        let σ := if x ≤ φ then (x * W + (jA : ℚ) * φ) / ((W : ℚ) + jA)
                 else (x * W + (jB : ℚ) * φ) / ((W : ℚ) + jB)
        o * (σ - x - ((jA : ℚ) * φ^2 - (jB : ℚ) * (1 - φ)^2) / (2 * W))
  | .cutR runs uG =>
      let φ := fracQ (uG * W)
      if φ = 0 then 0 else
        (runs.map fun (r : ℚ × ℤ × ℚ × ℤ) =>
          let (xlo, jlo, xup, jup) := r
          let g := (xup - xlo) * φ + ((jup : ℚ) - jlo) * φ^2 / (2 * W)
          let slo := (xlo * W + (jlo : ℚ) * φ) / ((W : ℚ) + jlo)
          let sup := (xup * W + (jup : ℚ) * φ) / ((W : ℚ) + jup)
          max (min sup φ - max slo 0) 0 - g).sum
  | .cutL runs uG =>
      let φ := fracQ (uG * W)
      if φ = 0 then 0 else
        (runs.map fun (r : ℚ × ℤ × ℚ × ℤ) =>
          let (xlo, jlo, xup, jup) := r
          let g := (xup - xlo) * (1 - φ) - ((jup : ℚ) - jlo) * (1 - φ)^2 / (2 * W)
          let slo := (xlo * W + (jlo : ℚ) * φ) / ((W : ℚ) + jlo)
          let sup := (xup * W + (jup : ℚ) * φ) / ((W : ℚ) + jup)
          max (min sup 1 - max slo φ) 0 - g).sum

/-- The total predicted budget at integer time `W`. -/
def budget (segs : List Seg) (evs : List Ev) (W : ℤ) : ℚ :=
  (segs.map (segPQ W)).sum + (evs.map (evK W)).sum

/-- **THM-750, the master identity as a decidable spec.**  For the tail-family body at
integer time `W ≥ Wz`: `W (L − Area) = budget`.  The semantic witness is the exact-ℚ
referee (6/6 Fraction equalities, both shapes, lrc14_closed_budget_thm748_opus_S283.py);
each concrete instance is decidable by evaluation. -/
def ClosedBudgetSpec (segs : List Seg) (evs : List Ev) (W : ℤ) (L Area : ℚ) : Prop :=
  (W : ℚ) * (L - Area) = budget segs evs W

instance (segs : List Seg) (evs : List Ev) (W : ℤ) (L Area : ℚ) :
    Decidable (ClosedBudgetSpec segs evs W L Area) := by
  unfold ClosedBudgetSpec
  infer_instance

/-! ### THM-755: the capped-envelope summation kernel (opus-S289/S291)

The origin cap `|c l| ≤ A` spliced with the spoke envelope `l * |c l| ≤ B` at `m`:
`∑ (c l)² ≤ m A² + B²/m` — the abstract engine behind `disc_v ≤ 4 r |G'|/(π v) + 2|G'|²`
and the (H)-band edge `v* = r/(π |G'|)`. -/

/-- Telescoping tail (strengthened): for `1 ≤ m ≤ N`,
`∑_{l ∈ Ioc m N} 1/l² ≤ 1/m − 1/N`. -/
theorem tail_inv_sq_le_sub (m : ℕ) (hm : 1 ≤ m) :
    ∀ N, m ≤ N → ∑ l ∈ Finset.Ioc m N, (1 : ℝ) / (l : ℝ) ^ 2 ≤ 1 / m - 1 / N := by
  intro N hN
  induction N, hN using Nat.le_induction with
  | base => simp
  | succ n hn ih =>
      rw [Finset.sum_Ioc_succ_top hn]
      have h1n : 1 ≤ n := le_trans hm hn
      have hn0 : (0 : ℝ) < (n : ℝ) := by exact_mod_cast h1n
      have hstep : (1 : ℝ) / ((n : ℝ) + 1) ^ 2 ≤ 1 / n - 1 / ((n : ℝ) + 1) := by
        have hid : (1 : ℝ) / n - 1 / ((n : ℝ) + 1) = 1 / ((n : ℝ) * ((n : ℝ) + 1)) := by
          field_simp
          ring
        rw [hid]
        have hle : (n : ℝ) * ((n : ℝ) + 1) ≤ ((n : ℝ) + 1) ^ 2 := by nlinarith
        exact one_div_le_one_div_of_le (by positivity) hle
      have hgoal := add_le_add ih hstep
      push_cast
      push_cast at hgoal
      linarith

/-- **The capped-envelope kernel (THM-755).**  Origin cap + spoke envelope, spliced at `m`:
`∑_{l ∈ Icc 1 N} (c l)² ≤ m A² + B²/m`. -/
theorem capped_envelope_kernel (c : ℕ → ℝ) (A B : ℝ) (_hA : 0 ≤ A)
    (m N : ℕ) (hm : 1 ≤ m)
    (hcap : ∀ l ∈ Finset.Icc 1 N, |c l| ≤ A)
    (henv : ∀ l ∈ Finset.Icc 1 N, (l : ℝ) * |c l| ≤ B) :
    ∑ l ∈ Finset.Icc 1 N, (c l) ^ 2 ≤ (m : ℝ) * A ^ 2 + B ^ 2 / m := by
  have hm0 : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hBm : 0 ≤ B ^ 2 / m := by positivity
  have hsq : ∀ l ∈ Finset.Icc 1 N, (c l) ^ 2 ≤ A ^ 2 := by
    intro l hl
    have h := abs_le.mp (hcap l hl)
    exact sq_le_sq' h.1 h.2
  by_cases hNm : N ≤ m
  · have hhead : ∑ l ∈ Finset.Icc 1 N, (c l) ^ 2 ≤ (N : ℝ) * A ^ 2 := by
      calc ∑ l ∈ Finset.Icc 1 N, (c l) ^ 2 ≤ ∑ _l ∈ Finset.Icc 1 N, A ^ 2 :=
            Finset.sum_le_sum hsq
        _ = (N : ℝ) * A ^ 2 := by
            rw [Finset.sum_const, Nat.card_Icc]
            simp [nsmul_eq_mul]
    have hNA : (N : ℝ) * A ^ 2 ≤ (m : ℝ) * A ^ 2 := by
      have hc : (N : ℝ) ≤ (m : ℝ) := by exact_mod_cast hNm
      nlinarith [sq_nonneg A]
    linarith
  · have hNm' : m < N := by omega
    have hsplit : Finset.Icc 1 N = Finset.Icc 1 m ∪ Finset.Ioc m N := by
      ext x
      simp only [Finset.mem_Icc, Finset.mem_Ioc, Finset.mem_union]
      omega
    have hdisj : Disjoint (Finset.Icc 1 m) (Finset.Ioc m N) := by
      apply Finset.disjoint_left.mpr
      intro a ha hb
      simp only [Finset.mem_Icc] at ha
      simp only [Finset.mem_Ioc] at hb
      omega
    rw [hsplit, Finset.sum_union hdisj]
    have hhead : ∑ l ∈ Finset.Icc 1 m, (c l) ^ 2 ≤ (m : ℝ) * A ^ 2 := by
      have hsub : ∀ l ∈ Finset.Icc 1 m, (c l) ^ 2 ≤ A ^ 2 := by
        intro l hl
        simp only [Finset.mem_Icc] at hl
        exact hsq l (by simp only [Finset.mem_Icc]; omega)
      calc ∑ l ∈ Finset.Icc 1 m, (c l) ^ 2 ≤ ∑ _l ∈ Finset.Icc 1 m, A ^ 2 :=
            Finset.sum_le_sum hsub
        _ = (m : ℝ) * A ^ 2 := by
            rw [Finset.sum_const, Nat.card_Icc]
            simp [nsmul_eq_mul]
    have htail : ∑ l ∈ Finset.Ioc m N, (c l) ^ 2 ≤ B ^ 2 / m := by
      have hterm : ∀ l ∈ Finset.Ioc m N, (c l) ^ 2 ≤ B ^ 2 * (1 / (l : ℝ) ^ 2) := by
        intro l hl
        simp only [Finset.mem_Ioc] at hl
        have hlN : l ∈ Finset.Icc 1 N := by
          simp only [Finset.mem_Icc]
          omega
        have hl1 : 1 ≤ l := by omega
        have hl0 : (0 : ℝ) < (l : ℝ) := by exact_mod_cast hl1
        have habs : |c l| ≤ B / l := by
          rw [le_div_iff₀ hl0, mul_comm]
          exact henv l hlN
        have h2 := abs_le.mp habs
        have := sq_le_sq' h2.1 h2.2
        calc (c l) ^ 2 ≤ (B / l) ^ 2 := this
          _ = B ^ 2 * (1 / (l : ℝ) ^ 2) := by
              rw [div_pow]
              ring
      have hsum : ∑ l ∈ Finset.Ioc m N, (c l) ^ 2
          ≤ ∑ l ∈ Finset.Ioc m N, B ^ 2 * (1 / (l : ℝ) ^ 2) :=
        Finset.sum_le_sum hterm
      have htel : ∑ l ∈ Finset.Ioc m N, (1 : ℝ) / (l : ℝ) ^ 2 ≤ 1 / m := by
        have := tail_inv_sq_le_sub m hm N (le_of_lt hNm')
        have hN0 : (0 : ℝ) < (N : ℝ) := by
          have : 1 ≤ N := by omega
          exact_mod_cast this
        have : (0 : ℝ) ≤ 1 / (N : ℝ) := by positivity
        linarith
      calc ∑ l ∈ Finset.Ioc m N, (c l) ^ 2
          ≤ ∑ l ∈ Finset.Ioc m N, B ^ 2 * (1 / (l : ℝ) ^ 2) := hsum
        _ = B ^ 2 * ∑ l ∈ Finset.Ioc m N, (1 / (l : ℝ) ^ 2) := by
            rw [Finset.mul_sum]
        _ ≤ B ^ 2 * (1 / m) := by
            have hB2 : (0 : ℝ) ≤ B ^ 2 := sq_nonneg B
            exact mul_le_mul_of_nonneg_left htel hB2
        _ = B ^ 2 / m := by ring
    linarith

end ClosedBudget
end LRC14
end LonelyRunner
