/-
  TournamentH7.LRCArcCount — the ARC-COUNT PIGEONHOLE for good-period existence
  (opus-2026-07-09-S169).

  The finite-`Vmax` glue of the LRC(14) covering capstone (mac-mini-S58/S61, opus-S167/S168).
  A *good period* is an integer `j` with `j/V ∈ Good θ E` (some θ-arc empty of the orbit at the
  rational phase `j/V`).  The dissociated branch closes by EXISTENCE of a good period, obtained
  from the pigeonhole

      (#arcs of `Good θ E ∩ [0,1)`)  <  ρ*·V        ⟹        a good period exists,

  where `ρ* = vol(Good θ E ∩ [0,1))`.  Concretely: if the good set contains `N` intervals of
  total length `> N/V`, one of them is longer than `1/V`, hence catches a grid point `j/V`.

  This file proves that pigeonhole in full (kernel-pure).  The two inputs are supplied elsewhere:
    • `ρ* ≥ D3 ≥ bar > 0` — the density floor (THM-661 / klein finite `D3_c` table), and
    • `#arcs ≤ c·spread` with `c < ρ*` — the bounded-arc-count bound (mac-mini-S58; the trivial
      `O(k² spread)` bound is `~10³×` too loose, so the small-constant a-priori bound remains the
      one open analytic item; here it enters as the hypothesis `N/V < total-length`).

  The heart — `exists_gridpoint_Ico` (a long interval catches a grid point) and `exists_long_Ico`
  (of `N` intervals summing to `> N/V`, one exceeds `1/V`) — is elementary and unconditional.
-/
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace TailDiameter

open scoped ENNReal
open Set

/-! ### 1. A long interval catches a grid point

If `Ico a b` is longer than `1/V`, then some grid point `j/V` (`j : ℤ`) lies in it.
Take `j = ⌈a·V⌉`: then `a·V ≤ j < a·V + 1 < b·V`, so `a ≤ j/V < b`. -/

theorem exists_gridpoint_Ico {V : ℕ} (hV : 0 < V) {a b : ℝ}
    (hlen : 1 / (V : ℝ) < b - a) : ∃ j : ℤ, (j : ℝ) / V ∈ Ico a b := by
  have hVR : (0 : ℝ) < V := by exact_mod_cast hV
  refine ⟨⌈a * V⌉, ?_, ?_⟩
  · -- `a ≤ ⌈a·V⌉ / V`
    rw [le_div_iff₀ hVR]
    exact Int.le_ceil (a * V)
  · -- `⌈a·V⌉ / V < b`
    rw [div_lt_iff₀ hVR]
    -- from `1/V < b - a` : `V·a + 1 < V·b`
    have h1 : (1 : ℝ) < V * (b - a) := by
      have := (div_lt_iff₀ hVR).1 hlen        -- `1 < (b - a) * V`
      nlinarith [this]
    have hceil : (⌈a * V⌉ : ℝ) < a * V + 1 := Int.ceil_lt_add_one (a * V)
    nlinarith [hceil, h1]

/-! ### 2. Pigeonhole on interval lengths

Of `N` intervals whose lengths sum to more than `N/V`, at least one exceeds `1/V`. -/

theorem exists_long_Ico {V : ℕ} (hV : 0 < V) {N : ℕ} (a b : Fin N → ℝ)
    (hsum : (N : ℝ) / V < ∑ i, (b i - a i)) : ∃ i, 1 / (V : ℝ) < b i - a i := by
  by_contra h
  push_neg at h                                  -- `∀ i, b i - a i ≤ 1/V`
  have hle : (∑ i, (b i - a i)) ≤ ∑ _i : Fin N, (1 / (V : ℝ)) :=
    Finset.sum_le_sum (fun i _ => h i)
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hle
  -- `∑ ≤ N · (1/V) = N/V`, contradicting `N/V < ∑`
  rw [mul_one_div] at hle
  linarith [hsum, hle]

/-! ### 3. The arc-count pigeonhole: a good period exists

If the good set `Good θ E` contains `N` intervals of total length `> N/V`, then some grid
point `j/V` is good — a *good period*.  (Take `Good θ E ∩ [0,1)` to be the union of its `N`
maximal arcs; then `∑ lengths = ρ*` and the hypothesis is exactly `N < ρ*·V`.) -/

theorem good_period_of_arccount {θ : ℝ} {E : Finset ℤ} {V : ℕ} (hV : 0 < V) {N : ℕ}
    (a b : Fin N → ℝ)
    (hcover : (⋃ i, Ico (a i) (b i)) ⊆ Good θ E)
    (hmeas : (N : ℝ) / V < ∑ i, (b i - a i)) :
    ∃ j : ℤ, ((j : ℝ) / V) ∈ Good θ E := by
  obtain ⟨i, hi⟩ := exists_long_Ico hV a b hmeas
  obtain ⟨j, hj⟩ := exists_gridpoint_Ico hV hi
  exact ⟨j, hcover (mem_iUnion.2 ⟨i, hj⟩)⟩

/-- Contrapositive reading: **no good period ⟹ the good arcs are few** (`ρ*·V ≤ #arcs`).
This is the exact obstruction — realized only by the tight complete-residue AP at `V` prime
(opus-S164), the LRC extremal.  Any speed set beating `#arcs < ρ*·V` has a good period. -/
theorem arccount_le_of_no_good_period {θ : ℝ} {E : Finset ℤ} {V : ℕ} (hV : 0 < V) {N : ℕ}
    (a b : Fin N → ℝ)
    (hcover : (⋃ i, Ico (a i) (b i)) ⊆ Good θ E)
    (hno : ∀ j : ℤ, ((j : ℝ) / V) ∉ Good θ E) :
    (∑ i, (b i - a i)) ≤ (N : ℝ) / V := by
  by_contra h
  push_neg at h
  obtain ⟨j, hj⟩ := good_period_of_arccount hV a b hcover h
  exact hno j hj

/-! ### 4. The SMOOTH averaging route (opus-S170)

The sharp good-set indicator has Fourier decay `~1/m` (L¹-divergent resonant sum ⟹ existence needs
signed CANCELLATION, Mertens-hard — and the arc-count bound is vacuous on the extremal near-AP family,
klein MISTAKE-127).  The *smooth* surrogate `maxgap(x)` is continuous piecewise-linear, Fourier decay
`~1/m^α` with `α>1` (opus-S170: `α≈1.5–2.0`), so its resonant sum converges ABSOLUTELY — no
cancellation.  A good period is `j` with `maxgap(j/V) > 1/7`; by `max ≥ mean` it suffices that the
ruler-grid MEAN of `maxgap` exceeds `1/7`, and the grid mean equals the dilation-invariant continuous
mean `E_x[maxgap]` up to the (small, absolutely-bounded) resonant discrepancy `D`.  Since
`E_x[maxgap] > 1/7` with a uniform margin (opus-S170: `≥1.047×`, and `≈1.48×` for the extremal AP
`{1..13}`), `E_j[maxgap] > 1/7` a-priori once `D <` margin. -/

/-- **Smooth averaging existence.**  `g` is the (max-gap) surrogate on a nonempty ruler grid `J`;
`meanx` its dilation-invariant continuous mean; the grid mean sits within `D` of `meanx`
(`|mean_J g − meanx| ≤ D`, the resonant discrepancy — absolutely bounded since `α>1`); and the
a-priori margin `thr + D < meanx` holds.  Then some grid point beats the threshold: a good period. -/
theorem exists_good_of_smooth_mean {ι : Type*} (J : Finset ι) (hJ : J.Nonempty) (g : ι → ℝ)
    (thr meanx D : ℝ)
    (hdisc : |(∑ j ∈ J, g j) / (J.card : ℝ) - meanx| ≤ D)
    (hmargin : thr + D < meanx) :
    ∃ j ∈ J, thr < g j := by
  have hcard : (0 : ℝ) < (J.card : ℝ) := by exact_mod_cast hJ.card_pos
  rw [abs_le] at hdisc
  have hmean : thr < (∑ j ∈ J, g j) / (J.card : ℝ) := by linarith [hdisc.1]
  by_contra h
  push_neg at h                                     -- `∀ j ∈ J, g j ≤ thr`
  have hsum : (∑ j ∈ J, g j) ≤ ∑ _j ∈ J, thr := Finset.sum_le_sum (fun j hj => h j hj)
  rw [Finset.sum_const, nsmul_eq_mul] at hsum        -- `∑ ≤ card · thr`
  rw [lt_div_iff₀ hcard] at hmean                    -- `thr · card < ∑`
  nlinarith [hmean, hsum]

/-! ### 5. The residual is ABSOLUTELY bounded — no cancellation (opus-S171)

kps-S96's `E_grid[W]` route closes the good period once the resonance residual `R = Σ_{V|m} 𝒲̂(m)`
satisfies `|R| < (6/7)^k`.  The question was whether `R` is small only by SIGNED cancellation
(Mertens-hard) or by an ABSOLUTE bound.  opus-S171 (verified, adversarial): `R_abs := Σ_{n≥1}
2|𝒲̂(nV)| < (6/7)^k` for ALL clusters — dissociated (`≤0.40·main`) AND the 7-structured hard case
(`≤0.41·main`, MISTAKE-128) — because `W = Σ(gap−1/7)_+` is CONTINUOUS piecewise-linear so `𝒲̂(m)`
decays (`~1/m^α`, `α>1`; opus-S170), making the resonant sum converge ABSOLUTELY.  So `|R| ≤ R_abs <
main` needs no cancellation.  (For 7-structured `E` the arc-Fourier `b(m)=(1−e(m/7))/(2πim)=0` at
`7|m` additionally suppresses the balanced resonances — opus-S167.)  This lemma turns the signed
hypothesis `|R| < main` into the absolute one `Σ|t| < main`.

**opus-S172 CORRECTION.**  The S171 "no cancellation" reading was an OVER-CLAIM (parallel to kps-S97→
S98).  The absolutely-summable object is the *per-frequency* `Σ_m |𝒲̂₁(m)|` (`< main`), NOT the
*per-covering* `Σ_n |𝒲̂(n)|` (which DIVERGES — the LEM-011 shells grow).  And even the per-frequency
`R_abs < main` is NOT provable by the crude BV bound `|𝒲̂₁(m)| ≤ TV(W')/(2πm)²`, because
`TV(W') ~ spread²` (measured 2.03) makes that bound `≈ 8·main` — the true `R_abs` is `40–200×` below by
CANCELLATION (the Mertens wall).  Root cause: the LRC covering over-covers (`k/7 = 13/7 > 1`).  So
`abs_residual_lt` remains a true *reduction*, but its hypothesis `Σ|t| < main` is a cancellation fact,
not an a-priori absolute bound.  The a-priori `|R|<main` is NOT proof-critical (mac-mini-S64's `j=1`
wraparound + LEM-013 exhaustion close the good period). -/
theorem abs_residual_lt {ι : Type*} (S : Finset ι) (t : ι → ℝ) (main : ℝ)
    (habs : ∑ n ∈ S, |t n| < main) : |∑ n ∈ S, t n| < main :=
  lt_of_le_of_lt (Finset.abs_sum_le_sum_abs t S) habs

/-! ### 6. The resonant TAIL is a-priori (opus-S172)

What IS rigorous in the resonant sum is the high-frequency TAIL.  If the resonant terms decay as
`|a n| ≤ C/n²` (the `p=2` BV rate of the smooth `W`), the tail `Σ_{M<n≤N}|a n| ≤ C/M` — a-priori,
telescoping `1/n² ≤ 1/((n−1)n) = 1/(n−1)−1/n`.  So the resonant residual splits as `|R| ≤ (finite head
n ≤ M) + C/M`; the tail is the easy part and only the head carries the Mertens cancellation.  This
delimits exactly where the a-priori boundary lies. -/
theorem resonant_tail_le {a : ℕ → ℝ} {C : ℝ} (hC : 0 ≤ C) {M N : ℕ} (hM : 1 ≤ M)
    (hb : ∀ n, M + 1 ≤ n → |a n| ≤ C / (n : ℝ) ^ 2) :
    ∑ n ∈ Finset.Ico (M + 1) (N + 1), |a n| ≤ C / (M : ℝ) := by
  have hMpos : (0 : ℝ) < (M : ℝ) := by exact_mod_cast hM
  set F : ℕ → ℝ := fun k => C / (k : ℝ) with hF
  -- termwise: |a n| ≤ F(n-1) - F n on the tail range
  have step : ∀ n ∈ Finset.Ico (M + 1) (N + 1), |a n| ≤ F (n - 1) - F n := by
    intro n hn
    rw [Finset.mem_Ico] at hn
    have h1 : M + 1 ≤ n := hn.1
    have hn1 : 1 ≤ n := by omega
    have hn2 : (2 : ℝ) ≤ (n : ℝ) := by exact_mod_cast (show 2 ≤ n by omega)
    have hp : (0 : ℝ) < (n : ℝ) := by linarith
    have hp1 : (0 : ℝ) < (n : ℝ) - 1 := by linarith
    have hcast : ((n - 1 : ℕ) : ℝ) = (n : ℝ) - 1 := by
      rw [Nat.cast_sub hn1]; simp
    have hle : C / (n : ℝ) ^ 2 ≤ C / (((n : ℝ) - 1) * (n : ℝ)) := by
      gcongr
      nlinarith
    have heq : C / (((n : ℝ) - 1) * (n : ℝ)) = F (n - 1) - F n := by
      simp only [hF, hcast]; field_simp; ring
    exact le_trans (hb n h1) (le_trans hle (le_of_eq heq))
  refine le_trans (Finset.sum_le_sum step) ?_
  -- telescoping ∑ (F(n-1) - F n) = F M - F N ≤ F M = C/M
  have tele : ∀ P, M ≤ P →
      ∑ n ∈ Finset.Ico (M + 1) (P + 1), (F (n - 1) - F n) = F M - F P := by
    intro P hP
    induction P, hP using Nat.le_induction with
    | base => simp
    | succ Q hQ ih =>
      rw [Finset.sum_Ico_succ_top (by omega : M + 1 ≤ Q + 1), ih]
      simp only [Nat.add_sub_cancel]
      ring
  by_cases hMN : M ≤ N
  · rw [tele N hMN]
    have hFN : (0 : ℝ) ≤ F N := by
      simp only [hF]; positivity
    have hFM : F M = C / (M : ℝ) := rfl
    linarith
  · -- N < M ⟹ Ico (M+1) (N+1) empty
    rw [Finset.Ico_eq_empty (by omega), Finset.sum_empty]
    positivity

end TailDiameter
end LonelyRunner
