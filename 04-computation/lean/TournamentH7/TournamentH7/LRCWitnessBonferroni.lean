/-
  TournamentH7.LRCWitnessBonferroni -- the ELEMENTARY Bonferroni reduction of the
  LRC(14) global 1/7-witness floor (kind-pasteur-2026-06-22-S30).

  The remaining analytic crux of the witness route is a uniform positive floor for
  the global-witness density
      rho*_glob(P,E) = witnessG2(P,E)
                     = meas{ x in G_P : circular maxgap{frac(e x): e in E} > 1/7 }.
  THM-527 Part G called this "the genuine remaining crux" and framed it as a
  COMPACTNESS problem (continuity + closure on a bounded-spread shape space).

  THIS module records the structural advance that DISSOLVES the compactness
  obstruction.  Writing `GOOD(E) = {x : maxgap{frac(e x)} > 1/7}`, with
  `nu(E) = meas(GOOD(E))` (NO small part) and `measGP(P) = meas(G_P)`, the
  elementary union bound (Bonferroni) gives

      witnessG2 = meas(GOOD cap G_P) >= nu(E) + meas(G_P) - 1.            (BONFERRONI)

  Hence the floor reduces to two DECOUPLED universal lower bounds plus finite
  arithmetic:

      LEMMA A (pure three-distance):  nu(E)      >= nu_consec(k);
      LEMMA B (proved finite rational): meas(G_P) >= cap_k = min meas(G_P);
      ARITHMETIC: nu_consec(k) + cap_k - 1 > 0  (and >= m_P) for all k = 8..13.

  `nu` is SCALE-INVARIANT (`nu(cE)=nu(E)`), so Lemma A is a pure primitive-shape
  statement; it is VERIFIED (consecutive is the strict minimizer; wide shapes
  decorrelate to nu->1, the SAFE direction) and its rigorous closure is the
  [finite-core + decorrelation-tail] architecture (same as the proven single-far
  / Leg-C legs).  Lemma B is an EXACT finite computation (`lrc_capGP_exact_kps.py`)
  with min_{|P|=13-k} meas(G_P) = capRat(k) exactly (a duality with the p0 cap).

  WHAT IS PROVED HERE (sorry-free, axiom-free beyond Lean foundations +
  native_decide):
    * the exact-rational floor table `nu_consec(k) + cap_k - 1 > 0` and `>= m_P`
      (`bonferroni_floor_pos`, `bonferroni_floor_ge_mP`);
    * the STRUCTURAL reduction `witness_floor_from_bonferroni_nodes`: given the
      Bonferroni inequality (a measure union bound), Lemma A and Lemma B as
      hypotheses, the witness floor `m_P <= witnessG2 s` follows for every shape
      with cluster size 8..13.  This carries NO `sorry`: it isolates exactly the
      three obligations, two of which are essentially settled.
    * the top-level assembly wrappers
      `lrc14_from_bonferroni_split_nodes` and
      `lrc14_from_p0_wide_bound_split_nodes`, which feed those floor reductions
      into the skeleton's LRC(14) witness-route theorem.  They do not assert the
      remaining analytic nodes; they make the exact Lean boundary explicit.

  This replaces the single opaque `hfloor` obligation by Bonferroni (proved
  measure fact) + Lemma B (proved finite rational) + decidable arithmetic +
  Lemma A (the one isolated three-distance lemma).
-/

import Mathlib.Tactic
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14
namespace Bonferroni

open scoped BigOperators

/-! ## The exact rational floor table

`nuConsec k` = `meas{x : the Steinhaus orbit {0,x,...,(k-1)x} has maxgap > 1/7}`,
an exact rational (three-distance; good `x` live only near `a/b`, `b <= 6`).
`capRat k` (reused from the skeleton) = `min_{|P|=13-k} meas(G_P)` (the duality).
For `k <= 7`, `nu = 1` by pigeonhole (k points => maxgap >= 1/k >= 1/7). -/

/-- `nu_consec(k)` for `k = 8..13` (exact rationals, `lrc_nu_floor_and_tail_kps.py`);
`= 1` for `k <= 7` (pigeonhole). -/
def nuConsec : ℕ → ℚ
  | 8  => 691 / 735
  | 9  => 247 / 294
  | 10 => 38 / 49
  | 11 => 1381 / 2205
  | 12 => 13823 / 24255
  | 13 => 477 / 1078
  | _  => 1

/-- **The elementary Bonferroni floor is strictly positive** for every binding
cluster size `k = 8..13`: `nu_consec(k) + cap_k - 1 > 0`.  Finite exact-rational
fact (worst value `1891/5880 ≈ 0.3216` at `k = 8`). -/
theorem bonferroni_floor_pos :
    ∀ k ∈ [8, 9, 10, 11, 12, 13], 0 < nuConsec k + capRat k - 1 := by
  native_decide

/-- **The Bonferroni floor dominates the admissible witness floor** `m_P = witnessMP
= 14249/252252` for every `k = 8..13`.  Hence this route proves the skeleton's
`hfloor : witnessMP <= witnessG2` obligation a fortiori. -/
theorem bonferroni_floor_ge_mP :
    ∀ k ∈ [8, 9, 10, 11, 12, 13], witnessMP ≤ nuConsec k + capRat k - 1 := by
  native_decide

/-- The same floor lower bound, stated per `k` with order hypotheses (so it can be
specialized to `k = clusterSize s`), as a `ℚ` inequality. -/
theorem floor_ge_mP_of_mem (k : ℕ) (h8 : 8 ≤ k) (h13 : k ≤ 13) :
    witnessMP ≤ nuConsec k + capRat k - 1 := by
  interval_cases k <;> native_decide

/-! ## The structural reduction (sorry-free, axiom-free)

The witness floor for any shape with cluster size `8..13` follows from:
* `hbonf` -- the Bonferroni / union bound `nu + measGP - 1 <= witnessG2` (a
  Lebesgue-measure fact: `meas(A ∩ B) >= meas A + meas B - 1` for `A,B ⊆ [0,1]`);
* `hA` -- **Lemma A** `nu_consec(k) <= nu(s)` (the scale-invariant three-distance
  floor; VERIFIED, rigorous closure = finite core + decorrelation tail);
* `hB` -- **Lemma B** `cap_k <= measGP(s)` (PROVED finite rational, `min meas(G_P)`).

We carry `nuShape, measGP : Shape → ℝ` as parameters (the Lebesgue measures of
`GOOD(E)` and `G_P`); the theorem is then pure real-arithmetic glue. -/

/-- **The Bonferroni witness-floor reduction.**  Given the Bonferroni measure
inequality, Lemma A, and Lemma B, the admissible witness floor `witnessMP` lower-
bounds `witnessG2 s` for every shape whose cluster size lies in `8..13`. -/
theorem witness_floor_from_bonferroni_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (s : Shape) (h8 : 8 ≤ clusterSize s) (h13 : clusterSize s ≤ 13) :
    (witnessMP : ℝ) ≤ witnessG2 s := by
  set k := clusterSize s with hk
  -- (1) the finite arithmetic floor, cast to ℝ
  have hrat : witnessMP ≤ nuConsec k + capRat k - 1 := floor_ge_mP_of_mem k h8 h13
  have hcast : (witnessMP : ℝ) ≤ (nuConsec k : ℝ) + (capRat k : ℝ) - 1 := by
    have := (Rat.cast_le (K := ℝ)).mpr hrat
    push_cast at this
    linarith
  -- (2) Lemma A + Lemma B raise nu_consec + cap to the genuine measures
  have hAk : (nuConsec k : ℝ) ≤ nuShape s := hA s h8 h13
  have hBk : (capRat k : ℝ) ≤ measGP s := hB s h8 h13
  -- (3) Bonferroni
  have hbk : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  -- chain
  linarith

/-- The large-cluster (`8..13`) floor node in the exact shape expected by the
LRC14 skeleton's case-split theorem. -/
theorem large_witness_floor_from_bonferroni_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v _hv h8 h13
  exact witness_floor_from_bonferroni_nodes nuShape measGP hbonf hA hB
    (shapeOf v) h8 h13

/-! ## ★ The UNIFICATION with the p0 wide bound (the rigorous closure)

The witness floor is IMPLIED by the team's already-closed p0 wide bound, with NO
new analytic lemma.  Writing the 1/7-dense measure `D = 1 - nu`, the elementary
set inclusion `{1/7-dense} ⊆ {all 6 inner sectors hit}` gives `D(E) ≤ p0(E)`
(`p0 = measS7`, the repo cover atom).  Hence by Bonferroni

    witnessG2 ≥ nu + measGP - 1 = measGP - D ≥ measGP - p0,

and with `p0 ≤ cap_k - δ_k` (the wide bound, `δ_k > 0` the team's margin) and
`cap_k ≤ measGP` (Lemma B + the duality `cap_k = min meas(G_P)`), we get
`witnessG2 ≥ δ_k > 0`.  This carries `p0Shape` as a parameter; `hDp0` is the
elementary inclusion, `hp0cap` the team's wide bound, `hmeasGP` Lemma B. -/

/-- **The witness floor via the p0 wide bound (UNIFICATION).**  Given Bonferroni,
the elementary `D ≤ p0` inclusion, the wide bound `p0 ≤ cap - δ`, and the cap floor
`cap ≤ measGP`, the witness density is at least the wide-bound margin `δ > 0`. -/
theorem witness_floor_from_p0_wide_bound
    (nuShape measGP p0Shape : Shape → ℝ) (cap delta : ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)          -- D = 1 - nu ≤ p0
    (hp0cap : ∀ s, p0Shape s ≤ cap - delta)             -- team's wide bound, margin δ
    (hmeasGP : ∀ s, cap ≤ measGP s)                     -- Lemma B + duality
    (s : Shape) :
    delta ≤ witnessG2 s := by
  have h1 : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have h2 : (1 - nuShape s) ≤ p0Shape s := hDp0 s
  have h3 : p0Shape s ≤ cap - delta := hp0cap s
  have h4 : cap ≤ measGP s := hmeasGP s
  -- witnessG2 ≥ nu+measGP-1 = measGP-(1-nu) ≥ measGP-p0 ≥ cap-(cap-δ) = δ
  linarith

/-- The same conclusion as strict positivity of the witness density, when the
wide-bound margin `delta` is positive. -/
theorem witnessG2_pos_from_p0_wide_bound
    (nuShape measGP p0Shape : Shape → ℝ) (cap delta : ℝ) (hδ : 0 < delta)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap - delta)
    (hmeasGP : ∀ s, cap ≤ measGP s)
    (s : Shape) :
    0 < witnessG2 s :=
  lt_of_lt_of_le hδ
    (witness_floor_from_p0_wide_bound nuShape measGP p0Shape cap delta
      hbonf hDp0 hp0cap hmeasGP s)

/-- Shape-indexed version of `witness_floor_from_p0_wide_bound`, with the margin
allowed to depend on the cluster shape.  This is the exact form needed for the
`k`-dependent wide-bound margins in the LRC14 proof DAG. -/
theorem witness_floor_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hδm : ∀ s, (witnessMP : ℝ) ≤ delta s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (s : Shape) :
    (witnessMP : ℝ) ≤ witnessG2 s := by
  have h1 : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have h2 : (1 - nuShape s) ≤ p0Shape s := hDp0 s
  have h3 : p0Shape s ≤ cap s - delta s := hp0cap s
  have h4 : cap s ≤ measGP s := hmeasGP s
  have h5 : (witnessMP : ℝ) ≤ delta s := hδm s
  linarith

/-- Large-cluster p0-unification floor in the shape expected by the skeleton's
case-split theorem.  The hypotheses are restricted to `8..13`, matching the
current sector/witness proof split. -/
theorem large_witness_floor_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hδm : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (witnessMP : ℝ) ≤ delta s)
    (hbonf : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      cap s ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v _hv h8 h13
  have h1 : nuShape (shapeOf v) + measGP (shapeOf v) - 1 ≤
      witnessG2 (shapeOf v) := hbonf (shapeOf v) h8 h13
  have h2 : (1 - nuShape (shapeOf v)) ≤ p0Shape (shapeOf v) :=
    hDp0 (shapeOf v) h8 h13
  have h3 : p0Shape (shapeOf v) ≤ cap (shapeOf v) - delta (shapeOf v) :=
    hp0cap (shapeOf v) h8 h13
  have h4 : cap (shapeOf v) ≤ measGP (shapeOf v) :=
    hmeasGP (shapeOf v) h8 h13
  have h5 : (witnessMP : ℝ) ≤ delta (shapeOf v) := hδm (shapeOf v) h8 h13
  linarith

/-- Shape-indexed p0-unification margin, with no comparison to the placeholder
`witnessMP`.  This is the positive-floor form needed by Part A: Bonferroni,
`D <= p0`, the wide p0 margin, and the cap floor give
`delta s <= witnessG2 s`. -/
theorem witness_margin_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (s : Shape) :
    delta s ≤ witnessG2 s := by
  have h1 : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have h2 : (1 - nuShape s) ≤ p0Shape s := hDp0 s
  have h3 : p0Shape s ≤ cap s - delta s := hp0cap s
  have h4 : cap s ≤ measGP s := hmeasGP s
  linarith

/-- Large-cluster p0-unification positivity in the shape expected by the current
split proof DAG.  Unlike `large_witness_floor_from_p0_wide_bound_shapes`, this
requires only a positive p0 margin, not `witnessMP <= delta`. -/
theorem large_witness_pos_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hδpos : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 → 0 < delta s)
    (hbonf : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      cap s ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      0 < witnessG2 (shapeOf v) := by
  intro v _hv h8 h13
  have hmargin : delta (shapeOf v) ≤ witnessG2 (shapeOf v) := by
    have h1 : nuShape (shapeOf v) + measGP (shapeOf v) - 1 ≤
        witnessG2 (shapeOf v) := hbonf (shapeOf v) h8 h13
    have h2 : (1 - nuShape (shapeOf v)) ≤ p0Shape (shapeOf v) :=
      hDp0 (shapeOf v) h8 h13
    have h3 : p0Shape (shapeOf v) ≤ cap (shapeOf v) - delta (shapeOf v) :=
      hp0cap (shapeOf v) h8 h13
    have h4 : cap (shapeOf v) ≤ measGP (shapeOf v) :=
      hmeasGP (shapeOf v) h8 h13
    linarith
  exact lt_of_lt_of_le (hδpos (shapeOf v) h8 h13) hmargin

/-- Shape-wise positivity version of the p0-wide-bound route.  This is weaker
than `witness_floor_from_p0_wide_bound_shapes`, but it is the natural endpoint
when the wide-bound margin is known only as a positive shape-wise quantity. -/
theorem witnessG2_pos_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape capShape deltaShape : Shape → ℝ)
    (hδ : ∀ s, 0 < deltaShape s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ capShape s - deltaShape s)
    (hmeasGP : ∀ s, capShape s ≤ measGP s)
    (s : Shape) :
    0 < witnessG2 s := by
  have h1 : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have h2 : (1 - nuShape s) ≤ p0Shape s := hDp0 s
  have h3 : p0Shape s ≤ capShape s - deltaShape s := hp0cap s
  have h4 : capShape s ≤ measGP s := hmeasGP s
  have hδs : 0 < deltaShape s := hδ s
  linarith

/-- Top-level conditional LRC14 assembly through the positive p0/witness
unification route.  This formulation records the weakest endpoint needed by
Part A: `G2 > 0`, rather than the stronger uniform `m_P` floor. -/
theorem lrc14_from_p0_wide_bound_given_partA
    (nuShape measGP p0Shape capShape deltaShape : Shape → ℝ)
    (hδ : ∀ s, 0 < deltaShape s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ capShape s - deltaShape s)
    (hmeasGP : ∀ s, capShape s ≤ measGP s)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  intro v hv
  have hpos : 0 < witnessG2 (shapeOf v) :=
    witnessG2_pos_from_p0_wide_bound_shapes
      nuShape measGP p0Shape capShape deltaShape hδ hbonf hDp0 hp0cap hmeasGP
      (shapeOf v)
  exact lonely_of_Mreach_ge v hv (hpartA v hpos)

/-! ## The `k <= 7` pigeonhole leg (no Lemma A needed)

For `k <= 7` the cluster has `<= 7` phases, so `maxgap >= 1/k >= 1/7` for EVERY
`x`; thus `GOOD(E) = [0,1)` (full), `nu = 1`, and `witnessG2 = measGP >= cap_k > 0`.
We record the arithmetic that `cap_k > 0` and that `nu = 1` here is exact. -/

/-- For `k <= 7`, `nuConsec k = 1`. -/
theorem nuConsec_eq_one_of_le_seven (k : ℕ) (hk : k ≤ 7) : nuConsec k = 1 := by
  interval_cases k <;> rfl

/-- In the `k <= 7` pigeonhole leg, with `nu = 1` the Bonferroni bound is exactly
`witnessG2 >= measGP`, so any positive `cap_k` floor on `measGP` gives a positive
witness floor.  (Stated as the clean specialization.) -/
theorem witness_floor_pigeonhole_leg
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hB : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (s : Shape) (h7 : clusterSize s ≤ 7) :
    (capRat (clusterSize s) : ℝ) ≤ witnessG2 s := by
  have hbk : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have hnu : nuShape s = 1 := hnu1 s h7
  have hBk : (capRat (clusterSize s) : ℝ) ≤ measGP s := hB s h7
  -- nu = 1 ⟹ measGP ≤ witnessG2 ⟹ cap ≤ witnessG2
  have : measGP s ≤ witnessG2 s := by linarith
  linarith

/-- The admissible witness floor `m_P` is below the trivial small-cluster cap
`capRat k = 1` for `k <= 7`. -/
theorem witnessMP_le_capRat_of_le_seven (k : ℕ) (hk : k ≤ 7) :
    (witnessMP : ℝ) ≤ (capRat k : ℝ) := by
  interval_cases k <;> norm_num [witnessMP, capRat]

/-- The pigeonhole leg promoted to the skeleton's `m_P <= witnessG2` floor. -/
theorem witness_floor_pigeonhole_leg_ge_mP
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hB : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (s : Shape) (h7 : clusterSize s ≤ 7) :
    (witnessMP : ℝ) ≤ witnessG2 s :=
  le_trans (witnessMP_le_capRat_of_le_seven (clusterSize s) h7)
    (witness_floor_pigeonhole_leg nuShape measGP hbonf hnu1 hB s h7)

/-- Small-cluster (`k <= 7`) floor node in the exact shape expected by the
LRC14 skeleton's case-split theorem. -/
theorem small_witness_floor_from_pigeonhole_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hB : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v _hv h7
  exact witness_floor_pigeonhole_leg_ge_mP nuShape measGP hbonf hnu1 hB
    (shapeOf v) h7

/-- Compatibility name for the small-cluster Bonferroni feeder. -/
theorem witness_small_floor_from_bonferroni_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hB : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) :=
  small_witness_floor_from_pigeonhole_nodes nuShape measGP hbonf hnu1 hB

/-- Compatibility name for the large-cluster Bonferroni feeder. -/
theorem witness_large_floor_from_bonferroni_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) :=
  large_witness_floor_from_bonferroni_nodes nuShape measGP hbonf hA hB

/-- Full all-cluster witness-floor bridge from the small and large Bonferroni
nodes, in the monolithic floor shape expected by `lrc14_from_witness_floor`. -/
theorem witness_floor_from_bonferroni_nodes_all_clusters
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hBsmall : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hBlarge : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) :=
  witness_floor_from_cluster_cases
    (witness_small_floor_from_bonferroni_nodes nuShape measGP hbonf hnu1 hBsmall)
    (witness_large_floor_from_bonferroni_nodes nuShape measGP hbonf hA hBlarge)
    hsize

/-! ## LRC14 assembly from the Bonferroni / p0 floor nodes -/

/-- **Top-level LRC14 assembly from the Bonferroni split floor.**  The remaining
inputs are the actual measure facts: Bonferroni, the `k<=7` full-good
pigeonhole statement, Lemma A, Lemma B, the cluster-size bound, and the direct
`G2 > 0 -> Mreach >= 1/14` Part-A bridge. -/
theorem lrc14_from_bonferroni_split_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hBsmall : ∀ s, clusterSize s ≤ 7 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hBlarge : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor_cases_given_nodes
    (small_witness_floor_from_pigeonhole_nodes nuShape measGP hbonf hnu1 hBsmall)
    (large_witness_floor_from_bonferroni_nodes nuShape measGP hbonf hA hBlarge)
    hsize hpartA lonely_of_Mreach_ge

/-- Compatibility top-level assembly through the all-cluster Bonferroni floor
node.  This is logically the same route as `lrc14_from_bonferroni_split_nodes`,
but exposes the monolithic floor theorem for downstream audits. -/
theorem lrc14_from_bonferroni_nodes_given_partA
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hnu1 : ∀ s, clusterSize s ≤ 7 → nuShape s = 1)
    (hA : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (nuConsec (clusterSize s) : ℝ) ≤ nuShape s)
    (hBsmall : ∀ s, clusterSize s ≤ 7 → (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hBlarge : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor
    (witness_floor_from_bonferroni_nodes_all_clusters
      nuShape measGP hbonf hnu1 hA hBsmall hBlarge hsize)
    hpartA

/-- **Top-level LRC14 assembly from an all-shapes p0 wide-bound margin.**  This
is the cleanest formal reading of HYP-2832: a positive p0 margin at least `m_P`,
plus `D <= p0`, `p0 <= cap-delta`, and `cap <= measGP`, supplies the witness
floor needed by the skeleton. -/
theorem lrc14_from_p0_wide_bound_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hδm : ∀ s, (witnessMP : ℝ) ≤ delta s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor_given_nodes
    (fun v _hv =>
      witness_floor_from_p0_wide_bound_shapes nuShape measGP p0Shape cap delta
        hδm hbonf hDp0 hp0cap hmeasGP (shapeOf v))
    hpartA lonely_of_Mreach_ge

/-- **Top-level LRC14 assembly from the current split p0 route.**  Small clusters
are supplied separately (typically by pigeonhole); large clusters are discharged
by the p0 wide-bound margin. -/
theorem lrc14_from_p0_wide_bound_split_nodes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hδm : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (witnessMP : ℝ) ≤ delta s)
    (hbonf : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      cap s ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor_cases_given_nodes hsmall
    (large_witness_floor_from_p0_wide_bound_shapes
      nuShape measGP p0Shape cap delta hδm hbonf hDp0 hp0cap hmeasGP)
    hsize hpartA lonely_of_Mreach_ge

/-- **Top-level LRC14 assembly from the corrected positive-margin p0 route.**
The p0 route only needs to make `witnessG2` positive before invoking Part A; it
does not need the p0 margin to dominate the conservative placeholder
`witnessMP`.  Small clusters may still be supplied by the existing `m_P` floor,
while large clusters require only `0 < delta`. -/
theorem lrc14_from_p0_positive_wide_bound_split_nodes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hδpos : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 → 0 < delta s)
    (hbonf : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      cap s ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  intro v hv
  have hpos : 0 < witnessG2 (shapeOf v) := by
    by_cases h7 : clusterSize (shapeOf v) ≤ 7
    · exact witness_floor_positive (shapeOf v) (hsmall v hv h7)
    · have h8 : 8 ≤ clusterSize (shapeOf v) := by omega
      exact large_witness_pos_from_p0_wide_bound_shapes
        nuShape measGP p0Shape cap delta hδpos hbonf hDp0 hp0cap hmeasGP
        v hv h8 (hsize v hv)
  exact lonely_of_Mreach_ge v hv (hpartA v hpos)

/-! ## Axiom audit -/

#print axioms bonferroni_floor_pos
#print axioms bonferroni_floor_ge_mP
#print axioms floor_ge_mP_of_mem
#print axioms witness_floor_from_bonferroni_nodes
#print axioms large_witness_floor_from_bonferroni_nodes
#print axioms witness_floor_from_p0_wide_bound
#print axioms witnessG2_pos_from_p0_wide_bound
#print axioms witness_floor_from_p0_wide_bound_shapes
#print axioms large_witness_floor_from_p0_wide_bound_shapes
#print axioms witness_margin_from_p0_wide_bound_shapes
#print axioms large_witness_pos_from_p0_wide_bound_shapes
#print axioms witnessG2_pos_from_p0_wide_bound_shapes
#print axioms lrc14_from_p0_wide_bound_given_partA
#print axioms nuConsec_eq_one_of_le_seven
#print axioms witness_floor_pigeonhole_leg
#print axioms witnessMP_le_capRat_of_le_seven
#print axioms witness_floor_pigeonhole_leg_ge_mP
#print axioms small_witness_floor_from_pigeonhole_nodes
#print axioms witness_small_floor_from_bonferroni_nodes
#print axioms witness_large_floor_from_bonferroni_nodes
#print axioms witness_floor_from_bonferroni_nodes_all_clusters
#print axioms lrc14_from_bonferroni_split_nodes
#print axioms lrc14_from_bonferroni_nodes_given_partA
#print axioms lrc14_from_p0_wide_bound_shapes
#print axioms lrc14_from_p0_wide_bound_split_nodes
#print axioms lrc14_from_p0_positive_wide_bound_split_nodes

end Bonferroni
end LRC14
end LonelyRunner
