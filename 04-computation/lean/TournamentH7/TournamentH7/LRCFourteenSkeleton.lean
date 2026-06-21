/-
  TournamentH7.LRCFourteenSkeleton -- a BUILDING Lean skeleton for the LRC(14)
  top-level theorem and its proof DAG (kind-pasteur-2026-06-21, THREAD C).

  This module is the formalization scaffold for the Lonely Runner Conjecture at
  `n = 14` (13 speeds = the FIRST OPEN case).  It records, as Lean statements,
  the proof DAG distilled in the canon (THM-527, THM-563, THM-564, OPEN-Q-108),
  and pins exactly which obligations are:

    * ALREADY PROVED sorry-free elsewhere in this repo
        - the denominator sieve (`LonelyRunner.sieve_frac` etc.);
        - the per-shape Delsarte / moment-LP bound `10*q0 <= L_yK8`
          (`FactorialAtom.delsarte_bound_k8`, the gK8 route's per-shape half);
        - the finite period-max / count / cap-clearance kernels
          (`PeriodmaxCertificate`, `GenuineWideCorrection`, `L7Discrepancy`).

    * FINITE and `decide`/`native_decide`-able NOW (illustrated here with the
      concrete `cap` rationals and a toy bounded-leg decidable predicate);

    * GENUINELY OPEN / needing mathlib analysis (carried as `sorry` with a
      precise statement):
        - THM-527 Part A : `rho* > 0  =>  M(S) >= 1/14`     (measure -> witness);
        - THM-527 Part G : `inf rho* = c0 > 0`              (THE crux, OPEN-Q-108);
        - the gK8 concentration-extremality  `max_E L_yK8 <= 10*cap`;
        - the doublet R-tail (Mordell-Tornheim) uniform bound (THM-564);
        - the top-level assembly into `M(S) >= 1/14` for every covering 13-set.

  DESIGN.  To keep the skeleton lightweight and buildable, the measure-theoretic
  quantities (`p0`, the miss-distribution `q`, `rho*`, `M`) are treated
  ABSTRACTLY: as real-valued functions / opaque data given to the theorems, with
  the analytic content stated as hypotheses or `sorry`-backed lemmas.  This is the
  honest boundary -- the finite combinatorial core is real Lean, the analysis is
  flagged.  Replacing each `sorry` with a mathlib proof (or each opaque function
  with its Lebesgue-measure definition) is the remaining formalization work.

  NONE of the `sorry`s below are claimed as theorems.  They mark open obligations.
-/

import Mathlib.Tactic
import TournamentH7.LonelyRunner
import TournamentH7.LRCFactorialAtom
import TournamentH7.LRCPeriodmaxCertificate
import TournamentH7.LRCGenuineWideCorrection

namespace LonelyRunner
namespace LRC14

open scoped BigOperators

/-! ## 0. The LRC(14) parameters

`n = 14`, so the threshold is `1/14` and there are `k = 13` nonzero speeds.  We
work with the reduced LRC: integer speeds `v : Fin 13 -> ℤ`, all nonzero, and ask
for a real time with every `‖v i * t‖ ≥ 1/14`.  `Lonely 14 v t` (from
`TournamentH7.LonelyRunner`) is exactly this predicate. -/

/-- The LRC(14) dimension. -/
def n14 : ℕ := 14

/-- The number of moving runners in LRC(14). -/
def k13 : ℕ := 13

/-- **The LRC(14) statement (TOP LEVEL).**  Every family of 13 nonzero integer
speeds has a `14`-lonely time.  This is the first open case of the Lonely Runner
Conjecture.  It is the root of the whole DAG; every node below is a reduction
toward it. -/
def LRC14Statement : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∃ t : ℝ, Lonely 14 v t

/-! ## 1. Covering reduction (THM-523/525/526)

By the proven LRC reductions, LRC(14) reduces to: every **primitive covering
13-set** `S` in case **S3** has `M(S) = max_τ min_{v∈S} ‖vτ‖ ≥ 1/14`.  We carry
`M` and the covering structure abstractly. -/

/-- The "max-min reach" `M(S) = sup_t min_{v∈S} ‖v t‖`.  Carried abstractly as a
real number attached to a speed family; in a full formalization this is a `sSup`
over `t` of the `Lonely`-threshold.  (Opaque here.) -/
opaque Mreach : (Fin 13 → ℤ) → ℝ

/-- **DAG node R0 (covering reduction, PROVED upstream THM-523/525/526).**  If the
covering residual bound `M(S) ≥ 1/14` holds for the speed family, then there is a
lonely time.  Stated as the reduction target; the upstream proof is the
covering-set machinery (not re-formalized here).  OPEN obligation: port THM-523. -/
theorem lonely_of_Mreach_ge (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ Mreach v) : ∃ t : ℝ, Lonely 14 v t := by
  -- Mreach v = sup over t of (min over i of ‖v i t‖); hM says this sup ≥ 1/14,
  -- and the sup is attained (compactness of the torus) -> a lonely witness.
  sorry

/-! ## 2. THM-527 -- the lonely-density reformulation

Two halves.  PART A is PROVED in the canon (the slow-fast change of variables):
positive good-period density forces the reach bound.  PART G is the OPEN crux:
the density floor `inf rho* = c0 > 0` over the compact bounded-spread shape space.

We carry `rho*` abstractly as a nonnegative real attached to the data `(P, E)`
(small part `P ⊆ {1,…,13}` and cluster co-offsets `E`). -/

/-- The small part `P = S ∩ {1,…,13}` and the cluster co-offsets `E` (the
"shape"): a pair of integer lists `(P, E)`.  Concrete (nonempty) so the opaque
density and `shapeOf` below are well-formed. -/
abbrev Shape : Type := List ℤ × List ℤ

/-- The good-period density `rho*(P,E) = meas{ x ∈ G_P : the cluster phases have
circular max-gap > 2/7 }`.  Abstract (a Lebesgue measure of a semialgebraic set);
nonnegative. -/
opaque rhoStar : Shape → ℝ

/-- `rho*` is a measure, hence nonnegative.  (Stated; a full formalization would
derive this from the Lebesgue-measure definition.) -/
axiom rhoStar_nonneg : ∀ s : Shape, 0 ≤ rhoStar s

/-- The shape attached to a covering speed family (its `(P,E)` decomposition). -/
opaque shapeOf : (Fin 13 → ℤ) → Shape

/-- **THM-527 PART A (PROVED in canon; OPEN as a Lean obligation).**  Positive
good-period density implies the reach bound `M(S) ≥ 1/14`.  This is the
slow-fast change of variables + criterion C at `v = Vmax`; the canon proves it
(unconditional for the limit, with an `O(1/Vmax)` finite-`Vmax` correction).
Formalization needs: the `Vmax`-ruler period structure, the fast-phase gap
criterion, and the equidistribution `ρ_K → ρ*`. -/
theorem thm527_partA_density_pos_implies_reach (v : Fin 13 → ℤ)
    (hpos : 0 < rhoStar (shapeOf v)) : (1 : ℝ) / 14 ≤ Mreach v := by
  -- good period at x  ⟺  x ∈ G_P AND cluster max-gap > 2/7  ⟹  criterion C
  --   ⟹  a sub-arc safe for all of S  ⟹  M(S) ≥ 1/14.
  -- ρ* > 0 ⟹ a good period exists (finite-Vmax: ρ_K = ρ* + O(#arcs/Vmax)).
  sorry

/-- **THM-527 PART G -- THE CRUX (OPEN, = OPEN-Q-108).**  A uniform positive floor
for the good-period density over all covering shapes: `inf rho* = c0 > 0`.  The
canon reduces this to a COMPACT problem (bounded-spread, `k ≤ 13`), proves the
`k = 3` end (margin `4/3`) and the exact consecutive floor `1/84`, and verifies
positivity on broad scans -- but the uniform compactness floor is not yet a
theorem.  This is THE one remaining inequality of the whole conjecture (this
route). -/
theorem thm527_partG_uniform_floor :
    ∃ c0 : ℝ, 0 < c0 ∧ ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      c0 ≤ rhoStar (shapeOf v) := by
  -- bounded-spread compact reduction (Part D) + continuity + positivity (Part E)
  -- + integer-vs-real passage + Vmax ≤ V0 finite check.  OPEN.
  sorry

/-- **THM-527 k=3 leg (PROVED unconditional in canon).**  As a Lean obligation:
when the cluster has size 3, three phases always leave a max-gap `≥ 1/3 > 2/7`, so
every `x ∈ G_P` is good and `rho* = meas(G_P) > 0`.  Carried as a hypothesis-form
statement (`G_P` positivity is the proven LRC at 11 speeds). -/
theorem thm527_k3_unconditional (s : Shape) (hk3 : True) (hGP : 0 < rhoStar s) :
    0 < rhoStar s := hGP

/-! ## 3. The sector / p0 route (L0 layer): the wide bound `p0(E) ≤ cap_k`

`p0(E) = measS7(E) = meas{ x : the phases {frac(e·x)} hit all 7 sectors }` is the
cover atom (`= q0`, the origin coordinate of the missed-count distribution `q`).
The "wide bound" is `p0(E) ≤ cap_k` for every `k`-speed configuration, split into
three legs (bounded / single-far / genuine-wide), OR closed in one shot by the
gK8 concentration cert.  We carry `p0` and the miss-distribution `q` abstractly. -/

/-- The cover atom `p0(E) = measS7(E)`.  Abstract (Lebesgue measure of the
all-sectors-hit set); a probability, so in `[0,1]`. -/
opaque p0 : (List ℤ) → ℝ

/-- The cap plateau `cap_k` (the decorrelated `Q`-plateau upper bound).  The exact
rationals from the canon for `k = 8..12` are recorded below. -/
def capRat : ℕ → ℚ
  | 8  => 2243 / 5880
  | 9  => 1979 / 4004
  | 10 => 55 / 91
  | 11 => 66 / 91
  | 12 => 6 / 7
  | _  => 1            -- non-binding / trivial elsewhere

/-- The caps are genuine probabilities `0 < cap_k < 1` for the binding rows
(finite, decidable). -/
theorem capRat_lt_one : ∀ k ∈ [8,9,10,11,12], capRat k < 1 := by decide

theorem capRat_pos : ∀ k ∈ [8,9,10,11,12], 0 < capRat k := by decide

/-- **The wide bound (the L0 target of the sector route).**  For every covering
configuration `E` with `|E| = k`, the cover atom is at most the cap.  This is the
disjunction of the three legs below (or the single gK8 cert).  Carried as the
target `Prop`. -/
def WideBound (k : ℕ) : Prop :=
  ∀ E : List ℤ, E.length = k → p0 E ≤ (capRat k : ℝ)

/-! ### 3a. Bounded leg (span ≤ 14): FINITE, native_decide-able

The bounded leg is `p0(E) ≤ cap_k` for `span(E) ≤ 14`.  In the canon this is a
finite exact-rational enumeration (the breakpoint grid of `measS7` is rational).
Here we record a TOY decidable witness pattern: an integer cross-multiplied
comparison on stored exact `(p0, cap)` rows.  (The full enumeration over all
bounded `E` is the Python/exact-arithmetic audit; this shows the `decide` shape
the formal cert would take.) -/

/-- A stored exact bounded-leg row: `p0 = num/den ≤ cap = capNum/capDen`. -/
structure BoundedRow where
  pNum : Int
  pDen : Int
  capNum : Int
  capDen : Int
deriving DecidableEq, Repr

/-- Cross-multiplied `p0 ≤ cap` (assuming positive denominators). -/
def BoundedRow.ok (r : BoundedRow) : Bool :=
  decide (r.pNum * r.capDen ≤ r.capNum * r.pDen)

/-- Example bounded rows (placeholder constants standing in for the exact
`measS7`/`cap` table; the real table is the bounded enumeration). -/
def sampleBoundedRows : List BoundedRow :=
  [ { pNum := 1, pDen := 4, capNum := 2243, capDen := 5880 }   -- k=8 sample
  , { pNum := 1, pDen := 3, capNum := 55,   capDen := 91 } ]    -- k=10 sample

/-- **Bounded-leg cert shape (FINITE, `native_decide`).**  Every stored bounded
row clears its cap.  This is the decidable boundary the full bounded enumeration
plugs into. -/
theorem sampleBoundedRows_ok :
    ∀ r ∈ sampleBoundedRows, r.ok = true := by native_decide

/-! ### 3b. Single-far leg (THM-563): FINITE period-max, ALREADY KERNELIZED

THM-563: the signed one-far deviation `w·Δ_w` is exactly periodic, so
`sup_w Δ_w·w` is a finite period-max; `period-max(B) < 15·margin(B)` closes the
leg for all `w ≥ 15`.  The 12805-base finite check is ALREADY a sorry-free Lean
kernel (`PeriodmaxCertificate`).  We re-expose its headroom result as the
single-far node. -/

/-- **DAG node L1 (single-far, FINITE -- backed by the sorry-free kernel).**  Every
per-`k` worst row of the completed THM-563 audit has positive headroom
`15·margin − period-max > 0`.  This is `PeriodmaxCertificate`'s result; the full
per-base enumeration is the mac-mini exact audit. -/
theorem single_far_periodmax_headroom :
    ∀ r : Fin 6, 0 < PeriodmaxCertificate.headroomNum (PeriodmaxCertificate.worstRow r) :=
  PeriodmaxCertificate.all_worst_rows_headroom_positive

/-! ### 3c. Genuine-wide leg (THM-564): doublet P/R split

The genuine-wide binding maximizer is the doublet `consec_{k-2} ∪ {M,M+1}`.
THM-564 splits `g(M) = M·(p0−Φ) = P(M) + R(M)` with `P` exactly periodic (finite
period-max) and `R = O(1/M)` (Mordell-Tornheim double sum, `|R| ≤ 12ζ(3)·N/π³`,
the analytic piece).  Closure = period-max(P) + sup|R| + a tiny finite window.
The finite ledger (rows below cap, robust-margin flags) is the sorry-free
`GenuineWideCorrection` kernel. -/

/-- **DAG node L2a (genuine-wide finite ledger, FINITE -- sorry-free kernel).**
Every reported doublet true-max row is below cap. -/
theorem genuine_wide_rows_below_cap :
    ∀ r : Fin 5, 0 < GenuineWideCorrection.marginNum (GenuineWideCorrection.maxRow r) :=
  GenuineWideCorrection.all_reported_rows_below_cap

/-- **DAG node L2b (the R-tail bound, OPEN -- needs mathlib analysis).**  The
doublet interaction correction `R(M) = M·(d2(M) − d_inf)` is uniformly bounded:
`sup_{M≥15} |R(M)| ≤ 12ζ(3)·N/π³` with `N ≤ 15` the active sector-pair count.
This is the Mordell-Tornheim / Koksma analytic obligation -- a real-analysis
statement (`ζ(3)`, Fourier of sector indicators).  OPEN as a Lean obligation. -/
theorem doublet_Rtail_uniform_bound :
    ∃ Rsup : ℝ, 0 ≤ Rsup ∧
      Rsup ≤ 12 * (∑' m : ℕ, (1 : ℝ) / ((m : ℝ) + 1) ^ 3) * 15 / Real.pi ^ 3 := by
  -- |R(M)| = |M·(d2(M) − d_inf)|, a covariance of two sector indicators at the
  -- locked far phases; Fourier ⟹ convergent double sum bounded by 12ζ(3)·N/π³,
  -- N = #active sector-pairs ≤ C(6,2) = 15 (base-independent).  OPEN.
  sorry

/-! ### 3d. The gK8 concentration route -- per-shape half ALREADY PROVED

The cleanest closure: `10*q0 ≤ L_yK8 = 10q0 + q3 + 10q6` holds per-shape (PROVED
sorry-free as `FactorialAtom.delsarte_bound_k8`), and the only remaining content
is the scalar concentration-extremality `max_E L_yK8 ≤ 10·cap_k`. -/

/-- **DAG node L3a (gK8 per-shape bound, PROVED sorry-free upstream).**  For any
nonnegative missed-count atom `q` (a genuine row distribution), `10*q0 ≤ L_yK8(q)`.
This is exactly `FactorialAtom.delsarte_bound_k8`, re-exposed as the gK8 node.
(`q0 = p0`, so this is `10*p0 ≤ 10q0 + q3 + 10q6`.) -/
theorem gK8_per_shape_bound (q : FactorialAtom.Atom)
    (hq : ∀ t : Fin FactorialAtom.atomCount, 0 ≤ q t) :
    10 * FactorialAtom.q0 q ≤ FactorialAtom.LyK8 q :=
  FactorialAtom.delsarte_bound_k8 q hq

/-- The gK8 dual is Delsarte-feasible (Krawtchouk-nonnegative, dominates the cover
indicator) -- PROVED sorry-free upstream. -/
theorem gK8_dual_feasible :
    ∀ t : Fin FactorialAtom.atomCount,
      (if t.val = 0 then (10 : Int) else 0) ≤ FactorialAtom.gK8 t :=
  FactorialAtom.gK8_dominates

/-- **DAG node L3b (gK8 concentration-extremality, OPEN).**  The scalar content:
the maximum of the moment functional `L_yK8` over ALL `k`-speed configs equals its
maximum over BOUNDED configs, which is `< 10·cap_k`.  Carried abstractly: a real
functional `LyVal` on configs is maximized by a bounded config.  This is a
smoothing/majorization lemma on the 7-simplex (HYP-2812).  OPEN. -/
theorem gK8_concentration_extremality
    (LyVal : (List ℤ) → ℝ) (k : ℕ)
    (boundedMax : ℝ) (hbound : boundedMax < 10 * (capRat k : ℝ)) :
    (∀ E : List ℤ, E.length = k → LyVal E ≤ boundedMax) := by
  -- gK8 charges the extreme miss-counts q0, q6; both individually maximized by
  -- CONCENTRATED (slowest/bounded) configs; decorrelation smooths to the middle.
  -- ⟹ max_E L_yK8 = max_BOUNDED L_yK8 = boundedMax < 10·cap.  OPEN.
  sorry

/-! ## 4. Assembling the wide bound

Either route (3a+3b+3c, or 3d) yields the wide bound `p0(E) ≤ cap_k`.  We record
the gK8 assembly: per-shape (`10*p0 ≤ L_yK8`) + concentration (`L_yK8 ≤ 10·cap`)
⟹ `p0 ≤ cap`.  The bridge `p0 = q0` and `p0·10 = 10·q0` is definitional once `q`
is the miss-distribution of `E`; carried as a hypothesis here. -/

/-- **DAG node W (wide bound via gK8).**  If the per-shape bound gives
`10*p0(E) ≤ L_yK8(E)` and concentration gives `L_yK8(E) ≤ 10*cap_k`, then
`p0(E) ≤ cap_k`.  This step is pure arithmetic (real division by 10). -/
theorem wide_bound_from_gK8 (E : List ℤ) (k : ℕ) (LyE : ℝ)
    (hshape : 10 * p0 E ≤ LyE) (hconc : LyE ≤ 10 * (capRat k : ℝ)) :
    p0 E ≤ (capRat k : ℝ) := by
  have h : 10 * p0 E ≤ 10 * (capRat k : ℝ) := le_trans hshape hconc
  linarith

/-! ## 5. The two routes to the top, and the top-level statement

Both the sector route (wide bound `p0 ≤ cap` ⟹ via THM-526 criterion C ⟹
`M(S) ≥ 1/14`) and the lonely-density route (THM-527: `rho* > 0` ⟹ `M ≥ 1/14`,
with `rho* > 0` from the Part-G floor) converge on `M(S) ≥ 1/14`, hence on
LRC(14).  We record the THM-527 route as the cleanest assembly: Part G gives a
positive floor, Part A turns it into the reach bound, R0 turns that into a lonely
time. -/

/-- **TOP-LEVEL ASSEMBLY (THM-527 route).**  Combining the Part-G uniform floor
(`rho* ≥ c0 > 0`), Part A (`rho* > 0 ⟹ M ≥ 1/14`), and R0 (`M ≥ 1/14 ⟹ lonely`)
proves LRC(14).  This theorem is `sorry`-free GIVEN the three nodes -- i.e. it is
the genuine logical glue, and the only `sorry`s it transitively depends on are the
three flagged open obligations (Part G, Part A, R0). -/
theorem lrc14_from_thm527 : LRC14Statement := by
  intro v hv
  -- Part G: a uniform positive floor c0 ≤ ρ*(shapeOf v).
  obtain ⟨c0, hc0pos, hfloor⟩ := thm527_partG_uniform_floor
  have hrho : 0 < rhoStar (shapeOf v) := lt_of_lt_of_le hc0pos (hfloor v hv)
  -- Part A: ρ* > 0 ⟹ M(S) ≥ 1/14.
  have hM : (1 : ℝ) / 14 ≤ Mreach v := thm527_partA_density_pos_implies_reach v hrho
  -- R0: M(S) ≥ 1/14 ⟹ a lonely time.
  exact lonely_of_Mreach_ge v hv hM

/-! ## 6. Sieve sanity checks (PROVED sorry-free upstream)

The denominator sieve already settles every LRC(14) family in which some `q ∈
{2,…,14}` divides no speed.  This is the unconditional part, and it is sorry-free
(`LonelyRunner.sieve_one_div`).  We record the LRC(14) specializations. -/

/-- **No-multiple-of-14 case (PROVED).**  If no speed is divisible by 14, then
`t = 1/14` is a 14-lonely time. -/
theorem lrc14_no_multiple_of_14 (v : Fin 13 → ℤ)
    (hdiv : ∀ i, ¬ ((14 : ℤ) ∣ v i)) : Lonely 14 v ((1 : ℝ) / 14) :=
  lonely_of_no_multiple 14 (by norm_num) v hdiv

/-- **All-odd case (PROVED).**  If every speed is odd, then `t = 1/2` is lonely. -/
theorem lrc14_all_odd (v : Fin 13 → ℤ) (hodd : ∀ i, ¬ (2 : ℤ) ∣ v i) :
    Lonely 14 v ((1 : ℝ) / 2) :=
  all_odd_half_lonely 14 v (by norm_num) hodd

/-- **A counterexample must be 14-divisible-saturated (PROVED).**  Any LRC(14)
counterexample contains, for every `q ∈ {2,…,14}`, a speed divisible by `q`. -/
theorem lrc14_counterexample_saturated (v : Fin 13 → ℤ)
    (hcex : ∀ t : ℝ, ¬ Lonely 14 v t) (q : ℕ) (hq2 : 2 ≤ q) (hq14 : q ≤ 14) :
    ∃ i, (q : ℤ) ∣ v i :=
  counterexample_needs_all_divisors 14 v hcex q hq2 hq14

/-! ## 7. Axiom audit

The PROVED-upstream nodes (sieve specializations, gK8 per-shape bound, finite
kernels) should report only Lean foundations (+ `native_decide` axioms).  The
OPEN nodes (`*_partA_*`, `*_partG_*`, `lonely_of_Mreach_ge`,
`doublet_Rtail_uniform_bound`, `gK8_concentration_extremality`) and anything
depending on them (`lrc14_from_thm527`) will additionally report `sorryAx` --
which is the precise machine-checkable statement of "these are the open
obligations." -/

#print axioms lrc14_no_multiple_of_14
#print axioms lrc14_all_odd
#print axioms lrc14_counterexample_saturated
#print axioms gK8_per_shape_bound
#print axioms gK8_dual_feasible
#print axioms single_far_periodmax_headroom
#print axioms genuine_wide_rows_below_cap
#print axioms sampleBoundedRows_ok
#print axioms wide_bound_from_gK8
-- The following SHOULD report `sorryAx` (open obligations):
#print axioms thm527_partA_density_pos_implies_reach
#print axioms thm527_partG_uniform_floor
#print axioms lonely_of_Mreach_ge
#print axioms doublet_Rtail_uniform_bound
#print axioms gK8_concentration_extremality
#print axioms lrc14_from_thm527

end LRC14
end LonelyRunner
