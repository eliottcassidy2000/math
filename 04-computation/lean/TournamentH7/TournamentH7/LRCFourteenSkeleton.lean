/-
  TournamentH7.LRCFourteenSkeleton -- a BUILDING Lean skeleton for the LRC(14)
  top-level theorem and its proof DAG (kind-pasteur-2026-06-21, THREAD C).

  This module is the formalization scaffold for the Lonely Runner Conjecture at
  `n = 14` (13 speeds = the FIRST OPEN case).  It records, as Lean statements,
  the proof DAG distilled in the canon (THM-527, THM-563, THM-564, OPEN-Q-108),
  and pins exactly which obligations are:

    * ALREADY PROVED sorry-free elsewhere in this repo
        - the denominator sieve (`LonelyRunner.sieve_frac` etc.);
        - the concrete `Mreach` compactness handoff
          (`LRC14Concrete.lonely_of_Mreach_ge`);
        - the per-shape Delsarte / moment-LP bound `10*q0 <= L_yK8`
          (`FactorialAtom.delsarte_bound_k8`, the gK8 route's per-shape half);
        - the finite period-max / count / cap-clearance kernels
          (`PeriodmaxCertificate`, `GenuineWideCorrection`,
          `Gk8SingleFar`, `DoubletWitnessFloor`, `L7Discrepancy`).

    * FINITE and `decide`/`native_decide`-able NOW (illustrated here with the
      concrete `cap` rationals and a toy bounded-leg decidable predicate);

    * GENUINELY OPEN / needing mathlib analysis (carried as named `Prop`
      obligations, not as asserted theorems):
        - THM-527 Part A : positive witness density => M(S) >= 1/14
          (measure -> witness);
        - the 1/7 witness floor `G2 >= m_P` over the remaining k=8..13 shapes
          (the new easier witness route; k<=7 is pigeonhole-elementary);
        - the corrected global-witness/rho-glob floor route at threshold 1/7;
        - the gK8 concentration-extremality  `max_E L_yK8 <= 10*cap`;
        - the doublet R-tail (Mordell-Tornheim) uniform bound (THM-564);
        - the top-level assembly into `M(S) >= 1/14` for every covering 13-set.

  DESIGN.  To keep the skeleton lightweight and buildable, the remaining
  measure-theoretic quantities (the miss-distribution `q`, `rho*`, and the
  decorrelation residuals) are treated ABSTRACTLY: as real-valued functions /
  opaque data given to the theorems, with the analytic content stated as
  hypotheses or named obligations.  The cover atom `p0` is now concrete via
  `DenseCovers.p0 = (slowμ (coverSet E)).toReal`.  The max-min reach `M`
  is now concrete: it is the supremum of the finite min-reach function over
  `[0,1]`, imported from `LRCMreachConcrete`.  This is the honest boundary --
  the finite combinatorial core and compactness handoff are real Lean, while the
  remaining measure/equidistribution analysis is flagged.  Replacing each
  remaining obligation with a mathlib proof (or each remaining opaque function
  with its Lebesgue-measure definition) is the remaining formalization work.

  No declaration below asserts an unproved analytic fact.  The open analytic
  targets are `Prop` definitions consumed only as hypotheses by conditional
  assembly theorems.
-/

import Mathlib.Tactic
import TournamentH7.LonelyRunner
import TournamentH7.LRCMreachConcrete
import TournamentH7.LRCFactorialAtom
import TournamentH7.LRCPeriodmaxCertificate
import TournamentH7.LRCGenuineWideCorrection
import TournamentH7.LRCGk8SingleFar
import TournamentH7.LRCDoubletWitnessFloor
import TournamentH7.LRCP0Concrete

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
`M` concretely and the covering structure abstractly. -/

/-- The concrete "max-min reach" `M(S) = sup_t min_{v∈S} ‖v t‖`, formalized as
the supremum of the finite min-reach function over `[0,1]`. -/
noncomputable abbrev Mreach : (Fin 13 → ℤ) → ℝ := LRC14Concrete.Mreach

/-- **DAG node R0.**  If the concrete max-min reach `M(S) ≥ 1/14` holds for the
speed family, then there is a lonely time.  This is discharged by the compactness
theorem in `LRCMreachConcrete`: the continuous min-reach function attains its
supremum on `[0,1]`. -/
theorem lonely_of_Mreach_ge (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ Mreach v) : ∃ t : ℝ, Lonely 14 v t := by
  exact LRC14Concrete.lonely_of_Mreach_ge v hv hM

/-! ## 2. THM-527 -- density reformulations

There are now two related density routes:

* the older via-`Vmax` `rho*` route at the `2/7` max-gap threshold;
* the corrected global-witness/rho-glob route at the `1/7` max-gap threshold,
  represented here by `G2`.

The current repo canon treats the `G2` route as the sharper final target.  The
k<=7 binding cases are pigeonhole-elementary; the remaining `k=8..13` floor has
large computational slack.  KPS Thread A later found admissible zeros for the
old via-`Vmax` `rho*` object, so its uniform positive floor is recorded only as
an obsolete conditional proposition below, not as a theorem.  The corrected
object is the global 1/7 witness density.

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

/-- `rho*` is a measure, hence nonnegative.  A full formalization should derive
this from the Lebesgue-measure definition; here it is only recorded as an
obligation proposition. -/
def rhoStar_nonneg_obligation : Prop := ∀ s : Shape, 0 ≤ rhoStar s

/-- The shape attached to a covering speed family (its `(P,E)` decomposition). -/
opaque shapeOf : (Fin 13 → ℤ) → Shape

/-- **THM-527 PART A (PROVED in canon; OPEN as a Lean obligation).**  Positive
good-period density implies the reach bound `M(S) ≥ 1/14`.  This is the
slow-fast change of variables + criterion C at `v = Vmax`; the canon proves it
(unconditional for the limit, with an `O(1/Vmax)` finite-`Vmax` correction).
Formalization needs: the `Vmax`-ruler period structure, the fast-phase gap
criterion, and the equidistribution `ρ_K → ρ*`.

This is deliberately a `Prop` obligation, not a theorem: the old 2/7
via-`Vmax` object has admissible zero examples, so downstream code must pass an
actual Part-A proof if it wants to use this obsolete route. -/
def thm527_partA_density_pos_implies_reach : Prop :=
  ∀ v : Fin 13 → ℤ, 0 < rhoStar (shapeOf v) → (1 : ℝ) / 14 ≤ Mreach v

/-- **OBSOLETE THM-527 PART G candidate.**  A uniform positive floor for the
via-`Vmax` 2/7 good-period density.  KPS Thread A now gives admissible
counterexamples with this density equal to zero, so this proposition is retained
only to keep the failed route explicit.  It must not be used as an asserted
theorem in the current proof DAG. -/
def obsoleteRhoStarUniformFloor : Prop :=
    ∃ c0 : ℝ, 0 < c0 ∧ ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      c0 ≤ rhoStar (shapeOf v)

/-- **THM-527 k=3 leg (PROVED unconditional in canon).**  As a Lean obligation:
when the cluster has size 3, three phases always leave a max-gap `≥ 1/3 > 2/7`, so
every `x ∈ G_P` is good and `rho* = meas(G_P) > 0`.  Carried as a hypothesis-form
statement (`G_P` positivity is the proven LRC at 11 speeds). -/
theorem thm527_k3_unconditional (s : Shape) (_hk3 : True) (hGP : 0 < rhoStar s) :
    0 < rhoStar s := hGP

/-! ### 2a. Doublet rho*/witness scout checksum (finite arithmetic import)

The S4 genuine-wide doublet scout reports exact positive floors over the audited
doublet bank: `rho* >= 2/147` and witness `>= 1066/2205`.  The measurement engine
is Python/exact-rational; this Lean node records only the final rational
comparisons and keeps the coverage claim out of Lean. -/

/-- **DAG node D0 (doublet witness-floor checksum, FINITE arithmetic import).**
The audited doublet floors are positive, and the reported `rho*`/witness minima
are above the older comparison floors `1/84` and `14249/252252`. -/
theorem doublet_witness_floor_checksum :
    0 < DoubletWitnessFloor.rhoFloorNum ∧
      0 < DoubletWitnessFloor.witnessFloorNum ∧
      DoubletWitnessFloor.consecRhoFloorNum * DoubletWitnessFloor.rhoFloorDen <
        DoubletWitnessFloor.rhoFloorNum * DoubletWitnessFloor.consecRhoFloorDen ∧
      DoubletWitnessFloor.thm530WitnessFloorNum *
          DoubletWitnessFloor.witnessFloorDen <
        DoubletWitnessFloor.witnessFloorNum *
          DoubletWitnessFloor.thm530WitnessFloorDen ∧
      DoubletWitnessFloor.rhoFloorNum * DoubletWitnessFloor.witnessFloorDen <
        DoubletWitnessFloor.witnessFloorNum * DoubletWitnessFloor.rhoFloorDen :=
  DoubletWitnessFloor.audited_doublet_floors_positive_and_separated

/-! ### 2b. The 1/7 witness route (KPS HYP-2825/HYP-2826/HYP-2827)

`G2(P,E)` is the measure of good slow phases in `G_P` whose cluster phases leave a
circular gap `> 1/7`.  This is the direct witness scale: a positive `G2` produces
`M(S) >= 1/14`.  KPS HYP-2827 proves the binding `k<=7` cases by the elementary
pigeonhole bound `maxgap >= 1/k >= 1/7`; only k=8..13 remains, with large slack.
-/

/-- The 1/7-scale witness density
`G2(P,E)=meas{x in G_P : maxgap(frac(e_i*x)) > 1/7}`.  Abstract here; a full
formalization will define it as a Lebesgue measure. -/
opaque witnessG2 : Shape → ℝ

/-- Current cluster size parameter: the number of co-offsets in the shape. -/
def clusterSize (s : Shape) : ℕ := s.2.length

/-- Alias for KPS's corrected `rho*_glob` notation: the global-witness density is
the same 1/7 object called `witnessG2` in this skeleton. -/
abbrev rhoGlob : Shape → ℝ := witnessG2

/-- The exact admissible small-part floor from THM-530/KPS witness route:
`m_P = 14249/252252`. -/
def witnessMP : ℚ := 14249 / 252252

/-- Exact lower bounds for the corrected global-witness density `rhoGlob=G2` in
the k=8..13 branch, as reported by the KPS/Claude exact rational engines.  These
are not definitions of the measure; they are the Lean-facing arithmetic ledger
for the remaining witness-floor proof once the measure lower-bound theorem is
formalized. -/
def rhoGlobFloorRat : ℕ → ℚ
  | 8 => 8152 / 24255
  | 9 => 143 / 420
  | 10 => 19 / 49
  | 11 => 3193 / 8820
  | 12 => 10358 / 24255
  | 13 => 477 / 1078
  | _ => 0

/-- The proved floor is strictly positive, as rational arithmetic. -/
theorem witnessMP_pos_rat : (0 : ℚ) < witnessMP := by
  native_decide

/-- The same strict positivity after casting to `ℝ`. -/
theorem witnessMP_pos_real : (0 : ℝ) < (witnessMP : ℝ) := by
  norm_num [witnessMP]

/-- If a shape has the `m_P` witness floor, then its witness density is positive.
This is the tiny arithmetic bridge that turns the exact floor into the THM-527
Part-A input. -/
theorem witness_floor_positive (s : Shape)
    (hfloor : (witnessMP : ℝ) ≤ witnessG2 s) : 0 < witnessG2 s :=
  lt_of_lt_of_le witnessMP_pos_real hfloor

/-- The exact k=8..13 `rhoGlob` floor ledger is uniformly stronger than the
admissible `m_P` floor needed by the witness route. -/
theorem witnessMP_le_rhoGlobFloorRat_real
    (k : ℕ) (h8 : 8 ≤ k) (h13 : k ≤ 13) :
    (witnessMP : ℝ) ≤ (rhoGlobFloorRat k : ℝ) := by
  interval_cases k <;> norm_num [witnessMP, rhoGlobFloorRat]

/-- Pure arithmetic core of the KPS HYP-2827 pigeonhole step: for any nonempty
cluster of size `k<=7`, the pigeonhole lower bound `maxgap >= 1/k` is already at
the 1/7 witness threshold.  The geometric max-gap theorem is separate; this is
the exact threshold arithmetic. -/
theorem one_seventh_le_inv_of_pos_le_seven (k : ℕ) (hkpos : 0 < k) (hk : k ≤ 7) :
    (1 : ℝ) / 7 ≤ (1 : ℝ) / k := by
  interval_cases k <;> norm_num at hkpos ⊢

/-- No-sorry bookkeeping for the corrected witness-floor proof split.  If the
small cluster cases `k<=7` have the pigeonhole floor, and the remaining
`8<=k<=13` cases have the computational/analytic floor, then the global
`G2>=m_P` floor follows. -/
theorem witness_floor_from_cluster_cases
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hlarge : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v hv
  by_cases hle : clusterSize (shapeOf v) ≤ 7
  · exact hsmall v hv hle
  · have h8 : 8 ≤ clusterSize (shapeOf v) := by omega
    exact hlarge v hv h8 (hsize v hv)

/-- No-sorry bridge from the exact `rhoGlobFloorRat` ledger to the abstract
large-cluster floor obligation.  A future formal proof only has to show that
each k=8..13 shape has `rhoGlobFloorRat k <= rhoGlob`; the comparison with the
global `m_P` floor is already checked here by exact rational arithmetic. -/
theorem witness_large_floor_from_rhoGlob_lower_bound
    (hlower : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (rhoGlobFloorRat (clusterSize (shapeOf v)) : ℝ) ≤ rhoGlob (shapeOf v)) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v hv h8 h13
  exact le_trans
    (witnessMP_le_rhoGlobFloorRat_real (clusterSize (shapeOf v)) h8 h13)
    (hlower v hv h8 h13)

/-- Sorry-free logical glue for the new witness route, parameterized by the two
analysis nodes: the uniform `G2>=m_P` floor and the direct-witness implication
`G2>0 -> Mreach>=1/14`.  This is the route now preferred over the old `rho*`
floor. -/
theorem lrc14_from_witness_floor_given_nodes
    (hfloor : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v)
    (hR0 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (1 : ℝ) / 14 ≤ Mreach v → ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  intro v hv
  have hG2pos : 0 < witnessG2 (shapeOf v) :=
    witness_floor_positive (shapeOf v) (hfloor v hv)
  exact hR0 v hv (hpartA v hG2pos)

/-- **TOP-LEVEL ASSEMBLY (1/7 witness route).**  This specializes the previous
glue to the abstract `Mreach` reduction already carried in the skeleton.  It
still depends on the supplied direct-witness implication `G2>0 -> Mreach>=1/14`,
but the reach-to-lonely and witness-floor-to-top logical wiring are explicit and
checked. -/
theorem lrc14_from_witness_floor
    (hfloor : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_witness_floor_given_nodes hfloor hpartA lonely_of_Mreach_ge

/-- Same top-level witness-route glue, but with the floor split into the current
two proof obligations: `k<=7` pigeonhole and `8<=k<=13` floor. -/
theorem lrc14_from_witness_floor_cases_given_nodes
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hlarge : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v)
    (hR0 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (1 : ℝ) / 14 ≤ Mreach v → ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement :=
  lrc14_from_witness_floor_given_nodes
    (witness_floor_from_cluster_cases hsmall hlarge hsize) hpartA hR0

/-! ## 3. The sector / p0 route (L0 layer): the wide bound `p0(E) ≤ cap_k`

`p0(E) = measS7(E) = meas{ x : the phases {frac(e·x)} hit all 7 sectors }` is the
cover atom (`= q0`, the origin coordinate of the missed-count distribution `q`).
The "wide bound" is `p0(E) ≤ cap_k` for every `k`-speed configuration, split into
three legs (bounded / single-far / genuine-wide), OR closed in one shot by the
gK8 concentration cert.  We carry the miss-distribution `q` abstractly, while
`p0` is now the concrete `DenseCovers.p0`. -/

/-- The cover atom `p0(E) = measS7(E)`.  **Concrete** (kps-S31b): the Lebesgue
measure of the all-6-inner-sectors-hit event `coverSet E` on the slow-time
probability space, `(slowμ (coverSet E)).toReal` — a probability in `[0,1]`.  The
elementary cover-bound cores (`p0_nonneg`, `p0_le_one`, `p0_eq_zero_of_card_lt_six`,
`p0_mono`, `wideBoundConcrete_of_decorrelation`) live in `LRCP0Concrete`. -/
noncomputable def p0 : (List ℤ) → ℝ := DenseCovers.p0

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
theorem capRat_lt_one : ∀ k ∈ [8,9,10,11,12], capRat k < 1 := by native_decide

theorem capRat_pos : ∀ k ∈ [8,9,10,11,12], 0 < capRat k := by native_decide

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
This is the Mordell-Tornheim / Koksma analytic obligation, parameterized by the
actual tail function so the statement is not vacuous. -/
def doublet_Rtail_uniform_bound (Rtail : ℕ → ℝ) : Prop :=
  ∃ Rsup : ℝ, 0 ≤ Rsup ∧
    (∀ M : ℕ, 15 ≤ M → |Rtail M| ≤ Rsup) ∧
    Rsup ≤ 12 * (∑' m : ℕ, (1 : ℝ) / ((m : ℝ) + 1) ^ 3) * 15 / Real.pi ^ 3

/-! ### 3d. The gK8 concentration route -- per-shape half ALREADY PROVED

The cleanest closure: `10*q0 ≤ L_yK8 = 10q0 + q3 + 10q6` holds per-shape (PROVED
sorry-free as `FactorialAtom.delsarte_bound_k8`), and the only remaining content
is the scalar concentration-extremality `max_E L_yK8 ≤ 10·cap_k`.

Mac-mini HYP-2829 now refines this open node into a far-count reduction:
bounded `r=0` is finite, the binding wide case is single-far `r=1` and should
route through THM-563 periodicity, and `r>=2` has large decorrelation margin.
That reduction is not formalized here, but it is the current intended structure
of this obligation. -/

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

/-- **DAG node L3a' (HYP-2829 single-far arithmetic, PROVED sorry-free).**
The exact rational inequalities for the bounded and single-far gK8 binding rows
at k=8,9,10.  This is the arithmetic checksum under the planned period-max
proof of the single-far sup; it does not by itself prove the analytic
periodicity/window reduction. -/
theorem gK8_singlefar_arithmetic_kernel :
    ((2633 : Nat) * 588 < 2243 * 735) ∧
    ((3259 : Nat) * 2002 < 9895 * 735) ∧
    ((37 : Nat) * 7 < 40 * 7) ∧
    ((2323 : Nat) * 588 < 2243 * 980) ∧
    ((2876 : Nat) * 2002 < 9895 * 735) ∧
    ((62267 : Nat) * 7 < 40 * 12936) ∧
    ((2323 : Nat) * 735 < 2633 * 980) ∧
    ((2876 : Nat) < 3259) ∧
    ((62267 : Nat) * 7 < 37 * 12936) :=
  Gk8SingleFar.all_binding_checks

/-- **DAG node L3a' (HYP-2829 single-far arithmetic checksum, FINITE import).**
Alias kept for S81 handoffs.  The bounded binding values and the reported
single-far window suprema clear `10*cap` for k=8,9,10, and the single-far
values are below the bounded values. -/
theorem gK8_singlefar_arithmetic_checksum :
    (2633 : Nat) * 588 < 2243 * 735 ∧
      (3259 : Nat) * 2002 < 9895 * 735 ∧
      (37 : Nat) * 7 < 40 * 7 ∧
      (2323 : Nat) * 588 < 2243 * 980 ∧
      (2876 : Nat) * 2002 < 9895 * 735 ∧
      (62267 : Nat) * 7 < 40 * 12936 ∧
      (2323 : Nat) * 735 < 2633 * 980 ∧
      (2876 : Nat) < 3259 ∧
      (62267 : Nat) * 7 < 37 * 12936 :=
  gK8_singlefar_arithmetic_kernel

/-- **DAG node L3b (gK8 concentration-extremality, OPEN).**  The scalar content:
the maximum of the moment functional `L_yK8` over ALL `k`-speed configs is below
`10·cap_k`.  HYP-2829's preferred split is:

* `r=0` bounded configs: finite certificate;
* `r=1` single-far configs: THM-563 periodicity, apparently the binding wide
  case;
* `r>=2`: decorrelated/far-count margin.

Carried abstractly as a real functional `LyVal` on configs.  This is an
obligation `Prop`, not an asserted theorem, because arbitrary `LyVal` cannot be
bounded without the missing concentration/far-count hypotheses. -/
def gK8_concentration_extremality
    (LyVal : (List ℤ) → ℝ) (k : ℕ) (boundedMax : ℝ) : Prop :=
  boundedMax < 10 * (capRat k : ℝ) ∧
    ∀ E : List ℤ, E.length = k → LyVal E ≤ boundedMax

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

/-! ## 5. The routes to the top, and the top-level statement

The sector route (wide bound `p0 ≤ cap` ⟹ via THM-526 criterion C ⟹
`M(S) ≥ 1/14`) and the corrected global-witness route (`G2>0 ⟹ M ≥ 1/14`, with
`G2>0` from the 1/7 floor) converge on `M(S) ≥ 1/14`, hence on LRC(14).

For auditability we also keep the logical glue for the obsolete 2/7 `rhoStar`
route as a **conditional theorem**: if the now-refuted uniform floor were true,
then the old route would imply LRC(14).  No theorem below asserts that floor. -/

/-- **CONDITIONAL OLD ROUTE.**  If the obsolete 2/7 via-`Vmax` uniform floor and
the Part-A implication were available, Part A plus R0 would imply LRC(14).  The
incoming KPS zero probes show that the uniform-floor hypothesis is false for the
intended object, so the current proof should use
`lrc14_from_witness_floor_cases_given_nodes` instead. -/
theorem lrc14_from_obsolete_rhoStar_floor
    (hpartA : thm527_partA_density_pos_implies_reach)
    (hobsolete : obsoleteRhoStarUniformFloor) : LRC14Statement := by
  intro v hv
  obtain ⟨c0, hc0pos, hfloor⟩ := hobsolete
  have hrho : 0 < rhoStar (shapeOf v) := lt_of_lt_of_le hc0pos (hfloor v hv)
  have hM : (1 : ℝ) / 14 ≤ Mreach v := hpartA v hrho
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
kernels, concrete Mreach handoff, and conditional glue) should report only Lean
foundations (+ `native_decide` axioms).  The analytic frontiers are now named
`Prop` obligations rather than declarations using `sorry`; future submissions
should replace their abstract data by definitions and then prove these `Prop`s. -/

#print axioms lrc14_no_multiple_of_14
#print axioms lrc14_all_odd
#print axioms lrc14_counterexample_saturated
#print axioms gK8_per_shape_bound
#print axioms gK8_dual_feasible
#print axioms gK8_singlefar_arithmetic_kernel
#print axioms gK8_singlefar_arithmetic_checksum
#print axioms single_far_periodmax_headroom
#print axioms genuine_wide_rows_below_cap
#print axioms doublet_witness_floor_checksum
#print axioms sampleBoundedRows_ok
#print axioms wide_bound_from_gK8
#print axioms witnessMP_pos_rat
#print axioms witnessMP_pos_real
#print axioms witness_floor_positive
#print axioms witnessMP_le_rhoGlobFloorRat_real
#print axioms one_seventh_le_inv_of_pos_le_seven
#print axioms witness_floor_from_cluster_cases
#print axioms witness_large_floor_from_rhoGlob_lower_bound
#print axioms lrc14_from_witness_floor_given_nodes
#print axioms lrc14_from_witness_floor_cases_given_nodes
#print axioms lrc14_from_obsolete_rhoStar_floor
#print axioms lrc14_from_witness_floor
-- Open obligations, checked as propositions rather than asserted theorems:
#print axioms thm527_partA_density_pos_implies_reach
#print axioms lonely_of_Mreach_ge
#print axioms doublet_Rtail_uniform_bound
#print axioms gK8_concentration_extremality

end LRC14
end LonelyRunner
