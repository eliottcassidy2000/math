/-
  TournamentH7.LRCTailDiameter — THE TAIL-DIAMETER THEOREM (monad-explorer-S2, HYP-4817;
  monotonicity trick from kind-pasteur-S59, HYP-4797; formalized mac-mini-2026-07-07-S42).

  For a 13-family (k=13 / P=∅ shape) E ⊆ ℤ, the "good set" at threshold θ is
      Good θ E = {x : the circular point set {frac(e·x) : e ∈ E} has max-gap > θ},
  and the k=13 leg of the witness-floor node consumes μ_{1/7}(E) = vol(Good (1/7) E ∩ [0,1])
  ≥ m_P = 14249/252252 (THM-530 / `hlarge` in LRCFourteenSkeleton).

  THE THEOREM (three elementary ingredients):
   1. SUBSET ANTITONICITY.  E ⊆ F ⟹ Good θ F ⊆ Good θ E: an arc empty of the F-orbit is
      empty of the E-orbit.  (Removing points merges gaps; max-gap can only grow.)
   2. TRANSLATION INVARIANCE.  Good θ (E − m) = Good θ E: the (E−m)-orbit at slow phase x
      is the rigid rotation of the E-orbit by m·x, and "∃ an empty θ-arc" is
      rotation-invariant.  (Here per-x, with the explicit witness shift a ↦ a + m·x.)
   3. THE LEDGER LEAF.  μ_{1/7}(AP₇₆) = vol(Good (1/7) {0,…,75} ∩ [0,1])
      = 2314528732/40290957525 ≈ 0.0574452 ≥ m_P (exact Farey-cell computation,
      `lrc14_tail_diameter_theorem_monad_S2.py`; carried here as the certificate Prop
      `AP76Certificate`), and the rational comparison ≥ 14249/252252 is checked by norm_num.

  ASSEMBLY.  diam(E) = D ≤ 75 ⟹ E − min(E) ⊆ {0,…,D} ⊆ {0,…,75} ⟹ (by 1+2)
      Good (1/7) {0,…,75} ⊆ Good (1/7) E   ⟹   μ_{1/7}(E) ≥ μ_{1/7}(AP₇₆) ≥ m_P.
  So EVERY 13-family of diameter ≤ 75 clears the k=13 `hlarge` bar — the entire
  small-diameter regime (all known extremal families: records diam 20/22, min-E[U] diam 43,
  consecutive blocks, GW, parity/cascade shapes).  The open k=13 residual is diam ≥ 76,
  where adversarial floors are ~10× the bar (monad-S2 probe).

  DESIGN NOTES (for auditors).
  · `Good` is defined by the ∃-empty-arc form with STRICT inequalities:
       ∃ a, ∀ e ∈ E, θ < Int.fract (e·x − a),
    i.e. all orbit points avoid the CLOSED arc [a, a+θ] (mod 1).  For a finite orbit this
    is pointwise EQUIVALENT to max-gap > θ (center a closed θ-arc in a gap of length > θ;
    conversely an avoided closed θ-arc sits inside a gap of length > θ).  The equivalence
    is definitional repackaging, not used by the proofs below.
  · Monotonicity/translation are proved as SET inclusions/equalities; measure enters only
    through `measure_mono` (valid for arbitrary sets — no measurability side conditions).
  · The one non-elementary leaf (the exact value of μ_{1/7}(AP₇₆)) is an explicit
    certificate `Prop`, exactly like the skeleton's `rhoGlobFloorRat` ledger: the Lean
    theorem is conditional on it, and it is discharged by the exact rational engine
    outside Lean (Farey-cell affine decomposition, all cells asserted).
  · The bridge from `muGood` to the skeleton's opaque `witnessG2` cannot be a theorem
    here (`witnessG2` is opaque); see the BRIDGE NOTE at the end of the file.
-/
import Mathlib

namespace LonelyRunner
namespace TailDiameter

open MeasureTheory
open scoped ENNReal

/-- The closed θ-arc `[a, a+θ]` (mod 1) at slow phase `x` is empty of the `E`-orbit:
every orbit point `frac(e·x)` keeps fractional distance `> θ` past the arc's left end. -/
def EmptyArc (θ : ℝ) (E : Finset ℤ) (x a : ℝ) : Prop :=
  ∀ e ∈ E, θ < Int.fract ((e : ℝ) * x - a)

/-- The good set at threshold `θ`: some closed θ-arc is empty of the orbit.  For a finite
orbit this says exactly "circular max-gap > θ". -/
def Good (θ : ℝ) (E : Finset ℤ) : Set ℝ :=
  {x | ∃ a : ℝ, EmptyArc θ E x a}

/-- The (one-period) good-set measure `μ_θ(E) = vol(Good θ E ∩ [0,1])`. -/
noncomputable def muGood (θ : ℝ) (E : Finset ℤ) : ℝ≥0∞ :=
  volume (Good θ E ∩ Set.Icc (0 : ℝ) 1)

/-! ### 1. Subset antitonicity -/

/-- An arc empty of a superset's orbit is empty of the subset's orbit. -/
theorem emptyArc_anti {θ : ℝ} {E F : Finset ℤ} (hEF : E ⊆ F) {x a : ℝ}
    (h : EmptyArc θ F x a) : EmptyArc θ E x a :=
  fun e he => h e (hEF he)

/-- `E ⊆ F ⟹ Good θ F ⊆ Good θ E` : fewer points, bigger gaps. -/
theorem good_anti {θ : ℝ} {E F : Finset ℤ} (hEF : E ⊆ F) :
    Good θ F ⊆ Good θ E := by
  rintro x ⟨a, ha⟩
  exact ⟨a, emptyArc_anti hEF ha⟩

/-- Measure form of subset antitonicity. -/
theorem muGood_anti (θ : ℝ) {E F : Finset ℤ} (hEF : E ⊆ F) :
    muGood θ F ≤ muGood θ E :=
  measure_mono (Set.inter_subset_inter_left _ (good_anti hEF))

/-! ### 2. Translation invariance -/

/-- Orbit translation: the `(e−m)`-phase at witness `a` equals the `e`-phase at the
rotated witness `a + m·x` — the same real number inside `Int.fract`. -/
theorem emptyArc_translate {θ : ℝ} (E : Finset ℤ) (m : ℤ) (x a : ℝ) :
    EmptyArc θ (E.image (fun e => e - m)) x a ↔ EmptyArc θ E x (a + (m : ℝ) * x) := by
  constructor
  · intro h e he
    have hmem : e - m ∈ E.image (fun e => e - m) := Finset.mem_image_of_mem _ he
    have := h (e - m) hmem
    have harg : ((e - m : ℤ) : ℝ) * x - a = (e : ℝ) * x - (a + (m : ℝ) * x) := by
      push_cast; ring
    rwa [harg] at this
  · intro h f hf
    rcases Finset.mem_image.mp hf with ⟨e, he, rfl⟩
    have := h e he
    have harg : ((e - m : ℤ) : ℝ) * x - a = (e : ℝ) * x - (a + (m : ℝ) * x) := by
      push_cast; ring
    rwa [harg]

/-- `Good θ (E − m) = Good θ E` : existence of an empty arc is rotation-invariant. -/
theorem good_translate (θ : ℝ) (E : Finset ℤ) (m : ℤ) :
    Good θ (E.image (fun e => e - m)) = Good θ E := by
  ext x
  constructor
  · rintro ⟨a, ha⟩
    exact ⟨a + (m : ℝ) * x, (emptyArc_translate E m x a).mp ha⟩
  · rintro ⟨b, hb⟩
    refine ⟨b - (m : ℝ) * x, (emptyArc_translate E m x (b - (m : ℝ) * x)).mpr ?_⟩
    simpa using hb

/-- Measure form of translation invariance. -/
theorem muGood_translate (θ : ℝ) (E : Finset ℤ) (m : ℤ) :
    muGood θ (E.image (fun e => e - m)) = muGood θ E := by
  unfold muGood
  rw [good_translate]

/-! ### 3. The diameter chain -/

/-- **The tail-diameter chain.**  If `E` (translated by `m`) fits inside `{0,…,D}` with
`0 ≤ D ≤ 75`, then the AP₇₆ good set is contained (after translation) in `E`'s good set,
so `μ_θ(E) ≥ μ_θ(AP₇₆)` — for EVERY threshold `θ`. -/
theorem muGood_ge_AP76 (θ : ℝ) (E : Finset ℤ) (m D : ℤ) (hD : D ≤ 75)
    (hE : ∀ e ∈ E, e - m ∈ Finset.Icc (0 : ℤ) D) :
    muGood θ (Finset.Icc (0 : ℤ) 75) ≤ muGood θ E := by
  have hsub1 : E.image (fun e => e - m) ⊆ Finset.Icc (0 : ℤ) D := by
    intro f hf
    rcases Finset.mem_image.mp hf with ⟨e, he, rfl⟩
    exact hE e he
  have hsub2 : Finset.Icc (0 : ℤ) D ⊆ Finset.Icc (0 : ℤ) 75 := by
    intro j hj
    rcases Finset.mem_Icc.mp hj with ⟨h0, hD'⟩
    exact Finset.mem_Icc.mpr ⟨h0, le_trans hD' hD⟩
  calc muGood θ (Finset.Icc (0 : ℤ) 75)
      ≤ muGood θ (Finset.Icc (0 : ℤ) D) := muGood_anti θ hsub2
    _ ≤ muGood θ (E.image (fun e => e - m)) := muGood_anti θ hsub1
    _ = muGood θ E := muGood_translate θ E m

/-! ### 4. The ledger leaf (certificate) and the assembly -/

/-- The exact-engine certificate: `μ_{1/7}(AP₇₆) = 2314528732/40290957525` (Farey-cell
exact computation, `lrc14_tail_diameter_theorem_monad_S2.py`; this Prop asserts only the
`≥` direction, which is what the assembly consumes). -/
def AP76Certificate : Prop :=
  ENNReal.ofReal ((2314528732 : ℝ) / 40290957525) ≤
    muGood (1 / 7) (Finset.Icc (0 : ℤ) 75)

/-- The certified value clears the THM-530 bar `m_P = 14249/252252`: pure arithmetic. -/
theorem cert_value_ge_mP :
    (14249 : ℝ) / 252252 ≤ (2314528732 : ℝ) / 40290957525 := by
  norm_num

/-- **THE TAIL-DIAMETER THEOREM (assembly).**  Conditional only on the AP₇₆ ledger
certificate: every integer family that fits in a translated window `{m, …, m+D}` with
`D ≤ 75` has `μ_{1/7} ≥ m_P` — the k=13 `hlarge` bar. -/
theorem hlarge_floor_of_diam_le (hcert : AP76Certificate)
    (E : Finset ℤ) (m D : ℤ) (hD : D ≤ 75)
    (hE : ∀ e ∈ E, e - m ∈ Finset.Icc (0 : ℤ) D) :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) E := by
  refine le_trans ?_ (le_trans hcert (muGood_ge_AP76 (1 / 7) E m D hD hE))
  exact ENNReal.ofReal_le_ofReal cert_value_ge_mP

/- BRIDGE NOTE (documentation, not a def): on pure-cluster (k=13, P=∅) shapes the
skeleton's `witnessG2` is intended to BE this file's `muGood (1/7)`; since `witnessG2`
is an opaque constant, that identification cannot be a theorem here — it is discharged
by replacing the opaque constant with a definition when the skeleton is next rewired
(monad-S1's "sieve BEFORE floor" rewiring is the natural occasion). -/

end TailDiameter
end LonelyRunner
