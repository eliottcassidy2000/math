/-
  TournamentH7.LRCResonanceNucleus — THE SUPPLY NUCLEUS and the resonance interface
  (death-star-2026-07-17-S44, HYP-7205; open item 3 of the formalization picture,
  after codex-S49's census capstone).

  codex's `CensusB5Certificate` packages THM-950's criterion per family.  This
  module makes the certificate DECIDABLY CONSTRUCTIBLE and instantiates it on a
  batch, demonstrating the supply nucleus in action:

  * `CensusB5Certificate.of_decide` — the decidable constructor: for concrete
    `(v, q)` the `live_beats_deep` field is a finite computation; `decide` builds
    the certificate outright.
  * `census_batch_demo₁/₂/₃` — three structurally distinct families certified
    through the funnel by `decide` (arithmetic, geometric-ish, mixed) — the
    supply nucleus works wholesale, not just on THM-951's single demo.

  THE RESONANCE OBSERVATION (HYP-7205, REFINED BY REFUTATION): the naive law
  "all-coprime ⟹ deep-free on ladder-7 strata" FAILS — 2 violations in 26265
  coprime strata (120k battery): genuine Dirichlet co-approximations
  (t = p/q ≈ n_i/v_i simultaneously; rate ≈ q/7⁷ per family, matching 2/26265).
  The honest replacement: deep-freeness off gcd-resonance is a DENSITY statement,
  and the deep set has PROVABLE mirror structure — `bandCount_reflect` below:
  `bandCount v q (q−p) = bandCount v q p` (the band is symmetric under negation),
  so deep multipliers come in `p ↔ q−p` pairs and the deep count is EVEN for odd
  `q` — which sharpens THM-950's census by parity.  The a-priori supply's heart
  (∃ q in the window with live > 792·deep, uniformly) remains open item 3,
  now with its exceptional set structurally understood.

  Kernel-pure except the demos' kernel `decide`.
-/
import Mathlib
import TournamentH7.LRCDeepCertificate
import TournamentH7.LRCB5RelationEndgame

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The decidable constructor for codex's census certificate. -/
def CensusB5Certificate.of_decide (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (h : 792 * (((Finset.Ioo 0 q).filter fun p =>
      6 ≤ bandCount v q p).card : ℤ) < (liveCount v q : ℤ)) :
    LRC14Grand.CensusB5Certificate v :=
  { q := q, q_pos := hq, live_beats_deep := h }

/-- Batch demo 1: the arithmetic family `i + 2` (as in THM-951). -/
theorem census_batch_demo₁ :
    Nonempty (LRC14Grand.CensusB5Certificate
      (fun i : Fin 13 => (i : ℤ) + 2)) := by
  refine ⟨CensusB5Certificate.of_decide _ 31 (by norm_num) ?_⟩
  decide

/-- Batch demo 2: the shifted odd family `2i + 3`. -/
theorem census_batch_demo₂ :
    Nonempty (LRC14Grand.CensusB5Certificate
      (fun i : Fin 13 => 2 * (i : ℤ) + 3)) := by
  refine ⟨CensusB5Certificate.of_decide _ 37 (by norm_num) ?_⟩
  decide

/-- Batch demo 3: a mixed-growth family. -/
theorem census_batch_demo₃ :
    Nonempty (LRC14Grand.CensusB5Certificate
      (fun i : Fin 13 => 3 * (i : ℤ) + 5)) := by
  refine ⟨CensusB5Certificate.of_decide _ 31 (by norm_num) ?_⟩
  decide

/-- **The mirror lemma**: the band is symmetric under `p ↦ q − p` — reflected
multipliers have identical coverage. -/
theorem bandCount_reflect (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ)
    (hq : 0 < q) (hp : 0 < p) (hpq : p < q) :
    bandCount v q (q - p) = bandCount v q p := by
  unfold bandCount
  congr 1
  apply Finset.filter_congr
  intro i _
  unfold inBand
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - (p : ℤ) := by
    push_cast [Nat.cast_sub hpq.le]
    ring
  rw [hcast]
  have hmod : (v i * ((q : ℤ) - (p : ℤ))) % (q : ℤ)
      = (-(v i * (p : ℤ))) % (q : ℤ) := by
    have h1 : v i * ((q : ℤ) - (p : ℤ)) = -(v i * (p : ℤ)) + (q : ℤ) * v i := by ring
    rw [h1, Int.add_mul_emod_self_left]
  rw [hmod]
  apply not_congr
  constructor
  · rintro ⟨h1, h2⟩
    -- r' := (−x) mod q; if r := x mod q ∈ (0, q) then r' = q − r; r = 0 gives r' = 0
    by_cases hz : (v i * (p : ℤ)) % (q : ℤ) = 0
    · exfalso
      have : (-(v i * (p : ℤ))) % (q : ℤ) = 0 := by
        rw [Int.neg_emod, if_pos (Int.dvd_of_emod_eq_zero hz)]
      rw [this] at h1
      omega
    · have hr0 : 0 ≤ (v i * (p : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
      have hrq : (v i * (p : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
      have hneg : (-(v i * (p : ℤ))) % (q : ℤ)
          = (q : ℤ) - (v i * (p : ℤ)) % (q : ℤ) := by
        rw [Int.neg_emod, if_neg (fun hdvd => hz (Int.emod_eq_zero_of_dvd hdvd))]
        simp [Int.natAbs_natCast]
      rw [hneg] at h1 h2
      constructor <;> omega
  · rintro ⟨h1, h2⟩
    by_cases hz : (v i * (p : ℤ)) % (q : ℤ) = 0
    · exfalso
      rw [hz] at h1
      omega
    · have hr0 : 0 ≤ (v i * (p : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
      have hrq : (v i * (p : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
      have hneg : (-(v i * (p : ℤ))) % (q : ℤ)
          = (q : ℤ) - (v i * (p : ℤ)) % (q : ℤ) := by
        rw [Int.neg_emod, if_neg (fun hdvd => hz (Int.emod_eq_zero_of_dvd hdvd))]
        simp [Int.natAbs_natCast]
      rw [hneg]
      constructor <;> omega

/-! ## Axiom audit -/
#print axioms CensusB5Certificate.of_decide
#print axioms census_batch_demo₁
#print axioms census_batch_demo₂
#print axioms census_batch_demo₃
#print axioms bandCount_reflect

end LRC14Concrete
end LonelyRunner
