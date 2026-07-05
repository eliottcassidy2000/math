/-
  TournamentH7.LRCHcompSurface — THE hcomp RESIDUAL, PINNED AS TWO NAMED PREDICATES
  (klein-2026-07-05-S133, HYP-4096).

  Assembles the S132 trio (one-window peel, primitivity split, CRT free-rider) into the
  full reduction:  hcomp  ⟸  TightLooseDichotomy + CornerLonely,  hence

      LRC(14)  ≤  LRC(13) citation  +  TightLooseDichotomy  +  CornerLonely.

  * `tight_free_rider'` GENERALIZES the S132 free-rider: the base values need only be a
    SUBSET of `c·{1,…,12}` up to sign (each `|v i| = c·j`, no exhaustion, no ordering) —
    the per-runner residue argument never used more, and the sign is absorbed by the
    `∀ m` quantifier.  This removes ALL permutation/multiset glue from the tight side.
  * `gcd_killer_of_primitive`: `gcd(c, killer) = 1` is DERIVED from family primitivity.
  * `TightLooseDichotomy` — the n=12 rigidity dichotomy at the argmax peel (the open
    math: tight base values sit inside `c·{1,…,12}`, or the base has margin `2/25`).
  * `CornerLonely` — the sub-threshold corner (killer `≤ (25/3)·B`): the alignment-band
    lane (mac-mini THM-619 / S49 sweep).
  Kernel-pure.
-/
import TournamentH7.LRCTightAPFreeRider
import TournamentH7.LRCOneWindowPeel
import TournamentH7.LRCEndgameAssembly

namespace LonelyRunner
namespace HcompSurface

open LRC14 TightAPFreeRider OneWindowPeel

/-- **The generalized free-rider**: base values ANY subset of `c·{1,…,12}` up to sign,
killer attached with `gcd(c, v*) = 1`.  No multiset exhaustion, no ordering. -/
theorem tight_free_rider' (v : Fin 13 → ℤ) (istar : Fin 13) (c : ℤ) (hc : 2 ≤ c)
    (hbase : ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j)
    (hgcd : Int.gcd c (v istar) = 1) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨k, r, qq, hk13, hEqK, hrlo, hrhi⟩ := killer_target c (v istar) hc hgcd
  refine ⟨_, lonely14_of_ratio v k (13 * c) (by omega) ?_⟩
  intro i m
  by_cases hi : i = istar
  · subst hi
    exact residue_key k (13 * c) c m (v i) qq r (by omega) (by omega) hEqK hrlo (by omega)
  · obtain ⟨j, hj1, hj12, habs⟩ := hbase i hi
    have hcj : (0 : ℤ) ≤ c * j := by positivity
    rcases abs_eq hcj |>.mp habs with hpos | hneg
    · -- v i = c·j: the S132 base argument verbatim
      rw [hpos]
      have hp13 : Prime (13 : ℤ) := by
        rw [Int.prime_iff_natAbs_prime]; norm_num
      have hjk : ¬ ((13 : ℤ) ∣ j * k) := by
        intro hdvd
        rcases hp13.dvd_mul.mp hdvd with hdj | hdk
        · have := Int.le_of_dvd (by omega) hdj
          omega
        · exact hk13 hdk
      have hmod0 : (j * k) % 13 ≠ 0 := fun h => hjk (Int.dvd_of_emod_eq_zero h)
      have hmodnn : 0 ≤ (j * k) % 13 := Int.emod_nonneg _ (by norm_num)
      have hmodlt : (j * k) % 13 < 13 := Int.emod_lt_of_pos _ (by norm_num)
      have hde : 13 * ((j * k) / 13) + (j * k) % 13 = j * k := Int.ediv_add_emod _ _
      have hs1 : 1 ≤ (j * k) % 13 := by omega
      have hs12 : (j * k) % 13 ≤ 12 := by omega
      apply residue_key k (13 * c) c m (c * j) ((j * k) / 13) (c * ((j * k) % 13))
        (by omega) (by omega)
      · linear_combination (-c) * hde
      · nlinarith [hs1, hc]
      · nlinarith [hs12, hc]
    · -- v i = −c·j: absorb the sign into the m-quantifier
      rw [hneg]
      have hkey : ∀ m' : ℤ, 13 * c ≤ 14 * |(c * j) * k - m' * (13 * c)| := by
        intro m'
        have hp13 : Prime (13 : ℤ) := by
          rw [Int.prime_iff_natAbs_prime]; norm_num
        have hjk : ¬ ((13 : ℤ) ∣ j * k) := by
          intro hdvd
          rcases hp13.dvd_mul.mp hdvd with hdj | hdk
          · have := Int.le_of_dvd (by omega) hdj
            omega
          · exact hk13 hdk
        have hmod0 : (j * k) % 13 ≠ 0 := fun h => hjk (Int.dvd_of_emod_eq_zero h)
        have hmodnn : 0 ≤ (j * k) % 13 := Int.emod_nonneg _ (by norm_num)
        have hmodlt : (j * k) % 13 < 13 := Int.emod_lt_of_pos _ (by norm_num)
        have hde : 13 * ((j * k) / 13) + (j * k) % 13 = j * k := Int.ediv_add_emod _ _
        have hs1 : 1 ≤ (j * k) % 13 := by omega
        have hs12 : (j * k) % 13 ≤ 12 := by omega
        apply residue_key k (13 * c) c m' (c * j) ((j * k) / 13) (c * ((j * k) % 13))
          (by omega) (by omega)
        · linear_combination (-c) * hde
        · nlinarith [hs1, hc]
        · nlinarith [hs12, hc]
      have := hkey (-m)
      have hflip : |(c * j) * k - (-m) * (13 * c)| = |(-(c * j)) * k - m * (13 * c)| := by
        rw [show (-(c * j)) * k - m * (13 * c) = -((c * j) * k - (-m) * (13 * c)) by ring,
          abs_neg]
      rw [hflip] at this
      exact this

/-- `gcd(c, killer) = 1` from family primitivity: `c` divides every base entry, so any
common divisor of `c` and the killer divides the whole family's gcd. -/
theorem gcd_killer_of_primitive (v : Fin 13 → ℤ) (istar : Fin 13) (c : ℤ) (hc : 2 ≤ c)
    (hbase : ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j)
    (hprim : tupleGcd v = 1) :
    Int.gcd c (v istar) = 1 := by
  set g : ℕ := Int.gcd c (v istar) with hg
  have hgc : (g : ℤ) ∣ c := by
    rw [hg]; exact Int.gcd_dvd_left (a := c) (b := v istar)
  have hgk : (g : ℤ) ∣ v istar := by
    rw [hg]; exact Int.gcd_dvd_right (a := c) (b := v istar)
  have hall : ∀ i, g ∣ (v i).natAbs := by
    intro i
    by_cases hi : i = istar
    · subst hi
      exact Int.natCast_dvd_natCast.mp (Int.dvd_natAbs.mpr hgk)
    · obtain ⟨j, -, -, habs⟩ := hbase i hi
      have hcv : (g : ℤ) ∣ |v i| := by
        rw [habs]
        exact Dvd.dvd.mul_right hgc j
      rw [Int.abs_eq_natAbs] at hcv
      exact Int.natCast_dvd_natCast.mp hcv
  have hdvd : g ∣ tupleGcd v := by
    unfold tupleGcd
    exact Finset.dvd_gcd (fun i _ => hall i)
  rw [hprim] at hdvd
  exact Nat.dvd_one.mp hdvd

/-! ## The two named residual predicates -/

/-- **THE n=12 RIGIDITY DICHOTOMY at the argmax peel** (the open mathematics of the
tight side): for every primitive compressed covering family, the base at the argmax
peel either has all its values inside `c·{1,…,12}` (some dilation `c ≥ 2`) up to sign,
or carries margin `2/25` at some point.  Empirically TRUE with the 12-runner spectrum
`1/13 (AP) < 2/25 ({1..11,24}·c) < …` (mac-mini S46/S47); the classification is
MISTAKE-100-class open math until proved or citation-confirmed. -/
def TightLooseDichotomy : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
    ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
      (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2/25 ≤ |(v i : ℝ) * tstar - m|)

/-- **THE SUB-THRESHOLD CORNER** (the alignment-band lane, mac-mini THM-619/S49):
a loose base (margin `2/25`) whose killer sits BELOW the one-window threshold
`(25/3)·B`.  Everything above the threshold is closed by `lonely_of_window_margin`. -/
def CornerLonely : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
    ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
    ∀ B : ℤ, (∀ i, i ≠ istar → |v i| ≤ B) → (∃ i, i ≠ istar ∧ |v i| = B) →
      (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2/25 ≤ |(v i : ℝ) * tstar - m|) →
      3 * |v istar| ≤ 25 * B →
      ∃ t : ℝ, Lonely 14 v t

/-! ## The assembly -/

/-- **hcomp FROM THE TWO NAMED PREDICATES**: primitivity split → argmax peel →
dichotomy dispatch → {tight: generalized free-rider} / {loose: one-window lemma above
the threshold, corner hypothesis below}. -/
theorem hcomp_of_dichotomy_and_corner
    (hdich : TightLooseDichotomy) (hcorner : CornerLonely) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t := by
  apply hcomp_of_primitive
  intro v hv hcov hcomp hprim
  obtain ⟨istar, -, hmaxmem⟩ :=
    Finset.exists_max_image (Finset.univ : Finset (Fin 13)) (fun i => |v i|)
      ⟨0, Finset.mem_univ 0⟩
  have hmax : ∀ i, |v i| ≤ |v istar| := fun i => hmaxmem i (Finset.mem_univ i)
  rcases hdich v hv hcov hcomp hprim istar hmax with ⟨c, hc, hbase⟩ | ⟨tstar, hloose⟩
  · exact tight_free_rider' v istar c hc hbase
      (gcd_killer_of_primitive v istar c hc hbase hprim)
  · -- loose: build the base max B
    have hne : (Finset.univ.erase istar).Nonempty := by
      rw [← Finset.card_pos, Finset.card_erase_of_mem (Finset.mem_univ istar)]
      simp
    obtain ⟨ibase, hibmem, hibmax⟩ :=
      Finset.exists_max_image (Finset.univ.erase istar) (fun i => |v i|) hne
    have hibne : ibase ≠ istar := Finset.ne_of_mem_erase hibmem
    set B : ℤ := |v ibase| with hB
    have hBbound : ∀ i, i ≠ istar → |v i| ≤ B := by
      intro i hi
      exact hibmax i (Finset.mem_erase.mpr ⟨hi, Finset.mem_univ i⟩)
    have hBpos : 0 < B := by
      rw [hB]
      exact abs_pos.mpr (hv ibase)
    by_cases hthr : 25 * B < 3 * |v istar|
    · -- above threshold: the one-window lemma at β = 2/25
      refine lonely_of_window_margin v istar tstar (2/25) ((B : ℤ) : ℝ)
        (by norm_num) hloose ?_ (by exact_mod_cast hBpos) (hv istar) ?_
      · intro i hi
        have := hBbound i hi
        calc |(v i : ℝ)| = ((|v i| : ℤ) : ℝ) := by push_cast; rfl
          _ ≤ ((B : ℤ) : ℝ) := by exact_mod_cast this
      · -- (B:ℝ) < 14·(2/25 − 1/14)·|(v istar : ℝ)| = (3/25)·|v istar|
        have h1 : (25 * B : ℤ) < 3 * |v istar| := hthr
        have h2 : ((25 * B : ℤ) : ℝ) < ((3 * |v istar| : ℤ) : ℝ) := by exact_mod_cast h1
        push_cast at h2
        have habs : |(v istar : ℝ)| = ((|v istar| : ℤ) : ℝ) := by push_cast; rfl
        rw [show (14 : ℝ) * (2/25 - 1/14) * |(v istar : ℝ)|
            = (3/25) * |(v istar : ℝ)| by ring, habs]
        push_cast
        linarith
    · exact hcorner v hv hcov hcomp hprim istar hmax B hBbound
        ⟨ibase, hibne, rfl⟩ ⟨tstar, hloose⟩ (by omega)

/-- **THE PINNED SURFACE**: LRC(14) from the citation node and the two named residual
predicates.  After S132+S133, the ENTIRE remaining mathematics of LRC(14) is
`TightLooseDichotomy` (n=12 rigidity) and `CornerLonely` (sub-threshold corner). -/
theorem lrc14_of_dichotomy_and_corner (cite : LonelyRunner.LRCUpTo13)
    (hdich : TightLooseDichotomy) (hcorner : CornerLonely) :
    LRC14Statement :=
  lrc14_of_compressed cite (hcomp_of_dichotomy_and_corner hdich hcorner)

/-! ## The β-parametric surface

The margin 2/25 above is one instantiation; the assembly works for ANY loose margin
β > 1/14, trading dichotomy strength against corner size (threshold `B/(14(β−1/14))`).
The parametric form lets the rigidity lane land at whatever margin it can prove. -/

/-- The dichotomy at loose margin `β`. -/
def TightLooseDichotomyAt (β : ℝ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
    ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
      (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|)

/-- The corner at loose margin `β`: killer below the one-window threshold. -/
def CornerLonelyAt (β : ℝ) : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
    ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
    ∀ B : ℤ, (∀ i, i ≠ istar → |v i| ≤ B) → (∃ i, i ≠ istar ∧ |v i| = B) →
      (∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * tstar - m|) →
      14 * (β - 1/14) * |(v istar : ℝ)| ≤ (B : ℝ) →
      ∃ t : ℝ, Lonely 14 v t

/-- **The β-parametric assembly**: `hcomp` from the dichotomy and corner at ANY
loose margin `β > 1/14`. -/
theorem hcomp_of_dichotomy_and_corner_at (β : ℝ) (hβ : 1/14 < β)
    (hdich : TightLooseDichotomyAt β) (hcorner : CornerLonelyAt β) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t := by
  apply hcomp_of_primitive
  intro v hv hcov hcomp hprim
  obtain ⟨istar, -, hmaxmem⟩ :=
    Finset.exists_max_image (Finset.univ : Finset (Fin 13)) (fun i => |v i|)
      ⟨0, Finset.mem_univ 0⟩
  have hmax : ∀ i, |v i| ≤ |v istar| := fun i => hmaxmem i (Finset.mem_univ i)
  rcases hdich v hv hcov hcomp hprim istar hmax with ⟨c, hc, hbase⟩ | ⟨tstar, hloose⟩
  · exact tight_free_rider' v istar c hc hbase
      (gcd_killer_of_primitive v istar c hc hbase hprim)
  · have hne : (Finset.univ.erase istar).Nonempty := by
      rw [← Finset.card_pos, Finset.card_erase_of_mem (Finset.mem_univ istar)]
      simp
    obtain ⟨ibase, hibmem, hibmax⟩ :=
      Finset.exists_max_image (Finset.univ.erase istar) (fun i => |v i|) hne
    have hibne : ibase ≠ istar := Finset.ne_of_mem_erase hibmem
    set B : ℤ := |v ibase| with hB
    have hBbound : ∀ i, i ≠ istar → |v i| ≤ B := by
      intro i hi
      exact hibmax i (Finset.mem_erase.mpr ⟨hi, Finset.mem_univ i⟩)
    have hBpos : 0 < B := by
      rw [hB]
      exact abs_pos.mpr (hv ibase)
    by_cases hthr : (B : ℝ) < 14 * (β - 1/14) * |(v istar : ℝ)|
    · refine lonely_of_window_margin v istar tstar β ((B : ℤ) : ℝ)
        hβ hloose ?_ (by exact_mod_cast hBpos) (hv istar) hthr
      intro i hi
      have := hBbound i hi
      calc |(v i : ℝ)| = ((|v i| : ℤ) : ℝ) := by push_cast; rfl
        _ ≤ ((B : ℤ) : ℝ) := by exact_mod_cast this
    · exact hcorner v hv hcov hcomp hprim istar hmax B hBbound
        ⟨ibase, hibne, rfl⟩ ⟨tstar, hloose⟩ (le_of_not_gt hthr)

/-- **The β-parametric pinned surface**. -/
theorem lrc14_of_dichotomy_and_corner_at (cite : LonelyRunner.LRCUpTo13)
    (β : ℝ) (hβ : 1/14 < β)
    (hdich : TightLooseDichotomyAt β) (hcorner : CornerLonelyAt β) :
    LRC14Statement :=
  lrc14_of_compressed cite (hcomp_of_dichotomy_and_corner_at β hβ hdich hcorner)

#print axioms tight_free_rider'
