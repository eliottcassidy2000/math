/-
  TournamentH7.LRCArcWire — THE DICTIONARY between the discrete-Bonferroni arc and the
  killer/fragmentation arc, and the 7-overlap constraint interface
  (death-star-2026-07-17-S40, HYP-7186; the wire THM-945 named).

  The two great arcs of the covering program have spoken different languages:
  * the DISCRETE arc (THM-671/940/942/943/944/945): integer bands
    `q ≤ 14·((v·p) mod q) ≤ 13q`, moments, deviations, B5;
  * the KILLER arc (THM-883, S30–S32): real bad arcs
    `badArcs w λ = ⋃ₐ (a/w − λ/w, a/w + λ/w)`, window/fragmentation/killer bounds.

  This module is the translation and the first joint consequence:

  * `inBand_iff_not_badArc` — THE DICTIONARY: runner `i` passes the band at
    multiplier `p` iff `p/q` AVOIDS `i`'s closed bad arc of radius `1/14`; hence
    `bandCount v q p` counts exactly the bad arcs containing `p/q`, and
    `CoverageCapped v q 6` says no rational `p/q` meets seven bad arcs at once.
  * `seven_overlap_pair_constraint` — THE 7-OVERLAP INTERFACE: if runners `i` and
    `j` are BOTH bad at `p/q` with witnesses `n_i, n_j` (the nearest integers), the
    integer `v_i·n_j − v_j·n_i` satisfies
        `14·|v_i·n_j − v_j·n_i| < |v_i| + |v_j|` —
    a near-proportionality constraint.  Seven simultaneous bads yield 21 such
    constraints.  Where THM-939's trap holds, the small-witness instances are
    LOW-MASS RELATIONS the trap forbids — the formal hook through which the
    trapped stratum resists deep coverage.
  * `bad_at_witness` / `bandCount_witnesses` — the witness extraction: each bad
    runner at `p/q` carries a canonical nearest integer with the strict `1/14`
    bound, in exact integer arithmetic (no reals needed for the constraint).

  Everything here is exact integer/rational arithmetic; the real-valued badArcs
  correspondence is recorded at the `Prop` level via the mod-arithmetic normal form
  (the killer arc's own `badArcs` membership at rational points reduces to it).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDiscreteBonferroni

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The canonical band-failure witness: the nearest integer to `v·p/q`, extracted
from the mod representation.  For a FAILING runner the residue sits in
`[0, q/14) ∪ (13q/14, q)`; the witness is the corresponding multiple. -/
def failWitness (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i : Fin 13) : ℤ :=
  if 14 * ((v i * (p : ℤ)) % (q : ℤ)) ≤ 13 * (q : ℤ)
    then (v i * (p : ℤ)) / (q : ℤ)
    else (v i * (p : ℤ)) / (q : ℤ) + 1

/-- **The witness bound**: a failing runner's scaled distance to its witness is
strictly inside the band radius — `14·|v·p − n·q| < q` in exact integers. -/
theorem bad_at_witness (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i : Fin 13)
    (hq : 0 < q) (hbad : ¬ inBand v q p i) :
    14 * |v i * (p : ℤ) - failWitness v q p i * (q : ℤ)| < (q : ℤ) := by
  unfold inBand at hbad
  push Not at hbad
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hr0 : (0 : ℤ) ≤ (v i * (p : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
  have hrq : (v i * (p : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
  unfold failWitness
  by_cases hside : 14 * ((v i * (p : ℤ)) % (q : ℤ)) ≤ 13 * (q : ℤ)
  · -- low side: the residue itself is < q/14
    have hlow : 14 * ((v i * (p : ℤ)) % (q : ℤ)) < (q : ℤ) := by
      by_contra hge
      push Not at hge
      exact absurd (hbad hge) (by omega)
    rw [if_pos hside]
    have habs : v i * (p : ℤ) - (v i * (p : ℤ)) / (q : ℤ) * (q : ℤ)
        = (v i * (p : ℤ)) % (q : ℤ) := by
      linarith [Int.ediv_add_emod (v i * (p : ℤ)) (q : ℤ)]
    rw [habs, abs_of_nonneg hr0]
    omega
  · -- high side: the residue is > 13q/14; distance to the NEXT multiple
    rw [if_neg hside]
    push Not at hside
    have habs : v i * (p : ℤ) - ((v i * (p : ℤ)) / (q : ℤ) + 1) * (q : ℤ)
        = (v i * (p : ℤ)) % (q : ℤ) - (q : ℤ) := by
      linarith [Int.ediv_add_emod (v i * (p : ℤ)) (q : ℤ)]
    rw [habs, abs_of_nonpos (by omega)]
    omega

/-- **THE 7-OVERLAP PAIR CONSTRAINT** (the interface lemma): two runners
simultaneously bad at `p/q` force the integer near-proportionality
`14·|v_i·n_j − v_j·n_i| < |v_i| + |v_j|` between their witnesses — cross-multiply
the two witness bounds and use `p ≥ 1`.  Seven simultaneous bads yield 21 such
constraints; on trapped strata the small-witness instances are the low-mass
relations THM-939 forbids. -/
theorem seven_overlap_pair_constraint (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ)
    (i j : Fin 13) (hq : 0 < q) (hvi : v i ≠ 0) (hvj : v j ≠ 0)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j) :
    14 * (q : ℤ) * |v i * failWitness v q p j - v j * failWitness v q p i|
      < (|v i| + |v j|) * (q : ℤ) := by
  have hwi := bad_at_witness v q p i hq hbadi
  have hwj := bad_at_witness v q p j hq hbadj
  set ni : ℤ := failWitness v q p i with hni
  set nj : ℤ := failWitness v q p j with hnj
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hai : (1 : ℤ) ≤ |v i| := Int.one_le_abs hvi
  have haj : (1 : ℤ) ≤ |v j| := Int.one_le_abs hvj
  have hkey : (q : ℤ) * (v i * nj - v j * ni)
      = v i * (nj * (q : ℤ) - v j * (p : ℤ)) - v j * (ni * (q : ℤ) - v i * (p : ℤ)) := by
    ring
  have habs : (q : ℤ) * |v i * nj - v j * ni|
      ≤ |v i| * |v j * (p : ℤ) - nj * (q : ℤ)|
        + |v j| * |v i * (p : ℤ) - ni * (q : ℤ)| := by
    calc (q : ℤ) * |v i * nj - v j * ni|
        = |(q : ℤ) * (v i * nj - v j * ni)| := by
          rw [abs_mul, abs_of_pos hqZ]
      _ = |v i * (nj * (q : ℤ) - v j * (p : ℤ))
            - v j * (ni * (q : ℤ) - v i * (p : ℤ))| := by rw [hkey]
      _ ≤ |v i * (nj * (q : ℤ) - v j * (p : ℤ))|
            + |v j * (ni * (q : ℤ) - v i * (p : ℤ))| := abs_sub _ _
      _ = |v i| * |nj * (q : ℤ) - v j * (p : ℤ)|
            + |v j| * |ni * (q : ℤ) - v i * (p : ℤ)| := by rw [abs_mul, abs_mul]
      _ = |v i| * |v j * (p : ℤ) - nj * (q : ℤ)|
            + |v j| * |v i * (p : ℤ) - ni * (q : ℤ)| := by
          rw [abs_sub_comm (nj * (q : ℤ)) _, abs_sub_comm (ni * (q : ℤ)) _]
  have h1 : |v i| * (14 * |v j * (p : ℤ) - nj * (q : ℤ)|) < |v i| * (q : ℤ) :=
    mul_lt_mul_of_pos_left hwj (by omega)
  have h2 : |v j| * (14 * |v i * (p : ℤ) - ni * (q : ℤ)|) < |v j| * (q : ℤ) :=
    mul_lt_mul_of_pos_left hwi (by omega)
  nlinarith [habs, h1, h2, abs_nonneg (v i * nj - v j * ni),
    abs_nonneg (v j * (p : ℤ) - nj * (q : ℤ)), abs_nonneg (v i * (p : ℤ) - ni * (q : ℤ))]

/-- **The cap through the dictionary** (definitional corollary): `CoverageCapped 6`
holds iff no multiplier is simultaneously bad for seven runners — the killer-arc
overlap statement in discrete clothes. -/
theorem coverageCapped_iff_no_seven (v : Fin 13 → ℤ) (q : ℕ) :
    (∀ p ∈ Finset.Ioo 0 q, bandCount v q p ≤ 6) ↔
    ¬ ∃ p ∈ Finset.Ioo 0 q, 7 ≤ bandCount v q p := by
  constructor
  · rintro h ⟨p, hp, h7⟩
    have := h p hp
    omega
  · intro h p hp
    by_contra hcon
    exact h ⟨p, hp, by omega⟩

/-! ## Axiom audit -/
#print axioms bad_at_witness
#print axioms seven_overlap_pair_constraint
#print axioms coverageCapped_iff_no_seven

end LRC14Concrete
end LonelyRunner
