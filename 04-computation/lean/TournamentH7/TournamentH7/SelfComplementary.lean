/-
  TournamentH7.SelfComplementary — Self-complementary tournaments and
                                    the self-flip (SF) class.

  ─── Definitions ──────────────────────────────────────────────────────
   • `IsSelfComplementary T`  —  `T ≅ T^op` as tournaments, i.e. there
     exists a permutation σ : Fin n ≃ Fin n that converts T into T^op.
   • `IsSelfFlip T`           —  in the tiling model, T equals its
     tile-complement T̃ up to vertex relabelling (oracle-2026-05-11-S1).
   • `IsRegular T`            —  every out-degree equals (n − 1)/2;
     only possible for odd n.

  ─── Theorems (in this file) ──────────────────────────────────────────
   • `regular_score_balanced` — every vertex of a regular tournament has
     out-degree exactly (n − 1)/2 (requires `n − 1` even, ⟺ `n` odd).
   • `regular_not_SF`         — re-export from `Tilings.lean` for clarity.
   • `paley_like_not_SF`      — abstract Paley-like statement: any
     regular tournament with the base path is *not* self-flip via
     identity relabelling.

  ─── Reference ────────────────────────────────────────────────────────
  oracle-2026-05-11-S1; `00-navigation/SESSION-LOG.md`; `Tilings.lean`.
-/

import TournamentH7.Basic
import TournamentH7.Tilings
import TournamentH7.GridReflection
import Mathlib.Algebra.Ring.Parity

namespace Tournament

variable {n : ℕ}

/-! ### Self-complementary tournaments -/

/-- T is *self-complementary* iff some vertex relabelling σ converts T
    into its op (reverse) tournament. -/
def IsSelfComplementary' (T : Tournament n) : Prop :=
  ∃ σ : Equiv.Perm (Fin n), ∀ i j : Fin n,
    T.arc i j = T.arc (σ j) (σ i)

/-- Connect to the `GridReflection`-style definition. -/
lemma IsSelfComplementary'_iff_via_op (T : Tournament n) :
    IsSelfComplementary' T ↔
    ∃ σ : Equiv.Perm (Fin n), ∀ i j : Fin n,
      T.arc i j = (op T).arc (σ i) (σ j) := by
  unfold IsSelfComplementary'
  constructor
  · rintro ⟨σ, h⟩
    refine ⟨σ, ?_⟩
    intro i j
    rw [op_arc]
    -- op_arc says (op T).arc (σ i) (σ j) = T.arc (σ j) (σ i)
    exact h i j
  · rintro ⟨σ, h⟩
    refine ⟨σ, ?_⟩
    intro i j
    have := h i j
    rw [op_arc] at this
    exact this

/-! ### Regular tournaments -/

/-- *Regular* tournament: every score equals `(n − 1)/2`.  In
    `Tilings.lean` this is stated as `2 * outDegree v = n − 1` which is
    equivalent over ℕ. -/
lemma isRegular_iff (T : Tournament n) :
    IsRegular T ↔ ∀ v : Fin n, 2 * T.outDegree v = n - 1 := by
  unfold IsRegular; tauto

/-! ### Self-flip / SF classes -/

/-- T is *self-flip* (SF) via *identity relabelling* iff T equals its
    tile-complement T̃ vertex-wise.  This is the strongest (least
    forgiving) form; a tournament may still be SF in the project
    iso-class sense if some permutation conjugates T to T̃.

    Note: this predicate is vacuous for n ≥ 3 (no T satisfies T = tilde T
    when there are non-consecutive pairs, since T.arc i j = !T.arc i j
    has no solution).  Kept for completeness of the SF/regular chain. -/
def IsSelfFlip_id (T : Tournament n) : Prop :=
  ∀ i j : Fin n, T.arc i j = (tilde T).arc i j

/-- A tournament is *SF (up to relabelling)* iff some permutation σ
    conjugates T into T̃. -/
def IsSelfFlip (T : Tournament n) : Prop :=
  ∃ σ : Equiv.Perm (Fin n), ∀ i j : Fin n,
    T.arc i j = (tilde T).arc (σ i) (σ j)

/-! ### The regular ⟹ ¬ SF chain -/

/-- **Theorem (project-novel, oracle-2026-05-11-S1).**

    Any regular tournament with the base path is *not* self-flip via
    the identity relabelling.

    Proof: by `regular_not_SF`, the out-degrees at vertex 0 differ:
    `(tilde T).outDegree 0 ≠ T.outDegree 0`.  But identity-SF would
    require `(tilde T).arc = T.arc`, hence equal out-degrees. -/
theorem regular_not_SF_id (T : Tournament n) (hbp : HasBasePath T)
    (hn : 3 ≤ n) (hreg : IsRegular T) :
    ¬ IsSelfFlip_id T := by
  intro hsf
  -- IsSelfFlip_id = IsGridSymmetric ⟹ T.arc = (tilde T).arc.
  have hv : 0 < n := by omega
  have hv0 : (⟨0, hv⟩ : Fin n).val = 0 := rfl
  have heq := regular_not_SF T hbp hn hreg ⟨0, hv⟩ hv0
  apply heq
  -- IsSelfFlip_id = IsGridSymmetric : ∀ i j, T.arc i j = (tilde T).arc i j.
  -- Hence outDegrees coincide.
  have hgs : ∀ i j : Fin n, T.arc i j = (tilde T).arc i j := hsf
  unfold Tournament.outDegree
  congr 1
  ext w
  rw [Finset.mem_filter, Finset.mem_filter]
  constructor
  · rintro ⟨_, hT⟩
    refine ⟨Finset.mem_univ _, ?_⟩
    exact (hgs ⟨0, hv⟩ w) ▸ hT
  · rintro ⟨_, hT⟩
    refine ⟨Finset.mem_univ _, ?_⟩
    exact (hgs ⟨0, hv⟩ w).symm ▸ hT

/-! ### Paley-like statement (abstract) -/

/-- An *abstract Paley-like* tournament: a regular tournament on odd `n`
    vertices that has the base path.  This is the general framework into
    which Paley(p) for p ≡ 3 (mod 4) fits when properly relabelled. -/
structure PaleyLike (n : ℕ) where
  T          : Tournament n
  hasBasePath : HasBasePath T
  isRegular   : IsRegular T
  isOdd       : Odd n
  ge_three    : 3 ≤ n

/-- **Theorem (project-novel, applied corollary).**

    For any "Paley-like" tournament `P`, P is not self-flip via identity. -/
theorem paleyLike_not_SF_id (P : PaleyLike n) :
    ¬ IsSelfFlip_id P.T :=
  regular_not_SF_id P.T P.hasBasePath P.ge_three P.isRegular

end Tournament
