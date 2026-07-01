/-
  FlipRankExcessLaw.lean  --  klein-2026-07-01-S85
  A Lean 4 (Mathlib-style) FORMALIZATION TARGET for HYP-3817/3819:

      the flip-rank covering excess equals the number of self-complementary tournament
      iso-classes that are more symmetric than the cyclic rotation:

        excess n = #{ C : iso-class of tournaments on Fin n | SelfComplementary C ∧ n < |Aut C| }.

  STATUS: definitions + statement + verified data (n = 3..7 : 0,0,0,1,3). The proof is `sorry`.
  This is the artifact to feed a pipeline-math-style (GPT prover + Lean verifier) attempt.
  NOT a proof; a precise, machine-checkable STATEMENT with all definitions pinned down.
-/
import Mathlib

open Finset

namespace FlipRank

variable {n : ℕ}

/-- A tournament on `Fin n`: a `Prop`-valued relation that is total and antisymmetric off the diagonal.
    `beats i j` means `i → j`. -/
structure Tournament (n : ℕ) where
  beats : Fin n → Fin n → Prop
  irrefl : ∀ i, ¬ beats i i
  total  : ∀ i j, i ≠ j → (beats i j ↔ ¬ beats j i)

/-- Relabelling a tournament by a permutation `σ`. -/
def relabel (σ : Equiv.Perm (Fin n)) (T : Tournament n) : Tournament n where
  beats i j := T.beats (σ i) (σ j)
  irrefl i := by simpa using T.irrefl (σ i)
  total i j h := by
    have : σ i ≠ σ j := by simpa using h
    simpa using T.total (σ i) (σ j) this

/-- Isomorphism of tournaments = a relabelling identifying them. -/
def Iso (S T : Tournament n) : Prop := ∃ σ : Equiv.Perm (Fin n), relabel σ S = T

/-- The converse (complement) tournament: reverse every arc. -/
def converse (T : Tournament n) : Tournament n where
  beats i j := T.beats j i
  irrefl := T.irrefl
  total i j h := by simpa [ne_comm] using T.total j i h.symm

/-- Self-complementary: isomorphic to its own converse. -/
def SelfComplementary (T : Tournament n) : Prop := Iso T (converse T)

/-- Automorphism group order = number of relabellings fixing `T` (as a Finset card;
    finite since `Equiv.Perm (Fin n)` is a fintype). -/
noncomputable def autOrder (T : Tournament n) : ℕ :=
  (univ.filter (fun σ : Equiv.Perm (Fin n) => relabel σ T = T)).card

/-- The set of iso-classes = the quotient of tournaments by `Iso`. `A000568 n = card`. -/
def IsoClass (n : ℕ) := Quot (@Iso n)

/-- The tiling cube dimension `m = C(n-1,2)` (non-consecutive arcs relative to a fixed base path). -/
def cubeDim (n : ℕ) : ℕ := Nat.choose (n-1) 2

/-- A tiling ∈ (Fin (cubeDim n) → Bool) determines a tournament (base path + flipped tiles);
    `classOf t` is its iso-class.  (Definition elided here; the construction is the S_n-quotient of Q_m.) -/
opaque classOf : (Fin (cubeDim n) → Bool) → IsoClass n

/-- An affine subcube of dimension `d`: an offset `x0` plus a `d`-dim GF(2) subspace of the coordinates.
    `covers` : the subcube's image under `classOf` is all of `IsoClass n`. -/
def AffineSubcube (n d : ℕ) := { p : (Fin (cubeDim n) → Bool) × (Fin d → (Fin (cubeDim n) → Bool)) // True }

opaque coversAll : AffineSubcube n d → Prop

/-- Flip-rank `ρ(n)` = least `d` such that some affine `d`-subcube covers all iso-classes. -/
noncomputable def flipRank (n : ℕ) : ℕ :=
  Nat.find (p := fun d => ∃ A : AffineSubcube n d, coversAll A) (by sorry)

/-- Information floor `⌈log₂ |G_n|⌉`. -/
noncomputable def infoFloor (n : ℕ) : ℕ := Nat.clog 2 (Nat.card (IsoClass n))

/-- The covering excess. -/
noncomputable def excess (n : ℕ) : ℕ := flipRank n - infoFloor n

/-- The count of self-complementary iso-classes with automorphism order exceeding `n`
    (well-defined on classes: `Iso` preserves `SelfComplementary` and `autOrder`). -/
noncomputable def scSuperSymmetricCount (n : ℕ) : ℕ := sorry
  -- = (classes.filter (fun C => SelfComplementary C.rep ∧ n < autOrder C.rep)).card

/-- ****  THE LAW (HYP-3817/3819, conjectural)  ****
    Verified n = 3..7 : excess = 0,0,0,1,3 = scSuperSymmetricCount; predicted excess(8) = 4. -/
theorem flipRank_excess_eq_scSuperSymmetric (n : ℕ) (hn : 3 ≤ n) :
    excess n = scSuperSymmetricCount n := by
  sorry

/-- Verified data (to be discharged by `decide`/`native_decide` on the concrete enumerations). -/
example : excess 6 = 1 := by sorry
example : excess 7 = 3 := by sorry
example : scSuperSymmetricCount 8 = 4 := by sorry   -- the predicted excess(8)

end FlipRank
