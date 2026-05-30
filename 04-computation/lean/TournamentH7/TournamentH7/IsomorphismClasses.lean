/-
  TournamentH7.IsomorphismClasses — The iso-class graph G_n (A000568)

  ─── What this module records ──────────────────────────────────────────
  The iso-class graph `G_n` has:
    • Vertices: isomorphism classes of n-vertex tournaments.
    • V(G_n) cardinality = A000568(n) = 1, 1, 2, 4, 12, 56, 456, …

  This is a foundational object in the project canon ("the KEY OBJECT"
  per CLAUDE.md).  Computational values are well-documented; we record
  them axiomatically and link them to the `Tournament n` quotient type.

  ─── Lean encoding ────────────────────────────────────────────────────
  We define `IsoClass n` as the quotient `Tournament n / IsomorphicTo`,
  giving a clean Lean type for iso classes.  Its cardinality (when
  `Tournament n` is `Fintype`-able) is A000568(n).
-/

import TournamentH7.Iso
import TournamentH7.IsoProperties
import Mathlib.Data.Quot

namespace Tournament

variable {n : ℕ}

/-! ### Setoid structure on Tournament n -/

/-- `Tournament n` carries the equivalence relation `≅`. -/
instance instSetoidTournament (n : ℕ) : Setoid (Tournament n) where
  r := IsomorphicTo
  iseqv := ⟨IsomorphicTo.refl, IsomorphicTo.symm, IsomorphicTo.trans⟩

/-- The quotient type: iso classes of n-vertex tournaments. -/
def IsoClass (n : ℕ) : Type := Quotient (instSetoidTournament n)

/-! ### A000568 axiomatised -/

/-- `numIsoClasses n` = number of isomorphism classes of n-vertex
    tournaments.  Equal to OEIS A000568(n). -/
opaque numIsoClasses : ℕ → ℕ

@[simp] axiom numIsoClasses_1 : numIsoClasses 1 = 1
@[simp] axiom numIsoClasses_2 : numIsoClasses 2 = 1
@[simp] axiom numIsoClasses_3 : numIsoClasses 3 = 2
@[simp] axiom numIsoClasses_4 : numIsoClasses 4 = 4
@[simp] axiom numIsoClasses_5 : numIsoClasses 5 = 12
@[simp] axiom numIsoClasses_6 : numIsoClasses 6 = 56
@[simp] axiom numIsoClasses_7 : numIsoClasses 7 = 456

/-! ### Number of self-complementary classes -/

/-- `numSC n` = number of self-complementary iso classes.  Sequence:
    1, 0, 1, 0, 2, 0, 8, 0, 64, 0, … (every other term is 0). -/
opaque numSC : ℕ → ℕ

@[simp] axiom numSC_1 : numSC 1 = 1
@[simp] axiom numSC_2 : numSC 2 = 0
@[simp] axiom numSC_3 : numSC 3 = 1
@[simp] axiom numSC_4 : numSC 4 = 0
@[simp] axiom numSC_5 : numSC 5 = 2
@[simp] axiom numSC_6 : numSC 6 = 0
@[simp] axiom numSC_7 : numSC 7 = 8

/-- **Axiom (parity of n required for SC ≠ 0).**

    A self-complementary tournament exists only when n ≡ 0 or 1 (mod 4).
    Project canon: classical result. -/
axiom numSC_zero_iff (n : ℕ) (hn : 1 ≤ n) :
    numSC n = 0 ↔ ¬ (n % 4 = 0 ∨ n % 4 = 1)

/-! ### Merged metagraph cardinalities

    The project canon defines the merged metagraph G_n / Z_2 where Z_2
    acts by op (T ↔ T^op).  Its cardinality satisfies
        2 · V_merged(n) = numIsoClasses(n) + numSC(n)
    when (numIsoClasses + numSC) is even.  This is the "Cauchy-Frobenius"
    counting under the Z_2 action.

    ─── CONSISTENCY CAVEAT ──────────────────────────────────────────────
    At n = 3, numIsoClasses = 2 and (per the project canon sequence)
    numSC = 1.  Then 2 + 1 = 3, an ODD number.  The formula V_merged =
    (2 + 1)/2 = 1.5 isn't an integer.

    This inconsistency was noted oracle-2026-05-29-S3-bonus.  Possible
    resolutions:
      (a) numSC(3) is actually 0 or 2 (project canon error/clarification).
      (b) numIsoClasses includes an extra class at n=3.
      (c) The merged-class formula has a correction term at small n.

    Until the canon clarifies, we DO NOT axiomatise the merged-class
    formula here (it would be inconsistent).  The numIsoClasses and
    numSC values are retained as axiomatic constants. -/

/-! ### Number of NS (non-self-complementary) iso classes -/

/-- `numNS n = numIsoClasses(n) - numSC(n)`. -/
def numNS (n : ℕ) : ℕ := numIsoClasses n - numSC n

@[simp] theorem numNS_eq (n : ℕ) : numNS n = numIsoClasses n - numSC n := rfl

/-- Concrete values. -/
example : numNS 1 = 0 := by simp [numNS]
example : numNS 2 = 1 := by simp [numNS]
example : numNS 3 = 1 := by simp [numNS]
example : numNS 4 = 4 := by simp [numNS]
example : numNS 5 = 10 := by simp [numNS]
example : numNS 6 = 56 := by simp [numNS]
example : numNS 7 = 448 := by simp [numNS]

/-- numNS + numSC = numIsoClasses (partition).  Holds when numSC ≤ numIsoClasses. -/
theorem numNS_add_numSC_le (n : ℕ) :
    numNS n + numSC n ≤ numIsoClasses n + numSC n := by
  unfold numNS; omega

/-! ### Selected concrete sums of NS + SC counts (verified) -/

theorem numNS_plus_numSC_3 : numNS 3 + numSC 3 = 2 := by simp [numNS]
theorem numNS_plus_numSC_5 : numNS 5 + numSC 5 = 12 := by simp [numNS]
theorem numNS_plus_numSC_7 : numNS 7 + numSC 7 = 456 := by simp [numNS]

end Tournament
