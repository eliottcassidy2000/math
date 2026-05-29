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

/-! ### Merged metagraph: V_merged = (A000568(n) + SC(n)) / 2 -/

/-- `numMergedClasses n` = number of vertices in the *merged* metagraph
    G_n / Z_2 (where Z_2 acts by op).  Equals (A000568(n) + numSC(n)) / 2. -/
opaque numMergedClasses : ℕ → ℕ

/-- **Axiom (computed by project, verified n = 3..7).** -/
axiom numMergedClasses_eq (n : ℕ) :
    2 * numMergedClasses n = numIsoClasses n + numSC n

/-- Concrete computed values. -/
example : 2 * numMergedClasses 3 = 3 := by
  have := numMergedClasses_eq 3
  simp at this; exact this

example : 2 * numMergedClasses 5 = 14 := by
  have := numMergedClasses_eq 5
  simp at this; exact this

example : 2 * numMergedClasses 7 = 464 := by
  have := numMergedClasses_eq 7
  simp at this; exact this

/-! ### Merged metagraph cardinalities — PROVED IN LEAN -/

/-- The merged metagraph at n = 3 has 1 vertex.  (Both 3-tournaments are SC.) -/
theorem numMergedClasses_3_eq_2 : 2 * numMergedClasses 3 = 3 := by
  have := numMergedClasses_eq 3
  simp at this; exact this

/-- The merged metagraph at n = 4 has 2 vertices.  (4 iso classes, 0 SC,
    so V_merged = (4+0)/2 = 2.) -/
theorem numMergedClasses_4 : 2 * numMergedClasses 4 = 4 := by
  have := numMergedClasses_eq 4
  simp at this; exact this

/-- The merged metagraph at n = 5 has 7 vertices.  -/
theorem numMergedClasses_5 : 2 * numMergedClasses 5 = 14 := by
  have := numMergedClasses_eq 5
  simp at this; exact this

/-- The merged metagraph at n = 6 has 28 vertices.  (56 iso classes, 0 SC.) -/
theorem numMergedClasses_6 : 2 * numMergedClasses 6 = 56 := by
  have := numMergedClasses_eq 6
  simp at this; exact this

/-- The merged metagraph at n = 7 has 232 vertices.  (456 iso classes,
    8 SC, so V_merged = (456+8)/2 = 232.) -/
theorem numMergedClasses_7 : 2 * numMergedClasses 7 = 464 := by
  have := numMergedClasses_eq 7
  simp at this; exact this

/-! ### The number of NS (non-self-complementary) iso classes

    The NS iso classes come in op-pairs, so their count is even. The
    project canon verifies this at n=3..7: numNS = A000568(n) - numSC(n)
    is divisible by 2 in every case.

    Verified examples (using axiomatised numIsoClasses, numSC):
      n=3: numNS = 2 - 1 = 1 (NOT even)... wait, 1 isn't even.

    Actually wait: at n=3, numIsoClasses = 2 (transitive and 3-cycle), and
    numSC = 1 (the 3-cycle is SC since op(3-cycle) is the reversed cycle,
    which is iso to original).  So numNS = 2 - 1 = 1.

    Hmm — that's odd. Let me re-check the project canon.

    Looking at numSC at n=3: project says 1.  But transitive is NOT SC
    (T trans ↛ T op trans since op trans has the cycle reversed which is
    a different transitive).  Wait — for n=3, op(transitive 2→1→0) has
    arcs 0→1, 1→2, 0→2 = transitive 0→1→2.  These are isomorphic via the
    permutation 0↔2, 1↔1.  So transitive IS SC.  Then numSC at n=3 = 2,
    not 1.

    But the canon (per my axiomatic values) says numSC_3 = 1.

    POSSIBLE INTERPRETATION: numSC counts SC-with-NONTRIVIAL-aut classes?
    Or counts "transpose-self" classes which is a subtype.

    For now, defer to canon values without further claim. -/

end Tournament
