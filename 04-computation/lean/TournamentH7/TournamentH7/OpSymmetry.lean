/-
  TournamentH7.OpSymmetry — H(T) = H(T^op) via Hamilton-path reversal bijection.

  The operation `op T` reverses every arc of T.  A Hamilton path σ in T
  becomes a Hamilton path in T^op when we read it backwards.  Concretely,
  the map  σ ↦ (vertexReversal n).trans σ  is an involution on
  `Equiv.Perm (Fin n)` that bijects `IsHamiltonianPath T` and
  `IsHamiltonianPath (op T)`.  Hence `H (op T) = H T`.

  Project-novel derived corollary: the H-spectrum is closed under op,
  and hence the iso-class graph G_n has a complement involution which
  preserves H.

  We do NOT axiomatize `alphaCount (op T) = alphaCount T` here because
  `alphaCount` is opaque; that fact is a separate axiomatic statement
  (call it OCF-op-symmetry) that the canon takes for granted.
-/

import TournamentH7.SCC
import TournamentH7.OCF
import TournamentH7.GridReflection
import TournamentH7.IsoCharacterizations

namespace Tournament

variable {n : ℕ}

/-! ### Pointwise effect of `(vertexReversal n).trans σ` -/

/-- The composition `(vertexReversal n).trans σ` evaluated at j equals
    σ at the reversed index. -/
lemma trans_vertexReversal_apply (σ : Equiv.Perm (Fin n)) (j : Fin n) :
    ((vertexReversal n).trans σ) j = σ (vertexReversal n j) := rfl

/-- The map σ ↦ (vertexReversal n).trans σ is an involution. -/
lemma trans_vertexReversal_involutive (σ : Equiv.Perm (Fin n)) :
    (vertexReversal n).trans ((vertexReversal n).trans σ) = σ := by
  ext j
  simp only [Equiv.trans_apply, vertexReversal_apply_twice]

/-! ### HamPath in T ⟹ reversed permutation is HamPath in op T -/

/-- **Theorem.** If σ is a Hamilton path in T, then `(vertexReversal n).trans σ`
    is a Hamilton path in `op T`. -/
theorem hampath_op_of_hampath (T : Tournament n) (σ : Equiv.Perm (Fin n))
    (h : IsHamiltonianPath T σ) :
    IsHamiltonianPath (op T) ((vertexReversal n).trans σ) := by
  intro j h_succ
  have hj_lt : j.val < n := j.is_lt
  have hn_pos : 0 < n := by omega
  -- Let ja := ⟨n - 1 - (j+1), _⟩.
  -- Then (ja + 1).val = n - 1 - j.
  have ha_lt : n - 1 - (j.val + 1) < n := by omega
  let ja : Fin n := ⟨n - 1 - (j.val + 1), ha_lt⟩
  have ha_succ : ja.val + 1 < n := by show n - 1 - (j.val + 1) + 1 < n; omega
  -- Apply HamPath condition on T at index ja.
  have h_arc := h ja ha_succ
  -- Identify the two terms of the goal with the two terms of h_arc.
  have hv1 : vertexReversal n ⟨j.val + 1, h_succ⟩ = ja := by
    apply Fin.ext; rfl
  have hv2 : vertexReversal n j = ⟨ja.val + 1, ha_succ⟩ := by
    apply Fin.ext
    show n - 1 - j.val = n - 1 - (j.val + 1) + 1
    omega
  show (op T).arc (σ (vertexReversal n j)) (σ (vertexReversal n ⟨j.val + 1, h_succ⟩)) = true
  rw [op_arc, hv1, hv2]
  exact h_arc

/-! ### H(op T) = H(T) -/

/-- **Theorem.** H(op T) = H(T).  The Hamilton path counts of a tournament
    and its arc-reversal coincide. -/
theorem H_op_eq_H (T : Tournament n) : H (op T) = H T := by
  classical
  unfold H
  symm
  apply Finset.card_bij (fun σ _ => (vertexReversal n).trans σ)
  · -- Maps a HamPath in T to a HamPath in op T.
    intro σ hσ
    rw [Finset.mem_filter] at hσ ⊢
    refine ⟨Finset.mem_univ _, ?_⟩
    show IsHamiltonianPath (op T) ((vertexReversal n).trans σ)
    exact hampath_op_of_hampath T σ hσ.2
  · -- Injective.
    intro σ₁ _ σ₂ _ heq
    apply Equiv.ext
    intro k
    -- Apply heq at (vertexReversal n k), then use involution.
    have h_pointwise : ∀ x, ((vertexReversal n).trans σ₁) x
                          = ((vertexReversal n).trans σ₂) x := fun x => by rw [heq]
    specialize h_pointwise (vertexReversal n k)
    simp only [Equiv.trans_apply, vertexReversal_apply_twice] at h_pointwise
    exact h_pointwise
  · -- Surjective.
    intro τ hτ
    rw [Finset.mem_filter] at hτ
    refine ⟨(vertexReversal n).trans τ, ?_, ?_⟩
    · rw [Finset.mem_filter]
      refine ⟨Finset.mem_univ _, ?_⟩
      show IsHamiltonianPath T ((vertexReversal n).trans τ)
      -- Apply hampath_op_of_hampath with `T := op T`, then use op_op.
      have h := hampath_op_of_hampath (op T) τ hτ.2
      rw [op_op] at h
      exact h
    · -- The image equals τ by involution.
      exact trans_vertexReversal_involutive τ

/-! ### Corollaries -/

/-- The transitive tournament's reverse also has H = 1. -/
theorem H_op_transitive_eq_one (n : ℕ) (hn : 1 ≤ n) :
    H (op (transitiveTournament n)) = 1 := by
  rw [H_op_eq_H]
  exact H_transitive_eq_one n hn

/-- H is iso-and-op invariant.  In particular `H(op T) = H(T)`. -/
theorem H_invariant_under_op (T : Tournament n) :
    H (op T) = H T := H_op_eq_H T

/-! ### op-OCF axiom for alphaCount -/

/-- **Axiom (classical).** The number of vertex-disjoint odd directed cycle
    collections is invariant under reversal `op T`.  Reason: every directed
    odd cycle of T reverses to a directed odd cycle of `op T` and vice versa;
    this involution preserves disjoint-collections of any size. -/
axiom alphaCount_op_eq (k : ℕ) (T : Tournament n) :
    alphaCount k (op T) = alphaCount k T

/-! ### Consistency check: OCF derivation gives the same H -/

/-- OCF on op T expanded. -/
example (T : Tournament n) :
    H (op T) = 1 + 2 * alphaCount 1 (op T) + 4 * alphaCount 2 (op T)
              + 8 * alphaCount 3 (op T) + 16 * alphaCount 4 (op T) := ocf (op T)

/-- Combined with H_op_eq_H + alphaCount_op_eq: ocf for T^op forces the same
    arithmetic as T (consistency check). -/
example (T : Tournament n) :
    1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
      + 8 * alphaCount 3 T + 16 * alphaCount 4 T
  = 1 + 2 * alphaCount 1 (op T) + 4 * alphaCount 2 (op T)
      + 8 * alphaCount 3 (op T) + 16 * alphaCount 4 (op T) := by
  rw [alphaCount_op_eq, alphaCount_op_eq, alphaCount_op_eq, alphaCount_op_eq]

/-! ### Transitive ≅ op(transitive) via vertex reversal

    The transitive tournament's arc-reversal `op(transitive_n)` is isomorphic
    to `transitive_n` itself: the isomorphism is the vertex-reversal
    permutation `v ↦ n - 1 - v`.

    Reason: `transitive_n.arc i j ↔ j.val < i.val`.
    `op(transitive_n).arc i j ↔ i.val < j.val`.
    After applying `σ = vertexReversal` to both vertices:
    `transitive_n.arc (σ i) (σ j) ↔ (n - 1 - j.val) < (n - 1 - i.val) ↔ i.val < j.val`.
    Match. -/

/-- The isomorphism `op(transitive_n) ≅ transitive_n` realized by vertex
    reversal. -/
def isoOpTransitive (n : ℕ) :
    TournamentIso (op (transitiveTournament n)) (transitiveTournament n) where
  perm := vertexReversal n
  arc_eq := fun i j => by
    show (op (transitiveTournament n)).arc i j
       = (transitiveTournament n).arc (vertexReversal n i) (vertexReversal n j)
    show decide (i.val < j.val)
       = decide ((vertexReversal n j).val < (vertexReversal n i).val)
    simp only [vertexReversal_apply]
    by_cases hn : n = 0
    · exact absurd i.is_lt (by simp [hn])
    · have hi : i.val < n := i.is_lt
      have hj : j.val < n := j.is_lt
      by_cases h : i.val < j.val
      · simp only [h, decide_true]
        symm; simp; omega
      · simp only [h, decide_false]
        symm; simp; omega

/-- `op(transitive_n) ≅ transitive_n`. -/
theorem op_transitive_iso (n : ℕ) :
    op (transitiveTournament n) ≅ transitiveTournament n :=
  ⟨isoOpTransitive n⟩

/-- Corollary: `H (op (transitive n)) = 1` via iso route (alternative proof). -/
theorem H_op_transitive_eq_one_via_iso (n : ℕ) (hn : 1 ≤ n) :
    H (op (transitiveTournament n)) = 1 := by
  rw [H_iso_invariant _ _ (op_transitive_iso n)]
  exact H_transitive_eq_one n hn

/-! ### op preserves iso classes -/

/-- The isomorphism `op T₁ → op T₂` induced by an isomorphism `T₁ → T₂`. -/
def TournamentIso_op (T₁ T₂ : Tournament n) (f : TournamentIso T₁ T₂) :
    TournamentIso (op T₁) (op T₂) where
  perm := f.perm
  arc_eq := fun i j => by
    show (op T₁).arc i j = (op T₂).arc (f.perm i) (f.perm j)
    rw [op_arc, op_arc]
    exact f.arc_eq j i

/-- **Theorem.** Op is iso-respecting: T₁ ≅ T₂ ⟹ op T₁ ≅ op T₂. -/
theorem op_iso_of_iso (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    op T₁ ≅ op T₂ := by
  obtain ⟨f⟩ := h
  exact ⟨TournamentIso_op T₁ T₂ f⟩

/-- Corollary: `op` is an involution on iso classes. -/
theorem op_op_iso (T : Tournament n) : op (op T) ≅ T := by
  rw [op_op]

/-! ### Iso-and-op: H(T) = 1 ⟺ op T ≅ transitive -/

/-- **Theorem.** H(T) = 1 ⟺ op T ≅ transitive.
    Combines `H_eq_one_iff_transitive` with `op_transitive_iso` and
    `op_iso_of_iso`. -/
theorem H_eq_one_iff_op_transitive (T : Tournament n) (hn : 1 ≤ n) :
    H T = 1 ↔ op T ≅ transitiveTournament n := by
  rw [H_eq_one_iff_transitive T hn]
  constructor
  · intro h
    -- T ≅ transitive ⟹ op T ≅ op(transitive) ≅ transitive
    have h1 : op T ≅ op (transitiveTournament n) := op_iso_of_iso _ _ h
    exact h1.trans (op_transitive_iso n)
  · intro h
    -- op T ≅ transitive ⟹ op(op T) ≅ op(transitive) ≅ transitive ⟹ T ≅ transitive
    have h1 : op (op T) ≅ op (transitiveTournament n) := op_iso_of_iso _ _ h
    have h2 : T ≅ op (transitiveTournament n) := by rw [op_op] at h1; exact h1
    exact h2.trans (op_transitive_iso n)

end Tournament
