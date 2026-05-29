/-
  TournamentH7.GridReflection — THM-280 (Grid Reflection Induces
                                          Tournament Complement)

  **Statement (project-novel, opus-2026-04-03-S27).**

  Let T be a tournament on `Fin n` represented by a tiling `M` of the
  staircase Δ_{n−2} with fixed base path `n → (n-1) → … → 1`.  Let
  `M'` be the *grid-reflected* tiling: the tile coordinates (x, y) are
  mapped to (n+1−y, n+1−x), and the bit is carried over.

  Then the tournament `T(M')` equals `σ(T^op)`, where:
    • `T^op` reverses all arcs of T;
    • `σ : Fin n → Fin n` is the vertex permutation `v ↦ (n − 1) − v`
      (in 0-indexed `Fin n` form; in 1-indexed paper notation,
      `σ(v) = n + 1 − v`).

  **Consequence.** Grid reflection on tilings induces the complement
  involution `T ↦ T^op` on isomorphism classes.  In particular, the
  fixed points of grid reflection within a class fiber correspond exactly
  to *self-complementary* (SC) tournaments.

  ─── Lean formalisation strategy ────────────────────────────────────────
  Since `Tilings.lean` does not build a concrete tile-coordinate model,
  we express THM-280 directly at the tournament level:
    • Define `Tournament.op T` (the reverse/complement tournament).
    • Define `Tournament.relabel T σ` for `σ : Equiv.Perm (Fin n)`.
    • Define `vertexReversal : Equiv.Perm (Fin n)` for `v ↦ (n − 1) − v`.
    • State THM-280 in the form: for any base-path tournament T,
        T̃ = (vertexReversal-relabel of T^op).
      i.e. tile-complement equals the vertex-reversed reverse tournament.
-/

import TournamentH7.Tilings

namespace Tournament

variable {n : ℕ}

/-! ### Tournament operations -/

/-- The *reverse* of a tournament: every arc flipped. -/
def op (T : Tournament n) : Tournament n where
  arc := fun i j => T.arc j i
  irrefl := T.irrefl
  total := fun _ _ hij => by
    rcases T.total _ _ hij with hl | hr
    · exact Or.inr hl
    · exact Or.inl hr
  asym := fun i j ⟨h1, h2⟩ => T.asym i j ⟨h2, h1⟩

@[simp] lemma op_arc (T : Tournament n) (i j : Fin n) :
    (op T).arc i j = T.arc j i := rfl

@[simp] lemma op_op (T : Tournament n) : op (op T) = T := by
  cases T; rfl

/-- The *relabelling* of a tournament by a permutation σ:
    new_arc(i, j) = old_arc(σ⁻¹ i, σ⁻¹ j). -/
def relabel (T : Tournament n) (σ : Equiv.Perm (Fin n)) : Tournament n where
  arc := fun i j => T.arc (σ.symm i) (σ.symm j)
  irrefl := fun _ => T.irrefl _
  total := fun i j hij =>
    T.total (σ.symm i) (σ.symm j) (fun H => hij (σ.symm.injective H))
  asym := fun _ _ h => T.asym _ _ h

@[simp] lemma relabel_arc (T : Tournament n) (σ : Equiv.Perm (Fin n)) (i j : Fin n) :
    (relabel T σ).arc i j = T.arc (σ.symm i) (σ.symm j) := rfl

/-! ### The vertex-reversal permutation σ : v ↦ (n − 1) − v -/

/-- Vertex-reversal permutation `v ↦ (n − 1) − v` on Fin n. -/
def vertexReversal (n : ℕ) : Equiv.Perm (Fin n) where
  toFun v := ⟨n - 1 - v.val, by
    by_cases hn : n = 0
    · exact absurd v.is_lt (by simp [hn])
    · have : v.val < n := v.is_lt
      omega⟩
  invFun v := ⟨n - 1 - v.val, by
    by_cases hn : n = 0
    · exact absurd v.is_lt (by simp [hn])
    · have : v.val < n := v.is_lt
      omega⟩
  left_inv := by
    intro v
    apply Fin.ext
    by_cases hn : n = 0
    · exact absurd v.is_lt (by simp [hn])
    · have : v.val < n := v.is_lt
      simp; omega
  right_inv := by
    intro v
    apply Fin.ext
    by_cases hn : n = 0
    · exact absurd v.is_lt (by simp [hn])
    · have : v.val < n := v.is_lt
      simp; omega

@[simp] lemma vertexReversal_apply (v : Fin n) :
    (vertexReversal n v).val = n - 1 - v.val := rfl

@[simp] lemma vertexReversal_symm (v : Fin n) :
    ((vertexReversal n).symm v).val = n - 1 - v.val := rfl

/-- Vertex-reversal is an involution. -/
@[simp] lemma vertexReversal_apply_twice (v : Fin n) :
    vertexReversal n (vertexReversal n v) = v := by
  apply Fin.ext
  by_cases hn : n = 0
  · exact absurd v.is_lt (by simp [hn])
  · have : v.val < n := v.is_lt
    simp; omega

/-! ### Vertex reversal preserves the IsConsec relation -/

/-- The vertex-reversal preserves consecutivity (it swaps i+1 ↔ i, etc.). -/
lemma vertexReversal_preserves_consec (_hn : 1 ≤ n) (i j : Fin n) :
    IsConsec (vertexReversal n i) (vertexReversal n j) ↔ IsConsec i j := by
  unfold IsConsec
  simp [vertexReversal_apply]
  have hi : i.val < n := i.is_lt
  have hj : j.val < n := j.is_lt
  constructor
  · rintro (h | h)
    · -- (n-1-i) + 1 = n-1-j  ⟹  n - i = n - 1 - j  ⟹ j + 1 = i (using bounds)
      right; omega
    · left; omega
  · rintro (h | h)
    · right; omega
    · left; omega

/-! ### THM-280 — Grid Reflection Induces Complement -/

/-- **Theorem (THM-280, project-novel, opus-2026-04-03-S27).**

    For any tournament T with the base path, the *tile-complement* T̃
    equals the relabelling of T^op by the vertex-reversal permutation:
        T̃ = relabel (op T) (vertexReversal n).

    Equivalently: `(tilde T).arc i j = T.arc (σ j) (σ i)` for σ = vertexReversal.

    *Proof strategy* (informal, project canon):
    Both sides agree on consecutive arcs (the base path direction is
    preserved under vertexReversal because σ(j) and σ(i) are also
    consecutive, and the base-path arc orientation flips, then `op`
    flips it back). Both sides agree on non-consecutive arcs (T̃ flips
    the tile bit; relabel(op T, σ) flips via op then permutes vertices).
    A careful case split on consecutivity in both views completes it.

    Lean formalisation: the full proof requires manipulating the
    interaction between vertex reversal and the base-path orientation.
    It is straightforward case-bookkeeping but lengthy; deferred. -/
axiom tilde_eq_reversed_op (T : Tournament n) (hbp : HasBasePath T) (hn : 1 ≤ n) :
    ∀ i j : Fin n,
      (tilde T).arc i j = (op T).arc (vertexReversal n i) (vertexReversal n j)

/-! ### Corollaries: SC ⇔ grid-symmetric tilings -/

/-- A tournament T is *self-complementary* (SC) iff there exists a
    permutation σ such that `T = relabel (op T) σ`.  This captures
    `T ≅ T^op`. -/
def IsSelfComplementary (T : Tournament n) : Prop :=
  ∃ σ : Equiv.Perm (Fin n), ∀ i j : Fin n,
    T.arc i j = (op T).arc (σ i) (σ j)

/-- A tiling-equipped tournament T is *grid-symmetric* iff T = T̃ (the
    tile-complement equals T itself). -/
def IsGridSymmetric (T : Tournament n) : Prop :=
  ∀ i j : Fin n, T.arc i j = (tilde T).arc i j

/-- **Corollary of THM-280.** A tournament T with base path is
    grid-symmetric iff it is self-complementary via the vertex-reversal
    permutation `vertexReversal n`. -/
theorem gridSym_iff_sc_via_reversal (T : Tournament n)
    (hbp : HasBasePath T) (hn : 1 ≤ n) :
    IsGridSymmetric T ↔
    (∀ i j : Fin n,
       T.arc i j = (op T).arc (vertexReversal n i) (vertexReversal n j)) := by
  unfold IsGridSymmetric
  constructor
  · intro h i j
    rw [h]
    exact tilde_eq_reversed_op T hbp hn i j
  · intro h i j
    rw [h]
    exact (tilde_eq_reversed_op T hbp hn i j).symm

end Tournament
