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

/-! ### THM-280 — Grid Reflection ≠ Tile-Complement

    ─── Correction (oracle-2026-05-29-S3-final2) ───────────────────────
    The project canon distinguishes TWO different involutions on tilings:
     • *Tile-complement* (`tilde T` in this Lean library):
       flip every non-base-path arc bit.
     • *Grid reflection*: apply the tile-coordinate map
       (x, y) ↦ (n+1−y, n+1−x), permuting WHICH tile gets which bit.

    These are different operations.  Concretely, for n = 4 and tiling
    b = (b_{3,1}, b_{4,1}, b_{4,2}) = (1, 0, 0):
       tile-complement(b) = (0, 1, 1)
       grid-reflection(b) = (0, 0, 1)
    are distinct.

    THM-280 is about the GRID REFLECTION, not the tile-complement:
       T(grid-reflection(b)) ≅ σ(T(b)^op)
    where σ = vertex-reversal.

    Since this Lean library represents tilings only at the *tournament
    level* (no explicit tile coordinates), we cannot directly
    formalise THM-280 here.  A future module `StaircaseTileModel.lean`
    would introduce a concrete `Tile (n)` type and a tile-indexed
    `Tiling n := Tile n → Bool`; THM-280 lives there.

    ─── What this section provides ──────────────────────────────────────
    The corrected statement: vertex reversal on the tournament side
    interacts nicely with `op T` but NOT directly with `tilde T`.  We
    record the trivial lemma that `op` and `relabel` distribute, and
    leave THM-280 proper for future formalisation. -/

/-- Identity: (relabel (op T) σ).arc i j = T.arc (σ.symm j) (σ.symm i). -/
lemma relabel_op_arc (T : Tournament n) (σ : Equiv.Perm (Fin n)) (i j : Fin n) :
    (relabel (op T) σ).arc i j = T.arc (σ.symm j) (σ.symm i) := by
  show (op T).arc (σ.symm i) (σ.symm j) = T.arc (σ.symm j) (σ.symm i)
  exact op_arc T (σ.symm i) (σ.symm j)

/-! ### Self-complementary predicate (kept; was correct) -/

/-- A tournament T is *self-complementary* (SC) iff there exists a
    permutation σ such that `T = relabel (op T) σ`.  This captures
    `T ≅ T^op`. -/
def IsSelfComplementary (T : Tournament n) : Prop :=
  ∃ σ : Equiv.Perm (Fin n), ∀ i j : Fin n,
    T.arc i j = (op T).arc (σ i) (σ j)

end Tournament
