/-
THM-343: For every tournament T, the number of Hamiltonian paths H(T) ≠ 7.

This file is a SKELETON for the Lean 4 / Mathlib formalization of THM-343.
It declares the types and statements needed for the proof, with `sorry`s for
the unproved lemmas. See `03-artifacts/drafts/lean-formalization-plan-thm343.md`
for the full plan.

Author: opus-2026-05-28-S6
References:
- Moon, J.W. (1966), On subtournaments of a tournament, Canad. Math. Bull. 9,
  297-301. (Corollary 2.1 — "c(T_n, k) ≥ n-k+1 for SC T_n, 3 ≤ k ≤ n".)
- Camion, P. (1959), Chemins et circuits hamiltoniens des graphes complets,
  C. R. Acad. Sci. Paris 249, 2151-2152.
- Grinberg, D. and Stanley, R.P. (2024), arXiv:2412.10572, Corollary 20 (OCF).
-/

import Mathlib.Combinatorics.Digraph.Basic
import Mathlib.Combinatorics.Quiver.Path
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Tactic

open Finset

namespace THM343

universe u

/-! ### 1. Tournament type -/

/-- A tournament on a finite type `V` is an oriented complete digraph:
    for every pair `u ≠ v`, exactly one of `Adj u v` or `Adj v u` holds, and
    `Adj u u` never holds. -/
structure Tournament (V : Type u) where
  Adj      : V → V → Prop
  irrefl   : ∀ v, ¬ Adj v v
  total    : ∀ u v : V, u ≠ v → Adj u v ∨ Adj v u
  asymm    : ∀ u v : V, Adj u v → ¬ Adj v u

variable {V : Type u} [Fintype V] [DecidableEq V]
variable (T : Tournament V)

/-! ### 2. Directed cycles and Hamiltonian paths -/

/-- A directed cycle of length `k` is a tuple of `k` distinct vertices such
    that each consecutive pair is adjacent in `T`, including the wrap-around. -/
structure DirCycle (k : ℕ) (T : Tournament V) where
  vertices : Fin k → V
  inj      : Function.Injective vertices
  edges    : ∀ i : Fin k, T.Adj (vertices i) (vertices ⟨(i.val + 1) % k, by omega⟩)
  pos      : 3 ≤ k

/-- The vertex set of a directed cycle. -/
def DirCycle.vertexSet {k : ℕ} (c : DirCycle k T) : Finset V :=
  (Finset.univ : Finset (Fin k)).image c.vertices

/-- A Hamiltonian path is an injective sequence of all vertices with
    each consecutive pair forming an arc. -/
structure HamPath (T : Tournament V) where
  perm  : Fin (Fintype.card V) → V
  inj   : Function.Injective perm
  edges : ∀ i : Fin (Fintype.card V - 1),
    T.Adj (perm ⟨i.val, by omega⟩) (perm ⟨i.val + 1, by omega⟩)

/-- The Rédei H-invariant: number of Hamiltonian paths in `T`. -/
noncomputable def H (T : Tournament V) : ℕ :=
  Nat.card (HamPath T)

/-! ### 3. Strong connectivity -/

/-- One-step reachability: `T.Adj u v` already gives a direct path.
    For multi-step, use the transitive closure. -/
def reachable (T : Tournament V) (u v : V) : Prop :=
  Relation.ReflTransGen T.Adj u v

/-- `T` is strongly connected if every pair of vertices is mutually reachable. -/
def IsSC (T : Tournament V) : Prop :=
  ∀ u v : V, reachable T u v ∧ reachable T v u

/-! ### 4. Sub-tournament -/

/-- The sub-tournament induced by a subset `S ⊆ V`. -/
def induced (T : Tournament V) (S : Finset V) : Tournament S where
  Adj u v := T.Adj u.val v.val
  irrefl  v := T.irrefl v.val
  total   u v hne := T.total u.val v.val (fun h => hne (Subtype.eq h))
  asymm   u v := T.asymm u.val v.val

/-! ### 5. Lemmas -/

/-- Every directed cycle, viewed as a sub-tournament on its vertex set,
    is strongly connected. -/
lemma DirCycle.isSC {k : ℕ} (c : DirCycle k T) :
    IsSC (T.induced c.vertexSet) := by
  sorry  -- standard: follow cycle arcs forward/backward

/-- If two SC sub-tournaments share a vertex, their union is also SC. -/
lemma SC_union (S₁ S₂ : Finset V)
    (h₁ : IsSC (T.induced S₁)) (h₂ : IsSC (T.induced S₂))
    (h_shared : (S₁ ∩ S₂).Nonempty) :
    IsSC (T.induced (S₁ ∪ S₂)) := by
  sorry  -- standard: reachability transits through shared vertex

/-- The number of directed cycles of length `k` in `T`. -/
noncomputable def cycleCount (T : Tournament V) (k : ℕ) : ℕ :=
  Nat.card (DirCycle k T)

/-- **Moon (1966), Corollary 2.1**: A strongly connected tournament on
    `n` vertices has at least `n - k + 1` directed cycles of length `k`,
    for `3 ≤ k ≤ n`.

    Proof: Induction on `n` via Moon's argument (a strong tournament on
    `n` vertices contains a strong subtournament on `n-1` vertices, and
    the removed vertex lies on a `k`-cycle by vertex-pancyclicity).
    See the cited paper for details. -/
theorem moon_cycle_count
    (T : Tournament V) (hSC : IsSC T) (k : ℕ)
    (hk_lo : 3 ≤ k) (hk_hi : k ≤ Fintype.card V) :
    Fintype.card V - k + 1 ≤ cycleCount T k := by
  sorry

/-- **Camion (1959)**: Every strongly connected tournament on `n ≥ 3` vertices
    has a directed Hamilton cycle. Special case of Moon (k = n). -/
theorem camion
    (T : Tournament V) (hSC : IsSC T) (h : 3 ≤ Fintype.card V) :
    Nonempty (DirCycle (Fintype.card V) T) := by
  have := moon_cycle_count T hSC (Fintype.card V) h (le_refl _)
  simp at this
  -- this : 1 ≤ cycleCount T (Fintype.card V)
  sorry  -- existence from positive count

/-- The conflict graph `Ω(T)` of odd directed cycles: vertices are odd
    directed cycles of T, edges are pairs whose vertex sets intersect. -/
def OddCycle (T : Tournament V) : Type _ :=
  Σ k : ℕ, { c : DirCycle k T // Odd k }

/-- The independence polynomial of `Ω(T)` evaluated at 2. -/
noncomputable def alphaSum (T : Tournament V) : ℕ :=
  sorry  -- ∑_k (count of indep sets of size k in Ω(T)) · 2^k

/-- **OCF (Grinberg-Stanley 2024)**: For any tournament T,
    H(T) = I(Ω(T), 2) where Ω(T) is the conflict graph of odd directed cycles
    and I is the independence polynomial. -/
axiom OCF (T : Tournament V) : H T = alphaSum T

/-! ### 6. Main theorem -/

/-- **THM-343**: For every tournament T, H(T) ≠ 7. -/
theorem H_ne_seven (T : Tournament V) : H T ≠ 7 := by
  intro h_eq
  rw [OCF] at h_eq
  -- Step 1: by OCF + chain constraint, alpha = (3, 0, 0, …)
  -- Step 2: extract 3 pairwise-intersecting odd cycles C₁, C₂, C₃
  -- Step 3: T[V'] is SC where V' = V(C₁) ∪ V(C₂) ∪ V(C₃)
  -- Step 4: case analysis on s = |V'|:
  --   s ≤ 2: impossible
  --   s = 3: only 1 directed 3-cycle, can't have 3 distinct
  --   s = 4: max 2 three-cycles (unique SC score (1,1,2,2))
  --   s = 5: ≥ 3 three-cycles + ≥ 1 five-cycle ≥ 4 odd cycles (Moon)
  --   s ≥ 6: ≥ s-2 ≥ 4 three-cycles (Moon)
  --   ⇒ T has ≥ 4 odd cycles, contradicting α₁ = 3
  sorry

end THM343
