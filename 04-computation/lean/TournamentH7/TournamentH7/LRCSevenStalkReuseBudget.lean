import TournamentH7.LRCSevenOverlapRelations

/-!
# Exact reuse budgets for rooted seven-stalk aggregation

The colored dense-core ledger sums local charges over six-element lower faces
attached to one root.  This module supplies the missing finite double count:
every `k`-face is reused in exactly `choose (m-k) (6-k)` such lower faces.
For at most twelve lower vertices, a root spoke is therefore reused at most
`choose 11 5 = 462` times and a lower pair at most `choose 10 4 = 210` times.

These are transport budgets, not a payment theorem: determinant magnitude and
zero-color gluing must still beat the stated reuse factors.  The proof-relevant
vertices are lower faces, the binary observable is containment, and the tie
Hamiltonian path is lexicographic inclusion order.  A runner tournament loses
face multiplicity; the face-incidence bipartite graph preserves it exactly.
-/

namespace LonelyRunner
namespace LRCSevenStalkReuseBudget

open Finset
open LRC14Concrete
open scoped BigOperators

/-- Exact face-incidence double count.  Summing a `k`-face weight through all
`r`-faces multiplies its global mass by `choose (|S|-k) (r-k)`. -/
theorem sum_powersetCard_sum_powersetCard
    {α : Type*} [DecidableEq α] (vertices : Finset α)
    (weight : Finset α → ℕ) (k r : ℕ) (hkr : k ≤ r) :
    (∑ face ∈ vertices.powersetCard r,
        ∑ subface ∈ face.powersetCard k, weight subface) =
      (vertices.card - k).choose (r - k) *
        ∑ subface ∈ vertices.powersetCard k, weight subface := by
  calc
    (∑ face ∈ vertices.powersetCard r,
        ∑ subface ∈ face.powersetCard k, weight subface) =
        ∑ face ∈ vertices.powersetCard r,
          ∑ subface ∈ vertices.powersetCard k,
            if subface ⊆ face then weight subface else 0 := by
      apply Finset.sum_congr rfl
      intro face hface
      have hfaceVertices : face ⊆ vertices :=
        (Finset.mem_powersetCard.mp hface).1
      have hsub : face.powersetCard k ⊆ vertices.powersetCard k := by
        intro subface hsubface
        rw [Finset.mem_powersetCard] at hsubface ⊢
        exact ⟨hsubface.1.trans hfaceVertices, hsubface.2⟩
      exact Finset.sum_subset_zero_on_sdiff hsub
        (by
          intro subface houtside
          have hsubfaceVertices := (Finset.mem_sdiff.mp houtside).1
          have hsubfaceOutside := (Finset.mem_sdiff.mp houtside).2
          by_cases hcontained : subface ⊆ face
          · exact False.elim (hsubfaceOutside
              (Finset.mem_powersetCard.mpr
                ⟨hcontained,
                  (Finset.mem_powersetCard.mp hsubfaceVertices).2⟩))
          · simp [hcontained])
        (by
          intro subface hsubface
          simp [(Finset.mem_powersetCard.mp hsubface).1])
    _ = ∑ subface ∈ vertices.powersetCard k,
          ∑ face ∈ vertices.powersetCard r,
            if subface ⊆ face then weight subface else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ subface ∈ vertices.powersetCard k,
          (vertices.card - k).choose (r - k) * weight subface := by
      apply Finset.sum_congr rfl
      intro subface hsubface
      rw [← Finset.sum_filter]
      have hsubfaceVertices : subface ⊆ vertices :=
        (Finset.mem_powersetCard.mp hsubface).1
      have hsubfaceCard : subface.card = k :=
        (Finset.mem_powersetCard.mp hsubface).2
      have hcard :
          ((vertices.powersetCard r).filter fun face =>
              subface ⊆ face).card =
            (vertices.card - k).choose (r - k) := by
        simpa [hsubfaceCard] using
          Finset.card_filter_powersetCard_subset subface vertices r
            hsubfaceVertices (by simpa [hsubfaceCard])
      rw [Finset.sum_const_nat (m := weight subface) (fun _ _ => rfl),
        hcard]
    _ = (vertices.card - k).choose (r - k) *
          ∑ subface ∈ vertices.powersetCard k, weight subface := by
      rw [Finset.mul_sum]

/-- Exact spoke reuse through rooted six-faces. -/
theorem sum_sixFace_spokeMass
    {α : Type*} [DecidableEq α] (vertices : Finset α)
    (weight : α → ℕ) :
    (∑ face ∈ vertices.powersetCard 6,
        ∑ vertex ∈ face, weight vertex) =
      (vertices.card - 1).choose 5 * ∑ vertex ∈ vertices, weight vertex := by
  calc
    (∑ face ∈ vertices.powersetCard 6,
        ∑ vertex ∈ face, weight vertex) =
        ∑ face ∈ vertices.powersetCard 6,
          ∑ vertex ∈ vertices,
            if vertex ∈ face then weight vertex else 0 := by
      apply Finset.sum_congr rfl
      intro face hface
      have hsub : face ⊆ vertices :=
        (Finset.mem_powersetCard.mp hface).1
      exact Finset.sum_subset_zero_on_sdiff hsub
        (by
          intro vertex hvertex
          simp [(Finset.mem_sdiff.mp hvertex).2])
        (by
          intro vertex hvertex
          simp [hvertex])
    _ = ∑ vertex ∈ vertices,
          ∑ face ∈ vertices.powersetCard 6,
            if vertex ∈ face then weight vertex else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ vertex ∈ vertices,
          (vertices.card - 1).choose 5 * weight vertex := by
      apply Finset.sum_congr rfl
      intro vertex hvertex
      rw [← Finset.sum_filter]
      have hcard :
          ((vertices.powersetCard 6).filter fun face =>
              vertex ∈ face).card =
            (vertices.card - 1).choose 5 := by
        have hsingleton : ({vertex} : Finset α) ⊆ vertices := by
          simpa using hvertex
        simpa using
          Finset.card_filter_powersetCard_subset
            ({vertex} : Finset α) vertices 6 hsingleton (by norm_num)
      rw [Finset.sum_const_nat (m := weight vertex) (fun _ _ => rfl),
        hcard]
    _ = (vertices.card - 1).choose 5 *
          ∑ vertex ∈ vertices, weight vertex := by
      rw [Finset.mul_sum]

/-- Exact lower-pair reuse through rooted six-faces. -/
theorem sum_sixFace_pairMass
    {α : Type*} [DecidableEq α] (vertices : Finset α)
    (pairWeight : Finset α → ℕ) :
    (∑ face ∈ vertices.powersetCard 6,
        ∑ pair ∈ face.powersetCard 2, pairWeight pair) =
      (vertices.card - 2).choose 4 *
        ∑ pair ∈ vertices.powersetCard 2, pairWeight pair := by
  exact sum_powersetCard_sum_powersetCard vertices pairWeight 2 6 (by omega)

/-- With at most twelve lower vertices, one spoke can occur in at most 462
rooted six-faces. -/
theorem spokeReuseFactor_le_462 (vertexCount : ℕ)
    (hvertexCount : vertexCount ≤ 12) :
    (vertexCount - 1).choose 5 ≤ 462 := by
  interval_cases vertexCount <;> decide

/-- With at most twelve lower vertices, one lower pair can occur in at most
210 rooted six-faces. -/
theorem pairReuseFactor_le_210 (vertexCount : ℕ)
    (hvertexCount : vertexCount ≤ 12) :
    (vertexCount - 2).choose 4 ≤ 210 := by
  interval_cases vertexCount <;> decide

/-- Global spoke transport upper budget on the thirteen-runner carrier. -/
theorem sum_sixFace_spokeMass_le_462
    {α : Type*} [DecidableEq α] (vertices : Finset α)
    (hvertices : vertices.card ≤ 12) (weight : α → ℕ) :
    (∑ face ∈ vertices.powersetCard 6,
        ∑ vertex ∈ face, weight vertex) ≤
      462 * ∑ vertex ∈ vertices, weight vertex := by
  rw [sum_sixFace_spokeMass]
  exact Nat.mul_le_mul_right _
    (spokeReuseFactor_le_462 vertices.card hvertices)

/-- Global lower-pair transport upper budget on the thirteen-runner carrier. -/
theorem sum_sixFace_pairMass_le_210
    {α : Type*} [DecidableEq α] (vertices : Finset α)
    (hvertices : vertices.card ≤ 12) (pairWeight : Finset α → ℕ) :
    (∑ face ∈ vertices.powersetCard 6,
        ∑ pair ∈ face.powersetCard 2, pairWeight pair) ≤
      210 * ∑ pair ∈ vertices.powersetCard 2, pairWeight pair := by
  rw [sum_sixFace_pairMass]
  exact Nat.mul_le_mul_right _
    (pairReuseFactor_le_210 vertices.card hvertices)

/-! ## Concrete colored-wall carrier -/

/-- Bad lower vertices at one multiplier event, with the root deleted. -/
def badLowerVertices (v : Fin 13 → ℤ) (q p : ℕ)
    (root : Fin 13) : Finset (Fin 13) :=
  (Finset.univ.erase root).filter fun index => ¬ inBand v q p index

/-- Absolute determinant mass on one root spoke. -/
def rootSpokeMass (v : Fin 13 → ℤ) (q p : ℕ)
    (root index : Fin 13) : ℕ :=
  (overlapDet v q p index root).natAbs

theorem badLowerVertices_card_le_twelve
    (v : Fin 13 → ℤ) (q p : ℕ) (root : Fin 13) :
    (badLowerVertices v q p root).card ≤ 12 := by
  calc
    (badLowerVertices v q p root).card ≤
        (Finset.univ.erase root).card := by
      exact Finset.card_le_card (Finset.filter_subset _ _)
    _ = 12 := by simp

/-- Exact reuse identity on the actual rooted colored-wall carrier. -/
theorem badLower_sixFace_spokeMass_eq
    (v : Fin 13 → ℤ) (q p : ℕ) (root : Fin 13) :
    (∑ face ∈ (badLowerVertices v q p root).powersetCard 6,
        ∑ index ∈ face, rootSpokeMass v q p root index) =
      ((badLowerVertices v q p root).card - 1).choose 5 *
        ∑ index ∈ badLowerVertices v q p root,
          rootSpokeMass v q p root index := by
  exact sum_sixFace_spokeMass
    (badLowerVertices v q p root) (rootSpokeMass v q p root)

/-- Concrete global transport budget: after summing through every rooted
six-face at one multiplier, no spoke determinant is charged more than 462
times. -/
theorem badLower_sixFace_spokeMass_le_462
    (v : Fin 13 → ℤ) (q p : ℕ) (root : Fin 13) :
    (∑ face ∈ (badLowerVertices v q p root).powersetCard 6,
        ∑ index ∈ face, rootSpokeMass v q p root index) ≤
      462 * ∑ index ∈ badLowerVertices v q p root,
        rootSpokeMass v q p root index := by
  exact sum_sixFace_spokeMass_le_462
    (badLowerVertices v q p root)
    (badLowerVertices_card_le_twelve v q p root)
    (rootSpokeMass v q p root)

#print axioms sum_powersetCard_sum_powersetCard
#print axioms sum_sixFace_spokeMass
#print axioms sum_sixFace_pairMass
#print axioms sum_sixFace_spokeMass_le_462
#print axioms sum_sixFace_pairMass_le_210
#print axioms badLower_sixFace_spokeMass_eq
#print axioms badLower_sixFace_spokeMass_le_462

end LRCSevenStalkReuseBudget
end LonelyRunner
