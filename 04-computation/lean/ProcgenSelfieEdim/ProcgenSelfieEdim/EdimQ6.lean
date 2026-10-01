import ProcgenSelfieEdim.Hypercube

set_option autoImplicit false

/-!
# THM-4525 (a): the paper's 15-set resolves `Q_6`

Allikvere's certificate for `edim_m(Q_6) ≤ 15` (arXiv:2608.09983, Table 1) is the
vertex mask `0x02283022a042a00a`, i.e. the vertices
`{1, 3, 13, 15, 17, 22, 29, 31, 33, 37, 44, 45, 51, 53, 57}` (bit `v` of the mask
is vertex `v`). The kernel evaluates the checker `checkResolving`, whose soundness
is `resolving_of_check`.
-/

namespace ProcgenSelfieEdim

/-- The paper's resolving 15-set of `Q_6`. -/
def paperSet : List Nat := [1, 3, 13, 15, 17, 22, 29, 31, 33, 37, 44, 45, 51, 53, 57]

/-- The list is a set of 15 vertices of `Q_6`, and it is the paper's mask. -/
theorem paperSet_is_vertex_set :
    paperSet.length = 15 ∧ paperSet.Nodup ∧ (∀ s ∈ paperSet, s < 2 ^ 6) ∧
      paperSet = (List.range 64).filter (fun v => Nat.testBit 0x02283022a042a00a v) := by
  decide

theorem check_paperSet : checkResolving 6 paperSet = true := by
  decide +kernel

/-- `Q_6` has `6 · 2^5 = 192` edges, listed without repetition. -/
theorem edgeList_six : (edgeList 6).length = 192 ∧ (edgeList 6).Nodup :=
  ⟨by decide +kernel, nodup_of_map_nodup _ _ (nodup_of_distinctR _ check_paperSet)⟩

/-- **THM-4525 (a).** The paper's 15-set is edge-multiset resolving in `Q_6`. -/
theorem paperSet_resolving : Resolving 6 paperSet :=
  resolving_of_check 6 paperSet check_paperSet

end ProcgenSelfieEdim
