import Mathlib

/-!
# The owner-local pre-nerve projection principle

The exact common-scale certificates first allow unit choices to vary with the
owner.  A global covariant unit word therefore maps to a feasible choice in
every owner projection, but the converse is false in general.  This file
kernel-checks both the one-way implication used by the c=14,...,20 terminal
certificates and a two-owner counterexample to the converse.
-/

namespace LRCPreNerveProjection

universe uOwner uGlobal uLocal

variable {Owner : Type uOwner}
variable {Global : Type uGlobal}
variable {Local : Owner → Type uLocal}

/-- One shared global word meets every owner obligation after projection. -/
def GloballyFeasible
    (project : (owner : Owner) → Global → Local owner)
    (good : (owner : Owner) → Local owner → Prop) : Prop :=
  ∃ word, ∀ owner, good owner (project owner word)

/-- Each owner has a witness, with no requirement that the witnesses glue. -/
def OwnerLocallyFeasible
    (good : (owner : Owner) → Local owner → Prop) : Prop :=
  ∀ owner, ∃ localWord, good owner localWord

/-- A global covariant word supplies an owner-local witness at every owner. -/
theorem ownerLocallyFeasible_of_globallyFeasible
    (project : (owner : Owner) → Global → Local owner)
    (good : (owner : Owner) → Local owner → Prop)
    (hglobal : GloballyFeasible project good) :
    OwnerLocallyFeasible good := by
  rcases hglobal with ⟨word, hword⟩
  intro owner
  exact ⟨project owner word, hword owner⟩

/-- One empty owner projection is an immediate global obstruction. -/
theorem not_globallyFeasible_of_empty_owner
    (project : (owner : Owner) → Global → Local owner)
    (good : (owner : Owner) → Local owner → Prop)
    (owner : Owner)
    (hempty : ∀ localWord, ¬good owner localWord) :
    ¬GloballyFeasible project good := by
  rintro ⟨word, hword⟩
  exact hempty (project owner word) (hword owner)

/-- The global obligation at one owner is the preimage of its local one. -/
theorem global_obligation_eq_preimage
    (project : (owner : Owner) → Global → Local owner)
    (good : (owner : Owner) → Local owner → Prop)
    (owner : Owner) :
    {word : Global | good owner (project owner word)} =
      project owner ⁻¹' {localWord : Local owner | good owner localWord} := by
  rfl

/-- Empty local obligation sets have empty global preimages. -/
theorem global_obligation_empty_of_local_empty
    (project : (owner : Owner) → Global → Local owner)
    (good : (owner : Owner) → Local owner → Prop)
    (owner : Owner)
    (hempty : {localWord : Local owner | good owner localWord} = ∅) :
    {word : Global | good owner (project owner word)} = ∅ := by
  rw [global_obligation_eq_preimage, hempty]
  exact Set.preimage_empty

section ConverseFailure

/-- Both owners see the same global bit. -/
def twoOwnerProject (_owner : Bool) (word : Bool) : Bool := word

/-- Owner `false` demands `false`; owner `true` demands `true`. -/
def twoOwnerGood (owner localWord : Bool) : Prop := localWord = owner

/-- Each owner projection is nonempty when its choice may vary by owner. -/
theorem twoOwner_locally_feasible :
    OwnerLocallyFeasible (Local := fun _ : Bool => Bool) twoOwnerGood := by
  intro owner
  exact ⟨owner, rfl⟩

/-- The two local witnesses do not glue to one shared global word. -/
theorem twoOwner_not_globally_feasible :
    ¬GloballyFeasible (Local := fun _ : Bool => Bool)
      twoOwnerProject twoOwnerGood := by
  rintro ⟨word, hword⟩
  have hfalse : word = false := hword false
  have htrue : word = true := hword true
  cases word <;> simp_all [twoOwnerProject, twoOwnerGood]

/-- Hence owner-local feasibility is necessary but not sufficient globally. -/
theorem ownerLocal_converse_fails :
    OwnerLocallyFeasible (Local := fun _ : Bool => Bool) twoOwnerGood ∧
      ¬GloballyFeasible (Local := fun _ : Bool => Bool)
        twoOwnerProject twoOwnerGood :=
  ⟨twoOwner_locally_feasible, twoOwner_not_globally_feasible⟩

end ConverseFailure

#print axioms ownerLocallyFeasible_of_globallyFeasible
#print axioms not_globallyFeasible_of_empty_owner
#print axioms global_obligation_empty_of_local_empty
#print axioms ownerLocal_converse_fails

end LRCPreNerveProjection
