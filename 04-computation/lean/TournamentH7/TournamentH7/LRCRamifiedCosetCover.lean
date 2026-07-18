import TournamentH7.LRCPreNerveProjection

/-!
# Ramified coset-cover projection and fibre bounds

This file isolates two proof principles used by the finite common-scale
certificates.

* A single shared assignment which covers every owner projects to a covering
  assignment in every owner-local bank.  This is only a necessary condition:
  the converse is deliberately not asserted.
* A union of whole fibres together with pieces meeting each fibre at most once
  is bounded by a saturated fibre score.  Consequently a full cover forces
  that score to equal the number of sheets.

The second statement is the abstract consumer behind the scale-27
`Z/27 -> Z/9` flag obstruction.  A concrete scale instantiation must still
prove that its lower-order masks are pullbacks of fibre flags, that every
order-27 mask meets each fibre in at most one sheet, and that every enumerated
owner key has score below the target.  No raw-bank completeness or numerical
table is asserted here.
-/

namespace LonelyRunner
namespace RamifiedCosetCover

open Finset

universe uOwner uGlobal uLocal uSheet uProvider uChoice uFibre

/-! ## Shared assignments project to owner-local covers -/

section OwnerProjection

variable {Owner : Type uOwner}
variable {Global : Type uGlobal}
variable {Local : Owner → Type uLocal}
variable {Sheet : Type uSheet} [Fintype Sheet] [DecidableEq Sheet]

/-- A local word covers every sheet in its owner chart. -/
def FullCoverGood
    (covered : (owner : Owner) → Local owner → Finset Sheet)
    (owner : Owner) (word : Local owner) : Prop :=
  covered owner word = Finset.univ

/-- One shared word covers every owner after applying the owner projection. -/
def SharedFullCover
    (project : (owner : Owner) → Global → Local owner)
    (covered : (owner : Owner) → Local owner → Finset Sheet) : Prop :=
  LRCPreNerveProjection.GloballyFeasible project (FullCoverGood covered)

/-- Every owner bank is nonempty when its word may vary with the owner. -/
def EveryOwnerLocallyCoverable
    (covered : (owner : Owner) → Local owner → Finset Sheet) : Prop :=
  LRCPreNerveProjection.OwnerLocallyFeasible (FullCoverGood covered)

/-- Predicate-preserving bridge from a shared assignment to all owner banks. -/
theorem everyOwnerLocallyCoverable_of_sharedFullCover
    (project : (owner : Owner) → Global → Local owner)
    (covered : (owner : Owner) → Local owner → Finset Sheet)
    (hshared : SharedFullCover project covered) :
    EveryOwnerLocallyCoverable covered :=
  LRCPreNerveProjection.ownerLocallyFeasible_of_globallyFeasible
    project (FullCoverGood covered) hshared

/-- A single owner with no local covering word obstructs every shared word. -/
theorem not_sharedFullCover_of_empty_owner
    (project : (owner : Owner) → Global → Local owner)
    (covered : (owner : Owner) → Local owner → Finset Sheet)
    (owner : Owner)
    (hempty : ∀ word, covered owner word ≠ Finset.univ) :
    ¬SharedFullCover project covered :=
  LRCPreNerveProjection.not_globallyFeasible_of_empty_owner
    project (FullCoverGood covered) owner hempty

/-- A uniform strict cardinality deficit at one owner is enough to obstruct
every shared assignment.  This is the direct consumer of the exact maximum-
union and sound upper-relaxation tables used by the common-scale certificates. -/
theorem not_sharedFullCover_of_owner_card_lt
    (project : (owner : Owner) → Global → Local owner)
    (covered : (owner : Owner) → Local owner → Finset Sheet)
    (owner : Owner)
    (hdeficit : ∀ word, (covered owner word).card < Fintype.card Sheet) :
    ¬SharedFullCover project covered := by
  apply not_sharedFullCover_of_empty_owner project covered owner
  intro word hcover
  have hlt := hdeficit word
  rw [hcover, Finset.card_univ] at hlt
  omega

end OwnerProjection

/-! ## Anchor union plus independent remainder caps -/

section AnchorBound

variable {Sheet : Type uSheet} [DecidableEq Sheet]
variable {Provider : Type uProvider} [DecidableEq Provider]

/-- Add a finite family of provider pieces to an already fixed anchor union. -/
def extendUnion (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet) : Finset Sheet :=
  base ∪ providers.biUnion piece

/-- Only the part of a new piece outside the anchor union can enlarge it. -/
theorem extendUnion_eq_remainders
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet) :
    extendUnion base providers piece =
      base ∪ providers.biUnion (fun i ↦ piece i \ base) := by
  ext sheet
  simp only [extendUnion, mem_union, mem_biUnion, mem_sdiff]
  constructor
  · rintro (hsheet | ⟨i, hi, hpiece⟩)
    · exact Or.inl hsheet
    · by_cases hbase : sheet ∈ base
      · exact Or.inl hbase
      · exact Or.inr ⟨i, hi, hpiece, hbase⟩
  · rintro (hsheet | ⟨i, hi, hpiece, _hbase⟩)
    · exact Or.inl hsheet
    · exact Or.inr ⟨i, hi, hpiece⟩

/-- Sound nested-fibre relaxation: each nonanchor provider may independently
use its best certified contribution outside the fixed anchor union. -/
theorem extendUnion_card_le_base_add_caps
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet) (cap : Provider → ℕ)
    (hcap : ∀ i ∈ providers, (piece i \ base).card ≤ cap i) :
    (extendUnion base providers piece).card ≤
      base.card + ∑ i ∈ providers, cap i := by
  rw [extendUnion_eq_remainders]
  calc
    (base ∪ providers.biUnion (fun i ↦ piece i \ base)).card ≤
        base.card + (providers.biUnion (fun i ↦ piece i \ base)).card :=
      Finset.card_union_le _ _
    _ ≤ base.card + ∑ i ∈ providers, (piece i \ base).card := by
      exact Nat.add_le_add_left Finset.card_biUnion_le _
    _ ≤ base.card + ∑ i ∈ providers, cap i := by
      exact Nat.add_le_add_left (Finset.sum_le_sum hcap) _

/-- Assignment form of the preceding bound.  The cap may be proved uniformly
over all choices, while the selected choices remain globally coupled. -/
theorem assignment_extendUnion_card_le
    {Choice : Provider → Type uChoice}
    (base : Finset Sheet) (providers : Finset Provider)
    (mask : (i : Provider) → Choice i → Finset Sheet)
    (word : (i : Provider) → Choice i) (cap : Provider → ℕ)
    (hcap : ∀ i ∈ providers, ∀ choice, (mask i choice \ base).card ≤ cap i) :
    (extendUnion base providers (fun i ↦ mask i (word i))).card ≤
      base.card + ∑ i ∈ providers, cap i := by
  exact extendUnion_card_le_base_add_caps base providers
    (fun i ↦ mask i (word i)) cap (fun i hi ↦ hcap i hi (word i))

/-- If the relaxed anchor-plus-cap total is smaller than a target, then the
literal selected union cannot cover that target. -/
theorem not_covers_of_base_add_caps_lt
    (target base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet) (cap : Provider → ℕ)
    (hcap : ∀ i ∈ providers, (piece i \ base).card ≤ cap i)
    (hdeficit : base.card + ∑ i ∈ providers, cap i < target.card) :
    ¬ target ⊆ extendUnion base providers piece := by
  intro hcover
  have hlower : target.card ≤ (extendUnion base providers piece).card :=
    Finset.card_le_card hcover
  have hupper :=
    extendUnion_card_le_base_add_caps base providers piece cap hcap
  omega

end AnchorBound

/-! ## Saturated fibres and transversal pieces -/

section FibreFlag

variable {Sheet : Type uSheet} [Fintype Sheet] [DecidableEq Sheet]
variable {Provider : Type uProvider} [DecidableEq Provider]
variable {Fibre : Type uFibre} [Fintype Fibre] [DecidableEq Fibre]

/-- The part of a sheet set lying over one quotient fibre. -/
def fibreSlice (fibre : Sheet → Fibre) (j : Fibre)
    (sheets : Finset Sheet) : Finset Sheet :=
  sheets.filter fun sheet ↦ fibre sheet = j

/-- A base is a union of complete fibres, with `complete` recording exactly
which fibres occur. -/
def IsCompleteFibreUnion (fibre : Sheet → Fibre)
    (complete : Finset Fibre) (base : Finset Sheet) : Prop :=
  ∀ sheet, sheet ∈ base ↔ fibre sheet ∈ complete

/-- Fibre-constant membership is exactly what is needed to turn a periodic
mask into a union of complete quotient fibres. -/
theorem completeFibreUnion_image_of_fibre_constant
    (fibre : Sheet → Fibre) (base : Finset Sheet)
    (hconstant : ∀ s t, fibre s = fibre t → (s ∈ base ↔ t ∈ base)) :
    IsCompleteFibreUnion fibre (base.image fibre) base := by
  intro sheet
  constructor
  · intro hsheet
    exact Finset.mem_image.mpr ⟨sheet, hsheet, rfl⟩
  · intro himage
    rcases Finset.mem_image.mp himage with ⟨source, hsource, heq⟩
    exact (hconstant source sheet heq).mp hsource

/-- Conversely, a recorded complete-fibre presentation makes membership
constant on every fibre. -/
theorem fibre_constant_of_completeFibreUnion
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    {s t : Sheet} (hst : fibre s = fibre t) :
    s ∈ base ↔ t ∈ base := by
  constructor
  · intro hs
    apply (hbase t).mpr
    rw [← hst]
    exact (hbase s).mp hs
  · intro ht
    apply (hbase s).mpr
    rw [hst]
    exact (hbase t).mp ht

/-- Number of selected transversal pieces which actually meet a fibre. -/
def fibreHitCount (fibre : Sheet → Fibre) (j : Fibre)
    (providers : Finset Provider) (piece : Provider → Finset Sheet) : ℕ :=
  (providers.filter fun i ↦ (fibreSlice fibre j (piece i)).Nonempty).card

/-- Per-fibre saturated flag.  Complete base fibres contribute their full
size.  On every other fibre, transversal contributions are capped both by
the fibre size and by the number of pieces which meet it. -/
def fibreFlagTerm (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (providers : Finset Provider) (piece : Provider → Finset Sheet)
    (j : Fibre) : ℕ :=
  if j ∈ complete then (fibreSlice fibre j Finset.univ).card
  else min (fibreSlice fibre j Finset.univ).card
    (fibreHitCount fibre j providers piece)

/-- Sum of the saturated flags over all quotient fibres. -/
def fibreFlagScore (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (providers : Finset Provider) (piece : Provider → Finset Sheet) : ℕ :=
  ∑ j : Fibre, fibreFlagTerm fibre complete providers piece j

theorem biUnion_eq_biUnion_nonempty
    (providers : Finset Provider) (piece : Provider → Finset Sheet) :
    providers.biUnion piece =
      (providers.filter fun i ↦ (piece i).Nonempty).biUnion piece := by
  ext sheet
  simp only [mem_biUnion, mem_filter]
  constructor
  · rintro ⟨i, hi, hsheet⟩
    exact ⟨i, ⟨hi, ⟨sheet, hsheet⟩⟩, hsheet⟩
  · rintro ⟨i, ⟨hi, -⟩, hsheet⟩
    exact ⟨i, hi, hsheet⟩

/-- A union of pieces of size at most one has cardinality at most the number
of nonempty pieces. -/
theorem card_biUnion_le_nonempty_count
    (providers : Finset Provider) (piece : Provider → Finset Sheet)
    (hone : ∀ i ∈ providers, (piece i).card ≤ 1) :
    (providers.biUnion piece).card ≤
      (providers.filter fun i ↦ (piece i).Nonempty).card := by
  rw [biUnion_eq_biUnion_nonempty]
  calc
    ((providers.filter fun i ↦ (piece i).Nonempty).biUnion piece).card ≤
        ∑ i ∈ providers.filter (fun i ↦ (piece i).Nonempty),
          (piece i).card := Finset.card_biUnion_le
    _ ≤ ∑ _i ∈ providers.filter (fun i ↦ (piece i).Nonempty), 1 := by
      apply Finset.sum_le_sum
      intro i hi
      exact hone i (Finset.mem_of_mem_filter i hi)
    _ = (providers.filter fun i ↦ (piece i).Nonempty).card := by simp

theorem fibreSlice_extendUnion_of_incomplete
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (j : Fibre) (hj : j ∉ complete) :
    fibreSlice fibre j (extendUnion base providers piece) =
      providers.biUnion fun i ↦ fibreSlice fibre j (piece i) := by
  ext sheet
  simp only [fibreSlice, extendUnion, mem_filter, mem_union, mem_biUnion]
  constructor
  · rintro ⟨hbaseOrPiece, hfibre⟩
    rcases hbaseOrPiece with hsheetBase | ⟨i, hi, hsheet⟩
    · have : j ∈ complete := by
        rw [← hfibre]
        exact (hbase sheet).mp hsheetBase
      exact (hj this).elim
    · exact ⟨i, hi, hsheet, hfibre⟩
  · rintro ⟨i, hi, hsheet, hfibre⟩
    exact ⟨Or.inr ⟨i, hi, hsheet⟩, hfibre⟩

theorem fibreSlice_card_le_univ
    (fibre : Sheet → Fibre) (j : Fibre) (sheets : Finset Sheet) :
    (fibreSlice fibre j sheets).card ≤
      (fibreSlice fibre j Finset.univ).card := by
  apply Finset.card_le_card
  intro sheet hsheet
  rw [fibreSlice, Finset.mem_filter] at hsheet ⊢
  exact ⟨Finset.mem_univ _, hsheet.2⟩

/-- Outside a complete base fibre, selected transversal pieces cover at most
one sheet for every piece which meets that fibre. -/
theorem incomplete_fibre_card_le_hitCount
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (htransversal : ∀ i ∈ providers, ∀ j,
      (fibreSlice fibre j (piece i)).card ≤ 1)
    (j : Fibre) (hj : j ∉ complete) :
    (fibreSlice fibre j (extendUnion base providers piece)).card ≤
      fibreHitCount fibre j providers piece := by
  rw [fibreSlice_extendUnion_of_incomplete fibre complete base providers piece
    hbase j hj]
  exact card_biUnion_le_nonempty_count providers
    (fun i ↦ fibreSlice fibre j (piece i)) (fun i hi ↦ htransversal i hi j)

theorem covered_fibre_card_le_flagTerm
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (htransversal : ∀ i ∈ providers, ∀ j,
      (fibreSlice fibre j (piece i)).card ≤ 1)
    (j : Fibre) :
    (fibreSlice fibre j (extendUnion base providers piece)).card ≤
      fibreFlagTerm fibre complete providers piece j := by
  by_cases hj : j ∈ complete
  · rw [fibreFlagTerm, if_pos hj]
    exact fibreSlice_card_le_univ fibre j _
  · rw [fibreFlagTerm, if_neg hj]
    exact le_min (fibreSlice_card_le_univ fibre j _)
      (incomplete_fibre_card_le_hitCount fibre complete base providers piece
        hbase htransversal j hj)

/-- The saturated fibre score is a sound upper bound for every literal union
represented by the complete-fibre base and transversal pieces. -/
theorem extendUnion_card_le_fibreFlagScore
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (htransversal : ∀ i ∈ providers, ∀ j,
      (fibreSlice fibre j (piece i)).card ≤ 1) :
    (extendUnion base providers piece).card ≤
      fibreFlagScore fibre complete providers piece := by
  have hsplit :
      (extendUnion base providers piece).card =
        ∑ j : Fibre,
          (fibreSlice fibre j (extendUnion base providers piece)).card := by
    simpa [fibreSlice] using
      (Finset.card_eq_sum_card_fiberwise
        (f := fibre) (s := extendUnion base providers piece)
        (t := (Finset.univ : Finset Fibre)) (by simp))
  rw [hsplit, fibreFlagScore]
  exact Finset.sum_le_sum fun j _ ↦
    covered_fibre_card_le_flagTerm fibre complete base providers piece
      hbase htransversal j

/-- A saturated flag never exceeds the number of literal sheets. -/
theorem fibreFlagScore_le_card
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (providers : Finset Provider) (piece : Provider → Finset Sheet) :
    fibreFlagScore fibre complete providers piece ≤ Fintype.card Sheet := by
  have hsplit :
      Fintype.card Sheet =
        ∑ j : Fibre, (fibreSlice fibre j Finset.univ).card := by
    simpa [fibreSlice] using
      (Finset.card_eq_sum_card_fiberwise
        (f := fibre) (s := (Finset.univ : Finset Sheet))
        (t := (Finset.univ : Finset Fibre)) (by simp))
  rw [fibreFlagScore, hsplit]
  apply Finset.sum_le_sum
  intro j _
  by_cases hj : j ∈ complete
  · simp [fibreFlagTerm, hj]
  · simp only [fibreFlagTerm, if_neg hj]
    exact min_le_left _ _

/-- Precise c27-facing implication: a literal full cover forces the saturated
flag score to be exactly the number of sheets. -/
theorem fibreFlagScore_eq_card_of_fullCover
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (htransversal : ∀ i ∈ providers, ∀ j,
      (fibreSlice fibre j (piece i)).card ≤ 1)
    (hcover : extendUnion base providers piece = Finset.univ) :
    fibreFlagScore fibre complete providers piece = Fintype.card Sheet := by
  apply le_antisymm
  · exact fibreFlagScore_le_card fibre complete providers piece
  · have hbound := extendUnion_card_le_fibreFlagScore fibre complete
      base providers piece hbase htransversal
    rw [hcover, Finset.card_univ] at hbound
    exact hbound

/-- Terminal no-cover consumer: any strict saturated-fibre deficit rules out
the literal cover. -/
theorem not_fullCover_of_fibreFlagScore_lt
    (fibre : Sheet → Fibre) (complete : Finset Fibre)
    (base : Finset Sheet) (providers : Finset Provider)
    (piece : Provider → Finset Sheet)
    (hbase : IsCompleteFibreUnion fibre complete base)
    (htransversal : ∀ i ∈ providers, ∀ j,
      (fibreSlice fibre j (piece i)).card ≤ 1)
    (hdeficit : fibreFlagScore fibre complete providers piece < Fintype.card Sheet) :
    extendUnion base providers piece ≠ Finset.univ := by
  intro hcover
  have heq := fibreFlagScore_eq_card_of_fullCover fibre complete base
    providers piece hbase htransversal hcover
  omega

#print axioms everyOwnerLocallyCoverable_of_sharedFullCover
#print axioms not_sharedFullCover_of_empty_owner
#print axioms not_sharedFullCover_of_owner_card_lt
#print axioms assignment_extendUnion_card_le
#print axioms not_covers_of_base_add_caps_lt
#print axioms completeFibreUnion_image_of_fibre_constant
#print axioms fibre_constant_of_completeFibreUnion
#print axioms extendUnion_card_le_fibreFlagScore
#print axioms fibreFlagScore_eq_card_of_fullCover
#print axioms not_fullCover_of_fibreFlagScore_lt

end FibreFlag

end RamifiedCosetCover
end LonelyRunner
