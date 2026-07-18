/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-c12-lean (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# The scale-twelve owner-orthogonality obstruction

This file formalizes the terminal symbolic quotient in THM-976.  After the
exact hereditary-order, scalar-capacity, and owner-local reductions, every
remaining support chooses one sign above each of the six projective classes
in `F_13^*/{+-1}`, and every provider has effective order twelve.  Thus a
support is a six-bit sign word and a global unit assignment is a six-letter
word in the units `1,5,7,11` modulo twelve.

Every literal order-twelve provider mask has two sheets.  Hence full coverage
at one owner, by six providers on twelve sheets, forces the six masks to be
pairwise disjoint.  A kernel-checked two-provider core then proves that owner
obligations zero and four are disjoint.  Common sign negation preserves all
ratios; eight tables of only `4 * 4^2 = 64` rows cover the lower 32 sign
sections, and negation supplies the upper 32.  This one disjoint owner pair is
already enough to make the sixfold owner intersection empty.

The coordinates are indexed by projective class rather than by increasing
selected label; this is an exact reindexing of every fibre.  `sheetMask` is
the frozen literal twelve-sheet CRT table from THM-976.  This module does
**not** formalize the preceding reduction of 2,413,458,432 raw labelled
contexts to the 64 sign-transversal fibres, nor the reported exact owner size
48 or the other fourteen pairwise-disjoint owner pairs.  There is no `sorry`,
no `native_decide`, and no added axiom.
-/

namespace LonelyRunner
namespace ScaleTwelveOwnerOrthogonality

-- Keep all declaration elaboration and kernel checking sequential.  With the
-- Lean 4.30 command-line default, even the small independent finite lemmas in
-- this module would otherwise be retained concurrently.
set_option Elab.async false

abbrev ProjectiveClass := Fin 6
abbrev SignWord := Fin 64
abbrev UnitDigit := Fin 4
abbrev UnitWord := ProjectiveClass → UnitDigit

def representative (i : ProjectiveClass) : Nat :=
  ![1, 2, 3, 4, 5, 6] i

def signAt (a : SignWord) (i : ProjectiveClass) : Bool :=
  (a.val / 2 ^ i.val) % 2 == 1

def signedLabel (a : SignWord) (i : ProjectiveClass) : Nat :=
  if signAt a i then 13 - representative i else representative i

def inverse13 : Nat → Nat
  | 1 => 1
  | 2 => 7
  | 3 => 9
  | 4 => 10
  | 5 => 8
  | 6 => 11
  | 7 => 2
  | 8 => 5
  | 9 => 3
  | 10 => 4
  | 11 => 6
  | 12 => 12
  | _ => 0

def ratio (a : SignWord) (provider owner : ProjectiveClass) : Fin 13 :=
  ⟨signedLabel a provider * inverse13 (signedLabel a owner) % 13,
    Nat.mod_lt _ (by omega)⟩

/-- Common negation of all selected labels. -/
def globalNeg (a : SignWord) : SignWord :=
  ⟨63 - a.val, by omega⟩

abbrev HalfSignWord := Fin 32

def halfSign (h : HalfSignWord) : SignWord :=
  ⟨h.val, by omega⟩

/-- Every sign section is in the lower half or is the common negative of one
in the lower half. -/
theorem sign_half_cases (a : SignWord) :
    ∃ h : HalfSignWord, a = halfSign h ∨ a = globalNeg (halfSign h) := by
  revert a
  decide

/-- Common sign negation cancels from every provider/owner ratio. -/
theorem ratio_globalNeg :
    ∀ (a : SignWord) (provider owner : ProjectiveClass),
      ratio (globalNeg a) provider owner = ratio a provider owner := by
  decide

/-- Literal order-twelve sheet masks at normalized owner one.  Rows are
ratios `1,...,12`; columns are the units `1,5,7,11`. -/
def sheetMask (q : Nat) : UnitDigit → Nat :=
  match q with
  | 1 => ![0xc00, 0x840, 0x810, 0x801]
  | 2 => ![0x600, 0x042, 0x210, 0x003]
  | 3 => ![0x300, 0x102, 0x204, 0x006]
  | 4 => ![0x180, 0x108, 0x084, 0x00c]
  | 5 => ![0x0c0, 0x408, 0x081, 0x018]
  | 6 => ![0x060, 0x420, 0x021, 0x030]
  | 7 => ![0x030, 0x021, 0x420, 0x060]
  | 8 => ![0x018, 0x081, 0x408, 0x0c0]
  | 9 => ![0x00c, 0x084, 0x108, 0x180]
  | 10 => ![0x006, 0x204, 0x102, 0x300]
  | 11 => ![0x003, 0x210, 0x042, 0x600]
  | 12 => ![0x001, 0x010, 0x040, 0x400]
  | _ => ![0, 0, 0, 0]

def providerMask (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) : Nat :=
  sheetMask (ratio a provider owner).val (x provider)

def maskSheets (m : Nat) : Finset (Fin 12) :=
  Finset.univ.filter fun ell ↦ m.testBit ell.val

def providerSheets (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) : Finset (Fin 12) :=
  maskSheets (providerMask a x provider owner)

def ownerSheets (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) : Finset (Fin 12) :=
  Finset.univ.biUnion fun provider ↦ providerSheets a x provider owner

def ownerSat (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) : Prop :=
  ownerSheets a x owner = Finset.univ

theorem ratio_ne_zero :
    ∀ (a : SignWord) (provider owner : ProjectiveClass),
      ratio a provider owner ≠ 0 := by
  decide

/-- Ratio `-1` would require two selected labels in the same projective
class.  A sign transversal has only one, including when provider=owner. -/
theorem ratio_ne_twelve :
    ∀ (a : SignWord) (provider owner : ProjectiveClass),
      ratio a provider owner ≠ 12 := by
  decide

/-- Every admissible ratio row has two sheets.  The excluded `q=12` row is
the one-sheet negative-label row and cannot occur in a sign transversal. -/
theorem maskSheets_card_admissible :
    ∀ (q : Fin 13) (d : UnitDigit), q ≠ 0 → q ≠ 12 →
      (maskSheets (sheetMask q.val d)).card = 2 := by
  decide

theorem providerSheets_card (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) :
    (providerSheets a x provider owner).card = 2 := by
  exact maskSheets_card_admissible (ratio a provider owner) (x provider)
    (ratio_ne_zero a provider owner) (ratio_ne_twelve a provider owner)

theorem provider_capacity_sum (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) :
    ∑ provider : ProjectiveClass,
      (providerSheets a x provider owner).card = 12 := by
  calc
    ∑ provider : ProjectiveClass,
        (providerSheets a x provider owner).card =
        ∑ _provider : ProjectiveClass, 2 := by
          apply Finset.sum_congr rfl
          intro provider _
          exact providerSheets_card a x provider owner
    _ = 12 := by decide

/-- A full union with exactly twelve total incidences is a partition. -/
theorem pairwiseDisjoint_of_full_capacity
    (t : ProjectiveClass → Finset (Fin 12))
    (hfull : Finset.univ.biUnion t = Finset.univ)
    (hcard : ∑ i : ProjectiveClass, (t i).card = 12) :
    (Set.univ : Set ProjectiveClass).PairwiseDisjoint t := by
  let Incidence := Σ i : ProjectiveClass, {s // s ∈ t i}
  let sheet : Incidence → Fin 12 := fun z ↦ z.2.1
  have hsurj : Function.Surjective sheet := by
    intro s
    have hs : s ∈ Finset.univ.biUnion t := by simpa [hfull]
    rw [Finset.mem_biUnion] at hs
    obtain ⟨i, -, hsi⟩ := hs
    exact ⟨⟨i, ⟨s, hsi⟩⟩, rfl⟩
  have hincidence : Fintype.card Incidence = Fintype.card (Fin 12) := by
    simp only [Incidence, Fintype.card_sigma, Fintype.card_coe]
    simpa using hcard
  have hinj : Function.Injective sheet :=
    ((Fintype.bijective_iff_surjective_and_card sheet).2
      ⟨hsurj, hincidence⟩).1
  rintro i - j - hij
  change Disjoint (t i) (t j)
  rw [Finset.disjoint_left]
  intro s hsi hsj
  have heq : (⟨i, ⟨s, hsi⟩⟩ : Incidence) = ⟨j, ⟨s, hsj⟩⟩ :=
    hinj rfl
  exact hij (congrArg Sigma.fst heq)

theorem ownerSat_pairwise (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) (howner : ownerSat a x owner) :
    (Set.univ : Set ProjectiveClass).PairwiseDisjoint
      (fun provider ↦ providerSheets a x provider owner) := by
  exact pairwiseDisjoint_of_full_capacity
    (fun provider ↦ providerSheets a x provider owner)
    howner (provider_capacity_sum a x owner)

def maskDisjoint (m n : Nat) : Bool :=
  (m &&& n) == 0

abbrev literalDisjointStatement (q r : Fin 13) (d e : UnitDigit) : Prop :=
  Disjoint (maskSheets (sheetMask q.val d))
      (maskSheets (sheetMask r.val e)) ↔
    maskDisjoint (sheetMask q.val d) (sheetMask r.val e) = true

/-- Split the `13^2 * 4^2` representation audit by its first ratio so no one
kernel declaration retains the full finite proof term. -/
theorem literalDisjoint_0 : ∀ r d e, literalDisjointStatement 0 r d e := by decide
theorem literalDisjoint_1 : ∀ r d e, literalDisjointStatement 1 r d e := by decide
theorem literalDisjoint_2 : ∀ r d e, literalDisjointStatement 2 r d e := by decide
theorem literalDisjoint_3 : ∀ r d e, literalDisjointStatement 3 r d e := by decide
theorem literalDisjoint_4 : ∀ r d e, literalDisjointStatement 4 r d e := by decide
theorem literalDisjoint_5 : ∀ r d e, literalDisjointStatement 5 r d e := by decide
theorem literalDisjoint_6 : ∀ r d e, literalDisjointStatement 6 r d e := by decide
theorem literalDisjoint_7 : ∀ r d e, literalDisjointStatement 7 r d e := by decide
theorem literalDisjoint_8 : ∀ r d e, literalDisjointStatement 8 r d e := by decide
theorem literalDisjoint_9 : ∀ r d e, literalDisjointStatement 9 r d e := by decide
theorem literalDisjoint_10 : ∀ r d e, literalDisjointStatement 10 r d e := by decide
theorem literalDisjoint_11 : ∀ r d e, literalDisjointStatement 11 r d e := by decide
theorem literalDisjoint_12 : ∀ r d e, literalDisjointStatement 12 r d e := by decide

theorem literalMask_disjoint_iff :
    ∀ (q r : Fin 13) (d e : UnitDigit), literalDisjointStatement q r d e := by
  intro q
  fin_cases q
  · exact literalDisjoint_0
  · exact literalDisjoint_1
  · exact literalDisjoint_2
  · exact literalDisjoint_3
  · exact literalDisjoint_4
  · exact literalDisjoint_5
  · exact literalDisjoint_6
  · exact literalDisjoint_7
  · exact literalDisjoint_8
  · exact literalDisjoint_9
  · exact literalDisjoint_10
  · exact literalDisjoint_11
  · exact literalDisjoint_12

theorem providerPairDisjoint_of_ownerSat
    (a : SignWord) (x : UnitWord) (owner provider₁ provider₂ : ProjectiveClass)
    (hneq : provider₁ ≠ provider₂) (howner : ownerSat a x owner) :
    maskDisjoint (providerMask a x provider₁ owner)
      (providerMask a x provider₂ owner) = true := by
  have hp := ownerSat_pairwise a x owner howner
  have hdisjoint :
      Disjoint (providerSheets a x provider₁ owner)
        (providerSheets a x provider₂ owner) :=
    hp (by simp) (by simp) hneq
  rw [providerSheets, providerMask] at hdisjoint
  exact (literalMask_disjoint_iff
    (ratio a provider₁ owner) (ratio a provider₂ owner)
    (x provider₁) (x provider₂)).1 hdisjoint

/-! ## Small two-owner provider cores -/

def projectedUnits {k : Nat} (x : UnitWord)
    (providers : Fin k → ProjectiveClass) : Fin k → UnitDigit :=
  fun i ↦ x (providers i)

def coreProviderMask {k : Nat} (a : SignWord) (u : Fin k → UnitDigit)
    (providers : Fin k → ProjectiveClass) (i : Fin k)
    (owner : ProjectiveClass) : Nat :=
  sheetMask (ratio a (providers i) owner).val (u i)

theorem coreProviderMask_projected {k : Nat} (a : SignWord) (x : UnitWord)
    (providers : Fin k → ProjectiveClass) (i : Fin k)
    (owner : ProjectiveClass) :
    coreProviderMask a (projectedUnits x providers) providers i owner =
      providerMask a x (providers i) owner := by
  rfl

def corePartition2 (a : SignWord) (u : Fin 2 → UnitDigit)
    (providers : Fin 2 → ProjectiveClass) (owner : ProjectiveClass) : Bool :=
  maskDisjoint (coreProviderMask a u providers 0 owner)
    (coreProviderMask a u providers 1 owner)

def corePair2 (a : SignWord) (u : Fin 2 → UnitDigit)
    (providers : Fin 2 → ProjectiveClass) (left right : ProjectiveClass) : Bool :=
  corePartition2 a u providers left && corePartition2 a u providers right

theorem coreProviderMask_globalNeg {k : Nat} (a : SignWord)
    (u : Fin k → UnitDigit) (providers : Fin k → ProjectiveClass)
    (i : Fin k) (owner : ProjectiveClass) :
    coreProviderMask (globalNeg a) u providers i owner =
      coreProviderMask a u providers i owner := by
  rw [coreProviderMask, coreProviderMask, ratio_globalNeg]

theorem corePair2_globalNeg (a : SignWord) (u : Fin 2 → UnitDigit)
    (providers : Fin 2 → ProjectiveClass) (left right : ProjectiveClass) :
    corePair2 (globalNeg a) u providers left right =
      corePair2 a u providers left right := by
  simp only [corePair2, corePartition2, coreProviderMask_globalNeg]

theorem corePair2_false_of_half
    (providers : Fin 2 → ProjectiveClass) (left right : ProjectiveClass)
    (hhalf : ∀ (h : HalfSignWord) (u : Fin 2 → UnitDigit),
      corePair2 (halfSign h) u providers left right = false) :
    ∀ (a : SignWord) (u : Fin 2 → UnitDigit),
      corePair2 a u providers left right = false := by
  intro a u
  obtain ⟨h, rfl | rfl⟩ := sign_half_cases a
  · exact hhalf h u
  · rw [corePair2_globalNeg]
    exact hhalf h u

theorem corePartition2_of_ownerSat (a : SignWord) (x : UnitWord)
    (providers : Fin 2 → ProjectiveClass) (hinj : Function.Injective providers)
    (owner : ProjectiveClass) (howner : ownerSat a x owner) :
    corePartition2 a (projectedUnits x providers) providers owner = true := by
  rw [corePartition2, coreProviderMask_projected, coreProviderMask_projected]
  exact providerPairDisjoint_of_ownerSat a x owner (providers 0) (providers 1)
    (hinj.ne (by decide)) howner

theorem ownerPair_false_of_core2
    (providers : Fin 2 → ProjectiveClass) (hinj : Function.Injective providers)
    (left right : ProjectiveClass)
    (hcore : ∀ (a : SignWord) (u : Fin 2 → UnitDigit),
      corePair2 a u providers left right = false)
    (a : SignWord) (x : UnitWord)
    (hleft : ownerSat a x left) (hright : ownerSat a x right) : False := by
  have htrue : corePair2 a (projectedUnits x providers) providers left right = true := by
    simp only [corePair2, Bool.and_eq_true]
    exact ⟨corePartition2_of_ownerSat a x providers hinj left hleft,
      corePartition2_of_ownerSat a x providers hinj right hright⟩
  rw [hcore a (projectedUnits x providers)] at htrue
  contradiction

/-! ## A 512-row exact terminal core -/

def providers04 : Fin 2 → ProjectiveClass := ![1, 2]

abbrev HalfBlock := Fin 2
abbrev HalfOffset := Fin 16
abbrev OffsetBlock := Fin 4
abbrev OffsetDigit := Fin 4

def blockedHalfSign (block : HalfBlock) (offset : HalfOffset) : HalfSignWord :=
  ⟨16 * block.val + offset.val, by omega⟩

def halfBlock (h : HalfSignWord) : HalfBlock :=
  ⟨h.val / 16, by omega⟩

def halfOffset (h : HalfSignWord) : HalfOffset :=
  ⟨h.val % 16, Nat.mod_lt _ (by omega)⟩

theorem blockedHalfSign_reconstruct (h : HalfSignWord) :
    blockedHalfSign (halfBlock h) (halfOffset h) = h := by
  apply Fin.ext
  simp [blockedHalfSign, halfBlock, halfOffset]
  omega

def blockedOffset (block : OffsetBlock) (digit : OffsetDigit) : HalfOffset :=
  ⟨4 * block.val + digit.val, by omega⟩

def offsetBlock (offset : HalfOffset) : OffsetBlock :=
  ⟨offset.val / 4, by omega⟩

def offsetDigit (offset : HalfOffset) : OffsetDigit :=
  ⟨offset.val % 4, Nat.mod_lt _ (by omega)⟩

theorem blockedOffset_reconstruct (offset : HalfOffset) :
    blockedOffset (offsetBlock offset) (offsetDigit offset) = offset := by
  apply Fin.ext
  simp [blockedOffset, offsetBlock, offsetDigit]
  omega

-- Source-level serialization prevents Lean 4.30 from retaining simultaneous
-- kernel checks.  Each literal table has only `4 * 4^2 = 64` rows.
set_option maxRecDepth 100000
set_option maxHeartbeats 500000

theorem core04_b0_o0_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 0 (blockedOffset 0 d)))
      u providers04 0 4 = false := by decide
theorem core04_b0_o1_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 0 (blockedOffset 1 d)))
      u providers04 0 4 = false := by decide
theorem core04_b0_o2_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 0 (blockedOffset 2 d)))
      u providers04 0 4 = false := by decide
theorem core04_b0_o3_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 0 (blockedOffset 3 d)))
      u providers04 0 4 = false := by decide
theorem core04_b1_o0_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 1 (blockedOffset 0 d)))
      u providers04 0 4 = false := by decide
theorem core04_b1_o1_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 1 (blockedOffset 1 d)))
      u providers04 0 4 = false := by decide
theorem core04_b1_o2_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 1 (blockedOffset 2 d)))
      u providers04 0 4 = false := by decide
theorem core04_b1_o3_false : ∀ d u,
    corePair2 (halfSign (blockedHalfSign 1 (blockedOffset 3 d)))
      u providers04 0 4 = false := by decide

theorem core04_tiny_false (block : HalfBlock) (oblock : OffsetBlock)
    (d : OffsetDigit) (u : Fin 2 → UnitDigit) :
    corePair2 (halfSign (blockedHalfSign block (blockedOffset oblock d)))
      u providers04 0 4 = false := by
  fin_cases block <;> fin_cases oblock
  · exact core04_b0_o0_false d u
  · exact core04_b0_o1_false d u
  · exact core04_b0_o2_false d u
  · exact core04_b0_o3_false d u
  · exact core04_b1_o0_false d u
  · exact core04_b1_o1_false d u
  · exact core04_b1_o2_false d u
  · exact core04_b1_o3_false d u

theorem core04_block_false (block : HalfBlock) (offset : HalfOffset)
    (u : Fin 2 → UnitDigit) :
    corePair2 (halfSign (blockedHalfSign block offset)) u providers04 0 4 = false := by
  have htiny := core04_tiny_false block (offsetBlock offset) (offsetDigit offset) u
  rw [blockedOffset_reconstruct] at htiny
  exact htiny

theorem core04_half_false :
    ∀ h u, corePair2 (halfSign h) u providers04 0 4 = false := by
  intro h u
  have hblock := core04_block_false (halfBlock h) (halfOffset h) u
  rw [blockedHalfSign_reconstruct] at hblock
  exact hblock

theorem core04_false :
    ∀ a u, corePair2 a u providers04 0 4 = false :=
  corePair2_false_of_half providers04 0 4 core04_half_false

theorem ownerPair04_false (a : SignWord) (x : UnitWord)
    (h0 : ownerSat a x 0) (h4 : ownerSat a x 4) : False :=
  ownerPair_false_of_core2 providers04 (by decide) 0 4 core04_false a x h0 h4

theorem ownerPair04_sets_disjoint (a : SignWord) :
    Disjoint {x : UnitWord | ownerSat a x 0}
      {x : UnitWord | ownerSat a x 4} := by
  rw [Set.disjoint_left]
  intro x h0 h4
  exact ownerPair04_false a x h0 h4

theorem no_common_unit_word (a : SignWord) (x : UnitWord)
    (hall : ∀ i : ProjectiveClass, ownerSat a x i) : False := by
  exact ownerPair04_false a x (hall 0) (hall 4)

#print axioms core04_false
#print axioms ownerPair04_false
#print axioms ownerPair04_sets_disjoint
#print axioms no_common_unit_word

end ScaleTwelveOwnerOrthogonality
end LonelyRunner
