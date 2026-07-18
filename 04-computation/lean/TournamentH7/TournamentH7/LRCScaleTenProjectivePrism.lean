/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-c10-lean (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# The scale-ten projective-prism obstruction

This file formalizes the terminal symbolic quotient in THM-970.  After the
exact order and owner-local reductions, a support is a sign section of the six
projective classes in `F_13^*/{+-1}`, every effective order is ten, and a
global unit word has six letters in `{1,3,7,9}`.

The compact predicates below replay only that final `2^6 * 4^6` fibre.  The
literal one- or two-sheet masks are the exact owner-one CRT masks from THM-970,
equations (8)--(10), transported by provider/owner ratios.  Ordinary kernel
`decide` verifies that owner obligations joined by multiplication by two on
the projective quotient are disjoint.  Multiplication by two is the six-cycle

`[1] -> [2] -> [4] -> [5] -> [3] -> [6] -> [1]`.

Thus the complementary zero-intersection graph is `C_6`; equivalently, the
possible pair-intersection graph is contained in `K_6 \ C_6`, the triangular
prism.  A single adjacent pair already excludes a common six-owner unit word.

This module does **not** formalize the preceding 821,620,800-context reduction
to the sign-section fibre, nor does it claim the positive intersection size
four on nonadjacent pairs.  There is no `sorry`, no `native_decide`, and no
added axiom.
-/

namespace LonelyRunner
namespace ScaleTenProjectivePrism

/- The command-line frontend otherwise checks independent theorem declarations
asynchronously.  Serialize this audit so every finite kernel proof is released
before the next starts.  This is a scheduling option only; it does not alter
the trusted kernel or the proof terms. -/
set_option Elab.async false

/-- The six projective classes, indexed by their representatives `1,...,6`. -/
abbrev ProjectiveClass := Fin 6

/-- A choice of one sign above each projective class. -/
abbrev SignWord := Fin 64

/-- One of the four units `1,3,7,9` modulo ten. -/
abbrev UnitDigit := Fin 4

/-- A six-letter word in the four effective order-ten units. -/
abbrev UnitWord := Fin 4096

/-- Positive representative in `{1,...,6}` of a projective class. -/
def representative (i : ProjectiveClass) : Nat :=
  ![1, 2, 3, 4, 5, 6] i

/-- Read a sign from the compact six-bit encoding. -/
def signAt (a : SignWord) (i : ProjectiveClass) : Bool :=
  (a.val / 2 ^ i.val) % 2 == 1

/-- The selected nonzero residue modulo thirteen, represented in `1,...,12`. -/
def signedLabel (a : SignWord) (i : ProjectiveClass) : Nat :=
  if signAt a i then 13 - representative i else representative i

/-- Read a base-four digit from the compact unit-word encoding. -/
def unitAt (x : UnitWord) (i : ProjectiveClass) : UnitDigit :=
  ⟨(x.val / 4 ^ i.val) % 4, by omega⟩

/-- Multiplicative inverse modulo thirteen on representatives `1,...,12`. -/
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

/-- Provider/owner ratio in `F_13^*`, represented in `Fin 13`. -/
def ratio (a : SignWord) (provider owner : ProjectiveClass) : Fin 13 :=
  ⟨signedLabel a provider * inverse13 (signedLabel a owner) % 13,
    Nat.mod_lt _ (by omega)⟩

/-- Exact order-ten sheet mask for a normalized provider ratio and one of the
four units `1,3,7,9`.  Bit `ell` records coverage of sheet `ell`.

These are the literal masks reconstructed from the CRT formula in THM-970,
not merely their cardinalities. -/
def sheetMask (q : Nat) : UnitDigit → Nat :=
  match q with
  | 1 => ![12, 72, 9, 24]
  | 2 => ![6, 576, 129, 48]
  | 3 => ![3, 516, 144, 96]
  | 4 => ![1, 4, 16, 64]
  | 5 => ![512, 32, 2, 128]
  | 6 => ![768, 288, 258, 384]
  | 7 => ![384, 258, 288, 768]
  | 8 => ![128, 2, 32, 512]
  | 9 => ![64, 16, 4, 1]
  | 10 => ![96, 144, 516, 3]
  | 11 => ![48, 129, 576, 6]
  | 12 => ![16, 1, 64, 4]
  | _ => ![0, 0, 0, 0]

/-- The literal sheet mask contributed by one provider at one owner. -/
def providerMask (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) : Nat :=
  sheetMask (ratio a provider owner).val (unitAt x provider)

/-- The finite set encoded by a ten-sheet bit mask. -/
def maskSheets (m : Nat) : Finset (Fin 10) :=
  Finset.univ.filter fun ell ↦ m.testBit ell.val

/-- The sheet set contributed by one provider at one owner. -/
def providerSheets (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) : Finset (Fin 10) :=
  maskSheets (providerMask a x provider owner)

/-- Union of all six provider sheet sets at an owner. -/
def ownerSheets (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) : Finset (Fin 10) :=
  Finset.univ.biUnion fun provider ↦ providerSheets a x provider owner

/-- The owner obligation: all ten sheets are covered at this owner. -/
def ownerSat (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) : Prop :=
  ownerSheets a x owner = Finset.univ

/-- Cardinality of an order-ten mask, as a function of its normalized ratio.
The four projective ratio classes `[1],[2],[3],[6]` have two sheets and
`[4],[5]` have one. -/
def maskCapacity (q : Fin 13) : Nat :=
  ![0, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1] q

/-- The frozen literal mask table has the advertised ratio-dependent
cardinality, independently of its unit digit. -/
theorem maskSheets_card :
    ∀ (q : Fin 13) (d : UnitDigit),
      (maskSheets (sheetMask q.val d)).card = maskCapacity q := by
  decide

/-- On every sign section and at every owner, the six mask cardinalities add
to exactly ten.  This is the exact-capacity feature that turns coverage into
a partition. -/
theorem ratio_capacity_sum :
    ∀ (a : SignWord) (owner : ProjectiveClass),
      ∑ provider : ProjectiveClass, maskCapacity (ratio a provider owner) = 10 := by
  decide

theorem providerSheets_card (a : SignWord) (x : UnitWord)
    (provider owner : ProjectiveClass) :
    (providerSheets a x provider owner).card =
      maskCapacity (ratio a provider owner) := by
  exact maskSheets_card (ratio a provider owner) (unitAt x provider)

theorem provider_capacity_sum (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) :
    ∑ provider : ProjectiveClass, (providerSheets a x provider owner).card = 10 := by
  calc
    ∑ provider : ProjectiveClass, (providerSheets a x provider owner).card =
        ∑ provider : ProjectiveClass, maskCapacity (ratio a provider owner) := by
          apply Finset.sum_congr rfl
          intro provider _
          exact providerSheets_card a x provider owner
    _ = 10 := ratio_capacity_sum a owner

/-- A full union with no excess total cardinality is a partition.  The proof
uses the incidence type `Σ i, {s // s ∈ t i}`: its projection onto sheets is
surjective and has equal finite cardinality, hence is injective. -/
theorem pairwiseDisjoint_of_full_capacity
    (t : ProjectiveClass → Finset (Fin 10))
    (hfull : Finset.univ.biUnion t = Finset.univ)
    (hcard : ∑ i : ProjectiveClass, (t i).card = 10) :
    (Set.univ : Set ProjectiveClass).PairwiseDisjoint t := by
  let Incidence := Σ i : ProjectiveClass, {s // s ∈ t i}
  let sheet : Incidence → Fin 10 := fun z ↦ z.2.1
  have hsurj : Function.Surjective sheet := by
    intro s
    have hs : s ∈ Finset.univ.biUnion t := by simpa [hfull]
    rw [Finset.mem_biUnion] at hs
    obtain ⟨i, -, hsi⟩ := hs
    exact ⟨⟨i, ⟨s, hsi⟩⟩, rfl⟩
  have hincidence : Fintype.card Incidence = Fintype.card (Fin 10) := by
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

/-- Therefore literal full coverage at an owner forces all six provider masks
to be pairwise disjoint. -/
theorem ownerSat_pairwise (a : SignWord) (x : UnitWord)
    (owner : ProjectiveClass) (howner : ownerSat a x owner) :
    (Set.univ : Set ProjectiveClass).PairwiseDisjoint
      (fun provider ↦ providerSheets a x provider owner) := by
  exact pairwiseDisjoint_of_full_capacity
    (fun provider ↦ providerSheets a x provider owner)
    howner (provider_capacity_sum a x owner)

/-- Multiplication by two on `F_13^*/{+-1}`. -/
def doubleClass (i : ProjectiveClass) : ProjectiveClass :=
  ![1, 3, 5, 4, 2, 0] i

/-- Cyclic coordinate on the projective doubling orbit
`[1],[2],[4],[5],[3],[6]`. -/
def projectiveCycle (i : Fin 6) : ProjectiveClass :=
  ![0, 1, 3, 4, 2, 5] i

/-- Successor in cyclic coordinates. -/
def cycleSuccessor (i : Fin 6) : Fin 6 :=
  ⟨(i.val + 1) % 6, by omega⟩

/-- The explicit projective quotient really is the multiplication-by-two
six-cycle displayed in THM-970. -/
theorem doubleClass_projectiveCycle :
    ∀ i : Fin 6,
      doubleClass (projectiveCycle i) =
        projectiveCycle (cycleSuccessor i) := by
  decide

/-- The multiplication-by-two action is a permutation of all six classes. -/
theorem doubleClass_bijective : Function.Bijective doubleClass := by
  decide

/-- After six doublings every projective class returns to itself. -/
theorem doubleClass_six (i : ProjectiveClass) :
    doubleClass (doubleClass (doubleClass (doubleClass
      (doubleClass (doubleClass i))))) = i := by
  fin_cases i <;> rfl

/-- Five providers already carry the adjacent-owner contradiction.  The base
edge keeps projective classes `[1],[2],[3],[4],[6]` and omits `[5]`; subsequent
rows are its images under projective doubling. -/
abbrev CoreClass := Fin 5

def selectedClass (edge : Fin 6) (i : CoreClass) : ProjectiveClass :=
  ![![0, 1, 2, 3, 5],
    ![1, 3, 5, 4, 0],
    ![3, 4, 0, 2, 1],
    ![4, 2, 1, 5, 3],
    ![2, 5, 3, 0, 4],
    ![5, 0, 4, 1, 2]] edge i

theorem selectedClass_injective (edge : Fin 6) :
    Function.Injective (selectedClass edge) := by
  fin_cases edge <;> decide

/-- The first two selected classes are the endpoints of the corresponding
doubling edge. -/
theorem selected_owner_zero (edge : Fin 6) :
    selectedClass edge 0 = projectiveCycle edge := by
  fin_cases edge <;> rfl

theorem selected_owner_one (edge : Fin 6) :
    selectedClass edge 1 = doubleClass (projectiveCycle edge) := by
  fin_cases edge <;> rfl

/-- Sign and unit assignments on the five-provider core. -/
abbrev CoreSigns := CoreClass → Bool
abbrev CoreUnits := CoreClass → UnitDigit

def coreSignedLabel (edge : Fin 6) (b : CoreSigns) (i : CoreClass) : Nat :=
  if b i then 13 - representative (selectedClass edge i)
  else representative (selectedClass edge i)

def coreRatio (edge : Fin 6) (b : CoreSigns)
    (provider owner : CoreClass) : Fin 13 :=
  ⟨coreSignedLabel edge b provider *
      inverse13 (coreSignedLabel edge b owner) % 13,
    Nat.mod_lt _ (by omega)⟩

def coreProviderMask (edge : Fin 6) (b : CoreSigns) (u : CoreUnits)
    (provider owner : CoreClass) : Nat :=
  sheetMask (coreRatio edge b provider owner).val (u provider)

/-- Bit-mask disjointness. -/
def maskDisjoint (m n : Nat) : Bool :=
  (m &&& n) == 0

/-- For literal scale-ten masks, bit-mask disjointness agrees with disjointness
of the represented sheet finsets.  Only `13^2 * 4^2` local rows occur. -/
theorem literalMask_disjoint_iff :
    ∀ (q r : Fin 13) (d e : UnitDigit),
      Disjoint (maskSheets (sheetMask q.val d))
          (maskSheets (sheetMask r.val e)) ↔
        maskDisjoint (sheetMask q.val d) (sheetMask r.val e) = true := by
  decide

/-- Pairwise disjointness of five masks, written as the ten edges of `K_5`. -/
def fivePairwise (m0 m1 m2 m3 m4 : Nat) : Bool :=
  maskDisjoint m0 m1 && maskDisjoint m0 m2 &&
  maskDisjoint m0 m3 && maskDisjoint m0 m4 &&
  maskDisjoint m1 m2 && maskDisjoint m1 m3 &&
  maskDisjoint m1 m4 && maskDisjoint m2 m3 &&
  maskDisjoint m2 m4 && maskDisjoint m3 m4

def coreOwnerPartition (edge : Fin 6) (b : CoreSigns) (u : CoreUnits)
    (owner : CoreClass) : Bool :=
  fivePairwise
    (coreProviderMask edge b u 0 owner)
    (coreProviderMask edge b u 1 owner)
    (coreProviderMask edge b u 2 owner)
    (coreProviderMask edge b u 3 owner)
    (coreProviderMask edge b u 4 owner)

def adjacentCore (edge : Fin 6) (b : CoreSigns) (u : CoreUnits) : Bool :=
  coreOwnerPartition edge b u 0 && coreOwnerPartition edge b u 1

/-- Simultaneously reversing all five selected signs does not change a
provider/owner ratio.  Gauge the first selected sign to `false`. -/
def gaugeSigns (b : CoreSigns) : CoreSigns :=
  fun i ↦ b i != b 0

theorem coreRatio_gauge :
    ∀ (b : CoreSigns) (provider owner : CoreClass),
      coreRatio 0 b provider owner = coreRatio 0 (gaugeSigns b) provider owner := by
  decide

theorem coreProviderMask_gauge (b : CoreSigns) (u : CoreUnits)
    (provider owner : CoreClass) :
    coreProviderMask 0 b u provider owner =
      coreProviderMask 0 (gaugeSigns b) u provider owner := by
  unfold coreProviderMask
  rw [coreRatio_gauge b provider owner]

theorem coreOwnerPartition_gauge (b : CoreSigns) (u : CoreUnits)
    (owner : CoreClass) :
    coreOwnerPartition 0 b u owner =
      coreOwnerPartition 0 (gaugeSigns b) u owner := by
  unfold coreOwnerPartition
  rw [coreProviderMask_gauge b u 0 owner,
    coreProviderMask_gauge b u 1 owner,
    coreProviderMask_gauge b u 2 owner,
    coreProviderMask_gauge b u 3 owner,
    coreProviderMask_gauge b u 4 owner]

theorem adjacentCore_gauge (b : CoreSigns) (u : CoreUnits) :
    adjacentCore 0 b u = adjacentCore 0 (gaugeSigns b) u := by
  unfold adjacentCore
  rw [coreOwnerPartition_gauge b u 0, coreOwnerPartition_gauge b u 1]

/-- Combined disjointness of one provider pair at the two adjacent owners. -/
def literalPairAllowed (b : CoreSigns) (provider₁ provider₂ : CoreClass)
    (u₁ u₂ : UnitDigit) : Bool :=
  maskDisjoint
      (sheetMask (coreRatio 0 b provider₁ 0).val u₁)
      (sheetMask (coreRatio 0 b provider₂ 0).val u₂) &&
    maskDisjoint
      (sheetMask (coreRatio 0 b provider₁ 1).val u₁)
      (sheetMask (coreRatio 0 b provider₂ 1).val u₂)

/-- A 16-bit binary relation on two order-ten unit digits. -/
def frozenRelation (code : Nat) (u v : UnitDigit) : Bool :=
  code.testBit (4 * u.val + v.val)

/-- The eight-edge CSP core common to all sign blocks.  The constants are the
exact two-owner disjointness relations induced by the literal CRT masks. -/
def frozenCore (b₁ b₂ b₃ b₄ : Bool)
    (u₀ u₁ u₂ u₃ u₄ : UnitDigit) : Bool :=
  frozenRelation (if b₁ then 60855 else 31710) u₀ u₁ &&
  frozenRelation (if b₂ then 50595 else 14940) u₀ u₂ &&
  frozenRelation 38505 u₀ u₃ &&
  frozenRelation (if b₄ then 48765 else 55275) u₀ u₄ &&
  frozenRelation (if b₁ != b₂ then 50595 else 14940) u₁ u₂ &&
  frozenRelation 38505 u₁ u₄ &&
  frozenRelation (if b₂ != b₃ then 44085 else 21450) u₂ u₃ &&
  frozenRelation (if b₃ != b₄ then 50595 else 14940) u₃ u₄

/-- Reversing the unit-digit order, the permutation
`0 <-> 3, 1 <-> 2`. -/
def flipUnit (b : Bool) (u : UnitDigit) : UnitDigit :=
  if b then ⟨3 - u.val, by omega⟩ else u

/- The apparent sixteen sign-dependent CSPs are all the same object: changing
one selected sign reverses that vertex's unit alphabet.  These eight local
identities expose that switching symmetry without any global search. -/
theorem frozenRelation_01_normalized :
    ∀ (b : Bool) (u v : UnitDigit),
      frozenRelation (if b then 60855 else 31710) u v =
        frozenRelation 31710 u (flipUnit b v) := by
  intro b
  cases b <;> decide

theorem frozenRelation_02_normalized :
    ∀ (b : Bool) (u v : UnitDigit),
      frozenRelation (if b then 50595 else 14940) u v =
        frozenRelation 14940 u (flipUnit b v) := by
  intro b
  cases b <;> decide

theorem frozenRelation_03_normalized :
    ∀ (b : Bool) (u v : UnitDigit),
      frozenRelation 38505 u v =
        frozenRelation 38505 u (flipUnit b v) := by
  intro b
  cases b <;> decide

theorem frozenRelation_04_normalized :
    ∀ (b : Bool) (u v : UnitDigit),
      frozenRelation (if b then 48765 else 55275) u v =
        frozenRelation 55275 u (flipUnit b v) := by
  intro b
  cases b <;> decide

theorem frozenRelation_12_normalized :
    ∀ (b c : Bool) (u v : UnitDigit),
      frozenRelation (if b != c then 50595 else 14940) u v =
        frozenRelation 14940 (flipUnit b u) (flipUnit c v) := by
  intro b c
  cases b <;> cases c <;> decide

theorem frozenRelation_14_normalized :
    ∀ (b c : Bool) (u v : UnitDigit),
      frozenRelation 38505 u v =
        frozenRelation 38505 (flipUnit b u) (flipUnit c v) := by
  intro b c
  cases b <;> cases c <;> decide

theorem frozenRelation_23_normalized :
    ∀ (b c : Bool) (u v : UnitDigit),
      frozenRelation (if b != c then 44085 else 21450) u v =
        frozenRelation 21450 (flipUnit b u) (flipUnit c v) := by
  intro b c
  cases b <;> cases c <;> decide

theorem frozenRelation_34_normalized :
    ∀ (b c : Bool) (u v : UnitDigit),
      frozenRelation (if b != c then 50595 else 14940) u v =
        frozenRelation 14940 (flipUnit b u) (flipUnit c v) := by
  intro b c
  cases b <;> cases c <;> decide

/-- Vertex switching removes all four residual sign bits from the frozen CSP. -/
theorem frozenCore_sign_normalized
    (b₁ b₂ b₃ b₄ : Bool) (u₀ u₁ u₂ u₃ u₄ : UnitDigit) :
    frozenCore b₁ b₂ b₃ b₄ u₀ u₁ u₂ u₃ u₄ =
      frozenCore false false false false u₀
        (flipUnit b₁ u₁) (flipUnit b₂ u₂)
        (flipUnit b₃ u₃) (flipUnit b₄ u₄) := by
  unfold frozenCore
  rw [frozenRelation_01_normalized b₁ u₀ u₁,
    frozenRelation_02_normalized b₂ u₀ u₂,
    frozenRelation_03_normalized b₃ u₀ u₃,
    frozenRelation_04_normalized b₄ u₀ u₄,
    frozenRelation_12_normalized b₁ b₂ u₁ u₂,
    frozenRelation_14_normalized b₁ b₄ u₁ u₄,
    frozenRelation_23_normalized b₂ b₃ u₂ u₃,
    frozenRelation_34_normalized b₃ b₄ u₃ u₄]
  simp

/-- Each lookup below has only sixteen unit pairs after its sign cases are
split.  Together they certify every constant in `frozenCore` directly from the
literal CRT masks. -/
theorem literalPairAllowed_01 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 0 1 u v =
        frozenRelation (if b₁ then 60855 else 31710) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_02 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 0 2 u v =
        frozenRelation (if b₂ then 50595 else 14940) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_03 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 0 3 u v =
        frozenRelation 38505 u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_04 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 0 4 u v =
        frozenRelation (if b₄ then 48765 else 55275) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_12 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 1 2 u v =
        frozenRelation (if b₁ != b₂ then 50595 else 14940) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_14 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 1 4 u v =
        frozenRelation 38505 u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_23 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 2 3 u v =
        frozenRelation (if b₂ != b₃ then 44085 else 21450) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

theorem literalPairAllowed_34 :
    ∀ (b₁ b₂ b₃ b₄ : Bool) (u v : UnitDigit),
      literalPairAllowed ![false, b₁, b₂, b₃, b₄] 3 4 u v =
        frozenRelation (if b₃ != b₄ then 50595 else 14940) u v := by
  intro b₁ b₂ b₃ b₄
  cases b₁ <;> cases b₂ <;> cases b₃ <;> cases b₄ <;> decide

/- The command-line frontend otherwise checks independent theorem declarations
asynchronously.  The option above serializes this finite certificate bank so
each 256-row kernel proof is released before the next one starts. -/

/-- The frozen terminal CSP is split into 64 sign/first-unit blocks.
Each ordinary kernel `decide` proof checks exactly `4^4 = 256` unit words. -/
theorem adjacentCore_false_0_0000_u0_0 :
    ∀ (u1 u2 u3 u4 : UnitDigit),
      frozenCore false false false false
        0 u1 u2 u3 u4 = false := by
  decide

theorem adjacentCore_false_0_0000_u0_1 :
    ∀ (u1 u2 u3 u4 : UnitDigit),
      frozenCore false false false false
        1 u1 u2 u3 u4 = false := by
  decide

theorem adjacentCore_false_0_0000_u0_2 :
    ∀ (u1 u2 u3 u4 : UnitDigit),
      frozenCore false false false false
        2 u1 u2 u3 u4 = false := by
  decide

theorem adjacentCore_false_0_0000_u0_3 :
    ∀ (u1 u2 u3 u4 : UnitDigit),
      frozenCore false false false false
        3 u1 u2 u3 u4 = false := by
  decide

theorem adjacentCore_false_0_0000 :
    ∀ (u0 u1 u2 u3 u4 : UnitDigit),
      frozenCore false false false false
        u0 u1 u2 u3 u4 = false := by
  intro u0
  fin_cases u0
  · exact adjacentCore_false_0_0000_u0_0
  · exact adjacentCore_false_0_0000_u0_1
  · exact adjacentCore_false_0_0000_u0_2
  · exact adjacentCore_false_0_0000_u0_3

/-- Transport the single `4^5 = 1,024`-row certificate across the vertex
switching symmetry, thereby covering all `2^4 * 4^5 = 16,384` gauge-fixed
rows. -/
theorem frozenCore_false_gauged :
    ∀ (b1 b2 b3 b4 : Bool) (u0 u1 u2 u3 u4 : UnitDigit),
      frozenCore b1 b2 b3 b4 u0 u1 u2 u3 u4 = false := by
  intro b1 b2 b3 b4 u0 u1 u2 u3 u4
  rw [frozenCore_sign_normalized]
  exact adjacentCore_false_0_0000 u0
    (flipUnit b1 u1) (flipUnit b2 u2) (flipUnit b3 u3) (flipUnit b4 u4)

/-- The literal five-provider partition predicate implies all eight frozen
binary relations.  No exhaustive search occurs here: the proof only projects
eight pair constraints from the two `K_5` partition conditions and applies the
small lookup lemmas above. -/
theorem adjacentCore_implies_frozen
    (b1 b2 b3 b4 : Bool) (u0 u1 u2 u3 u4 : UnitDigit)
    (h : adjacentCore 0 ![false, b1, b2, b3, b4]
      ![u0, u1, u2, u3, u4] = true) :
    frozenCore b1 b2 b3 b4 u0 u1 u2 u3 u4 = true := by
  have hparts :
      coreOwnerPartition 0 ![false, b1, b2, b3, b4]
          ![u0, u1, u2, u3, u4] 0 = true ∧
        coreOwnerPartition 0 ![false, b1, b2, b3, b4]
          ![u0, u1, u2, u3, u4] 1 = true := by
    simpa [adjacentCore] using h
  have h0 : coreOwnerPartition 0 ![false, b1, b2, b3, b4]
      ![u0, u1, u2, u3, u4] 0 = true := hparts.1
  have h1 : coreOwnerPartition 0 ![false, b1, b2, b3, b4]
      ![u0, u1, u2, u3, u4] 1 = true := hparts.2
  simp only [coreOwnerPartition, fivePairwise, Bool.and_eq_true,
    coreProviderMask] at h0 h1
  simp only [frozenCore, Bool.and_eq_true]
  rw [← literalPairAllowed_01 b1 b2 b3 b4 u0 u1,
    ← literalPairAllowed_02 b1 b2 b3 b4 u0 u2,
    ← literalPairAllowed_03 b1 b2 b3 b4 u0 u3,
    ← literalPairAllowed_04 b1 b2 b3 b4 u0 u4,
    ← literalPairAllowed_12 b1 b2 b3 b4 u1 u2,
    ← literalPairAllowed_14 b1 b2 b3 b4 u1 u4,
    ← literalPairAllowed_23 b1 b2 b3 b4 u2 u3,
    ← literalPairAllowed_34 b1 b2 b3 b4 u3 u4]
  simp only [literalPairAllowed, Bool.and_eq_true]
  aesop

/-- The literal adjacent-owner core is empty after fixing the global sign
gauge. -/
theorem adjacentCore_false_0_gauged :
    ∀ (b1 b2 b3 b4 : Bool) (u0 u1 u2 u3 u4 : UnitDigit),
      adjacentCore 0 ![false, b1, b2, b3, b4]
        ![u0, u1, u2, u3, u4] = false := by
  intro b1 b2 b3 b4 u0 u1 u2 u3 u4
  cases hcore : adjacentCore 0 ![false, b1, b2, b3, b4]
      ![u0, u1, u2, u3, u4] with
  | false => rfl
  | true =>
      have hfrozen := adjacentCore_implies_frozen
        b1 b2 b3 b4 u0 u1 u2 u3 u4 hcore
      have hempty := frozenCore_false_gauged
        b1 b2 b3 b4 u0 u1 u2 u3 u4
      simp [hempty] at hfrozen

theorem coreUnits_eta (u : CoreUnits) :
    (![u 0, u 1, u 2, u 3, u 4] : CoreUnits) = u := by
  funext i
  fin_cases i <;> rfl

theorem gaugeSigns_eta (b : CoreSigns) :
    (![false, gaugeSigns b 1, gaugeSigns b 2, gaugeSigns b 3, gaugeSigns b 4] :
      CoreSigns) = gaugeSigns b := by
  funext i
  fin_cases i <;> simp [gaugeSigns]

theorem adjacentCore_false_0 (b : CoreSigns) (u : CoreUnits) :
    adjacentCore 0 b u = false := by
  rw [adjacentCore_gauge]
  rw [← gaugeSigns_eta b]
  rw [← coreUnits_eta u]
  exact adjacentCore_false_0_gauged
    (gaugeSigns b 1) (gaugeSigns b 2) (gaugeSigns b 3) (gaugeSigns b 4)
    (u 0) (u 1) (u 2) (u 3) (u 4)

/-- Sign correction incurred when multiplication by `2^edge` crosses the
chosen positive representative of a projective class. -/
def edgeSignCorrection (edge : Fin 6) (i : CoreClass) : Bool :=
  ![![false, false, false, false, false],
    ![false, false, false, true,  true],
    ![false, true,  true,  false, true],
    ![true,  false, true,  false, true],
    ![false, false, true,  true,  false],
    ![false, true,  false, true,  true]] edge i

def normalizedSigns (edge : Fin 6) (b : CoreSigns) : CoreSigns :=
  fun i ↦ b i != edgeSignCorrection edge i

/-- Multiplying every selected label by the same power of two cancels in all
provider/owner ratios.  Only the representative sign correction remains. -/
theorem coreRatio_normalized :
    ∀ (edge : Fin 6) (b : CoreSigns) (provider owner : CoreClass),
      coreRatio edge b provider owner =
        coreRatio 0 (normalizedSigns edge b) provider owner := by
  decide

theorem coreProviderMask_normalized (edge : Fin 6) (b : CoreSigns)
    (u : CoreUnits) (provider owner : CoreClass) :
    coreProviderMask edge b u provider owner =
      coreProviderMask 0 (normalizedSigns edge b) u provider owner := by
  unfold coreProviderMask
  rw [coreRatio_normalized edge b provider owner]

theorem coreOwnerPartition_normalized (edge : Fin 6) (b : CoreSigns)
    (u : CoreUnits) (owner : CoreClass) :
    coreOwnerPartition edge b u owner =
      coreOwnerPartition 0 (normalizedSigns edge b) u owner := by
  unfold coreOwnerPartition
  rw [coreProviderMask_normalized edge b u 0 owner,
    coreProviderMask_normalized edge b u 1 owner,
    coreProviderMask_normalized edge b u 2 owner,
    coreProviderMask_normalized edge b u 3 owner,
    coreProviderMask_normalized edge b u 4 owner]

theorem adjacentCore_normalized (edge : Fin 6) (b : CoreSigns) (u : CoreUnits) :
    adjacentCore edge b u = adjacentCore 0 (normalizedSigns edge b) u := by
  unfold adjacentCore
  rw [coreOwnerPartition_normalized edge b u 0,
    coreOwnerPartition_normalized edge b u 1]

theorem adjacentCore_false (edge : Fin 6) (b : CoreSigns) (u : CoreUnits) :
    adjacentCore edge b u = false := by
  rw [adjacentCore_normalized]
  exact adjacentCore_false_0 (normalizedSigns edge b) u

def projectedSigns (edge : Fin 6) (a : SignWord) : CoreSigns :=
  fun i ↦ signAt a (selectedClass edge i)

def projectedUnits (edge : Fin 6) (x : UnitWord) : CoreUnits :=
  fun i ↦ unitAt x (selectedClass edge i)

theorem coreProviderMask_projected (edge : Fin 6) (a : SignWord)
    (x : UnitWord) (provider owner : CoreClass) :
    coreProviderMask edge (projectedSigns edge a) (projectedUnits edge x)
        provider owner =
      providerMask a x (selectedClass edge provider) (selectedClass edge owner) := by
  rfl

theorem corePair_disjoint_of_ownerSat (edge : Fin 6) (a : SignWord)
    (x : UnitWord) (owner provider₁ provider₂ : CoreClass)
    (hneq : provider₁ ≠ provider₂)
    (howner : ownerSat a x (selectedClass edge owner)) :
    maskDisjoint
        (coreProviderMask edge (projectedSigns edge a) (projectedUnits edge x)
          provider₁ owner)
        (coreProviderMask edge (projectedSigns edge a) (projectedUnits edge x)
          provider₂ owner) = true := by
  have hp := ownerSat_pairwise a x (selectedClass edge owner) howner
  have hselected : selectedClass edge provider₁ ≠ selectedClass edge provider₂ :=
    (selectedClass_injective edge).ne hneq
  have hdisjoint :
      Disjoint
        (providerSheets a x (selectedClass edge provider₁) (selectedClass edge owner))
        (providerSheets a x (selectedClass edge provider₂) (selectedClass edge owner)) :=
    hp (by simp) (by simp) hselected
  rw [providerSheets, providerMask] at hdisjoint
  have hmask := (literalMask_disjoint_iff
    (ratio a (selectedClass edge provider₁) (selectedClass edge owner))
    (ratio a (selectedClass edge provider₂) (selectedClass edge owner))
    (unitAt x (selectedClass edge provider₁))
    (unitAt x (selectedClass edge provider₂))).1 hdisjoint
  simpa [coreProviderMask_projected] using hmask

theorem ownerSat_implies_corePartition (edge : Fin 6) (a : SignWord)
    (x : UnitWord) (owner : CoreClass)
    (howner : ownerSat a x (selectedClass edge owner)) :
    coreOwnerPartition edge (projectedSigns edge a) (projectedUnits edge x)
      owner = true := by
  simp only [coreOwnerPartition, fivePairwise, Bool.and_eq_true]
  repeat' apply And.intro
  all_goals
    apply corePair_disjoint_of_ownerSat edge a x owner <;>
      first | decide | assumption

/-- Cyclic coordinate inverse. -/
def cycleIndex (i : ProjectiveClass) : Fin 6 :=
  ![0, 1, 4, 2, 3, 5] i

theorem selected_cycleIndex_zero (i : ProjectiveClass) :
    selectedClass (cycleIndex i) 0 = i := by
  fin_cases i <;> rfl

theorem selected_cycleIndex_one (i : ProjectiveClass) :
    selectedClass (cycleIndex i) 1 = doubleClass i := by
  fin_cases i <;> rfl

/-- Exact structural audit: adjacent owner obligations on the projective
doubling cycle are disjoint, uniformly over all sign sections and unit words. -/
theorem adjacent_owner_disjoint
    (a : SignWord) (x : UnitWord) (i : ProjectiveClass)
    (hi : ownerSat a x i) : ¬ ownerSat a x (doubleClass i) := by
  intro hj
  let edge := cycleIndex i
  have h0 : coreOwnerPartition edge (projectedSigns edge a)
      (projectedUnits edge x) 0 = true := by
    apply ownerSat_implies_corePartition edge a x 0
    simpa [edge, selected_cycleIndex_zero] using hi
  have h1 : coreOwnerPartition edge (projectedSigns edge a)
      (projectedUnits edge x) 1 = true := by
    apply ownerSat_implies_corePartition edge a x 1
    simpa [edge, selected_cycleIndex_one] using hj
  have hcore : adjacentCore edge (projectedSigns edge a)
      (projectedUnits edge x) = true := by
    simp [adjacentCore, h0, h1]
  rw [adjacentCore_false] at hcore
  contradiction

/-- No global unit word lies in all six owner obligations. -/
theorem no_common_unit_word
    (a : SignWord) (x : UnitWord)
    (hall : ∀ i : ProjectiveClass, ownerSat a x i) : False := by
  exact adjacent_owner_disjoint a x 0 (hall 0) (hall (doubleClass 0))

/-- Set-theoretic form: each edge of the projective doubling cycle is absent
from the owner-obligation nerve. -/
theorem obligation_sets_disjoint (a : SignWord) (i : ProjectiveClass) :
    Disjoint {x : UnitWord | ownerSat a x i}
      {x : UnitWord | ownerSat a x (doubleClass i)} := by
  rw [Set.disjoint_left]
  intro x hi hj
  exact adjacent_owner_disjoint a x i hi hj

#print axioms doubleClass_projectiveCycle
#print axioms adjacent_owner_disjoint
#print axioms no_common_unit_word
#print axioms obligation_sets_disjoint

end ScaleTenProjectivePrism
end LonelyRunner
