import Mathlib

/-!
# Mixed-period beat masks and the Hunter-tree ledger (THM-1217)

At a sum/difference beat denominator `q`, the danger predicate for a speed
with reduced period `Q` lifts to a master clock `L` in `L / Q` identical
fibres.  The first part below kernel-checks this product-cardinality model.

The main formal result is the five-mask rooted-tree ledger.  When masks are
adjoined in a root-first ordering, intersection with any earlier parent mask
is credit against the naive sum of cardinalities.  This is precisely the
Hunter spanning-tree estimate used by THM-1217.  The common-zero estimate is
included separately, as is the strict cardinality consumer that turns either
ledger into a block escape.

The number-theoretic identification `L = q / gcd(q,b₁,...,b₆) = lcmᵢ Qᵢ`,
the strict-radius mask-count formula, and the interval-to-consecutive-block
bridge remain analytic provider lemmas.  The final section kernel-checks the
complete arithmetic and exact four-point master-clock certificate for the corrected
`a = 79`, `q = 280`, `p = 41` packet.
-/

namespace LRC14.MixedPeriodBeatMaskTree

open scoped Classical

/-- Exact strict `1/14` window count on a reduced clock of size `Q`. -/
def windowCount (Q : ℕ) : ℕ :=
  2 * ((Q + 13) / 14) - 1

/-- The abstract lift of a reduced mask through `R` identical fibres.  The
analytic residue lift to a master clock is canonically bijective to this
product when `L = Q * R`. -/
def fibreLift {Q R : ℕ} (base : Finset (Fin Q)) :
    Finset (Fin Q × Fin R) :=
  base ×ˢ Finset.univ

theorem fibreLift_card {Q R : ℕ} (base : Finset (Fin Q)) :
    (fibreLift (R := R) base).card = base.card * R := by
  simp [fibreLift]

theorem strictMask_fibreLift_card
    {Q R : ℕ} (base : Finset (Fin Q))
    (hbase : base.card = windowCount Q) :
    (fibreLift (R := R) base).card = R * windowCount Q := by
  rw [fibreLift_card, hbase]
  exact Nat.mul_comm _ _

/-- Five masks sharing one residue save four units from the naive sum. -/
theorem five_union_common_point_ledger
    {α : Type*} [DecidableEq α]
    (x : α) (s0 s1 s2 s3 s4 : Finset α)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (hx3 : x ∈ s3) (hx4 : x ∈ s4) :
    (s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4).card + 4 ≤
      s0.card + s1.card + s2.card + s3.card + s4.card := by
  have h01 := Finset.card_union_add_card_inter s0 s1
  have h012 := Finset.card_union_add_card_inter (s0 ∪ s1) s2
  have h0123 := Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2) s3
  have h01234 :=
    Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2 ∪ s3) s4
  have hi01 : 1 ≤ (s0 ∩ s1).card := by
    exact Finset.card_pos.mpr ⟨x, by simp [hx0, hx1]⟩
  have hi012 : 1 ≤ ((s0 ∪ s1) ∩ s2).card := by
    exact Finset.card_pos.mpr ⟨x, by simp [hx0, hx2]⟩
  have hi0123 : 1 ≤ ((s0 ∪ s1 ∪ s2) ∩ s3).card := by
    exact Finset.card_pos.mpr ⟨x, by simp [hx0, hx3]⟩
  have hi01234 : 1 ≤ ((s0 ∪ s1 ∪ s2 ∪ s3) ∩ s4).card := by
    exact Finset.card_pos.mpr ⟨x, by simp [hx0, hx4]⟩
  omega

/-- Rooted-tree Hunter ledger for five masks.  The `pᵢ` are parent masks in
a root-first ordering: each is only required to lie in the union already
exposed when `sᵢ` is adjoined.  Choosing an actual earlier mask as `pᵢ`
specializes this to every rooted spanning tree, after relabelling vertices in
root-first order. -/
theorem five_union_rooted_tree_ledger
    {α : Type*} [DecidableEq α]
    (s0 s1 s2 s3 s4 p1 p2 p3 p4 : Finset α)
    (hp1 : p1 ⊆ s0)
    (hp2 : p2 ⊆ s0 ∪ s1)
    (hp3 : p3 ⊆ s0 ∪ s1 ∪ s2)
    (hp4 : p4 ⊆ s0 ∪ s1 ∪ s2 ∪ s3) :
    (s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4).card
        + (p1 ∩ s1).card + (p2 ∩ s2).card
        + (p3 ∩ s3).card + (p4 ∩ s4).card ≤
      s0.card + s1.card + s2.card + s3.card + s4.card := by
  have hi1 : p1 ∩ s1 ⊆ s0 ∩ s1 := by
    intro x hx
    simp only [Finset.mem_inter] at hx ⊢
    exact ⟨hp1 hx.1, hx.2⟩
  have hi2 : p2 ∩ s2 ⊆ (s0 ∪ s1) ∩ s2 := by
    intro x hx
    simp only [Finset.mem_inter] at hx ⊢
    exact ⟨hp2 hx.1, hx.2⟩
  have hi3 : p3 ∩ s3 ⊆ (s0 ∪ s1 ∪ s2) ∩ s3 := by
    intro x hx
    simp only [Finset.mem_inter] at hx ⊢
    exact ⟨hp3 hx.1, hx.2⟩
  have hi4 : p4 ∩ s4 ⊆ (s0 ∪ s1 ∪ s2 ∪ s3) ∩ s4 := by
    intro x hx
    simp only [Finset.mem_inter] at hx ⊢
    exact ⟨hp4 hx.1, hx.2⟩
  have hc1 := Finset.card_le_card hi1
  have hc2 := Finset.card_le_card hi2
  have hc3 := Finset.card_le_card hi3
  have hc4 := Finset.card_le_card hi4
  have h01 := Finset.card_union_add_card_inter s0 s1
  have h012 := Finset.card_union_add_card_inter (s0 ∪ s1) s2
  have h0123 := Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2) s3
  have h01234 :=
    Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2 ∪ s3) s4
  omega

/-- Subtraction-free cardinality consumer.  If a block were covered by the
five masks, its cardinality plus the tree intersection credit could not
exceed the sum of mask cardinalities. -/
theorem tree_block_cover_impossible
    {α : Type*} [DecidableEq α]
    (block union : Finset α) (credit maskSum : ℕ)
    (hledger : union.card + credit ≤ maskSum)
    (hlarge : maskSum < block.card + credit)
    (hcover : block ⊆ union) : False := by
  have hcard := Finset.card_le_card hcover
  omega

/-- Common-zero specialization of the same block consumer. -/
theorem common_point_block_cover_impossible
    {α : Type*} [DecidableEq α]
    (block union : Finset α) (maskSum : ℕ)
    (hledger : union.card + 4 ≤ maskSum)
    (hlarge : maskSum < block.card + 4)
    (hcover : block ⊆ union) : False := by
  exact tree_block_cover_impossible block union 4 maskSum hledger hlarge hcover

/-- `everyBlockEscapes U step B` records the exact run supplier without
pretending that a full cyclic mask has a finite escape run. -/
def everyBlockEscapes
    {α : Type*} [DecidableEq α]
    (U : Finset α) (step : α → ℕ → α) (B : ℕ) : Prop :=
  ∀ x, ∃ j < B, step x j ∉ U

/-- An exact cyclic-run supplier forces properness.  This makes the guard
missing from several informal run-length formulations explicit. -/
theorem proper_of_everyBlockEscapes
    {α : Type*} [Fintype α] [DecidableEq α] [Nonempty α]
    (U : Finset α) (step : α → ℕ → α) (B : ℕ)
    (hescape : everyBlockEscapes U step B) :
    U ≠ Finset.univ := by
  intro hfull
  obtain ⟨x⟩ := ‹Nonempty α›
  obtain ⟨j, _hj, hout⟩ := hescape x
  rw [hfull] at hout
  simp at hout

/-- Least circular residue modulo `q`; used for the exact packet audit. -/
def circularResidue (q n : ℕ) : ℕ :=
  min (n % q) ((q - n % q) % q)

/-- Strict danger mask on the normalized master clock. -/
def dangerMask (L r : ℕ) : Finset (Fin L) :=
  Finset.univ.filter fun p => 14 * circularResidue L (r * p.val) < L

theorem corrected_master_masks :
    dangerMask 4 2 = {0, 2} ∧
    dangerMask 4 3 = {0} ∧
    dangerMask 4 5 = {0} ∧
    dangerMask 4 6 = {0, 2} ∧
    dangerMask 4 7 = {0} ∧
    dangerMask 4 11 = dangerMask 4 7 := by
  decide

theorem corrected_master_union :
    dangerMask 4 2 ∪ dangerMask 4 3 ∪ dangerMask 4 5 ∪
        dangerMask 4 6 ∪ dangerMask 4 7 = {0, 2} := by
  decide

theorem corrected_witness_is_hole :
    (1 : Fin 4) ∉
      dangerMask 4 2 ∪ dangerMask 4 3 ∪ dangerMask 4 5 ∪
        dangerMask 4 6 ∪ dangerMask 4 7 := by
  decide

/-- A realizable full mixed beat clock with no period-one mask.  This is a
guardrail for the method, not an LRC covering packet. -/
theorem full_clock_without_period_one :
    dangerMask 16 17 ∪ dangerMask 16 35 ∪ dangerMask 16 53 ∪
        dangerMask 16 71 ∪ dangerMask 16 88 = Finset.univ := by
  decide

theorem full_clock_period_data :
    Nat.gcd 17 16 = 1 ∧ Nat.gcd 35 16 = 1 ∧
    Nat.gcd 53 16 = 1 ∧ Nat.gcd 71 16 = 1 ∧
    Nat.gcd 88 16 = 8 ∧ Nat.gcd 104 16 = 8 ∧
    104 - 88 = 16 := by
  norm_num

theorem corrected_period_data :
    Nat.gcd 140 280 = 140 ∧ Nat.gcd 210 280 = 70 ∧
    Nat.gcd 350 280 = 70 ∧ Nat.gcd 420 280 = 140 ∧
    Nat.gcd 490 280 = 70 ∧ Nat.gcd 770 280 = 70 ∧
    Nat.lcm 2 (Nat.lcm 4 (Nat.lcm 4 (Nat.lcm 2 (Nat.lcm 4 4)))) = 4 := by
  norm_num

theorem corrected_beat_block :
    280 * (14 * 11 + 1) ≤ 14 * 79 * 41 ∧
    14 * 79 * 41 ≤ 280 * (14 * 11 + 13) := by
  norm_num

theorem corrected_packet_residues :
    [1, 2, 3, 4, 5, 6, 79, 140, 210, 350, 420, 490, 770].map
        (fun b => circularResidue 280 (b * 41)) =
      [41, 82, 123, 116, 75, 34, 121, 140, 70, 70, 140, 70, 70] := by
  norm_num [circularResidue]

theorem corrected_packet_safe :
    [41, 82, 123, 116, 75, 34, 121, 140, 70, 70, 140, 70, 70].all
      (fun residue => decide (280 ≤ 14 * residue)) = true := by
  decide

theorem corrected_frontier_arithmetic :
    (140 : ℚ) / 79 < 13 / 6 ∧
    (70 : ℚ) / 79 < 397 / 432 ∧
    (79 : ℚ) * (1 / 140 + 1 / 210 + 1 / 350 +
      1 / 420 + 1 / 490 + 1 / 770) = 21804 / 13475 ∧
    (21804 : ℚ) / 13475 > 1 := by
  norm_num

#print axioms fibreLift_card
#print axioms strictMask_fibreLift_card
#print axioms five_union_common_point_ledger
#print axioms five_union_rooted_tree_ledger
#print axioms tree_block_cover_impossible
#print axioms common_point_block_cover_impossible
#print axioms proper_of_everyBlockEscapes
#print axioms corrected_master_masks
#print axioms corrected_master_union
#print axioms corrected_witness_is_hole
#print axioms full_clock_without_period_one
#print axioms full_clock_period_data
#print axioms corrected_period_data
#print axioms corrected_beat_block
#print axioms corrected_packet_residues
#print axioms corrected_packet_safe
#print axioms corrected_frontier_arithmetic

end LRC14.MixedPeriodBeatMaskTree
