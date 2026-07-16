/-
  TournamentH7.CompositionDefect — klein-2026-07-16-S314 (cont.9)

  Formalizes the composition-defect identity and the reversal (frame) invariance lemmas
  of the relation-as-object reflection (HYP-6948 / cont.8 anchors):

   • `Tournament.op`, `op_op`           : the reversal Tᵒᵖ (perspective flip); an involution.
   • `outDegree_op_add`                 : s_op(v) + s(v) = n − 1  (scores reversal-COVARIANT).
   • `dev_op`, `xLevel_op`              : dev flips sign; x = Σ dev² is reversal-INVARIANT
                                           (the second moment is the first frame-free functional).
   • `card_cyc_op`                      : the composition-defect count is reversal-invariant.
   • `card_cyc_eq_three_mul_anchored`   : #cyc = 3 · c3 — transitivity R∘R ⊆ R fails by
                                           exactly 3·c3 ordered witnesses (HYP-6948 core).
-/

import Mathlib
import TournamentH7.Basic

open Finset

namespace Tournament

variable {n : ℕ}

/-- The reversal (converse) tournament: every arc flipped. -/
def op (T : Tournament n) : Tournament n where
  arc i j := T.arc j i
  irrefl i := T.irrefl i
  total i j h := (T.total i j h).symm
  asym i j := T.asym j i

@[simp] theorem op_arc (T : Tournament n) (i j : Fin n) : (T.op).arc i j = T.arc j i := rfl

theorem op_op (T : Tournament n) : T.op.op = T := rfl

/-- Scores are reversal-covariant: s_op(v) + s(v) = n − 1. -/
theorem outDegree_op_add (T : Tournament n) (v : Fin n) :
    (T.op).outDegree v + T.outDegree v = n - 1 := by
  classical
  have hunion :
      (univ.filter (fun w => T.arc w v = true)) ∪ (univ.filter (fun w => T.arc v w = true))
        = univ.erase v := by
    ext w
    simp only [mem_union, mem_filter, mem_univ, true_and, mem_erase, and_true]
    constructor
    · rintro (h | h)
      · intro he; subst he; simp [T.irrefl w] at h
      · intro he; subst he; simp [T.irrefl w] at h
    · intro hw
      rcases T.total w v hw with h | h
      · exact Or.inl h
      · exact Or.inr h
  have hdisj :
      Disjoint (univ.filter (fun w => T.arc w v = true))
               (univ.filter (fun w => T.arc v w = true)) := by
    rw [Finset.disjoint_left]
    intro w hw hw'
    simp only [mem_filter, mem_univ, true_and] at hw hw'
    exact T.asym v w ⟨hw', hw⟩
  have hcard := Finset.card_union_of_disjoint hdisj
  rw [hunion] at hcard
  have herase : (univ.erase v).card = n - 1 := by
    rw [Finset.card_erase_of_mem (mem_univ v), Finset.card_univ, Fintype.card_fin]
  have h1 : (T.op).outDegree v = (univ.filter (fun w => T.arc w v = true)).card := rfl
  have h2 : T.outDegree v = (univ.filter (fun w => T.arc v w = true)).card := rfl
  omega

/-- Centered doubled score, over ℤ. -/
def dev (T : Tournament n) (v : Fin n) : ℤ := 2 * (T.outDegree v : ℤ) - ((n : ℤ) - 1)

theorem dev_op (T : Tournament n) (v : Fin n) : (T.op).dev v = - T.dev v := by
  have h := outDegree_op_add T v
  have hn : 1 ≤ n := v.pos
  have h' : ((T.op).outDegree v : ℤ) + (T.outDegree v : ℤ) = (n : ℤ) - 1 := by
    have hc := congrArg (fun m : ℕ => (m : ℤ)) h
    push_cast at hc
    omega
  unfold dev
  omega

/-- The axis coordinate x = Σ dev². -/
def xLevel (T : Tournament n) : ℤ := ∑ v : Fin n, (T.dev v) ^ 2

/-- The second moment is reversal-invariant: the first frame-free functional. -/
theorem xLevel_op (T : Tournament n) : (T.op).xLevel = T.xLevel := by
  unfold xLevel
  refine Finset.sum_congr rfl (fun v _ => ?_)
  rw [dev_op]
  ring

/-- Ordered cyclic triples: witnesses of R∘R ∩ Rᵒᵖ — the composition defect. -/
def cyc (T : Tournament n) : Finset (Fin n × Fin n × Fin n) :=
  univ.filter (fun p => T.arc p.1 p.2.1 = true ∧ T.arc p.2.1 p.2.2 = true ∧
                        T.arc p.2.2 p.1 = true)

theorem mem_cyc_iff (T : Tournament n) (p : Fin n × Fin n × Fin n) :
    p ∈ T.cyc ↔ T.arc p.1 p.2.1 = true ∧ T.arc p.2.1 p.2.2 = true ∧
                 T.arc p.2.2 p.1 = true := by
  unfold cyc; simp

/-- Entries of a cyclic triple are pairwise distinct. -/
theorem mem_cyc_distinct (T : Tournament n) {p : Fin n × Fin n × Fin n} (hp : p ∈ T.cyc) :
    p.1 ≠ p.2.1 ∧ p.2.1 ≠ p.2.2 ∧ p.2.2 ≠ p.1 := by
  rw [mem_cyc_iff] at hp
  obtain ⟨h1, h2, h3⟩ := hp
  refine ⟨?_, ?_, ?_⟩
  · intro he; rw [he] at h1; simp [T.irrefl] at h1
  · intro he; rw [he] at h2; simp [T.irrefl] at h2
  · intro he; rw [he] at h3; simp [T.irrefl] at h3

/-- Reversal invariance of the defect count: the bijection (a,b,c) ↦ (c,b,a). -/
theorem card_cyc_op (T : Tournament n) : (T.op).cyc.card = T.cyc.card := by
  classical
  refine Finset.card_bij (fun p _ => (p.2.2, p.2.1, p.1)) ?_ ?_ ?_
  · intro p hp
    rw [mem_cyc_iff] at hp ⊢
    simp only [op_arc] at hp
    exact ⟨hp.2.1, hp.1, hp.2.2⟩
  · intro p1 h1 p2 h2 h
    have e1 := congrArg (fun q : Fin n × Fin n × Fin n => q.2.2) h
    have e2 := congrArg (fun q : Fin n × Fin n × Fin n => q.2.1) h
    have e3 := congrArg (fun q : Fin n × Fin n × Fin n => q.1) h
    simp only at e1 e2 e3
    exact Prod.ext e1 (Prod.ext e2 e3)
  · intro q hq
    refine ⟨(q.2.2, q.2.1, q.1), ?_, rfl⟩
    rw [mem_cyc_iff] at hq
    rw [mem_cyc_iff]
    simp only [op_arc]
    exact ⟨hq.2.1, hq.1, hq.2.2⟩

/-- Rotation of a triple. -/
def rot (p : Fin n × Fin n × Fin n) : Fin n × Fin n × Fin n := (p.2.1, p.2.2, p.1)

@[simp] theorem rot_rot_rot (p : Fin n × Fin n × Fin n) : rot (rot (rot p)) = p := rfl

theorem rot_mem_cyc (T : Tournament n) {p : Fin n × Fin n × Fin n} (hp : p ∈ T.cyc) :
    rot p ∈ T.cyc := by
  rw [mem_cyc_iff] at hp ⊢
  exact ⟨hp.2.1, hp.2.2, hp.1⟩

/-- Min-anchored cyclic triples: canonical rotation representatives (c3). -/
def cycAnchored (T : Tournament n) : Finset (Fin n × Fin n × Fin n) :=
  T.cyc.filter (fun p => p.1 < p.2.1 ∧ p.1 < p.2.2)

theorem mem_cycAnchored_iff (T : Tournament n) (p : Fin n × Fin n × Fin n) :
    p ∈ T.cycAnchored ↔ p ∈ T.cyc ∧ p.1 < p.2.1 ∧ p.1 < p.2.2 := by
  unfold cycAnchored; simp

/-- Two different rotations of one triple cannot both be anchored. -/
theorem anchored_unique₁ (T : Tournament n) (p : Fin n × Fin n × Fin n)
    (h1 : p ∈ T.cycAnchored) (h2 : rot p ∈ T.cycAnchored) : False := by
  rw [mem_cycAnchored_iff] at h1 h2
  have a := h1.2.1
  have b := h2.2.2
  exact absurd a (not_lt.mpr (le_of_lt b))

theorem anchored_unique₂ (T : Tournament n) (p : Fin n × Fin n × Fin n)
    (h1 : p ∈ T.cycAnchored) (h2 : rot (rot p) ∈ T.cycAnchored) : False := by
  rw [mem_cycAnchored_iff] at h1 h2
  have a := h1.2.2
  have b := h2.2.1
  exact absurd a (not_lt.mpr (le_of_lt b))

theorem anchored_unique₃ (T : Tournament n) (p : Fin n × Fin n × Fin n)
    (h1 : rot p ∈ T.cycAnchored) (h2 : rot (rot p) ∈ T.cycAnchored) : False := by
  rw [mem_cycAnchored_iff] at h1 h2
  have a := h1.2.1
  have b := h2.2.2
  exact absurd a (not_lt.mpr (le_of_lt b))

/-- Some rotation of a cyclic triple is anchored (the least entry can lead). -/
theorem anchored_exists (T : Tournament n) {p : Fin n × Fin n × Fin n} (hp : p ∈ T.cyc) :
    p ∈ T.cycAnchored ∨ rot p ∈ T.cycAnchored ∨ rot (rot p) ∈ T.cycAnchored := by
  obtain ⟨d1, d2, d3⟩ := mem_cyc_distinct T hp
  have h1 : rot p ∈ T.cyc := rot_mem_cyc T hp
  have h2 : rot (rot p) ∈ T.cyc := rot_mem_cyc T h1
  have e1 : p.1.val ≠ p.2.1.val := fun h => d1 (Fin.ext h)
  have e2 : p.2.1.val ≠ p.2.2.val := fun h => d2 (Fin.ext h)
  have e3 : p.2.2.val ≠ p.1.val := fun h => d3 (Fin.ext h)
  have harith : (p.1.val < p.2.1.val ∧ p.1.val < p.2.2.val) ∨
      (p.2.1.val < p.2.2.val ∧ p.2.1.val < p.1.val) ∨
      (p.2.2.val < p.1.val ∧ p.2.2.val < p.2.1.val) := by omega
  rcases harith with ⟨a, b⟩ | ⟨a, b⟩ | ⟨a, b⟩
  · exact Or.inl ((mem_cycAnchored_iff T p).mpr
      ⟨hp, Fin.lt_def.mpr a, Fin.lt_def.mpr b⟩)
  · exact Or.inr (Or.inl ((mem_cycAnchored_iff T (rot p)).mpr
      ⟨h1, Fin.lt_def.mpr a, Fin.lt_def.mpr b⟩))
  · exact Or.inr (Or.inr ((mem_cycAnchored_iff T (rot (rot p))).mpr
      ⟨h2, Fin.lt_def.mpr a, Fin.lt_def.mpr b⟩))

/-- THE COMPOSITION-DEFECT IDENTITY: #cyc = 3 · c3. -/
theorem card_cyc_eq_three_mul_anchored (T : Tournament n) :
    T.cyc.card = 3 * (T.cycAnchored).card := by
  classical
  have hcover : T.cyc = (T.cycAnchored).biUnion (fun q => {q, rot q, rot (rot q)}) := by
    ext p
    simp only [Finset.mem_biUnion, Finset.mem_insert, Finset.mem_singleton]
    constructor
    · intro hp
      rcases anchored_exists T hp with h | h | h
      · exact ⟨p, h, Or.inl rfl⟩
      · exact ⟨rot p, h, Or.inr (Or.inr (rot_rot_rot p).symm)⟩
      · exact ⟨rot (rot p), h, Or.inr (Or.inl (rot_rot_rot p).symm)⟩
    · rintro ⟨q, hq, h | h | h⟩
      · rw [h]; exact ((mem_cycAnchored_iff T q).mp hq).1
      · subst h; exact rot_mem_cyc T ((mem_cycAnchored_iff T q).mp hq).1
      · subst h; exact rot_mem_cyc T (rot_mem_cyc T ((mem_cycAnchored_iff T q).mp hq).1)
  rw [hcover, Finset.card_biUnion]
  · have hblock : ∀ q ∈ T.cycAnchored, ({q, rot q, rot (rot q)} : Finset _).card = 3 := by
      intro q hq
      have hqc : q ∈ T.cyc := ((mem_cycAnchored_iff T q).mp hq).1
      obtain ⟨d1, d2, d3⟩ := mem_cyc_distinct T hqc
      have ne1 : q ≠ rot q := fun h => d1 (congrArg Prod.fst h)
      have ne2 : q ≠ rot (rot q) := fun h => d3 (congrArg Prod.fst h).symm
      have ne3 : rot q ≠ rot (rot q) := fun h => d2 (congrArg Prod.fst h)
      rw [Finset.card_insert_of_notMem (by simp [ne1, ne2]),
          Finset.card_insert_of_notMem (by simp [ne3]), Finset.card_singleton]
    rw [Finset.sum_congr rfl hblock, Finset.sum_const, smul_eq_mul, Nat.mul_comm]
  · intro q1 hq1 q2 hq2 hne
    simp only [Finset.disjoint_left, Finset.mem_insert, Finset.mem_singleton]
    rintro p (rfl | rfl | rfl) hp2
    · rcases hp2 with h | h | h
      · exact hne h
      · exact anchored_unique₁ T q2 hq2 (h ▸ hq1)
      · exact anchored_unique₂ T q2 hq2 (h ▸ hq1)
    · rcases hp2 with h | h | h
      · exact anchored_unique₁ T q1 hq1 (h.symm ▸ hq2)
      · refine hne ?_
        have hc := congrArg (fun x => rot (rot x)) h
        simpa using hc
      · have h'' : q1 = rot q2 := by
          have hc := congrArg (fun x => rot (rot x)) h
          simpa using hc
        exact anchored_unique₁ T q2 hq2 (h'' ▸ hq1)
    · rcases hp2 with h | h | h
      · exact anchored_unique₂ T q1 hq1 (h.symm ▸ hq2)
      · have h'' : q2 = rot q1 := by
          have hc := congrArg (fun x => rot (rot x)) h.symm
          simpa using hc
        exact anchored_unique₁ T q1 hq1 (h'' ▸ hq2)
      · refine hne ?_
        have hc := congrArg (fun x => rot x) h
        simpa using hc

end Tournament
