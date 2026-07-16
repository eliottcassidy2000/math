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



/-! ### The Kendall–Babington Smith score formula (klein-S315):
     c3 + Σ_v C(s_v, 2) = C(n, 3). -/

theorem arc_of_not_arc (T : Tournament n) {u v : Fin n} (h : u ≠ v)
    (hf : ¬ T.arc u v = true) : T.arc v u = true := by
  rcases T.total u v h with h1 | h1
  · exact absurd h1 hf
  · exact h1

/-- Ordered transitive triples (u beats v and w; v beats w): one per transitive 3-set. -/
def orderedTrans (T : Tournament n) : Finset (Fin n × Fin n × Fin n) :=
  univ.filter (fun p => T.arc p.1 p.2.1 = true ∧ T.arc p.2.1 p.2.2 = true ∧
                        T.arc p.1 p.2.2 = true)

theorem mem_orderedTrans_iff (T : Tournament n) (p : Fin n × Fin n × Fin n) :
    p ∈ T.orderedTrans ↔ T.arc p.1 p.2.1 = true ∧ T.arc p.2.1 p.2.2 = true ∧
                          T.arc p.1 p.2.2 = true := by
  unfold orderedTrans; simp

theorem mem_orderedTrans_distinct (T : Tournament n) {p : Fin n × Fin n × Fin n}
    (hp : p ∈ T.orderedTrans) : p.1 ≠ p.2.1 ∧ p.2.1 ≠ p.2.2 ∧ p.1 ≠ p.2.2 := by
  rw [mem_orderedTrans_iff] at hp
  obtain ⟨h1, h2, h3⟩ := hp
  refine ⟨?_, ?_, ?_⟩
  · intro he; rw [he] at h1; simp [T.irrefl] at h1
  · intro he; rw [he] at h2; simp [T.irrefl] at h2
  · intro he; rw [he] at h3; simp [T.irrefl] at h3

/-- The support 3-set of a triple. -/
def supp (p : Fin n × Fin n × Fin n) : Finset (Fin n) := {p.1, p.2.1, p.2.2}

theorem mem_supp_iff (p : Fin n × Fin n × Fin n) (y : Fin n) :
    y ∈ supp p ↔ y = p.1 ∨ y = p.2.1 ∨ y = p.2.2 := by
  unfold supp; simp

theorem supp_card_of_distinct {p : Fin n × Fin n × Fin n}
    (d1 : p.1 ≠ p.2.1) (d2 : p.2.1 ≠ p.2.2) (d3 : p.1 ≠ p.2.2) :
    (supp p).card = 3 := by
  unfold supp
  rw [Finset.card_insert_of_notMem (by simp [d1, d3]),
      Finset.card_insert_of_notMem (by simp [d2]), Finset.card_singleton]

theorem supp_rot (p : Fin n × Fin n × Fin n) : supp (rot p) = supp p := by
  unfold supp rot
  ext y
  simp only [Finset.mem_insert, Finset.mem_singleton]
  tauto

/-- No 3-set carries both a cyclic and a transitive orientation. -/
theorem cyc_trans_supp_ne (T : Tournament n) {p q : Fin n × Fin n × Fin n}
    (hp : p ∈ T.cyc) (hq : q ∈ T.orderedTrans) : supp p ≠ supp q := by
  intro hsupp
  rw [mem_cyc_iff] at hp
  rw [mem_orderedTrans_iff] at hq
  obtain ⟨c1, c2, c3⟩ := hp
  obtain ⟨t1, t2, t3⟩ := hq
  have hq1 : q.1 ∈ supp p := by rw [hsupp, mem_supp_iff]; left; rfl
  rw [mem_supp_iff] at hq1
  have key : ∃ z, z ∈ supp p ∧ T.arc z q.1 = true := by
    rcases hq1 with h | h | h
    · exact ⟨p.2.2, (mem_supp_iff p _).mpr (Or.inr (Or.inr rfl)), h ▸ c3⟩
    · exact ⟨p.1, (mem_supp_iff p _).mpr (Or.inl rfl), h ▸ c1⟩
    · exact ⟨p.2.1, (mem_supp_iff p _).mpr (Or.inr (Or.inl rfl)), h ▸ c2⟩
  obtain ⟨z, hz, hzq⟩ := key
  have hz' : z ∈ supp q := by rw [← hsupp]; exact hz
  rw [mem_supp_iff] at hz'
  rcases hz' with h | h | h
  · subst h; simp [T.irrefl] at hzq
  · subst h; exact T.asym q.1 q.2.1 ⟨t1, hzq⟩
  · subst h; exact T.asym q.1 q.2.2 ⟨t3, hzq⟩

/-- Same-support cyclic triples are rotations of one another. -/
theorem cyc_supp_rot (T : Tournament n) {p q : Fin n × Fin n × Fin n}
    (hp : p ∈ T.cyc) (hq : q ∈ T.cyc) (hsupp : supp p = supp q) :
    q = p ∨ q = rot p ∨ q = rot (rot p) := by
  have c1 := ((mem_cyc_iff T p).mp hp).1
  have c2 := ((mem_cyc_iff T p).mp hp).2.1
  have c3 := ((mem_cyc_iff T p).mp hp).2.2
  have e1 := ((mem_cyc_iff T q).mp hq).1
  have e2 := ((mem_cyc_iff T q).mp hq).2.1
  have e3 := ((mem_cyc_iff T q).mp hq).2.2
  obtain ⟨dq1, dq2, dq3⟩ := mem_cyc_distinct T hq
  have hmem : ∀ y, y ∈ supp q → y = p.1 ∨ y = p.2.1 ∨ y = p.2.2 := by
    intro y hy; rw [← hsupp, mem_supp_iff] at *; exact hy
  have h1 := hmem q.1 ((mem_supp_iff q _).mpr (Or.inl rfl))
  have h2 := hmem q.2.1 ((mem_supp_iff q _).mpr (Or.inr (Or.inl rfl)))
  have h3 := hmem q.2.2 ((mem_supp_iff q _).mpr (Or.inr (Or.inr rfl)))
  rcases h1 with k1 | k1 | k1
  · left
    have k2 : q.2.1 = p.2.1 := by
      rcases h2 with m | m | m
      · exact absurd (k1.trans m.symm) dq1
      · exact m
      · exfalso; rw [k1, m] at e1; exact T.asym p.2.2 p.1 ⟨c3, e1⟩
    have k3 : q.2.2 = p.2.2 := by
      rcases h3 with m | m | m
      · exact absurd (m.trans k1.symm) dq3
      · exact absurd (k2.trans m.symm) dq2
      · exact m
    exact Prod.ext k1 (Prod.ext k2 k3)
  · right; left
    have k2 : q.2.1 = p.2.2 := by
      rcases h2 with m | m | m
      · exfalso; rw [k1, m] at e1; exact T.asym p.1 p.2.1 ⟨c1, e1⟩
      · exact absurd (k1.trans m.symm) dq1
      · exact m
    have k3 : q.2.2 = p.1 := by
      rcases h3 with m | m | m
      · exact m
      · exact absurd (m.trans k1.symm) dq3
      · exact absurd (k2.trans m.symm) dq2
    exact Prod.ext k1 (Prod.ext k2 k3)
  · right; right
    have k2 : q.2.1 = p.1 := by
      rcases h2 with m | m | m
      · exact m
      · exfalso; rw [k1, m] at e1; exact T.asym p.2.1 p.2.2 ⟨c2, e1⟩
      · exact absurd (k1.trans m.symm) dq1
    have k3 : q.2.2 = p.2.1 := by
      rcases h3 with m | m | m
      · exact absurd (k2.trans m.symm) dq2
      · exact m
      · exact absurd (m.trans k1.symm) dq3
    exact Prod.ext k1 (Prod.ext k2 k3)

/-- Same-support ordered transitive triples are equal: source and order are forced. -/
theorem trans_supp_eq (T : Tournament n) {p q : Fin n × Fin n × Fin n}
    (hp : p ∈ T.orderedTrans) (hq : q ∈ T.orderedTrans) (hsupp : supp p = supp q) :
    q = p := by
  have t1 := ((mem_orderedTrans_iff T p).mp hp).1
  have t2 := ((mem_orderedTrans_iff T p).mp hp).2.1
  have t3 := ((mem_orderedTrans_iff T p).mp hp).2.2
  have u1 := ((mem_orderedTrans_iff T q).mp hq).1
  have u2 := ((mem_orderedTrans_iff T q).mp hq).2.1
  have u3 := ((mem_orderedTrans_iff T q).mp hq).2.2
  obtain ⟨dp1, dp2, dp3⟩ := mem_orderedTrans_distinct T hp
  obtain ⟨dq1, dq2, dq3⟩ := mem_orderedTrans_distinct T hq
  have hmem : ∀ y, y ∈ supp q → y = p.1 ∨ y = p.2.1 ∨ y = p.2.2 := by
    intro y hy
    rw [← hsupp] at hy
    exact (mem_supp_iff p y).mp hy
  have h1 := hmem q.1 ((mem_supp_iff q q.1).mpr (Or.inl rfl))
  have h2 := hmem q.2.1 ((mem_supp_iff q q.2.1).mpr (Or.inr (Or.inl rfl)))
  have h3 := hmem q.2.2 ((mem_supp_iff q q.2.2).mpr (Or.inr (Or.inr rfl)))
  have hp1 : p.1 = q.1 ∨ p.1 = q.2.1 ∨ p.1 = q.2.2 := by
    have hin : p.1 ∈ supp q := by
      rw [← hsupp]; exact (mem_supp_iff p p.1).mpr (Or.inl rfl)
    exact (mem_supp_iff q p.1).mp hin
  have k1 : q.1 = p.1 := by
    rcases h1 with m | m | m
    · exact m
    · exfalso
      rcases hp1 with w | w | w
      · exact dp1 (w.trans m)
      · have hu := u1; rw [m, ← w] at hu; exact T.asym p.1 p.2.1 ⟨t1, hu⟩
      · have hu := u3; rw [m, ← w] at hu; exact T.asym p.1 p.2.1 ⟨t1, hu⟩
    · exfalso
      rcases hp1 with w | w | w
      · exact dp3 (w.trans m)
      · have hu := u1; rw [m, ← w] at hu; exact T.asym p.1 p.2.2 ⟨t3, hu⟩
      · have hu := u3; rw [m, ← w] at hu; exact T.asym p.1 p.2.2 ⟨t3, hu⟩
  have k2 : q.2.1 = p.2.1 := by
    rcases h2 with m | m | m
    · exact absurd (k1.trans m.symm) dq1
    · exact m
    · exfalso
      have k3' : q.2.2 = p.2.1 := by
        rcases h3 with m3 | m3 | m3
        · exact absurd (k1.trans m3.symm) dq3
        · exact m3
        · exact absurd (m.trans m3.symm) dq2
      rw [m, k3'] at u2
      exact T.asym p.2.1 p.2.2 ⟨t2, u2⟩
  have k3 : q.2.2 = p.2.2 := by
    rcases h3 with m | m | m
    · exact absurd (k1.trans m.symm) dq3
    · exact absurd (k2.trans m.symm) dq2
    · exact m
  exact Prod.ext k1 (Prod.ext k2 k3)

/-- Every 3-set of vertices is the support of an anchored-cyclic or transitive triple. -/
theorem exists_rep (T : Tournament n) {a b c : Fin n}
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c) :
    ∃ p ∈ T.cycAnchored ∪ T.orderedTrans, supp p = ({a, b, c} : Finset (Fin n)) := by
  classical
  have setperm : ∀ (x y z : Fin n), ({x, y, z} : Finset (Fin n)) = {a, b, c} →
      True := fun _ _ _ _ => trivial
  by_cases h1 : T.arc a b = true <;> by_cases h2 : T.arc b c = true <;>
    by_cases h3 : T.arc a c = true
  · exact ⟨(a, b, c), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h1, h2, h3⟩)), rfl⟩
  · have h3' : T.arc c a = true := arc_of_not_arc T hac h3
    have hcy : (a, b, c) ∈ T.cyc := (mem_cyc_iff T _).mpr ⟨h1, h2, h3'⟩
    rcases anchored_exists T hcy with h | h | h
    · exact ⟨(a, b, c), Finset.mem_union.mpr (Or.inl h), rfl⟩
    · exact ⟨rot (a, b, c), Finset.mem_union.mpr (Or.inl h), by rw [supp_rot]; rfl⟩
    · exact ⟨rot (rot (a, b, c)), Finset.mem_union.mpr (Or.inl h), by
        rw [supp_rot, supp_rot]; rfl⟩
  · have h2' : T.arc c b = true := arc_of_not_arc T hbc h2
    refine ⟨(a, c, b), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h3, h2', h1⟩)), ?_⟩
    unfold supp; ext y; simp; tauto
  · have h2' : T.arc c b = true := arc_of_not_arc T hbc h2
    have h3' : T.arc c a = true := arc_of_not_arc T hac h3
    refine ⟨(c, a, b), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h3', h1, h2'⟩)), ?_⟩
    unfold supp; ext y; simp; tauto
  · have h1' : T.arc b a = true := arc_of_not_arc T hab h1
    refine ⟨(b, a, c), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h1', h3, h2⟩)), ?_⟩
    unfold supp; ext y; simp; tauto
  · have h1' : T.arc b a = true := arc_of_not_arc T hab h1
    have h3' : T.arc c a = true := arc_of_not_arc T hac h3
    refine ⟨(b, c, a), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h2, h3', h1'⟩)), ?_⟩
    unfold supp; ext y; simp; tauto
  · have h1' : T.arc b a = true := arc_of_not_arc T hab h1
    have h2' : T.arc c b = true := arc_of_not_arc T hbc h2
    have hcy : (a, c, b) ∈ T.cyc := (mem_cyc_iff T _).mpr ⟨h3, h2', h1'⟩
    rcases anchored_exists T hcy with h | h | h
    · refine ⟨(a, c, b), Finset.mem_union.mpr (Or.inl h), ?_⟩
      unfold supp; ext y; simp; tauto
    · refine ⟨rot (a, c, b), Finset.mem_union.mpr (Or.inl h), ?_⟩
      rw [supp_rot]; unfold supp; ext y; simp; tauto
    · refine ⟨rot (rot (a, c, b)), Finset.mem_union.mpr (Or.inl h), ?_⟩
      rw [supp_rot, supp_rot]; unfold supp; ext y; simp; tauto
  · have h1' : T.arc b a = true := arc_of_not_arc T hab h1
    have h2' : T.arc c b = true := arc_of_not_arc T hbc h2
    have h3' : T.arc c a = true := arc_of_not_arc T hac h3
    refine ⟨(c, b, a), Finset.mem_union.mpr (Or.inr ((mem_orderedTrans_iff T _).mpr
      ⟨h2', h1', h3'⟩)), ?_⟩
    unfold supp; ext y; simp; tauto

/-- THE 3-SET PARTITION: c3 + #transitive-3-sets = C(n, 3). -/
theorem cycAnchored_card_add_orderedTrans_card (T : Tournament n) :
    (T.cycAnchored).card + (T.orderedTrans).card = n.choose 3 := by
  classical
  have hdisj : Disjoint T.cycAnchored T.orderedTrans := by
    rw [Finset.disjoint_left]
    intro p hp hq
    exact cyc_trans_supp_ne T ((mem_cycAnchored_iff T p).mp hp).1 hq rfl
  rw [← Finset.card_union_of_disjoint hdisj]
  have hpow : ((univ : Finset (Fin n)).powersetCard 3).card = n.choose 3 := by
    rw [Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]
  rw [← hpow]
  refine Finset.card_bij (fun p _ => supp p) ?_ ?_ ?_
  · intro p hp
    rw [Finset.mem_powersetCard]
    rcases Finset.mem_union.mp hp with h | h
    · obtain ⟨d1, d2, d3⟩ := mem_cyc_distinct T ((mem_cycAnchored_iff T p).mp h).1
      exact ⟨Finset.subset_univ _, supp_card_of_distinct d1 d2 d3.symm⟩
    · obtain ⟨d1, d2, d3⟩ := mem_orderedTrans_distinct T h
      exact ⟨Finset.subset_univ _, supp_card_of_distinct d1 d2 d3⟩
  · intro p1 h1 p2 h2 heq
    rcases Finset.mem_union.mp h1 with a1 | a1 <;> rcases Finset.mem_union.mp h2 with a2 | a2
    · have heq' : supp p1 = supp p2 := heq
      rcases cyc_supp_rot T ((mem_cycAnchored_iff T p1).mp a1).1
        ((mem_cycAnchored_iff T p2).mp a2).1 heq' with h | h | h
      · exact h.symm
      · exfalso; rw [h] at a2; exact anchored_unique₁ T p1 a1 a2
      · exfalso; rw [h] at a2; exact anchored_unique₂ T p1 a1 a2
    · exact absurd heq (cyc_trans_supp_ne T ((mem_cycAnchored_iff T p1).mp a1).1 a2)
    · exact absurd heq.symm (cyc_trans_supp_ne T ((mem_cycAnchored_iff T p2).mp a2).1 a1)
    · have heq' : supp p2 = supp p1 := heq.symm
      exact trans_supp_eq T a2 a1 heq'
  · intro s hs
    rw [Finset.mem_powersetCard] at hs
    obtain ⟨a, b, c, hab, hac, hbc, rfl⟩ := Finset.card_eq_three.mp hs.2
    obtain ⟨p, hp, hsupp⟩ := exists_rep T hab hac hbc
    exact ⟨p, hp, hsupp⟩

/-- The out-pair count: #orderedTrans = Σ_v C(s_v, 2). -/
theorem orderedTrans_card_eq_sum_choose (T : Tournament n) :
    (T.orderedTrans).card = ∑ v : Fin n, (T.outDegree v).choose 2 := by
  classical
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun p : Fin n × Fin n × Fin n => p.1) (t := univ) (fun p _ => mem_univ _)]
  refine Finset.sum_congr rfl (fun v _ => ?_)
  -- fiber over v ≃ 2-subsets of the out-set of v
  have : (T.orderedTrans.filter (fun p => p.1 = v)).card
      = ((univ.filter (fun w => T.arc v w = true)).powersetCard 2).card := by
    refine Finset.card_bij (fun p _ => ({p.2.1, p.2.2} : Finset (Fin n))) ?_ ?_ ?_
    · intro p hp
      simp only [Finset.mem_filter] at hp
      obtain ⟨hpt, hpv⟩ := hp
      have t1 := ((mem_orderedTrans_iff T p).mp hpt).1
      have t2 := ((mem_orderedTrans_iff T p).mp hpt).2.1
      have t3 := ((mem_orderedTrans_iff T p).mp hpt).2.2
      obtain ⟨d1, d2, d3⟩ := mem_orderedTrans_distinct T hpt
      rw [Finset.mem_powersetCard]
      constructor
      · intro y hy
        simp only [Finset.mem_insert, Finset.mem_singleton] at hy
        rcases hy with rfl | rfl
        · exact Finset.mem_filter.mpr ⟨mem_univ _, hpv ▸ t1⟩
        · exact Finset.mem_filter.mpr ⟨mem_univ _, hpv ▸ t3⟩
      · rw [Finset.card_insert_of_notMem (by simp [d2]), Finset.card_singleton]
    · intro p1 h1 p2 h2 heq
      simp only [Finset.mem_filter] at h1 h2
      obtain ⟨hp1, hv1⟩ := h1
      obtain ⟨hp2, hv2⟩ := h2
      have s1 := ((mem_orderedTrans_iff T p1).mp hp1).2.1
      have s2 := ((mem_orderedTrans_iff T p2).mp hp2).2.1
      have heqb : ({p1.2.1, p1.2.2} : Finset (Fin n)) = {p2.2.1, p2.2.2} := heq
      have hm : p2.2.1 = p1.2.1 ∨ p2.2.1 = p1.2.2 := by
        have hin0 : p2.2.1 ∈ ({p1.2.1, p1.2.2} : Finset (Fin n)) := by
          rw [heqb]; simp
        simpa using hin0
      rcases hm with m | m
      · have hm2 : p2.2.2 = p1.2.2 := by
          have hin : p2.2.2 ∈ ({p1.2.1, p1.2.2} : Finset (Fin n)) := by
            rw [heqb]; simp
          simp only [Finset.mem_insert, Finset.mem_singleton] at hin
          rcases hin with m2 | m2
          · exact absurd (m.trans m2.symm) (mem_orderedTrans_distinct T hp2).2.1
          · exact m2
        have h11 : p2.1 = p1.1 := hv2.trans hv1.symm
        exact (Prod.ext h11 (Prod.ext m hm2)).symm
      · exfalso
        have hm2 : p2.2.2 = p1.2.1 := by
          have hin : p2.2.2 ∈ ({p1.2.1, p1.2.2} : Finset (Fin n)) := by
            rw [heqb]; simp
          simp only [Finset.mem_insert, Finset.mem_singleton] at hin
          rcases hin with m2 | m2
          · exact m2
          · exact absurd (m.trans m2.symm) (mem_orderedTrans_distinct T hp2).2.1
        rw [m, hm2] at s2
        exact T.asym p1.2.1 p1.2.2 ⟨s1, s2⟩
    · intro s hs
      rw [Finset.mem_powersetCard] at hs
      obtain ⟨x, y, hxy, rfl⟩ := Finset.card_eq_two.mp hs.2
      have hx : T.arc v x = true := by
        have := hs.1 (by simp : x ∈ ({x, y} : Finset (Fin n)))
        exact (Finset.mem_filter.mp this).2
      have hy : T.arc v y = true := by
        have := hs.1 (by simp : y ∈ ({x, y} : Finset (Fin n)))
        exact (Finset.mem_filter.mp this).2
      rcases T.total x y hxy with hxy' | hxy'
      · refine ⟨(v, x, y), Finset.mem_filter.mpr
          ⟨(mem_orderedTrans_iff T _).mpr ⟨hx, hxy', hy⟩, rfl⟩, rfl⟩
      · refine ⟨(v, y, x), Finset.mem_filter.mpr
          ⟨(mem_orderedTrans_iff T _).mpr ⟨hy, hxy', hx⟩, rfl⟩, ?_⟩
        ext z; simp; tauto
  rw [this, Finset.card_powersetCard]
  rfl

/-- KENDALL–BABINGTON SMITH: c3 + Σ_v C(s_v, 2) = C(n, 3). -/
theorem kendall (T : Tournament n) :
    (T.cycAnchored).card + ∑ v : Fin n, (T.outDegree v).choose 2 = n.choose 3 := by
  rw [← orderedTrans_card_eq_sum_choose]
  exact cycAnchored_card_add_orderedTrans_card T

end Tournament
