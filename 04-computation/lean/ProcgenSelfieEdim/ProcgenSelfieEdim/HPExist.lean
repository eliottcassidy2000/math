import ProcgenSelfieEdim.Circulant

set_option autoImplicit false

/-!
# Existence of Hamiltonian paths, and HP-covered tournaments (THM-4524 §3)

* `exists_path_on`: on any duplicate-free vertex list `W` on which `T` is a tournament
  there is a path of `T` through exactly the vertices of `W` (insertion).
* `exists_hamPath`, `hpCount_pos`: every tournament on `n` vertices has a Hamiltonian
  path (the existence half of Rédei's theorem).
* `no_covered_universal_arc` (THM-4524, "no tournament on `n ≥ 3` vertices has every arc on
  some HP and some arc on every HP").
-/

namespace ProcgenSelfieEdim

/-- Insert `w` into the path `Q`, before the first vertex that `w` beats. -/
def insertV (T : Nat → Nat → Bool) (w : Nat) : List Nat → List Nat
  | [] => [w]
  | q :: Q => if T w q = true then w :: q :: Q else q :: insertV T w Q

theorem mem_insertV (T : Nat → Nat → Bool) (w x : Nat) :
    ∀ (Q : List Nat), x ∈ insertV T w Q ↔ x = w ∨ x ∈ Q
  | [] => by
    show x ∈ [w] ↔ _
    rw [List.mem_singleton]
    exact ⟨Or.inl, fun h => h.elim id (fun h' => absurd h' List.not_mem_nil)⟩
  | q :: Q => by
    unfold insertV
    by_cases h : T w q = true
    · rw [if_pos h, List.mem_cons]
    · rw [if_neg h, List.mem_cons, mem_insertV T w x Q, List.mem_cons]
      constructor
      · rintro (h1 | h1 | h1)
        · exact Or.inr (Or.inl h1)
        · exact Or.inl h1
        · exact Or.inr (Or.inr h1)
      · rintro (h1 | h1 | h1)
        · exact Or.inr (Or.inl h1)
        · exact Or.inl h1
        · exact Or.inr (Or.inr h1)

theorem nodup_insertV (T : Nat → Nat → Bool) (w : Nat) :
    ∀ (Q : List Nat), Q.Nodup → w ∉ Q → (insertV T w Q).Nodup
  | [], _, _ => List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩
  | q :: Q, hnd, hw => by
    unfold insertV
    by_cases h : T w q = true
    · rw [if_pos h]
      exact List.nodup_cons.2 ⟨hw, hnd⟩
    · rw [if_neg h]
      rw [List.nodup_cons] at hnd
      refine List.nodup_cons.2 ⟨fun hq => ?_, nodup_insertV T w Q hnd.2
        (fun h' => hw (List.mem_cons_of_mem q h'))⟩
      rcases (mem_insertV T w q Q).1 hq with h' | h'
      · exact hw (h' ▸ List.mem_cons_self)
      · exact hnd.1 h'

/-- The first vertex of `insertV T w Q` is `w` or the first vertex of `Q`. -/
theorem head_insertV (T : Nat → Nat → Bool) (w : Nat) :
    ∀ (Q : List Nat), ∃ rest, insertV T w Q = w :: rest ∨ ∃ q Q', Q = q :: Q' ∧ insertV T w Q = q :: rest
  | [] => ⟨[], Or.inl rfl⟩
  | q :: Q => by
    unfold insertV
    by_cases h : T w q = true
    · rw [if_pos h]
      exact ⟨q :: Q, Or.inl rfl⟩
    · rw [if_neg h]
      exact ⟨insertV T w Q, Or.inr ⟨q, Q, rfl, rfl⟩⟩

theorem isPath_insertV (T : Nat → Nat → Bool) (w : Nat) :
    ∀ (Q : List Nat), IsPath T Q → (∀ q ∈ Q, T w q = false → T q w = true) →
      IsPath T (insertV T w Q)
  | [], _, _ => trivial
  | q :: Q, hp, hw => by
    unfold insertV
    by_cases h : T w q = true
    · rw [if_pos h]
      exact ⟨h, hp⟩
    · rw [if_neg h]
      have hqw : T q w = true := hw q List.mem_cons_self (by cases hc : T w q <;> simp_all)
      have hrest : IsPath T (insertV T w Q) := isPath_insertV T w Q
        (by cases Q with
          | nil => trivial
          | cons _ _ => exact hp.2)
        (fun x hx => hw x (List.mem_cons_of_mem q hx))
      obtain ⟨rest, hr | ⟨q', Q', hQ, hr⟩⟩ := head_insertV T w Q
      · rw [hr] at hrest ⊢
        exact ⟨hqw, hrest⟩
      · rw [hr] at hrest ⊢
        subst hQ
        exact ⟨hp.1, hrest⟩

/-- **Paths through any vertex set of a tournament** (by insertion). -/
theorem exists_path_on (T : Nat → Nat → Bool) :
    ∀ (W : List Nat), W.Nodup → (∀ a ∈ W, ∀ b ∈ W, a ≠ b → T b a = !T a b) →
      ∃ Q, IsPath T Q ∧ Q.Nodup ∧ ∀ x, x ∈ Q ↔ x ∈ W
  | [], _, _ => ⟨[], trivial, List.Pairwise.nil, fun _ => Iff.rfl⟩
  | w :: W, hnd, hT => by
    rw [List.nodup_cons] at hnd
    obtain ⟨Q, hpath, hQnd, hmem⟩ := exists_path_on T W hnd.2
      (fun a ha b hb hab => hT a (List.mem_cons_of_mem w ha) b (List.mem_cons_of_mem w hb) hab)
    have hwQ : w ∉ Q := fun h => hnd.1 ((hmem w).1 h)
    refine ⟨insertV T w Q, isPath_insertV T w Q hpath ?_, nodup_insertV T w Q hQnd hwQ, ?_⟩
    · intro q hq hwq
      have hqW := (hmem q).1 hq
      have hne : w ≠ q := fun h => hwQ (h ▸ hq)
      rw [hT w List.mem_cons_self q (List.mem_cons_of_mem w hqW) hne, hwq]
      rfl
    · intro x
      rw [mem_insertV, List.mem_cons, hmem]

/-- **Every tournament has a Hamiltonian path.** -/
theorem exists_hamPath (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    ∃ p, IsHamPath n T p := by
  obtain ⟨Q, hpath, hnd, hmem⟩ := exists_path_on T (downR n) (nodup_downR n)
    (fun a ha b hb hab => hT a b ((mem_downR n a).1 ha) ((mem_downR n b).1 hb) hab)
  refine ⟨Q, ?_, hnd, fun x hx => (mem_downR n x).1 ((hmem x).1 hx), hpath⟩
  rw [← length_downR n]
  exact length_eq_of_same_members Q (downR n) hnd (nodup_downR n) hmem

theorem hpCount_pos (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) : 0 < hpCount n T := by
  obtain ⟨p, hp⟩ := exists_hamPath n T hT
  exact List.length_pos_of_mem ((mem_hamPathsR n T p).2 hp)

/-! ## No HP-covered tournament has an arc on every HP -/

/-- In a duplicate-free list a vertex has at most one successor. -/
theorem usesArc_succ_unique (u w v : Nat) (p : List Nat) (hnd : p.Nodup)
    (h1 : usesArc u w p = true) (h2 : usesArc u v p = true) : w = v := by
  obtain ⟨l1, l2, rfl⟩ := (usesArc_iff u w p).1 h1
  obtain ⟨l1', l2', h⟩ := (usesArc_iff u v _).1 h2
  have hl := pos_unique u l1 (w :: l2) l1' (v :: l2') hnd h
  have := (List.append_inj h hl).2
  exact (List.cons.inj (List.cons.inj this).2).1

/-- In a duplicate-free list a vertex has at most one predecessor. -/
theorem usesArc_pred_unique (u w v : Nat) (p : List Nat) (hnd : p.Nodup)
    (h1 : usesArc w v p = true) (h2 : usesArc u v p = true) : w = u := by
  obtain ⟨l1, l2, rfl⟩ := (usesArc_iff w v p).1 h1
  obtain ⟨l1', l2', h⟩ := (usesArc_iff u v _).1 h2
  have h' : (l1 ++ [w]) ++ v :: l2 = (l1' ++ [u]) ++ v :: l2' := by
    simp only [List.append_assoc, List.cons_append, List.nil_append]
    exact h
  have hnd' : ((l1 ++ [w]) ++ v :: l2).Nodup := by
    simp only [List.append_assoc, List.cons_append, List.nil_append]
    exact hnd
  have hl := pos_unique v (l1 ++ [w]) l2 (l1' ++ [u]) l2' hnd' h'
  have := (List.append_inj h' hl).1
  have hl2 : l1.length = l1'.length := by
    simp only [List.length_append, List.length_singleton] at hl
    omega
  exact (List.cons.inj (List.append_inj this hl2).2).1

/-- The last vertex of a duplicate-free list is never followed by anything. -/
theorem usesArc_last (u v : Nat) (Q : List Nat) (hQ : Q.Nodup) (hu : u ∉ Q) :
    usesArc u v (Q ++ [u]) = false := by
  cases h : usesArc u v (Q ++ [u])
  · rfl
  · exfalso
    obtain ⟨l1, l2, hl⟩ := (usesArc_iff u v _).1 h
    have hnd : (Q ++ u :: []).Nodup := by
      rw [List.nodup_append]
      refine ⟨hQ, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩, ?_⟩
      intro x hx y hy hxy
      rw [List.mem_singleton] at hy
      subst hy
      subst hxy
      exact hu hx
    have hpos := pos_unique u Q [] l1 (v :: l2) hnd hl
    have hlen := congrArg List.length hl
    rw [List.length_append, List.length_singleton, length_append_two] at hlen
    omega

/-- **THM-4524 (§3).** No tournament on `n ≥ 3` vertices has every arc on some Hamiltonian
path and some arc on every Hamiltonian path. -/
theorem no_covered_universal_arc (n : Nat) (hn : 3 ≤ n) (T : Nat → Nat → Bool)
    (hT : IsTournament n T)
    (hcover : ∀ a b, a < n → b < n → a ≠ b → T a b = true →
      ∃ p, IsHamPath n T p ∧ usesArc a b p = true) :
    ¬ ∃ u v, u < n ∧ v < n ∧ u ≠ v ∧ T u v = true ∧
      ∀ p, IsHamPath n T p → usesArc u v p = true := by
  rintro ⟨u, v, hu, hv, huv, harc, hall⟩
  -- `v` is `u`'s only out-neighbour and `u` is `v`'s only in-neighbour
  have hout : ∀ w, w < n → w ≠ u → w ≠ v → T w u = true ∧ T v w = true := by
    intro w hw hwu hwv
    constructor
    · cases h : T u w
      · rw [hT u w hu hw (Ne.symm hwu), h]
        rfl
      · obtain ⟨p, hp, hpu⟩ := hcover u w hu hw (Ne.symm hwu) h
        exact absurd (usesArc_succ_unique u w v p hp.2.1 hpu (hall p hp)) hwv
    · cases h : T w v
      · rw [hT w v hw hv hwv, h]
        rfl
      · obtain ⟨p, hp, hpv⟩ := hcover w v hw hv hwv h
        exact absurd (usesArc_pred_unique u w v p hp.2.1 hpv (hall p hp)) hwu
  -- a path through all vertices but `u`
  let W := (downR n).filter (fun x => decide (x ≠ u))
  have hWmem : ∀ x, x ∈ W ↔ x < n ∧ x ≠ u := by
    intro x
    show x ∈ (downR n).filter _ ↔ _
    rw [List.mem_filter, mem_downR, decide_eq_true_iff]
  obtain ⟨Q, hpath, hQnd, hQmem⟩ := exists_path_on T W
    (List.Nodup.sublist List.filter_sublist (nodup_downR n))
    (fun a ha b hb hab => hT a b ((hWmem a).1 ha).1 ((hWmem b).1 hb).1 hab)
  have hQlen : Q.length + 1 = n := by
    have h1 := length_eq_of_same_members Q W hQnd
      (List.Nodup.sublist List.filter_sublist (nodup_downR n)) hQmem
    have h2 := length_filter_ne u (downR n) ((mem_downR n u).2 hu) (nodup_downR n)
    rw [length_downR] at h2
    show Q.length + 1 = n
    rw [h1]
    exact h2
  -- `Q` starts at `v` (a source of `T - u`) and ends at some `w ≠ v`
  have hvQ : v ∈ Q := (hQmem v).2 ((hWmem v).2 ⟨hv, Ne.symm huv⟩)
  obtain ⟨A, B, hAB⟩ := List.append_of_mem hvQ
  have hA : A = [] := by
    rcases List.eq_nil_or_concat A with h0 | ⟨A', p, h0⟩
    · exact h0
    · exfalso
      rw [List.concat_eq_append] at h0
      have hQ' : Q = A' ++ p :: v :: B := by
        rw [hAB, h0, List.append_assoc]
        rfl
      have hpQ : p ∈ Q := by
        rw [hQ']
        exact List.mem_append_right _ List.mem_cons_self
      have hpv : T p v = true := (isPath_iff_pairs T Q).1 hpath A' B p v hQ'
      have hpW := (hWmem p).1 ((hQmem p).1 hpQ)
      have hpne : p ≠ v := by
        have hnd2 : (p :: v :: B).Nodup :=
          List.Nodup.sublist (List.sublist_append_right A' _) (hQ' ▸ hQnd)
        exact fun h => (List.nodup_cons.1 hnd2).1 (h ▸ List.mem_cons_self)
      have hvp := (hout p hpW.1 hpW.2 hpne).2
      rw [hT v p hv hpW.1 (Ne.symm hpne), hvp] at hpv
      exact Bool.noConfusion hpv
  subst hA
  -- `B` is nonempty; its last vertex beats `u`
  have hB : B ≠ [] := by
    intro hB
    subst hB
    rw [hAB] at hQlen
    simp only [List.nil_append, List.length_singleton] at hQlen
    omega
  obtain ⟨mid, hmid⟩ := eq_append_lastOr B v hB
  have hlast : lastOr v B ∈ Q := by
    rw [hAB]
    exact List.mem_cons_of_mem v (lastOr_mem B v hB)
  have hlastW := (hWmem _).1 ((hQmem _).1 hlast)
  have hlastv : lastOr v B ≠ v := by
    intro h
    have := hQnd
    rw [hAB, List.nil_append, List.nodup_cons] at this
    exact this.1 (h ▸ lastOr_mem B v hB)
  have hwu := (hout _ hlastW.1 hlastW.2 hlastv).1
  -- the Hamiltonian path `Q ++ [u]` ends at `u`, so it does not use `u → v`
  have hP : IsHamPath n T (Q ++ [u]) := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · rw [List.length_append, List.length_singleton, hQlen]
    · rw [List.nodup_append]
      refine ⟨hQnd, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩, ?_⟩
      intro x hx y hy hxy
      rw [List.mem_singleton] at hy
      subst hy
      subst hxy
      exact ((hWmem x).1 ((hQmem x).1 hx)).2 rfl
    · intro x hx
      rw [List.mem_append, List.mem_singleton] at hx
      rcases hx with hx | rfl
      · exact ((hWmem x).1 ((hQmem x).1 hx)).1
      · exact hu
    · have hpath' : IsPath T (v :: B) := by
        rw [hAB, List.nil_append] at hpath
        exact hpath
      rw [hAB, List.nil_append]
      exact isPath_snoc T v B u hpath' hwu
  have hnot := usesArc_last u v Q hQnd (fun h => ((hWmem u).1 ((hQmem u).1 h)).2 rfl)
  rw [hall _ hP] at hnot
  exact Bool.noConfusion hnot

end ProcgenSelfieEdim
