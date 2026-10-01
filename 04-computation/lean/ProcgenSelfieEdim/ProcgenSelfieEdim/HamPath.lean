import ProcgenSelfieEdim.ListLemmas
import ProcgenSelfieEdim.Tournament

set_option autoImplicit false

/-!
# Hamiltonian paths: a verified enumerator and arc counts

A Hamiltonian path (HP) of `T` on `{0, …, n-1}` is a list containing every vertex
exactly once, with consecutive vertices joined by arcs of `T` (`IsHamPath`). The
enumerator `hamPathsR` is a depth-first search that grows a path backwards, only
along arcs of `T`; `mem_hamPathsR` and `nodup_hamPathsR` show that it lists every
HP exactly once. Counts are lengths of such lists:

* `hpCount n T = H(T)`, the number of HPs;
* `arcCount n T u v = c(u → v)`, the number of HPs using the arc `u → v`.

`hpCount_unique` / `arcCount_unique` show these numbers do not depend on the
enumerator: any duplicate-free list of exactly those paths has the same length.
-/

namespace ProcgenSelfieEdim

/-! ## Definitions -/

/-- Consecutive vertices of `p` are joined by arcs of `T`. -/
def IsPath (T : Nat → Nat → Bool) : List Nat → Prop
  | a :: b :: rest => T a b = true ∧ IsPath T (b :: rest)
  | _ => True

/-- `p` is a Hamiltonian path of `T` on `{0, …, n-1}`. -/
def IsHamPath (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) : Prop :=
  p.length = n ∧ p.Nodup ∧ (∀ x ∈ p, x < n) ∧ IsPath T p

/-- `usesArc u v p = true` iff `u` is immediately followed by `v` in `p`
(see `usesArc_iff`). -/
def usesArc (u v : Nat) : List Nat → Bool
  | a :: b :: rest => (Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)
  | _ => false

/-- A duplicate-free list of vertices below `n` that is a path of `T`. -/
def Valid (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) : Prop :=
  p.Nodup ∧ (∀ x ∈ p, x < n) ∧ IsPath T p

theorem isPath_cons_cons (T : Nat → Nat → Bool) (a b : Nat) (rest : List Nat) :
    IsPath T (a :: b :: rest) ↔ T a b = true ∧ IsPath T (b :: rest) := Iff.rfl

theorem isPath_append_right (T : Nat → Nat → Bool) :
    ∀ (l1 l2 : List Nat), IsPath T (l1 ++ l2) → IsPath T l2
  | [], _, h => h
  | [_], [], _ => trivial
  | [_], _ :: _, h => h.2
  | _ :: y :: t, l2, h => isPath_append_right T (y :: t) l2 h.2

theorem length_append_two (l1 l2 : List Nat) (u v : Nat) :
    (l1 ++ u :: v :: l2).length = l1.length + l2.length + 2 := by
  rw [List.length_append, List.length_cons, List.length_cons]
  omega

theorem usesArc_iff (u v : Nat) :
    ∀ (p : List Nat), usesArc u v p = true ↔ ∃ l1 l2, p = l1 ++ u :: v :: l2
  | [] => ⟨fun h => Bool.noConfusion h, fun ⟨l1, l2, h⟩ => by
      have := congrArg List.length h
      rw [length_append_two, List.length_nil] at this
      omega⟩
  | [a] => ⟨fun h => Bool.noConfusion h, fun ⟨l1, l2, h⟩ => by
      have := congrArg List.length h
      rw [length_append_two] at this
      simp only [List.length_singleton] at this
      omega⟩
  | a :: b :: rest => by
    show ((Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)) = true ↔ _
    rw [Bool.or_eq_true, Bool.and_eq_true, Nat.beq_eq, Nat.beq_eq, usesArc_iff u v (b :: rest)]
    constructor
    · rintro (⟨rfl, rfl⟩ | ⟨l1, l2, h⟩)
      · exact ⟨[], rest, rfl⟩
      · exact ⟨a :: l1, l2, by rw [h]; rfl⟩
    · rintro ⟨l1, l2, h⟩
      cases l1 with
      | nil =>
        rw [List.nil_append, List.cons.injEq, List.cons.injEq] at h
        exact Or.inl ⟨h.1, h.2.1⟩
      | cons x t =>
        rw [List.cons_append, List.cons.injEq] at h
        exact Or.inr ⟨t, l2, h.2⟩

theorem usesArc_mem (u v : Nat) :
    ∀ (p : List Nat), usesArc u v p = true → u ∈ p ∧ v ∈ p := by
  intro p h
  obtain ⟨l1, l2, rfl⟩ := (usesArc_iff u v p).1 h
  exact ⟨List.mem_append_right _ List.mem_cons_self,
    List.mem_append_right _ (List.mem_cons_of_mem _ List.mem_cons_self)⟩

/-- A path only uses arcs of `T`. -/
theorem usesArc_arc (T : Nat → Nat → Bool) (u v : Nat) :
    ∀ (p : List Nat), IsPath T p → usesArc u v p = true → T u v = true
  | [], _, h => Bool.noConfusion h
  | [_], _, h => Bool.noConfusion h
  | a :: b :: rest, hp, h => by
    have h' : ((Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)) = true := h
    rw [Bool.or_eq_true, Bool.and_eq_true, Nat.beq_eq, Nat.beq_eq] at h'
    rcases h' with ⟨rfl, rfl⟩ | h''
    · exact hp.1
    · exact usesArc_arc T u v (b :: rest) hp.2 h''

/-- A duplicate-free list never uses a loop `u → u`. -/
theorem usesArc_self (u : Nat) (p : List Nat) (hp : p.Nodup) : usesArc u u p = false := by
  cases h : usesArc u u p
  · rfl
  · obtain ⟨l1, l2, rfl⟩ := (usesArc_iff u u p).1 h
    have := (List.nodup_cons.1 (List.Nodup.sublist (List.sublist_append_right l1 _) hp)).1
    exact absurd List.mem_cons_self this

/-! ## The enumerator -/

/-- Depth-first extension backwards: all ways to prepend `k` further vertices to `p`,
each an in-neighbour (in `T`) of the current first vertex and not yet used. -/
noncomputable def extR (n : Nat) (T : Nat → Nat → Bool) (k : Nat) : List Nat → List (List Nat) :=
  Nat.rec (motive := fun _ => List Nat → List (List Nat))
    (fun p => [p])
    (fun _ ih p => List.casesOn (motive := fun _ => List (List Nat)) p []
      (fun h _ => flatR (fun w => cond (T w h && !memR w p) (ih (w :: p)) []) (downR n)))
    k

theorem extR_zero (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) : extR n T 0 p = [p] := rfl

theorem extR_succ_cons (n : Nat) (T : Nat → Nat → Bool) (k h : Nat) (rest : List Nat) :
    extR n T (k + 1) (h :: rest) =
      flatR (fun w => cond (T w h && !memR w (h :: rest)) (extR n T k (w :: h :: rest)) [])
        (downR n) := rfl

/-- All Hamiltonian paths of `T` on `{0, …, n-1}`. -/
noncomputable def hamPathsR (n : Nat) (T : Nat → Nat → Bool) : List (List Nat) :=
  Nat.casesOn (motive := fun _ => List (List Nat)) n [[]]
    (fun m => flatR (fun v => extR (m + 1) T m [v]) (downR (m + 1)))

theorem hamPathsR_zero (T : Nat → Nat → Bool) : hamPathsR 0 T = [[]] := rfl

theorem hamPathsR_succ (m : Nat) (T : Nat → Nat → Bool) :
    hamPathsR (m + 1) T = flatR (fun v => extR (m + 1) T m [v]) (downR (m + 1)) := rfl

/-- Every list produced by `extR n T k p` is `p` with `k` vertices prepended. -/
theorem extR_shape (n : Nat) (T : Nat → Nat → Bool) :
    ∀ (k : Nat) (p q : List Nat), q ∈ extR n T k p → ∃ e, q = e ++ p ∧ e.length = k
  | 0, p, q, hq => by
    rw [extR_zero, List.mem_singleton] at hq
    exact ⟨[], hq, rfl⟩
  | _ + 1, [], _, hq => absurd hq List.not_mem_nil
  | k + 1, h :: rest, q, hq => by
    rw [extR_succ_cons, flatR_eq, List.mem_flatMap] at hq
    obtain ⟨w, _, hw⟩ := hq
    cases hc : (T w h && !memR w (h :: rest))
    · rw [hc, cond_false] at hw
      exact absurd hw List.not_mem_nil
    · rw [hc, cond_true] at hw
      obtain ⟨e, he, hl⟩ := extR_shape n T k (w :: h :: rest) q hw
      refine ⟨e ++ [w], ?_, ?_⟩
      · rw [he, List.append_assoc]
        rfl
      · rw [List.length_append, hl]
        rfl

theorem valid_cons (n : Nat) (T : Nat → Nat → Bool) (w h : Nat) (rest : List Nat)
    (hp : Valid n T (h :: rest)) (hw : w < n) (hT : T w h = true) (hwp : w ∉ h :: rest) :
    Valid n T (w :: h :: rest) := by
  refine ⟨List.nodup_cons.2 ⟨hwp, hp.1⟩, ?_, ⟨hT, hp.2.2⟩⟩
  intro x hx
  rcases List.mem_cons.1 hx with rfl | hx'
  · exact hw
  · exact hp.2.1 x hx'

/-- Specification of `extR`: it lists exactly the valid extensions of `p` by `k`
prepended vertices. -/
theorem mem_extR (n : Nat) (T : Nat → Nat → Bool) :
    ∀ (k : Nat) (p : List Nat), p ≠ [] → Valid n T p →
      ∀ q, q ∈ extR n T k p ↔ (∃ e, q = e ++ p ∧ e.length = k) ∧ Valid n T q
  | 0, p, _, hp, q => by
    rw [extR_zero, List.mem_singleton]
    constructor
    · rintro rfl
      exact ⟨⟨[], rfl, rfl⟩, hp⟩
    · rintro ⟨⟨e, rfl, he⟩, _⟩
      cases e with
      | nil => rfl
      | cons _ _ =>
        rw [List.length_cons] at he
        omega
  | _ + 1, [], hne, _, _ => absurd rfl hne
  | k + 1, h :: rest, _, hp, q => by
    rw [extR_succ_cons, flatR_eq, List.mem_flatMap]
    constructor
    · rintro ⟨w, hw, hq⟩
      rw [mem_downR] at hw
      cases hc : (T w h && !memR w (h :: rest))
      · rw [hc, cond_false] at hq
        exact absurd hq List.not_mem_nil
      · rw [hc, cond_true] at hq
        rw [Bool.and_eq_true, Bool.not_eq_true'] at hc
        have hwp : w ∉ h :: rest := fun hm => by
          rw [(memR_iff w _).2 hm] at hc
          exact Bool.noConfusion hc.2
        have hv := valid_cons n T w h rest hp hw hc.1 hwp
        obtain ⟨⟨e, he, hl⟩, hq'⟩ := (mem_extR n T k (w :: h :: rest) (List.cons_ne_nil _ _) hv q).1 hq
        refine ⟨⟨e ++ [w], ?_, ?_⟩, hq'⟩
        · rw [he, List.append_assoc]
          rfl
        · rw [List.length_append, hl]
          rfl
    · rintro ⟨⟨e, he, hl⟩, hq⟩
      rcases List.eq_nil_or_concat e with h0 | ⟨e0, w, h0⟩
      · subst h0
        rw [List.length_nil] at hl
        omega
      · rw [List.concat_eq_append] at h0
        subst h0
        have hq2 : q = e0 ++ (w :: h :: rest) := by
          rw [he, List.append_assoc]
          rfl
        have hwq : w ∈ q := by
          rw [hq2]
          exact List.mem_append_right _ List.mem_cons_self
        have hw : w < n := hq.2.1 w hwq
        have hT : T w h = true := (isPath_append_right T e0 (w :: h :: rest) (hq2 ▸ hq.2.2)).1
        have hnd : (w :: h :: rest).Nodup :=
          List.Nodup.sublist (List.sublist_append_right e0 _) (hq2 ▸ hq.1)
        have hwp : w ∉ h :: rest := (List.nodup_cons.1 hnd).1
        have hc : (T w h && !memR w (h :: rest)) = true := by
          rw [hT, Bool.true_and]
          cases hm : memR w (h :: rest)
          · rfl
          · exact absurd ((memR_iff w _).1 hm) hwp
        refine ⟨w, (mem_downR n w).2 hw, ?_⟩
        rw [hc, cond_true]
        have hv := valid_cons n T w h rest hp hw hT hwp
        refine (mem_extR n T k (w :: h :: rest) (List.cons_ne_nil _ _) hv q).2 ⟨⟨e0, hq2, ?_⟩, hq⟩
        rw [List.length_append] at hl
        simp only [List.length_singleton] at hl
        omega

/-- `extR` produces no repetitions. -/
theorem nodup_extR (n : Nat) (T : Nat → Nat → Bool) :
    ∀ (k : Nat) (p : List Nat), (extR n T k p).Nodup
  | 0, p => by
    rw [extR_zero]
    exact List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩
  | _ + 1, [] => List.Pairwise.nil
  | k + 1, h :: rest => by
    rw [extR_succ_cons, flatR_eq]
    apply nodup_flatMap _ _ (nodup_downR n)
    · intro w _
      cases (T w h && !memR w (h :: rest))
      · exact List.Pairwise.nil
      · exact nodup_extR n T k (w :: h :: rest)
    · intro w w' _ _ hww' q hq hq'
      cases hc : (T w h && !memR w (h :: rest))
      · rw [hc, cond_false] at hq
        exact absurd hq List.not_mem_nil
      · cases hc' : (T w' h && !memR w' (h :: rest))
        · rw [hc', cond_false] at hq'
          exact absurd hq' List.not_mem_nil
        · rw [hc, cond_true] at hq
          rw [hc', cond_true] at hq'
          obtain ⟨e, he, hl⟩ := extR_shape n T k _ q hq
          obtain ⟨e', he', hl'⟩ := extR_shape n T k _ q hq'
          have := List.append_inj (he.symm.trans he') (hl.trans hl'.symm)
          exact hww' (List.cons.inj this.2).1

/-- **Specification of the enumerator.** `hamPathsR n T` lists exactly the
Hamiltonian paths of `T` on `{0, …, n-1}`. -/
theorem mem_hamPathsR (n : Nat) (T : Nat → Nat → Bool) (q : List Nat) :
    q ∈ hamPathsR n T ↔ IsHamPath n T q := by
  cases n with
  | zero =>
    rw [hamPathsR_zero, List.mem_singleton]
    constructor
    · rintro rfl
      exact ⟨rfl, List.Pairwise.nil, fun x hx => absurd hx List.not_mem_nil, trivial⟩
    · rintro ⟨hl, -, -, -⟩
      cases q with
      | nil => rfl
      | cons _ _ =>
        rw [List.length_cons] at hl
        omega
  | succ m =>
    rw [hamPathsR_succ, flatR_eq, List.mem_flatMap]
    constructor
    · rintro ⟨v, hv, hq⟩
      rw [mem_downR] at hv
      have hval : Valid (m + 1) T [v] :=
        ⟨List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩,
          fun x hx => by rw [List.mem_singleton] at hx; omega, trivial⟩
      obtain ⟨⟨e, he, hl⟩, hq'⟩ :=
        (mem_extR (m + 1) T m [v] (List.cons_ne_nil _ _) hval q).1 hq
      refine ⟨?_, hq'.1, hq'.2.1, hq'.2.2⟩
      rw [he, List.length_append, hl]
      rfl
    · rintro ⟨hl, hnd, hlt, hpath⟩
      rcases List.eq_nil_or_concat q with h0 | ⟨e, v, h0⟩
      · subst h0
        rw [List.length_nil] at hl
        omega
      · rw [List.concat_eq_append] at h0
        subst h0
        have hv : v < m + 1 := hlt v (List.mem_append_right _ List.mem_cons_self)
        refine ⟨v, (mem_downR (m + 1) v).2 hv, ?_⟩
        have hval : Valid (m + 1) T [v] :=
          ⟨List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩,
            fun x hx => by rw [List.mem_singleton] at hx; omega, trivial⟩
        refine (mem_extR (m + 1) T m [v] (List.cons_ne_nil _ _) hval _).2
          ⟨⟨e, rfl, ?_⟩, hnd, hlt, hpath⟩
        rw [List.length_append] at hl
        simp only [List.length_singleton] at hl
        omega

/-- **The enumerator lists each Hamiltonian path once.** -/
theorem nodup_hamPathsR (n : Nat) (T : Nat → Nat → Bool) : (hamPathsR n T).Nodup := by
  cases n with
  | zero =>
    rw [hamPathsR_zero]
    exact List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩
  | succ m =>
    rw [hamPathsR_succ, flatR_eq]
    apply nodup_flatMap _ _ (nodup_downR (m + 1))
    · intro v _
      exact nodup_extR (m + 1) T m [v]
    · intro v v' _ _ hvv' q hq hq'
      obtain ⟨e, he, hl⟩ := extR_shape (m + 1) T m _ q hq
      obtain ⟨e', he', hl'⟩ := extR_shape (m + 1) T m _ q hq'
      have := List.append_inj (he.symm.trans he') (hl.trans hl'.symm)
      exact hvv' (List.cons.inj this.2).1

/-! ## Counts -/

/-- `H(T)`: the number of Hamiltonian paths of `T` on `{0, …, n-1}`. -/
noncomputable def hpCount (n : Nat) (T : Nat → Nat → Bool) : Nat := (hamPathsR n T).length

/-- `c(u → v)`: the number of Hamiltonian paths of `T` using the arc `u → v`. -/
noncomputable def arcCount (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) : Nat :=
  ((hamPathsR n T).filter (usesArc u v)).length

/-- `H(T)` is the length of any duplicate-free list of exactly the HPs. -/
theorem hpCount_unique (n : Nat) (T : Nat → Nat → Bool) (L : List (List Nat)) (hL : L.Nodup)
    (hmem : ∀ p, p ∈ L ↔ IsHamPath n T p) : L.length = hpCount n T :=
  length_eq_of_same_members L _ hL (nodup_hamPathsR n T)
    (fun p => (hmem p).trans (mem_hamPathsR n T p).symm)

/-- `c(u → v)` is the length of any duplicate-free list of exactly the HPs using `u → v`. -/
theorem arcCount_unique (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) (L : List (List Nat))
    (hL : L.Nodup) (hmem : ∀ p, p ∈ L ↔ IsHamPath n T p ∧ usesArc u v p = true) :
    L.length = arcCount n T u v := by
  apply length_eq_of_same_members L _ hL (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T))
  intro p
  rw [hmem p, List.mem_filter, mem_hamPathsR]

end ProcgenSelfieEdim
