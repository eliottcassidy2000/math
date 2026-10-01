import ProcgenSelfieEdim.Shaved

set_option autoImplicit false

/-!
# Copies of `H_n` and Hamiltonian cycles: `#copies(H_n) = H(T) - n · hc(T)`

For a tournament on `n ≥ 2` vertices, every Hamiltonian path either has its first vertex
beating its last (a copy of `H_n`; `copiesH`) or closes into a Hamiltonian cycle
(`closes`), and never both. A directed Hamiltonian cycle has exactly one rotation
starting at vertex `0`, so the number of Hamiltonian cycles is `hcCount`, the number of
closing Hamiltonian paths starting at `0`. Each cycle yields exactly `n` closing paths
(its `n` rotations; `closing_count`), so

  `copiesH n T + n · hcCount n T = H(T)`   (`copiesH_add`).

`Aut(H_n)` is trivial (S15 note), so copies and embeddings coincide.
-/

namespace ProcgenSelfieEdim

/-- The last vertex of `p` beats its first: `p` closes into a Hamiltonian cycle. -/
noncomputable def closes (T : Nat → Nat → Bool) (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a rest => T (lastOr a rest) a)

/-- `p` starts at vertex `0`. -/
noncomputable def startsAtZero (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a _ => Nat.beq a 0)

/-- The number of copies of `H_n`: Hamiltonian paths whose first vertex beats the last. -/
noncomputable def copiesH (n : Nat) (T : Nat → Nat → Bool) : Nat :=
  ((hamPathsR n T).filter (headBeatsLast T)).length

/-- The number of directed Hamiltonian cycles: closing Hamiltonian paths starting at `0`. -/
noncomputable def hcCount (n : Nat) (T : Nat → Nat → Bool) : Nat :=
  ((hamPathsR n T).filter (fun p => startsAtZero p && closes T p)).length

theorem lastOr_mem (l : List Nat) (a : Nat) (h : l ≠ []) : lastOr a l ∈ l := by
  obtain ⟨mid, hmid⟩ := eq_append_lastOr l a h
  have : lastOr a l ∈ mid ++ [lastOr a l] := List.mem_append_right _ List.mem_cons_self
  rw [← hmid] at this
  exact this

/-- `T` contains `H_n` iff it has at least one copy. -/
theorem containsH_iff_copiesH_pos (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) :
    ContainsH n T ↔ 0 < copiesH n T := by
  constructor
  · intro h
    obtain ⟨a, b, mid, hham, hab⟩ := h
    have hmem : (a :: (mid ++ [b])) ∈ (hamPathsR n T).filter (headBeatsLast T) := by
      rw [List.mem_filter, mem_hamPathsR]
      refine ⟨hham, ?_⟩
      show T a (lastOr a (mid ++ [b])) = true
      rw [lastOr_append]
      exact hab
    exact List.length_pos_of_mem hmem
  · intro h
    obtain ⟨p, hp⟩ : ∃ p, p ∈ (hamPathsR n T).filter (headBeatsLast T) := by
      unfold copiesH at h
      cases hl : (hamPathsR n T).filter (headBeatsLast T) with
      | nil =>
        rw [hl, List.length_nil] at h
        exact absurd h (Nat.lt_irrefl 0)
      | cons p _ => exact ⟨p, List.mem_cons_self⟩
    rw [List.mem_filter] at hp
    exact containsH_of_mem n hn T p hp.1 hp.2

/-! ## Every Hamiltonian path either closes or is a copy of `H_n` -/

theorem closes_eq_not (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (hT : IsTournament n T)
    (p : List Nat) (hp : IsHamPath n T p) : closes T p = !headBeatsLast T p := by
  cases p with
  | nil =>
    exfalso
    have := hp.1
    rw [List.length_nil] at this
    omega
  | cons a rest =>
    have hrest : rest ≠ [] := by
      intro hr
      subst hr
      have := hp.1
      rw [List.length_singleton] at this
      omega
    show T (lastOr a rest) a = !T a (lastOr a rest)
    have hb := lastOr_mem rest a hrest
    have hne : a ≠ lastOr a rest := fun h => (List.nodup_cons.1 hp.2.1).1 (h ▸ hb)
    exact hT a (lastOr a rest) (hp.2.2.1 a List.mem_cons_self)
      (hp.2.2.1 _ (List.mem_cons_of_mem a hb)) hne

theorem split_closes (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    ((hamPathsR n T).filter (closes T)).length + copiesH n T = hpCount n T := by
  have h := length_filter_add_not (headBeatsLast T) (hamPathsR n T)
  have hc : (hamPathsR n T).filter (closes T) =
      (hamPathsR n T).filter (fun p => !headBeatsLast T p) :=
    List.filter_congr (fun p hp => closes_eq_not n hn T hT p ((mem_hamPathsR n T p).1 hp))
  unfold copiesH hpCount
  rw [hc]
  omega

/-! ## Rotations -/

/-- Rotate left by one. -/
def rot1 : List Nat → List Nat
  | [] => []
  | a :: rest => rest ++ [a]

/-- Rotate left by `k`. -/
def rotN : Nat → List Nat → List Nat
  | 0, p => p
  | k + 1, p => rotN k (rot1 p)

theorem rotN_append : ∀ (A l : List Nat), rotN A.length (A ++ l) = l ++ A
  | [], l => by rw [List.length_nil, List.nil_append, List.append_nil]; rfl
  | a :: A, l => by
    rw [List.length_cons, List.cons_append]
    show rotN A.length ((A ++ l) ++ [a]) = l ++ a :: A
    rw [List.append_assoc, rotN_append A (l ++ [a]), List.append_assoc]
    rfl

theorem rotN_eq (c : List Nat) (k : Nat) (hk : k ≤ c.length) : rotN k c = c.drop k ++ c.take k := by
  have h := rotN_append (c.take k) (c.drop k)
  rw [List.take_append_drop, List.length_take, Nat.min_eq_left hk] at h
  exact h

theorem isPath_snoc (T : Nat → Nat → Bool) : ∀ (c : Nat) (l : List Nat) (x : Nat),
    IsPath T (c :: l) → T (lastOr c l) x = true → IsPath T (c :: l ++ [x])
  | _, [], _, _, h => ⟨h, trivial⟩
  | _, d :: l, x, hp, h => ⟨hp.1, isPath_snoc T d l x hp.2 h⟩

/-- A closing Hamiltonian path. -/
def IsClosingHP (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) : Prop :=
  IsHamPath n T p ∧ closes T p = true

theorem rot1_closing (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (p : List Nat)
    (h : IsClosingHP n T p) : IsClosingHP n T (rot1 p) := by
  obtain ⟨⟨hl, hnd, hlt, hpath⟩, hcl⟩ := h
  rcases p with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · exfalso
    rw [List.length_nil] at hl
    omega
  · exfalso
    rw [List.length_singleton] at hl
    omega
  · have hcl' : T (lastOr b rest) a = true := hcl
    rw [List.nodup_cons] at hnd
    show IsClosingHP n T (b :: rest ++ [a])
    refine ⟨⟨?_, ?_, ?_, isPath_snoc T b rest a hpath.2 hcl'⟩, ?_⟩
    · rw [← hl]
      simp only [List.length_cons, List.length_append, List.length_nil]
    · rw [List.nodup_append]
      refine ⟨hnd.2, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩, ?_⟩
      intro x hx y hy hxy
      rw [List.mem_singleton] at hy
      subst hy
      subst hxy
      exact hnd.1 hx
    · intro x hx
      rw [List.mem_append, List.mem_singleton] at hx
      rcases hx with hx | rfl
      · exact hlt x (List.mem_cons_of_mem a hx)
      · exact hlt _ List.mem_cons_self
    · show T (lastOr b (rest ++ [a])) b = true
      rw [lastOr_append]
      exact hpath.1

theorem rotN_closing (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) :
    ∀ (k : Nat) (p : List Nat), IsClosingHP n T p → IsClosingHP n T (rotN k p)
  | 0, _, h => h
  | k + 1, p, h => rotN_closing n hn T k (rot1 p) (rot1_closing n hn T p h)

/-- The position of an element in a duplicate-free list is unique. -/
theorem pos_unique (x : Nat) : ∀ (X Y X' Y' : List Nat), (X ++ x :: Y).Nodup →
    X ++ x :: Y = X' ++ x :: Y' → X.length = X'.length
  | [], _, [], _, _, _ => rfl
  | [], Y, z :: X', Y', hnd, h => by
    exfalso
    rw [List.nil_append, List.cons_append, List.cons.injEq] at h
    rw [List.nil_append, List.nodup_cons] at hnd
    apply hnd.1
    rw [h.2]
    exact List.mem_append_right _ List.mem_cons_self
  | z :: X, Y, [], Y', hnd, h => by
    exfalso
    rw [List.nil_append, List.cons_append, List.cons.injEq] at h
    rw [List.cons_append, List.nodup_cons] at hnd
    apply hnd.1
    rw [h.1]
    exact List.mem_append_right _ List.mem_cons_self
  | z :: X, Y, z' :: X', Y', hnd, h => by
    rw [List.cons_append, List.cons_append, List.cons.injEq] at h
    rw [List.cons_append, List.nodup_cons] at hnd
    rw [List.length_cons, List.length_cons, pos_unique x X Y X' Y' hnd.2 h.2]

theorem rotN_zero_head (B : List Nat) (k : Nat) (hk : 1 ≤ k) (hkB : k ≤ B.length) :
    rotN k (0 :: B) = B.drop (k - 1) ++ 0 :: B.take (k - 1) := by
  rw [rotN_eq (0 :: B) k (by rw [List.length_cons]; omega)]
  obtain ⟨j, rfl⟩ : ∃ j, k = j + 1 := ⟨k - 1, by omega⟩
  rw [List.drop_succ_cons, List.take_succ_cons, Nat.add_sub_cancel]

/-- Where `0` sits in the `k`-th rotation of a list starting at `0`. -/
theorem rotN_zero_split (n : Nat) (B : List Nat) (hB : B.length + 1 = n) (k : Nat) (hk : k < n) :
    ∃ X Y, rotN k (0 :: B) = X ++ 0 :: Y ∧ X.length = (if k = 0 then 0 else n - k) := by
  by_cases hk0 : k = 0
  · subst hk0
    exact ⟨[], B, rfl, by rw [if_pos rfl]; rfl⟩
  · refine ⟨B.drop (k - 1), B.take (k - 1), rotN_zero_head B k (by omega) (by omega), ?_⟩
    rw [if_neg hk0, List.length_drop]
    omega

theorem rotN_nodup : ∀ (k : Nat) (p : List Nat), p.Nodup → (rotN k p).Nodup
  | 0, _, h => h
  | k + 1, [], h => rotN_nodup k [] h
  | k + 1, a :: rest, h => by
    apply rotN_nodup k (rest ++ [a])
    rw [List.nodup_cons] at h
    rw [List.nodup_append]
    refine ⟨h.2, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩, ?_⟩
    intro x hx y hy hxy
    rw [List.mem_singleton] at hy
    subst hy
    subst hxy
    exact h.1 hx

/-- The rotation amount is determined by the position of `0`. -/
theorem rot_index_unique (n : Nat) (B B' : List Nat) (hB : B.length + 1 = n)
    (hB' : B'.length + 1 = n) (hnd : (0 :: B).Nodup) (k k' : Nat) (hk : k < n) (hk' : k' < n)
    (h : rotN k (0 :: B) = rotN k' (0 :: B')) : k = k' := by
  obtain ⟨X, Y, hX, hXl⟩ := rotN_zero_split n B hB k hk
  obtain ⟨X', Y', hX', hXl'⟩ := rotN_zero_split n B' hB' k' hk'
  have hnd' : (X ++ 0 :: Y).Nodup := hX ▸ rotN_nodup k _ hnd
  have := pos_unique 0 X Y X' Y' hnd' (hX.symm.trans (h.trans hX'))
  rw [hXl, hXl'] at this
  by_cases h0 : k = 0 <;> by_cases h0' : k' = 0
  · omega
  · rw [if_pos h0, if_neg h0'] at this
    omega
  · rw [if_neg h0, if_pos h0'] at this
    omega
  · rw [if_neg h0, if_neg h0'] at this
    omega

theorem rotN_inj (k : Nat) (c c' : List Nat) (hl : c.length = c'.length) (hk : k ≤ c.length)
    (h : rotN k c = rotN k c') : c = c' := by
  rw [rotN_eq c k hk, rotN_eq c' k (hl ▸ hk)] at h
  have hd : (c.drop k).length = (c'.drop k).length := by
    rw [List.length_drop, List.length_drop, hl]
  obtain ⟨h1, h2⟩ := List.append_inj h hd
  rw [← List.take_append_drop k c, ← List.take_append_drop k c', h1, h2]

/-! ## Counting the closing paths -/

theorem startsAtZero_iff (p : List Nat) : startsAtZero p = true ↔ ∃ B, p = 0 :: B := by
  cases p with
  | nil => exact ⟨fun h => Bool.noConfusion h, fun ⟨_, h⟩ => nomatch h⟩
  | cons a B =>
    show Nat.beq a 0 = true ↔ _
    rw [Nat.beq_eq]
    constructor
    · rintro rfl
      exact ⟨B, rfl⟩
    · rintro ⟨B', h⟩
      exact (List.cons.inj h).1

theorem zero_mem_of_hamPath (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) (p : List Nat)
    (hp : IsHamPath n T p) : 0 ∈ p := by
  by_cases h : 0 ∈ p
  · exact h
  · exfalso
    have h1 := pigeonhole p ((downR n).filter (fun x => decide (x ≠ 0))) hp.2.1 (fun x hx => by
      rw [List.mem_filter]
      exact ⟨(mem_downR n x).2 (hp.2.2.1 x hx), decide_eq_true (fun h0 => h (h0 ▸ hx))⟩)
    have h2 := length_filter_lt (fun x => decide (x ≠ 0)) 0 (downR n)
      ((mem_downR n 0).2 (by omega)) rfl
    rw [length_downR] at h2
    rw [hp.1] at h1
    omega

/-- **Each Hamiltonian cycle yields exactly `n` closing Hamiltonian paths**:
`#{closing HPs} = n · hc(T)` (any relation `T`, `n ≥ 2`). -/
theorem closing_count (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) :
    ((hamPathsR n T).filter (closes T)).length = n * hcCount n T := by
  let C0 := (hamPathsR n T).filter (fun p => startsAtZero p && closes T p)
  let L := C0.flatMap (fun c => (downR n).map (fun k => rotN k c))
  have hC0 : ∀ c, c ∈ C0 ↔ IsClosingHP n T c ∧ ∃ B, c = 0 :: B := by
    intro c
    show c ∈ (hamPathsR n T).filter _ ↔ _
    rw [List.mem_filter, mem_hamPathsR, Bool.and_eq_true, startsAtZero_iff]
    constructor
    · rintro ⟨h1, h2, h3⟩
      exact ⟨⟨h1, h3⟩, h2⟩
    · rintro ⟨⟨h1, h3⟩, h2⟩
      exact ⟨h1, h2, h3⟩
  have hmem : ∀ q, q ∈ L ↔ q ∈ (hamPathsR n T).filter (closes T) := by
    intro q
    show q ∈ C0.flatMap _ ↔ _
    rw [List.mem_flatMap, List.mem_filter, mem_hamPathsR]
    constructor
    · rintro ⟨c, hc, hq⟩
      rw [List.mem_map] at hq
      obtain ⟨k, _, rfl⟩ := hq
      exact rotN_closing n hn T k c ((hC0 c).1 hc).1
    · rintro ⟨hq1, hq2⟩
      obtain ⟨A, B, hAB⟩ := List.append_of_mem (zero_mem_of_hamPath n (by omega) T q hq1)
      have hcq : IsClosingHP n T (rotN A.length q) := rotN_closing n hn T _ q ⟨hq1, hq2⟩
      have hrot : rotN A.length q = 0 :: (B ++ A) := by
        rw [hAB, rotN_append]
        rfl
      refine ⟨rotN A.length q, (hC0 _).2 ⟨hcq, B ++ A, hrot⟩, ?_⟩
      rw [List.mem_map]
      have hlen : A.length + B.length + 1 = n := by
        rw [← hq1.1, hAB, List.length_append, List.length_cons]
        omega
      by_cases hA : A = []
      · subst hA
        refine ⟨0, (mem_downR n 0).2 (by omega), ?_⟩
        rw [hrot, hAB, List.append_nil, List.nil_append]
        rfl
      · have hApos : 0 < A.length := by
          cases A with
          | nil => exact absurd rfl hA
          | cons _ _ => exact Nat.succ_pos _
        refine ⟨B.length + 1, (mem_downR n _).2 (by omega), ?_⟩
        rw [hrot]
        have := rotN_append (0 :: B) A
        rw [List.length_cons] at this
        rw [show 0 :: (B ++ A) = (0 :: B) ++ A from rfl, this, hAB]
  have hnd : L.Nodup := by
    apply nodup_flatMap _ _ (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T))
    · intro c hc
      obtain ⟨⟨hham, _⟩, B, rfl⟩ := (hC0 c).1 hc
      have hB : B.length + 1 = n := by rw [← hham.1, List.length_cons]
      apply nodup_map_of_inj_on _ _ (nodup_downR n)
      intro k hk k' hk' hkk
      exact rot_index_unique n B B hB hB hham.2.1 k k' ((mem_downR n k).1 hk)
        ((mem_downR n k').1 hk') hkk
    · intro c c' hc hc' hne q hq hq'
      obtain ⟨⟨hham, _⟩, B, rfl⟩ := (hC0 c).1 hc
      obtain ⟨⟨hham', _⟩, B', rfl⟩ := (hC0 c').1 hc'
      rw [List.mem_map] at hq hq'
      obtain ⟨k, hk, rfl⟩ := hq
      obtain ⟨k', hk', hkk⟩ := hq'
      have hB : B.length + 1 = n := by rw [← hham.1, List.length_cons]
      have hB' : B'.length + 1 = n := by rw [← hham'.1, List.length_cons]
      have heq := rot_index_unique n B B' hB hB' hham.2.1 k k' ((mem_downR n k).1 hk)
        ((mem_downR n k').1 hk') hkk.symm
      subst heq
      apply hne
      apply rotN_inj k _ _ (by rw [List.length_cons, List.length_cons]; omega)
        (by rw [List.length_cons]; have := (mem_downR n k).1 hk; omega) hkk.symm
  have hlenL : L.length = n * C0.length := by
    show (C0.flatMap _).length = _
    rw [List.length_flatMap]
    rw [List.map_congr_left (g := fun _ => n) (fun c _ => by rw [List.length_map, length_downR])]
    exact sum_map_const n C0
  have := length_eq_of_same_members L _ hnd
    (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T)) hmem
  unfold hcCount
  rw [← this, hlenL]

/-- **`#copies(H_n) = H(T) - n · hc(T)`** for every tournament on `n ≥ 2` vertices. -/
theorem copiesH_add (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    copiesH n T + n * hcCount n T = hpCount n T := by
  have h1 := split_closes n hn T hT
  have h2 := closing_count n hn T
  omega

/-- Example: `C3[1, C3, 1]` has `H = 15` Hamiltonian paths, `hc = 3` Hamiltonian cycles and
no copy of `H_5`: all 15 paths close (`15 = 0 + 5 · 3`). -/
theorem c3c3_counts : hpCount 5 c3c3 = 15 ∧ hcCount 5 c3c3 = 3 ∧ copiesH 5 c3c3 = 0 := by
  decide +kernel

end ProcgenSelfieEdim
