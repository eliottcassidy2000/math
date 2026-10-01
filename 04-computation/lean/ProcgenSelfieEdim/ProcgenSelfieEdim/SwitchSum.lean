import ProcgenSelfieEdim.HPExist

set_option autoImplicit false

/-!
# THM-4524 A2: the loop-set sum `Σ_L H(switch_L T) = 2 · N!`

* `hpCount_complete`: there are `N!` orderings of `{0, …, N-1}` (Hamiltonian paths of the
  complete relation), counted with the verified enumerator.
* `two_loop_sets`: for a tournament `T` and an ordering `π`, the loop sets `L` making `π` a
  directed path of `switch L T` are exactly `L₀` and its complement (A2, first claim).
* `switch_sum`: `Σ_{L ⊆ {0,…,N-1}} H(switch_L T) = 2 · N!` (A2), the loop sets running over
  the `2^N` bit lists `allBools N`.
-/

namespace ProcgenSelfieEdim

/-! ## Counting orderings -/

/-- `n!`. -/
def fact : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fact n

/-- The falling factorial `a (a-1) ⋯ (a-k+1)`. -/
def ffact : Nat → Nat → Nat
  | _, 0 => 1
  | a, k + 1 => a * ffact (a - 1) k

theorem ffact_self : ∀ (a : Nat), ffact a a = fact a
  | 0 => rfl
  | a + 1 => by
    show (a + 1) * ffact (a + 1 - 1) a = (a + 1) * fact a
    rw [Nat.add_sub_cancel, ffact_self a]

theorem sum_map_ite {α : Type} (q : α → Bool) (c : Nat) :
    ∀ (l : List α), (l.map (fun x => if q x = true then c else 0)).sum = c * (l.filter q).length
  | [] => by rw [List.map_nil, List.sum_nil, List.filter_nil, List.length_nil, Nat.mul_zero]
  | x :: t => by
    rw [List.map_cons, List.sum_cons, sum_map_ite q c t]
    cases hx : q x
    · rw [filter_cons_false q x t hx, if_neg (by decide)]
      omega
    · rw [filter_cons_true q x t hx, if_pos rfl, List.length_cons, Nat.mul_succ]
      omega

/-- The vertices below `n` not on a duplicate-free list `p ⊆ {0, …, n-1}` number `n - |p|`. -/
theorem count_not_mem (n : Nat) (p : List Nat) (hnd : p.Nodup) (hlt : ∀ x ∈ p, x < n) :
    ((downR n).filter (fun w => !memR w p)).length + p.length = n := by
  have h1 := length_filter_add_not (fun w => memR w p) (downR n)
  have h2 : ((downR n).filter (fun w => memR w p)).length = p.length := by
    apply length_eq_of_same_members _ _ (List.Nodup.sublist List.filter_sublist (nodup_downR n)) hnd
    intro x
    rw [List.mem_filter, memR_iff, mem_downR]
    exact ⟨fun h => h.2, fun h => ⟨hlt x h, h⟩⟩
  rw [length_downR] at h1
  omega

/-- The depth-first extension in the complete relation has `(n - |p|)(n - |p| - 1)⋯` branches. -/
theorem length_extR_complete (n : Nat) :
    ∀ (k : Nat) (p : List Nat), p ≠ [] → p.Nodup → (∀ x ∈ p, x < n) →
      (extR n complete k p).length = ffact (n - p.length) k
  | 0, _, _, _, _ => rfl
  | _ + 1, [], hne, _, _ => absurd rfl hne
  | k + 1, h :: rest, _, hnd, hlt => by
    rw [extR_succ_cons, flatR_eq, List.length_flatMap]
    have hblock : ∀ w ∈ downR n,
        (cond (complete w h && !memR w (h :: rest)) (extR n complete k (w :: h :: rest)) []).length =
          if (!memR w (h :: rest)) = true then ffact (n - (h :: rest).length - 1) k else 0 := by
      intro w hw
      have hwn := (mem_downR n w).1 hw
      cases hm : memR w (h :: rest)
      · rw [if_pos (show (!false) = true from rfl)]
        show (cond (true && true) _ []).length = _
        rw [Bool.and_true, cond_true]
        have hwp : w ∉ h :: rest := fun hm' => by rw [(memR_iff w _).2 hm'] at hm; exact Bool.noConfusion hm
        rw [length_extR_complete n k (w :: h :: rest) (List.cons_ne_nil _ _)
          (List.nodup_cons.2 ⟨hwp, hnd⟩) (fun x hx => by
            rcases List.mem_cons.1 hx with rfl | hx'
            · exact hwn
            · exact hlt x hx')]
        rw [List.length_cons, List.length_cons]
        congr 1
      · rw [if_neg (show ¬ (!true) = true by decide)]
        show (cond (true && false) _ []).length = 0
        rfl
    rw [List.map_congr_left hblock, sum_map_ite, ffact]
    have hc : ((downR n).filter (fun w => !memR w (h :: rest))).length = n - (h :: rest).length := by
      have := count_not_mem n (h :: rest) hnd hlt
      omega
    rw [hc]
    exact Nat.mul_comm _ _

/-- **There are `N!` orderings of `{0, …, N-1}`.** -/
theorem hpCount_complete (N : Nat) : hpCount N complete = fact N := by
  cases N with
  | zero => rfl
  | succ m =>
    unfold hpCount
    rw [hamPathsR_succ, flatR_eq, List.length_flatMap]
    rw [List.map_congr_left (g := fun _ => fact m) (fun v hv => by
      rw [length_extR_complete (m + 1) m [v] (List.cons_ne_nil _ _)
        (List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩)
        (fun x hx => by rw [List.mem_singleton] at hx; rw [hx]; exact (mem_downR _ v).1 hv),
        List.length_singleton, Nat.add_sub_cancel, ffact_self])]
    rw [sum_map_const, length_downR]
    show fact m * (m + 1) = (m + 1) * fact m
    exact Nat.mul_comm _ _

/-! ## Loop sets as bit lists -/

/-- All bit lists of length `k`; `nthB bs` is the loop set `{i < k : bit i set}`. -/
def allBools : Nat → List (List Bool)
  | 0 => [[]]
  | k + 1 => (allBools k).flatMap (fun l => [l ++ [false], l ++ [true]])

theorem mem_allBools : ∀ (k : Nat) (l : List Bool), l ∈ allBools k ↔ l.length = k
  | 0, l => by
    show l ∈ [[]] ↔ _
    rw [List.mem_singleton]
    constructor
    · rintro rfl
      rfl
    · intro h
      cases l with
      | nil => rfl
      | cons _ _ =>
        rw [List.length_cons] at h
        omega
  | k + 1, l => by
    show l ∈ (allBools k).flatMap _ ↔ _
    rw [List.mem_flatMap]
    constructor
    · rintro ⟨l', hl', hl⟩
      have := (mem_allBools k l').1 hl'
      rcases List.mem_cons.1 hl with rfl | hl2
      · rw [List.length_append, this]
        rfl
      · rw [List.mem_singleton] at hl2
        subst hl2
        rw [List.length_append, this]
        rfl
    · intro h
      rcases List.eq_nil_or_concat l with h0 | ⟨l', b, h0⟩
      · subst h0
        rw [List.length_nil] at h
        omega
      · rw [List.concat_eq_append] at h0
        subst h0
        rw [List.length_append, List.length_singleton] at h
        refine ⟨l', (mem_allBools k l').2 (by omega), ?_⟩
        cases b
        · exact List.mem_cons_self
        · exact List.mem_cons_of_mem _ List.mem_cons_self

theorem nodup_allBools : ∀ (k : Nat), (allBools k).Nodup
  | 0 => List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩
  | k + 1 => by
    show ((allBools k).flatMap _).Nodup
    apply nodup_flatMap _ _ (nodup_allBools k)
    · intro l _
      refine List.nodup_cons.2 ⟨fun h => ?_, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩⟩
      rw [List.mem_singleton] at h
      have := (List.append_inj h rfl).2
      exact Bool.noConfusion (List.cons.inj this).1
    · intro l l' hl hl' hne x hx hx'
      have h1 := (mem_allBools k l).1 hl
      have h2 := (mem_allBools k l').1 hl'
      have e1 : ∃ b, x = l ++ [b] := by
        rcases List.mem_cons.1 hx with h | h
        · exact ⟨false, h⟩
        · exact ⟨true, List.mem_singleton.1 h⟩
      have e2 : ∃ b, x = l' ++ [b] := by
        rcases List.mem_cons.1 hx' with h | h
        · exact ⟨false, h⟩
        · exact ⟨true, List.mem_singleton.1 h⟩
      obtain ⟨b, rfl⟩ := e1
      obtain ⟨b', hb'⟩ := e2
      exact hne (List.append_inj hb' (h1.trans h2.symm)).1

theorem length_allBools : ∀ (k : Nat), (allBools k).length = 2 ^ k
  | 0 => rfl
  | k + 1 => by
    show ((allBools k).flatMap _).length = _
    rw [List.length_flatMap, List.map_congr_left (f := fun a => [a ++ [false], a ++ [true]].length)
      (g := fun _ => 2) (fun l _ => rfl), sum_map_const,
      length_allBools k, Nat.pow_succ, Nat.mul_comm]

theorem nthB_ext : ∀ (l1 l2 : List Bool), l1.length = l2.length →
    (∀ i, i < l1.length → nthB l1 i = nthB l2 i) → l1 = l2
  | [], [], _, _ => rfl
  | [], _ :: _, h, _ => by rw [List.length_nil, List.length_cons] at h; omega
  | _ :: _, [], h, _ => by rw [List.length_nil, List.length_cons] at h; omega
  | x :: l1, y :: l2, h, hi => by
    have h0 : x = y := hi 0 (by rw [List.length_cons]; omega)
    rw [h0, nthB_ext l1 l2 (by rw [List.length_cons, List.length_cons] at h; omega)
      (fun i hi' => hi (i + 1) (by rw [List.length_cons]; omega))]

/-! ## Paths of switched tournaments -/

/-- A Boolean path test. -/
def isPathB (T : Nat → Nat → Bool) : List Nat → Bool
  | a :: b :: rest => T a b && isPathB T (b :: rest)
  | _ => true

theorem isPathB_iff (T : Nat → Nat → Bool) : ∀ (p : List Nat), isPathB T p = true ↔ IsPath T p
  | [] => ⟨fun _ => trivial, fun _ => rfl⟩
  | [_] => ⟨fun _ => trivial, fun _ => rfl⟩
  | a :: b :: rest => by
    show (T a b && isPathB T (b :: rest)) = true ↔ T a b = true ∧ IsPath T (b :: rest)
    rw [Bool.and_eq_true, isPathB_iff T (b :: rest)]

theorem switch_pair (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) (L : Nat → Bool)
    (a b : Nat) (ha : a < N) (hb : b < N) (hab : a ≠ b) :
    switch L T a b = true ↔ (L a = L b ↔ T a b = true) := by
  unfold switch
  by_cases hl : L a = L b
  · rw [if_pos hl]
    exact ⟨fun h => ⟨fun _ => h, fun _ => hl⟩, fun h => h.1 hl⟩
  · rw [if_neg hl, hT a b ha hb hab]
    cases T a b
    · exact ⟨fun _ => ⟨fun h => absurd h hl, fun h => Bool.noConfusion h⟩, fun _ => rfl⟩
    · exact ⟨fun h => Bool.noConfusion h, fun h => absurd (h.2 rfl) hl⟩

theorem bool_flip : ∀ (p q p' q' t : Bool), (p = q ↔ t = true) → (p' = q' ↔ t = true) →
    (q = q' ↔ p = p') := by
  intro p q p' q' t
  cases p <;> cases q <;> cases p' <;> cases q' <;> cases t <;> decide

theorem isPath_switch_congr (T : Nat → Nat → Bool) (L L' : Nat → Bool) (p : List Nat)
    (h : ∀ x ∈ p, L x = L' x) : IsPath (switch L T) p ↔ IsPath (switch L' T) p := by
  rw [isPath_iff_pairs, isPath_iff_pairs]
  have key : ∀ l1 l2 a b, p = l1 ++ a :: b :: l2 → switch L T a b = switch L' T a b := by
    intro l1 l2 a b hp
    have ha : a ∈ p := by rw [hp]; exact List.mem_append_right _ List.mem_cons_self
    have hb : b ∈ p := by
      rw [hp]
      exact List.mem_append_right _ (List.mem_cons_of_mem _ List.mem_cons_self)
    unfold switch
    rw [h a ha, h b hb]
  constructor
  · intro h1 l1 l2 a b hp
    rw [← key l1 l2 a b hp]
    exact h1 l1 l2 a b hp
  · intro h1 l1 l2 a b hp
    rw [key l1 l2 a b hp]
    exact h1 l1 l2 a b hp

/-- Two loop sets that make the same ordering a path agree or are complementary along it. -/
theorem loops_along (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) (L L' : Nat → Bool) :
    ∀ (rest : List Nat) (a : Nat), (a :: rest).Nodup → (∀ x ∈ a :: rest, x < N) →
      IsPath (switch L T) (a :: rest) → IsPath (switch L' T) (a :: rest) →
      ∀ x ∈ a :: rest, (L x = L' x ↔ L a = L' a)
  | [], a, _, _, _, _, x, hx => by
    rw [List.mem_singleton] at hx
    rw [hx]
  | b :: rest, a, hnd, hlt, hp, hp', x, hx => by
    have hab : a ≠ b := fun h => (List.nodup_cons.1 hnd).1 (h ▸ List.mem_cons_self)
    have ha := hlt a List.mem_cons_self
    have hb := hlt b (List.mem_cons_of_mem a List.mem_cons_self)
    have e1 := (switch_pair N T hT L a b ha hb hab).1 hp.1
    have e2 := (switch_pair N T hT L' a b ha hb hab).1 hp'.1
    have hflip := bool_flip (L a) (L b) (L' a) (L' b) (T a b) e1 e2
    rcases List.mem_cons.1 hx with rfl | hx'
    · exact Iff.rfl
    · rw [loops_along N T hT L L' rest b (List.nodup_cons.1 hnd).2
        (fun y hy => hlt y (List.mem_cons_of_mem a hy)) hp.2 hp'.2 x hx', hflip]

/-- The loop set built along an ordering. -/
def loopsAlong (T : Nat → Nat → Bool) : List Nat → Nat → Bool
  | [] => fun _ => false
  | [_] => fun _ => false
  | a :: b :: rest => fun x =>
    if x = a then (if T a b = true then loopsAlong T (b :: rest) b else !loopsAlong T (b :: rest) b)
    else loopsAlong T (b :: rest) x

theorem loopsAlong_path (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) :
    ∀ (p : List Nat), p.Nodup → (∀ x ∈ p, x < N) → IsPath (switch (loopsAlong T p) T) p
  | [], _, _ => trivial
  | [_], _, _ => trivial
  | a :: b :: rest, hnd, hlt => by
    have hab : a ≠ b := fun h => (List.nodup_cons.1 hnd).1 (h ▸ List.mem_cons_self)
    have ha := hlt a List.mem_cons_self
    have hb := hlt b (List.mem_cons_of_mem a List.mem_cons_self)
    have hrest := loopsAlong_path N T hT (b :: rest) (List.nodup_cons.1 hnd).2
      (fun y hy => hlt y (List.mem_cons_of_mem a hy))
    have hagree : ∀ x ∈ b :: rest, loopsAlong T (a :: b :: rest) x = loopsAlong T (b :: rest) x := by
      intro x hx
      have hxa : x ≠ a := fun h => (List.nodup_cons.1 hnd).1 (h ▸ hx)
      show (if x = a then _ else _) = _
      rw [if_neg hxa]
    refine ⟨?_, (isPath_switch_congr T _ _ (b :: rest) hagree).2 hrest⟩
    apply (switch_pair N T hT _ a b ha hb hab).2
    rw [hagree b List.mem_cons_self]
    show (if a = a then _ else _) = _ ↔ _
    rw [if_pos rfl]
    cases T a b
    · rw [if_neg (by decide)]
      constructor
      · intro h
        cases hc : loopsAlong T (b :: rest) b <;> rw [hc] at h <;> exact absurd h (by decide)
      · intro h
        exact Bool.noConfusion h
    · rw [if_pos rfl]
      exact ⟨fun _ => rfl, fun _ => rfl⟩

theorem mem_of_hamPath (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) (hp : IsHamPath n T p)
    (x : Nat) (hx : x < n) : x ∈ p := by
  by_cases h : x ∈ p
  · exact h
  · exfalso
    have h1 := pigeonhole p ((downR n).filter (fun y => decide (y ≠ x))) hp.2.1 (fun y hy => by
      rw [List.mem_filter]
      exact ⟨(mem_downR n y).2 (hp.2.2.1 y hy), decide_eq_true (fun hyx => h (hyx ▸ hy))⟩)
    have h2 := length_filter_lt (fun y => decide (y ≠ x)) x (downR n) ((mem_downR n x).2 hx)
      (by simp)
    rw [length_downR] at h2
    rw [hp.1] at h1
    omega

/-- **A2, first claim.** For a tournament `T` and an ordering `π` of `{0, …, N-1}`, the loop
sets making `π` a directed path of `switch_L T` are exactly `L₀` and its complement. -/
theorem two_loop_sets (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) (π : List Nat)
    (hπ : IsHamPath N complete π) :
    ∃ L₀ : Nat → Bool, ∀ L : Nat → Bool,
      IsPath (switch L T) π ↔ (SameLoops N L L₀ ∨ SameLoops N L (compl L₀)) := by
  refine ⟨loopsAlong T π, fun L => ?_⟩
  have h0 := loopsAlong_path N T hT π hπ.2.1 hπ.2.2.1
  constructor
  · intro hL
    cases π with
    | nil =>
      refine Or.inl (fun x hx => ?_)
      have := hπ.1
      rw [List.length_nil] at this
      omega
    | cons a rest =>
      have hrel := loops_along N T hT L (loopsAlong T (a :: rest)) rest a hπ.2.1 hπ.2.2.1 hL h0
      by_cases ha : L a = loopsAlong T (a :: rest) a
      · exact Or.inl (fun x hx => (hrel x (mem_of_hamPath N complete _ hπ x hx)).2 ha)
      · refine Or.inr (fun x hx => ?_)
        have hne : ¬ L x = loopsAlong T (a :: rest) x :=
          fun h => ha ((hrel x (mem_of_hamPath N complete _ hπ x hx)).1 h)
        unfold compl
        cases h1 : L x <;> cases h2 : loopsAlong T (a :: rest) x <;> rw [h1, h2] at hne <;>
          first | rfl | exact absurd rfl hne
  · rintro (h | h)
    · exact (isPath_switch_congr T L _ π (fun x hx => h x (hπ.2.2.1 x hx))).2 h0
    · have h1 := (isPath_switch_congr T L _ π (fun x hx => h x (hπ.2.2.1 x hx))).2
      apply h1
      rw [isPath_iff_pairs] at h0 ⊢
      intro l1 l2 a b hp
      rw [switch_compl]
      exact h0 l1 l2 a b hp

/-- Exactly two of the `2^N` loop sets make a given ordering a path. -/
theorem count_loop_sets (N : Nat) (hN : 1 ≤ N) (T : Nat → Nat → Bool) (hT : IsTournament N T)
    (π : List Nat) (hπ : IsHamPath N complete π) :
    ((allBools N).filter (fun bs => isPathB (switch (nthB bs) T) π)).length = 2 := by
  obtain ⟨L₀, hL₀⟩ := two_loop_sets N T hT π hπ
  have hdist : ascList L₀ N ≠ ascList (compl L₀) N := by
    intro h
    have := congrArg (fun l => nthB l 0) h
    simp only at this
    rw [nthB_ascList _ N 0 (by omega), nthB_ascList _ N 0 (by omega)] at this
    unfold compl at this
    cases hc : L₀ 0 <;> rw [hc] at this <;> exact Bool.noConfusion this
  have hpair : [ascList L₀ N, ascList (compl L₀) N].Nodup :=
    List.nodup_cons.2 ⟨fun h => hdist (List.mem_singleton.1 h),
      List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩⟩
  have := length_eq_of_same_members _ [ascList L₀ N, ascList (compl L₀) N]
    (List.Nodup.sublist List.filter_sublist (nodup_allBools N)) hpair (fun bs => by
      rw [List.mem_filter, mem_allBools, isPathB_iff, hL₀, List.mem_cons, List.mem_singleton]
      constructor
      · rintro ⟨hl, h | h⟩
        · exact Or.inl (nthB_ext _ _ (by rw [ascList_length]; exact hl) (fun i hi =>
            by rw [nthB_ascList _ N i (by omega)]; exact h i (by omega)))
        · exact Or.inr (nthB_ext _ _ (by rw [ascList_length]; exact hl) (fun i hi =>
            by rw [nthB_ascList _ N i (by omega)]; exact h i (by omega)))
      · rintro (rfl | rfl)
        · exact ⟨ascList_length _ N, Or.inl (fun x hx => nthB_ascList _ N x hx)⟩
        · exact ⟨ascList_length _ N, Or.inr (fun x hx => nthB_ascList _ N x hx)⟩)
  rw [this]
  rfl

/-! ## The sum -/

theorem length_filter_cons {α : Type} (p : α → Bool) (x : α) (l : List α) :
    ((x :: l).filter p).length = ind (p x) + (l.filter p).length := by
  cases hx : p x
  · rw [filter_cons_false p x l hx]
    show _ = (if false = true then 1 else 0) + _
    rw [if_neg (by decide), Nat.zero_add]
  · rw [filter_cons_true p x l hx, List.length_cons]
    show _ = (if true = true then 1 else 0) + _
    rw [if_pos rfl]
    omega

theorem sum_map_add {α : Type} (f g : α → Nat) :
    ∀ (l : List α), (l.map (fun x => f x + g x)).sum = (l.map f).sum + (l.map g).sum
  | [] => rfl
  | x :: t => by
    rw [List.map_cons, List.map_cons, List.map_cons, List.sum_cons, List.sum_cons, List.sum_cons,
      sum_map_add f g t]
    omega

/-- Double counting of the pairs `(x, y)` with `R x y`. -/
theorem sum_filter_swap {α β : Type} (R : α → β → Bool) (Y : List β) :
    ∀ (X : List α), (X.map (fun x => (Y.filter (R x)).length)).sum =
      (Y.map (fun y => (X.filter (fun x => R x y)).length)).sum
  | [] => by
    rw [List.map_nil, List.sum_nil]
    rw [List.map_congr_left (f := fun y => (List.filter (fun x => R x y) []).length) (g := fun _ => 0)
      (fun y _ => rfl), sum_map_const, Nat.zero_mul]
  | x :: X => by
    rw [List.map_cons, List.sum_cons, sum_filter_swap R Y X]
    rw [List.map_congr_left (g := fun y => ind (R x y) + (X.filter (fun x' => R x' y)).length)
      (fun y _ => length_filter_cons (fun x' => R x' y) x X)]
    rw [sum_map_add, length_filter_eq_sum]

/-- **THM-4524 A2.** For every tournament `T` on `N ≥ 1` vertices,
`Σ_{L ⊆ {0,…,N-1}} H(switch_L T) = 2 · N!` (the sum runs over the `2^N` loop sets). -/
theorem switch_sum (N : Nat) (hN : 1 ≤ N) (T : Nat → Nat → Bool) (hT : IsTournament N T) :
    ((allBools N).map (fun bs => hpCount N (switch (nthB bs) T))).sum = 2 * fact N := by
  have hH : ∀ L : Nat → Bool, hpCount N (switch L T) =
      ((hamPathsR N complete).filter (isPathB (switch L T))).length := by
    intro L
    apply (hpCount_unique N (switch L T) _
      (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR N complete)) _).symm
    intro p
    rw [List.mem_filter, mem_hamPathsR, isPathB_iff]
    constructor
    · rintro ⟨⟨h1, h2, h3, _⟩, h5⟩
      exact ⟨h1, h2, h3, h5⟩
    · rintro ⟨h1, h2, h3, h5⟩
      exact ⟨⟨h1, h2, h3, isPath_complete p⟩, h5⟩
  rw [List.map_congr_left (fun bs _ => hH (nthB bs))]
  rw [sum_filter_swap (fun bs π => isPathB (switch (nthB bs) T) π) (hamPathsR N complete)]
  rw [List.map_congr_left (g := fun _ => 2) (fun π hπ =>
    count_loop_sets N hN T hT π ((mem_hamPathsR N complete π).1 hπ))]
  rw [sum_map_const]
  have := hpCount_complete N
  unfold hpCount at this
  rw [this]

end ProcgenSelfieEdim
