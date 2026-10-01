set_option autoImplicit false

/-!
# Small list lemmas

Several core lemmas (`List.nodup_range`, `List.length_erase_of_mem`,
`List.mem_erase_of_ne`) depend on `Classical.choice` in Lean 4.30.0. The
versions below are proved directly so that the package stays inside
`{propext, Quot.sound}`.

The counting convention of the package: "the number of objects with property
`P`" is the length of a duplicate-free list whose members are exactly the
objects with property `P`. `length_eq_of_same_members` shows that this number
does not depend on the list chosen.
-/

namespace ProcgenSelfieEdim

theorem nodup_range (n : Nat) : (List.range n).Nodup := by
  induction n with
  | zero => exact List.Pairwise.nil
  | succ n ih =>
    rw [List.range_succ, List.nodup_append]
    refine ⟨ih, List.Pairwise.cons (fun _ h => absurd h (List.not_mem_nil)) List.Pairwise.nil, ?_⟩
    intro a ha b hb
    rw [List.mem_range] at ha
    rw [List.mem_singleton] at hb
    omega

theorem length_filter_lt {α : Type} (p : α → Bool) (a : α) :
    ∀ (l : List α), a ∈ l → p a = false → (l.filter p).length < l.length
  | [], ha, _ => absurd ha List.not_mem_nil
  | x :: t, ha, hp => by
    rcases List.mem_cons.1 ha with h | h
    · subst h
      rw [List.filter_cons_of_neg (by simp [hp])]
      exact Nat.lt_succ_of_le (List.length_filter_le p t)
    · cases hx : p x
      · rw [List.filter_cons_of_neg (by simp [hx])]
        exact Nat.lt_succ_of_le (List.length_filter_le p t)
      · rw [List.filter_cons_of_pos hx]
        exact Nat.succ_lt_succ (length_filter_lt p a t h hp)

/-- Pigeonhole: a duplicate-free list contained in `l2` is no longer than `l2`. -/
theorem pigeonhole {α : Type} [DecidableEq α] :
    ∀ (l1 l2 : List α), l1.Nodup → (∀ x ∈ l1, x ∈ l2) → l1.length ≤ l2.length
  | [], _, _, _ => Nat.zero_le _
  | a :: t, l2, hnd, hsub => by
    have ha : a ∈ l2 := hsub a List.mem_cons_self
    rw [List.nodup_cons] at hnd
    have ht : ∀ x ∈ t, x ∈ l2.filter (fun y => decide (y ≠ a)) := by
      intro x hx
      rw [List.mem_filter]
      refine ⟨hsub x (List.mem_cons_of_mem a hx), ?_⟩
      apply decide_eq_true
      intro hxa
      subst hxa
      exact hnd.1 hx
    have h1 := pigeonhole t (l2.filter (fun y => decide (y ≠ a))) hnd.2 ht
    have h2 := length_filter_lt (fun y => decide (y ≠ a)) a l2 ha (by simp)
    rw [List.length_cons]
    omega

/-- Two duplicate-free lists with the same members have the same length. -/
theorem length_eq_of_same_members {α : Type} [DecidableEq α] (l1 l2 : List α)
    (h1 : l1.Nodup) (h2 : l2.Nodup) (h : ∀ x, x ∈ l1 ↔ x ∈ l2) :
    l1.length = l2.length :=
  Nat.le_antisymm (pigeonhole l1 l2 h1 (fun x hx => (h x).1 hx))
    (pigeonhole l2 l1 h2 (fun x hx => (h x).2 hx))

/-- A function that is injective on a duplicate-free image list. -/
theorem eq_of_map_nodup {α β : Type} (f : α → β) :
    ∀ (l : List α), (l.map f).Nodup → ∀ a b, a ∈ l → b ∈ l → f a = f b → a = b
  | [], _, _, _, ha, _, _ => absurd ha List.not_mem_nil
  | x :: t, hnd, a, b, ha, hb, hab => by
    rw [List.map_cons, List.nodup_cons] at hnd
    rcases List.mem_cons.1 ha with ha' | ha' <;> rcases List.mem_cons.1 hb with hb' | hb'
    · rw [ha', hb']
    · subst ha'
      exact absurd (hab ▸ List.mem_map_of_mem hb') hnd.1
    · subst hb'
      exact absurd (hab ▸ List.mem_map_of_mem ha') hnd.1
    · exact eq_of_map_nodup f t hnd.2 a b ha' hb' hab

/-- A list whose image is duplicate-free is duplicate-free. -/
theorem nodup_of_map_nodup {α β : Type} (f : α → β) :
    ∀ (l : List α), (l.map f).Nodup → l.Nodup
  | [], _ => List.Pairwise.nil
  | x :: t, hnd => by
    rw [List.map_cons, List.nodup_cons] at hnd
    rw [List.nodup_cons]
    exact ⟨fun hx => hnd.1 (List.mem_map_of_mem hx), nodup_of_map_nodup f t hnd.2⟩

/-- Duplicate-freeness of a `flatMap` whose blocks are duplicate-free and disjoint. -/
theorem nodup_flatMap {α β : Type} (f : α → List β) :
    ∀ (l : List α), l.Nodup → (∀ a ∈ l, (f a).Nodup) →
      (∀ a b, a ∈ l → b ∈ l → a ≠ b → ∀ x, x ∈ f a → x ∈ f b → False) →
      (l.flatMap f).Nodup
  | [], _, _, _ => List.Pairwise.nil
  | a :: t, hnd, hblock, hdisj => by
    rw [List.flatMap_cons, List.nodup_append]
    rw [List.nodup_cons] at hnd
    refine ⟨hblock a List.mem_cons_self,
      nodup_flatMap f t hnd.2 (fun b hb => hblock b (List.mem_cons_of_mem a hb))
        (fun b c hb hc hbc => hdisj b c (List.mem_cons_of_mem a hb) (List.mem_cons_of_mem a hc) hbc),
      ?_⟩
    intro x hx y hy hxy
    subst hxy
    rw [List.mem_flatMap] at hy
    obtain ⟨b, hb, hxb⟩ := hy
    have hab : a ≠ b := by
      intro h
      subst h
      exact hnd.1 hb
    exact hdisj a b List.mem_cons_self (List.mem_cons_of_mem a hb) hab x hx hxb

/-- Boolean duplicate test, cheap to evaluate in the kernel. -/
def allDistinct {α : Type} [BEq α] : List α → Bool
  | [] => true
  | x :: xs => !(xs.contains x) && allDistinct xs

theorem nodup_of_allDistinct {α : Type} [BEq α] [LawfulBEq α] :
    ∀ (l : List α), allDistinct l = true → l.Nodup
  | [], _ => List.Pairwise.nil
  | x :: xs, h => by
    simp only [allDistinct, Bool.and_eq_true, Bool.not_eq_true'] at h
    rw [List.nodup_cons]
    refine ⟨?_, nodup_of_allDistinct xs h.2⟩
    intro hx
    have : xs.contains x = true := by simpa using hx
    rw [this] at h
    exact Bool.noConfusion h.1

/-- Every element of a list satisfies a Boolean test. -/
theorem all_of_all_eq_true {α : Type} (p : α → Bool) (l : List α) (h : l.all p = true) :
    ∀ x ∈ l, p x = true := by
  intro x hx
  exact (List.all_eq_true.1 h) x hx

/-- Exhaustive check over `0, …, n-1`. -/
theorem forall_lt_of_all_range (p : Nat → Bool) (n : Nat) (h : (List.range n).all p = true) :
    ∀ x, x < n → p x = true := by
  intro x hx
  exact all_of_all_eq_true p (List.range n) h x (List.mem_range.2 hx)

/-! ## Filter-length lemmas -/

theorem filter_cons_true {α : Type} (p : α → Bool) (x : α) (l : List α) (h : p x = true) :
    (x :: l).filter p = x :: l.filter p := List.filter_cons_of_pos h

theorem filter_cons_false {α : Type} (p : α → Bool) (x : α) (l : List α) (h : p x = false) :
    (x :: l).filter p = l.filter p :=
  List.filter_cons_of_neg (fun h' => Bool.noConfusion (h.symm.trans h'))

theorem exists_of_length_filter_ne_zero {α : Type} (p : α → Bool) :
    ∀ (l : List α), (l.filter p).length ≠ 0 → ∃ x, x ∈ l ∧ p x = true
  | [], h => absurd rfl h
  | x :: t, h => by
    cases hx : p x
    · rw [filter_cons_false p x t hx] at h
      obtain ⟨y, hy, hpy⟩ := exists_of_length_filter_ne_zero p t h
      exact ⟨y, List.mem_cons_of_mem x hy, hpy⟩
    · exact ⟨x, List.mem_cons_self, hx⟩

theorem length_filter_or_le {α : Type} (p q : α → Bool) :
    ∀ (l : List α), (l.filter (fun x => p x || q x)).length ≤ (l.filter p).length + (l.filter q).length
  | [] => Nat.le_refl 0
  | x :: t => by
    have ih := length_filter_or_le p q t
    cases hp : p x <;> cases hq : q x
    · rw [filter_cons_false _ x t (by rw [hp, hq]; rfl), filter_cons_false p x t hp,
        filter_cons_false q x t hq]
      exact ih
    · rw [filter_cons_true _ x t (by rw [hp, hq]; rfl), filter_cons_false p x t hp,
        filter_cons_true q x t hq]
      simp only [List.length_cons]
      omega
    · rw [filter_cons_true _ x t (by rw [hp, hq]; rfl), filter_cons_true p x t hp,
        filter_cons_false q x t hq]
      simp only [List.length_cons]
      omega
    · rw [filter_cons_true _ x t (by rw [hp, hq]; rfl), filter_cons_true p x t hp,
        filter_cons_true q x t hq]
      simp only [List.length_cons]
      omega

theorem length_filter_mono {α : Type} (p q : α → Bool) :
    ∀ (l : List α), (∀ x ∈ l, p x = true → q x = true) → (l.filter p).length ≤ (l.filter q).length
  | [], _ => Nat.le_refl 0
  | x :: t, h => by
    have ih := length_filter_mono p q t (fun y hy => h y (List.mem_cons_of_mem x hy))
    cases hp : p x
    · rw [filter_cons_false p x t hp]
      cases hq : q x
      · rw [filter_cons_false q x t hq]
        exact ih
      · rw [filter_cons_true q x t hq, List.length_cons]
        omega
    · rw [filter_cons_true p x t hp, filter_cons_true q x t (h x List.mem_cons_self hp),
        List.length_cons, List.length_cons]
      omega

theorem length_filter_add_not {α : Type} (p : α → Bool) :
    ∀ (l : List α), (l.filter p).length + (l.filter (fun x => !p x)).length = l.length
  | [] => rfl
  | x :: t => by
    have ih := length_filter_add_not p t
    cases hp : p x
    · rw [filter_cons_false p x t hp, filter_cons_true (fun x => !p x) x t (by show (!p x) = true; rw [hp]; rfl),
        List.length_cons, List.length_cons]
      omega
    · rw [filter_cons_true p x t hp, filter_cons_false (fun x => !p x) x t (by show (!p x) = false; rw [hp]; rfl),
        List.length_cons, List.length_cons]
      omega

theorem length_filter_false {α : Type} : ∀ (l : List α), (l.filter (fun _ => false)).length = 0
  | [] => rfl
  | x :: t => by rw [filter_cons_false _ x t rfl, length_filter_false t]

/-- Union bound: an element hit by some member of `S` is counted by that member. -/
theorem length_filter_any_le {α β : Type} (q : β → α → Bool) (E : List α) :
    ∀ (S : List β), (E.filter (fun e => S.any (fun s => q s e))).length ≤
      (S.map (fun s => (E.filter (q s)).length)).sum
  | [] => by
    show (E.filter (fun _ => false)).length ≤ 0
    rw [length_filter_false]
    exact Nat.le_refl 0
  | s :: t => by
    have h1 := length_filter_or_le (q s) (fun e => t.any (fun s' => q s' e)) E
    have h2 := length_filter_any_le q E t
    show (E.filter (fun e => q s e || t.any (fun s' => q s' e))).length ≤ _
    rw [List.map_cons, List.sum_cons]
    omega

theorem sum_le_mul {β : Type} (f : β → Nat) (c : Nat) :
    ∀ (S : List β), (∀ s ∈ S, f s ≤ c) → (S.map f).sum ≤ c * S.length
  | [], _ => Nat.zero_le _
  | s :: t, h => by
    have h1 := h s List.mem_cons_self
    have h2 := sum_le_mul f c t (fun s' hs' => h s' (List.mem_cons_of_mem s hs'))
    rw [List.map_cons, List.sum_cons, List.length_cons, Nat.mul_succ]
    omega

/-- A duplicate-free list mapped by a function injective on it is duplicate-free. -/
theorem nodup_map_of_inj_on {α β : Type} (f : α → β) :
    ∀ (l : List α), l.Nodup → (∀ a ∈ l, ∀ b ∈ l, f a = f b → a = b) → (l.map f).Nodup
  | [], _, _ => List.Pairwise.nil
  | x :: t, hnd, hinj => by
    rw [List.nodup_cons] at hnd
    rw [List.map_cons, List.nodup_cons]
    refine ⟨fun hx => ?_, nodup_map_of_inj_on f t hnd.2
      (fun a ha b hb hab => hinj a (List.mem_cons_of_mem x ha) b (List.mem_cons_of_mem x hb) hab)⟩
    rw [List.mem_map] at hx
    obtain ⟨y, hy, hxy⟩ := hx
    have := hinj y (List.mem_cons_of_mem x hy) x List.mem_cons_self hxy
    subst this
    exact hnd.1 hy


/-- Membership test through `Nat.beq`. -/
noncomputable def memR (x : Nat) (l : List Nat) : Bool :=
  List.rec (motive := fun _ => Bool) false (fun y _ ih => Nat.beq x y || ih) l

/-- Duplicate test through `memR`. -/
noncomputable def distinctR (l : List Nat) : Bool :=
  List.rec (motive := fun _ => Bool) true (fun x xs ih => !(memR x xs) && ih) l

theorem memR_iff (x : Nat) : ∀ (l : List Nat), memR x l = true ↔ x ∈ l
  | [] => by
    show false = true ↔ x ∈ []
    simp
  | y :: ys => by
    show (Nat.beq x y || memR x ys) = true ↔ x ∈ y :: ys
    rw [Bool.or_eq_true, memR_iff x ys, List.mem_cons, Nat.beq_eq]

theorem nodup_of_distinctR : ∀ (l : List Nat), distinctR l = true → l.Nodup
  | [], _ => List.Pairwise.nil
  | x :: xs, h => by
    have h' : (!(memR x xs) && distinctR xs) = true := h
    rw [Bool.and_eq_true, Bool.not_eq_true'] at h'
    rw [List.nodup_cons]
    refine ⟨fun hx => ?_, nodup_of_distinctR xs h'.2⟩
    rw [(memR_iff x xs).2 hx] at h'
    exact Bool.noConfusion h'.1


/-! ## Kernel-friendly list functions

Definitions by structural recursion compile to `brecOn`, which allocates a tuple at
each step of kernel reduction. The functions below are written with the recursors
`Nat.rec` / `List.rec` and are proved equal to their core counterparts; they are
used only inside finite checks evaluated by `decide +kernel`. -/

/-- Append, by `List.rec`. -/
noncomputable def appR {α : Type} (l1 l2 : List α) : List α :=
  List.rec (motive := fun _ => List α) l2 (fun x _ ih => x :: ih) l1

theorem appR_eq {α : Type} : ∀ (l1 l2 : List α), appR l1 l2 = l1 ++ l2
  | [], _ => rfl
  | x :: t, l2 => by
    show x :: appR t l2 = x :: (t ++ l2)
    rw [appR_eq t l2]

/-- `flatMap`, by `List.rec`. -/
noncomputable def flatR {α β : Type} (f : α → List β) (l : List α) : List β :=
  List.rec (motive := fun _ => List β) [] (fun x _ ih => appR (f x) ih) l

theorem flatR_eq {α β : Type} (f : α → List β) : ∀ (l : List α), flatR f l = l.flatMap f
  | [] => rfl
  | x :: t => by
    show appR (f x) (flatR f t) = (x :: t).flatMap f
    rw [appR_eq, flatR_eq f t, List.flatMap_cons]

/-- `[n-1, …, 1, 0]`, by `Nat.rec`. -/
noncomputable def downR (n : Nat) : List Nat :=
  Nat.rec (motive := fun _ => List Nat) [] (fun k ih => k :: ih) n

theorem downR_succ (n : Nat) : downR (n + 1) = n :: downR n := rfl

theorem mem_downR : ∀ (n x : Nat), x ∈ downR n ↔ x < n
  | 0, x => ⟨fun h => absurd h List.not_mem_nil, fun h => absurd h (Nat.not_lt_zero x)⟩
  | n + 1, x => by
    rw [downR_succ, List.mem_cons, mem_downR n x]
    constructor
    · rintro (h | h)
      · omega
      · omega
    · intro h
      by_cases hx : x = n
      · exact Or.inl hx
      · exact Or.inr (by omega)

theorem length_downR : ∀ (n : Nat), (downR n).length = n
  | 0 => rfl
  | n + 1 => by rw [downR_succ, List.length_cons, length_downR n]

theorem nodup_downR : ∀ (n : Nat), (downR n).Nodup
  | 0 => List.Pairwise.nil
  | n + 1 => by
    rw [downR_succ, List.nodup_cons]
    exact ⟨fun h => by rw [mem_downR] at h; omega, nodup_downR n⟩

/-- Number of list members satisfying `p`, by `List.rec`. -/
noncomputable def countR {α : Type} (p : α → Bool) (l : List α) : Nat :=
  List.rec (motive := fun _ => Nat) 0 (fun x _ ih => cond (p x) (ih + 1) ih) l

theorem countR_eq {α : Type} (p : α → Bool) : ∀ (l : List α), countR p l = (l.filter p).length
  | [] => rfl
  | x :: t => by
    show cond (p x) (countR p t + 1) (countR p t) = ((x :: t).filter p).length
    cases hp : p x
    · rw [filter_cons_false p x t hp, cond_false, countR_eq p t]
    · rw [filter_cons_true p x t hp, cond_true, countR_eq p t, List.length_cons]

/-- `l.all p`, by `List.rec`. -/
noncomputable def allR {α : Type} (p : α → Bool) (l : List α) : Bool :=
  List.rec (motive := fun _ => Bool) true (fun x _ ih => p x && ih) l

theorem allR_eq_true {α : Type} (p : α → Bool) : ∀ (l : List α), allR p l = true → ∀ x ∈ l, p x = true
  | [], _, x, hx => absurd hx List.not_mem_nil
  | y :: t, h, x, hx => by
    have h' : (p y && allR p t) = true := h
    rw [Bool.and_eq_true] at h'
    rcases List.mem_cons.1 hx with rfl | hx'
    · exact h'.1
    · exact allR_eq_true p t h'.2 x hx'

/-- `l.any p`, by `List.rec`. -/
noncomputable def anyR {α : Type} (p : α → Bool) (l : List α) : Bool :=
  List.rec (motive := fun _ => Bool) false (fun x _ ih => p x || ih) l

theorem anyR_eq_true {α : Type} (p : α → Bool) : ∀ (l : List α), anyR p l = true → ∃ x ∈ l, p x = true
  | [], h => Bool.noConfusion h
  | y :: t, h => by
    have h' : (p y || anyR p t) = true := h
    rw [Bool.or_eq_true] at h'
    rcases h' with h1 | h1
    · exact ⟨y, List.mem_cons_self, h1⟩
    · obtain ⟨x, hx, hpx⟩ := anyR_eq_true p t h1
      exact ⟨x, List.mem_cons_of_mem y hx, hpx⟩

/-! ## Finite sums -/

/-- `rsum n f = f 0 + f 1 + ⋯ + f (n-1)`. -/
def rsum : Nat → (Nat → Nat) → Nat
  | 0, _ => 0
  | n + 1, f => rsum n f + f n

theorem rsum_congr : ∀ (n : Nat) (f g : Nat → Nat), (∀ r, r < n → f r = g r) → rsum n f = rsum n g
  | 0, _, _, _ => rfl
  | n + 1, f, g, h => by
    show rsum n f + f n = rsum n g + g n
    rw [rsum_congr n f g (fun r hr => h r (by omega)), h n (by omega)]

theorem rsum_add : ∀ (n : Nat) (f g : Nat → Nat), rsum n (fun r => f r + g r) = rsum n f + rsum n g
  | 0, _, _ => rfl
  | n + 1, f, g => by
    show rsum n (fun r => f r + g r) + (f n + g n) = (rsum n f + f n) + (rsum n g + g n)
    rw [rsum_add n f g]
    omega

theorem rsum_zero : ∀ (n : Nat), rsum n (fun _ => 0) = 0
  | 0 => rfl
  | n + 1 => by
    show rsum n (fun _ => 0) + 0 = 0
    rw [rsum_zero n]

/-- Only the term `r = x` survives. -/
theorem rsum_single : ∀ (n x : Nat) (g : Nat → Nat), x < n →
    rsum n (fun r => if x = r then g r else 0) = g x
  | 0, _, _, hx => absurd hx (Nat.not_lt_zero _)
  | n + 1, x, g, hx => by
    show rsum n (fun r => if x = r then g r else 0) + (if x = n then g n else 0) = g x
    by_cases hxn : x = n
    · subst hxn
      rw [if_pos rfl, rsum_congr x _ (fun _ => 0) (fun r hr => if_neg (by omega)), rsum_zero]
      omega
    · rw [if_neg hxn, rsum_single n x g (by omega)]
      omega

end ProcgenSelfieEdim
