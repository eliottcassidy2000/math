import ProcgenSelfieEdim.SwitchSum

set_option autoImplicit false

/-!
# Rédei's theorem: every tournament has an odd number of Hamiltonian paths

Proof by induction on `n`, deleting the top vertex `x = n` of a tournament `T` on
`{0, …, n}`. Every Hamiltonian path of `T` has `x` first, last, or between two vertices
`a → x → b`; hence (`hp_decomp`, `count_head`, `count_last`, `count_interior`)

  `H(T) = #{Q : x → first Q} + #{Q : last Q → x} + Σ_{a,b} [a → x → b] · c_{T_ab}(a → b)`,

where `Q` runs over the Hamiltonian paths of `T' = T - x` and `T_ab` is `T'` with the pair
`{a, b}` oriented `a → b`. By the induction hypothesis and single-arc shaving,
`c_{T_ab}(a → b) ≡ c'(a → b) + c'(b → a) (mod 2)`, so each `Q` contributes
`[x → first] + [last → x] + #{consecutive pairs of Q separated by x}`, which is odd
(the side of `x` changes an odd number of times along `0, s_1, …, s_n, 1`). Therefore
`H(T) ≡ H(T') ≡ 1 (mod 2)`.
-/

namespace ProcgenSelfieEdim

/-! ## Path helpers -/

theorem isPath_append_cons (T : Nat → Nat → Bool) : ∀ (l1 : List Nat) (c : Nat) (l2 : List Nat),
    IsPath T (l1 ++ c :: l2) ↔ IsPath T (l1 ++ [c]) ∧ IsPath T (c :: l2)
  | [], _, _ => ⟨fun h => ⟨trivial, h⟩, fun h => h.2⟩
  | [_], _, _ => ⟨fun h => ⟨⟨h.1, trivial⟩, h.2⟩, fun h => ⟨h.1.1, h.2⟩⟩
  | y :: z :: t, c, l2 => by
    show (T y z = true ∧ IsPath T (z :: t ++ c :: l2)) ↔
      (T y z = true ∧ IsPath T (z :: t ++ [c])) ∧ IsPath T (c :: l2)
    rw [isPath_append_cons T (z :: t) c l2]
    constructor
    · rintro ⟨h1, h2, h3⟩
      exact ⟨⟨h1, h2⟩, h3⟩
    · rintro ⟨⟨h1, h2⟩, h3⟩
      exact ⟨h1, h2, h3⟩

theorem isPath_prefix (T : Nat → Nat → Bool) : ∀ (l1 l2 : List Nat), IsPath T (l1 ++ l2) → IsPath T l1
  | [], _, _ => trivial
  | [_], _, _ => trivial
  | _ :: z :: t, l2, h => ⟨h.1, isPath_prefix T (z :: t) l2 h.2⟩

theorem isPath_snoc_iff (T : Nat → Nat → Bool) : ∀ (c : Nat) (l : List Nat) (x : Nat),
    IsPath T (c :: l ++ [x]) ↔ IsPath T (c :: l) ∧ T (lastOr c l) x = true
  | _, [], _ => ⟨fun h => ⟨trivial, h.1⟩, fun h => ⟨h.2, trivial⟩⟩
  | c, d :: l, x => by
    show (T c d = true ∧ IsPath T (d :: l ++ [x])) ↔
      (T c d = true ∧ IsPath T (d :: l)) ∧ T (lastOr d l) x = true
    rw [isPath_snoc_iff T d l x]
    constructor
    · rintro ⟨h1, h2, h3⟩
      exact ⟨⟨h1, h2⟩, h3⟩
    · rintro ⟨⟨h1, h2⟩, h3⟩
      exact ⟨h1, h2, h3⟩

/-- Relations agreeing on the consecutive pairs of `l` give the same paths. -/
theorem isPath_congr_pairs (T T' : Nat → Nat → Bool) (l : List Nat)
    (h : ∀ l1 l2 a b, l = l1 ++ a :: b :: l2 → T a b = T' a b) : IsPath T l ↔ IsPath T' l := by
  rw [isPath_iff_pairs, isPath_iff_pairs]
  constructor
  · intro h1 l1 l2 a b hl
    rw [← h l1 l2 a b hl]
    exact h1 l1 l2 a b hl
  · intro h1 l1 l2 a b hl
    rw [h l1 l2 a b hl]
    exact h1 l1 l2 a b hl

/-- An HP of `T` on `{0, …, n}` minus `x = n`: membership bookkeeping. -/
theorem lt_of_mem_ne (n : Nat) (p : List Nat) (hlt : ∀ z ∈ p, z < n + 1) (y : Nat) (hy : y ∈ p)
    (hne : y ≠ n) : y < n := by
  have := hlt y hy
  omega

/-! ## The first and the last vertex -/

/-- The first vertex of `p` is `x`. -/
noncomputable def headIs (x : Nat) (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a _ => Nat.beq a x)

/-- The last vertex of `p` is `x`. -/
noncomputable def lastIs (x : Nat) (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a rest => Nat.beq (lastOr a rest) x)

/-- `x` beats the first vertex of `Q`. -/
noncomputable def beatsHead (T : Nat → Nat → Bool) (x : Nat) (Q : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) Q false (fun q _ => T x q)

/-- The last vertex of `Q` beats `x`. -/
noncomputable def lastBeats (T : Nat → Nat → Bool) (x : Nat) (Q : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) Q false (fun q rest => T (lastOr q rest) x)

/-- Hamiltonian paths of `T` on `{0, …, n}` starting at `n` ↔ Hamiltonian paths `Q` of `T` on
`{0, …, n-1}` whose first vertex `n` beats. -/
theorem count_head (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) :
    ((hamPathsR (n + 1) T).filter (headIs n)).length =
      ((hamPathsR n T).filter (beatsHead T n)).length := by
  have hL2 : ((hamPathsR n T).filter (beatsHead T n)).Nodup :=
    List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T)
  have hmap : (((hamPathsR n T).filter (beatsHead T n)).map (fun Q => n :: Q)).Nodup :=
    nodup_map_of_inj_on _ _ hL2 (fun a _ b _ h => (List.cons.inj h).2)
  rw [length_eq_of_same_members _ _ (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR _ T))
    hmap (fun p => ?_), List.length_map]
  rw [List.mem_filter, mem_hamPathsR, List.mem_map]
  constructor
  · rintro ⟨⟨hl, hnd, hlt, hpath⟩, hh⟩
    cases p with
    | nil => exact Bool.noConfusion hh
    | cons a Q =>
      have ha : a = n := Nat.eq_of_beq_eq_true hh
      subst ha
      rw [List.nodup_cons] at hnd
      have hQlt : ∀ y ∈ Q, y < a := fun y hy =>
        lt_of_mem_ne a (a :: Q) hlt y (List.mem_cons_of_mem a hy) (fun h => hnd.1 (h ▸ hy))
      cases Q with
      | nil =>
        exfalso
        rw [List.length_singleton] at hl
        omega
      | cons q rest =>
        refine ⟨q :: rest, ?_, rfl⟩
        rw [List.mem_filter, mem_hamPathsR]
        refine ⟨⟨by simp only [List.length_cons] at hl ⊢; omega, hnd.2, hQlt, hpath.2⟩, hpath.1⟩
  · rintro ⟨Q, hQ, rfl⟩
    rw [List.mem_filter, mem_hamPathsR] at hQ
    obtain ⟨⟨hl, hnd, hlt, hpath⟩, hh⟩ := hQ
    cases Q with
    | nil => exact Bool.noConfusion hh
    | cons q rest =>
      refine ⟨⟨by rw [List.length_cons, hl], ?_, ?_, ⟨hh, hpath⟩⟩, Nat.beq_refl n⟩
      · exact List.nodup_cons.2 ⟨fun h => Nat.lt_irrefl n (hlt n h), hnd⟩
      · intro y hy
        rcases List.mem_cons.1 hy with rfl | hy'
        · exact Nat.lt_succ_self _
        · exact Nat.lt_succ_of_lt (hlt y hy')

/-- Hamiltonian paths ending at `n` ↔ Hamiltonian paths `Q` on `{0, …, n-1}` whose last
vertex beats `n`. -/
theorem count_last (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) :
    ((hamPathsR (n + 1) T).filter (lastIs n)).length =
      ((hamPathsR n T).filter (lastBeats T n)).length := by
  have hL2 : ((hamPathsR n T).filter (lastBeats T n)).Nodup :=
    List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T)
  have hmap : (((hamPathsR n T).filter (lastBeats T n)).map (fun Q => Q ++ [n])).Nodup :=
    nodup_map_of_inj_on _ _ hL2 (fun a _ b _ h => List.append_cancel_right h)
  rw [length_eq_of_same_members _ _ (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR _ T))
    hmap (fun p => ?_), List.length_map]
  rw [List.mem_filter, mem_hamPathsR, List.mem_map]
  constructor
  · rintro ⟨⟨hl, hnd, hlt, hpath⟩, hh⟩
    cases p with
    | nil => exact Bool.noConfusion hh
    | cons a rest =>
      have hlast : lastOr a rest = n := Nat.eq_of_beq_eq_true hh
      cases rest with
      | nil =>
        exfalso
        rw [List.length_singleton] at hl
        omega
      | cons r rest' =>
        obtain ⟨mid, hmid⟩ := eq_append_lastOr (r :: rest') a (List.cons_ne_nil _ _)
        rw [hlast] at hmid
        have hp : a :: r :: rest' = (a :: mid) ++ [n] := by rw [hmid]; rfl
        rw [hp] at hnd hlt hpath hl
        have hnd' := hnd
        rw [List.nodup_append] at hnd'
        have hnQ : n ∉ a :: mid := fun h => hnd'.2.2 n h n List.mem_cons_self rfl
        refine ⟨a :: mid, ?_, hp.symm⟩
        rw [List.mem_filter, mem_hamPathsR]
        refine ⟨⟨?_, hnd'.1, fun y hy => lt_of_mem_ne n _ hlt y (List.mem_append_left _ hy)
          (fun h => hnQ (h ▸ hy)), (isPath_snoc_iff T a mid n).1 hpath |>.1⟩, ?_⟩
        · rw [List.length_append, List.length_singleton] at hl
          omega
        · exact ((isPath_snoc_iff T a mid n).1 hpath).2
  · rintro ⟨Q, hQ, rfl⟩
    rw [List.mem_filter, mem_hamPathsR] at hQ
    obtain ⟨⟨hl, hnd, hlt, hpath⟩, hh⟩ := hQ
    cases Q with
    | nil => exact Bool.noConfusion hh
    | cons q rest =>
      refine ⟨⟨by rw [List.length_append, List.length_singleton, hl], ?_, ?_,
        (isPath_snoc_iff T q rest n).2 ⟨hpath, hh⟩⟩, ?_⟩
      · rw [List.nodup_append]
        refine ⟨hnd, List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩, ?_⟩
        intro y hy z hz hyz
        rw [List.mem_singleton] at hz
        subst hz
        subst hyz
        exact Nat.lt_irrefl _ (hlt _ hy)
      · intro y hy
        rw [List.mem_append, List.mem_singleton] at hy
        rcases hy with hy | rfl
        · exact Nat.lt_succ_of_lt (hlt y hy)
        · exact Nat.lt_succ_self _
      · show Nat.beq (lastOr q (rest ++ [n])) n = true
        rw [lastOr_append, Nat.beq_refl]

/-! ## Forcing one pair -/

/-- `T` with the pair `{a, b}` oriented `a → b`. -/
def forceArc (T : Nat → Nat → Bool) (a b : Nat) (u v : Nat) : Bool :=
  if u = a ∧ v = b then true else if u = b ∧ v = a then false else T u v

theorem forceArc_ab (T : Nat → Nat → Bool) (a b : Nat) : forceArc T a b a b = true := by
  unfold forceArc
  rw [if_pos ⟨rfl, rfl⟩]

theorem forceArc_ba (T : Nat → Nat → Bool) (a b : Nat) (hab : a ≠ b) : forceArc T a b b a = false := by
  unfold forceArc
  rw [if_neg (fun h => hab h.2), if_pos ⟨rfl, rfl⟩]

theorem forceArc_other (T : Nat → Nat → Bool) (a b u v : Nat) (h1 : ¬ (u = a ∧ v = b))
    (h2 : ¬ (u = b ∧ v = a)) : forceArc T a b u v = T u v := by
  unfold forceArc
  rw [if_neg h1, if_neg h2]

theorem forceArc_isTournament (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) (a b : Nat)
    (hab : a ≠ b) : IsTournament n (forceArc T a b) := by
  intro u v hu hv huv
  by_cases h1 : u = a ∧ v = b
  · obtain ⟨rfl, rfl⟩ := h1
    rw [forceArc_ab, forceArc_ba T u v hab]
    rfl
  · by_cases h2 : u = b ∧ v = a
    · obtain ⟨rfl, rfl⟩ := h2
      rw [forceArc_ab, forceArc_ba T v u hab]
      rfl
    · rw [forceArc_other T a b v u (fun h => h2 ⟨h.2, h.1⟩) (fun h => h1 ⟨h.2, h.1⟩),
        forceArc_other T a b u v h1 h2]
      exact hT u v hu hv huv

theorem forceArc_sameOn (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) (a b : Nat)
    (ha : a < n) (hb : b < n) (hab : a ≠ b) (h : T a b = true) : SameOn n (forceArc T a b) T := by
  intro u v _ _ _
  by_cases h1 : u = a ∧ v = b
  · obtain ⟨rfl, rfl⟩ := h1
    rw [forceArc_ab, h]
  · by_cases h2 : u = b ∧ v = a
    · obtain ⟨rfl, rfl⟩ := h2
      rw [forceArc_ba T v u hab, hT v u ha hb hab, h]
      rfl
    · exact forceArc_other T a b u v h1 h2

/-! ## Vertices between `a` and `b` -/

/-- Remove `x` from a list. -/
def eraseN (x : Nat) (l : List Nat) : List Nat := l.filter (fun y => !Nat.beq y x)

theorem eraseN_of_not_mem (x : Nat) (l : List Nat) (h : x ∉ l) : eraseN x l = l := by
  unfold eraseN
  apply List.filter_eq_self.2
  intro y hy
  have : y ≠ x := fun e => h (e ▸ hy)
  rw [beq_false_of_ne y x this]
  rfl

theorem eraseN_cons_self (x : Nat) (l : List Nat) : eraseN x (x :: l) = eraseN x l := by
  unfold eraseN
  rw [filter_cons_false _ x l (by rw [Nat.beq_refl]; rfl)]

theorem eraseN_cons_ne (x y : Nat) (l : List Nat) (h : y ≠ x) : eraseN x (y :: l) = y :: eraseN x l := by
  unfold eraseN
  rw [filter_cons_true _ y l (by rw [beq_false_of_ne y x h]; rfl)]

theorem eraseN_append (x : Nat) (l1 l2 : List Nat) : eraseN x (l1 ++ l2) = eraseN x l1 ++ eraseN x l2 := by
  unfold eraseN
  rw [List.filter_append]

/-- A vertex with both a predecessor `a` and a successor `b` sits as `… a x b …`. -/
theorem split_around (x a b : Nat) (p : List Nat) (hnd : p.Nodup) (h1 : usesArc a x p = true)
    (h2 : usesArc x b p = true) : ∃ L R, p = L ++ a :: x :: b :: R := by
  obtain ⟨L1, R1, hp1⟩ := (usesArc_iff a x p).1 h1
  obtain ⟨L2, R2, hp2⟩ := (usesArc_iff x b p).1 h2
  have e : (L1 ++ [a]) ++ x :: R1 = L2 ++ x :: b :: R2 := by
    rw [List.append_assoc, ← hp2]
    exact hp1.symm
  have hnd' : ((L1 ++ [a]) ++ x :: R1).Nodup := by
    rw [List.append_assoc]
    exact hp1 ▸ hnd
  have hl := pos_unique x (L1 ++ [a]) R1 L2 (b :: R2) hnd' e
  have hR := (List.append_inj e hl).2
  rw [List.cons.injEq] at hR
  refine ⟨L1, R2, ?_⟩
  rw [hp1, hR.2]

theorem nodup_insert_middle (L R : List Nat) (a x b : Nat) (h : (L ++ a :: b :: R).Nodup)
    (hx : x ∉ L ++ a :: b :: R) : (L ++ a :: x :: b :: R).Nodup := by
  have hperm : (L ++ a :: x :: b :: R).Perm (x :: (L ++ a :: b :: R)) := by
    have := List.perm_middle (a := x) (l₁ := L ++ [a]) (l₂ := b :: R)
    simp only [List.append_assoc, List.cons_append, List.nil_append] at this
    exact this
  exact hperm.nodup_iff.2 (List.nodup_cons.2 ⟨hx, h⟩)

/-- In `L ++ [a] = l1 ++ u :: v :: l2`, the vertex `u` lies in `L`. -/
theorem mem_of_snoc_eq : ∀ (L : List Nat) (a u v : Nat) (l1 l2 : List Nat),
    L ++ [a] = l1 ++ u :: v :: l2 → u ∈ L
  | [], _, _, _, l1, l2, h => by
    exfalso
    have := congrArg List.length h
    rw [List.nil_append, List.length_singleton, length_append_two] at this
    omega
  | y :: L, a, u, v, [], l2, h => by
    rw [List.cons_append, List.nil_append, List.cons.injEq] at h
    rw [h.1]
    exact List.mem_cons_self
  | y :: L, a, u, v, z :: l1, l2, h => by
    rw [List.cons_append, List.cons_append, List.cons.injEq] at h
    exact List.mem_cons_of_mem y (mem_of_snoc_eq L a u v l1 l2 h.2)

/-- Paths through `a → b` in `forceArc T a b` and paths through `a → x → b` in `T`. -/
theorem isPath_splice (T : Nat → Nat → Bool) (L R : List Nat) (a x b : Nat)
    (hnd : (L ++ a :: b :: R).Nodup) :
    IsPath T (L ++ a :: x :: b :: R) ↔
      IsPath (forceArc T a b) (L ++ a :: b :: R) ∧ T a x = true ∧ T x b = true := by
  have haL : a ∉ L := fun h => by
    have := (List.nodup_append.1 hnd).2.2 a h a List.mem_cons_self
    exact this rfl
  have hbL : b ∉ L := fun h => by
    have := (List.nodup_append.1 hnd).2.2 b h b (List.mem_cons_of_mem a List.mem_cons_self)
    exact this rfl
  have hR : (a :: b :: R).Nodup := List.Nodup.sublist (List.sublist_append_right L _) hnd
  have haR : a ∉ R := fun h => (List.nodup_cons.1 hR).1 (List.mem_cons_of_mem b h)
  have hbR : b ∉ R := fun h => (List.nodup_cons.1 (List.nodup_cons.1 hR).2).1 h
  -- the two pieces avoid the pair `{a, b}`
  have hleft : IsPath (forceArc T a b) (L ++ [a]) ↔ IsPath T (L ++ [a]) := by
    apply isPath_congr_pairs
    intro l1 l2 u v hl
    have hu : u ∈ L := mem_of_snoc_eq L a u v l1 l2 hl
    exact forceArc_other T a b u v (fun h => haL (h.1 ▸ hu)) (fun h => hbL (h.1 ▸ hu))
  have hright : IsPath (forceArc T a b) (b :: R) ↔ IsPath T (b :: R) := by
    apply isPath_congr_pairs
    intro l1 l2 u v hl
    have hv : v ∈ R := by
      cases l1 with
      | nil =>
        rw [List.nil_append, List.cons.injEq] at hl
        rw [hl.2]
        exact List.mem_cons_self
      | cons y l1' =>
        rw [List.cons_append, List.cons.injEq] at hl
        rw [hl.2]
        exact List.mem_append_right _ (List.mem_cons_of_mem _ List.mem_cons_self)
    exact forceArc_other T a b u v (fun h => hbR (h.2 ▸ hv)) (fun h => haR (h.2 ▸ hv))
  rw [isPath_append_cons, isPath_append_cons (forceArc T a b)]
  show IsPath T (L ++ [a]) ∧ (T a x = true ∧ T x b = true ∧ IsPath T (b :: R)) ↔
    (IsPath (forceArc T a b) (L ++ [a]) ∧ (forceArc T a b a b = true ∧ IsPath (forceArc T a b) (b :: R))) ∧
      T a x = true ∧ T x b = true
  rw [hleft, hright, forceArc_ab]
  constructor
  · rintro ⟨h1, h2, h3, h4⟩
    exact ⟨⟨h1, rfl, h4⟩, h2, h3⟩
  · rintro ⟨⟨h1, _, h4⟩, h2, h3⟩
    exact ⟨h1, h2, h3, h4⟩

theorem eraseN_splice (L R : List Nat) (a x b : Nat) (hnd : (L ++ a :: x :: b :: R).Nodup)
    (hax : a ≠ x) (hbx : b ≠ x) : eraseN x (L ++ a :: x :: b :: R) = L ++ a :: b :: R := by
  have hxL : x ∉ L := fun h => by
    have := (List.nodup_append.1 hnd).2.2 x h x (List.mem_cons_of_mem a List.mem_cons_self)
    exact this rfl
  have hsuf : (x :: b :: R).Nodup :=
    (List.nodup_cons.1 (List.Nodup.sublist (List.sublist_append_right L _) hnd)).2
  have hxR : x ∉ R := fun h => (List.nodup_cons.1 hsuf).1 (List.mem_cons_of_mem b h)
  rw [eraseN_append, eraseN_of_not_mem x L hxL, eraseN_cons_ne x a _ hax, eraseN_cons_self,
    eraseN_cons_ne x b _ hbx, eraseN_of_not_mem x R hxR]

/-- Hamiltonian paths of `T` on `{0, …, n}` through `a → n → b` ↔ Hamiltonian paths of
`forceArc T a b` on `{0, …, n-1}` through `a → b`, when `a → n → b`. -/
theorem count_interior (n : Nat) (T : Nat → Nat → Bool) (a b : Nat) (ha : a < n) (hb : b < n)
    (hab : a ≠ b) :
    ((hamPathsR (n + 1) T).filter (fun p => usesArc a n p && usesArc n b p)).length =
      ind (T a n && T n b) * arcCount n (forceArc T a b) a b := by
  cases hx : (T a n && T n b)
  · rw [show ind false = 0 from rfl, Nat.zero_mul]
    have := length_filter_mono (fun p => usesArc a n p && usesArc n b p) (fun _ => false)
      (hamPathsR (n + 1) T) (fun p hp hpp => by
        rw [Bool.and_eq_true] at hpp
        have hpath := ((mem_hamPathsR _ T p).1 hp).2.2.2
        rw [usesArc_arc T a n p hpath hpp.1, usesArc_arc T n b p hpath hpp.2] at hx
        exact Bool.noConfusion hx)
    rw [length_filter_false] at this
    omega
  · rw [show ind true = 1 from rfl, Nat.one_mul]
    rw [Bool.and_eq_true] at hx
    unfold arcCount
    have hL1 : ((hamPathsR (n + 1) T).filter (fun p => usesArc a n p && usesArc n b p)).Nodup :=
      List.Nodup.sublist List.filter_sublist (nodup_hamPathsR _ T)
    have hL2 : ((hamPathsR n (forceArc T a b)).filter (usesArc a b)).Nodup :=
      List.Nodup.sublist List.filter_sublist (nodup_hamPathsR _ _)
    have hshape : ∀ p ∈ (hamPathsR (n + 1) T).filter (fun p => usesArc a n p && usesArc n b p),
        ∃ L R, p = L ++ a :: n :: b :: R ∧ IsHamPath (n + 1) T p := by
      intro p hp
      rw [List.mem_filter, mem_hamPathsR, Bool.and_eq_true] at hp
      obtain ⟨L, R, hLR⟩ := split_around n a b p hp.1.2.1 hp.2.1 hp.2.2
      exact ⟨L, R, hLR, hp.1⟩
    have hnot : a ≠ n ∧ b ≠ n := ⟨by omega, by omega⟩
    have hmapnd : (((hamPathsR (n + 1) T).filter
        (fun p => usesArc a n p && usesArc n b p)).map (eraseN n)).Nodup := by
      apply nodup_map_of_inj_on _ _ hL1
      intro p hp p' hp' heq
      obtain ⟨L, R, rfl, hham⟩ := hshape p hp
      obtain ⟨L', R', rfl, hham'⟩ := hshape p' hp'
      rw [eraseN_splice L R a n b hham.2.1 hnot.1 hnot.2,
        eraseN_splice L' R' a n b hham'.2.1 hnot.1 hnot.2] at heq
      have hndq : (L ++ a :: b :: R).Nodup := by
        rw [← eraseN_splice L R a n b hham.2.1 hnot.1 hnot.2]
        exact List.Nodup.sublist List.filter_sublist hham.2.1
      have hl := pos_unique a L (b :: R) L' (b :: R') hndq heq
      obtain ⟨h1, h2⟩ := List.append_inj heq hl
      rw [List.cons.injEq, List.cons.injEq] at h2
      rw [h1, h2.2.2]
    rw [← List.length_map (f := eraseN n)]
    apply length_eq_of_same_members _ _ hmapnd hL2
    intro q
    rw [List.mem_map, List.mem_filter, mem_hamPathsR]
    constructor
    · rintro ⟨p, hp, rfl⟩
      obtain ⟨L, R, rfl, hham⟩ := hshape p hp
      have hsub : (eraseN n (L ++ a :: n :: b :: R)).Sublist (L ++ a :: n :: b :: R) :=
        List.filter_sublist
      have hndq : (eraseN n (L ++ a :: n :: b :: R)).Nodup := List.Nodup.sublist hsub hham.2.1
      have hmemq : ∀ y ∈ eraseN n (L ++ a :: n :: b :: R), y < n := by
        intro y hy
        have hy' := hy
        unfold eraseN at hy'
        rw [List.mem_filter] at hy'
        have hyn : y ≠ n := fun h => by
          rw [h, Nat.beq_refl] at hy'
          exact Bool.noConfusion hy'.2
        exact lt_of_mem_ne n _ hham.2.2.1 y hy'.1 hyn
      rw [eraseN_splice L R a n b hham.2.1 hnot.1 hnot.2] at hndq hmemq ⊢
      refine ⟨⟨?_, hndq, hmemq, ((isPath_splice T L R a n b hndq).1 hham.2.2.2).1⟩, ?_⟩
      · have := hham.1
        simp only [List.length_append, List.length_cons] at this ⊢
        omega
      · exact (usesArc_iff a b _).2 ⟨L, R, rfl⟩
    · rintro ⟨⟨hl, hnd, hlt, hpath⟩, huse⟩
      obtain ⟨L, R, rfl⟩ := (usesArc_iff a b q).1 huse
      have hnq : n ∉ L ++ a :: b :: R := fun h => Nat.lt_irrefl n (hlt n h)
      have hndp := nodup_insert_middle L R a n b hnd hnq
      refine ⟨L ++ a :: n :: b :: R, ?_, eraseN_splice L R a n b hndp hnot.1 hnot.2⟩
      rw [List.mem_filter, mem_hamPathsR, Bool.and_eq_true]
      refine ⟨⟨?_, hndp, ?_, (isPath_splice T L R a n b hnd).2 ⟨hpath, hx.1, hx.2⟩⟩, ?_, ?_⟩
      · simp only [List.length_append, List.length_cons] at hl ⊢
        omega
      · intro y hy
        rw [List.mem_append, List.mem_cons, List.mem_cons, List.mem_cons] at hy
        rcases hy with hy | rfl | rfl | rfl | hy
        · exact Nat.lt_succ_of_lt (hlt y (List.mem_append_left _ hy))
        · exact Nat.lt_succ_of_lt ha
        · exact Nat.lt_succ_self _
        · exact Nat.lt_succ_of_lt hb
        · exact Nat.lt_succ_of_lt (hlt y (List.mem_append_right _
            (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ hy))))
      · exact (usesArc_iff a n _).2 ⟨L, b :: R, rfl⟩
      · exact (usesArc_iff n b _).2 ⟨L ++ [a], R, by simp only [List.append_assoc, List.cons_append,
          List.nil_append]⟩

/-! ## Where the top vertex sits -/

theorem lastIs_snoc (x z : Nat) : ∀ (l : List Nat), lastIs x (l ++ [z]) = Nat.beq z x
  | [] => rfl
  | y :: l => by
    show Nat.beq (lastOr y (l ++ [z])) x = Nat.beq z x
    rw [lastOr_append]

/-- In a duplicate-free `L ++ x :: R`, `u` precedes `x` iff `L` ends with `u`. -/
theorem usesArc_into (x u : Nat) (L R : List Nat) (hnd : (L ++ x :: R).Nodup) :
    usesArc u x (L ++ x :: R) = true ↔ ∃ L', L = L' ++ [u] := by
  constructor
  · intro h
    obtain ⟨l1, l2, hl⟩ := (usesArc_iff u x _).1 h
    have e : L ++ x :: R = (l1 ++ [u]) ++ x :: l2 := by
      rw [hl, List.append_assoc]
      rfl
    have hpos := pos_unique x L R (l1 ++ [u]) l2 hnd e
    exact ⟨l1, (List.append_inj e hpos).1⟩
  · rintro ⟨L', rfl⟩
    exact (usesArc_iff u x _).2 ⟨L', R, by rw [List.append_assoc]; rfl⟩

/-- In a duplicate-free `L ++ x :: R`, `w` follows `x` iff `R` starts with `w`. -/
theorem usesArc_from (x w : Nat) (L R : List Nat) (hnd : (L ++ x :: R).Nodup) :
    usesArc x w (L ++ x :: R) = true ↔ ∃ R', R = w :: R' := by
  constructor
  · intro h
    obtain ⟨l1, l2, hl⟩ := (usesArc_iff x w _).1 h
    have hpos := pos_unique x L R l1 (w :: l2) hnd hl
    have h2 := (List.append_inj hl hpos).2
    rw [List.cons.injEq] at h2
    exact ⟨l2, h2.2⟩
  · rintro ⟨R', rfl⟩
    exact (usesArc_iff x w _).2 ⟨L, R', rfl⟩

theorem ind_and_false_left (b : Bool) : ind (false && b) = 0 := rfl

theorem ind_and_false_right (a : Bool) : ind (a && false) = 0 := by
  cases a <;> rfl

/-- Each Hamiltonian path of `T` on `{0, …, n}` has `n` first, last, or between exactly one
pair `a, b`. -/
theorem hp_split (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) (p : List Nat)
    (hp : IsHamPath (n + 1) T p) :
    ind (headIs n p) + ind (lastIs n p) +
      rsum n (fun a => rsum n (fun b => ind (usesArc a n p && usesArc n b p))) = 1 := by
  obtain ⟨L, R, rfl⟩ := List.append_of_mem (mem_of_hamPath (n + 1) T p hp n (Nat.lt_succ_self n))
  have hnd := hp.2.1
  have hlt := hp.2.2.1
  have hnL : n ∉ L := fun h => by
    have := (List.nodup_append.1 hnd).2.2 n h n List.mem_cons_self
    exact this rfl
  have hnR : n ∉ R := fun h => (List.nodup_cons.1
    (List.Nodup.sublist (List.sublist_append_right L _) hnd)).1 h
  have hzero : ∀ (f : Nat → Nat → Nat), (∀ a b, f a b = 0) → rsum n (fun a => rsum n (fun b => f a b)) = 0 :=
    fun f hf => (rsum_congr n _ (fun _ => 0) (fun a _ =>
      (rsum_congr n _ (fun _ => 0) (fun b _ => hf a b)).trans (rsum_zero n))).trans (rsum_zero n)
  rcases List.eq_nil_or_concat L with hL | ⟨L', a0, hL⟩
  · -- `n` first
    subst hL
    cases R with
    | nil =>
      exfalso
      have := hp.1
      simp only [List.nil_append, List.length_singleton] at this
      omega
    | cons r R' =>
      have hhead : headIs n ([] ++ n :: r :: R') = true := Nat.beq_refl n
      obtain ⟨R'', z, hz⟩ : ∃ R'' z, r :: R' = R'' ++ [z] := by
        rcases List.eq_nil_or_concat (r :: R') with h | ⟨R'', z, h⟩
        · exact absurd h (List.cons_ne_nil _ _)
        · exact ⟨R'', z, by rw [h, List.concat_eq_append]⟩
      have hzR : z ∈ r :: R' := by rw [hz]; exact List.mem_append_right _ List.mem_cons_self
      have hlast : lastIs n ([] ++ n :: r :: R') = false := by
        rw [hz, show [] ++ n :: (R'' ++ [z]) = (n :: R'') ++ [z] from rfl, lastIs_snoc]
        exact beq_false_of_ne z n (fun h => hnR (h ▸ hzR))
      rw [hhead, hlast, hzero _ (fun a b => by
        have : usesArc a n ([] ++ n :: r :: R') = false := by
          cases h : usesArc a n ([] ++ n :: r :: R')
          · rfl
          · obtain ⟨L'', hL''⟩ := (usesArc_into n a [] (r :: R') hnd).1 h
            exact absurd (congrArg List.length hL'') (by
              rw [List.length_nil, List.length_append, List.length_singleton]
              omega)
        rw [this, ind_and_false_left])]
      rfl
  · rw [List.concat_eq_append] at hL
    subst hL
    have hhead : headIs n ((L' ++ [a0]) ++ n :: R) = false := by
      cases L' with
      | nil => exact beq_false_of_ne a0 n (fun h => hnL (h ▸ List.mem_cons_self))
      | cons y L'' =>
        exact beq_false_of_ne y n (fun h => hnL (h ▸ List.mem_cons_self))
    cases R with
    | nil =>
      have hlast : lastIs n ((L' ++ [a0]) ++ [n]) = true := by
        rw [lastIs_snoc, Nat.beq_refl]
      rw [hhead, hlast, hzero _ (fun a b => by
        have : usesArc n b ((L' ++ [a0]) ++ [n]) = false := by
          cases h : usesArc n b ((L' ++ [a0]) ++ [n])
          · rfl
          · obtain ⟨R'', hR''⟩ := (usesArc_from n b (L' ++ [a0]) [] hnd).1 h
            exact nomatch hR''
        rw [this, ind_and_false_right])]
      rfl
    | cons b0 R' =>
      obtain ⟨R'', z, hz⟩ : ∃ R'' z, b0 :: R' = R'' ++ [z] := by
        rcases List.eq_nil_or_concat (b0 :: R') with h | ⟨R'', z, h⟩
        · exact absurd h (List.cons_ne_nil _ _)
        · exact ⟨R'', z, by rw [h, List.concat_eq_append]⟩
      have hzR : z ∈ b0 :: R' := by rw [hz]; exact List.mem_append_right _ List.mem_cons_self
      have hlast : lastIs n ((L' ++ [a0]) ++ n :: b0 :: R') = false := by
        rw [hz, show (L' ++ [a0]) ++ n :: (R'' ++ [z]) = ((L' ++ [a0]) ++ n :: R'') ++ [z] by
          simp only [List.append_assoc, List.cons_append], lastIs_snoc]
        exact beq_false_of_ne z n (fun h => hnR (h ▸ hzR))
      have ha0 : a0 < n := lt_of_mem_ne n _ hlt a0
        (List.mem_append_left _ (List.mem_append_right _ List.mem_cons_self))
        (fun h => hnL (h ▸ List.mem_append_right _ List.mem_cons_self))
      have hb0 : b0 < n := lt_of_mem_ne n _ hlt b0
        (List.mem_append_right _ (List.mem_cons_of_mem _ List.mem_cons_self))
        (fun h => hnR (h ▸ List.mem_cons_self))
      have hterm : ∀ a b, ind (usesArc a n ((L' ++ [a0]) ++ n :: b0 :: R') &&
          usesArc n b ((L' ++ [a0]) ++ n :: b0 :: R')) =
          if a0 = a then (if b0 = b then 1 else 0) else 0 := by
        intro a b
        have e1 : usesArc a n ((L' ++ [a0]) ++ n :: b0 :: R') = true ↔ a0 = a := by
          rw [usesArc_into n a _ _ hnd]
          constructor
          · rintro ⟨L'', hL''⟩
            have := (List.append_inj' hL'' rfl).2
            exact (List.cons.inj this).1
          · rintro rfl
            exact ⟨L', rfl⟩
        have e2 : usesArc n b ((L' ++ [a0]) ++ n :: b0 :: R') = true ↔ b0 = b := by
          rw [usesArc_from n b _ _ hnd]
          constructor
          · rintro ⟨R'', hR''⟩
            exact (List.cons.inj hR'').1
          · rintro rfl
            exact ⟨R', rfl⟩
        by_cases h1 : a0 = a
        · rw [if_pos h1, (e1.2 h1)]
          by_cases h2 : b0 = b
          · rw [if_pos h2, e2.2 h2]
            rfl
          · rw [if_neg h2]
            cases h : usesArc n b ((L' ++ [a0]) ++ n :: b0 :: R')
            · rfl
            · exact absurd (e2.1 h) h2
        · rw [if_neg h1]
          cases h : usesArc a n ((L' ++ [a0]) ++ n :: b0 :: R')
          · rfl
          · exact absurd (e1.1 h) h1
      rw [hhead, hlast, rsum_congr n _ _ (fun a _ => rsum_congr n _ _ (fun b _ => hterm a b)),
        rsum_rsum_point n a0 b0 ha0 hb0]
      rfl

theorem sum_map_one' {α : Type} : ∀ (l : List α), (l.map (fun _ => 1)).sum = l.length
  | [] => rfl
  | _ :: t => by rw [List.map_cons, List.sum_cons, sum_map_one' t, List.length_cons]; omega

/-- **The decomposition of `H(T)` by the position of the top vertex.** -/
theorem hp_decomp (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) :
    hpCount (n + 1) T = ((hamPathsR (n + 1) T).filter (headIs n)).length +
      ((hamPathsR (n + 1) T).filter (lastIs n)).length +
      rsum n (fun a => rsum n (fun b =>
        ((hamPathsR (n + 1) T).filter (fun p => usesArc a n p && usesArc n b p)).length)) := by
  unfold hpCount
  rw [← sum_map_one' (hamPathsR (n + 1) T)]
  rw [List.map_congr_left (g := fun p => (ind (headIs n p) + ind (lastIs n p)) +
      rsum n (fun a => rsum n (fun b => ind (usesArc a n p && usesArc n b p))))
    (fun p hp => (hp_split n hn T p ((mem_hamPathsR _ T p).1 hp)).symm)]
  rw [sum_map_add, sum_map_add, length_filter_eq_sum, length_filter_eq_sum]
  congr 1
  rw [← rsum_listSum n (fun a p => rsum n (fun b => ind (usesArc a n p && usesArc n b p)))]
  apply rsum_congr
  intro a _
  rw [← rsum_listSum n (fun b p => ind (usesArc a n p && usesArc n b p))]
  apply rsum_congr
  intro b _
  rw [length_filter_eq_sum]

/-! ## Sums and parities -/

theorem rsum_comm (n : Nat) (f : Nat → Nat → Nat) :
    ∀ (m : Nat), rsum m (fun a => rsum n (fun b => f a b)) = rsum n (fun b => rsum m (fun a => f a b))
  | 0 => ((rsum_congr n _ (fun _ => 0) (fun _ _ => rfl)).trans (rsum_zero n)).symm
  | m + 1 => by
    show rsum m (fun a => rsum n (fun b => f a b)) + rsum n (fun b => f m b) =
      rsum n (fun b => rsum m (fun a => f a b) + f m b)
    rw [rsum_comm n f m, rsum_add]

theorem mul_list_sum {α : Type} (c : Nat) (f : α → Nat) :
    ∀ (l : List α), c * (l.map f).sum = (l.map (fun x => c * f x)).sum
  | [] => by rw [List.map_nil, List.map_nil, List.sum_nil, Nat.mul_zero]
  | x :: t => by rw [List.map_cons, List.map_cons, List.sum_cons, List.sum_cons, Nat.mul_add,
      mul_list_sum c f t]

theorem list_sum_mod_two {α : Type} (f g : α → Nat) :
    ∀ (l : List α), (∀ x ∈ l, f x % 2 = g x % 2) → (l.map f).sum % 2 = (l.map g).sum % 2
  | [], _ => rfl
  | x :: t, h => by
    rw [List.map_cons, List.map_cons, List.sum_cons, List.sum_cons]
    have h1 := h x List.mem_cons_self
    have h2 := list_sum_mod_two f g t (fun y hy => h y (List.mem_cons_of_mem x hy))
    omega

/-- Sum of `g` over the consecutive pairs of a list. -/
def pairSum (g : Nat → Nat → Nat) : List Nat → Nat
  | a :: b :: rest => g a b + pairSum g (b :: rest)
  | _ => 0

theorem rsum_rsum_weighted_point (n a b : Nat) (g : Nat → Nat → Nat) (ha : a < n) (hb : b < n) :
    rsum n (fun u => rsum n (fun v => g u v * (if a = u then (if b = v then 1 else 0) else 0))) = g a b := by
  have inner : ∀ u, rsum n (fun v => g u v * (if a = u then (if b = v then 1 else 0) else 0)) =
      if a = u then g u b else 0 := by
    intro u
    by_cases h : a = u
    · rw [if_pos h]
      rw [rsum_congr n _ (fun v => if b = v then g u v else 0) (fun v _ => by
        show g u v * (if a = u then (if b = v then 1 else 0) else 0) = (if b = v then g u v else 0)
        rw [if_pos h]
        by_cases h' : b = v
        · rw [if_pos h', if_pos h', Nat.mul_one]
        · rw [if_neg h', if_neg h', Nat.mul_zero])]
      exact rsum_single n b (fun v => g u v) hb
    · rw [if_neg h]
      exact (rsum_congr n _ (fun _ => 0) (fun v _ => by
        show g u v * (if a = u then (if b = v then 1 else 0) else 0) = 0
        rw [if_neg h, Nat.mul_zero])).trans (rsum_zero n)
  rw [rsum_congr n _ _ (fun u _ => inner u)]
  exact rsum_single n a (fun u => g u b) ha

/-- Weighted version of `arcs_of_path`. -/
theorem weighted_arcs (n : Nat) (g : Nat → Nat → Nat) : ∀ (p : List Nat), p.Nodup → (∀ x ∈ p, x < n) →
    rsum n (fun u => rsum n (fun v => g u v * ind (usesArc u v p))) = pairSum g p
  | [], _, _ => (rsum_congr n _ (fun _ => 0) (fun u _ =>
      (rsum_congr n _ (fun _ => 0) (fun v _ => Nat.mul_zero _)).trans (rsum_zero n))).trans (rsum_zero n)
  | [_], _, _ => (rsum_congr n _ (fun _ => 0) (fun u _ =>
      (rsum_congr n _ (fun _ => 0) (fun v _ => Nat.mul_zero _)).trans (rsum_zero n))).trans (rsum_zero n)
  | a :: b :: rest, hnd, hlt => by
    have ha : a ∉ b :: rest := (List.nodup_cons.1 hnd).1
    have ih := weighted_arcs n g (b :: rest) (List.nodup_cons.1 hnd).2
      (fun x hx => hlt x (List.mem_cons_of_mem a hx))
    show _ = g a b + pairSum g (b :: rest)
    rw [rsum_congr n _ (fun u => rsum n (fun v => g u v * (if a = u then (if b = v then 1 else 0) else 0)) +
        rsum n (fun v => g u v * ind (usesArc u v (b :: rest)))) (fun u _ => by
          dsimp only
          rw [← rsum_add]
          exact rsum_congr n _ _ (fun v _ => by rw [ind_usesArc_cons_cons u v a b rest ha, Nat.mul_add]))]
    rw [rsum_add, rsum_rsum_weighted_point n a b g (hlt a List.mem_cons_self)
      (hlt b (List.mem_cons_of_mem a List.mem_cons_self)), ih]

/-- `x` separates `a` and `b`: one of them beats `x`, which beats the other. -/
def cross (T : Nat → Nat → Bool) (x a b : Nat) : Nat := ind (T a x && T x b) + ind (T b x && T x a)

/-- Along a path of `T - x`, `[x → first] + #(separated consecutive pairs) + [last → x]` is odd. -/
theorem odd_contribution (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament (n + 1) T) :
    ∀ (rest : List Nat) (q : Nat), (∀ y ∈ q :: rest, y < n) →
      (ind (T n q) + pairSum (cross T n) (q :: rest) + ind (T (lastOr q rest) n)) % 2 = 1
  | [], q, hlt => by
    have hq := hlt q List.mem_cons_self
    have := hT n q (Nat.lt_succ_self n) (Nat.lt_succ_of_lt hq) (by omega)
    show (ind (T n q) + 0 + ind (T q n)) % 2 = 1
    rw [this]
    cases T n q <;> rfl
  | r :: rest, q, hlt => by
    have hq := hlt q List.mem_cons_self
    have hr := hlt r (List.mem_cons_of_mem q List.mem_cons_self)
    have ih := odd_contribution n T hT rest r (fun y hy => hlt y (List.mem_cons_of_mem q hy))
    show (ind (T n q) + (cross T n q r + pairSum (cross T n) (r :: rest)) +
      ind (T (lastOr r rest) n)) % 2 = 1
    have e : (ind (T n q) + cross T n q r) % 2 = ind (T n r) % 2 := by
      unfold cross
      rw [hT n q (Nat.lt_succ_self n) (Nat.lt_succ_of_lt hq) (by omega),
        hT n r (Nat.lt_succ_self n) (Nat.lt_succ_of_lt hr) (by omega)]
      cases T n q <;> cases T n r <;> rfl
    omega

/-- A pair cannot both precede and follow the same vertex in a duplicate-free list. -/
theorem not_pred_and_succ (a x : Nat) (p : List Nat) (hnd : p.Nodup) :
    (usesArc a x p && usesArc x a p) = false := by
  cases h1 : usesArc a x p
  · rfl
  · cases h2 : usesArc x a p
    · rfl
    · exfalso
      obtain ⟨L, R, hLR⟩ := split_around x a a p hnd h1 h2
      rw [hLR] at hnd
      have := List.Nodup.sublist (List.sublist_append_right L _) hnd
      rw [List.nodup_cons] at this
      exact this.1 (List.mem_cons_of_mem x List.mem_cons_self)

theorem deleteArc_force_sameOn (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) (a b : Nat)
    (ha : a < n) (hb : b < n) (hab : a ≠ b) (hba : T b a = true) :
    SameOn n (deleteArc (forceArc T a b) a b) (deleteArc T b a) := by
  intro u v hu hv huv
  unfold deleteArc
  by_cases h1 : u = a ∧ v = b
  · obtain ⟨rfl, rfl⟩ := h1
    rw [forceArc_ab, Nat.beq_refl, Nat.beq_refl, beq_false_of_ne u v hab]
    have : T u v = false := by rw [hT v u hb ha (Ne.symm hab), hba]; rfl
    rw [this]
    rfl
  · by_cases h2 : u = b ∧ v = a
    · obtain ⟨rfl, rfl⟩ := h2
      rw [forceArc_ba T v u hab, Nat.beq_refl, Nat.beq_refl]
      show false = (T u v && !(true && true))
      cases T u v <;> rfl
    · rw [forceArc_other T a b u v h1 h2]
      have e1 : (Nat.beq u a && Nat.beq v b) = false := by
        cases hua : Nat.beq u a
        · rfl
        · cases hvb : Nat.beq v b
          · rfl
          · exact absurd ⟨Nat.eq_of_beq_eq_true hua, Nat.eq_of_beq_eq_true hvb⟩ h1
      have e2 : (Nat.beq u b && Nat.beq v a) = false := by
        cases hub : Nat.beq u b
        · rfl
        · cases hva : Nat.beq v a
          · rfl
          · exact absurd ⟨Nat.eq_of_beq_eq_true hub, Nat.eq_of_beq_eq_true hva⟩ h2
      rw [e1, e2]

/-! ## Rédei's theorem -/

/-- **Rédei's theorem (1934).** Every tournament has an odd number of Hamiltonian paths. -/
theorem redei : ∀ (n : Nat) (T : Nat → Nat → Bool), IsTournament n T → hpCount n T % 2 = 1
  | 0, _, _ => rfl
  | 1, _, _ => rfl
  | m + 2, T, hT => by
    -- `n = m + 1 ≥ 1` remaining vertices, top vertex `x = m + 1`
    have hn : 1 ≤ m + 1 := Nat.succ_pos m
    have hT' : IsTournament (m + 1) T := fun a b ha hb hab => hT a b (by omega) (by omega) hab
    have ih := redei (m + 1) T hT'
    -- the interior terms, mod 2
    have hint : ∀ a b, a < m + 1 → b < m + 1 →
        ((hamPathsR (m + 2) T).filter (fun p => usesArc a (m + 1) p && usesArc (m + 1) b p)).length % 2 =
          (ind (T a (m + 1) && T (m + 1) b) *
            (arcCount (m + 1) T a b + arcCount (m + 1) T b a)) % 2 := by
      intro a b ha hb
      by_cases hab : a = b
      · subst hab
        have h0 : ((hamPathsR (m + 2) T).filter
            (fun p => usesArc a (m + 1) p && usesArc (m + 1) a p)).length = 0 := by
          have := length_filter_mono (fun p => usesArc a (m + 1) p && usesArc (m + 1) a p)
            (fun _ => false) (hamPathsR (m + 2) T) (fun p hp hpp => by
              have hpp' : (usesArc a (m + 1) p && usesArc (m + 1) a p) = true := hpp
              rw [not_pred_and_succ a (m + 1) p ((mem_hamPathsR _ T p).1 hp).2.1] at hpp'
              exact hpp')
          rw [length_filter_false] at this
          omega
        have hx : (T a (m + 1) && T (m + 1) a) = false := by
          rw [hT a (m + 1) (by omega) (Nat.lt_succ_self _) (by omega)]
          cases T a (m + 1) <;> rfl
        rw [h0, hx, show ind false = 0 from rfl, Nat.zero_mul]
      · rw [count_interior (m + 1) T a b ha hb hab]
        have hF : IsTournament (m + 1) (forceArc T a b) := forceArc_isTournament (m + 1) T hT' a b hab
        cases hTab : T a b
        · -- `b → a`: compare with reversing the pair
          have hba : T b a = true := by rw [hT' a b ha hb hab, hTab]; rfl
          have hzero : arcCount (m + 1) T a b = 0 := arcCount_eq_zero_of_not (m + 1) T a b (fun p hp => by
            cases h : usesArc a b p
            · rfl
            · rw [usesArc_arc T a b p hp.2.2.2 h] at hTab
              exact Bool.noConfusion hTab)
          have s1 := shave_arc (m + 1) (forceArc T a b) a b
          have s2 := shave_arc (m + 1) T b a
          have hc := hpCount_congr (m + 1) _ _ (deleteArc_force_sameOn (m + 1) T hT' a b ha hb hab hba)
          have o1 := redei (m + 1) (forceArc T a b) hF
          rw [hzero, Nat.zero_add]
          have key : arcCount (m + 1) (forceArc T a b) a b % 2 = arcCount (m + 1) T b a % 2 := by
            omega
          rw [Nat.mul_mod, key, ← Nat.mul_mod]
        · -- `a → b`: forcing changes nothing
          have hzero : arcCount (m + 1) T b a = 0 := arcCount_eq_zero_of_not (m + 1) T b a (fun p hp => by
            cases h : usesArc b a p
            · rfl
            · have := usesArc_arc T b a p hp.2.2.2 h
              rw [hT' a b ha hb hab, hTab] at this
              exact Bool.noConfusion this)
          rw [arcCount_congr (m + 1) _ _ (forceArc_sameOn (m + 1) T hT' a b ha hb hab hTab), hzero]
          rfl
    -- assemble
    rw [show m + 2 = (m + 1) + 1 from rfl, hp_decomp (m + 1) hn T, count_head (m + 1) hn T,
      count_last (m + 1) hn T]
    have hsum : rsum (m + 1) (fun a => rsum (m + 1) (fun b =>
        ((hamPathsR (m + 1 + 1) T).filter (fun p => usesArc a (m + 1) p && usesArc (m + 1) b p)).length)) % 2 =
        rsum (m + 1) (fun a => rsum (m + 1) (fun b =>
          cross T (m + 1) a b * arcCount (m + 1) T a b)) % 2 := by
      rw [rsum_mod_two_congr (m + 1) _ (fun a => rsum (m + 1) (fun b =>
        ind (T a (m + 1) && T (m + 1) b) * (arcCount (m + 1) T a b + arcCount (m + 1) T b a)))
        (fun a ha => rsum_mod_two_congr (m + 1) _ _ (fun b hb => hint a b ha hb))]
      congr 1
      rw [rsum_congr (m + 1) _ (fun a => rsum (m + 1) (fun b =>
          ind (T a (m + 1) && T (m + 1) b) * arcCount (m + 1) T a b) +
          rsum (m + 1) (fun b => ind (T a (m + 1) && T (m + 1) b) * arcCount (m + 1) T b a))
        (fun a _ => by dsimp only; rw [← rsum_add]; exact rsum_congr _ _ _ (fun b _ => Nat.mul_add _ _ _))]
      rw [rsum_add, rsum_comm (m + 1) (fun a b => ind (T a (m + 1) && T (m + 1) b) * arcCount (m + 1) T b a)]
      rw [← rsum_add]
      apply rsum_congr
      intro a _
      dsimp only
      rw [← rsum_add]
      apply rsum_congr
      intro b _
      unfold cross
      rw [Nat.add_mul]
    -- the cross sum as a sum over the Hamiltonian paths of `T - x`
    have hcross : rsum (m + 1) (fun a => rsum (m + 1) (fun b => cross T (m + 1) a b * arcCount (m + 1) T a b)) =
        ((hamPathsR (m + 1) T).map (fun Q => rsum (m + 1) (fun a => rsum (m + 1) (fun b =>
          cross T (m + 1) a b * ind (usesArc a b Q))))).sum := by
      rw [← rsum_listSum (m + 1) (fun a Q => rsum (m + 1) (fun b => cross T (m + 1) a b * ind (usesArc a b Q)))]
      apply rsum_congr
      intro a _
      rw [← rsum_listSum (m + 1) (fun b Q => cross T (m + 1) a b * ind (usesArc a b Q))]
      apply rsum_congr
      intro b _
      unfold arcCount
      rw [length_filter_eq_sum, mul_list_sum]
    have htotal : (((hamPathsR (m + 1) T).filter (beatsHead T (m + 1))).length +
        ((hamPathsR (m + 1) T).filter (lastBeats T (m + 1))).length +
        rsum (m + 1) (fun a => rsum (m + 1) (fun b => cross T (m + 1) a b * arcCount (m + 1) T a b))) % 2 =
        hpCount (m + 1) T % 2 := by
      rw [hcross, length_filter_eq_sum, length_filter_eq_sum, ← sum_map_add, ← sum_map_add]
      unfold hpCount
      rw [← sum_map_one' (hamPathsR (m + 1) T)]
      apply list_sum_mod_two
      intro Q hQ
      have hQ' := (mem_hamPathsR _ T Q).1 hQ
      cases Q with
      | nil =>
        exfalso
        have := hQ'.1
        rw [List.length_nil] at this
        omega
      | cons q rest =>
        rw [weighted_arcs (m + 1) (cross T (m + 1)) (q :: rest) hQ'.2.1 hQ'.2.2.1]
        show (ind (T (m + 1) q) + ind (T (lastOr q rest) (m + 1)) + pairSum (cross T (m + 1)) (q :: rest)) % 2 = 1 % 2
        have := odd_contribution (m + 1) T hT rest q hQ'.2.2.1
        omega
    omega

/-! ## Consequences -/

/-- **THM-4524 (C), unconditional.** If every arc of a tournament lies on an odd number of
Hamiltonian paths, then `n ≡ 1, 2 (mod 4)`. -/
theorem allOdd_mod_four_of_tournament (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool)
    (hT : IsTournament n T)
    (hodd : ∀ u v, u < n → v < n → u ≠ v → T u v = true → arcCount n T u v % 2 = 1) :
    n % 4 = 1 ∨ n % 4 = 2 :=
  allOdd_mod_four n hn T hT hodd (redei n T hT)

/-- **THM-4524 (C), unconditional.** If every arc lies on an even number of Hamiltonian
paths, then `n` is odd. -/
theorem allEven_odd_of_tournament (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool)
    (hT : IsTournament n T) (heven : ∀ u v, u < n → v < n → arcCount n T u v % 2 = 0) :
    n % 2 = 1 :=
  allEven_odd n hn T heven (redei n T hT)

/-- **THM-4524 (C), unconditional.** The number of arcs on an odd number of Hamiltonian
paths is `≡ n - 1 (mod 2)`. -/
theorem oddArcs_parity_of_tournament (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    rsum n (fun u => rsum n (fun v => arcCount n T u v % 2)) % 2 = (n - 1) % 2 :=
  oddArcs_parity n T (redei n T hT)

/-- For even `n ≥ 2`, some arc lies on an odd number of Hamiltonian paths. -/
theorem exists_odd_arc_of_even (n : Nat) (hn : 2 ≤ n) (heven : n % 2 = 0) (T : Nat → Nat → Bool)
    (hT : IsTournament n T) : ∃ u v, u < n ∧ v < n ∧ arcCount n T u v % 2 = 1 := by
  have h := oddArcs_parity_of_tournament n T hT
  by_cases hall : ∀ u, u < n → ∀ v, v < n → arcCount n T u v % 2 = 0
  · exfalso
    rw [rsum_congr n _ (fun _ => 0) (fun u hu =>
      (rsum_congr n _ (fun _ => 0) (fun v hv => hall u hu v hv)).trans (rsum_zero n)), rsum_zero] at h
    omega
  · -- a pair with an odd count exists (decidable search over `u, v < n`)
    have : ∃ u, u < n ∧ ∃ v, v < n ∧ arcCount n T u v % 2 = 1 := by
      apply Decidable.byContradiction
      intro hno
      apply hall
      intro u hu v hv
      have : ¬ arcCount n T u v % 2 = 1 := fun h' => hno ⟨u, hu, v, hv, h'⟩
      omega
    obtain ⟨u, hu, v, hv, huv⟩ := this
    exact ⟨u, v, hu, hv, huv⟩

/-- **Rédei shaving (THM-4524 §3).** Deleting the arc `u → v` keeps the number of
Hamiltonian paths odd iff `c(u → v)` is even. -/
theorem shave_keeps_odd_iff (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) (u v : Nat) :
    hpCount n (deleteArc T u v) % 2 = 1 ↔ arcCount n T u v % 2 = 0 := by
  have h1 := shave_arc n T u v
  have h2 := redei n T hT
  constructor
  · intro h
    omega
  · intro h
    omega

/-- **THM-4526, Theorem A for even `n`.** A tournament on an even number `n ≥ 2` of vertices
contains an odd number of copies of `H_n` (Hamiltonian path plus first→last arc). -/
theorem copiesH_odd_of_even (n : Nat) (hn : 2 ≤ n) (heven : n % 2 = 0) (T : Nat → Nat → Bool)
    (hT : IsTournament n T) : copiesH n T % 2 = 1 := by
  have h1 := copiesH_add n hn T hT
  have h2 := redei n T hT
  have h3 : (n * hcCount n T) % 2 = 0 := by
    rw [Nat.mul_mod, heven, Nat.zero_mul]
  omega

/-- In particular every tournament on an even number `n ≥ 2` of vertices contains `H_n`. -/
theorem containsH_of_even (n : Nat) (hn : 2 ≤ n) (heven : n % 2 = 0) (T : Nat → Nat → Bool)
    (hT : IsTournament n T) : ContainsH n T := by
  apply (containsH_iff_copiesH_pos n hn T).2
  have := copiesH_odd_of_even n hn heven T hT
  omega

end ProcgenSelfieEdim
