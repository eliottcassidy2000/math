import ProcgenSelfieEdim.CycleCount

set_option autoImplicit false

/-!
# THM-4524 C1 for cyclic groups: circulants of odd order are all-even

For a circulant digraph `circ n D` on `Z_n` (`a → b` iff `b - a mod n ∈ D`) with `n` odd,
every pair `(u, v)` lies on an even number of Hamiltonian paths as the arc `u → v`
(`circulant_arcCount_even`). This is THM-4524's C1 (Cayley tournaments of abelian groups
of odd order are all-even) for the cyclic groups `Z_n`; `D` need not even define a
tournament.

Proof (the note's): `φ(x) = u + v - x` is an involution of `Z_n` swapping `u` and `v` with
`T(φ a, φ b) = T(b, a)`. Then `P ↦ reverse (φ P)` is an involution on the Hamiltonian
paths through `u → v`; a fixed path would have `u → v` in the middle, impossible for `n`
odd. A fixed-point-free involution on a finite set has even size.
-/

namespace ProcgenSelfieEdim

/-- `(b - a) mod n` for `a, b < n`. -/
def cdist (n a b : Nat) : Nat := if a ≤ b then b - a else b + n - a

/-- The circulant digraph on `Z_n` with connection set `D`: `a → b` iff `(b - a) mod n ∈ D`. -/
def circ (n : Nat) (D : Nat → Bool) (a b : Nat) : Bool := D (cdist n a b)

/-- The reflection `x ↦ s - x (mod n)` for `x < n`, `s < 2n`. -/
def refl (n s x : Nat) : Nat := if x ≤ s then (if s - x < n then s - x else s - x - n) else s + n - x

theorem refl_spec (n s x : Nat) :
    (x ≤ s ∧ s - x < n ∧ refl n s x = s - x) ∨ (x ≤ s ∧ n ≤ s - x ∧ refl n s x = s - x - n) ∨
      (s < x ∧ refl n s x = s + n - x) := by
  unfold refl
  by_cases h1 : x ≤ s
  · rw [if_pos h1]
    by_cases h2 : s - x < n
    · rw [if_pos h2]
      exact Or.inl ⟨h1, h2, rfl⟩
    · rw [if_neg h2]
      exact Or.inr (Or.inl ⟨h1, by omega, rfl⟩)
  · rw [if_neg h1]
    exact Or.inr (Or.inr ⟨by omega, rfl⟩)

theorem cdist_spec (n a b : Nat) :
    (a ≤ b ∧ cdist n a b = b - a) ∨ (b < a ∧ cdist n a b = b + n - a) := by
  unfold cdist
  by_cases h : a ≤ b
  · rw [if_pos h]
    exact Or.inl ⟨h, rfl⟩
  · rw [if_neg h]
    exact Or.inr ⟨by omega, rfl⟩

theorem refl_lt (n s x : Nat) (hs : s < 2 * n) (hx : x < n) : refl n s x < n := by
  have := refl_spec n s x
  omega

theorem refl_refl (n s x : Nat) (hs : s < 2 * n) (hx : x < n) : refl n s (refl n s x) = x := by
  have h1 := refl_spec n s x
  have h2 := refl_spec n s (refl n s x)
  omega

theorem cdist_refl (n s a b : Nat) (hs : s < 2 * n) (ha : a < n) (hb : b < n) :
    cdist n (refl n s a) (refl n s b) = cdist n b a := by
  have h1 := refl_spec n s a
  have h2 := refl_spec n s b
  have h3 := cdist_spec n (refl n s a) (refl n s b)
  have h4 := cdist_spec n b a
  omega

/-- `φ` is an anti-automorphism of every circulant digraph. -/
theorem circ_refl (n : Nat) (D : Nat → Bool) (s a b : Nat) (hs : s < 2 * n) (ha : a < n) (hb : b < n) :
    circ n D (refl n s a) (refl n s b) = circ n D b a := by
  unfold circ
  rw [cdist_refl n s a b hs ha hb]

/-! ## Paths and reversal -/

/-- Paths through their consecutive pairs. -/
theorem isPath_iff_pairs (T : Nat → Nat → Bool) : ∀ (p : List Nat),
    IsPath T p ↔ ∀ (l1 l2 : List Nat) (a b : Nat), p = l1 ++ a :: b :: l2 → T a b = true
  | [] => ⟨fun _ l1 _ _ _ h => by
      have := congrArg List.length h
      rw [length_append_two, List.length_nil] at this
      omega, fun _ => trivial⟩
  | [x] => ⟨fun _ l1 _ _ _ h => by
      have := congrArg List.length h
      rw [length_append_two, List.length_singleton] at this
      omega, fun _ => trivial⟩
  | x :: y :: rest => by
    rw [isPath_cons_cons, isPath_iff_pairs T (y :: rest)]
    constructor
    · rintro ⟨hxy, hrest⟩ l1 l2 a b h
      cases l1 with
      | nil =>
        rw [List.nil_append, List.cons.injEq, List.cons.injEq] at h
        rw [← h.1, ← h.2.1]
        exact hxy
      | cons z l1 =>
        rw [List.cons_append, List.cons.injEq] at h
        exact hrest l1 l2 a b h.2
    · intro h
      exact ⟨h [] rest x y rfl, fun l1 l2 a b h' => h (x :: l1) l2 a b (by rw [h']; rfl)⟩

/-- The reverse of a path of the opposite relation is a path. -/
theorem isPath_reverse (T : Nat → Nat → Bool) (p : List Nat) (h : IsPath (fun x y => T y x) p) :
    IsPath T p.reverse := by
  rw [isPath_iff_pairs] at h ⊢
  intro l1 l2 a b hp
  have : p = l2.reverse ++ b :: a :: l1.reverse := by
    rw [← List.reverse_reverse p, hp]
    simp only [List.reverse_append, List.reverse_cons, List.append_assoc, List.cons_append,
      List.nil_append]
  exact h l2.reverse l1.reverse b a this

theorem nodup_reverse_of (l : List Nat) (h : l.Nodup) : l.reverse.Nodup := by
  unfold List.Nodup
  rw [List.pairwise_reverse]
  exact h.imp (fun h' => Ne.symm h')

/-! ## Fixed-point-free involutions -/

theorem length_filter_ne {α : Type} [DecidableEq α] (y : α) : ∀ (t : List α), y ∈ t → t.Nodup →
    (t.filter (fun z => decide (z ≠ y))).length + 1 = t.length
  | [], hy, _ => absurd hy List.not_mem_nil
  | a :: t, hy, hnd => by
    rw [List.nodup_cons] at hnd
    by_cases hay : a = y
    · subst hay
      rw [filter_cons_false _ a t (by simp)]
      have : t.filter (fun z => decide (z ≠ a)) = t := by
        apply List.filter_eq_self.2
        intro z hz
        exact decide_eq_true (fun hza => hnd.1 (hza ▸ hz))
      rw [this, List.length_cons]
    · have hy' : y ∈ t := by
        rcases List.mem_cons.1 hy with h | h
        · exact absurd h.symm hay
        · exact h
      rw [filter_cons_true _ a t (decide_eq_true hay), List.length_cons, List.length_cons,
        length_filter_ne y t hy' hnd.2]

/-- A duplicate-free list closed under a fixed-point-free involution has even length. -/
theorem even_of_involution (ι : List Nat → List Nat) : ∀ (k : Nat) (L : List (List Nat)),
    L.length ≤ k → L.Nodup → (∀ x ∈ L, ι x ∈ L) → (∀ x ∈ L, ι (ι x) = x) →
    (∀ x ∈ L, ι x ≠ x) → L.length % 2 = 0
  | _, [], _, _, _, _, _ => rfl
  | 0, _ :: _, hk, _, _, _, _ => by
    rw [List.length_cons] at hk
    omega
  | k + 1, x :: t, hk, hnd, hcl, hinv, hfix => by
    rw [List.nodup_cons] at hnd
    have hy : ι x ∈ t := by
      rcases List.mem_cons.1 (hcl x List.mem_cons_self) with h | h
      · exact absurd h (hfix x List.mem_cons_self)
      · exact h
    let L' := t.filter (fun z => decide (z ≠ ι x))
    have hlen : L'.length + 1 = t.length := length_filter_ne (ι x) t hy hnd.2
    have hL' : L'.length % 2 = 0 := by
      apply even_of_involution ι k L' (by rw [List.length_cons] at hk; omega)
        (List.Nodup.sublist List.filter_sublist hnd.2)
      · intro z hz
        rw [List.mem_filter] at hz
        have hzt := hz.1
        have hzy : z ≠ ι x := of_decide_eq_true hz.2
        have hzx : z ≠ x := fun h => hnd.1 (h ▸ hzt)
        have hiz := hcl z (List.mem_cons_of_mem x hzt)
        rw [List.mem_filter]
        rcases List.mem_cons.1 hiz with h | h
        · exfalso
          apply hzy
          rw [← hinv z (List.mem_cons_of_mem x hzt), h]
        · refine ⟨h, decide_eq_true (fun h' => hzx ?_)⟩
          rw [← hinv z (List.mem_cons_of_mem x hzt), h', hinv x List.mem_cons_self]
      · intro z hz
        exact hinv z (List.mem_cons_of_mem x (List.mem_filter.1 hz).1)
      · intro z hz
        exact hfix z (List.mem_cons_of_mem x (List.mem_filter.1 hz).1)
    rw [List.length_cons]
    omega

/-! ## The theorem -/

/-- **THM-4524 C1 for `Z_n`.** In every circulant digraph on `Z_n` with `n` odd, every
pair `(u, v)` is used as an arc `u → v` by an even number of Hamiltonian paths. -/
theorem circulant_arcCount_even (n : Nat) (hn : n % 2 = 1) (D : Nat → Bool) (u v : Nat)
    (hu : u < n) (hv : v < n) : arcCount n (circ n D) u v % 2 = 0 := by
  let T := circ n D
  let φ := refl n (u + v)
  have hs : u + v < 2 * n := by omega
  have hφu : φ u = v := by
    show refl n (u + v) u = v
    unfold refl
    rw [if_pos (by omega), if_pos (by omega)]
    omega
  have hφv : φ v = u := by
    show refl n (u + v) v = u
    unfold refl
    rw [if_pos (by omega), if_pos (by omega)]
    omega
  let ι : List Nat → List Nat := fun p => (p.map φ).reverse
  let L := (hamPathsR n T).filter (usesArc u v)
  have hmemL : ∀ p, p ∈ L ↔ IsHamPath n T p ∧ usesArc u v p = true := by
    intro p
    show p ∈ (hamPathsR n T).filter _ ↔ _
    rw [List.mem_filter, mem_hamPathsR]
  -- `ι` maps paths through `u → v` to paths through `u → v`
  have hmaps : ∀ p, IsHamPath n T p → usesArc u v p = true →
      IsHamPath n T (ι p) ∧ usesArc u v (ι p) = true := by
    intro p hp huv
    obtain ⟨hl, hnd, hlt, hpath⟩ := hp
    refine ⟨⟨?_, ?_, ?_, ?_⟩, ?_⟩
    · show (p.map φ).reverse.length = n
      rw [List.length_reverse, List.length_map, hl]
    · show (p.map φ).reverse.Nodup
      apply nodup_reverse_of
      exact nodup_map_of_inj_on φ p hnd (fun a ha b hb hab => by
        rw [← refl_refl n (u + v) a hs (hlt a ha), ← refl_refl n (u + v) b hs (hlt b hb)]
        exact congrArg (refl n (u + v)) hab)
    · intro x hx
      have hx' : x ∈ p.map φ := List.mem_reverse.1 hx
      rw [List.mem_map] at hx'
      obtain ⟨y, hy, rfl⟩ := hx'
      exact refl_lt n (u + v) y hs (hlt y hy)
    · apply isPath_reverse
      exact isPath_map n T (fun x y => T y x) φ
        (fun i j hi hj _ => circ_refl n D (u + v) j i hs hj hi) p hnd hlt hpath
    · obtain ⟨l1, l2, rfl⟩ := (usesArc_iff u v p).1 huv
      apply (usesArc_iff u v _).2
      refine ⟨(l2.map φ).reverse, (l1.map φ).reverse, ?_⟩
      show ((l1 ++ u :: v :: l2).map φ).reverse = _
      rw [List.map_append, List.map_cons, List.map_cons, hφu, hφv]
      simp only [List.reverse_append, List.reverse_cons, List.append_assoc, List.cons_append,
        List.nil_append]
  have hinv : ∀ p, IsHamPath n T p → ι (ι p) = p := by
    intro p hp
    show ((p.map φ).reverse.map φ).reverse = p
    rw [List.map_reverse, List.reverse_reverse, List.map_map]
    conv => rhs; rw [← List.map_id p]
    apply List.map_congr_left
    intro x hx
    exact refl_refl n (u + v) x hs (hp.2.2.1 x hx)
  have hfix : ∀ p, IsHamPath n T p → usesArc u v p = true → ι p ≠ p := by
    intro p hp huv heq
    obtain ⟨l1, l2, hp12⟩ := (usesArc_iff u v p).1 huv
    have hι : ι p = (l2.map φ).reverse ++ u :: v :: (l1.map φ).reverse := by
      show (p.map φ).reverse = _
      rw [hp12, List.map_append, List.map_cons, List.map_cons, hφu, hφv]
      simp only [List.reverse_append, List.reverse_cons, List.append_assoc, List.cons_append,
        List.nil_append]
    have hnd : (l1 ++ u :: v :: l2).Nodup := hp12 ▸ hp.2.1
    have hpos := pos_unique u l1 (v :: l2) (l2.map φ).reverse (v :: (l1.map φ).reverse) hnd
      (by rw [← hp12, ← heq, hι])
    rw [List.length_reverse, List.length_map] at hpos
    have hl := hp.1
    rw [hp12, length_append_two] at hl
    omega
  have := even_of_involution ι L.length L (Nat.le_refl _)
    (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T))
    (fun p hp => by
      obtain ⟨h1, h2⟩ := (hmemL p).1 hp
      exact (hmemL _).2 (hmaps p h1 h2))
    (fun p hp => hinv p ((hmemL p).1 hp).1)
    (fun p hp => hfix p ((hmemL p).1 hp).1 ((hmemL p).1 hp).2)
  exact this

/-- `QR_7` is the circulant on `Z_7` with connection set `{1, 2, 4}`. -/
theorem qr7_eq_circ : ∀ a, a < 7 → ∀ b, b < 7 →
    qr7 a b = circ 7 (fun k => Nat.beq k 1 || Nat.beq k 2 || Nat.beq k 4) a b := by
  decide

/-- Every arc of `QR_7` lies on an even number of Hamiltonian paths, as a consequence of the
circulant theorem (no enumeration; compare `qr7_counts`, where the count is 54). -/
theorem qr7_allEven (u v : Nat) (hu : u < 7) (hv : v < 7) : arcCount 7 qr7 u v % 2 = 0 := by
  have hsame : SameOn 7 qr7 (circ 7 (fun k => Nat.beq k 1 || Nat.beq k 2 || Nat.beq k 4)) :=
    fun a b ha hb _ => qr7_eq_circ a ha b hb
  rw [arcCount_congr 7 _ _ hsame]
  exact circulant_arcCount_even 7 rfl _ u v hu hv

end ProcgenSelfieEdim
