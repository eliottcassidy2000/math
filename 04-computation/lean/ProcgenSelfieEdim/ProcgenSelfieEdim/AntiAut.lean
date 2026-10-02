import ProcgenSelfieEdim.Redei

set_option autoImplicit false

/-!
# Anti-automorphisms and arc counts (THM-4529 infrastructure, and THM-4529 (c))

An *anti-automorphism* of a digraph `T` on `{0, …, n-1}` is a bijection `π` with
`T(π a, π b) = T(b, a)`: it maps `T` onto its reverse.

* `anti_maps_hp`: `P ↦ reverse(π P)` maps the Hamiltonian paths through `a → b` to the
  Hamiltonian paths through `π b → π a`.
* `arcCount_anti`: hence `c(a → b) = c(π b → π a)` when `π` has an inverse `π'`.
* `reflection_arcCount_even`: if moreover `π` is an involution swapping `u` and `v` and `n`
  is odd, then `c(u → v)` is even (no path is fixed: the arc would sit in the middle). This is
  the proof of THM-4529 Theorem 5.1 / THM-4524 C1, abstracted from the group.
* `anti_cycle_mod_four` (**THM-4529 §6.3, the necessity**): a tournament on `N ≥ 2` vertices
  with an anti-automorphism that is a single `N`-cycle has `N ≡ 2 (mod 4)`. For `N` odd,
  `π^N = id` would reverse every arc; for `4 ∣ N`, `π^(N/2)` is an automorphism exchanging
  `0` and `π^(N/2)(0) ≠ 0`, so it would reverse the arc between them.
-/

namespace ProcgenSelfieEdim

/-- An anti-automorphism maps the Hamiltonian paths through `a → b` to Hamiltonian paths
through `π b → π a`. -/
theorem anti_maps_hp (n : Nat) (T : Nat → Nat → Bool) (π : Nat → Nat)
    (hπ : ∀ x, x < n → π x < n) (hinj : ∀ x y, x < n → y < n → π x = π y → x = y)
    (hanti : ∀ a b, a < n → b < n → a ≠ b → T (π a) (π b) = T b a)
    (a b : Nat) (p : List Nat) (hp : IsHamPath n T p) (hab : usesArc a b p = true) :
    IsHamPath n T (p.map π).reverse ∧ usesArc (π b) (π a) (p.map π).reverse = true := by
  obtain ⟨hl, hnd, hlt, hpath⟩ := hp
  refine ⟨⟨?_, ?_, ?_, ?_⟩, ?_⟩
  · rw [List.length_reverse, List.length_map, hl]
  · apply nodup_reverse_of
    exact nodup_map_of_inj_on π p hnd (fun x hx y hy hxy => hinj x y (hlt x hx) (hlt y hy) hxy)
  · intro x hx
    have hx' : x ∈ p.map π := List.mem_reverse.1 hx
    rw [List.mem_map] at hx'
    obtain ⟨y, hy, rfl⟩ := hx'
    exact hπ y (hlt y hy)
  · apply isPath_reverse
    exact isPath_map n T (fun x y => T y x) π
      (fun i j hi hj hij => hanti j i hj hi (Ne.symm hij)) p hnd hlt hpath
  · obtain ⟨l1, l2, rfl⟩ := (usesArc_iff a b p).1 hab
    apply (usesArc_iff _ _ _).2
    refine ⟨(l2.map π).reverse, (l1.map π).reverse, ?_⟩
    rw [List.map_append, List.map_cons, List.map_cons]
    simp only [List.reverse_append, List.reverse_cons, List.append_assoc, List.cons_append,
      List.nil_append]

/-- Applying `π'` then `π` (entrywise, with a reversal each time) gives back the path. -/
theorem map_reverse_map (π π' : Nat → Nat) (n : Nat) (p : List Nat) (hlt : ∀ x ∈ p, x < n)
    (h : ∀ x, x < n → π (π' x) = x) : (((p.map π').reverse).map π).reverse = p := by
  rw [List.map_reverse, List.reverse_reverse, List.map_map]
  conv => rhs; rw [← List.map_id p]
  apply List.map_congr_left
  intro x hx
  exact h x (hlt x hx)

/-- **Arc counts are invariant under anti-automorphisms**: `c(a → b) = c(π b → π a)`. -/
theorem arcCount_anti (n : Nat) (T : Nat → Nat → Bool) (π π' : Nat → Nat)
    (hπ : ∀ x, x < n → π x < n) (hπ' : ∀ x, x < n → π' x < n)
    (h1 : ∀ x, x < n → π' (π x) = x) (h2 : ∀ x, x < n → π (π' x) = x)
    (hanti : ∀ a b, a < n → b < n → a ≠ b → T (π a) (π b) = T b a)
    (a b : Nat) (ha : a < n) (hb : b < n) : arcCount n T a b = arcCount n T (π b) (π a) := by
  have hinj : ∀ x y, x < n → y < n → π x = π y → x = y := fun x y hx hy hxy => by
    rw [← h1 x hx, ← h1 y hy, hxy]
  have hinj' : ∀ x y, x < n → y < n → π' x = π' y → x = y := fun x y hx hy hxy => by
    rw [← h2 x hx, ← h2 y hy, hxy]
  have hanti' : ∀ x y, x < n → y < n → x ≠ y → T (π' x) (π' y) = T y x := by
    intro x y hx hy hxy
    have hne : π' y ≠ π' x := fun h => hxy (hinj' y x hy hx h).symm
    have := hanti (π' y) (π' x) (hπ' y hy) (hπ' x hx) hne
    rw [h2 y hy, h2 x hx] at this
    exact this.symm
  let ι : List Nat → List Nat := fun p => (p.map π).reverse
  let L1 := (hamPathsR n T).filter (usesArc a b)
  have hmem1 : ∀ p, p ∈ L1 ↔ IsHamPath n T p ∧ usesArc a b p = true := by
    intro p
    show p ∈ (hamPathsR n T).filter _ ↔ _
    rw [List.mem_filter, mem_hamPathsR]
  have hnd1 : L1.Nodup := List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T)
  have hlen : (L1.map ι).length = arcCount n T (π b) (π a) := by
    apply arcCount_unique
    · apply nodup_map_of_inj_on ι L1 hnd1
      intro p hp q hq hpq
      have hp' := ((hmem1 p).1 hp).1
      have hq' := ((hmem1 q).1 hq).1
      rw [← map_reverse_map π' π n p hp'.2.2.1 h1, ← map_reverse_map π' π n q hq'.2.2.1 h1]
      show ((ι p).map π').reverse = ((ι q).map π').reverse
      rw [hpq]
    · intro q
      rw [List.mem_map]
      constructor
      · rintro ⟨p, hp, rfl⟩
        obtain ⟨hp1, hp2⟩ := (hmem1 p).1 hp
        exact anti_maps_hp n T π hπ hinj hanti a b p hp1 hp2
      · rintro ⟨hq1, hq2⟩
        refine ⟨(q.map π').reverse, (hmem1 _).2 ?_, map_reverse_map π π' n q hq1.2.2.1 h2⟩
        have := anti_maps_hp n T π' hπ' hinj' hanti' (π b) (π a) q hq1 hq2
        rw [h1 a ha, h1 b hb] at this
        exact this
  rw [← hlen, List.length_map]
  rfl

/-- **Reflection lemma.** If `φ` is an involution of `{0, …, n-1}` reversing `T` and swapping
`u` and `v`, and `n` is odd, then the arc `u → v` lies on an even number of Hamiltonian
paths. -/
theorem reflection_arcCount_even (n : Nat) (hn : n % 2 = 1) (T : Nat → Nat → Bool) (φ : Nat → Nat)
    (hφ : ∀ x, x < n → φ x < n) (hinv : ∀ x, x < n → φ (φ x) = x)
    (hanti : ∀ a b, a < n → b < n → a ≠ b → T (φ a) (φ b) = T b a)
    (u v : Nat) (hu : φ u = v) (hv : φ v = u) : arcCount n T u v % 2 = 0 := by
  have hinj : ∀ x y, x < n → y < n → φ x = φ y → x = y := fun x y hx hy hxy => by
    rw [← hinv x hx, ← hinv y hy, hxy]
  let ι : List Nat → List Nat := fun p => (p.map φ).reverse
  let L := (hamPathsR n T).filter (usesArc u v)
  have hmemL : ∀ p, p ∈ L ↔ IsHamPath n T p ∧ usesArc u v p = true := by
    intro p
    show p ∈ (hamPathsR n T).filter _ ↔ _
    rw [List.mem_filter, mem_hamPathsR]
  have hmaps : ∀ p, IsHamPath n T p → usesArc u v p = true →
      IsHamPath n T (ι p) ∧ usesArc u v (ι p) = true := by
    intro p hp huv
    have := anti_maps_hp n T φ hφ hinj hanti u v p hp huv
    rw [hu, hv] at this
    exact this
  have hfix : ∀ p, IsHamPath n T p → usesArc u v p = true → ι p ≠ p := by
    intro p hp huv heq
    obtain ⟨l1, l2, hp12⟩ := (usesArc_iff u v p).1 huv
    have hι : ι p = (l2.map φ).reverse ++ u :: v :: (l1.map φ).reverse := by
      show (p.map φ).reverse = _
      rw [hp12, List.map_append, List.map_cons, List.map_cons, hu, hv]
      simp only [List.reverse_append, List.reverse_cons, List.append_assoc, List.cons_append,
        List.nil_append]
    have hnd : (l1 ++ u :: v :: l2).Nodup := hp12 ▸ hp.2.1
    have hpos := pos_unique u l1 (v :: l2) (l2.map φ).reverse (v :: (l1.map φ).reverse) hnd
      (by rw [← hp12, ← heq, hι])
    rw [List.length_reverse, List.length_map] at hpos
    have hl := hp.1
    rw [hp12, length_append_two] at hl
    omega
  exact even_of_involution ι L.length L (Nat.le_refl _)
    (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T))
    (fun p hp => by
      obtain ⟨h1, h2⟩ := (hmemL p).1 hp
      exact (hmemL _).2 (hmaps p h1 h2))
    (fun p hp => map_reverse_map φ φ n p ((hmemL p).1 hp).1.2.2.1 hinv)
    (fun p hp => hfix p ((hmemL p).1 hp).1 ((hmemL p).1 hp).2)

/-! ## THM-4529 §6.3: a cyclic anti-automorphism forces `N ≡ 2 (mod 4)` -/

/-- `f` iterated `k` times. -/
def iterN (f : Nat → Nat) : Nat → Nat → Nat
  | 0, x => x
  | k + 1, x => f (iterN f k x)

theorem iterN_succ' (f : Nat → Nat) : ∀ (j y : Nat), f (iterN f j y) = iterN f j (f y)
  | 0, _ => rfl
  | j + 1, y => by
    show f (f (iterN f j y)) = f (iterN f j (f y))
    rw [iterN_succ' f j y]

theorem iterN_add (f : Nat → Nat) (j : Nat) : ∀ (k x : Nat), iterN f (j + k) x = iterN f j (iterN f k x)
  | 0, _ => rfl
  | k + 1, x => by
    show f (iterN f (j + k) x) = iterN f j (f (iterN f k x))
    rw [iterN_add f j k x, iterN_succ' f j (iterN f k x)]

theorem iterN_lt (N : Nat) (f : Nat → Nat) (hf : ∀ a, a < N → f a < N) :
    ∀ (k a : Nat), a < N → iterN f k a < N
  | 0, _, h => h
  | k + 1, a, h => hf _ (iterN_lt N f hf k a h)

/-- **THM-4529 §6.3.** A tournament on `N ≥ 2` vertices with an anti-automorphism `π` that is
a single `N`-cycle (`π^N = id`, and `π^k(0) ≠ 0` for `0 < k < N`) has `N ≡ 2 (mod 4)`. -/
theorem anti_cycle_mod_four (N : Nat) (hN : 2 ≤ N) (T : Nat → Nat → Bool) (hT : IsTournament N T)
    (π : Nat → Nat) (hπ : ∀ a, a < N → π a < N)
    (hanti : ∀ a b, a < N → b < N → a ≠ b → T (π a) (π b) = T b a)
    (hcyc : ∀ a, a < N → iterN π N a = a)
    (hfree : ∀ k, 0 < k → k < N → iterN π k 0 ≠ 0) : N % 4 = 2 := by
  -- `π` is injective: otherwise `T(b, a) = T(π a, π a) = T(a, b)`
  have hinj : ∀ a b, a < N → b < N → π a = π b → a = b := by
    intro a b ha hb hab
    by_cases h : a = b
    · exact h
    · exfalso
      have e1 := hanti a b ha hb h
      have e2 := hanti b a hb ha (Ne.symm h)
      rw [hab] at e1 e2
      rw [e1] at e2
      have := hT a b ha hb h
      rw [e2] at this
      cases hc : T a b <;> rw [hc] at this <;> exact Bool.noConfusion this
  have hinjk : ∀ k a b, a < N → b < N → a ≠ b → iterN π k a ≠ iterN π k b := by
    intro k
    induction k with
    | zero => intro a b _ _ h; exact h
    | succ k ih =>
      intro a b ha hb h e
      exact ih a b ha hb h (hinj _ _ (iterN_lt N π hπ k a ha) (iterN_lt N π hπ k b hb) e)
  -- iterates: even powers preserve arcs, odd powers reverse them
  have hiter : ∀ k a b, a < N → b < N → a ≠ b →
      T (iterN π k a) (iterN π k b) = (if k % 2 = 0 then T a b else T b a) := by
    intro k
    induction k with
    | zero => intro a b _ _ _; rfl
    | succ k ih =>
      intro a b ha hb h
      show T (π (iterN π k a)) (π (iterN π k b)) = _
      rw [hanti _ _ (iterN_lt N π hπ k a ha) (iterN_lt N π hπ k b hb) (hinjk k a b ha hb h),
        ih b a hb ha (Ne.symm h)]
      by_cases hk : k % 2 = 0
      · rw [if_pos hk, if_neg (by omega)]
      · rw [if_neg hk, if_pos (by omega)]
  have hrev : ∀ a b, a < N → b < N → a ≠ b → T a b = T b a → False := by
    intro a b ha hb h e
    have := hT a b ha hb h
    rw [← e] at this
    cases hc : T a b <;> rw [hc] at this <;> exact Bool.noConfusion this
  -- `N` is even
  have heven : N % 2 = 0 := by
    by_cases hodd : N % 2 = 0
    · exact hodd
    · exfalso
      have := hiter N 0 1 (by omega) (by omega) (by decide)
      rw [hcyc 0 (by omega), hcyc 1 (by omega), if_neg hodd] at this
      exact hrev 0 1 (by omega) (by omega) (by decide) this
  -- `4 ∤ N`
  by_cases h4 : N % 4 = 0
  · exfalso
    obtain ⟨j, hj⟩ : ∃ j, N = j + j := ⟨N / 2, by omega⟩
    have hb := hfree j (by omega) (by omega)
    have hbN := iterN_lt N π hπ j 0 (by omega)
    have hback : iterN π j (iterN π j 0) = 0 := by
      rw [← iterN_add, ← hj]
      exact hcyc 0 (by omega)
    have := hiter j 0 (iterN π j 0) (by omega) hbN (Ne.symm hb)
    rw [hback, if_pos (by omega)] at this
    exact hrev 0 (iterN π j 0) (by omega) hbN (Ne.symm hb) this.symm
  · omega

end ProcgenSelfieEdim
