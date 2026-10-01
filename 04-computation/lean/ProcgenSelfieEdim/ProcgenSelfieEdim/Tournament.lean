set_option autoImplicit false

/-!
# Tournaments, switching, and the selfie gauge theorem (THM-4524, A1)

A tournament on `{0, …, N-1}` is a Boolean relation `T` (`T a b = true` means the
arc `a → b`) in which exactly one of `T a b`, `T b a` holds for distinct `a, b < N`.
Values on the diagonal and outside `{0, …, N-1}` play no role: two relations are the
same labeled tournament when they agree on all ordered pairs of distinct vertices
(`SameOn`).

**Tiling model.** The base path is `N-1 → N-2 → ⋯ → 0` (arcs `i+1 → i`; this is the
note's `N → ⋯ → 1` shifted down by one). The `C(N-1, 2)` non-path pairs `a + 2 ≤ b`
carry free tile bits: `t a b = true` means `a → b`.

**Selfie model.** A selfie tournament is a pair `(t, L)`: tile bits and a loop set
`L ⊆ {0, …, N-1}` (here `L : Nat → Bool`). Its labeled tournament is
`switch L (tiling t)`, where switching reverses every arc with exactly one end in `L`.

**Theorem A1** (proved below): the base-path arc `i+1 → i` is reversed iff exactly
one of `i, i+1` is looped; the map `(t, L) ↦ switch L (tiling t)` is onto the labeled
tournaments; and its fibres are exactly `{(t, L), (t, Lᶜ)}`, two distinct points.
-/

namespace ProcgenSelfieEdim

/-- `T` is a tournament on `{0, …, N-1}`. -/
def IsTournament (N : Nat) (T : Nat → Nat → Bool) : Prop :=
  ∀ a b, a < N → b < N → a ≠ b → T b a = !T a b

/-- `T` and `T'` are the same labeled tournament on `{0, …, N-1}`. -/
def SameOn (N : Nat) (T T' : Nat → Nat → Bool) : Prop :=
  ∀ a b, a < N → b < N → a ≠ b → T a b = T' a b

/-- Switching at `L` reverses every arc with exactly one end in `L`. -/
def switch (L : Nat → Bool) (T : Nat → Nat → Bool) (a b : Nat) : Bool :=
  if L a = L b then T a b else T b a

/-- The tiling tournament: base path arcs `i+1 → i`; for `a + 2 ≤ b` the tile bit
`t a b` says whether `a → b`. -/
def tiling (t : Nat → Nat → Bool) (a b : Nat) : Bool :=
  if a = b + 1 then true
  else if b = a + 1 then false
  else if a < b then t a b
  else !t b a

/-- Two tilings with the same tile bits. -/
def SameTiles (N : Nat) (t t' : Nat → Nat → Bool) : Prop :=
  ∀ a b, a + 2 ≤ b → b < N → t a b = t' a b

/-- Two loop sets with the same members in `{0, …, N-1}`. -/
def SameLoops (N : Nat) (L L' : Nat → Bool) : Prop := ∀ a, a < N → L a = L' a

/-- The complementary loop set. -/
def compl (L : Nat → Bool) (a : Nat) : Bool := !L a

/-! ## Basic facts -/

theorem tiling_down (t : Nat → Nat → Bool) (i : Nat) : tiling t (i + 1) i = true := by
  unfold tiling
  rw [if_pos rfl]

theorem tiling_up (t : Nat → Nat → Bool) (i : Nat) : tiling t i (i + 1) = false := by
  unfold tiling
  rw [if_neg (by omega), if_pos rfl]

theorem tiling_lt (t : Nat → Nat → Bool) (a b : Nat) (h : a + 2 ≤ b) : tiling t a b = t a b := by
  unfold tiling
  rw [if_neg (by omega), if_neg (by omega), if_pos (by omega)]

theorem tiling_gt (t : Nat → Nat → Bool) (a b : Nat) (h : b + 2 ≤ a) : tiling t a b = !t b a := by
  unfold tiling
  rw [if_neg (by omega), if_neg (by omega), if_neg (by omega)]

/-- Distinct `a, b` are in one of four positions. -/
theorem four_cases (a b : Nat) (hab : a ≠ b) :
    a = b + 1 ∨ b = a + 1 ∨ a + 2 ≤ b ∨ b + 2 ≤ a := by
  by_cases h1 : a = b + 1
  · exact Or.inl h1
  · by_cases h2 : b = a + 1
    · exact Or.inr (Or.inl h2)
    · by_cases h3 : a + 2 ≤ b
      · exact Or.inr (Or.inr (Or.inl h3))
      · exact Or.inr (Or.inr (Or.inr (by omega)))

theorem tiling_isTournament (N : Nat) (t : Nat → Nat → Bool) : IsTournament N (tiling t) := by
  intro a b _ _ hab
  rcases four_cases a b hab with h | h | h | h
  · subst h
    rw [tiling_down, tiling_up]
    rfl
  · subst h
    rw [tiling_down, tiling_up]
    rfl
  · rw [tiling_gt t b a h, tiling_lt t a b h]
  · rw [tiling_lt t b a h, tiling_gt t a b h, Bool.not_not]

theorem switch_isTournament (N : Nat) (L : Nat → Bool) (T : Nat → Nat → Bool)
    (hT : IsTournament N T) : IsTournament N (switch L T) := by
  intro a b ha hb hab
  unfold switch
  by_cases h : L a = L b
  · rw [if_pos h.symm, if_pos h]
    exact hT a b ha hb hab
  · rw [if_neg (Ne.symm h), if_neg h, hT b a hb ha (Ne.symm hab)]

/-- Switching twice at the same set is the identity. -/
theorem switch_switch (L : Nat → Bool) (T : Nat → Nat → Bool) (a b : Nat) :
    switch L (switch L T) a b = T a b := by
  unfold switch
  by_cases h : L a = L b
  · rw [if_pos h, if_pos h]
  · rw [if_neg h, if_neg (Ne.symm h)]

theorem not_eq_not (x y : Bool) : ((!x) = !y) ↔ x = y := by
  cases x <;> cases y <;> decide

/-- Switching at `L` and at its complement agree. -/
theorem switch_compl (L : Nat → Bool) (T : Nat → Nat → Bool) (a b : Nat) :
    switch (compl L) T a b = switch L T a b := by
  unfold switch compl
  by_cases h : L a = L b
  · rw [if_pos h, if_pos ((not_eq_not _ _).2 h)]
  · rw [if_neg h, if_neg (fun h' => h ((not_eq_not _ _).1 h'))]

theorem switch_congr (N : Nat) (L : Nat → Bool) (T T' : Nat → Nat → Bool) (h : SameOn N T T') :
    SameOn N (switch L T) (switch L T') := by
  intro a b ha hb hab
  unfold switch
  by_cases hl : L a = L b
  · rw [if_pos hl, if_pos hl]
    exact h a b ha hb hab
  · rw [if_neg hl, if_neg hl]
    exact h b a hb ha (Ne.symm hab)

theorem switch_loops_congr (N : Nat) (L L' : Nat → Bool) (T : Nat → Nat → Bool)
    (h : SameLoops N L L') : SameOn N (switch L T) (switch L' T) := by
  intro a b ha hb _
  unfold switch
  rw [h a ha, h b hb]

/-- Tiling tournaments are determined by their tile bits. -/
theorem tiling_congr (N : Nat) (t t' : Nat → Nat → Bool) (h : SameTiles N t t') :
    SameOn N (tiling t) (tiling t') := by
  intro a b ha hb hab
  rcases four_cases a b hab with h' | h' | h' | h'
  · subst h'
    rw [tiling_down, tiling_down]
  · subst h'
    rw [tiling_up, tiling_up]
  · rw [tiling_lt t a b h', tiling_lt t' a b h']
    exact h a b h' hb
  · rw [tiling_gt t a b h', tiling_gt t' a b h', h b a h' ha]

/-! ## Theorem A1 -/

/-- **A1, derivative law.** In `switch L (tiling t)`, the base-path arc `i+1 → i` is
present (not reversed) iff `i` and `i+1` are both looped or both unlooped. -/
theorem path_arc (L : Nat → Bool) (t : Nat → Nat → Bool) (i : Nat) :
    switch L (tiling t) (i + 1) i = true ↔ L i = L (i + 1) := by
  unfold switch
  by_cases h : L (i + 1) = L i
  · rw [if_pos h, tiling_down]
    exact ⟨fun _ => h.symm, fun _ => rfl⟩
  · rw [if_neg h, tiling_up]
    exact ⟨fun h' => Bool.noConfusion h', fun h' => absurd h'.symm h⟩

/-- Recovering the tiles from a switched tiling tournament with known loops. -/
theorem tiles_of_switch_eq (N : Nat) (L : Nat → Bool) (t t' : Nat → Nat → Bool)
    (h : SameOn N (switch L (tiling t)) (switch L (tiling t'))) : SameTiles N t t' := by
  intro a b hab hb
  have hT : SameOn N (tiling t) (tiling t') := by
    intro x y hx hy hxy
    have h1 := h x y hx hy hxy
    unfold switch at h1
    by_cases hl : L x = L y
    · rw [if_pos hl, if_pos hl] at h1
      exact h1
    · rw [if_neg hl, if_neg hl] at h1
      rw [tiling_isTournament N t x y hx hy hxy, tiling_isTournament N t' x y hx hy hxy] at h1
      exact (not_eq_not _ _).1 h1
  have := hT a b (by omega) hb (by omega)
  rw [tiling_lt t a b hab, tiling_lt t' a b hab] at this
  exact this

/-- `x xor y` is preserved along the base path when both loop sets give the same
path orientations. -/
theorem loops_relation (N : Nat) (L L' : Nat → Bool)
    (hpath : ∀ i, i + 1 < N → (L i = L (i + 1) ↔ L' i = L' (i + 1))) :
    ∀ a, a < N → (L a = L' a ↔ L 0 = L' 0)
  | 0, _ => Iff.rfl
  | a + 1, ha => by
    have ih := loops_relation N L L' hpath a (by omega)
    have hp := hpath a ha
    rw [← ih]
    cases h1 : L a <;> cases h2 : L (a + 1) <;> cases h3 : L' a <;> cases h4 : L' (a + 1) <;>
      rw [h1, h2, h3, h4] at hp <;> simp_all

/-- **A1, fibres.** Two selfies give the same labeled tournament iff they
have the same tiles and the same loop set or complementary loop sets. -/
theorem selfie_fibre (N : Nat) (t t' : Nat → Nat → Bool) (L L' : Nat → Bool) :
    SameOn N (switch L (tiling t)) (switch L' (tiling t')) ↔
      SameTiles N t t' ∧ (SameLoops N L L' ∨ SameLoops N L (compl L')) := by
  constructor
  · intro h
    have hpath : ∀ i, i + 1 < N → (L i = L (i + 1) ↔ L' i = L' (i + 1)) := by
      intro i hi
      rw [← path_arc L t i, ← path_arc L' t' i, h (i + 1) i hi (by omega) (by omega)]
    have hrel := loops_relation N L L' hpath
    by_cases h0 : L 0 = L' 0
    · have hL : SameLoops N L L' := fun a ha => (hrel a ha).2 h0
      refine ⟨tiles_of_switch_eq N L t t' ?_, Or.inl hL⟩
      intro a b ha hb hab
      rw [h a b ha hb hab]
      exact (switch_loops_congr N L L' (tiling t') hL a b ha hb hab).symm
    · have hL : SameLoops N L (compl L') := by
        intro a ha
        have : ¬ L a = L' a := fun h' => h0 ((hrel a ha).1 h')
        unfold compl
        cases h1 : L a <;> cases h2 : L' a <;> rw [h1, h2] at this <;> simp_all
      refine ⟨tiles_of_switch_eq N L t t' ?_, Or.inr hL⟩
      intro a b ha hb hab
      rw [h a b ha hb hab, ← switch_compl L' (tiling t') a b]
      exact (switch_loops_congr N L (compl L') (tiling t') hL a b ha hb hab).symm
  · rintro ⟨ht, hL | hL⟩
    · intro a b ha hb hab
      rw [switch_loops_congr N L L' (tiling t) hL a b ha hb hab]
      exact switch_congr N L' (tiling t) (tiling t') (tiling_congr N t t' ht) a b ha hb hab
    · intro a b ha hb hab
      rw [switch_loops_congr N L (compl L') (tiling t) hL a b ha hb hab, switch_compl]
      exact switch_congr N L' (tiling t) (tiling t') (tiling_congr N t t' ht) a b ha hb hab

/-- The two points of a fibre are distinct: a loop set is never its own complement
on a nonempty vertex set. -/
theorem loops_ne_compl (N : Nat) (hN : 1 ≤ N) (L : Nat → Bool) : ¬ SameLoops N L (compl L) := by
  intro h
  have := h 0 (by omega)
  unfold compl at this
  cases h1 : L 0 <;> rw [h1] at this <;> exact Bool.noConfusion this

/-- The gauge: loop bits that make every base-path arc point down. -/
def gauge (T : Nat → Nat → Bool) : Nat → Bool
  | 0 => false
  | i + 1 => if T (i + 1) i = true then gauge T i else !gauge T i

theorem gauge_path (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) (i : Nat)
    (hi : i + 1 < N) : switch (gauge T) T (i + 1) i = true := by
  unfold switch
  show (if (if T (i + 1) i = true then gauge T i else !gauge T i) = gauge T i
    then T (i + 1) i else T i (i + 1)) = true
  by_cases h : T (i + 1) i = true
  · rw [if_pos h, if_pos rfl]
    exact h
  · rw [if_neg h]
    have hne : (!gauge T i) ≠ gauge T i := by cases gauge T i <;> decide
    rw [if_neg hne, hT (i + 1) i hi (by omega) (by omega)]
    cases h2 : T (i + 1) i
    · rfl
    · exact absurd h2 h

/-- **A1, surjectivity.** Every labeled tournament is `switch L (tiling t)` for some
selfie `(t, L)`. -/
theorem selfie_onto (N : Nat) (T : Nat → Nat → Bool) (hT : IsTournament N T) :
    ∃ (t : Nat → Nat → Bool) (L : Nat → Bool), SameOn N (switch L (tiling t)) T := by
  let L := gauge T
  let T' := switch L T
  have hT' : IsTournament N T' := switch_isTournament N L T hT
  have hdown : ∀ i, i + 1 < N → T' (i + 1) i = true := fun i hi => gauge_path N T hT i hi
  refine ⟨T', L, ?_⟩
  have htil : SameOn N (tiling T') T' := by
    intro a b ha hb hab
    rcases four_cases a b hab with h | h | h | h
    · subst h
      rw [tiling_down, hdown b ha]
    · subst h
      rw [tiling_up, hT' (a + 1) a hb ha (by omega), hdown a hb]
      rfl
    · rw [tiling_lt T' a b h]
    · rw [tiling_gt T' a b h, hT' a b ha hb hab, Bool.not_not]
  intro a b ha hb hab
  rw [switch_congr N L _ _ htil a b ha hb hab]
  exact switch_switch L T a b

/-- Loops versus path arcs: the number of reversed base-path arcs has the parity of
`[0 ∈ L] + [N-1 ∈ L]` (the discrete derivative telescopes). -/
def reversedCount (L : Nat → Bool) : Nat → Nat
  | 0 => 0
  | i + 1 => reversedCount L i + (if L i = L (i + 1) then 0 else 1)

theorem reversedCount_parity (L : Nat → Bool) :
    ∀ n, reversedCount L n % 2 = (if L 0 = L n then 0 else 1)
  | 0 => by simp [reversedCount]
  | n + 1 => by
    have ih := reversedCount_parity L n
    unfold reversedCount
    cases h0 : L 0 <;> cases h1 : L n <;> cases h2 : L (n + 1) <;> rw [h0, h1] at ih <;>
      simp only [h1, h2, h0] at * <;> simp at * <;> omega

/-- The parameter count `C(N-1, 2) + N = C(N, 2) + 1` (with `pairs n = C(n, 2)`
counted as `0 + 1 + ⋯ + (n-1)`): the one extra bit is the global complement. -/
def pairs : Nat → Nat
  | 0 => 0
  | n + 1 => pairs n + n

theorem selfie_parameter_count (N : Nat) (hN : 1 ≤ N) : pairs (N - 1) + N = pairs N + 1 := by
  obtain ⟨m, rfl⟩ : ∃ m, N = m + 1 := ⟨N - 1, by omega⟩
  show pairs m + (m + 1) = (pairs m + m) + 1
  omega

theorem pairs_eq : ∀ (n : Nat), 2 * pairs n = n * (n - 1)
  | 0 => rfl
  | n + 1 => by
    have ih := pairs_eq n
    show 2 * (pairs n + n) = (n + 1) * (n + 1 - 1)
    rw [Nat.mul_add, ih, Nat.add_sub_cancel]
    cases n with
    | zero => rfl
    | succ k =>
      rw [Nat.add_sub_cancel]
      simp only [Nat.add_mul, Nat.mul_add, Nat.one_mul, Nat.mul_one]
      omega

end ProcgenSelfieEdim
