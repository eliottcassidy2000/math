import ProcgenSelfieEdim.SelfieFinite

set_option autoImplicit false

/-!
# Shaved tournaments: `H_4` and `H_5` (opus S15, THM-4526)

`H_n` is the oriented graph "Hamiltonian path `v_1 → ⋯ → v_n` plus the arc
`v_1 → v_n`". A copy of `H_n` in `T` is therefore a Hamiltonian path of `T` whose first
vertex beats its last (`ContainsH`).

* `every_four_contains_H4`: every tournament on 4 vertices contains `H_4`, i.e. vertices
  `A, B, C, D` with `A → B`, `B → C`, `C → D`, `A → D`.
* `avoids_H5_iff`: a tournament on 5 vertices contains no copy of `H_5` iff it is
  isomorphic to `C3[1, C3, 1]` (`a` beats a 3-cycle, the 3-cycle beats `b`, `b → a`).

Both are exhaustive over labeled tournaments (codes, `TournamentCode`); the `H_5`
check is split into 8 chunks of 128 codes.
-/

namespace ProcgenSelfieEdim

/-- `T` contains a copy of `H_n`: a Hamiltonian path `a → ⋯ → b` with the arc `a → b`. -/
def ContainsH (n : Nat) (T : Nat → Nat → Bool) : Prop :=
  ∃ (a b : Nat) (mid : List Nat), IsHamPath n T (a :: (mid ++ [b])) ∧ T a b = true

/-- `T` is isomorphic to `S` on `{0, …, n-1}`: an injective self-map `π` of
`{0, …, n-1}` with `T i j = S (π i) (π j)` for all distinct `i, j`. -/
def IsoTo (n : Nat) (T S : Nat → Nat → Bool) : Prop :=
  ∃ π : Nat → Nat, (∀ i, i < n → π i < n) ∧ (∀ i j, i < n → j < n → π i = π j → i = j) ∧
    ∀ i j, i < n → j < n → i ≠ j → S (π i) (π j) = T i j

/-- The labeled `C3[1, C3, 1]`: `0` beats the 3-cycle `1 → 2 → 3 → 1`, which beats `4`,
and `4 → 0`. -/
def c3c3 (x y : Nat) : Bool :=
  [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4), (4, 0)].contains (x, y)

theorem c3c3_isTournament_aux : ∀ a, a < 5 → ∀ b, b < 5 → a ≠ b → c3c3 b a = !c3c3 a b := by
  decide

theorem c3c3_isTournament : IsTournament 5 c3c3 :=
  fun a b ha hb hab => c3c3_isTournament_aux a ha b hb hab

/-! ## Last elements and the `H_n` test -/

/-- The last element of `a :: l`. -/
noncomputable def lastOr (a : Nat) (l : List Nat) : Nat :=
  List.rec (motive := fun _ => Nat → Nat) (fun x => x) (fun y _ ih _ => ih y) l a

theorem lastOr_nil (a : Nat) : lastOr a [] = a := rfl

theorem lastOr_cons (a y : Nat) (l : List Nat) : lastOr a (y :: l) = lastOr y l := rfl

theorem eq_append_lastOr : ∀ (l : List Nat) (a : Nat), l ≠ [] → ∃ mid, l = mid ++ [lastOr a l]
  | [], _, h => absurd rfl h
  | [y], a, _ => ⟨[], by rw [lastOr_cons, lastOr_nil]; rfl⟩
  | y :: z :: t, a, _ => by
    obtain ⟨mid, hmid⟩ := eq_append_lastOr (z :: t) y (List.cons_ne_nil _ _)
    refine ⟨y :: mid, ?_⟩
    rw [lastOr_cons]
    exact congrArg (List.cons y) hmid

theorem lastOr_append (a b : Nat) : ∀ (mid : List Nat), lastOr a (mid ++ [b]) = b
  | [] => rfl
  | y :: t => by
    rw [List.cons_append, lastOr_cons]
    exact lastOr_append y b t

/-- The first vertex of `p` beats its last vertex. -/
noncomputable def headBeatsLast (T : Nat → Nat → Bool) (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a rest => T a (lastOr a rest))

/-- Some Hamiltonian path of `tourC c` on `n` vertices is a copy of `H_n`. -/
noncomputable def hasH (n c : Nat) : Bool := anyR (headBeatsLast (tourC c)) (hamPathsR n (tourC c))

theorem containsH_of_mem (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (p : List Nat)
    (hp : p ∈ hamPathsR n T) (h : headBeatsLast T p = true) : ContainsH n T := by
  have hham := (mem_hamPathsR n T p).1 hp
  cases p with
  | nil =>
    exfalso
    have := hham.1
    rw [List.length_nil] at this
    omega
  | cons a rest =>
    have hrest : rest ≠ [] := by
      intro hr
      subst hr
      have := hham.1
      rw [List.length_singleton] at this
      omega
    obtain ⟨mid, hmid⟩ := eq_append_lastOr rest a hrest
    refine ⟨a, lastOr a rest, mid, ?_, h⟩
    rw [← hmid]
    exact hham

theorem hasH_sound (n : Nat) (hn : 2 ≤ n) (c : Nat) (h : hasH n c = true) : ContainsH n (tourC c) := by
  obtain ⟨p, hp, hpH⟩ := anyR_eq_true _ _ h
  exact containsH_of_mem n hn _ p hp hpH

theorem anyR_of_mem {α : Type} (p : α → Bool) : ∀ (l : List α) (x : α), x ∈ l → p x = true → anyR p l = true
  | [], _, hx, _ => absurd hx List.not_mem_nil
  | y :: t, x, hx, hpx => by
    show (p y || anyR p t) = true
    rcases List.mem_cons.1 hx with rfl | hx'
    · rw [hpx, Bool.true_or]
    · rw [anyR_of_mem p t x hx' hpx, Bool.or_true]

/-- Conversely, a copy of `H_n` is found by the test. -/
theorem anyR_headBeatsLast_of_containsH (n : Nat) (T : Nat → Nat → Bool) (h : ContainsH n T) :
    anyR (headBeatsLast T) (hamPathsR n T) = true := by
  obtain ⟨a, b, mid, hham, hab⟩ := h
  apply anyR_of_mem _ _ (a :: (mid ++ [b])) ((mem_hamPathsR n T _).2 hham)
  show T a (lastOr a (mid ++ [b])) = true
  rw [lastOr_append]
  exact hab

theorem containsH_congr (n : Nat) (T T' : Nat → Nat → Bool) (hs : SameOn n T T') (h : ContainsH n T) :
    ContainsH n T' := by
  obtain ⟨a, b, mid, hham, hab⟩ := h
  refine ⟨a, b, mid, (isHamPath_congr n T T' hs _).1 hham, ?_⟩
  have hb : b ∈ a :: (mid ++ [b]) :=
    List.mem_cons_of_mem a (List.mem_append_right _ List.mem_cons_self)
  have hne : a ≠ b := by
    intro h'
    subst h'
    exact (List.nodup_cons.1 hham.2.1).1 (List.mem_append_right _ List.mem_cons_self)
  rw [← hs a b (hham.2.2.1 a List.mem_cons_self) (hham.2.2.1 b hb) hne]
  exact hab

/-! ## `H_4` -/

theorem four_H_codes : checkFrom (hasH 4) 0 64 = true := by
  decide +kernel

/-- **Every tournament on 4 vertices contains `H_4`**: there are distinct `A, B, C, D`
with `A → B`, `B → C`, `C → D` and `A → D`. -/
theorem every_four_contains_H4 (T : Nat → Nat → Bool) (hT : IsTournament 4 T) :
    ∃ A B C D, A < 4 ∧ B < 4 ∧ C < 4 ∧ D < 4 ∧ [A, B, C, D].Nodup ∧
      T A B = true ∧ T B C = true ∧ T C D = true ∧ T A D = true := by
  obtain ⟨hlt, hsame⟩ := tourC_code 4 T hT
  have h64 : 2 ^ pairs 4 = 64 := by decide
  have hc := checkFrom_sound _ 0 64 four_H_codes (bitsToNat (codeList T 4)) (Nat.zero_le _)
    (by omega)
  have hH := containsH_congr 4 _ T (fun a b ha hb hab => (hsame a b ha hb hab).symm)
    (hasH_sound 4 (by decide) _ hc)
  obtain ⟨a, d, mid, ⟨hl, hnd, hrange, hpath⟩, had⟩ := hH
  rcases mid with _ | ⟨b, _ | ⟨c, _ | ⟨e, t⟩⟩⟩
  · exfalso
    simp only [List.nil_append, List.length_cons, List.length_nil] at hl
    omega
  · exfalso
    simp only [List.cons_append, List.nil_append, List.length_cons, List.length_nil] at hl
    omega
  · have hm : ∀ x, x ∈ [a, b, c, d] → x < 4 := hrange
    exact ⟨a, b, c, d, hm a List.mem_cons_self,
      hm b (List.mem_cons_of_mem _ List.mem_cons_self),
      hm c (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ List.mem_cons_self)),
      hm d (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ List.mem_cons_self))),
      hnd, hpath.1, hpath.2.1, hpath.2.2.1, had⟩
  · exfalso
    simp only [List.cons_append, List.length_cons, List.length_append] at hl
    omega

/-! ## `H_5` and `C3[1, C3, 1]` -/

/-- Indexing, by `List.rec` (default `0`). -/
noncomputable def nthR (l : List Nat) : Nat → Nat :=
  List.rec (motive := fun _ => Nat → Nat) (fun _ => 0)
    (fun x _ ih i => Nat.casesOn (motive := fun _ => Nat) i x (fun j => ih j)) l

theorem nthR_mem : ∀ (l : List Nat) (i : Nat), i < l.length → nthR l i ∈ l
  | [], _, h => absurd h (Nat.not_lt_zero _)
  | x :: _, 0, _ => List.mem_cons_self
  | x :: t, i + 1, h => by
    show nthR t i ∈ x :: t
    exact List.mem_cons_of_mem x (nthR_mem t i (by simp only [List.length_cons] at h; omega))

theorem nthR_inj : ∀ (l : List Nat), l.Nodup → ∀ (i j : Nat), i < l.length → j < l.length →
    nthR l i = nthR l j → i = j
  | [], _, _, _, hi, _, _ => absurd hi (Nat.not_lt_zero _)
  | x :: t, hnd, i, j, hi, hj, h => by
    rw [List.nodup_cons] at hnd
    cases i with
    | zero =>
      cases j with
      | zero => rfl
      | succ j =>
        have h' : x = nthR t j := h
        exact absurd (h' ▸ nthR_mem t j (by simp only [List.length_cons] at hj; omega)) hnd.1
    | succ i =>
      cases j with
      | zero =>
        have h' : nthR t i = x := h
        exact absurd (h' ▸ nthR_mem t i (by simp only [List.length_cons] at hi; omega)) hnd.1
      | succ j =>
        have h' : nthR t i = nthR t j := h
        have := nthR_inj t hnd.2 i j (by simp only [List.length_cons] at hi; omega)
          (by simp only [List.length_cons] at hj; omega) h'
        rw [this]

/-- The complete relation: its Hamiltonian paths are all permutations. -/
def complete (_ _ : Nat) : Bool := true

theorem isPath_complete : ∀ (p : List Nat), IsPath complete p
  | [] => trivial
  | [_] => trivial
  | _ :: b :: rest => ⟨rfl, isPath_complete (b :: rest)⟩

/-- `π` (a permutation list) maps `tourC c` onto `c3c3`. -/
noncomputable def isoCheck (c : Nat) (π : List Nat) : Bool :=
  allR (fun i => allR (fun j =>
    Nat.beq i j || (c3c3 (nthR π i) (nthR π j) == tourC c i j)) (downR 5)) (downR 5)

/-- Some permutation maps `tourC c` onto `c3c3`. -/
noncomputable def findIso (c : Nat) : Bool := anyR (isoCheck c) (hamPathsR 5 complete)

theorem findIso_sound (c : Nat) (h : findIso c = true) : IsoTo 5 (tourC c) c3c3 := by
  obtain ⟨π, hπ, hcheck⟩ := anyR_eq_true _ _ h
  obtain ⟨hl, hnd, hlt, _⟩ := (mem_hamPathsR 5 complete π).1 hπ
  refine ⟨nthR π, fun i hi => hlt _ (nthR_mem π i (by omega)),
    fun i j hi hj hij => nthR_inj π hnd i j (by omega) (by omega) hij, ?_⟩
  intro i j hi hj hij
  have h1 := allR_eq_true _ _ hcheck i ((mem_downR 5 i).2 hi)
  have h2 := allR_eq_true _ _ h1 j ((mem_downR 5 j).2 hj)
  rw [beq_false_of_ne i j hij, Bool.false_or] at h2
  exact beq_iff_eq.1 h2

theorem five_H_chunk_0 : checkFrom (fun c => hasH 5 c || findIso c) 0 128 = true := by decide +kernel

theorem five_H_chunk_1 : checkFrom (fun c => hasH 5 c || findIso c) 128 128 = true := by decide +kernel

theorem five_H_chunk_2 : checkFrom (fun c => hasH 5 c || findIso c) 256 128 = true := by decide +kernel

theorem five_H_chunk_3 : checkFrom (fun c => hasH 5 c || findIso c) 384 128 = true := by decide +kernel

theorem five_H_chunk_4 : checkFrom (fun c => hasH 5 c || findIso c) 512 128 = true := by decide +kernel

theorem five_H_chunk_5 : checkFrom (fun c => hasH 5 c || findIso c) 640 128 = true := by decide +kernel

theorem five_H_chunk_6 : checkFrom (fun c => hasH 5 c || findIso c) 768 128 = true := by decide +kernel

theorem five_H_chunk_7 : checkFrom (fun c => hasH 5 c || findIso c) 896 128 = true := by decide +kernel

theorem five_H_codes (c : Nat) (hc : c < 2 ^ pairs 5) : (hasH 5 c || findIso c) = true := by
  have hc' : c < 1024 := by
    have : 2 ^ pairs 5 = 1024 := by decide
    omega
  by_cases h0 : c < 128
  · exact checkFrom_sound _ 0 128 five_H_chunk_0 c (by omega) (by omega)
  by_cases h1 : c < 256
  · exact checkFrom_sound _ 128 128 five_H_chunk_1 c (by omega) (by omega)
  by_cases h2 : c < 384
  · exact checkFrom_sound _ 256 128 five_H_chunk_2 c (by omega) (by omega)
  by_cases h3 : c < 512
  · exact checkFrom_sound _ 384 128 five_H_chunk_3 c (by omega) (by omega)
  by_cases h4 : c < 640
  · exact checkFrom_sound _ 512 128 five_H_chunk_4 c (by omega) (by omega)
  by_cases h5 : c < 768
  · exact checkFrom_sound _ 640 128 five_H_chunk_5 c (by omega) (by omega)
  by_cases h6 : c < 896
  · exact checkFrom_sound _ 768 128 five_H_chunk_6 c (by omega) (by omega)
  · exact checkFrom_sound _ 896 128 five_H_chunk_7 c (by omega) (by omega)

theorem c3c3_avoids_check : anyR (headBeatsLast c3c3) (hamPathsR 5 c3c3) = false := by
  decide +kernel

/-- `C3[1, C3, 1]` contains no copy of `H_5`: all 15 of its Hamiltonian paths close up. -/
theorem c3c3_avoids_H5 : ¬ ContainsH 5 c3c3 := by
  intro h
  have := anyR_headBeatsLast_of_containsH 5 c3c3 h
  rw [c3c3_avoids_check] at this
  exact Bool.noConfusion this

/-- An isomorphism carries Hamiltonian paths to Hamiltonian paths. -/
theorem isPath_map (n : Nat) (T S : Nat → Nat → Bool) (π : Nat → Nat)
    (harc : ∀ i j, i < n → j < n → i ≠ j → S (π i) (π j) = T i j) :
    ∀ (p : List Nat), p.Nodup → (∀ x ∈ p, x < n) → IsPath T p → IsPath S (p.map π)
  | [], _, _, _ => trivial
  | [_], _, _, _ => trivial
  | a :: b :: rest, hnd, hlt, hp => by
    have hab : a ≠ b := fun h => by
      subst h
      exact (List.nodup_cons.1 hnd).1 List.mem_cons_self
    refine ⟨?_, isPath_map n T S π harc (b :: rest) (List.nodup_cons.1 hnd).2
      (fun x hx => hlt x (List.mem_cons_of_mem a hx)) hp.2⟩
    rw [harc a b (hlt a List.mem_cons_self) (hlt b (List.mem_cons_of_mem a List.mem_cons_self)) hab]
    exact hp.1

theorem containsH_of_iso (n : Nat) (T S : Nat → Nat → Bool) (hiso : IsoTo n T S)
    (h : ContainsH n T) : ContainsH n S := by
  obtain ⟨π, hπlt, hπinj, harc⟩ := hiso
  obtain ⟨a, b, mid, ⟨hl, hnd, hlt, hp⟩, hab⟩ := h
  refine ⟨π a, π b, mid.map π, ⟨?_, ?_, ?_, ?_⟩, ?_⟩
  · rw [← hl]
    simp only [List.length_cons, List.length_append, List.length_map]
  · have := nodup_map_of_inj_on π (a :: (mid ++ [b])) hnd
      (fun x hx y hy hxy => hπinj x y (hlt x hx) (hlt y hy) hxy)
    simp only [List.map_cons, List.map_append, List.map_cons, List.map_nil] at this
    exact this
  · intro x hx
    have hx' : x ∈ (a :: (mid ++ [b])).map π := by
      simp only [List.map_cons, List.map_append, List.map_nil]
      exact hx
    rw [List.mem_map] at hx'
    obtain ⟨y, hy, rfl⟩ := hx'
    exact hπlt y (hlt y hy)
  · have := isPath_map n T S π harc _ hnd hlt hp
    simp only [List.map_cons, List.map_append, List.map_nil] at this
    exact this
  · have hb : b ∈ a :: (mid ++ [b]) :=
      List.mem_cons_of_mem a (List.mem_append_right _ List.mem_cons_self)
    have hne : a ≠ b := by
      intro h'
      subst h'
      exact (List.nodup_cons.1 hnd).1 (List.mem_append_right _ List.mem_cons_self)
    rw [harc a b (hlt a List.mem_cons_self) (hlt b hb) hne]
    exact hab

/-- **THM-4526, `n = 5`.** A tournament on 5 vertices contains no copy of `H_5` iff it is
isomorphic to `C3[1, C3, 1]`. -/
theorem avoids_H5_iff (T : Nat → Nat → Bool) (hT : IsTournament 5 T) :
    ¬ ContainsH 5 T ↔ IsoTo 5 T c3c3 := by
  constructor
  · intro hno
    obtain ⟨hlt, hsame⟩ := tourC_code 5 T hT
    have hc := five_H_codes _ hlt
    rw [Bool.or_eq_true] at hc
    rcases hc with hH | hI
    · exact absurd (containsH_congr 5 _ T (fun a b ha hb hab => (hsame a b ha hb hab).symm)
        (hasH_sound 5 (by decide) _ hH)) hno
    · obtain ⟨π, h1, h2, h3⟩ := findIso_sound _ hI
      exact ⟨π, h1, h2, fun i j hi hj hij => by rw [h3 i j hi hj hij, hsame i j hi hj hij]⟩
  · intro hiso hH
    exact c3c3_avoids_H5 (containsH_of_iso 5 T c3c3 hiso hH)

end ProcgenSelfieEdim
