import ProcgenSelfieEdim.ArcParity
import ProcgenSelfieEdim.TournamentCode

set_option autoImplicit false

/-!
# THM-4524 (b): finite facts about arc-HP counts

* `qr7_counts`: `H(QR_7) = 189` and every arc of `QR_7` lies on exactly 54 HPs.
* `qr7del_counts`: for every vertex `z` of `QR_7`, the tournament `QR_7 - z` has
  `H = 45`; an arc `u → v` lies on 23 HPs if `u + v ≡ 2z (mod 7)` (the three arcs
  between `z`-antipodal pairs) and on 13 HPs otherwise. In particular every arc is
  on an odd number of HPs (`qr7del_allOdd`).
* `no_allOdd_five`: every tournament on 5 vertices has an arc on an even number of
  HPs (exhaustive over the `2^10` labeled tournaments, in 16 chunks).
-/

namespace ProcgenSelfieEdim

/-! ## Counts depend only on the labeled tournament -/

theorem isPath_congr (n : Nat) (T T' : Nat → Nat → Bool) (h : SameOn n T T') :
    ∀ (p : List Nat), p.Nodup → (∀ x ∈ p, x < n) → (IsPath T p ↔ IsPath T' p)
  | [], _, _ => Iff.rfl
  | [_], _, _ => Iff.rfl
  | a :: b :: rest, hnd, hlt => by
    rw [isPath_cons_cons, isPath_cons_cons]
    have hab : a ≠ b := fun hab => by
      subst hab
      exact (List.nodup_cons.1 hnd).1 List.mem_cons_self
    rw [h a b (hlt a List.mem_cons_self) (hlt b (List.mem_cons_of_mem a List.mem_cons_self)) hab,
      isPath_congr n T T' h (b :: rest) (List.nodup_cons.1 hnd).2
        (fun x hx => hlt x (List.mem_cons_of_mem a hx))]

theorem isHamPath_congr (n : Nat) (T T' : Nat → Nat → Bool) (h : SameOn n T T') (p : List Nat) :
    IsHamPath n T p ↔ IsHamPath n T' p := by
  unfold IsHamPath
  constructor
  · rintro ⟨h1, h2, h3, h4⟩
    exact ⟨h1, h2, h3, (isPath_congr n T T' h p h2 h3).1 h4⟩
  · rintro ⟨h1, h2, h3, h4⟩
    exact ⟨h1, h2, h3, (isPath_congr n T T' h p h2 h3).2 h4⟩

theorem hpCount_congr (n : Nat) (T T' : Nat → Nat → Bool) (h : SameOn n T T') :
    hpCount n T = hpCount n T' :=
  hpCount_unique n T' (hamPathsR n T) (nodup_hamPathsR n T)
    (fun p => (mem_hamPathsR n T p).trans (isHamPath_congr n T T' h p))

theorem arcCount_congr (n : Nat) (T T' : Nat → Nat → Bool) (h : SameOn n T T') (u v : Nat) :
    arcCount n T u v = arcCount n T' u v :=
  arcCount_unique n T' u v _ (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T))
    (fun p => by rw [List.mem_filter, mem_hamPathsR, isHamPath_congr n T T' h p])

/-! ## Kernel-friendly counting -/

/-- `usesArc u v (prev :: l)`, by `List.rec`. -/
noncomputable def usesArcFrom (u v : Nat) (l : List Nat) : Nat → Bool :=
  List.rec (motive := fun _ => Nat → Bool) (fun _ => false)
    (fun x _ ih prev => (Nat.beq prev u && Nat.beq x v) || ih x) l

/-- `usesArc`, by `List.rec`. -/
noncomputable def usesArcR (u v : Nat) (p : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) p false (fun a rest => usesArcFrom u v rest a)

theorem usesArcFrom_eq (u v : Nat) :
    ∀ (l : List Nat) (a : Nat), usesArcFrom u v l a = usesArc u v (a :: l)
  | [], _ => rfl
  | x :: t, a => by
    show ((Nat.beq a u && Nat.beq x v) || usesArcFrom u v t x) = usesArc u v (a :: x :: t)
    rw [usesArcFrom_eq u v t x]
    rfl

theorem usesArcR_eq (u v : Nat) : ∀ (p : List Nat), usesArcR u v p = usesArc u v p
  | [] => rfl
  | a :: rest => usesArcFrom_eq u v rest a

/-- `c(u → v)`, evaluated with recursor-based functions. -/
noncomputable def arcCountR (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) : Nat :=
  countR (usesArcR u v) (hamPathsR n T)

theorem arcCountR_eq (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) :
    arcCountR n T u v = arcCount n T u v := by
  unfold arcCountR arcCount
  rw [countR_eq, List.filter_congr (fun p _ => usesArcR_eq u v p)]

/-- Length, by `List.rec`. -/
noncomputable def lenR {α : Type} (l : List α) : Nat :=
  List.rec (motive := fun _ => Nat) 0 (fun _ _ ih => ih + 1) l

theorem lenR_eq {α : Type} : ∀ (l : List α), lenR l = l.length
  | [] => rfl
  | _ :: t => by
    show lenR t + 1 = _
    rw [lenR_eq t, List.length_cons]

/-- Check `H(T) = h` and `ok u v c(u → v)` for every arc `u → v` with `u, v < n`. -/
noncomputable def checkArcs (n : Nat) (T : Nat → Nat → Bool) (h : Nat)
    (ok : Nat → Nat → Nat → Bool) : Bool :=
  Nat.beq (lenR (hamPathsR n T)) h &&
    allR (fun u => allR (fun v => !(T u v) || ok u v (arcCountR n T u v)) (downR n)) (downR n)

theorem checkArcs_sound (n : Nat) (T : Nat → Nat → Bool) (h : Nat) (ok : Nat → Nat → Nat → Bool)
    (hc : checkArcs n T h ok = true) :
    hpCount n T = h ∧ ∀ u v, u < n → v < n → T u v = true → ok u v (arcCount n T u v) = true := by
  unfold checkArcs at hc
  rw [Bool.and_eq_true, Nat.beq_eq] at hc
  refine ⟨by unfold hpCount; rw [← lenR_eq]; exact hc.1, ?_⟩
  intro u v hu hv hT
  have h1 := allR_eq_true _ _ hc.2 u ((mem_downR n u).2 hu)
  have h2 := allR_eq_true _ _ h1 v ((mem_downR n v).2 hv)
  rw [hT, Bool.not_true, Bool.false_or, arcCountR_eq] at h2
  exact h2

/-! ## The Paley tournament `QR_7` -/

/-- `QR_7` on `Z_7 = {0, …, 6}`: `a → b` iff `b - a ∈ {1, 2, 4} (mod 7)`. -/
def qr7 (a b : Nat) : Bool :=
  Nat.beq ((b + 7 - a) % 7) 1 || Nat.beq ((b + 7 - a) % 7) 2 || Nat.beq ((b + 7 - a) % 7) 4

theorem qr7_isTournament_aux :
    ∀ a, a < 7 → ∀ b, b < 7 → a ≠ b → qr7 b a = !qr7 a b := by
  decide

theorem qr7_isTournament : IsTournament 7 qr7 :=
  fun a b ha hb hab => qr7_isTournament_aux a ha b hb hab

theorem qr7_check : checkArcs 7 qr7 189 (fun _ _ c => Nat.beq c 54) = true := by
  decide +kernel

/-- **`QR_7`.** `H(QR_7) = 189`, and every arc lies on exactly 54 Hamiltonian paths. -/
theorem qr7_counts : hpCount 7 qr7 = 189 ∧
    ∀ u v, u < 7 → v < 7 → qr7 u v = true → arcCount 7 qr7 u v = 54 := by
  obtain ⟨h1, h2⟩ := checkArcs_sound 7 qr7 189 _ qr7_check
  exact ⟨h1, fun u v hu hv huv => Nat.eq_of_beq_eq_true (h2 u v hu hv huv)⟩

/-- Relabel `{0, …, 5}` as the vertices of `QR_7` other than `z`, in increasing order. -/
def lift (z a : Nat) : Nat := if a < z then a else a + 1

/-- `QR_7 - z` on `{0, …, 5}`. -/
def qr7del (z a b : Nat) : Bool := qr7 (lift z a) (lift z b)

theorem lift_lt (z a : Nat) (ha : a < 6) : lift z a < 7 := by
  unfold lift
  split <;> omega

theorem lift_ne (z a : Nat) : lift z a ≠ z := by
  unfold lift
  split <;> omega

theorem lift_inj (z a b : Nat) (h : lift z a = lift z b) : a = b := by
  unfold lift at h
  split at h <;> split at h <;> omega

theorem qr7del_isTournament (z : Nat) : IsTournament 6 (qr7del z) := by
  intro a b ha hb hab
  exact qr7_isTournament (lift z a) (lift z b) (lift_lt z a ha) (lift_lt z b hb)
    (fun h => hab (lift_inj z a b h))

/-- `u → v` joins a `z`-antipodal pair: `u + v ≡ 2 z (mod 7)` in `QR_7` labels. -/
def antipodalTo (z u v : Nat) : Bool := Nat.beq ((lift z u + lift z v) % 7) (2 * z % 7)

/-- The arc-count rule for `QR_7 - z`. -/
def qr7delRule (z u v c : Nat) : Bool :=
  cond (antipodalTo z u v) (Nat.beq c 23) (Nat.beq c 13)

theorem qr7del_check : allR (fun z => checkArcs 6 (qr7del z) 45 (qr7delRule z)) (downR 7) = true := by
  decide +kernel

/-- **`QR_7` minus a vertex.** For every vertex `z` of `QR_7`: `H(QR_7 - z) = 45`, and an
arc `u → v` of `QR_7 - z` lies on 23 Hamiltonian paths if it joins a `z`-antipodal
pair (`u + v ≡ 2z mod 7` in `QR_7` labels) and on 13 otherwise. -/
theorem qr7del_counts (z : Nat) (hz : z < 7) : hpCount 6 (qr7del z) = 45 ∧
    ∀ u v, u < 6 → v < 6 → qr7del z u v = true →
      (antipodalTo z u v = true ∧ arcCount 6 (qr7del z) u v = 23) ∨
      (antipodalTo z u v = false ∧ arcCount 6 (qr7del z) u v = 13) := by
  have hc := allR_eq_true _ _ qr7del_check z ((mem_downR 7 z).2 hz)
  obtain ⟨h1, h2⟩ := checkArcs_sound 6 (qr7del z) 45 _ hc
  refine ⟨h1, fun u v hu hv huv => ?_⟩
  have h3 := h2 u v hu hv huv
  unfold qr7delRule at h3
  cases ha : antipodalTo z u v
  · rw [ha, cond_false] at h3
    exact Or.inr ⟨rfl, Nat.eq_of_beq_eq_true h3⟩
  · rw [ha, cond_true] at h3
    exact Or.inl ⟨rfl, Nat.eq_of_beq_eq_true h3⟩

/-- Exactly three arcs of `QR_7 - z` are `z`-antipodal (so 3 arcs have count 23 and the
other 12 have count 13). -/
theorem qr7del_antipodal_arcs : allR (fun z => Nat.beq (rsum 6 (fun u => rsum 6 (fun v =>
    ind (qr7del z u v && antipodalTo z u v)))) 3) (downR 7) = true := by
  decide +kernel

/-- **Every arc of `QR_7 - z` lies on an odd number of Hamiltonian paths.** -/
theorem qr7del_allOdd (z : Nat) (hz : z < 7) :
    ∀ u v, u < 6 → v < 6 → qr7del z u v = true → arcCount 6 (qr7del z) u v % 2 = 1 := by
  intro u v hu hv huv
  rcases (qr7del_counts z hz).2 u v hu hv huv with ⟨_, h⟩ | ⟨_, h⟩ <;> rw [h]

/-! ## No tournament on five vertices is all-odd -/

/-- The coded tournament `tourC c` on `n` vertices has an arc on an even number of HPs. -/
noncomputable def hasEvenArc (n c : Nat) : Bool :=
  anyR (fun u => anyR (fun v =>
    tourC c u v && Nat.beq (arcCountR n (tourC c) u v % 2) 0) (downR n)) (downR n)

theorem hasEvenArc_sound (n c : Nat) (h : hasEvenArc n c = true) :
    ∃ u v, u < n ∧ v < n ∧ tourC c u v = true ∧ arcCount n (tourC c) u v % 2 = 0 := by
  obtain ⟨u, hu, h1⟩ := anyR_eq_true _ _ h
  obtain ⟨v, hv, h2⟩ := anyR_eq_true _ _ h1
  rw [Bool.and_eq_true, Nat.beq_eq, arcCountR_eq] at h2
  exact ⟨u, v, (mem_downR n u).1 hu, (mem_downR n v).1 hv, h2.1, h2.2⟩

/-- Check `p (lo + i)` for `i < len`. -/
noncomputable def checkFrom (p : Nat → Bool) (lo len : Nat) : Bool :=
  allR (fun i => p (lo + i)) (downR len)

theorem checkFrom_sound (p : Nat → Bool) (lo len : Nat) (h : checkFrom p lo len = true) :
    ∀ c, lo ≤ c → c < lo + len → p c = true := by
  intro c h1 h2
  have := allR_eq_true _ _ h (c - lo) ((mem_downR len (c - lo)).2 (by omega))
  rw [show lo + (c - lo) = c by omega] at this
  exact this

/-- From a check of every code to a statement about every tournament. -/
theorem exists_even_arc_of_codes (n : Nat) (hcodes : ∀ c, c < 2 ^ pairs n → hasEvenArc n c = true)
    (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    ∃ u v, u < n ∧ v < n ∧ u ≠ v ∧ T u v = true ∧ arcCount n T u v % 2 = 0 := by
  obtain ⟨hlt, hsame⟩ := tourC_code n T hT
  obtain ⟨u, v, hu, hv, h1, h2⟩ := hasEvenArc_sound n _ (hcodes _ hlt)
  have huv : u ≠ v := fun h => by
    subst h
    rw [tourC_self] at h1
    exact Bool.noConfusion h1
  refine ⟨u, v, hu, hv, huv, ?_, ?_⟩
  · rw [hsame u v hu hv huv]
    exact h1
  · rw [arcCount_congr n T _ hsame]
    exact h2

theorem three_codes : checkFrom (hasEvenArc 3) 0 8 = true := by decide +kernel

theorem four_codes : checkFrom (hasEvenArc 4) 0 64 = true := by decide +kernel

theorem five_chunk_0 : checkFrom (hasEvenArc 5) 0 128 = true := by decide +kernel

theorem five_chunk_1 : checkFrom (hasEvenArc 5) 128 128 = true := by decide +kernel

theorem five_chunk_2 : checkFrom (hasEvenArc 5) 256 128 = true := by decide +kernel

theorem five_chunk_3 : checkFrom (hasEvenArc 5) 384 128 = true := by decide +kernel

theorem five_chunk_4 : checkFrom (hasEvenArc 5) 512 128 = true := by decide +kernel

theorem five_chunk_5 : checkFrom (hasEvenArc 5) 640 128 = true := by decide +kernel

theorem five_chunk_6 : checkFrom (hasEvenArc 5) 768 128 = true := by decide +kernel

theorem five_chunk_7 : checkFrom (hasEvenArc 5) 896 128 = true := by decide +kernel

theorem five_codes (c : Nat) (hc : c < 2 ^ pairs 5) : hasEvenArc 5 c = true := by
  have hc' : c < 1024 := by
    have : 2 ^ pairs 5 = 1024 := by decide
    omega
  by_cases h0 : c < 128
  · exact checkFrom_sound _ 0 128 five_chunk_0 c (by omega) (by omega)
  by_cases h1 : c < 256
  · exact checkFrom_sound _ 128 128 five_chunk_1 c (by omega) (by omega)
  by_cases h2 : c < 384
  · exact checkFrom_sound _ 256 128 five_chunk_2 c (by omega) (by omega)
  by_cases h3 : c < 512
  · exact checkFrom_sound _ 384 128 five_chunk_3 c (by omega) (by omega)
  by_cases h4 : c < 640
  · exact checkFrom_sound _ 512 128 five_chunk_4 c (by omega) (by omega)
  by_cases h5 : c < 768
  · exact checkFrom_sound _ 640 128 five_chunk_5 c (by omega) (by omega)
  by_cases h6 : c < 896
  · exact checkFrom_sound _ 768 128 five_chunk_6 c (by omega) (by omega)
  · exact checkFrom_sound _ 896 128 five_chunk_7 c (by omega) (by omega)

/-- **No all-odd tournament on 3, 4 or 5 vertices**: every tournament on `{0, …, n-1}`
with `3 ≤ n ≤ 5` has an arc lying on an even number of Hamiltonian paths
(exhaustive over the `2^C(n,2)` labeled tournaments). -/
theorem no_allOdd_three_to_five (n : Nat) (hn : 3 ≤ n ∧ n ≤ 5) (T : Nat → Nat → Bool)
    (hT : IsTournament n T) :
    ∃ u v, u < n ∧ v < n ∧ u ≠ v ∧ T u v = true ∧ arcCount n T u v % 2 = 0 := by
  have e3 : 2 ^ pairs 3 = 8 := by decide
  have e4 : 2 ^ pairs 4 = 64 := by decide
  have c3 : ∀ c, c < 2 ^ pairs 3 → hasEvenArc 3 c = true := fun c hc =>
    checkFrom_sound _ 0 8 three_codes c (Nat.zero_le c) (by omega)
  have c4 : ∀ c, c < 2 ^ pairs 4 → hasEvenArc 4 c = true := fun c hc =>
    checkFrom_sound _ 0 64 four_codes c (Nat.zero_le c) (by omega)
  obtain ⟨h3, h5⟩ := hn
  by_cases hn3 : n = 3
  · subst hn3
    exact exists_even_arc_of_codes 3 c3 T hT
  by_cases hn4 : n = 4
  · subst hn4
    exact exists_even_arc_of_codes 4 c4 T hT
  have hn5 : n = 5 := by omega
  subst hn5
  exact exists_even_arc_of_codes 5 five_codes T hT

end ProcgenSelfieEdim
