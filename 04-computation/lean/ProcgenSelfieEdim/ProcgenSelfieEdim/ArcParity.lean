import ProcgenSelfieEdim.HamPath

set_option autoImplicit false

/-!
# Double counting of arc-HP incidences and the arc-parity constraints (THM-4524 (C))

* `sum_arcCount`: `Σ_{u,v} c(u → v) = (n - 1) · H(T)` for every relation `T`
  (each Hamiltonian path uses exactly `n - 1` arcs).
* `arcTotal`: a tournament on `n` vertices has `C(n, 2) = pairs n` arcs.
* With `H(T)` odd (Rédei's theorem; taken here as a hypothesis, see the note):
  - `oddArcs_parity`: the number of arcs on an odd number of HPs is `≡ n - 1 (mod 2)`;
  - `allOdd_mod_four`: if every arc is on an odd number of HPs, then `n ≡ 1, 2 (mod 4)`;
  - `allEven_odd`: if every arc is on an even number of HPs, then `n` is odd.
* `shave_arc`: deleting the arc `e` removes exactly the `c(e)` paths through it.
-/

namespace ProcgenSelfieEdim

/-- The indicator of a Boolean. -/
def ind (b : Bool) : Nat := if b = true then 1 else 0

theorem beq_false_of_ne (a b : Nat) (h : a ≠ b) : Nat.beq a b = false := by
  cases hb : Nat.beq a b
  · rfl
  · exact absurd (Nat.eq_of_beq_eq_true hb) h

theorem length_filter_eq_sum {α : Type} (f : α → Bool) :
    ∀ (P : List α), (P.filter f).length = (P.map (fun p => ind (f p))).sum
  | [] => rfl
  | x :: t => by
    rw [List.map_cons, List.sum_cons, ← length_filter_eq_sum f t]
    cases hx : f x
    · rw [filter_cons_false f x t hx]
      show _ = (if false = true then 1 else 0) + _
      rw [if_neg (by decide), Nat.zero_add]
    · rw [filter_cons_true f x t hx, List.length_cons]
      show _ = (if true = true then 1 else 0) + _
      rw [if_pos rfl]
      omega

/-- Exchange a finite sum over `0, …, n-1` with a sum over a list. -/
theorem rsum_listSum {α : Type} (n : Nat) (g : Nat → α → Nat) :
    ∀ (P : List α), rsum n (fun u => (P.map (g u)).sum) = (P.map (fun p => rsum n (fun u => g u p))).sum
  | [] => by
    rw [List.map_nil, List.sum_nil]
    exact (rsum_congr n _ (fun _ => 0) (fun u _ => rfl)).trans (rsum_zero n)
  | x :: t => by
    rw [List.map_cons, List.sum_cons, ← rsum_listSum n g t, ← rsum_add]
    apply rsum_congr
    intro u _
    rw [List.map_cons, List.sum_cons]

theorem sum_map_const {α : Type} (c : Nat) : ∀ (P : List α), (P.map (fun _ => c)).sum = c * P.length
  | [] => by rw [List.map_nil, List.sum_nil, List.length_nil, Nat.mul_zero]
  | x :: t => by
    rw [List.map_cons, List.sum_cons, sum_map_const c t, List.length_cons, Nat.mul_succ]
    omega

theorem rsum_ite_out (n : Nat) (c : Prop) [Decidable c] (f : Nat → Nat) :
    rsum n (fun v => if c then f v else 0) = if c then rsum n f else 0 := by
  by_cases hc : c
  · rw [if_pos hc]
    exact rsum_congr n _ _ (fun v _ => if_pos hc)
  · rw [if_neg hc]
    exact (rsum_congr n _ (fun _ => 0) (fun v _ => if_neg hc)).trans (rsum_zero n)

/-- The point mass at `(a, b)` sums to one. -/
theorem rsum_rsum_point (n a b : Nat) (ha : a < n) (hb : b < n) :
    rsum n (fun u => rsum n (fun v => if a = u then (if b = v then 1 else 0) else 0)) = 1 := by
  have hin : rsum n (fun v => if b = v then 1 else 0) = 1 := by
    have := rsum_single n b (fun _ => 1) hb
    exact this
  rw [rsum_congr n _ (fun u => if a = u then rsum n (fun v => if b = v then 1 else 0) else 0)
    (fun u _ => rsum_ite_out n (a = u) (fun v => if b = v then 1 else 0))]
  rw [hin]
  exact rsum_single n a (fun _ => 1) ha

theorem ind_usesArc_cons_cons (u v a b : Nat) (rest : List Nat) (ha : a ∉ b :: rest) :
    ind (usesArc u v (a :: b :: rest)) =
      (if a = u then (if b = v then 1 else 0) else 0) + ind (usesArc u v (b :: rest)) := by
  show ind ((Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)) = _
  by_cases hau : a = u
  · subst hau
    have h2 : usesArc a v (b :: rest) = false := by
      cases h : usesArc a v (b :: rest)
      · rfl
      · exact absurd (usesArc_mem a v _ h).1 ha
    rw [h2, Nat.beq_refl, Bool.true_and, Bool.or_false, if_pos rfl]
    by_cases hbv : b = v
    · subst hbv
      rw [Nat.beq_refl, if_pos rfl]
      rfl
    · rw [beq_false_of_ne b v hbv, if_neg hbv]
      rfl
  · rw [beq_false_of_ne a u hau, Bool.false_and, Bool.false_or, if_neg hau, Nat.zero_add]

/-- A duplicate-free path on `{0, …, n-1}` uses exactly `length - 1` ordered pairs. -/
theorem arcs_of_path (n : Nat) : ∀ (p : List Nat), p.Nodup → (∀ x ∈ p, x < n) →
    rsum n (fun u => rsum n (fun v => ind (usesArc u v p))) = p.length - 1
  | [], _, _ =>
    (rsum_congr n _ (fun _ => 0) (fun u _ =>
      (rsum_congr n _ (fun _ => 0) (fun v _ => rfl)).trans (rsum_zero n))).trans (rsum_zero n)
  | [_], _, _ =>
    (rsum_congr n _ (fun _ => 0) (fun u _ =>
      (rsum_congr n _ (fun _ => 0) (fun v _ => rfl)).trans (rsum_zero n))).trans (rsum_zero n)
  | a :: b :: rest, hnd, hlt => by
    have ha : a ∉ b :: rest := (List.nodup_cons.1 hnd).1
    have ih := arcs_of_path n (b :: rest) (List.nodup_cons.1 hnd).2
      (fun x hx => hlt x (List.mem_cons_of_mem a hx))
    rw [rsum_congr n _ (fun u => rsum n (fun v => if a = u then (if b = v then 1 else 0) else 0) +
        rsum n (fun v => ind (usesArc u v (b :: rest)))) (fun u _ => by
          dsimp only
          rw [← rsum_add]
          exact rsum_congr n _ _ (fun v _ => ind_usesArc_cons_cons u v a b rest ha))]
    rw [rsum_add, rsum_rsum_point n a b (hlt a List.mem_cons_self)
      (hlt b (List.mem_cons_of_mem a List.mem_cons_self)), ih]
    simp only [List.length_cons]
    omega

/-- **Double counting.** `Σ_{u,v} c(u → v) = (n - 1) · H(T)`. -/
theorem sum_arcCount (n : Nat) (T : Nat → Nat → Bool) :
    rsum n (fun u => rsum n (fun v => arcCount n T u v)) = (n - 1) * hpCount n T := by
  unfold arcCount hpCount
  have h1 : ∀ u, rsum n (fun v => ((hamPathsR n T).filter (usesArc u v)).length) =
      ((hamPathsR n T).map (fun p => rsum n (fun v => ind (usesArc u v p)))).sum := by
    intro u
    rw [rsum_congr n _ (fun v => ((hamPathsR n T).map (fun p => ind (usesArc u v p))).sum)
      (fun v _ => length_filter_eq_sum (usesArc u v) _)]
    exact rsum_listSum n (fun v p => ind (usesArc u v p)) _
  rw [rsum_congr n _ _ (fun u _ => h1 u)]
  rw [rsum_listSum n (fun u p => rsum n (fun v => ind (usesArc u v p))) (hamPathsR n T)]
  rw [List.map_congr_left (g := fun _ => n - 1) (fun p hp => by
    obtain ⟨hl, hnd, hlt, _⟩ := (mem_hamPathsR n T p).1 hp
    rw [arcs_of_path n p hnd hlt, hl])]
  exact sum_map_const (n - 1) _

/-! ## Arcs of a tournament -/

/-- `u → v` is an arc (with `u ≠ v`). -/
def isArc (T : Nat → Nat → Bool) (u v : Nat) : Bool := !Nat.beq u v && T u v

theorem rsum_one (n : Nat) : rsum n (fun _ => 1) = n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show rsum n (fun _ => 1) + 1 = n + 1
    rw [ih]

/-- A tournament on `n` vertices has `pairs n = C(n, 2)` arcs. -/
theorem arcTotal : ∀ (n : Nat) (T : Nat → Nat → Bool), IsTournament n T →
    rsum n (fun u => rsum n (fun v => ind (isArc T u v))) = pairs n
  | 0, _, _ => rfl
  | n + 1, T, hT => by
    have hT' : IsTournament n T := fun a b ha hb hab => hT a b (by omega) (by omega) hab
    have ih := arcTotal n T hT'
    show rsum n (fun u => rsum n (fun v => ind (isArc T u v)) + ind (isArc T u n)) +
      (rsum n (fun v => ind (isArc T n v)) + ind (isArc T n n)) = pairs n + n
    have hnn : ind (isArc T n n) = 0 := by
      unfold isArc
      rw [Nat.beq_refl]
      rfl
    rw [rsum_add, ih, hnn, Nat.add_zero, Nat.add_assoc, ← rsum_add]
    rw [rsum_congr n _ (fun _ => 1) (fun u hu => by
      have hne : u ≠ n := by omega
      unfold isArc
      rw [beq_false_of_ne u n hne, beq_false_of_ne n u (Ne.symm hne), Bool.not_false,
        Bool.true_and, Bool.true_and, hT u n (by omega) (by omega) hne]
      cases T u n <;> rfl), rsum_one]

/-! ## Parity -/

theorem rsum_mod_two_congr (n : Nat) (f g : Nat → Nat) (h : ∀ r, r < n → f r % 2 = g r % 2) :
    rsum n f % 2 = rsum n g % 2 := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show (rsum n f + f n) % 2 = (rsum n g + g n) % 2
    have h1 := ih (fun r hr => h r (by omega))
    have h2 := h n (by omega)
    omega

theorem arcCount_eq_zero_of_not (n : Nat) (T : Nat → Nat → Bool) (u v : Nat)
    (h : ∀ p, IsHamPath n T p → usesArc u v p = false) : arcCount n T u v = 0 := by
  unfold arcCount
  have := length_filter_mono (usesArc u v) (fun _ => false) (hamPathsR n T) (fun p hp hu => by
    rw [h p ((mem_hamPathsR n T p).1 hp)] at hu
    exact hu)
  rw [length_filter_false] at this
  omega

/-- `c(u → v) mod 2` is the indicator of an arc whenever every arc count is odd. -/
theorem arcCount_mod_two_of_allOdd (n : Nat) (T : Nat → Nat → Bool)
    (hodd : ∀ u v, u < n → v < n → u ≠ v → T u v = true → arcCount n T u v % 2 = 1)
    (u v : Nat) (hu : u < n) (hv : v < n) : arcCount n T u v % 2 = ind (isArc T u v) % 2 := by
  by_cases huv : u = v
  · subst huv
    rw [arcCount_eq_zero_of_not n T u u (fun p hp => usesArc_self u p hp.2.1)]
    unfold isArc
    rw [Nat.beq_refl]
    rfl
  · cases hT : T u v
    · rw [arcCount_eq_zero_of_not n T u v (fun p hp => by
        cases h : usesArc u v p
        · rfl
        · rw [usesArc_arc T u v p hp.2.2.2 h] at hT
          exact Bool.noConfusion hT)]
      unfold isArc
      rw [hT, Bool.and_false]
      rfl
    · rw [hodd u v hu hv huv hT]
      unfold isArc
      rw [hT, beq_false_of_ne u v huv]
      rfl

theorem pairs_parity : ∀ (n : Nat), pairs n % 2 = (if n % 4 = 2 ∨ n % 4 = 3 then 1 else 0)
  | 0 => rfl
  | n + 1 => by
    have ih := pairs_parity n
    show (pairs n + n) % 2 = _
    by_cases h : n % 4 = 2 ∨ n % 4 = 3
    · rw [if_pos h] at ih
      by_cases h' : (n + 1) % 4 = 2 ∨ (n + 1) % 4 = 3
      · rw [if_pos h']
        omega
      · rw [if_neg h']
        omega
    · rw [if_neg h] at ih
      by_cases h' : (n + 1) % 4 = 2 ∨ (n + 1) % 4 = 3
      · rw [if_pos h']
        omega
      · rw [if_neg h']
        omega

/-- With `H(T)` odd, the sum of all arc counts is `≡ n - 1 (mod 2)`. -/
theorem sum_arcCount_mod_two (n : Nat) (T : Nat → Nat → Bool) (hH : hpCount n T % 2 = 1) :
    rsum n (fun u => rsum n (fun v => arcCount n T u v)) % 2 = (n - 1) % 2 := by
  rw [sum_arcCount, Nat.mul_mod, hH, Nat.mul_one, Nat.mod_mod]

/-- **#odd arcs.** If `H(T)` is odd, the number of ordered pairs `(u, v)` with
`c(u → v)` odd is `≡ n - 1 (mod 2)`. -/
theorem oddArcs_parity (n : Nat) (T : Nat → Nat → Bool) (hH : hpCount n T % 2 = 1) :
    rsum n (fun u => rsum n (fun v => arcCount n T u v % 2)) % 2 = (n - 1) % 2 := by
  rw [← sum_arcCount_mod_two n T hH]
  apply rsum_mod_two_congr
  intro u _
  exact rsum_mod_two_congr n _ _ (fun v _ => Nat.mod_mod _ _)

/-- **All arcs odd.** If every arc of a tournament lies on an odd number of
Hamiltonian paths and `H(T)` is odd, then `n ≡ 1, 2 (mod 4)`. -/
theorem allOdd_mod_four (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool) (hT : IsTournament n T)
    (hodd : ∀ u v, u < n → v < n → u ≠ v → T u v = true → arcCount n T u v % 2 = 1)
    (hH : hpCount n T % 2 = 1) : n % 4 = 1 ∨ n % 4 = 2 := by
  have h1 := sum_arcCount_mod_two n T hH
  have h2 : rsum n (fun u => rsum n (fun v => arcCount n T u v)) % 2 =
      rsum n (fun u => rsum n (fun v => ind (isArc T u v))) % 2 := by
    apply rsum_mod_two_congr
    intro u hu
    exact rsum_mod_two_congr n _ _ (fun v hv => arcCount_mod_two_of_allOdd n T hodd u v hu hv)
  rw [h2, arcTotal n T hT, pairs_parity] at h1
  by_cases h : n % 4 = 2 ∨ n % 4 = 3
  · rw [if_pos h] at h1
    omega
  · rw [if_neg h] at h1
    omega

/-- **All arcs even.** If every arc lies on an even number of Hamiltonian paths and
`H(T)` is odd, then `n` is odd. -/
theorem allEven_odd (n : Nat) (hn : 1 ≤ n) (T : Nat → Nat → Bool)
    (heven : ∀ u v, u < n → v < n → arcCount n T u v % 2 = 0)
    (hH : hpCount n T % 2 = 1) : n % 2 = 1 := by
  have h1 := sum_arcCount_mod_two n T hH
  have h2 : rsum n (fun u => rsum n (fun v => arcCount n T u v)) % 2 =
      rsum n (fun u => rsum n (fun _ => 0)) % 2 := by
    apply rsum_mod_two_congr
    intro u hu
    exact rsum_mod_two_congr n _ _ (fun v hv => by rw [heven u v hu hv])
  rw [h2, rsum_congr n _ (fun _ => 0) (fun u _ => rsum_zero n), rsum_zero] at h1
  omega

/-! ## Shaving one arc -/

/-- Delete the arc `u → v`. -/
def deleteArc (T : Nat → Nat → Bool) (u v : Nat) (a b : Nat) : Bool :=
  T a b && !(Nat.beq a u && Nat.beq b v)

theorem isPath_deleteArc (T : Nat → Nat → Bool) (u v : Nat) :
    ∀ (p : List Nat), IsPath (deleteArc T u v) p ↔ IsPath T p ∧ usesArc u v p = false
  | [] => ⟨fun _ => ⟨trivial, rfl⟩, fun _ => trivial⟩
  | [_] => ⟨fun _ => ⟨trivial, rfl⟩, fun _ => trivial⟩
  | a :: b :: rest => by
    rw [isPath_cons_cons, isPath_cons_cons, isPath_deleteArc T u v (b :: rest)]
    show (deleteArc T u v a b = true ∧ IsPath T (b :: rest) ∧ usesArc u v (b :: rest) = false) ↔
      (T a b = true ∧ IsPath T (b :: rest)) ∧
        ((Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)) = false
    unfold deleteArc
    rw [Bool.and_eq_true, Bool.or_eq_false_iff, Bool.not_eq_true']
    constructor
    · rintro ⟨⟨h1, h2⟩, h3, h4⟩
      exact ⟨⟨h1, h3⟩, h2, h4⟩
    · rintro ⟨⟨h1, h3⟩, h2, h4⟩
      exact ⟨⟨h1, h2⟩, h3, h4⟩

/-- **Shaving one arc** (THM-4524 §3): deleting `u → v` from `T` leaves exactly the
`H(T) - c(u → v)` Hamiltonian paths not through it. -/
theorem shave_arc (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) :
    hpCount n (deleteArc T u v) + arcCount n T u v = hpCount n T := by
  have hsplit := length_filter_add_not (usesArc u v) (hamPathsR n T)
  have heq : hpCount n (deleteArc T u v) =
      ((hamPathsR n T).filter (fun p => !usesArc u v p)).length := by
    apply (hpCount_unique n (deleteArc T u v) _
      (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n T)) _).symm
    intro p
    rw [List.mem_filter, mem_hamPathsR]
    unfold IsHamPath
    rw [isPath_deleteArc, Bool.not_eq_true']
    constructor
    · rintro ⟨⟨h1, h2, h3, h4⟩, h5⟩
      exact ⟨h1, h2, h3, h4, h5⟩
    · rintro ⟨h1, h2, h3, h4, h5⟩
      exact ⟨⟨h1, h2, h3, h4⟩, h5⟩
  unfold arcCount hpCount at *
  rw [heq]
  omega

end ProcgenSelfieEdim
