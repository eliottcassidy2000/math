import ProcgenSelfieEdim.Redei

set_option autoImplicit false

/-!
# THM-4524 §3: the parity-break number of an all-even tournament is 2

`ρ(T)` is the least number of arcs whose deletion makes `H` even. Single deletions:
`H(T - e)` is even iff `c(e)` is odd (`shave_keeps_odd_iff`), so `ρ = 1` iff some arc is odd.
For an all-even tournament (`n ≥ 2`) no single deletion works, and two do
(`parity_break_two`): `end(v)` is odd, so some `u → v` ends an odd number of Hamiltonian
paths; then `Σ_w c(u → v → w)` is odd, some `c(u → v → w)` is odd, and deleting
`e = u → v`, `f = v → w` leaves `H - c(e) - c(f) + c(e, f)` paths, an even number
(`delete_two`).
-/

namespace ProcgenSelfieEdim

/-- Hamiltonian paths ending at `v`. -/
noncomputable def endCount (n : Nat) (T : Nat → Nat → Bool) (v : Nat) : Nat :=
  ((hamPathsR n T).filter (lastIs v)).length

/-- Hamiltonian paths ending with the arc `u → v`. -/
noncomputable def endArcCount (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) : Nat :=
  ((hamPathsR n T).filter (fun p => usesArc u v p && lastIs v p)).length

/-- Hamiltonian paths using `u → v → w`. -/
noncomputable def pairCount (n : Nat) (T : Nat → Nat → Bool) (u v w : Nat) : Nat :=
  ((hamPathsR n T).filter (fun p => usesArc u v p && usesArc v w p)).length

/-- A vertex of a duplicate-free path is last or has exactly one successor. -/
theorem succ_or_last (n : Nat) (v : Nat) (p : List Nat) (hnd : p.Nodup) (hlt : ∀ x ∈ p, x < n)
    (hv : v ∈ p) : rsum n (fun w => ind (usesArc v w p)) + ind (lastIs v p) = 1 := by
  obtain ⟨L, R, rfl⟩ := List.append_of_mem hv
  have hvR : v ∉ R := fun h => (List.nodup_cons.1
    (List.Nodup.sublist (List.sublist_append_right L _) hnd)).1 h
  cases R with
  | nil =>
    rw [rsum_congr n _ (fun _ => 0) (fun w _ => by
      cases h : usesArc v w (L ++ [v])
      · rfl
      · obtain ⟨R', hR'⟩ := (usesArc_from v w L [] hnd).1 h
        exact nomatch hR'), rsum_zero, lastIs_snoc, Nat.beq_refl]
    rfl
  | cons w0 R' =>
    have hw0 : w0 < n := hlt w0 (List.mem_append_right _ (List.mem_cons_of_mem _ List.mem_cons_self))
    rw [rsum_congr n _ (fun w => if w0 = w then 1 else 0) (fun w _ => by
      dsimp only
      have e := usesArc_from v w L (w0 :: R') hnd
      by_cases hw : w0 = w
      · rw [if_pos hw, (e.2 ⟨R', by rw [hw]⟩)]
        rfl
      · rw [if_neg hw]
        cases h : usesArc v w (L ++ v :: w0 :: R')
        · rfl
        · obtain ⟨R'', hR''⟩ := e.1 h
          exact absurd (List.cons.inj hR'').1 hw)]
    rw [rsum_single n w0 (fun _ => 1) hw0]
    obtain ⟨R'', z, hz⟩ : ∃ R'' z, w0 :: R' = R'' ++ [z] := by
      rcases List.eq_nil_or_concat (w0 :: R') with h | ⟨R'', z, h⟩
      · exact absurd h (List.cons_ne_nil _ _)
      · exact ⟨R'', z, by rw [h, List.concat_eq_append]⟩
    have hzR : z ∈ w0 :: R' := by rw [hz]; exact List.mem_append_right _ List.mem_cons_self
    rw [hz, show L ++ v :: (R'' ++ [z]) = (L ++ v :: R'') ++ [z] by
      simp only [List.append_assoc, List.cons_append], lastIs_snoc,
      beq_false_of_ne z v (fun h => hvR (h ▸ hzR))]
    rfl

/-- In a Hamiltonian path on `n ≥ 2` vertices, a last vertex has exactly one predecessor. -/
theorem pred_of_last (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (v : Nat) (p : List Nat)
    (hp : IsHamPath n T p) :
    rsum n (fun u => ind (usesArc u v p && lastIs v p)) = ind (lastIs v p) := by
  cases hl : lastIs v p
  · exact (rsum_congr n _ (fun _ => 0) (fun u _ => by rw [Bool.and_false]; rfl)).trans (rsum_zero n)
  · -- `p = L ++ [v]` with `L` nonempty
    obtain ⟨L, hL⟩ : ∃ L, p = L ++ [v] := by
      cases p with
      | nil => exact Bool.noConfusion hl
      | cons a rest =>
        have hlast : lastOr a rest = v := Nat.eq_of_beq_eq_true hl
        cases rest with
        | nil => exact ⟨[], by rw [List.nil_append]; exact congrArg (fun x => [x]) hlast⟩
        | cons r rest' =>
          obtain ⟨mid, hmid⟩ := eq_append_lastOr (r :: rest') a (List.cons_ne_nil _ _)
          rw [lastOr_cons, show lastOr r rest' = lastOr a (r :: rest') from rfl, hlast] at hmid
          exact ⟨a :: mid, by rw [hmid]; rfl⟩
    subst hL
    have hnd := hp.2.1
    rcases List.eq_nil_or_concat L with h0 | ⟨L', u0, h0⟩
    · exfalso
      subst h0
      have := hp.1
      rw [List.nil_append, List.length_singleton] at this
      omega
    · rw [List.concat_eq_append] at h0
      subst h0
      have hu0 : u0 < n := hp.2.2.1 u0 (List.mem_append_left _ (List.mem_append_right _ List.mem_cons_self))
      rw [rsum_congr n _ (fun u => if u0 = u then 1 else 0) (fun u _ => by
        dsimp only
        rw [Bool.and_true]
        have e := usesArc_into v u (L' ++ [u0]) [] hnd
        by_cases hu : u0 = u
        · rw [if_pos hu, e.2 ⟨L', by rw [hu]⟩]
          rfl
        · rw [if_neg hu]
          cases h : usesArc u v ((L' ++ [u0]) ++ [v])
          · rfl
          · obtain ⟨L'', hL''⟩ := e.1 h
            exact absurd (List.cons.inj (List.append_inj' hL'' rfl).2).1 hu)]
      rw [rsum_single n u0 (fun _ => 1) hu0]
      rfl

/-- `Σ_w c(v → w) + end(v) = H`. -/
theorem out_arcs_sum (n : Nat) (T : Nat → Nat → Bool) (v : Nat) (hv : v < n) :
    rsum n (fun w => arcCount n T v w) + endCount n T v = hpCount n T := by
  unfold arcCount endCount hpCount
  rw [rsum_congr n _ _ (fun w _ => length_filter_eq_sum (usesArc v w) _),
    rsum_listSum n (fun w p => ind (usesArc v w p)), length_filter_eq_sum, ← sum_map_add,
    ← sum_map_one' (hamPathsR n T)]
  apply congrArg List.sum
  apply List.map_congr_left
  intro p hp
  have hham := (mem_hamPathsR n T p).1 hp
  exact succ_or_last n v p hham.2.1 hham.2.2.1 (mem_of_hamPath n T p hham v hv)

/-- `end(v) = Σ_u #(paths ending with u → v)` for `n ≥ 2`. -/
theorem end_sum (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (v : Nat) :
    endCount n T v = rsum n (fun u => endArcCount n T u v) := by
  unfold endCount endArcCount
  rw [rsum_congr n _ _ (fun u _ => length_filter_eq_sum _ _),
    rsum_listSum n (fun u p => ind (usesArc u v p && lastIs v p)), length_filter_eq_sum]
  apply congrArg List.sum
  apply List.map_congr_left
  intro p hp
  exact (pred_of_last n hn T v p ((mem_hamPathsR n T p).1 hp)).symm

/-- `c(u → v) = #(paths ending with u → v) + Σ_w c(u → v → w)`. -/
theorem arc_split (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) :
    arcCount n T u v = endArcCount n T u v + rsum n (fun w => pairCount n T u v w) := by
  unfold arcCount endArcCount pairCount
  rw [rsum_congr n _ _ (fun w _ => length_filter_eq_sum _ _),
    rsum_listSum n (fun w p => ind (usesArc u v p && usesArc v w p)), length_filter_eq_sum,
    length_filter_eq_sum, ← sum_map_add]
  apply congrArg List.sum
  apply List.map_congr_left
  intro p hp
  have hham := (mem_hamPathsR n T p).1 hp
  cases huv : usesArc u v p
  · rw [show ind false = 0 from rfl,
      rsum_congr n (fun w => ind (false && usesArc v w p)) (fun _ => 0) (fun w _ => rfl), rsum_zero]
    rfl
  · have hvp : v ∈ p := (usesArc_mem u v p huv).2
    have := succ_or_last n v p hham.2.1 hham.2.2.1 hvp
    rw [show ind true = 1 from rfl]
    rw [rsum_congr n (fun w => ind (true && usesArc v w p)) (fun w => ind (usesArc v w p))
      (fun w _ => rfl), Bool.true_and]
    omega

/-- **Deleting two consecutive arcs** (inclusion-exclusion):
`H(T - e - f) + c(e) + c(f) = H(T) + c(e, f)` for `e = u → v`, `f = v → w`. -/
theorem delete_two (n : Nat) (T : Nat → Nat → Bool) (u v w : Nat) :
    hpCount n (deleteArc (deleteArc T u v) v w) + arcCount n T u v + arcCount n T v w =
      hpCount n T + pairCount n T u v w := by
  have s1 := shave_arc n T u v
  have s2 := shave_arc n (deleteArc T u v) v w
  -- paths of `T - e` through `f` are the paths of `T` through `f` avoiding `e`
  have hP : pairCount n T u v w =
      (((hamPathsR n T).filter (usesArc v w)).filter (fun p => usesArc u v p)).length := by
    unfold pairCount
    rw [List.filter_filter]
  have hA : arcCount n T v w = ((hamPathsR n T).filter (usesArc v w)).length := rfl
  have hsplit : (((hamPathsR n T).filter (usesArc v w)).filter (fun p => usesArc u v p)).length +
      (((hamPathsR n T).filter (usesArc v w)).filter (fun p => !usesArc u v p)).length =
      ((hamPathsR n T).filter (usesArc v w)).length :=
    length_filter_add_not (fun p => usesArc u v p) _
  have h2 : (((hamPathsR n T).filter (usesArc v w)).filter (fun p => !usesArc u v p)).length =
      arcCount n (deleteArc T u v) v w := by
    unfold arcCount
    apply length_eq_of_same_members _ _
      (List.Nodup.sublist (List.Sublist.trans List.filter_sublist List.filter_sublist)
        (nodup_hamPathsR n T))
      (List.Nodup.sublist List.filter_sublist (nodup_hamPathsR n _))
    intro p
    rw [List.mem_filter, List.mem_filter, List.mem_filter, mem_hamPathsR, mem_hamPathsR]
    unfold IsHamPath
    rw [isPath_deleteArc, Bool.not_eq_true']
    constructor
    · rintro ⟨⟨⟨h1, h2, h3, h4⟩, h5⟩, h6⟩
      exact ⟨⟨h1, h2, h3, h4, h6⟩, h5⟩
    · rintro ⟨⟨h1, h2, h3, h4, h6⟩, h5⟩
      exact ⟨⟨⟨h1, h2, h3, h4⟩, h5⟩, h6⟩
  omega

/-- **THM-4524 §3: `ρ = 2` for all-even tournaments.** If every arc of a tournament on
`n ≥ 2` vertices lies on an even number of Hamiltonian paths, then no single arc deletion
makes `H` even, but deleting two consecutive arcs `u → v`, `v → w` does. -/
theorem parity_break_two (n : Nat) (hn : 2 ≤ n) (T : Nat → Nat → Bool) (hT : IsTournament n T)
    (heven : ∀ u v, u < n → v < n → arcCount n T u v % 2 = 0) :
    (∀ u v, u < n → v < n → hpCount n (deleteArc T u v) % 2 = 1) ∧
      ∃ u v w, u < n ∧ v < n ∧ w < n ∧ T u v = true ∧ T v w = true ∧
        hpCount n (deleteArc (deleteArc T u v) v w) % 2 = 0 := by
  refine ⟨fun u v hu hv => (shave_keeps_odd_iff n T hT u v).2 (heven u v hu hv), ?_⟩
  have hH := redei n T hT
  -- `end(0)` is odd
  have hend : endCount n T 0 % 2 = 1 := by
    have h := out_arcs_sum n T 0 (by omega)
    have h2 : rsum n (fun w => arcCount n T 0 w) % 2 = 0 := by
      rw [rsum_mod_two_congr n _ (fun _ => 0) (fun w hw => by rw [heven 0 w (by omega) hw]; rfl),
        rsum_zero]
    omega
  -- some `u → 0` ends an odd number of paths
  have hu : ∃ u, u < n ∧ endArcCount n T u 0 % 2 = 1 := by
    apply Decidable.byContradiction
    intro hno
    have hall : ∀ u, u < n → endArcCount n T u 0 % 2 = 0 := by
      intro u hu
      have : ¬ endArcCount n T u 0 % 2 = 1 := fun h => hno ⟨u, hu, h⟩
      omega
    rw [end_sum n hn T 0, rsum_mod_two_congr n _ (fun _ => 0) (fun u hu => by rw [hall u hu]; rfl),
      rsum_zero] at hend
    exact absurd hend (by decide)
  obtain ⟨u, hun, huodd⟩ := hu
  -- then some `u → 0 → w` lies on an odd number of paths
  have hw : ∃ w, w < n ∧ pairCount n T u 0 w % 2 = 1 := by
    apply Decidable.byContradiction
    intro hno
    have hall : ∀ w, w < n → pairCount n T u 0 w % 2 = 0 := by
      intro w hw
      have : ¬ pairCount n T u 0 w % 2 = 1 := fun h => hno ⟨w, hw, h⟩
      omega
    have hs := arc_split n T u 0
    have h0 : rsum n (fun w => pairCount n T u 0 w) % 2 = 0 := by
      rw [rsum_mod_two_congr n _ (fun _ => 0) (fun w hw => by rw [hall w hw]; rfl), rsum_zero]
    have := heven u 0 hun (by omega)
    omega
  obtain ⟨w, hwn, hwodd⟩ := hw
  -- the two arcs exist
  have harc1 : T u 0 = true := by
    have hpos : 0 < endArcCount n T u 0 := by omega
    unfold endArcCount at hpos
    obtain ⟨p, hp⟩ : ∃ p, p ∈ (hamPathsR n T).filter (fun p => usesArc u 0 p && lastIs 0 p) := by
      cases h : (hamPathsR n T).filter (fun p => usesArc u 0 p && lastIs 0 p) with
      | nil => rw [h, List.length_nil] at hpos; exact absurd hpos (Nat.lt_irrefl 0)
      | cons p _ => exact ⟨p, List.mem_cons_self⟩
    rw [List.mem_filter, mem_hamPathsR, Bool.and_eq_true] at hp
    exact usesArc_arc T u 0 p hp.1.2.2.2 hp.2.1
  have harc2 : T 0 w = true := by
    have hpos : 0 < pairCount n T u 0 w := by omega
    unfold pairCount at hpos
    obtain ⟨p, hp⟩ : ∃ p, p ∈ (hamPathsR n T).filter (fun p => usesArc u 0 p && usesArc 0 w p) := by
      cases h : (hamPathsR n T).filter (fun p => usesArc u 0 p && usesArc 0 w p) with
      | nil => rw [h, List.length_nil] at hpos; exact absurd hpos (Nat.lt_irrefl 0)
      | cons p _ => exact ⟨p, List.mem_cons_self⟩
    rw [List.mem_filter, mem_hamPathsR, Bool.and_eq_true] at hp
    exact usesArc_arc T 0 w p hp.1.2.2.2 hp.2.2
  refine ⟨u, 0, w, hun, by omega, hwn, harc1, harc2, ?_⟩
  have hd := delete_two n T u 0 w
  have e1 := heven u 0 hun (by omega)
  have e2 := heven 0 w (by omega) hwn
  omega

end ProcgenSelfieEdim
