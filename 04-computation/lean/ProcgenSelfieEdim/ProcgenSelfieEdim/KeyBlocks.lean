import ProcgenSelfieEdim.Hypercube

set_option autoImplicit false

/-!
# Certificates for larger hypercubes: keys in blocks, distinctness by sorting

For `Q_7, Q_8, …` a single-declaration check is too expensive for the kernel. A resolving
certificate is checked in three steps:

1. the one-pass edge keys `Σ_{s ∈ S} base^{d(e,s)}` (any `base`; soundness only needs the key
   to be a function of the histogram, `resolving_of_edgeKeys`) are compared, 64 edges at a
   time, with an explicit list of numbers (one declaration per block, `eq_of_blocks`);
2. the explicit list is shown duplicate-free: a verified merge sort (`msortR_perm`) and a
   strict-increase check (`nodup_of_msortR`);
3. hence the keys are pairwise distinct and the set is resolving.
-/

namespace ProcgenSelfieEdim

/-- The one-pass key `Σ_{s ∈ S} base^{d(uv, s)}`. -/
noncomputable def keyBR (base d : Nat) (S : List Nat) (u v : Nat) : Nat :=
  List.rec (motive := fun _ => Nat) 0
    (fun s _ ih => base ^ fmin (hamR d u s) (hamR d v s) + ih) S

theorem keyBR_eq (base d : Nat) (u v : Nat) :
    ∀ (S : List Nat), keyBR base d S u v = (S.map (fun s => base ^ edgeDist d u v s)).sum
  | [] => rfl
  | s :: t => by
    show base ^ fmin (hamR d u s) (hamR d v s) + keyBR base d t u v = _
    rw [keyBR_eq base d u v t, List.map_cons, List.sum_cons, fmin_eq, hamR_eq, hamR_eq]
    rfl

/-- The keys of the canonical edges of `Q_d`. -/
noncomputable def edgeKeys (base d : Nat) (S : List Nat) : List Nat :=
  (edgeList d).map (fun e => keyBR base d S e.1 e.2)

/-- Distinct keys certify resolvability (for any base). -/
theorem resolving_of_edgeKeys (base d : Nat) (S : List Nat) (h : (edgeKeys base d S).Nodup) :
    Resolving d S := by
  refine resolving_of_keys d S (fun e => keyBR base d S e.1 e.2) ?_ h
  intro a b a' b' he he' hh
  show keyBR base d S a b = keyBR base d S a' b'
  rw [keyBR_eq, keyBR_eq, weighted_sum_eq d a b he.2.2 (fun r => base ^ r) S,
    weighted_sum_eq d a' b' he'.2.2 (fun r => base ^ r) S]
  exact rsum_congr d _ _ (fun r _ => by rw [hh r])

/-! ## Blocks -/

/-- The `i`-th block of length `B`. -/
def blk (B i : Nat) (l : List Nat) : List Nat := (l.drop (i * B)).take B

theorem blk_succ (B i : Nat) (l : List Nat) : blk B (i + 1) l = blk B i (l.drop B) := by
  unfold blk
  rw [List.drop_drop, Nat.succ_mul, Nat.add_comm]

/-- Two lists with the same blocks are equal. -/
theorem eq_of_blocks (B : Nat) : ∀ (k : Nat) (l1 l2 : List Nat), l1.length ≤ k * B →
    l2.length ≤ k * B → (∀ i, i < k → blk B i l1 = blk B i l2) → l1 = l2
  | 0, l1, l2, h1, h2, _ => by
    rw [Nat.zero_mul] at h1 h2
    rw [List.length_eq_zero_iff.1 (Nat.le_zero.1 h1), List.length_eq_zero_iff.1 (Nat.le_zero.1 h2)]
  | k + 1, l1, l2, h1, h2, hb => by
    have h0 : l1.take B = l2.take B := by
      have := hb 0 (Nat.succ_pos k)
      unfold blk at this
      rw [Nat.zero_mul, List.drop_zero, List.drop_zero] at this
      exact this
    have hrest : l1.drop B = l2.drop B := by
      apply eq_of_blocks B k
      · rw [List.length_drop, Nat.succ_mul] at *
        omega
      · rw [List.length_drop, Nat.succ_mul] at *
        omega
      · intro i hi
        rw [← blk_succ, ← blk_succ]
        exact hb (i + 1) (by omega)
    rw [← List.take_append_drop B l1, ← List.take_append_drop B l2, h0, hrest]

/-- Boolean list equality through `Nat.beq`. -/
noncomputable def eqListR (l1 : List Nat) : List Nat → Bool :=
  List.rec (motive := fun _ => List Nat → Bool)
    (fun l2 => List.casesOn (motive := fun _ => Bool) l2 true (fun _ _ => false))
    (fun x _ ih l2 => List.casesOn (motive := fun _ => Bool) l2 false (fun y t => Nat.beq x y && ih t)) l1

theorem eq_of_eqListR : ∀ (l1 l2 : List Nat), eqListR l1 l2 = true → l1 = l2
  | [], [], _ => rfl
  | [], _ :: _, h => Bool.noConfusion h
  | _ :: _, [], h => Bool.noConfusion h
  | x :: t1, y :: t2, h => by
    have h' : (Nat.beq x y && eqListR t1 t2) = true := h
    rw [Bool.and_eq_true, Nat.beq_eq] at h'
    rw [h'.1, eq_of_eqListR t1 t2 h'.2]

/-! ## `O(n log n)` distinctness through a merge sort

Only "the output is a permutation of the input" is proved (`msortR_perm`); that the output
is strictly increasing is checked by computation (`strictIncR`). A wrong sort could only
make the check fail. -/

/-- Merge (fuel-bounded; with no fuel left it appends, which is still a permutation). -/
noncomputable def mergeR (fuel : Nat) : List Nat → List Nat → List Nat :=
  Nat.rec (motive := fun _ => List Nat → List Nat → List Nat)
    (fun l1 l2 => appR l1 l2)
    (fun _ ih l1 l2 => List.casesOn (motive := fun _ => List Nat) l1 l2 (fun x xs =>
      List.casesOn (motive := fun _ => List Nat) l2 (x :: xs) (fun y ys =>
        cond (Nat.ble x y) (x :: ih xs (y :: ys)) (y :: ih (x :: xs) ys))))
    fuel

theorem mergeR_perm : ∀ (fuel : Nat) (l1 l2 : List Nat), (mergeR fuel l1 l2).Perm (l1 ++ l2)
  | 0, l1, l2 => by
    show (appR l1 l2).Perm (l1 ++ l2)
    rw [appR_eq]
  | _ + 1, [], l2 => List.Perm.refl l2
  | _ + 1, x :: xs, [] => by
    show (x :: xs).Perm (x :: xs ++ [])
    rw [List.append_nil]
  | f + 1, x :: xs, y :: ys => by
    show (cond (Nat.ble x y) (x :: mergeR f xs (y :: ys)) (y :: mergeR f (x :: xs) ys)).Perm _
    cases Nat.ble x y
    · show (y :: mergeR f (x :: xs) ys).Perm (x :: xs ++ y :: ys)
      exact ((mergeR_perm f (x :: xs) ys).cons y).trans List.perm_middle.symm
    · show (x :: mergeR f xs (y :: ys)).Perm (x :: (xs ++ y :: ys))
      exact (mergeR_perm f xs (y :: ys)).cons x

/-- Merge adjacent lists (`pending` holds an unmatched list). -/
noncomputable def mergePairsR (fuel : Nat) (Ls : List (List Nat)) : Option (List Nat) → List (List Nat) :=
  List.rec (motive := fun _ => Option (List Nat) → List (List Nat))
    (fun o => Option.casesOn (motive := fun _ => List (List Nat)) o [] (fun a => [a]))
    (fun x _ ih o => Option.casesOn (motive := fun _ => List (List Nat)) o (ih (some x))
      (fun a => mergeR fuel a x :: ih none)) Ls

theorem mergePairsR_perm (fuel : Nat) : ∀ (Ls : List (List Nat)) (o : Option (List Nat)),
    (mergePairsR fuel Ls o).flatten.Perm ((o.getD []) ++ Ls.flatten)
  | [], none => List.Perm.refl _
  | [], some a => by
    show [a].flatten.Perm (a ++ [])
    rw [List.flatten_cons, List.flatten_nil]
  | x :: rest, none => by
    show (mergePairsR fuel rest (some x)).flatten.Perm ([] ++ (x :: rest).flatten)
    rw [List.nil_append, List.flatten_cons]
    exact mergePairsR_perm fuel rest (some x)
  | x :: rest, some a => by
    show (mergeR fuel a x :: mergePairsR fuel rest none).flatten.Perm (a ++ (x :: rest).flatten)
    rw [List.flatten_cons, List.flatten_cons, ← List.append_assoc]
    exact (mergeR_perm fuel a x).append (by
      have := mergePairsR_perm fuel rest none
      rw [Option.getD_none, List.nil_append] at this
      exact this)

/-- `rounds` rounds of pairwise merging, starting from singletons. -/
noncomputable def msortR (rounds : Nat) (l : List Nat) : List Nat :=
  (Nat.rec (motive := fun _ => List (List Nat)) (l.map (fun x => [x]))
    (fun _ ih => mergePairsR l.length ih none) rounds).flatten

theorem flatten_singletons : ∀ (l : List Nat), (l.map (fun x => [x])).flatten = l
  | [] => rfl
  | x :: t => by rw [List.map_cons, List.flatten_cons, flatten_singletons t]; rfl

theorem msortR_perm (rounds : Nat) (l : List Nat) : (msortR rounds l).Perm l := by
  unfold msortR
  induction rounds with
  | zero =>
    show (l.map (fun x => [x])).flatten.Perm l
    rw [flatten_singletons]
  | succ r ih =>
    have := mergePairsR_perm l.length
      (Nat.rec (motive := fun _ => List (List Nat)) (l.map (fun x => [x]))
        (fun _ ih => mergePairsR l.length ih none) r) none
    rw [Option.getD_none, List.nil_append] at this
    exact this.trans ih

/-- Strict increase, checked on adjacent pairs. -/
noncomputable def strictFrom (l : List Nat) : Nat → Bool :=
  List.rec (motive := fun _ => Nat → Bool) (fun _ => true)
    (fun x _ ih prev => Nat.blt prev x && ih x) l

noncomputable def strictIncR (l : List Nat) : Bool :=
  List.casesOn (motive := fun _ => Bool) l true (fun a rest => strictFrom rest a)

theorem strictFrom_lt : ∀ (l : List Nat) (a : Nat), strictFrom l a = true → ∀ x ∈ l, a < x
  | [], _, _, x, hx => absurd hx List.not_mem_nil
  | y :: t, a, h, x, hx => by
    have h' : (Nat.blt a y && strictFrom t y) = true := h
    rw [Bool.and_eq_true, Nat.blt_eq] at h'
    rcases List.mem_cons.1 hx with rfl | hx'
    · exact h'.1
    · exact Nat.lt_trans h'.1 (strictFrom_lt t y h'.2 x hx')

theorem nodup_of_strictFrom : ∀ (l : List Nat) (a : Nat), strictFrom l a = true → (a :: l).Nodup
  | [], _, _ => List.nodup_cons.2 ⟨List.not_mem_nil, List.Pairwise.nil⟩
  | y :: t, a, h => by
    have h' : (Nat.blt a y && strictFrom t y) = true := h
    rw [Bool.and_eq_true, Nat.blt_eq] at h'
    refine List.nodup_cons.2 ⟨fun hm => ?_, nodup_of_strictFrom t y h'.2⟩
    have := strictFrom_lt (y :: t) a h a hm
    exact Nat.lt_irrefl a this

theorem nodup_of_strictIncR (l : List Nat) (h : strictIncR l = true) : l.Nodup := by
  cases l with
  | nil => exact List.Pairwise.nil
  | cons a rest => exact nodup_of_strictFrom rest a h

/-- **Distinctness by sorting.** -/
theorem nodup_of_msortR (rounds : Nat) (l : List Nat) (h : strictIncR (msortR rounds l) = true) :
    l.Nodup :=
  (msortR_perm rounds l).nodup_iff.1 (nodup_of_strictIncR _ h)

end ProcgenSelfieEdim
