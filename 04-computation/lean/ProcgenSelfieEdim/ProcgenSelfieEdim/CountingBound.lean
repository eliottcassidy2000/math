import ProcgenSelfieEdim.EdimQ6

set_option autoImplicit false

/-!
# THM-4525 (c): Lemma L4, no set of at most six vertices resolves `Q_6`

The counting argument of the lane note (section 3, L4):
* an edge with `H_e(0) ≥ 1` contains a landmark, and an edge with `H_e(5) ≥ 1`
  contains the antipode of a landmark (L2); each vertex `s` of `Q_6` lies in 6 edges
  and so does its antipode, so at most `12 m` edges are "bad";
* the remaining `≥ 192 - 12 m` edges have histograms `(0, a, b, c, d, 0)` with
  `a + b + c + d = m`, and there are only `C(m+3, 3)` of those;
* `192 - 12 m > C(m+3, 3)` for `m ≤ 6`, so two good edges share a histogram.

The landmark list may even contain repetitions (histograms then count with
multiplicity); in particular no set of at most six vertices resolves `Q_6`, i.e.
`edim_m(Q_6) ≥ 7`.
-/

namespace ProcgenSelfieEdim

/-! ## Vertices at distance 0 and `d` -/

theorem eq_antipode_of_hamming_eq (d a s : Nat) (ha : a < 2 ^ d) (hs : s < 2 ^ d)
    (h : hamming d a s = d) : a = antipode d s := by
  have h1 := hamming_antipode d s a hs
  rw [hamming_comm d s a, h, Nat.sub_self] at h1
  exact (eq_of_hamming_eq_zero d (antipode d s) a (antipode_lt d s hs) ha h1).symm

/-! ## The bad edges -/

/-- `s` or its antipode is an end of the edge `e` of `Q_6`. -/
def touches (s : Nat) (e : Nat × Nat) : Bool :=
  Nat.beq s e.1 || Nat.beq s e.2 || Nat.beq (antipode 6 s) e.1 || Nat.beq (antipode 6 s) e.2

/-- Each vertex of `Q_6`, together with its antipode, touches at most 12 edges
(finite check over the 64 vertices). -/
theorem touch_count_check :
    (List.range 64).all (fun s => Nat.ble ((edgeList 6).filter (touches s)).length 12) = true := by
  decide +kernel

theorem touch_count (s : Nat) (hs : s < 2 ^ 6) : ((edgeList 6).filter (touches s)).length ≤ 12 :=
  Nat.le_of_ble_eq_true (forall_lt_of_all_range _ 64 touch_count_check s hs)

/-- An edge with a landmark at level `0` or level `5` is touched by that landmark. -/
theorem touched_of_bad (S : List Nat) (hS : ∀ s ∈ S, s < 2 ^ 6) (a b : Nat)
    (hab : (a, b) ∈ edgeList 6) (hbad : hist 6 S a b 0 ≠ 0 ∨ hist 6 S a b 5 ≠ 0) :
    S.any (fun s => touches s (a, b)) = true := by
  obtain ⟨⟨ha, hb, h1⟩, _⟩ := isEdge_of_mem_edgeList 6 a b hab
  rw [List.any_eq_true]
  rcases hbad with h0 | h5
  · obtain ⟨s, hs, hd⟩ := exists_of_length_filter_ne_zero _ S h0
    refine ⟨s, hs, ?_⟩
    have hd' : edgeDist 6 a b s = 0 := beq_iff_eq.1 hd
    have hs64 := hS s hs
    unfold edgeDist at hd'
    by_cases hz : hamming 6 a s = 0
    · have := eq_of_hamming_eq_zero 6 a s ha hs64 hz
      subst this
      show (Nat.beq a a || _ || _ || _) = true
      rw [Nat.beq_refl]
      rfl
    · have hz' : hamming 6 b s = 0 := by omega
      have := eq_of_hamming_eq_zero 6 b s hb hs64 hz'
      subst this
      show (Nat.beq b a || Nat.beq b b || _ || _) = true
      rw [Nat.beq_refl, Bool.or_true]
      rfl
  · obtain ⟨s, hs, hd⟩ := exists_of_length_filter_ne_zero _ S h5
    refine ⟨s, hs, ?_⟩
    have hd' : edgeDist 6 a b s = 5 := beq_iff_eq.1 hd
    have hs64 := hS s hs
    have hadj := hamming_adjacent 6 a b s h1
    have hla := hamming_le 6 a s
    have hlb := hamming_le 6 b s
    unfold edgeDist at hd'
    by_cases h6 : hamming 6 a s = 6
    · have := eq_antipode_of_hamming_eq 6 a s ha hs64 h6
      subst this
      show (_ || Nat.beq (antipode 6 s) (antipode 6 s) || _) = true
      rw [Nat.beq_refl, Bool.or_true]
      rfl
    · have h6' : hamming 6 b s = 6 := by omega
      have := eq_antipode_of_hamming_eq 6 b s hb hs64 h6'
      subst this
      show (_ || Nat.beq (antipode 6 s) (antipode 6 s)) = true
      rw [Nat.beq_refl, Bool.or_true]

/-! ## The good edges and their histograms -/

/-- The histogram of `e` as the explicit list `[H(0), …, H(5)]`. -/
def histVec6 (S : List Nat) (e : Nat × Nat) : List Nat :=
  [hist 6 S e.1 e.2 0, hist 6 S e.1 e.2 1, hist 6 S e.1 e.2 2, hist 6 S e.1 e.2 3,
    hist 6 S e.1 e.2 4, hist 6 S e.1 e.2 5]

/-- An edge is good if no landmark is at level 0 or 5. -/
def good (S : List Nat) (e : Nat × Nat) : Bool :=
  Nat.beq (hist 6 S e.1 e.2 0) 0 && Nat.beq (hist 6 S e.1 e.2 5) 0

/-- The histograms `(0, a, b, c, d, 0)` with `a + b + c + d = m`. -/
def comps (m : Nat) : List (List Nat) :=
  (List.range (m + 1)).flatMap (fun a =>
    (List.range (m + 1 - a)).flatMap (fun b =>
      (List.range (m + 1 - a - b)).map (fun c => [0, a, b, c, m - a - b - c, 0])))

theorem mem_comps (m a b c d : Nat) (h : a + b + c + d = m) : [0, a, b, c, d, 0] ∈ comps m := by
  unfold comps
  rw [List.mem_flatMap]
  refine ⟨a, List.mem_range.2 (by omega), ?_⟩
  rw [List.mem_flatMap]
  refine ⟨b, List.mem_range.2 (by omega), ?_⟩
  rw [List.mem_map]
  refine ⟨c, List.mem_range.2 (by omega), ?_⟩
  have hd : m - a - b - c = d := by omega
  rw [hd]

/-- `|comps m| = C(m+3, 3)` is smaller than `192 - 12 m` for `m ≤ 6`. -/
theorem comps_small : (List.range 7).all (fun m => Nat.blt ((comps m).length + 12 * m) 192) = true := by
  decide +kernel

theorem rsum_six (f : Nat → Nat) : rsum 6 f = f 0 + f 1 + f 2 + f 3 + f 4 + f 5 := by
  simp only [rsum]
  omega

theorem histVec6_mem_comps (S : List Nat) (e : Nat × Nat) (he : e ∈ edgeList 6)
    (hg : good S e = true) : histVec6 S e ∈ comps S.length := by
  have hedge := (isEdge_of_mem_edgeList 6 e.1 e.2 he).1
  have hsum := rsum_hist 6 S e.1 e.2 hedge.2.2
  rw [rsum_six] at hsum
  unfold good at hg
  rw [Bool.and_eq_true, Nat.beq_eq, Nat.beq_eq] at hg
  unfold histVec6
  rw [hg.1, hg.2]
  exact mem_comps _ _ _ _ _ (by rw [hg.1, hg.2] at hsum; omega)

/-- **L4.** No list of at most six vertices (in particular no set of at most six
vertices) is edge-multiset resolving in `Q_6`. Hence `edim_m(Q_6) ≥ 7`. -/
theorem not_resolving_of_length_le_six (S : List Nat) (hS : ∀ s ∈ S, s < 2 ^ 6)
    (hm : S.length ≤ 6) : ¬ Resolving 6 S := by
  intro hres
  -- at most `12 m` bad edges
  have hbad : ((edgeList 6).filter (fun e => !good S e)).length ≤ 12 * S.length := by
    have h1 := length_filter_mono (fun e => !good S e) (fun e => S.any (fun s => touches s e))
      (edgeList 6) (fun e he hb => by
        rcases e with ⟨a, b⟩
        have hb' : (!good S (a, b)) = true := hb
        apply touched_of_bad S hS a b he
        by_cases h0 : hist 6 S a b 0 = 0
        · refine Or.inr (fun h5 => ?_)
          have hg : good S (a, b) = true := by
            unfold good
            rw [Bool.and_eq_true, Nat.beq_eq, Nat.beq_eq]
            exact ⟨h0, h5⟩
          rw [hg] at hb'
          exact Bool.noConfusion hb'
        · exact Or.inl h0)
    have h2 := length_filter_any_le touches (edgeList 6) S
    have h3 := sum_le_mul (fun s => ((edgeList 6).filter (touches s)).length) 12 S
      (fun s hs => touch_count s (hS s hs))
    omega
  have hsplit := length_filter_add_not (good S) (edgeList 6)
  have hlen : (edgeList 6).length = 192 := edgeList_six.1
  -- the good edges inject into `comps m`
  have hgood : ((edgeList 6).filter (good S)).Nodup :=
    List.Nodup.sublist List.filter_sublist edgeList_six.2
  have hinj : (((edgeList 6).filter (good S)).map (histVec6 S)).Nodup := by
    apply nodup_map_of_inj_on _ _ hgood
    intro e he e' he' heq
    rcases e with ⟨a, b⟩
    rcases e' with ⟨a', b'⟩
    rw [List.mem_filter] at he he'
    obtain ⟨hE, hlt⟩ := isEdge_of_mem_edgeList 6 a b he.1
    obtain ⟨hE', hlt'⟩ := isEdge_of_mem_edgeList 6 a' b' he'.1
    simp only [histVec6, List.cons.injEq] at heq
    obtain ⟨e0, e1, e2, e3, e4, e5, -⟩ := heq
    have hh : ∀ r, hist 6 S a b r = hist 6 S a' b' r := by
      intro r
      match r with
      | 0 => exact e0
      | 1 => exact e1
      | 2 => exact e2
      | 3 => exact e3
      | 4 => exact e4
      | 5 => exact e5
      | r + 6 =>
        rw [hist_eq_zero_of_le 6 a b (r + 6) hE.2.2 (by omega),
          hist_eq_zero_of_le 6 a' b' (r + 6) hE'.2.2 (by omega)]
    rcases hres a b a' b' hE hE' hh with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · rw [h1, h2]
    · exact absurd hlt (by omega)
  have hsub : ∀ x ∈ ((edgeList 6).filter (good S)).map (histVec6 S), x ∈ comps S.length := by
    intro x hx
    rw [List.mem_map] at hx
    obtain ⟨e, he, rfl⟩ := hx
    rw [List.mem_filter] at he
    exact histVec6_mem_comps S e he.1 he.2
  have hle := pigeonhole _ _ hinj hsub
  rw [List.length_map] at hle
  have hc := forall_lt_of_all_range _ 7 comps_small S.length (by omega)
  rw [Nat.blt_eq] at hc
  omega

end ProcgenSelfieEdim
