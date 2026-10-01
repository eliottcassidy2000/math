import ProcgenSelfieEdim.EdimQ6

set_option autoImplicit false

/-!
# THM-4525 (d): Lemma L1, resolving sets have trivial stabilizer

A graph automorphism of `Q_d` is a bijection of the vertex set `{0, …, 2^d - 1}` (given
with its inverse) that preserves adjacency in both directions (`IsAutomorphism`).

* `automorphism_isometry`: graph automorphisms preserve Hamming distance (the graph
  distance of `Q_d`): geodesics are followed one bit at a time.
* `hist_invariant`: an isometry fixing the landmark set `S` setwise preserves every
  edge histogram, `H_{σe} = H_e`.
* `trivial_stabilizer` (**L1**): for `d ≥ 2`, if `S` is edge-multiset resolving and the
  automorphism `σ` maps `S` into itself, then `σ` is the identity.
-/

namespace ProcgenSelfieEdim

/-- `σ` maps vertices of `Q_d` to vertices and preserves Hamming distance. -/
def IsIsometry (d : Nat) (σ : Nat → Nat) : Prop :=
  (∀ u, u < 2 ^ d → σ u < 2 ^ d) ∧
    ∀ u v, u < 2 ^ d → v < 2 ^ d → hamming d (σ u) (σ v) = hamming d u v

/-- `σ` is a graph automorphism of `Q_d`: a bijection of `{0, …, 2^d - 1}` with inverse
`τ`, preserving adjacency (Hamming distance one) in both directions. -/
def IsAutomorphism (d : Nat) (σ : Nat → Nat) : Prop :=
  ∃ τ : Nat → Nat, (∀ u, u < 2 ^ d → σ u < 2 ^ d) ∧ (∀ u, u < 2 ^ d → τ u < 2 ^ d) ∧
    (∀ u, u < 2 ^ d → τ (σ u) = u) ∧ (∀ u, u < 2 ^ d → σ (τ u) = u) ∧
    ∀ u v, u < 2 ^ d → v < 2 ^ d → (hamming d u v = 1 ↔ hamming d (σ u) (σ v) = 1)

/-! ## Geodesics and isometries -/

/-- One step along a geodesic: if `d(u, v) = k + 1` there is a neighbour `w` of `u`
with `d(w, v) = k`. -/
theorem geodesic_step : ∀ (d u v k : Nat), u < 2 ^ d → v < 2 ^ d → hamming d u v = k + 1 →
    ∃ w, w < 2 ^ d ∧ hamming d u w = 1 ∧ hamming d w v = k
  | 0, u, v, _, _, _, h => by
    exfalso
    have : hamming 0 u v = 0 := rfl
    omega
  | d + 1, u, v, k, hu, hv, h => by
    have hp := two_pow_succ d
    rw [hamming_succ] at h
    by_cases hm : u % 2 = v % 2
    · rw [if_pos hm] at h
      obtain ⟨w, hw, h1, h2⟩ :=
        geodesic_step d (u / 2) (v / 2) k (by omega) (by omega) (by omega)
      refine ⟨2 * w + u % 2, by omega, ?_, ?_⟩
      · rw [hamming_succ]
        have e1 : (2 * w + u % 2) % 2 = u % 2 := by omega
        have e2 : (2 * w + u % 2) / 2 = w := by omega
        rw [e1, e2, if_pos rfl, h1]
      · rw [hamming_succ]
        have e1 : (2 * w + u % 2) % 2 = u % 2 := by omega
        have e2 : (2 * w + u % 2) / 2 = w := by omega
        rw [e1, e2, if_pos hm, h2, Nat.zero_add]
    · rw [if_neg hm] at h
      refine ⟨2 * (u / 2) + v % 2, by omega, ?_, ?_⟩
      · rw [hamming_succ]
        have e1 : (2 * (u / 2) + v % 2) % 2 = v % 2 := by omega
        have e2 : (2 * (u / 2) + v % 2) / 2 = u / 2 := by omega
        rw [e1, e2, if_neg hm, hamming_self]
      · rw [hamming_succ]
        have e1 : (2 * (u / 2) + v % 2) % 2 = v % 2 := by omega
        have e2 : (2 * (u / 2) + v % 2) / 2 = u / 2 := by omega
        rw [e1, e2, if_pos rfl]
        omega

/-- An adjacency-preserving map does not increase Hamming distance. -/
theorem hamming_map_le (d : Nat) (σ : Nat → Nat) (hmap : ∀ u, u < 2 ^ d → σ u < 2 ^ d)
    (hadj : ∀ u v, u < 2 ^ d → v < 2 ^ d → hamming d u v = 1 → hamming d (σ u) (σ v) = 1) :
    ∀ (k u v : Nat), u < 2 ^ d → v < 2 ^ d → hamming d u v = k → hamming d (σ u) (σ v) ≤ k
  | 0, u, v, hu, hv, h => by
    rw [eq_of_hamming_eq_zero d u v hu hv h, hamming_self]
    exact Nat.le_refl 0
  | k + 1, u, v, hu, hv, h => by
    obtain ⟨w, hw, h1, h2⟩ := geodesic_step d u v k hu hv h
    have ih := hamming_map_le d σ hmap hadj k w v hw hv h2
    have hstep := hadj u w hu hw h1
    have htri := hamming_triangle d (σ u) (σ w) (σ v)
    omega

/-- **Graph automorphisms of `Q_d` are Hamming isometries.** -/
theorem automorphism_isometry (d : Nat) (σ : Nat → Nat) (h : IsAutomorphism d σ) :
    IsIsometry d σ := by
  obtain ⟨τ, hσ, hτ, hτσ, hστ, hadj⟩ := h
  refine ⟨hσ, fun u v hu hv => Nat.le_antisymm ?_ ?_⟩
  · exact hamming_map_le d σ hσ (fun x y hx hy hxy => (hadj x y hx hy).1 hxy) _ u v hu hv rfl
  · have hτadj : ∀ x y, x < 2 ^ d → y < 2 ^ d → hamming d x y = 1 → hamming d (τ x) (τ y) = 1 := by
      intro x y hx hy hxy
      apply (hadj (τ x) (τ y) (hτ x hx) (hτ y hy)).2
      rw [hστ x hx, hστ y hy]
      exact hxy
    have := hamming_map_le d τ hτ hτadj _ (σ u) (σ v) (hσ u hu) (hσ v hv) rfl
    rw [hτσ u hu, hτσ v hv] at this
    exact this

/-! ## Histograms are invariant under the stabilizer -/

/-- An injective self-map of a duplicate-free list is onto it. -/
theorem surj_of_inj (f : Nat → Nat) (l : List Nat) (hnd : l.Nodup)
    (hinj : ∀ a ∈ l, ∀ b ∈ l, f a = f b → a = b) (hmaps : ∀ a ∈ l, f a ∈ l) :
    ∀ y ∈ l, ∃ x ∈ l, f x = y := by
  intro y hy
  by_cases hm : y ∈ l.map f
  · rw [List.mem_map] at hm
    obtain ⟨x, hx, hfx⟩ := hm
    exact ⟨x, hx, hfx⟩
  · exfalso
    have h1 := pigeonhole (l.map f) (l.filter (fun z => decide (z ≠ y)))
      (nodup_map_of_inj_on f l hnd hinj) (fun z hz => by
        rw [List.mem_filter]
        refine ⟨?_, decide_eq_true (fun hzy => hm (hzy ▸ hz))⟩
        rw [List.mem_map] at hz
        obtain ⟨x, hx, rfl⟩ := hz
        exact hmaps x hx)
    have h2 := length_filter_lt (fun z => decide (z ≠ y)) y l hy (by simp)
    rw [List.length_map] at h1
    omega

/-- Counting through a bijection of `S`: `#{s ∈ S : P s} = #{t ∈ S : P (σ t)}`. -/
theorem length_filter_comp (σ : Nat → Nat) (S : List Nat) (hnd : S.Nodup)
    (hinj : ∀ a ∈ S, ∀ b ∈ S, σ a = σ b → a = b) (hmaps : ∀ a ∈ S, σ a ∈ S) (P : Nat → Bool) :
    (S.filter P).length = (S.filter (fun t => P (σ t))).length := by
  have hnd1 : (S.filter P).Nodup := List.Nodup.sublist List.filter_sublist hnd
  have hnd2 : (S.filter (fun t => P (σ t))).Nodup := List.Nodup.sublist List.filter_sublist hnd
  apply Nat.le_antisymm
  · -- every `s ∈ S` with `P s` is `σ t` with `t ∈ S`, `P (σ t)`
    have := pigeonhole (S.filter P) ((S.filter (fun t => P (σ t))).map σ) hnd1 (fun s hs => by
      rw [List.mem_filter] at hs
      obtain ⟨t, ht, rfl⟩ := surj_of_inj σ S hnd hinj hmaps s hs.1
      exact List.mem_map_of_mem (List.mem_filter.2 ⟨ht, hs.2⟩))
    rw [List.length_map] at this
    exact this
  · have := pigeonhole ((S.filter (fun t => P (σ t))).map σ) (S.filter P)
      (nodup_map_of_inj_on σ _ hnd2 (fun a ha b hb hab =>
        hinj a (List.mem_filter.1 ha).1 b (List.mem_filter.1 hb).1 hab))
      (fun s hs => by
        rw [List.mem_map] at hs
        obtain ⟨t, ht, rfl⟩ := hs
        rw [List.mem_filter] at ht ⊢
        exact ⟨hmaps t ht.1, ht.2⟩)
    rw [List.length_map] at this
    exact this

/-- An isometry mapping `S` into itself preserves the histogram of every edge. -/
theorem hist_invariant (d : Nat) (S : List Nat) (hS : ∀ s ∈ S, s < 2 ^ d) (hnd : S.Nodup)
    (σ : Nat → Nat) (hσ : IsIsometry d σ) (hfix : ∀ s ∈ S, σ s ∈ S) (u v : Nat)
    (hu : u < 2 ^ d) (hv : v < 2 ^ d) (r : Nat) :
    hist d S (σ u) (σ v) r = hist d S u v r := by
  have hinj : ∀ a ∈ S, ∀ b ∈ S, σ a = σ b → a = b := by
    intro a ha b hb hab
    have := hσ.2 a b (hS a ha) (hS b hb)
    rw [hab, hamming_self] at this
    exact eq_of_hamming_eq_zero d a b (hS a ha) (hS b hb) this.symm
  unfold hist
  rw [length_filter_comp σ S hnd hinj hfix]
  congr 1
  apply List.filter_congr
  intro t ht
  unfold edgeDist
  rw [hσ.2 u t hu (hS t ht), hσ.2 v t hv (hS t ht)]

/-! ## L1 -/

/-- Flip bit `i` of `u`. -/
def flipBit (u i : Nat) : Nat := if bit u i = 0 then u + 2 ^ i else u - 2 ^ i

theorem isEdge_flipBit (d u i : Nat) (hu : u < 2 ^ d) (hi : i < d) : IsEdge d u (flipBit u i) := by
  unfold flipBit
  by_cases hb : bit u i = 0
  · rw [if_pos hb]
    exact ⟨hu, flip_lt d i u hi hu hb, hamming_flip d i u hi hb⟩
  · rw [if_neg hb]
    have hpos := Nat.two_pow_pos i
    have hge : 2 ^ i ≤ u := by
      apply Nat.le_of_not_lt
      intro hlt
      apply hb
      unfold bit
      rw [Nat.div_eq_of_lt hlt]
    have hb' : bit (u - 2 ^ i) i = 0 := by
      unfold bit at hb ⊢
      have := Nat.sub_mul_div u (2 ^ i) 1
      rw [Nat.mul_one] at this
      rw [this]
      omega
    have h1 := hamming_flip d i (u - 2 ^ i) hi hb'
    rw [show u - 2 ^ i + 2 ^ i = u by omega, hamming_comm] at h1
    exact ⟨hu, by omega, h1⟩

theorem flipBit_dist (u i : Nat) : flipBit u i + 2 ^ i = u ∨ u + 2 ^ i = flipBit u i := by
  unfold flipBit
  by_cases hb : bit u i = 0
  · rw [if_pos hb]
    exact Or.inr rfl
  · rw [if_neg hb]
    have hge : 2 ^ i ≤ u := by
      apply Nat.le_of_not_lt
      intro hlt
      apply hb
      unfold bit
      rw [Nat.div_eq_of_lt hlt]
    exact Or.inl (by omega)

/-- **L1 (trivial stabilizer).** Let `d ≥ 2`, let `S` be a set of vertices of `Q_d` that
is edge-multiset resolving, and let `σ` be a graph automorphism of `Q_d` with
`σ(S) ⊆ S`. Then `σ` fixes every vertex. -/
theorem trivial_stabilizer (d : Nat) (hd : 2 ≤ d) (S : List Nat) (hS : ∀ s ∈ S, s < 2 ^ d)
    (hnd : S.Nodup) (hres : Resolving d S) (σ : Nat → Nat) (hσ : IsAutomorphism d σ)
    (hfix : ∀ s ∈ S, σ s ∈ S) : ∀ u, u < 2 ^ d → σ u = u := by
  have hiso := automorphism_isometry d σ hσ
  -- every edge is mapped to itself (as an unordered pair)
  have hedge : ∀ u v, IsEdge d u v → SameEdge (σ u) (σ v) u v := by
    intro u v he
    obtain ⟨hu, hv, h1⟩ := he
    have he' : IsEdge d (σ u) (σ v) :=
      ⟨hiso.1 u hu, hiso.1 v hv, by rw [hiso.2 u v hu hv]; exact h1⟩
    exact hres (σ u) (σ v) u v he' ⟨hu, hv, h1⟩
      (fun r => hist_invariant d S hS hnd σ hiso hfix u v hu hv r)
  intro u hu
  have e0 := hedge u (flipBit u 0) (isEdge_flipBit d u 0 hu (by omega))
  have e1 := hedge u (flipBit u 1) (isEdge_flipBit d u 1 hu (by omega))
  have d0 := flipBit_dist u 0
  have d1 := flipBit_dist u 1
  have p0 : (2 : Nat) ^ 0 = 1 := rfl
  have p1 : (2 : Nat) ^ 1 = 2 := rfl
  rw [p0] at d0
  rw [p1] at d1
  rcases e0 with ⟨h, _⟩ | ⟨h, _⟩
  · exact h
  · rcases e1 with ⟨h', _⟩ | ⟨h', _⟩
    · exact h'
    · exfalso
      rw [h] at h'
      omega

/-- The paper's resolving 15-set of `Q_6` has trivial stabilizer in `Aut(Q_6)`. -/
theorem paperSet_trivial_stabilizer (σ : Nat → Nat) (hσ : IsAutomorphism 6 σ)
    (hfix : ∀ s ∈ paperSet, σ s ∈ paperSet) : ∀ u, u < 2 ^ 6 → σ u = u :=
  trivial_stabilizer 6 (by decide) paperSet paperSet_is_vertex_set.2.2.1 paperSet_is_vertex_set.2.1
    paperSet_resolving σ hσ hfix

end ProcgenSelfieEdim
