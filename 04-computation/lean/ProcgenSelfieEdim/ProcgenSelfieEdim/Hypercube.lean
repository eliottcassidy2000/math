import ProcgenSelfieEdim.ListLemmas

set_option autoImplicit false

/-!
# The hypercube `Q_d` and edge-multiset resolving sets

Vertices of `Q_d` are the natural numbers `0, …, 2^d - 1`, read as binary words
(bit `i` of `v` is coordinate `i`). Two vertices are adjacent iff they differ in
exactly one bit, so the graph distance of `Q_d` is the Hamming distance
`hamming d u v` (number of differing bits among bits `0, …, d-1`).

Following Allikvere (arXiv:2608.09983), for an edge `e = uv` and a vertex `s`,
`d(e, s) = min (d(u, s), d(v, s))`, and the histogram of `e` with respect to a
landmark list `S` is `H_e(r) = #{s ∈ S : d(e, s) = r}`. The set `S` is
edge-multiset resolving if distinct edges have distinct histograms.
-/

namespace ProcgenSelfieEdim

/-! ## Definitions -/

/-- Hamming distance between the `d`-bit words `u` and `v` (bits `0, …, d-1`). -/
def hamming : Nat → Nat → Nat → Nat
  | 0, _, _ => 0
  | d + 1, u, v => (if u % 2 = v % 2 then 0 else 1) + hamming d (u / 2) (v / 2)

/-- `uv` is an edge of `Q_d`: both ends are vertices and they differ in one bit. -/
def IsEdge (d u v : Nat) : Prop := u < 2 ^ d ∧ v < 2 ^ d ∧ hamming d u v = 1

/-- The edge-vertex distance `d(uv, s) = min (d(u, s), d(v, s))`. -/
def edgeDist (d u v s : Nat) : Nat := min (hamming d u s) (hamming d v s)

/-- The histogram `H_{uv}(r) = #{s ∈ S : d(uv, s) = r}`. -/
def hist (d : Nat) (S : List Nat) (u v r : Nat) : Nat :=
  (S.filter (fun s => edgeDist d u v s == r)).length

/-- `{u, v} = {u', v'}` as unordered pairs. -/
def SameEdge (u v u' v' : Nat) : Prop := (u = u' ∧ v = v') ∨ (u = v' ∧ v = u')

/-- `S` is edge-multiset resolving in `Q_d`: two edges with the same histogram
(the same multiset of edge-vertex distances) are the same edge. -/
def Resolving (d : Nat) (S : List Nat) : Prop :=
  ∀ u v u' v', IsEdge d u v → IsEdge d u' v' →
    (∀ r, hist d S u v r = hist d S u' v' r) → SameEdge u v u' v'

/-- The antipode (bitwise complement) of a vertex of `Q_d`. -/
def antipode (d u : Nat) : Nat := 2 ^ d - 1 - u

/-- The `i`-th bit of `u`. -/
def bit (u i : Nat) : Nat := u / 2 ^ i % 2

/-! ## Hamming distance -/

theorem hamming_succ (d u v : Nat) :
    hamming (d + 1) u v = (if u % 2 = v % 2 then 0 else 1) + hamming d (u / 2) (v / 2) := rfl

theorem hamming_le : ∀ (d u v : Nat), hamming d u v ≤ d
  | 0, _, _ => Nat.le_refl 0
  | d + 1, u, v => by
    have := hamming_le d (u / 2) (v / 2)
    rw [hamming_succ]
    split <;> omega

theorem hamming_comm : ∀ (d u v : Nat), hamming d u v = hamming d v u
  | 0, _, _ => rfl
  | d + 1, u, v => by
    rw [hamming_succ, hamming_succ, hamming_comm d (u / 2) (v / 2)]
    by_cases h : u % 2 = v % 2
    · rw [if_pos h, if_pos h.symm]
    · rw [if_neg h, if_neg (Ne.symm h)]

theorem hamming_self : ∀ (d u : Nat), hamming d u u = 0
  | 0, _ => rfl
  | d + 1, u => by rw [hamming_succ, hamming_self d (u / 2), if_pos rfl]

theorem two_pow_succ (d : Nat) : 2 ^ (d + 1) = 2 ^ d * 2 := Nat.pow_succ 2 d

/-- Vertices at Hamming distance zero coincide. -/
theorem eq_of_hamming_eq_zero :
    ∀ (d u v : Nat), u < 2 ^ d → v < 2 ^ d → hamming d u v = 0 → u = v
  | 0, u, v, hu, hv, _ => by
    simp only [Nat.pow_zero] at hu hv
    omega
  | d + 1, u, v, hu, hv, h => by
    have hp := two_pow_succ d
    rw [hamming_succ] at h
    by_cases hm : u % 2 = v % 2
    · rw [if_pos hm] at h
      have := eq_of_hamming_eq_zero d (u / 2) (v / 2) (by omega) (by omega) (by omega)
      omega
    · rw [if_neg hm] at h
      omega

/-- Vertices at Hamming distance zero are at the same distance from every vertex. -/
theorem hamming_congr_of_zero :
    ∀ (d a b s : Nat), hamming d a b = 0 → hamming d a s = hamming d b s
  | 0, _, _, _, _ => rfl
  | d + 1, a, b, s, h => by
    rw [hamming_succ] at h
    rw [hamming_succ, hamming_succ]
    by_cases hm : a % 2 = b % 2
    · rw [if_pos hm] at h
      rw [hamming_congr_of_zero d (a / 2) (b / 2) (s / 2) (by omega), hm]
    · rw [if_neg hm] at h
      omega

/-- The two ends of an edge are at distances differing by exactly one from every word. -/
theorem hamming_adjacent : ∀ (d u v s : Nat), hamming d u v = 1 →
    hamming d u s + 1 = hamming d v s ∨ hamming d v s + 1 = hamming d u s
  | 0, _, _, _, h => by simp [hamming] at h
  | d + 1, u, v, s, h => by
    rw [hamming_succ] at h
    rw [hamming_succ, hamming_succ]
    by_cases hm : u % 2 = v % 2
    · rw [if_pos hm] at h
      have ih := hamming_adjacent d (u / 2) (v / 2) (s / 2) (by omega)
      rw [hm]
      omega
    · rw [if_neg hm] at h
      rw [hamming_congr_of_zero d (u / 2) (v / 2) (s / 2) (by omega)]
      by_cases hu : u % 2 = s % 2
      · rw [if_pos hu, if_neg (by omega)]
        omega
      · rw [if_neg hu, if_pos (by omega)]
        omega

/-- The triangle inequality for the Hamming distance. -/
theorem hamming_triangle : ∀ (d u v w : Nat), hamming d u w ≤ hamming d u v + hamming d v w
  | 0, _, _, _ => Nat.le_refl 0
  | d + 1, u, v, w => by
    have ih := hamming_triangle d (u / 2) (v / 2) (w / 2)
    rw [hamming_succ, hamming_succ, hamming_succ]
    by_cases h1 : u % 2 = v % 2 <;> by_cases h2 : v % 2 = w % 2 <;>
      by_cases h3 : u % 2 = w % 2 <;> simp only [h1, h2, h3, if_true, if_false] <;> omega

/-! ## Antipodes -/

theorem antipode_lt (d u : Nat) (hu : u < 2 ^ d) : antipode d u < 2 ^ d := by
  unfold antipode
  omega

theorem antipode_antipode (d u : Nat) (hu : u < 2 ^ d) : antipode d (antipode d u) = u := by
  unfold antipode
  omega

/-- The antipode is at distance `d - h` from a word at distance `h`. -/
theorem hamming_antipode :
    ∀ (d u s : Nat), u < 2 ^ d → hamming d (antipode d u) s = d - hamming d u s
  | 0, _, _, _ => rfl
  | d + 1, u, s, hu => by
    have hp := two_pow_succ d
    have hpos : 0 < 2 ^ d := Nat.two_pow_pos d
    have h1 : antipode (d + 1) u % 2 = 1 - u % 2 := by
      unfold antipode
      omega
    have h2 : antipode (d + 1) u / 2 = antipode d (u / 2) := by
      unfold antipode
      omega
    have ih := hamming_antipode d (u / 2) (s / 2) (by omega)
    have hle := hamming_le d (u / 2) (s / 2)
    rw [hamming_succ, hamming_succ, h1, h2, ih]
    by_cases hm : u % 2 = s % 2
    · rw [if_pos hm, if_neg (by omega)]
      omega
    · rw [if_neg hm, if_pos (by omega)]
      omega

theorem hamming_antipode_antipode (d u v : Nat) (hu : u < 2 ^ d) (hv : v < 2 ^ d) :
    hamming d (antipode d u) (antipode d v) = hamming d u v := by
  rw [hamming_antipode d u _ hu, hamming_comm d u, hamming_antipode d v u hv]
  have := hamming_le d v u
  rw [hamming_comm d v u] at this ⊢
  omega

/-- The antipodal edge of an edge is an edge. -/
theorem isEdge_antipode (d u v : Nat) (he : IsEdge d u v) :
    IsEdge d (antipode d u) (antipode d v) := by
  obtain ⟨hu, hv, h1⟩ := he
  exact ⟨antipode_lt d u hu, antipode_lt d v hv, by rw [hamming_antipode_antipode d u v hu hv, h1]⟩

/-! ## Edge-vertex distances and histograms -/

theorem edgeDist_comm (d u v s : Nat) : edgeDist d u v s = edgeDist d v u s := by
  unfold edgeDist
  omega

theorem hist_comm (d : Nat) (S : List Nat) (u v r : Nat) : hist d S u v r = hist d S v u r := by
  unfold hist
  rw [List.filter_congr (fun s _ => by rw [edgeDist_comm])]

/-- Every edge-vertex distance is at most `d - 1`. -/
theorem edgeDist_le (d u v s : Nat) (h : hamming d u v = 1) : edgeDist d u v s + 1 ≤ d := by
  have h1 := hamming_le d u s
  have h2 := hamming_le d v s
  have := hamming_adjacent d u v s h
  unfold edgeDist
  omega

/-- Antipodal reversal of edge-vertex distances: `d(ē, s) = d - 1 - d(e, s)`. -/
theorem edgeDist_antipode (d u v s : Nat) (hu : u < 2 ^ d) (hv : v < 2 ^ d)
    (h : hamming d u v = 1) :
    edgeDist d (antipode d u) (antipode d v) s = d - 1 - edgeDist d u v s := by
  have h1 := hamming_le d u s
  have h2 := hamming_le d v s
  have := hamming_adjacent d u v s h
  unfold edgeDist
  rw [hamming_antipode d u s hu, hamming_antipode d v s hv]
  omega

/-- **L2 (antipodal reversal).** For every edge `e` of `Q_d`, every landmark list `S`
and every level `0 ≤ r ≤ d - 1`, `H_{ē}(r) = H_e(d - 1 - r)`. -/
theorem hist_antipode (d : Nat) (S : List Nat) (u v : Nat) (he : IsEdge d u v)
    (r : Nat) (hr : r < d) :
    hist d S (antipode d u) (antipode d v) r = hist d S u v (d - 1 - r) := by
  obtain ⟨hu, hv, h1⟩ := he
  unfold hist
  congr 1
  apply List.filter_congr
  intro s _
  rw [edgeDist_antipode d u v s hu hv h1]
  have := edgeDist_le d u v s h1
  rw [Bool.eq_iff_iff, beq_iff_eq, beq_iff_eq]
  constructor <;> intro <;> omega

/-! ## Bits and the canonical list of edges -/

theorem bit_zero (u : Nat) : bit u 0 = u % 2 := by
  simp [bit]

theorem bit_succ (u i : Nat) : bit u (i + 1) = bit (u / 2) i := by
  unfold bit
  rw [Nat.pow_succ, Nat.mul_comm, ← Nat.div_div_eq_div_mul]

/-- Flipping a zero bit moves to an adjacent vertex. -/
theorem hamming_flip : ∀ (d i a : Nat), i < d → bit a i = 0 → hamming d a (a + 2 ^ i) = 1
  | 0, _, _, hi, _ => absurd hi (Nat.not_lt_zero _)
  | d + 1, 0, a, _, hb => by
    rw [bit_zero] at hb
    rw [hamming_succ, Nat.pow_zero]
    have h1 : (a + 1) / 2 = a / 2 := by omega
    rw [h1, hamming_self, if_neg (by omega)]
  | d + 1, i + 1, a, hi, hb => by
    rw [bit_succ] at hb
    rw [hamming_succ]
    have hp := two_pow_succ i
    have h1 : (a + 2 ^ (i + 1)) % 2 = a % 2 := by omega
    have h2 : (a + 2 ^ (i + 1)) / 2 = a / 2 + 2 ^ i := by omega
    rw [h1, h2, if_pos rfl, hamming_flip d i (a / 2) (by omega) hb]

/-- Flipping a zero bit stays inside the vertex set. -/
theorem flip_lt : ∀ (d i a : Nat), i < d → a < 2 ^ d → bit a i = 0 → a + 2 ^ i < 2 ^ d
  | 0, _, _, hi, _, _ => absurd hi (Nat.not_lt_zero _)
  | d + 1, 0, a, _, ha, hb => by
    rw [bit_zero] at hb
    have hp := two_pow_succ d
    rw [Nat.pow_zero]
    omega
  | d + 1, i + 1, a, hi, ha, hb => by
    rw [bit_succ] at hb
    have hp := two_pow_succ d
    have hq := two_pow_succ i
    have := flip_lt d i (a / 2) (by omega) (by omega) hb
    omega

/-- Every edge flips one bit, from `0` at the lower end to `1` at the upper end. -/
theorem edge_cases : ∀ (d u v : Nat), u < 2 ^ d → v < 2 ^ d → hamming d u v = 1 →
    ∃ i, i < d ∧ ((bit u i = 0 ∧ v = u + 2 ^ i) ∨ (bit v i = 0 ∧ u = v + 2 ^ i))
  | 0, _, _, _, _, h => by simp [hamming] at h
  | d + 1, u, v, hu, hv, h => by
    have hp := two_pow_succ d
    rw [hamming_succ] at h
    by_cases hm : u % 2 = v % 2
    · rw [if_pos hm] at h
      obtain ⟨i, hi, hc⟩ := edge_cases d (u / 2) (v / 2) (by omega) (by omega) (by omega)
      have hq := two_pow_succ i
      refine ⟨i + 1, by omega, ?_⟩
      rw [bit_succ, bit_succ]
      rcases hc with ⟨hb, he⟩ | ⟨hb, he⟩
      · exact Or.inl ⟨hb, by omega⟩
      · exact Or.inr ⟨hb, by omega⟩
    · rw [if_neg hm] at h
      have heq := eq_of_hamming_eq_zero d (u / 2) (v / 2) (by omega) (by omega) (by omega)
      refine ⟨0, by omega, ?_⟩
      rw [bit_zero, bit_zero, Nat.pow_zero]
      by_cases hu0 : u % 2 = 0
      · exact Or.inl ⟨hu0, by omega⟩
      · exact Or.inr ⟨by omega, by omega⟩

/-- The canonical list of edges `(u, u + 2^i)` with bit `i` of `u` equal to `0`. -/
def edgeList (d : Nat) : List (Nat × Nat) :=
  (List.range (2 ^ d)).flatMap (fun u =>
    (List.range d).filterMap (fun i => if bit u i = 0 then some (u, u + 2 ^ i) else none))

theorem mem_edgeList (d a b : Nat) :
    (a, b) ∈ edgeList d ↔ a < 2 ^ d ∧ ∃ i, i < d ∧ bit a i = 0 ∧ b = a + 2 ^ i := by
  unfold edgeList
  rw [List.mem_flatMap]
  constructor
  · rintro ⟨u, hu, hmem⟩
    rw [List.mem_filterMap] at hmem
    obtain ⟨i, hi, hsome⟩ := hmem
    rw [List.mem_range] at hu hi
    by_cases hb : bit u i = 0
    · rw [if_pos hb, Option.some.injEq, Prod.mk.injEq] at hsome
      obtain ⟨rfl, rfl⟩ := hsome
      exact ⟨hu, i, hi, hb, rfl⟩
    · rw [if_neg hb] at hsome
      cases hsome
  · rintro ⟨ha, i, hi, hb, rfl⟩
    refine ⟨a, List.mem_range.2 ha, ?_⟩
    rw [List.mem_filterMap]
    exact ⟨i, List.mem_range.2 hi, by rw [if_pos hb]⟩

theorem isEdge_of_mem_edgeList (d a b : Nat) (h : (a, b) ∈ edgeList d) :
    IsEdge d a b ∧ a < b := by
  rw [mem_edgeList] at h
  obtain ⟨ha, i, hi, hb, rfl⟩ := h
  have hpos := Nat.two_pow_pos i
  exact ⟨⟨ha, flip_lt d i a hi ha hb, hamming_flip d i a hi hb⟩, by omega⟩

/-- Every edge appears in `edgeList d`, in one of its two orientations. -/
theorem canonical_mem (d u v : Nat) (he : IsEdge d u v) :
    (u, v) ∈ edgeList d ∨ (v, u) ∈ edgeList d := by
  obtain ⟨hu, hv, h1⟩ := he
  obtain ⟨i, hi, hc⟩ := edge_cases d u v hu hv h1
  rcases hc with ⟨hb, he⟩ | ⟨hb, he⟩
  · exact Or.inl ((mem_edgeList d u v).2 ⟨hu, i, hi, hb, he⟩)
  · exact Or.inr ((mem_edgeList d v u).2 ⟨hv, i, hi, hb, he⟩)

/-! ## Histogram sums -/

theorem hist_nil (d u v r : Nat) : hist d [] u v r = 0 := rfl

theorem hist_cons (d s : Nat) (t : List Nat) (u v r : Nat) :
    hist d (s :: t) u v r = (if edgeDist d u v s = r then 1 else 0) + hist d t u v r := by
  unfold hist
  by_cases h : edgeDist d u v s = r
  · rw [List.filter_cons_of_pos (p := fun s => edgeDist d u v s == r) (beq_iff_eq.2 h), if_pos h,
      List.length_cons]
    omega
  · rw [List.filter_cons_of_neg (p := fun s => edgeDist d u v s == r) (fun hc => h (beq_iff_eq.1 hc)),
      if_neg h]
    omega

/-- Histograms vanish above level `d - 1`. -/
theorem hist_eq_zero_of_le (d : Nat) (u v r : Nat) (h : hamming d u v = 1) (hr : d ≤ r) :
    ∀ (S : List Nat), hist d S u v r = 0
  | [] => rfl
  | s :: t => by
    have := edgeDist_le d u v s h
    rw [hist_cons, if_neg (by omega), hist_eq_zero_of_le d u v r h hr t]

/-- A sum over the landmarks of a function of the edge-vertex distance is a
function of the histogram: `Σ_{s ∈ S} g(d(e,s)) = Σ_{r < d} g(r) · H_e(r)`. -/
theorem weighted_sum_eq (d : Nat) (u v : Nat) (h : hamming d u v = 1) (g : Nat → Nat) :
    ∀ (S : List Nat), (S.map (fun s => g (edgeDist d u v s))).sum =
      rsum d (fun r => g r * hist d S u v r)
  | [] => by
    rw [List.map_nil, List.sum_nil, rsum_congr d _ (fun _ => 0) (fun r _ => by rw [hist_nil, Nat.mul_zero]),
      rsum_zero]
  | s :: t => by
    rw [List.map_cons, List.sum_cons, weighted_sum_eq d u v h g t,
      rsum_congr d (fun r => g r * hist d (s :: t) u v r)
        (fun r => (if edgeDist d u v s = r then g r else 0) + g r * hist d t u v r)
        (fun r _ => by
          dsimp only
          rw [hist_cons, Nat.mul_add]
          by_cases hr : edgeDist d u v s = r
          · rw [if_pos hr, if_pos hr, Nat.mul_one]
          · rw [if_neg hr, if_neg hr, Nat.mul_zero]),
      rsum_add, rsum_single d (edgeDist d u v s) g (by have := edgeDist_le d u v s h; omega)]

theorem sum_map_one (S : List Nat) : (S.map (fun _ => 1)).sum = S.length := by
  induction S with
  | nil => rfl
  | cons s t ih => rw [List.map_cons, List.sum_cons, ih, List.length_cons]; omega

/-- The histogram of an edge has total mass `|S|`. -/
theorem rsum_hist (d : Nat) (S : List Nat) (u v : Nat) (h : hamming d u v = 1) :
    rsum d (fun r => hist d S u v r) = S.length := by
  rw [← sum_map_one S]
  have := weighted_sum_eq d u v h (fun _ => 1) S
  simp only [Nat.one_mul] at this
  exact this.symm

/-! ## Generic soundness of key-based checks -/

/-- If `key` is a function of the histogram on edges, and the keys of the canonical
edges are pairwise distinct, then `S` is resolving. -/
theorem resolving_of_keys {α : Type} (d : Nat) (S : List Nat) (key : Nat × Nat → α)
    (hkey : ∀ a b a' b', IsEdge d a b → IsEdge d a' b' →
      (∀ r, hist d S a b r = hist d S a' b' r) → key (a, b) = key (a', b'))
    (hnd : ((edgeList d).map key).Nodup) : Resolving d S := by
  intro u v u' v' he he' hh
  have hcan : ∀ a b a' b', (a, b) ∈ edgeList d → (a', b') ∈ edgeList d →
      (∀ r, hist d S a b r = hist d S a' b' r) → a = a' ∧ b = b' := by
    intro a b a' b' hm hm' hab
    have := eq_of_map_nodup key (edgeList d) hnd (a, b) (a', b') hm hm'
      (hkey a b a' b' (isEdge_of_mem_edgeList d a b hm).1 (isEdge_of_mem_edgeList d a' b' hm').1 hab)
    rw [Prod.mk.injEq] at this
    exact this
  rcases canonical_mem d u v he with h1 | h1 <;> rcases canonical_mem d u' v' he' with h2 | h2
  · exact Or.inl (hcan u v u' v' h1 h2 hh)
  · have := hcan u v v' u' h1 h2 (fun r => by rw [hh r, hist_comm])
    exact Or.inr ⟨this.1, this.2⟩
  · have := hcan v u u' v' h1 h2 (fun r => by rw [hist_comm, hh r])
    exact Or.inr ⟨this.2, this.1⟩
  · have := hcan v u v' u' h1 h2 (fun r => by rw [hist_comm, hh r, hist_comm])
    exact Or.inl ⟨this.2, this.1⟩

/-! ## A kernel-friendly checker

Definitions by structural recursion compile to `brecOn`, which allocates a tuple at
each step of kernel reduction. The checker is therefore written with the recursors
`Nat.rec` and `List.rec`, and each piece is proved equal to its readable
counterpart. -/

/-- `hamming`, written with `Nat.rec`. -/
noncomputable def hamR (d : Nat) : Nat → Nat → Nat :=
  Nat.rec (motive := fun _ => Nat → Nat → Nat) (fun _ _ => 0)
    (fun _ ih u v => (u + v) % 2 + ih (u / 2) (v / 2)) d

theorem hamR_eq : ∀ (d u v : Nat), hamR d u v = hamming d u v
  | 0, _, _ => rfl
  | d + 1, u, v => by
    show (u + v) % 2 + hamR d (u / 2) (v / 2) = hamming (d + 1) u v
    rw [hamR_eq d, hamming_succ]
    split <;> omega

/-- `min`, through `Nat.ble`. -/
def fmin (a b : Nat) : Nat := cond (Nat.ble a b) a b

theorem fmin_eq (a b : Nat) : fmin a b = min a b := by
  unfold fmin
  cases h : Nat.ble a b
  · have : ¬ a ≤ b := fun hab => by rw [Nat.ble_eq_true_of_le hab] at h; cases h
    show b = min a b
    omega
  · have := Nat.le_of_ble_eq_true h
    show a = min a b
    omega

/-- The one-pass key `Σ_{s ∈ S} 16^{d(uv, s)}`. -/
noncomputable def keyR (d : Nat) (S : List Nat) (u v : Nat) : Nat :=
  List.rec (motive := fun _ => Nat) 0
    (fun s _ ih => 16 ^ fmin (hamR d u s) (hamR d v s) + ih) S

theorem keyR_eq (d : Nat) (u v : Nat) :
    ∀ (S : List Nat), keyR d S u v = (S.map (fun s => 16 ^ edgeDist d u v s)).sum
  | [] => rfl
  | s :: t => by
    show 16 ^ fmin (hamR d u s) (hamR d v s) + keyR d t u v = _
    rw [keyR_eq d u v t, List.map_cons, List.sum_cons, fmin_eq, hamR_eq, hamR_eq]
    rfl

/-- The checker: the one-pass keys of the canonical edges are pairwise distinct. -/
noncomputable def checkResolving (d : Nat) (S : List Nat) : Bool :=
  distinctR ((edgeList d).map (fun e => keyR d S e.1 e.2))

/-- Soundness of the checker. -/
theorem resolving_of_check (d : Nat) (S : List Nat) (h : checkResolving d S = true) :
    Resolving d S := by
  refine resolving_of_keys d S (fun e => keyR d S e.1 e.2) ?_ (nodup_of_distinctR _ h)
  intro a b a' b' he he' hh
  show keyR d S a b = keyR d S a' b'
  rw [keyR_eq, keyR_eq, weighted_sum_eq d a b he.2.2 (fun r => 16 ^ r) S,
    weighted_sum_eq d a' b' he'.2.2 (fun r => 16 ^ r) S]
  exact rsum_congr d _ _ (fun r _ => by rw [hh r])

/-! ## L2, parity corollary -/

theorem rsum_shift : ∀ (n : Nat) (f : Nat → Nat), rsum (n + 1) f = f 0 + rsum n (fun r => f (r + 1))
  | 0, f => by
    show rsum 0 f + f 0 = f 0 + rsum 0 (fun r => f (r + 1))
    rw [show rsum 0 f = 0 from rfl, show rsum 0 (fun r => f (r + 1)) = 0 from rfl]
    omega
  | n + 1, f => by
    show rsum (n + 1) f + f (n + 1) = f 0 + (rsum n (fun r => f (r + 1)) + f (n + 1))
    rw [rsum_shift n f]
    omega

/-- A palindromic sequence of even length has an even sum. -/
theorem rsum_palindrome_even : ∀ (k : Nat) (f : Nat → Nat),
    (∀ r, r < 2 * k → f r = f (2 * k - 1 - r)) → rsum (2 * k) f % 2 = 0
  | 0, _, _ => rfl
  | k + 1, f, hpal => by
    have h1 : rsum (2 * (k + 1)) f = rsum (2 * k + 1) f + f (2 * k + 1) := by
      rw [show 2 * (k + 1) = (2 * k + 1) + 1 by omega]
      rfl
    have hg := rsum_palindrome_even k (fun r => f (r + 1)) (fun r hr => by
      dsimp only
      rw [hpal (r + 1) (by omega)]
      congr 1
      omega)
    have h0 : f (2 * k + 1) = f 0 := by
      rw [hpal (2 * k + 1) (by omega)]
      congr 1
      omega
    rw [h1, rsum_shift, h0]
    omega

/-- **L2, corollary.** If `d` is even and `|S|` is odd, no edge histogram is a palindrome. -/
theorem hist_not_palindrome (d : Nat) (hd : d % 2 = 0) (S : List Nat) (hS : S.length % 2 = 1)
    (u v : Nat) (he : IsEdge d u v) :
    ¬ (∀ r, r < d → hist d S u v r = hist d S u v (d - 1 - r)) := by
  intro hpal
  obtain ⟨k, rfl⟩ : ∃ k, d = 2 * k := ⟨d / 2, by omega⟩
  have h1 := rsum_hist (2 * k) S u v he.2.2
  have h2 := rsum_palindrome_even k (fun r => hist (2 * k) S u v r) hpal
  rw [h1] at h2
  omega

/-- **L2, corollary.** An edge collides with its antipodal edge only if its histogram is a
palindrome; so for even `d` and odd `|S|` an edge never collides with its antipode. -/
theorem no_antipodal_collision (d : Nat) (hd : d % 2 = 0) (S : List Nat) (hS : S.length % 2 = 1)
    (u v : Nat) (he : IsEdge d u v) :
    ¬ (∀ r, hist d S u v r = hist d S (antipode d u) (antipode d v) r) := by
  intro hcol
  apply hist_not_palindrome d hd S hS u v he
  intro r hr
  rw [hcol r, hist_antipode d S u v he r hr]

end ProcgenSelfieEdim
