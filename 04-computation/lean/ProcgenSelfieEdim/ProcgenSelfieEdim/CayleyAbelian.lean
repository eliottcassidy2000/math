import ProcgenSelfieEdim.AntiAut

set_option autoImplicit false

/-!
# THM-4529 Theorem 5.1: Cayley digraphs of abelian groups of odd order are arc-even

`FinAbGroup N` is a finite abelian group on the carrier `{0, …, N-1}` (operations on `Nat`
with closure, associativity, commutativity, a neutral element and inverses on the carrier).
For a connection set `S` (any predicate, `0 ∈ S` allowed), the Cayley digraph is
`a → b` iff `b - a ∈ S` (`cayley`).

* `cayley_arcCount_even` (**Theorem 5.1**): if `N` is odd, every arc of `Cay(Γ, S)` lies on an
  even number of Hamiltonian paths. The involution is `P ↦ reverse(φ P)` with
  `φ(x) = u + v - x` (`reflection_arcCount_even`).
* `cyclic n`: `Z/n` as a `FinAbGroup`; `prod G₁ G₂`: the direct product on
  `{0, …, N₁N₂ - 1}` (`x = i N₂ + j`). So `Z/m × Z/n` and every product of cyclic groups is
  covered (`cayley_cyclic_prod_arcCount_even`, e.g. the non-cyclic `Z/3 × Z/3`).
-/

namespace ProcgenSelfieEdim

/-- A finite abelian group on the carrier `{0, …, N-1}`. -/
structure FinAbGroup (N : Nat) where
  add : Nat → Nat → Nat
  neg : Nat → Nat
  zero : Nat
  add_lt : ∀ a b, a < N → b < N → add a b < N
  neg_lt : ∀ a, a < N → neg a < N
  zero_lt : zero < N
  add_assoc : ∀ a b c, a < N → b < N → c < N → add (add a b) c = add a (add b c)
  add_comm : ∀ a b, a < N → b < N → add a b = add b a
  zero_add : ∀ a, a < N → add zero a = a
  neg_add : ∀ a, a < N → add (neg a) a = zero

theorem FinAbGroup.add_zero {N : Nat} (G : FinAbGroup N)
    (a : Nat) (ha : a < N) : G.add a G.zero = a := by
  rw [G.add_comm a G.zero ha G.zero_lt, G.zero_add a ha]

theorem FinAbGroup.add_neg {N : Nat} (G : FinAbGroup N)
    (a : Nat) (ha : a < N) : G.add a (G.neg a) = G.zero := by
  rw [G.add_comm a (G.neg a) ha (G.neg_lt a ha), G.neg_add a ha]

/-- Inverses are unique. -/
theorem FinAbGroup.neg_unique {N : Nat} (G : FinAbGroup N)
    (a b : Nat) (ha : a < N) (hb : b < N) (h : G.add a b = G.zero) :
    b = G.neg a := by
  have e := G.add_assoc (G.neg a) a b (G.neg_lt a ha) ha hb
  rw [G.neg_add a ha, G.zero_add b hb, h, G.add_zero _ (G.neg_lt a ha)] at e
  exact e

/-- `-(s - x) = x - s`. -/
theorem FinAbGroup.neg_sub {N : Nat} (G : FinAbGroup N) (s x : Nat) (hs : s < N) (hx : x < N) :
    G.neg (G.add s (G.neg x)) = G.add x (G.neg s) := by
  have hnx := G.neg_lt x hx
  have hns := G.neg_lt s hs
  symm
  apply G.neg_unique _ _ (G.add_lt _ _ hs hnx) (G.add_lt _ _ hx hns)
  rw [G.add_assoc s (G.neg x) (G.add x (G.neg s)) hs hnx (G.add_lt _ _ hx hns),
    ← G.add_assoc (G.neg x) x (G.neg s) hnx hx hns, G.neg_add x hx, G.zero_add _ hns,
    G.add_neg s hs]

/-- The reflection `φ_s(x) = s - x`. -/
def FinAbGroup.refl {N : Nat} (G : FinAbGroup N) (s x : Nat) : Nat := G.add s (G.neg x)

theorem FinAbGroup.refl_lt {N : Nat} (G : FinAbGroup N)
    (s x : Nat) (hs : s < N) (hx : x < N) : G.refl s x < N :=
  G.add_lt _ _ hs (G.neg_lt x hx)

theorem FinAbGroup.refl_refl {N : Nat} (G : FinAbGroup N)
    (s x : Nat) (hs : s < N) (hx : x < N) : G.refl s (G.refl s x) = x := by
  unfold refl
  have hns := G.neg_lt s hs
  rw [G.neg_sub s x hs hx, ← G.add_assoc s x (G.neg s) hs hx hns, G.add_comm s x hs hx,
    G.add_assoc x s (G.neg s) hx hs hns, G.add_neg s hs, G.add_zero x hx]

/-- `φ_s` reverses differences: `φ b - φ a = a - b`. -/
theorem FinAbGroup.refl_diff {N : Nat} (G : FinAbGroup N)
    (s a b : Nat) (hs : s < N) (ha : a < N) (hb : b < N) :
    G.add (G.refl s b) (G.neg (G.refl s a)) = G.add a (G.neg b) := by
  unfold refl
  have hns := G.neg_lt s hs
  have hnb := G.neg_lt b hb
  rw [G.neg_sub s a hs ha, G.add_assoc s (G.neg b) (G.add a (G.neg s)) hs hnb (G.add_lt _ _ ha hns),
    ← G.add_assoc (G.neg b) a (G.neg s) hnb ha hns,
    G.add_comm (G.add (G.neg b) a) (G.neg s) (G.add_lt _ _ hnb ha) hns,
    ← G.add_assoc s (G.neg s) (G.add (G.neg b) a) hs hns (G.add_lt _ _ hnb ha), G.add_neg s hs,
    G.zero_add _ (G.add_lt _ _ hnb ha), G.add_comm (G.neg b) a hnb ha]

/-- `φ_{u+v}` swaps `u` and `v`. -/
theorem FinAbGroup.refl_swap {N : Nat} (G : FinAbGroup N) (u v : Nat) (hu : u < N) (hv : v < N) :
    G.refl (G.add u v) u = v ∧ G.refl (G.add u v) v = u := by
  unfold refl
  constructor
  · rw [G.add_comm u v hu hv, G.add_assoc v u (G.neg u) hv hu (G.neg_lt u hu), G.add_neg u hu,
      G.add_zero v hv]
  · rw [G.add_assoc u v (G.neg v) hu hv (G.neg_lt v hv), G.add_neg v hv, G.add_zero u hu]


/-- The Cayley digraph `Cay(Γ, S)`: `a → b` iff `b - a ∈ S`. -/
def cayley {N : Nat} (G : FinAbGroup N) (S : Nat → Bool) (a b : Nat) : Bool :=
  S (G.add b (G.neg a))

/-- **THM-4529 Theorem 5.1.** In a Cayley digraph of an abelian group of odd order, every arc
lies on an even number of Hamiltonian paths. -/
theorem cayley_arcCount_even (N : Nat) (hN : N % 2 = 1) (G : FinAbGroup N) (S : Nat → Bool)
    (u v : Nat) (hu : u < N) (hv : v < N) : arcCount N (cayley G S) u v % 2 = 0 := by
  have hs := G.add_lt u v hu hv
  obtain ⟨h1, h2⟩ := G.refl_swap u v hu hv
  refine reflection_arcCount_even N hN (cayley G S) (G.refl (G.add u v))
    (fun x hx => G.refl_lt _ x hs hx) (fun x hx => G.refl_refl _ x hs hx) ?_ u v h1 h2
  intro a b ha hb _
  unfold cayley
  rw [G.refl_diff _ a b hs ha hb]

/-! ## Cyclic groups and direct products -/

/-- Addition in `Z/n`. -/
def cadd (n a b : Nat) : Nat := if a + b < n then a + b else a + b - n

/-- Negation in `Z/n`. -/
def cneg (n a : Nat) : Nat := if a = 0 then 0 else n - a

/-- The cyclic group `Z/n`, `n ≥ 1`. -/
def cyclic (n : Nat) (hn : 0 < n) : FinAbGroup n where
  add := cadd n
  neg := cneg n
  zero := 0
  add_lt := by
    intro a b ha hb
    unfold cadd
    split <;> omega
  neg_lt := by
    intro a ha
    unfold cneg
    split <;> omega
  zero_lt := hn
  add_assoc := by
    intro a b c ha hb hc
    unfold cadd
    (repeat' split) <;> omega
  add_comm := by
    intro a b _ _
    unfold cadd
    (repeat' split) <;> omega
  zero_add := by
    intro a ha
    unfold cadd
    split <;> omega
  neg_add := by
    intro a ha
    unfold cadd cneg
    (repeat' split) <;> omega

theorem divmod_enc (M i j : Nat) (hj : j < M) : (i * M + j) / M = i ∧ (i * M + j) % M = j := by
  have hM : 0 < M := by omega
  constructor
  · rw [Nat.add_comm, Nat.add_mul_div_right j i hM, Nat.div_eq_of_lt hj, Nat.zero_add]
  · rw [Nat.add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt hj]

theorem enc_lt (N M i j : Nat) (hi : i < N) (hj : j < M) : i * M + j < N * M := by
  have : (i + 1) * M ≤ N * M := Nat.mul_le_mul_right M hi
  rw [Nat.add_mul, Nat.one_mul] at this
  omega

theorem dec_enc (M x : Nat) : x / M * M + x % M = x := by
  have := Nat.div_add_mod x M
  rw [Nat.mul_comm] at this
  exact this

theorem div_lt_of_lt_mul' (N M x : Nat) (hM : 0 < M) (hx : x < N * M) : x / M < N :=
  (Nat.div_lt_iff_lt_mul hM).2 hx

/-- The direct product of two finite abelian groups, on `{0, …, N₁N₂ - 1}` via `x = i N₂ + j`. -/
def FinAbGroup.prod {N1 N2 : Nat} (G1 : FinAbGroup N1) (G2 : FinAbGroup N2) :
    FinAbGroup (N1 * N2) where
  add x y := G1.add (x / N2) (y / N2) * N2 + G2.add (x % N2) (y % N2)
  neg x := G1.neg (x / N2) * N2 + G2.neg (x % N2)
  zero := G1.zero * N2 + G2.zero
  add_lt := by
    intro x y hx hy
    have hM := G2.zero_lt
    exact enc_lt N1 N2 _ _ (G1.add_lt _ _ (div_lt_of_lt_mul' N1 N2 x (by omega) hx)
      (div_lt_of_lt_mul' N1 N2 y (by omega) hy))
      (G2.add_lt _ _ (Nat.mod_lt x (by omega)) (Nat.mod_lt y (by omega)))
  neg_lt := by
    intro x hx
    have hM := G2.zero_lt
    exact enc_lt N1 N2 _ _ (G1.neg_lt _ (div_lt_of_lt_mul' N1 N2 x (by omega) hx))
      (G2.neg_lt _ (Nat.mod_lt x (by omega)))
  zero_lt := enc_lt N1 N2 _ _ G1.zero_lt G2.zero_lt
  add_assoc := by
    intro x y z hx hy hz
    have hM : 0 < N2 := by have := G2.zero_lt; omega
    have hx1 := div_lt_of_lt_mul' N1 N2 x hM hx
    have hy1 := div_lt_of_lt_mul' N1 N2 y hM hy
    have hz1 := div_lt_of_lt_mul' N1 N2 z hM hz
    have hx2 := Nat.mod_lt x hM
    have hy2 := Nat.mod_lt y hM
    have hz2 := Nat.mod_lt z hM
    obtain ⟨e1, e2⟩ := divmod_enc N2 (G1.add (x / N2) (y / N2)) _ (G2.add_lt _ _ hx2 hy2)
    obtain ⟨f1, f2⟩ := divmod_enc N2 (G1.add (y / N2) (z / N2)) _ (G2.add_lt _ _ hy2 hz2)
    rw [e1, e2, f1, f2, G1.add_assoc _ _ _ hx1 hy1 hz1, G2.add_assoc _ _ _ hx2 hy2 hz2]
  add_comm := by
    intro x y hx hy
    have hM : 0 < N2 := by have := G2.zero_lt; omega
    rw [G1.add_comm _ _ (div_lt_of_lt_mul' N1 N2 x hM hx) (div_lt_of_lt_mul' N1 N2 y hM hy),
      G2.add_comm _ _ (Nat.mod_lt x hM) (Nat.mod_lt y hM)]
  zero_add := by
    intro x hx
    have hM : 0 < N2 := by have := G2.zero_lt; omega
    obtain ⟨e1, e2⟩ := divmod_enc N2 G1.zero G2.zero G2.zero_lt
    show G1.add ((G1.zero * N2 + G2.zero) / N2) (x / N2) * N2 +
      G2.add ((G1.zero * N2 + G2.zero) % N2) (x % N2) = x
    rw [e1, e2, G1.zero_add _ (div_lt_of_lt_mul' N1 N2 x hM hx), G2.zero_add _ (Nat.mod_lt x hM),
      dec_enc]
  neg_add := by
    intro x hx
    have hM : 0 < N2 := by have := G2.zero_lt; omega
    have hx1 := div_lt_of_lt_mul' N1 N2 x hM hx
    have hx2 := Nat.mod_lt x hM
    obtain ⟨e1, e2⟩ := divmod_enc N2 (G1.neg (x / N2)) (G2.neg (x % N2)) (G2.neg_lt _ hx2)
    show G1.add ((G1.neg (x / N2) * N2 + G2.neg (x % N2)) / N2) (x / N2) * N2 +
      G2.add ((G1.neg (x / N2) * N2 + G2.neg (x % N2)) % N2) (x % N2) = G1.zero * N2 + G2.zero
    rw [e1, e2, G1.neg_add _ hx1, G2.neg_add _ hx2]

/-- **Theorem 5.1 for `Z/m × Z/n`** (`m`, `n` odd): every Cayley digraph of
`Z/m × Z/n` is arc-even; e.g. `Z/3 × Z/3`, which is not cyclic. -/
theorem cayley_cyclic_prod_arcCount_even (m n : Nat) (hm : m % 2 = 1) (hn : n % 2 = 1)
    (S : Nat → Bool) (u v : Nat) (hu : u < m * n) (hv : v < m * n) :
    arcCount (m * n) (cayley ((cyclic m (by omega)).prod (cyclic n (by omega))) S) u v % 2 = 0 :=
  cayley_arcCount_even (m * n) (by rw [Nat.mul_mod, hm, hn]) _ S u v hu hv

end ProcgenSelfieEdim
