import ProcgenSelfieEdim.AntiAut
import ProcgenSelfieEdim.Hypercube

set_option autoImplicit false

/-!
# THM-4529 Theorem 6.1: antipodal arcs of anti-circulant tournaments are odd

Let `m` be odd and `N = 2m`. An *anti-circulant* tournament on `Z/N` is one for which
`τ(x) = x + 1` is an anti-automorphism (THM-4529 §6.3 shows these are exactly the `T_s`:
`x → y` iff `s(y - x) = (-1)^x`, with `s(-d) = -(-1)^d s(d)`).

* `cyclic_anti_antipodal_odd` (**Theorem 6.1**, general form): if `T` is a tournament on
  `Z/2m`, `m` odd, with the anti-automorphism `x ↦ x + 1`, then for every `x` the antipodal
  pair `{x, x + m}` carries an arc lying on an odd number of Hamiltonian paths
  (`c(x → x+m) + c(x+m → x)` is odd; one of the two terms is `0`).
* `antiCirc m s` is `T_s`; `antiCirc_isTournament` and `antiCirc_anti` check the definition;
  `antiCirc_antipodal_odd` is Theorem 6.1 as stated in the note.

Proof (the note's): `σ: P ↦ reverse(τ P)` maps HPs to HPs and the arc `a → b` to
`τ b → τ a`, so `c(a → b) = c(τ b → τ a)` (`arcCount_anti`). Group the double sum
`Σ c(x → y) = (N - 1) H` by differences `d = y - x`: the class sums satisfy
`C(d) = C(N - d)` (`σ` maps difference `d` to `N - d`), so the sum is `≡ C(0) + C(m)
(mod 2)`, with `C(0) = 0`. The antipodal class is `C(m) = m · A`, where
`A(x) = c(x → x+m) + c(x+m → x)` is `σ`-invariant, hence constant. Rédei (`H` odd) and
`N - 1`, `m` odd give `A` odd.
-/

namespace ProcgenSelfieEdim

/-! ## Arithmetic on `{0, …, N-1}` -/

/-- `x + 1 (mod N)`. -/
def succMod (N x : Nat) : Nat := if x + 1 < N then x + 1 else 0

/-- `x - 1 (mod N)`. -/
def predMod (N x : Nat) : Nat := if x = 0 then N - 1 else x - 1

/-- `x + d (mod N)` for `x, d < N`. -/
def addMod (N x d : Nat) : Nat := if x + d < N then x + d else x + d - N

theorem succMod_lt (N x : Nat) (hx : x < N) : succMod N x < N := by
  unfold succMod
  split <;> omega

theorem predMod_lt (N x : Nat) (hx : x < N) : predMod N x < N := by
  unfold predMod
  split <;> omega

theorem pred_succ (N x : Nat) (hx : x < N) : predMod N (succMod N x) = x := by
  unfold succMod predMod
  (repeat' split) <;> omega

theorem succ_pred (N x : Nat) (hx : x < N) : succMod N (predMod N x) = x := by
  unfold succMod predMod
  (repeat' split) <;> omega

theorem addMod_lt (N x d : Nat) (hx : x < N) (hd : d < N) : addMod N x d < N := by
  unfold addMod
  split <;> omega

theorem addMod_zero (N x : Nat) (hx : x < N) : addMod N x 0 = x := by
  unfold addMod
  split <;> omega

theorem zero_addMod (N d : Nat) (hd : d < N) : addMod N 0 d = d := by
  unfold addMod
  split <;> omega

theorem addMod_comm (N x d : Nat) : addMod N x d = addMod N d x := by
  unfold addMod
  (repeat' split) <;> omega

theorem addMod_succ (N x d : Nat) (hx : x + 1 < N) (hd : d < N) :
    addMod N (x + 1) d = addMod N x (succMod N d) := by
  unfold addMod succMod
  (repeat' split) <;> omega

theorem succ_addMod (N x d : Nat) (hx : x < N) (hd : d < N) :
    succMod N (addMod N x d) = addMod N x (succMod N d) := by
  unfold addMod succMod
  (repeat' split) <;> omega

theorem succ_addMod' (N x d : Nat) (hx : x < N) (hd : d < N) :
    succMod N (addMod N x d) = addMod N (succMod N x) d := by
  unfold addMod succMod
  (repeat' split) <;> omega

theorem addMod_back (N x d : Nat) (hx : x < N) (hd0 : 0 < d) (hd : d < N) :
    addMod N (addMod N x (succMod N d)) (N - d) = succMod N x := by
  unfold addMod succMod
  (repeat' split) <;> omega

theorem addMod_antipode (m x : Nat) (hx : x < m) :
    addMod (2 * m) x m = x + m ∧ addMod (2 * m) (m + x) m = x := by
  unfold addMod
  constructor <;> split <;> omega

theorem cdist_succMod (N a b : Nat) (ha : a < N) (hb : b < N) :
    cdist N (succMod N a) (succMod N b) = cdist N a b := by
  unfold cdist succMod
  (repeat' split) <;> omega

theorem cdist_addMod (N x d : Nat) (hx : x < N) (hd : d < N) : cdist N x (addMod N x d) = d := by
  unfold cdist addMod
  (repeat' split) <;> omega

/-! ## Sums -/

theorem rsum_rot1 (M : Nat) (h : Nat → Nat) :
    rsum (M + 1) (fun d => h (succMod (M + 1) d)) = rsum (M + 1) h := by
  rw [rsum_shift M h]
  show rsum M (fun d => h (succMod (M + 1) d)) + h (succMod (M + 1) M) = _
  have e1 : succMod (M + 1) M = 0 := by
    unfold succMod
    rw [if_neg (by omega)]
  have e2 : rsum M (fun d => h (succMod (M + 1) d)) = rsum M (fun r => h (r + 1)) :=
    rsum_congr M _ _ (fun d hd => by
      show h (succMod (M + 1) d) = h (d + 1)
      unfold succMod
      rw [if_pos (by omega)])
  rw [e1, e2]
  omega

/-- Rotation invariance: `Σ_d g(x + d) = Σ_d g(d)` over `Z/N`. -/
theorem rsum_rot (N : Nat) (g : Nat → Nat) : ∀ x, x < N → rsum N (fun d => g (addMod N x d)) = rsum N g
  | 0, _ => rsum_congr N _ _ (fun d hd => by
      show g (addMod N 0 d) = g d
      rw [zero_addMod N d hd])
  | x + 1, hx => by
    have e : rsum N (fun d => g (addMod N (x + 1) d)) =
        rsum N (fun d => (fun e => g (addMod N x e)) (succMod N d)) :=
      rsum_congr N _ _ (fun d hd => by
        show g (addMod N (x + 1) d) = g (addMod N x (succMod N d))
        rw [addMod_succ N x d hx hd])
    rw [e]
    obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
    rw [rsum_rot1 M (fun e => g (addMod (M + 1) x e))]
    exact rsum_rot (M + 1) g x (by omega)

/-- Grouping a double sum over `Z/N × Z/N` by differences. -/
theorem total_by_diff (N : Nat) (c : Nat → Nat → Nat) :
    rsum N (fun x => rsum N (fun y => c x y)) =
      rsum N (fun d => rsum N (fun x => c x (addMod N x d))) := by
  rw [← rsum_comm N (fun x d => c x (addMod N x d)) N]
  exact rsum_congr N _ _ (fun x hx => (rsum_rot N (fun y => c x y) x hx).symm)

theorem rsum_split (f : Nat → Nat) (a : Nat) :
    ∀ b, rsum (a + b) f = rsum a f + rsum b (fun i => f (a + i))
  | 0 => by
    show rsum a f = rsum a f + 0
    omega
  | b + 1 => by
    show rsum (a + b) f + f (a + b) = rsum a f + (rsum b (fun i => f (a + i)) + f (a + b))
    rw [rsum_split f a b]
    omega

theorem rsum_reverse (g : Nat → Nat) : ∀ k, rsum k (fun i => g (k - 1 - i)) = rsum k g
  | 0 => rfl
  | k + 1 => by
    rw [rsum_shift k (fun i => g (k + 1 - 1 - i))]
    have h1 : rsum k (fun r => g (k + 1 - 1 - (r + 1))) = rsum k (fun r => g (k - 1 - r)) :=
      rsum_congr k _ _ (fun r _ => by
        show g (k + 1 - 1 - (r + 1)) = g (k - 1 - r)
        rw [show k + 1 - 1 - (r + 1) = k - 1 - r by omega])
    show g (k + 1 - 1 - 0) + rsum k (fun r => g (k + 1 - 1 - (r + 1))) = rsum k g + g k
    rw [h1, rsum_reverse g k, show k + 1 - 1 - 0 = k by omega]
    omega

/-- A sum over `Z/2m` of a function with `f(d) = f(2m - d)` is `≡ f(0) + f(m) (mod 2)`. -/
theorem rsum_sym_mod_two (m : Nat) (hm : 1 ≤ m) (f : Nat → Nat)
    (hf : ∀ d, 0 < d → d < 2 * m → f d = f (2 * m - d)) :
    rsum (2 * m) f % 2 = (f 0 + f m) % 2 := by
  obtain ⟨k, rfl⟩ : ∃ k, m = k + 1 := ⟨m - 1, by omega⟩
  have e0 : 2 * (k + 1) = (k + 2) + k := by omega
  rw [e0, rsum_split f (k + 2) k]
  -- the upper half, reflected
  have e1 : rsum k (fun i => f (k + 2 + i)) = rsum k (fun i => (fun j => f (j + 1)) (k - 1 - i)) :=
    rsum_congr k _ _ (fun i hi => by
      show f (k + 2 + i) = f (k - 1 - i + 1)
      rw [hf (k + 2 + i) (by omega) (by omega)]
      congr 1
      omega)
  rw [e1, rsum_reverse (fun j => f (j + 1)) k]
  -- the lower half
  have e2 : rsum (k + 2) f = f 0 + rsum k (fun r => f (r + 1)) + f (k + 1) := by
    show rsum (k + 1) f + f (k + 1) = _
    rw [rsum_shift k f]
  rw [e2]
  omega

theorem rsum_const (k : Nat) : ∀ n, rsum n (fun _ => k) = n * k
  | 0 => by rw [Nat.zero_mul]; rfl
  | n + 1 => by
    show rsum n (fun _ => k) + k = (n + 1) * k
    rw [rsum_const k n, Nat.add_mul, Nat.one_mul]

/-! ## Theorem 6.1, general form -/

theorem arcCount_zero_of_not_arc (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) (h : T u v = false) :
    arcCount n T u v = 0 := by
  apply arcCount_eq_zero_of_not
  intro p hp
  cases hc : usesArc u v p
  · rfl
  · rw [usesArc_arc T u v p hp.2.2.2 hc] at h
    exact Bool.noConfusion h

/-- **THM-4529 Theorem 6.1.** Let `T` be a tournament on `Z/2m` (`m` odd) for which
`x ↦ x + 1` is an anti-automorphism. Then for every `x`, the arc between `x` and `x + m`
lies on an odd number of Hamiltonian paths. -/
theorem cyclic_anti_antipodal_odd (m : Nat) (hm : m % 2 = 1) (T : Nat → Nat → Bool)
    (hT : IsTournament (2 * m) T)
    (hanti : ∀ a b, a < 2 * m → b < 2 * m → a ≠ b →
      T (succMod (2 * m) a) (succMod (2 * m) b) = T b a)
    (x : Nat) (hx : x < 2 * m) :
    (arcCount (2 * m) T x (addMod (2 * m) x m) + arcCount (2 * m) T (addMod (2 * m) x m) x) % 2 = 1 := by
  have hN : 0 < 2 * m := by omega
  -- `σ`-invariance: `c(a → b) = c(τ b → τ a)`
  have hσ : ∀ a b, a < 2 * m → b < 2 * m →
      arcCount (2 * m) T a b = arcCount (2 * m) T (succMod (2 * m) b) (succMod (2 * m) a) :=
    fun a b ha hb => arcCount_anti (2 * m) T (succMod (2 * m)) (predMod (2 * m))
      (succMod_lt (2 * m)) (predMod_lt (2 * m)) (pred_succ (2 * m)) (succ_pred (2 * m)) hanti
      a b ha hb
  -- the class sums `C(d)` and their symmetry
  let C : Nat → Nat := fun d => rsum (2 * m) (fun y => arcCount (2 * m) T y (addMod (2 * m) y d))
  have hsym : ∀ d, 0 < d → d < 2 * m → C d = C (2 * m - d) := by
    intro d hd0 hd
    have e1 : ∀ y, y < 2 * m → arcCount (2 * m) T y (addMod (2 * m) y d) =
        (fun z => arcCount (2 * m) T z (addMod (2 * m) z (2 * m - d)))
          (addMod (2 * m) (succMod (2 * m) d) y) := by
      intro y hy
      show _ = arcCount (2 * m) T (addMod (2 * m) (succMod (2 * m) d) y)
        (addMod (2 * m) (addMod (2 * m) (succMod (2 * m) d) y) (2 * m - d))
      rw [hσ y _ hy (addMod_lt _ _ _ hy hd), succ_addMod _ _ _ hy hd,
        addMod_comm (2 * m) (succMod (2 * m) d) y, addMod_back _ _ _ hy hd0 hd]
    show rsum (2 * m) _ = rsum (2 * m) _
    rw [rsum_congr _ _ _ e1]
    exact rsum_rot (2 * m) (fun z => arcCount (2 * m) T z (addMod (2 * m) z (2 * m - d)))
      (succMod (2 * m) d) (succMod_lt _ _ hd)
  -- `C(0) = 0`
  have hC0 : C 0 = 0 := by
    show rsum (2 * m) _ = 0
    rw [rsum_congr (2 * m) _ (fun _ => 0) (fun y hy => by
      show arcCount (2 * m) T y (addMod (2 * m) y 0) = 0
      rw [addMod_zero _ _ hy]
      exact arcCount_eq_zero_of_not (2 * m) T y y (fun p hp => usesArc_self y p hp.2.1))]
    exact rsum_zero _
  -- `A(y) = c(y → y+m) + c(y+m → y)` is constant
  let A : Nat → Nat := fun y =>
    arcCount (2 * m) T y (addMod (2 * m) y m) + arcCount (2 * m) T (addMod (2 * m) y m) y
  have hAτ : ∀ y, y < 2 * m → A (succMod (2 * m) y) = A y := by
    intro y hy
    have hym := addMod_lt (2 * m) y m hy (by omega)
    show arcCount (2 * m) T (succMod (2 * m) y) (addMod (2 * m) (succMod (2 * m) y) m) +
        arcCount (2 * m) T (addMod (2 * m) (succMod (2 * m) y) m) (succMod (2 * m) y) =
      arcCount (2 * m) T y (addMod (2 * m) y m) + arcCount (2 * m) T (addMod (2 * m) y m) y
    rw [hσ y _ hy hym, hσ _ y hym hy, ← succ_addMod' _ _ _ hy (by omega)]
    omega
  have hAconst : ∀ y, y < 2 * m → A y = A 0 := by
    intro y
    induction y with
    | zero => intro _; rfl
    | succ k ih =>
      intro hk
      have : succMod (2 * m) k = k + 1 := by
        unfold succMod
        rw [if_pos hk]
      rw [← this, hAτ k (by omega), ih (by omega)]
  -- `C(m) = m · A(0)`
  have hCm : C m = m * A 0 := by
    show rsum (2 * m) _ = _
    rw [show 2 * m = m + m by omega, rsum_split _ m m, ← rsum_add]
    rw [rsum_congr m _ (fun _ => A 0) (fun y hy => by
      obtain ⟨h1, h2⟩ := addMod_antipode m y hy
      show arcCount (m + m) T y (addMod (m + m) y m) +
          arcCount (m + m) T (m + y) (addMod (m + m) (m + y) m) = A 0
      rw [← hAconst y (by omega)]
      rw [show m + m = 2 * m by omega, h1, h2]
      show _ = arcCount (2 * m) T y (addMod (2 * m) y m) + arcCount (2 * m) T (addMod (2 * m) y m) y
      rw [h1, Nat.add_comm m y])]
    exact rsum_const (A 0) m
  -- the total is `(2m - 1) H`, odd by Rédei
  have htot := sum_arcCount (2 * m) T
  have hH := redei (2 * m) T hT
  rw [total_by_diff] at htot
  have hpar := rsum_sym_mod_two m (by omega) C hsym
  have hodd : ((2 * m - 1) * hpCount (2 * m) T) % 2 = 1 := by
    rw [Nat.mul_mod, hH, show (2 * m - 1) % 2 = 1 by omega]
  have hsum : rsum (2 * m) C % 2 = 1 := by
    show rsum (2 * m) (fun d => rsum (2 * m) (fun y => arcCount (2 * m) T y (addMod (2 * m) y d))) % 2 = 1
    rw [htot]
    exact hodd
  rw [hpar, hC0, hCm] at hsum
  have hA0 : A 0 % 2 = 1 := by
    rw [Nat.zero_add, Nat.mul_mod, hm, Nat.one_mul, Nat.mod_mod] at hsum
    exact hsum
  have := hAconst x hx
  show A x % 2 = 1
  rw [this]
  exact hA0

/-! ## The tournaments `T_s` -/

/-- The anti-circulant tournament `T_s` on `Z/2m`: `x → y` iff `s(y - x) = (-1)^x`
(`true` stands for `+1`). -/
def antiCirc (m : Nat) (s : Nat → Bool) (x y : Nat) : Bool :=
  if x % 2 = 0 then s (cdist (2 * m) x y) else !s (cdist (2 * m) x y)

/-- The sign condition `s(-d) = -(-1)^d s(d)`. -/
def AntiSign (m : Nat) (s : Nat → Bool) : Prop :=
  ∀ d, 0 < d → d < 2 * m → s (2 * m - d) = (if d % 2 = 0 then !s d else s d)

theorem antiCirc_isTournament (m : Nat) (s : Nat → Bool) (hs : AntiSign m s) :
    IsTournament (2 * m) (antiCirc m s) := by
  intro a b ha hb hab
  have h1 := cdist_spec (2 * m) a b
  have h2 := cdist_spec (2 * m) b a
  have e : cdist (2 * m) b a = 2 * m - cdist (2 * m) a b := by omega
  have hpar : b % 2 = (a + cdist (2 * m) a b) % 2 := by omega
  have hd0 : 0 < cdist (2 * m) a b := by omega
  have hd1 : cdist (2 * m) a b < 2 * m := by omega
  unfold antiCirc
  rw [e, hs _ hd0 hd1]
  generalize cdist (2 * m) a b = d at hpar hd0 hd1
  cases hsd : s d <;> (repeat' split) <;> first | rfl | (exfalso; omega)

/-- `x ↦ x + 1` is an anti-automorphism of `T_s`. -/
theorem antiCirc_anti (m : Nat) (s : Nat → Bool) (hs : AntiSign m s) :
    ∀ a b, a < 2 * m → b < 2 * m → a ≠ b →
      antiCirc m s (succMod (2 * m) a) (succMod (2 * m) b) = antiCirc m s b a := by
  intro a b ha hb hab
  rw [antiCirc_isTournament m s hs a b ha hb hab]
  unfold antiCirc
  rw [cdist_succMod _ _ _ ha hb]
  have hpar : succMod (2 * m) a % 2 = (a + 1) % 2 := by
    unfold succMod
    split <;> omega
  generalize cdist (2 * m) a b = d
  cases hsd : s d <;> (repeat' split) <;> first | rfl | (exfalso; omega)

/-- **THM-4529 Theorem 6.1, as stated.** In the anti-circulant tournament `T_s` on `2m`
vertices (`m` odd), every antipodal arc `x → x + m` lies on an odd number of Hamiltonian
paths. -/
theorem antiCirc_antipodal_odd (m : Nat) (hm : m % 2 = 1) (s : Nat → Bool) (hs : AntiSign m s)
    (x : Nat) (hx : x < 2 * m) (harc : antiCirc m s x (addMod (2 * m) x m) = true) :
    arcCount (2 * m) (antiCirc m s) x (addMod (2 * m) x m) % 2 = 1 := by
  have h := cyclic_anti_antipodal_odd m hm (antiCirc m s) (antiCirc_isTournament m s hs)
    (antiCirc_anti m s hs) x hx
  have hxm := addMod_lt (2 * m) x m hx (by omega)
  have hne : x ≠ addMod (2 * m) x m := by
    unfold addMod
    split <;> omega
  have hback : antiCirc m s (addMod (2 * m) x m) x = false := by
    rw [antiCirc_isTournament m s hs x _ hx hxm hne, harc]
    rfl
  rw [arcCount_zero_of_not_arc _ _ _ _ hback, Nat.add_zero] at h
  exact h

/-- The sign condition forces `m` odd (for `m ≥ 1`), matching `N = 2m ≡ 2 (mod 4)`
(`anti_cycle_mod_four`). -/
theorem antiSign_odd (m : Nat) (hm : 1 ≤ m) (s : Nat → Bool) (hs : AntiSign m s) : m % 2 = 1 := by
  have := hs m (by omega) (by omega)
  rw [show 2 * m - m = m by omega] at this
  by_cases h : m % 2 = 0
  · rw [if_pos h] at this
    cases hsm : s m <;> rw [hsm] at this <;> exact Bool.noConfusion this
  · omega

end ProcgenSelfieEdim
