set_option autoImplicit false

/-!
# The Syracuse drop and the owner's label map (opus S15, THM-4527 Theorems 1-2)

For odd `A > 0`: `S(A) = oddpart(3A + 1) = (3A + 1) / 2^v` with `v = v_2(3A + 1)`,
and the drop is `K(A) = (A - S(A)) / 2 ∈ Z`.

* `drop_identity` (Thm 2(a)): `6 K(A) + 1 = (2^v - 3) S(A)`.
* `drop_forward`, `drop_backward`, `drop_injective` (Thm 2(b)): for every `m ∈ Z`,
  `A ↦ v_2(3A + 1)` is a bijection from the odd `A > 0` with `K(A) = m` onto the
  `v ≥ 1` with `(2^v - 3) ∣ 6m + 1` and `(6m + 1) / (2^v - 3) > 0`.
* `admissible_of_neg`, `admissible_of_nonneg`: for `m < 0` only `v = 1` is admissible;
  for `m ≥ 0` exactly the `v ≥ 2` with `(2^v - 3) ∣ 6m + 1`. Descents
  (`A ≡ 1 mod 4`) are exactly the `A` with `v ≥ 2` (`mod_four_iff`).
* `labelF_even`, `labelF_four_j_one`, `labelF_microcosm` (Thm 1(b)): the label map
  `F(M) = (S(2M - 1) + 1) / 2` satisfies `F(2N) = 3N`, `F(4j + 1) = 3j + 1`,
  `F(4n - 1) = F(n)`; and `F(M) = M - K(2M - 1)` (`labelF_eq_sub_drop`).
-/

namespace ProcgenSelfieEdim

/-! ## 2-adic valuation and odd part -/

/-- Strip up to `k` factors of two, counting them. -/
def v2Aux : Nat → Nat → Nat
  | 0, _ => 0
  | k + 1, n => if n % 2 = 0 ∧ 0 < n then v2Aux k (n / 2) + 1 else 0

/-- Strip up to `k` factors of two. -/
def oddAux : Nat → Nat → Nat
  | 0, n => n
  | k + 1, n => if n % 2 = 0 ∧ 0 < n then oddAux k (n / 2) else n

/-- The 2-adic valuation `v_2(n)` of `n > 0`. -/
def v2 (n : Nat) : Nat := v2Aux n n

/-- The odd part of `n > 0`. -/
def oddPart (n : Nat) : Nat := oddAux n n

theorem aux_spec : ∀ (v k s : Nat), s % 2 = 1 → v < k →
    v2Aux k (2 ^ v * s) = v ∧ oddAux k (2 ^ v * s) = s
  | 0, 0, _, _, hk => absurd hk (Nat.lt_irrefl 0)
  | 0, k + 1, s, hs, _ => by
    rw [Nat.pow_zero, Nat.one_mul]
    have hc : ¬ (s % 2 = 0 ∧ 0 < s) := fun h => by omega
    exact ⟨by rw [v2Aux, if_neg hc], by rw [oddAux, if_neg hc]⟩
  | _ + 1, 0, _, _, hk => absurd hk (Nat.not_lt_zero _)
  | v + 1, k + 1, s, hs, hk => by
    have hpos : 0 < 2 ^ v * s := Nat.mul_pos (Nat.two_pow_pos v) (by omega)
    have heq : 2 ^ (v + 1) * s = 2 * (2 ^ v * s) := by
      rw [Nat.pow_succ, Nat.mul_comm (2 ^ v) 2, Nat.mul_assoc]
    have hc : 2 * (2 ^ v * s) % 2 = 0 ∧ 0 < 2 * (2 ^ v * s) := ⟨by omega, by omega⟩
    have hdiv : 2 * (2 ^ v * s) / 2 = 2 ^ v * s := by omega
    have ih := aux_spec v k s hs (by omega)
    rw [heq]
    exact ⟨by rw [v2Aux, if_pos hc, hdiv, ih.1], by rw [oddAux, if_pos hc, hdiv, ih.2]⟩

/-- **Uniqueness of the decomposition** `n = 2^v s` (`s` odd). -/
theorem v2_oddPart_of_eq (n v s : Nat) (h : n = 2 ^ v * s) (hs : s % 2 = 1) :
    v2 n = v ∧ oddPart n = s := by
  have hv : v < n := by
    have h1 : v < 2 ^ v := Nat.lt_two_pow_self
    have h2 : 2 ^ v ≤ 2 ^ v * s := Nat.le_mul_of_pos_right _ (by omega)
    omega
  subst h
  exact aux_spec v _ s hs hv

/-- **Existence of the decomposition.** -/
theorem exists_decomp : ∀ (k n : Nat), 0 < n → n < k → ∃ v s, n = 2 ^ v * s ∧ s % 2 = 1
  | 0, _, _, hk => absurd hk (Nat.not_lt_zero _)
  | k + 1, n, hn, hk => by
    by_cases he : n % 2 = 0
    · obtain ⟨v, s, h1, h2⟩ := exists_decomp k (n / 2) (by omega) (by omega)
      refine ⟨v + 1, s, ?_, h2⟩
      rw [Nat.pow_succ, Nat.mul_comm (2 ^ v) 2, Nat.mul_assoc, ← h1]
      omega
    · exact ⟨0, n, by rw [Nat.pow_zero, Nat.one_mul], by omega⟩

theorem decomp (n : Nat) (hn : 0 < n) : n = 2 ^ v2 n * oddPart n ∧ oddPart n % 2 = 1 := by
  obtain ⟨v, s, h1, h2⟩ := exists_decomp (n + 1) n hn (Nat.lt_succ_self n)
  obtain ⟨hv, hs⟩ := v2_oddPart_of_eq n v s h1 h2
  rw [hv, hs]
  exact ⟨h1, h2⟩

theorem oddPart_two_pow_mul (a x : Nat) (hx : 0 < x) : oddPart (2 ^ a * x) = oddPart x := by
  obtain ⟨hx1, hx2⟩ := decomp x hx
  have : 2 ^ a * x = 2 ^ (a + v2 x) * oddPart x := by
    rw [Nat.pow_add, Nat.mul_assoc, ← hx1]
  exact (v2_oddPart_of_eq _ _ _ this hx2).2

/-! ## The Syracuse map and the drop -/

/-- The Syracuse map `S(A) = oddpart(3A + 1)`. -/
def syr (A : Nat) : Nat := oddPart (3 * A + 1)

/-- The drop `K(A) = (A - S(A)) / 2` (an integer; negative on ascents). -/
def drop (A : Nat) : Int := ((A : Int) - (syr A : Int)) / 2

theorem syr_decomp (A : Nat) :
    3 * A + 1 = 2 ^ v2 (3 * A + 1) * syr A ∧ syr A % 2 = 1 :=
  decomp (3 * A + 1) (by omega)

/-- **Theorem 2(a).** `6 K(A) + 1 = (2^v - 3) S(A)` with `v = v_2(3A + 1)`, for odd `A`. -/
theorem drop_identity (A : Nat) (hA : A % 2 = 1) :
    6 * drop A + 1 = ((2 : Int) ^ v2 (3 * A + 1) - 3) * (syr A : Int) := by
  obtain ⟨h1, h2⟩ := syr_decomp A
  have hc : ((3 * A + 1 : Nat) : Int) = (2 : Int) ^ v2 (3 * A + 1) * (syr A : Int) := by
    have e := congrArg (fun x : Nat => (x : Int)) h1
    simp only [Int.natCast_mul, Int.natCast_pow] at e
    exact e
  unfold drop
  rw [Int.sub_mul]
  omega

theorem two_pow_sub_three_ne_zero (v : Nat) : (2 : Int) ^ v - 3 ≠ 0 := by
  intro h
  have h3 : (2 : Int) ^ v = 3 := by omega
  have h4 : ((2 ^ v : Nat) : Int) = ((3 : Nat) : Int) := by rw [Int.natCast_pow]; exact h3
  have h5 : 2 ^ v = 3 := Int.ofNat.inj h4
  cases v with
  | zero => rw [Nat.pow_zero] at h5; omega
  | succ w => rw [Nat.pow_succ] at h5; omega

/-- **Theorem 2(b), forward.** For odd `A`, `v = v_2(3A + 1) ≥ 1` is admissible for
`m = K(A)`: `(2^v - 3) ∣ 6m + 1` with quotient `S(A) > 0`. -/
theorem drop_forward (A : Nat) (hA : A % 2 = 1) :
    1 ≤ v2 (3 * A + 1) ∧ ((2 : Int) ^ v2 (3 * A + 1) - 3) ∣ (6 * drop A + 1) ∧
      (6 * drop A + 1) / ((2 : Int) ^ v2 (3 * A + 1) - 3) = (syr A : Int) ∧ 0 < syr A := by
  obtain ⟨h1, h2⟩ := syr_decomp A
  have hid := drop_identity A hA
  refine ⟨?_, ⟨syr A, hid⟩, ?_, by omega⟩
  · cases hv : v2 (3 * A + 1) with
    | zero =>
      rw [hv, Nat.pow_zero, Nat.one_mul] at h1
      omega
    | succ w => omega
  · rw [hid]
    exact Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero _)

/-- **Theorem 2(b), injectivity.** An odd `A` is determined by `K(A)` and `v_2(3A + 1)`. -/
theorem drop_injective (A A' : Nat) (hA : A % 2 = 1) (hA' : A' % 2 = 1) (hd : drop A = drop A')
    (hv : v2 (3 * A + 1) = v2 (3 * A' + 1)) : A = A' := by
  have h1 := drop_identity A hA
  have h2 := drop_identity A' hA'
  rw [← hd, ← hv] at h2
  have hs : (syr A : Int) = (syr A' : Int) := by
    have := h1.symm.trans h2
    have hne := two_pow_sub_three_ne_zero (v2 (3 * A + 1))
    have e1 := Int.mul_ediv_cancel_left (syr A : Int) hne
    have e2 := Int.mul_ediv_cancel_left (syr A' : Int) hne
    rw [this] at e1
    exact e1.symm.trans e2
  have hs' : syr A = syr A' := Int.ofNat.inj hs
  obtain ⟨d1, _⟩ := syr_decomp A
  obtain ⟨d2, _⟩ := syr_decomp A'
  rw [hv, hs'] at d1
  omega

/-- **Theorem 2(b), backward.** Every admissible `v ≥ 1` for `m` is `v_2(3A + 1)` for an
odd `A > 0` with `K(A) = m` (namely `A = 2m + s`, `s = (6m + 1)/(2^v - 3)`). -/
theorem drop_backward (m : Int) (v : Nat) (hv : 1 ≤ v) (hdvd : ((2 : Int) ^ v - 3) ∣ (6 * m + 1))
    (hpos : 0 < (6 * m + 1) / ((2 : Int) ^ v - 3)) :
    ∃ A : Nat, A % 2 = 1 ∧ 0 < A ∧ v2 (3 * A + 1) = v ∧ drop A = m := by
  -- the quotient `s` and the relation `6m + 1 = (2^v - 3) s`
  have hq := Int.mul_ediv_cancel' hdvd
  generalize hs : (6 * m + 1) / ((2 : Int) ^ v - 3) = s at hq hpos
  -- `2^v` is even
  obtain ⟨w, rfl⟩ : ∃ w, v = w + 1 := ⟨v - 1, by omega⟩
  have hpow : (2 : Int) ^ (w + 1) = 2 * (2 : Int) ^ w := by rw [Int.pow_succ, Int.mul_comm]
  rw [hpow] at hq
  -- `s` is odd
  have hsodd : s % 2 = 1 := by
    have : (2 * (2 : Int) ^ w - 3) * s = 2 * ((2 : Int) ^ w * s) - 3 * s := by
      rw [Int.sub_mul, Int.mul_assoc]
    rw [this] at hq
    omega
  -- the natural number `A = 2m + s`
  have hA0 : 0 < 2 * m + s := by
    have hP : 0 < (2 : Int) ^ w * s := Int.mul_pos (Int.pow_pos (by decide)) hpos
    have : (2 * (2 : Int) ^ w - 3) * s = 2 * ((2 : Int) ^ w * s) - 3 * s := by
      rw [Int.sub_mul, Int.mul_assoc]
    rw [this] at hq
    omega
  refine ⟨(2 * m + s).toNat, ?_, ?_, ?_⟩
  · have := Int.toNat_of_nonneg (Int.le_of_lt hA0)
    omega
  · have := Int.toNat_of_nonneg (Int.le_of_lt hA0)
    omega
  -- `3A + 1 = 2^(w+1) s` in `Nat`
  have hsN := Int.toNat_of_nonneg (Int.le_of_lt hpos)
  have hAN := Int.toNat_of_nonneg (Int.le_of_lt hA0)
  have key : 3 * (2 * m + s).toNat + 1 = 2 ^ (w + 1) * s.toNat := by
    apply Int.natCast_inj.mp
    have e1 : ((3 * (2 * m + s).toNat + 1 : Nat) : Int) = 3 * (2 * m + s) + 1 := by
      rw [Int.natCast_add, Int.natCast_mul, hAN]
      rfl
    have e2 : ((2 ^ (w + 1) * s.toNat : Nat) : Int) = 2 * ((2 : Int) ^ w * s) := by
      rw [Int.natCast_mul, Int.natCast_pow, hsN]
      show (2 : Int) ^ (w + 1) * s = _
      rw [hpow, Int.mul_assoc]
    rw [e1, e2]
    have : (2 * (2 : Int) ^ w - 3) * s = 2 * ((2 : Int) ^ w * s) - 3 * s := by
      rw [Int.sub_mul, Int.mul_assoc]
    rw [this] at hq
    omega
  have hsodd' : s.toNat % 2 = 1 := by omega
  obtain ⟨hv2, hodd⟩ := v2_oddPart_of_eq _ _ _ key hsodd'
  refine ⟨hv2, ?_⟩
  unfold drop syr
  rw [hodd, hsN, hAN]
  omega

/-- The admissibility condition of Theorem 2(b). -/
def Admissible (m : Int) (v : Nat) : Prop :=
  1 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * m + 1) ∧ 0 < (6 * m + 1) / ((2 : Int) ^ v - 3)

/-- **Theorem 2(b), the correspondence.** For odd `A`: `K(A) = m` iff `v_2(3A+1)` is
admissible for `m` and `A` is the unique preimage. Packaged: the admissible `v` are
exactly the valuations `v_2(3A + 1)` of the odd `A > 0` with `K(A) = m`. -/
theorem admissible_iff (m : Int) (v : Nat) :
    Admissible m v ↔ ∃ A : Nat, A % 2 = 1 ∧ 0 < A ∧ drop A = m ∧ v2 (3 * A + 1) = v := by
  constructor
  · rintro ⟨h1, h2, h3⟩
    obtain ⟨A, hA1, hA2, hA3, hA4⟩ := drop_backward m v h1 h2 h3
    exact ⟨A, hA1, hA2, hA4, hA3⟩
  · rintro ⟨A, hA1, _, rfl, rfl⟩
    obtain ⟨f1, f2, f3, f4⟩ := drop_forward A hA1
    refine ⟨f1, f2, ?_⟩
    rw [f3]
    exact Int.ofNat_lt.mpr f4

theorem pow_two_ge_four (v : Nat) (hv : 2 ≤ v) : (4 : Int) ≤ (2 : Int) ^ v := by
  obtain ⟨w, rfl⟩ : ∃ w, v = w + 2 := ⟨v - 2, by omega⟩
  rw [Int.pow_add]
  have : (0 : Int) < (2 : Int) ^ w := Int.pow_pos (by decide)
  have h4 : (2 : Int) ^ 2 = 4 := by decide
  rw [h4]
  omega

/-- For `m < 0`, the only admissible valuation is `v = 1` (one ascent of each size). -/
theorem admissible_of_neg (m : Int) (hm : m < 0) (v : Nat) : Admissible m v ↔ v = 1 := by
  constructor
  · rintro ⟨h1, ⟨q, hq⟩, h3⟩
    by_cases hv : v = 1
    · exact hv
    · have h4 := pow_two_ge_four v (by omega)
      rw [hq, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero v)] at h3
      have : 0 < ((2 : Int) ^ v - 3) * q := Int.mul_pos (by omega) h3
      omega
  · rintro rfl
    have h21 : (2 : Int) ^ 1 - 3 = -1 := by decide
    refine ⟨Nat.le_refl 1, ⟨-(6 * m + 1), by rw [h21]; omega⟩, ?_⟩
    rw [h21, Int.ediv_neg, Int.ediv_one]
    omega

/-- For `m ≥ 0`, the admissible valuations are the `v ≥ 2` with `(2^v - 3) ∣ 6m + 1`. -/
theorem admissible_of_nonneg (m : Int) (hm : 0 ≤ m) (v : Nat) :
    Admissible m v ↔ 2 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * m + 1) := by
  constructor
  · rintro ⟨h1, h2, h3⟩
    refine ⟨?_, h2⟩
    by_cases hv : v = 1
    · subst hv
      have h21 : (2 : Int) ^ 1 - 3 = -1 := by decide
      rw [h21, Int.ediv_neg, Int.ediv_one] at h3
      clear h2
      omega
    · omega
  · rintro ⟨h1, ⟨q, hq⟩⟩
    have h4 := pow_two_ge_four v h1
    refine ⟨by omega, ⟨q, hq⟩, ?_⟩
    rw [hq, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero v)]
    by_cases hq0 : 0 < q
    · exact hq0
    · have : ((2 : Int) ^ v - 3) * q ≤ 0 := Int.mul_nonpos_of_nonneg_of_nonpos (by omega) (by omega)
      omega

/-- Descents `A ≡ 1 (mod 4)` are exactly the odd `A` with `v_2(3A + 1) ≥ 2`. -/
theorem mod_four_iff (A : Nat) (hA : A % 2 = 1) : A % 4 = 1 ↔ 2 ≤ v2 (3 * A + 1) := by
  obtain ⟨h1, h2⟩ := syr_decomp A
  constructor
  · intro h4
    cases hv : v2 (3 * A + 1) with
    | zero =>
      rw [hv, Nat.pow_zero, Nat.one_mul] at h1
      omega
    | succ w =>
      cases w with
      | zero =>
        rw [hv, Nat.pow_one] at h1
        omega
      | succ _ => omega
  · intro hv
    obtain ⟨w, hw⟩ : ∃ w, v2 (3 * A + 1) = w + 2 := ⟨v2 (3 * A + 1) - 2, by omega⟩
    rw [hw, Nat.pow_add, Nat.mul_comm (2 ^ w) (2 ^ 2), Nat.mul_assoc] at h1
    have : (2 : Nat) ^ 2 = 4 := rfl
    rw [this] at h1
    omega

/-! ## The owner's label map -/

/-- The owner's map in labels `M = (A + 1) / 2`: `F(M) = (S(2M - 1) + 1) / 2`. -/
def labelF (M : Nat) : Nat := (syr (2 * M - 1) + 1) / 2

/-- **Theorem 1(b).** `F(2N) = 3N`. -/
theorem labelF_even (N : Nat) (hN : 1 ≤ N) : labelF (2 * N) = 3 * N := by
  unfold labelF syr
  have h : 3 * (2 * (2 * N) - 1) + 1 = 2 ^ 1 * (6 * N - 1) := by
    rw [Nat.pow_one]
    omega
  rw [h, (v2_oddPart_of_eq _ 1 (6 * N - 1) rfl (by omega)).2]
  omega

/-- **Theorem 1(b).** `F(4j + 1) = 3j + 1`. -/
theorem labelF_four_j_one (j : Nat) : labelF (4 * j + 1) = 3 * j + 1 := by
  unfold labelF syr
  have h : 3 * (2 * (4 * j + 1) - 1) + 1 = 2 ^ 2 * (6 * j + 1) := by
    show _ = 4 * (6 * j + 1)
    omega
  rw [h, (v2_oddPart_of_eq _ 2 (6 * j + 1) rfl (by omega)).2]
  omega

/-- **Theorem 1(b), the microcosm.** `F(4n - 1) = F(n)`. -/
theorem labelF_microcosm (n : Nat) (hn : 1 ≤ n) : labelF (4 * n - 1) = labelF n := by
  unfold labelF syr
  have h1 : 3 * (2 * (4 * n - 1) - 1) + 1 = 2 ^ 3 * (3 * n - 1) := by
    show _ = 8 * (3 * n - 1)
    omega
  have h2 : 3 * (2 * n - 1) + 1 = 2 ^ 1 * (3 * n - 1) := by
    rw [Nat.pow_one]
    omega
  rw [h1, h2, oddPart_two_pow_mul 3 _ (by omega), oddPart_two_pow_mul 1 _ (by omega)]

/-- **Theorem 1(a).** `F(M) = M - K(2M - 1)`; with `M = 2N - 1` this is the owner's
`F(2N - 1) = 2N - 1 - K_N`, `K_N = K(4N - 3)`. -/
theorem labelF_eq_sub_drop (M : Nat) (hM : 1 ≤ M) :
    (labelF M : Int) = (M : Int) - drop (2 * M - 1) := by
  obtain ⟨_, h2⟩ := syr_decomp (2 * M - 1)
  have hF : 2 * labelF M = syr (2 * M - 1) + 1 := by
    unfold labelF
    omega
  unfold drop
  have : ((2 * M - 1 : Nat) : Int) = 2 * (M : Int) - 1 := by omega
  rw [this]
  omega

/-- **Theorem 1(c).** `K_{2j+1} = j`, `K_{4m} = 5m - 1`, `K_{4m-2} = K_m + 6m - 4`, where
`K_N = K(4N - 3)`. -/
theorem drop_two_regular (j m : Nat) (hm : 1 ≤ m) :
    drop (4 * (2 * j + 1) - 3) = j ∧ drop (4 * (4 * m) - 3) = 5 * (m : Int) - 1 ∧
      drop (4 * (4 * m - 2) - 3) = drop (4 * m - 3) + 6 * (m : Int) - 4 := by
  refine ⟨?_, ?_, ?_⟩
  · unfold drop syr
    have h : 3 * (4 * (2 * j + 1) - 3) + 1 = 2 ^ 2 * (6 * j + 1) := by
      show _ = 4 * (6 * j + 1)
      omega
    rw [h, (v2_oddPart_of_eq _ 2 (6 * j + 1) rfl (by omega)).2]
    omega
  · unfold drop syr
    have h : 3 * (4 * (4 * m) - 3) + 1 = 2 ^ 3 * (6 * m - 1) := by
      show _ = 8 * (6 * m - 1)
      omega
    rw [h, (v2_oddPart_of_eq _ 3 (6 * m - 1) rfl (by omega)).2]
    omega
  · unfold drop syr
    have h1 : 3 * (4 * (4 * m - 2) - 3) + 1 = 2 ^ 4 * (3 * m - 2) := by
      show _ = 16 * (3 * m - 2)
      omega
    have h2 : 3 * (4 * m - 3) + 1 = 2 ^ 2 * (3 * m - 2) := by
      show _ = 4 * (3 * m - 2)
      omega
    rw [h1, h2, oddPart_two_pow_mul 4 _ (by omega), oddPart_two_pow_mul 2 _ (by omega)]
    have := (decomp (3 * m - 2) (by omega)).2
    omega

/-! ## Worked values from the note -/

/-- The owner's `K_N = K(4N - 3)` for `N = 1, …, 11`. -/
theorem owner_K_values :
    (List.range 11).map (fun i => drop (4 * (i + 1) - 3)) = [0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5] := by
  decide

/-- The owner's examples `S(1) = 1`, `S(5) = 1`, `S(9) = 7`, `S(13) = 5`. -/
theorem owner_S_values : syr 1 = 1 ∧ syr 5 = 1 ∧ syr 9 = 7 ∧ syr 13 = 5 := by
  decide

instance (m : Int) (v : Nat) : Decidable (Admissible m v) :=
  inferInstanceAs (Decidable (1 ≤ v ∧ _ ∧ _))

theorem two_pow_mono (k v : Nat) (h : k ≤ v) : (2 : Int) ^ k ≤ (2 : Int) ^ v := by
  have h1 : 2 ^ k ≤ 2 ^ v := Nat.pow_le_pow_right (by decide) h
  have h2 := Int.ofNat_le.mpr h1
  rw [Int.natCast_pow, Int.natCast_pow] at h2
  exact h2

/-- A positive divisor cannot exceed `6m + 1`. -/
theorem not_admissible_of_large (m : Int) (v : Nat) (hm : 0 ≤ m)
    (hv : 6 * m + 1 < (2 : Int) ^ v - 3) : ¬ Admissible m v := by
  rintro ⟨_, hd, _⟩
  have := Int.le_of_dvd (by omega) hd
  omega

/-- **Multiplicities** (the note's examples): `m = 1` occurs once (`v = 2`), `m = 2`
twice (`v = 2, 4`; `13 = 2^4 - 3`), `m = 24` three times (`v = 2, 3, 5`; `145 = 5 · 29`). -/
theorem multiplicity_examples :
    (∀ v, Admissible 1 v ↔ v = 2) ∧ (∀ v, Admissible 2 v ↔ v = 2 ∨ v = 4) ∧
      (∀ v, Admissible 24 v ↔ v = 2 ∨ v = 3 ∨ v = 5) := by
  have big : ∀ (m : Int) (v : Nat), 0 ≤ m → 8 ≤ v → 6 * m + 1 < 253 → ¬ Admissible m v := by
    intro m v hm hv hlt
    apply not_admissible_of_large m v hm
    have := two_pow_mono 8 v hv
    have h8 : (2 : Int) ^ 8 = 256 := by decide
    omega
  have small1 : ∀ v, v < 8 → (Admissible 1 v ↔ v = 2) := by decide
  have small2 : ∀ v, v < 8 → (Admissible 2 v ↔ v = 2 ∨ v = 4) := by decide
  have small24 : ∀ v, v < 8 → (Admissible 24 v ↔ v = 2 ∨ v = 3 ∨ v = 5) := by decide
  refine ⟨fun v => ?_, fun v => ?_, fun v => ?_⟩
  · by_cases hv : v < 8
    · exact small1 v hv
    · exact ⟨fun h => absurd h (big 1 v (by decide) (by omega) (by decide)), fun h => by omega⟩
  · by_cases hv : v < 8
    · exact small2 v hv
    · exact ⟨fun h => absurd h (big 2 v (by decide) (by omega) (by decide)),
        fun h => by rcases h with h | h <;> omega⟩
  · by_cases hv : v < 8
    · exact small24 v hv
    · exact ⟨fun h => absurd h (big 24 v (by decide) (by omega) (by decide)),
        fun h => by rcases h with h | h | h <;> omega⟩

end ProcgenSelfieEdim
