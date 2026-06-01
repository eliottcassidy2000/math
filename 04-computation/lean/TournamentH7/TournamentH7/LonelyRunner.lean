/-
  TournamentH7.LonelyRunner — Lean formalization of the Lonely Runner sieve theory.

  The Lonely Runner Conjecture (reduced form): for `k` nonzero integer speeds
  `v : ι → ℤ` (`|ι| = k`) and threshold `1/n` with `n = k+1`, there is a time
  `t ∈ ℝ` with `‖v i · t‖ ≥ 1/n` for every `i`, where `‖x‖` is the distance from
  `x` to the nearest integer.

  This module formalizes the *denominator-sieve* layer of the conjecture, which
  this repo had previously only verified computationally (Python, S356–S362,
  S388) and catalogued as THM-358 / THM-360 / THM-366.  This Lean module is
  catalogued as THM-369 (the machine-checked formalization of THM-366).
  The master result `sieve_frac` proves that
  for any reduced fraction `a/q` with `q ≤ n`, if no speed is divisible by `q`
  then `a/q` is a lonely time.  It uniformly subsumes:
    · THM-360  (unit endpoint divisibility filter) — the `q = n` case;
    · the positive direction of THM-358 (initial-segment unit witnesses);
    · the "all speeds odd ⇒ t = 1/2 is lonely" antipodal tool for even `n`;
    · the CRT folding tools `t = 1/p` for `p ∣ n` used at composite `n`.
  Consequently a counterexample must contain a speed divisible by every
  `q ∈ {2,…,n}` (`counterexample_needs_all_divisors`).

  oracle-2026-05-31-S18.
-/
import Mathlib.Tactic

namespace LonelyRunner

/-- `t` is `n`-lonely for the speed family `v` iff every `v i · t` is at distance
`≥ 1/n` from every integer.  This is the standard `‖v i · t‖ ≥ 1/n` condition on
the circle `ℝ/ℤ`, written in the proof-friendly "far from every integer" form. -/
def Lonely {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) : Prop :=
  ∀ i, ∀ m : ℤ, (1 : ℝ) / n ≤ |(v i : ℝ) * t - m|

/-- The observer has no *strict* LRC ties to the moving runners at time `t`.
Equivalently, no runner lies in the open forbidden observer window
`dist(v_i t, ℤ) < 1/n`. -/
def ObserverTieFree {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) : Prop :=
  ∀ i, ¬ ∃ m : ℤ, |(v i : ℝ) * t - m| < (1 : ℝ) / n

/-- Closed-threshold loneliness is exactly strict observer tie-freeness.  This
is the Lean core of the LRC trienerment convention: strict ties record open
nearness, while boundary contacts are still allowed LRC witnesses. -/
theorem lonely_iff_observerTieFree {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) :
    Lonely n v t ↔ ObserverTieFree n v t := by
  unfold Lonely ObserverTieFree
  constructor
  · intro h i htie
    rcases htie with ⟨m, hlt⟩
    exact not_lt_of_ge (h i m) hlt
  · intro h i m
    exact le_of_not_gt (fun hlt => h i ⟨m, hlt⟩)

section Sieve
variable {ι : Type*}

/-- **Master sieve lemma.**  Let `q ≤ n` with `0 < q`, let `a` be coprime to `q`,
and suppose no speed `v i` is divisible by `q`.  Then the reduced fraction `a/q`
is an `n`-lonely time. -/
theorem sieve_frac (n q : ℕ) (a : ℤ) (v : ι → ℤ)
    (hqn : q ≤ n) (hq0 : 0 < q) (hcop : IsCoprime (q : ℤ) a)
    (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i)) :
    Lonely n v ((a : ℝ) / q) := by
  intro i m
  have hq0' : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq0
  have hqne : (q : ℝ) ≠ 0 := ne_of_gt hq0'
  -- Put the value over the common denominator `q`.
  have key : (v i : ℝ) * ((a : ℝ) / q) - m
      = (((v i * a - m * q : ℤ)) : ℝ) / q := by
    field_simp
    push_cast
    ring
  rw [key, abs_div, abs_of_pos hq0']
  -- The numerator is a nonzero integer: `q ∤ v i · a`, since `q ∤ v i` and `q ⟂ a`.
  have hnotdvd : ¬ ((q : ℤ) ∣ v i * a) := fun h => hdiv i (hcop.dvd_of_dvd_mul_right h)
  have hne : (v i * a - m * q) ≠ 0 := by
    intro h
    apply hnotdvd
    have hva : v i * a = m * q := by linarith [sub_eq_zero.mp h]
    exact ⟨m, by rw [hva]; ring⟩
  -- A nonzero integer has absolute value `≥ 1`.
  have h1 : (1 : ℝ) ≤ |((v i * a - m * q : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]
    have : (1 : ℤ) ≤ |v i * a - m * q| := Int.one_le_abs hne
    exact_mod_cast this
  -- Chain `1/n ≤ 1/q ≤ |num|/q`.
  have hnq : (1 : ℝ) / n ≤ (1 : ℝ) / q :=
    one_div_le_one_div_of_le hq0' (by exact_mod_cast hqn)
  calc (1 : ℝ) / n ≤ (1 : ℝ) / q := hnq
    _ ≤ |((v i * a - m * q : ℤ) : ℝ)| / q := by gcongr

/-- The `a = 1` specialization: if no speed is divisible by `q` (`0 < q ≤ n`) then
`t = 1/q` is lonely.  This is the "denominator sieve". -/
theorem sieve_one_div (n q : ℕ) (v : ι → ℤ)
    (hqn : q ≤ n) (hq0 : 0 < q) (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i)) :
    Lonely n v ((1 : ℝ) / q) := by
  have h := sieve_frac n q 1 v hqn hq0 isCoprime_one_right hdiv
  simpa using h

/-- **Necessary condition for a counterexample.**  If `v` has *no* `n`-lonely
time, then for every `q` with `2 ≤ q ≤ n` some speed is divisible by `q`.  In
particular a counterexample must contain a multiple of `n` and of every prime
power `≤ n`. -/
theorem counterexample_needs_all_divisors (n : ℕ) (v : ι → ℤ)
    (hcex : ∀ t : ℝ, ¬ Lonely n v t)
    (q : ℕ) (hq2 : 2 ≤ q) (hqn : q ≤ n) :
    ∃ i, (q : ℤ) ∣ v i := by
  by_contra h
  push Not at h
  exact hcex ((1 : ℝ) / q) (sieve_one_div n q v hqn (by omega) h)

/-- **Antipodal tool (even-`n`).**  If every speed is odd and `2 ≤ n`, then
`t = 1/2` is lonely.  (Special case `q = 2`.) -/
theorem all_odd_half_lonely (n : ℕ) (v : ι → ℤ)
    (hn : 2 ≤ n) (hodd : ∀ i, ¬ (2 : ℤ) ∣ v i) :
    Lonely n v ((1 : ℝ) / 2) := by
  have h := sieve_one_div n 2 v hn (by norm_num) hodd
  simpa using h

end Sieve

/-- **Positive direction of THM-358.**  For the initial segment `v i = i+1`
(`i : Fin (n-1)`, i.e. speeds `1,…,n-1`) and any `a` coprime to `n`, the unit
fraction `a/n` is an `n`-lonely time.  Hence `{1,…,n-1}` is at best tight: it
always has the unit witnesses `a/n`. -/
theorem initial_segment_unit_lonely (n : ℕ) (hn : 1 ≤ n) (a : ℤ)
    (hcop : IsCoprime (n : ℤ) a) :
    Lonely n (fun i : Fin (n - 1) => (i : ℤ) + 1) ((a : ℝ) / n) := by
  apply sieve_frac n n a _ (le_refl n) (by omega) hcop
  intro i hdvd
  -- `n ∤ (i+1)` because `1 ≤ i+1 ≤ n-1 < n`.
  have hpos : (0 : ℤ) < (i : ℤ) + 1 := by positivity
  have hle : (n : ℤ) ≤ (i : ℤ) + 1 := Int.le_of_dvd hpos hdvd
  have hb : (i : ℕ) < n - 1 := i.isLt
  omega

/-! ### Grounding: `Lonely` is the standard central-box / nearest-integer condition

These connect the proof-friendly "far from every integer" form of `Lonely` to
`Int.fract`, hence to the **view-obstruction central box** `[1/n, 1-1/n]`, and
prove that loneliness is `1`-periodic in `t` for integer speeds. -/

section Grounding
variable {ι : Type*}

/-- `x` is at distance `≥ c` from every integer iff its fractional part lies in
the central band `[c, 1-c]`. -/
theorem far_iff_fract (c x : ℝ) :
    (∀ m : ℤ, c ≤ |x - m|) ↔ (c ≤ Int.fract x ∧ Int.fract x ≤ 1 - c) := by
  constructor
  · intro h
    refine ⟨?_, ?_⟩
    · have h0 := h ⌊x⌋
      have e : x - (⌊x⌋ : ℝ) = Int.fract x := by linarith [Int.floor_add_fract x]
      rw [e, abs_of_nonneg (Int.fract_nonneg x)] at h0
      exact h0
    · have h1 := h (⌊x⌋ + 1)
      have e : x - ((⌊x⌋ + 1 : ℤ) : ℝ) = Int.fract x - 1 := by
        push_cast; linarith [Int.floor_add_fract x]
      rw [e, abs_of_nonpos (by linarith [Int.fract_lt_one x] : Int.fract x - 1 ≤ 0)] at h1
      linarith
  · rintro ⟨h1, h2⟩ m
    have key : x - (m : ℝ) = Int.fract x - ((m - ⌊x⌋ : ℤ) : ℝ) := by
      push_cast; linarith [Int.floor_add_fract x]
    rw [key]
    obtain hk | hk : (m - ⌊x⌋ ≤ 0) ∨ (1 ≤ m - ⌊x⌋) := by omega
    · have hk' : ((m - ⌊x⌋ : ℤ) : ℝ) ≤ 0 := by exact_mod_cast hk
      exact le_abs.mpr (Or.inl (by linarith))
    · have hk' : (1 : ℝ) ≤ ((m - ⌊x⌋ : ℤ) : ℝ) := by exact_mod_cast hk
      exact le_abs.mpr (Or.inr (by linarith))

/-- `Lonely` ⟺ every `v i · t` has fractional part in the central box
`[1/n, 1-1/n]` — the view-obstruction formulation of the conjecture. -/
theorem lonely_iff_fract_mem (n : ℕ) (v : ι → ℤ) (t : ℝ) :
    Lonely n v t ↔
      ∀ i, (1 : ℝ) / n ≤ Int.fract ((v i : ℝ) * t)
            ∧ Int.fract ((v i : ℝ) * t) ≤ 1 - 1 / n := by
  unfold Lonely
  exact forall_congr' (fun i => far_iff_fract ((1 : ℝ) / n) ((v i : ℝ) * t))

/-- Loneliness for integer speeds is `1`-periodic in `t`. -/
theorem lonely_add_one (n : ℕ) (v : ι → ℤ) (t : ℝ) :
    Lonely n v (t + 1) ↔ Lonely n v t := by
  unfold Lonely
  constructor
  · intro h i m
    have hm := h i (m + v i)
    convert hm using 2
    push_cast; ring
  · intro h i m
    have hm := h i (m - v i)
    convert hm using 2
    push_cast; ring

end Grounding

/-! ### Proven LRC cases (as high `n` as the methodology gives unconditionally)

For *every* `n`, the `t = 1/n` (n-gon vertex) witness settles the case where no
speed is divisible by `n`; and `n = 2` is unconditional. -/

section Cases
variable {ι : Type*}

/-- **LRC for all `n`, no-multiple case.** If `0 < n` and no speed is divisible by
`n`, then `t = 1/n` is `n`-lonely. (The geodesic hits the `n`-gon vertex inside the
box.) -/
theorem lonely_of_no_multiple (n : ℕ) (hn : 0 < n) (v : ι → ℤ)
    (hdiv : ∀ i, ¬ ((n : ℤ) ∣ v i)) : Lonely n v ((1 : ℝ) / n) :=
  sieve_one_div n n v (le_refl n) hn hdiv

/-- **The 2-runner case (n = 2).** Any single nonzero speed has a 2-lonely time
(`t = 1/(2a)`, where `a·t = 1/2` is at distance `1/2` from every integer). -/
theorem lonely_two (a : ℤ) (ha : a ≠ 0) :
    ∃ t : ℝ, ∀ m : ℤ, (1 : ℝ) / 2 ≤ |(a : ℝ) * t - m| := by
  have ha' : (a : ℝ) ≠ 0 := by exact_mod_cast ha
  refine ⟨1 / (2 * a), fun m => ?_⟩
  have hval : (a : ℝ) * (1 / (2 * a)) = 1 / 2 := by field_simp
  rw [hval]
  have hf : Int.fract ((1 : ℝ) / 2) = 1 / 2 :=
    Int.fract_eq_self.mpr ⟨by norm_num, by norm_num⟩
  exact (far_iff_fract (1 / 2) (1 / 2)).mpr ⟨le_of_eq hf.symm, by rw [hf]; norm_num⟩ m

/-- **The 3-runner sieve coverage.**  For `n = 3`, the two unit witnesses `t = 1/3`
(`q = 3`) and `t = 1/2` (`q = 2`) settle every speed family *except* those in which
some speed is divisible by `3` AND some speed is divisible by `2`.  Concretely: if
either no speed is a multiple of `3`, or no speed is a multiple of `2`, then `v` has
a `3`-lonely time.  This isolates the genuine residual kernel of LRC@3 (speeds
entangled with `6`), which the denominator sieve alone cannot reach. -/
theorem three_lonely_sieve_cover (v : ι → ℤ)
    (h : (∀ i, ¬ (3 : ℤ) ∣ v i) ∨ (∀ i, ¬ (2 : ℤ) ∣ v i)) :
    ∃ t : ℝ, Lonely 3 v t := by
  rcases h with h3 | h2
  · exact ⟨(1 : ℝ) / 3, sieve_one_div 3 3 v (le_refl 3) (by norm_num) h3⟩
  · exact ⟨(1 : ℝ) / 2, sieve_one_div 3 2 v (by norm_num) (by norm_num) h2⟩

/-- **A checked residual family at `n = 3`.**  The sieve cover above cannot see
families containing a multiple of `3`; nevertheless each pair `{1, 3r}` has the
explicit witness `t = 1/3 + 1/(9r)`.  The first runner lies just inside the
central band and the second runner lands exactly at fractional part `1/3`. -/
theorem three_one_three_mul_lonely (r : ℕ) (hr : 0 < r) :
    ∃ t : ℝ, Lonely 3
      (fun i : Fin 2 => if i = 0 then (1 : ℤ) else (3 * (r : ℤ))) t := by
  let t : ℝ := (1 : ℝ) / 3 + 1 / (9 * (r : ℝ))
  refine ⟨t, ?_⟩
  rw [lonely_iff_fract_mem]
  intro i
  have hrR : (0 : ℝ) < (r : ℝ) := by exact_mod_cast hr
  have hrR1 : (1 : ℝ) ≤ (r : ℝ) := by exact_mod_cast hr
  have hden : (9 : ℝ) ≤ 9 * (r : ℝ) := by nlinarith
  have hsmall : (1 : ℝ) / (9 * (r : ℝ)) ≤ 1 / 9 := by
    exact one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 9) hden
  have ht_ge : (1 : ℝ) / 3 ≤ t := by
    dsimp [t]
    have hnonneg : (0 : ℝ) ≤ 1 / (9 * (r : ℝ)) := by positivity
    linarith
  have ht_le : t ≤ (2 : ℝ) / 3 := by
    dsimp [t]
    linarith
  have ht_nonneg : (0 : ℝ) ≤ t := by linarith
  have ht_lt_one : t < 1 := by linarith
  have hfract_t : Int.fract t = t := Int.fract_eq_self.mpr ⟨ht_nonneg, ht_lt_one⟩
  fin_cases i
  · norm_num
    rw [hfract_t]
    constructor
    · exact ht_ge
    · linarith
  · norm_num
    have hmul : (3 : ℝ) * (r : ℝ) * t = (r : ℝ) + (1 : ℝ) / 3 := by
      dsimp [t]
      field_simp [ne_of_gt hrR]
      ring_nf
    rw [hmul]
    have hfract : Int.fract ((r : ℝ) + (1 : ℝ) / 3) = (1 : ℝ) / 3 := by
      rw [Int.fract_natCast_add]
      norm_num
    rw [hfract]
    constructor <;> norm_num

/-- **Scaling invariance (repeated-addition reduction, S548).** Loneliness is
invariant under scaling all speeds by a nonzero `c` and the time by `1/c`, because
`(c·v_i)·(t/c) = v_i·t`.  Since a runner's position `v_i·t` is `t` *repeated-added*
`v_i` times, this is the statement that an arithmetic-progression family
`c·(v_i)` reduces to `(v_i)` at the scaled time `c·t`. -/
theorem lonely_scale (n : ℕ) (v : ι → ℤ) (t : ℝ) (c : ℤ) (hc : c ≠ 0) :
    Lonely n (fun i => c * v i) (t / (c : ℝ)) ↔ Lonely n v t := by
  have hc' : (c : ℝ) ≠ 0 := by exact_mod_cast hc
  unfold Lonely
  refine forall_congr' fun i => forall_congr' fun m => ?_
  have h : ((c * v i : ℤ) : ℝ) * (t / (c : ℝ)) = (v i : ℝ) * t := by
    push_cast; field_simp
  rw [h]

/-- **Doubled-prime / `n*=n/2` sieve (S546).** For the *doubled* dimension `n = 2p`
(`p > 0`), if no speed is divisible by `p` then `t = 1/p` (margin `2/n = 1/p`) is
lonely.  So a doubled-prime dimension `n = 2p` inherits the clean `mod p` channel via
`q = p ≤ n`. -/
theorem lonely_doubled (p : ℕ) (hp : 0 < p) (v : ι → ℤ)
    (hdiv : ∀ i, ¬ ((p : ℤ) ∣ v i)) : Lonely (2 * p) v ((1 : ℝ) / p) :=
  sieve_one_div (2 * p) p v (by omega) hp hdiv

/-- **Dirichlet box pigeonhole — "always a near pair" (S536/S539).** Among any
`n + 1` reals, two have fractional parts within `1/n`.  In the runner picture
(observer + runners) some pair is always within `1/n`, so the half-turn relation
always carries a *tie* (the trienerment is never a pure tournament) and some gap is
`≤ 1/n` (the apex pigeonhole). -/
theorem near_pair (n : ℕ) (hn : 0 < n) (x : Fin (n + 1) → ℝ) :
    ∃ i j : Fin (n + 1), i ≠ j ∧ |Int.fract (x i) - Int.fract (x j)| < 1 / n := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  -- box of a point = ⌊n · fract(x)⌋ (as a ℕ), landing in {0, …, n-1}; explicit (no `set`)
  have hmaps : ∀ i ∈ (Finset.univ : Finset (Fin (n + 1))),
      (⌊(n : ℝ) * Int.fract (x i)⌋).toNat ∈ Finset.range n := by
    intro i _
    have h0 : (0 : ℝ) ≤ (n : ℝ) * Int.fract (x i) :=
      mul_nonneg (le_of_lt hnR) (Int.fract_nonneg _)
    have h1 : (n : ℝ) * Int.fract (x i) < n := by
      nlinarith [Int.fract_lt_one (x i), Int.fract_nonneg (x i)]
    have hge : 0 ≤ ⌊(n : ℝ) * Int.fract (x i)⌋ := Int.floor_nonneg.mpr h0
    have hlt : ⌊(n : ℝ) * Int.fract (x i)⌋ < (n : ℤ) := by
      rw [Int.floor_lt]; exact_mod_cast h1
    rw [Finset.mem_range]
    omega
  have hcard : (Finset.range n).card < (Finset.univ : Finset (Fin (n + 1))).card := by
    simp
  obtain ⟨i, _, j, _, hij, hbox_eq⟩ :=
    Finset.exists_ne_map_eq_of_card_lt_of_maps_to hcard hmaps
  refine ⟨i, j, hij, ?_⟩
  -- hbox_eq : (⌊n·fract i⌋).toNat = (⌊n·fract j⌋).toNat; with both floors ≥ 0 ⇒ floors equal
  have hgei : 0 ≤ ⌊(n : ℝ) * Int.fract (x i)⌋ :=
    Int.floor_nonneg.mpr (mul_nonneg (le_of_lt hnR) (Int.fract_nonneg _))
  have hgej : 0 ≤ ⌊(n : ℝ) * Int.fract (x j)⌋ :=
    Int.floor_nonneg.mpr (mul_nonneg (le_of_lt hnR) (Int.fract_nonneg _))
  have hfloor : ⌊(n : ℝ) * Int.fract (x i)⌋ = ⌊(n : ℝ) * Int.fract (x j)⌋ := by
    omega
  have hlt1 : |(n : ℝ) * Int.fract (x i) - (n : ℝ) * Int.fract (x j)| < 1 :=
    Int.abs_sub_lt_one_of_floor_eq_floor hfloor
  have hexp : (n : ℝ) * Int.fract (x i) - (n : ℝ) * Int.fract (x j)
      = (n : ℝ) * (Int.fract (x i) - Int.fract (x j)) := by ring
  rw [hexp, abs_mul, abs_of_pos hnR] at hlt1
  -- hlt1 : n * |fract i - fract j| < 1
  rw [lt_div_iff₀ hnR]
  linarith [hlt1]

end Cases

/-! ### Axiom audit

Each `#print axioms` should report only the foundational
`[propext, Classical.choice, Quot.sound]` — i.e. no `sorry` and no
project-specific axioms underlie the LRC sieve theory. -/

#print axioms sieve_frac
#print axioms lonely_iff_observerTieFree
#print axioms sieve_one_div
#print axioms counterexample_needs_all_divisors
#print axioms all_odd_half_lonely
#print axioms initial_segment_unit_lonely
#print axioms far_iff_fract
#print axioms lonely_iff_fract_mem
#print axioms lonely_add_one
#print axioms lonely_of_no_multiple
#print axioms lonely_two
#print axioms three_lonely_sieve_cover
#print axioms three_one_three_mul_lonely
#print axioms lonely_scale
#print axioms lonely_doubled
#print axioms near_pair

end LonelyRunner
