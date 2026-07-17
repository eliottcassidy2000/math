/-
  TournamentH7.LRCRationalLock — THE RATIONAL-RATIO LOCK AND ITS EXACT COUNT
  (death-star-2026-07-17-S48, HYP-7240; tandem with boxeph LEM-042/044).

  THM-966's integer rung lock generalizes to REDUCED RATIONAL ratios: speeds
  `g·i'` and `g·j'` with the relation `j'·v_a = i'·v_b`.

  * `rational_lock` — if `i' + j' ≤ 13`, joint failure locks the witnesses on
    the Bézout ray: `j'·w_a = i'·w_b` EXACTLY (the exact identity
    `j'(v_a p − w_a q) − i'(v_b p − w_b q) = (i' w_b − j' w_a)·q` plus the
    strict triangle `< (i'+j')q/14 ≤ 13q/14`).  THM-966 is the case `i' = 1`.
  * `rational_branch_bound` — if `i' + j' ≤ 27`, at most THREE witness
    branches: `|j'·w_a − i'·w_b| ≤ 1`.  All 78 pairs of the canonical family
    `{1,…,13}` have `i' + j' ≤ 25`: the whole pair layer is ≤3-branch.
  * `witness_mod_bridge` — the two-sided normal form: for `0 ≤ r < q`,
    `∃s, 14M·|x − s·q| < q ⟺ 14M·r < q ∨ 14M·(q−r) < q` (where `r = x mod q`)
    — only `s ∈ {⌊x/q⌋, ⌊x/q⌋+1}` can win, giving the two ℕ-arcs that
    THM-961's `narrowFailN_count` counts.
  * `rational_pair_fail_iff` — locked case (`gcd(i',j') = 1`): both fail ⟺
    the GCD-speed `g` fails the `14·max(i',j')`-narrow band (coprimality
    turns the Bézout ray into `w_a = i'·s, w_b = j'·s`).
  * `rational_pair_count` — composed with THM-961's mod transport: at
    `gcd(g, q) = 1` the joint-failure count is EXACTLY
    `2·⌊(q−1)/(14·max(i',j'))⌋`.

  THE DISCRETE–CONTINUOUS BRIDGE (recon `rationallock_recon_deathstar_S48.out`):
  `N/(q−1) → 1/(7·max(i',j'))`, reproducing boxeph LEM-044's
  `μ(D_k ∩ D_{k+1}) = 1/49 + r(6−r)/(49k(k+1))` exactly for `k ≤ 6`
  (1/21 at k=2, 1/28 at k=3, …), with their zero-excess-iff-`7∣k` marking the
  sparse boundary — the same pair law seen from both sides of the bridge, now
  with floors valid at EVERY finite q (which is what the census consumes).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCLockedChainCount

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **THE RATIONAL-RATIO LOCK**: under the reduced relation `j'·v_a = i'·v_b`
with `i' + j' ≤ 13`, joint failure locks `j'·w_a = i'·w_b` exactly. -/
theorem rational_lock (va vb wa wb : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 13)
    (hrel : (j' : ℤ) * va = (i' : ℤ) * vb)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q) :
    (j' : ℤ) * wa = (i' : ℤ) * wb := by
  have hiZ : (1 : ℤ) ≤ (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (1 : ℤ) ≤ (j' : ℤ) := by exact_mod_cast hj
  have hsumZ : (i' : ℤ) + (j' : ℤ) ≤ 13 := by exact_mod_cast hsum
  -- the exact identity
  have hkey : ((i' : ℤ) * wb - (j' : ℤ) * wa) * q
      = (j' : ℤ) * (va * (p : ℤ) - wa * q) - (i' : ℤ) * (vb * (p : ℤ) - wb * q) := by
    have : (j' : ℤ) * (va * (p : ℤ)) = (i' : ℤ) * (vb * (p : ℤ)) := by
      calc (j' : ℤ) * (va * (p : ℤ)) = ((j' : ℤ) * va) * (p : ℤ) := by ring
        _ = ((i' : ℤ) * vb) * (p : ℤ) := by rw [hrel]
        _ = (i' : ℤ) * (vb * (p : ℤ)) := by ring
    linear_combination -this
  -- triangle with the strict band bounds
  have habs : |((i' : ℤ) * wb - (j' : ℤ) * wa) * q|
      ≤ (j' : ℤ) * |va * (p : ℤ) - wa * q| + (i' : ℤ) * |vb * (p : ℤ) - wb * q| := by
    rw [hkey]
    calc |(j' : ℤ) * (va * (p : ℤ) - wa * q) - (i' : ℤ) * (vb * (p : ℤ) - wb * q)|
        ≤ |(j' : ℤ) * (va * (p : ℤ) - wa * q)| + |(i' : ℤ) * (vb * (p : ℤ) - wb * q)| :=
          abs_sub _ _
      _ = (j' : ℤ) * |va * (p : ℤ) - wa * q| + (i' : ℤ) * |vb * (p : ℤ) - wb * q| := by
          rw [abs_mul, abs_mul, abs_of_pos (by linarith : (0:ℤ) < (j' : ℤ)),
            abs_of_pos (by linarith : (0:ℤ) < (i' : ℤ))]
  have hscale : 14 * |((i' : ℤ) * wb - (j' : ℤ) * wa) * q| < ((i' : ℤ) + j') * q := by
    have h1 : (j' : ℤ) * (14 * |va * (p : ℤ) - wa * q|) < (j' : ℤ) * q :=
      mul_lt_mul_of_pos_left ha (by linarith)
    have h2 : (i' : ℤ) * (14 * |vb * (p : ℤ) - wb * q|) < (i' : ℤ) * q :=
      mul_lt_mul_of_pos_left hb (by linarith)
    nlinarith [habs]
  -- |k|·q < 14q with k integer and the modulus factor positive ⟹ k = 0
  have hfact : |((i' : ℤ) * wb - (j' : ℤ) * wa) * q|
      = |(i' : ℤ) * wb - (j' : ℤ) * wa| * q := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))]
  rw [hfact] at hscale
  set k : ℤ := (i' : ℤ) * wb - (j' : ℤ) * wa with hk
  have hk0 : k = 0 := by
    by_contra hne
    have h1 : 1 ≤ |k| := Int.one_le_abs hne
    have hq1 : 1 ≤ (q : ℤ) := by
      rcases Nat.eq_zero_or_pos q with rfl | hq
      · simp at hscale
      · exact_mod_cast hq
    nlinarith [hscale, hsumZ]
  have : (i' : ℤ) * wb = (j' : ℤ) * wa := by linarith [hk]
  linarith [this]

/-- **The branch bound**: for `i' + j' ≤ 27` there are at most three witness
branches, `|j'·w_a − i'·w_b| ≤ 1`.  Every pair of `{1,…,13}` qualifies. -/
theorem rational_branch_bound (va vb wa wb : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 27) (hq : 0 < q)
    (hrel : (j' : ℤ) * va = (i' : ℤ) * vb)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q) :
    |(j' : ℤ) * wa - (i' : ℤ) * wb| ≤ 1 := by
  have hiZ : (1 : ℤ) ≤ (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (1 : ℤ) ≤ (j' : ℤ) := by exact_mod_cast hj
  have hsumZ : (i' : ℤ) + (j' : ℤ) ≤ 27 := by exact_mod_cast hsum
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hkey : ((i' : ℤ) * wb - (j' : ℤ) * wa) * q
      = (j' : ℤ) * (va * (p : ℤ) - wa * q) - (i' : ℤ) * (vb * (p : ℤ) - wb * q) := by
    have : (j' : ℤ) * (va * (p : ℤ)) = (i' : ℤ) * (vb * (p : ℤ)) := by
      calc (j' : ℤ) * (va * (p : ℤ)) = ((j' : ℤ) * va) * (p : ℤ) := by ring
        _ = ((i' : ℤ) * vb) * (p : ℤ) := by rw [hrel]
        _ = (i' : ℤ) * (vb * (p : ℤ)) := by ring
    linear_combination -this
  have habs : |((i' : ℤ) * wb - (j' : ℤ) * wa) * q|
      ≤ (j' : ℤ) * |va * (p : ℤ) - wa * q| + (i' : ℤ) * |vb * (p : ℤ) - wb * q| := by
    rw [hkey]
    calc |(j' : ℤ) * (va * (p : ℤ) - wa * q) - (i' : ℤ) * (vb * (p : ℤ) - wb * q)|
        ≤ |(j' : ℤ) * (va * (p : ℤ) - wa * q)| + |(i' : ℤ) * (vb * (p : ℤ) - wb * q)| :=
          abs_sub _ _
      _ = (j' : ℤ) * |va * (p : ℤ) - wa * q| + (i' : ℤ) * |vb * (p : ℤ) - wb * q| := by
          rw [abs_mul, abs_mul, abs_of_pos (by linarith : (0:ℤ) < (j' : ℤ)),
            abs_of_pos (by linarith : (0:ℤ) < (i' : ℤ))]
  have hscale : 14 * (|(i' : ℤ) * wb - (j' : ℤ) * wa| * q) < ((i' : ℤ) + j') * q := by
    have h1 : (j' : ℤ) * (14 * |va * (p : ℤ) - wa * q|) < (j' : ℤ) * q :=
      mul_lt_mul_of_pos_left ha (by linarith)
    have h2 : (i' : ℤ) * (14 * |vb * (p : ℤ) - wb * q|) < (i' : ℤ) * q :=
      mul_lt_mul_of_pos_left hb (by linarith)
    have hfact : |((i' : ℤ) * wb - (j' : ℤ) * wa) * q|
        = |(i' : ℤ) * wb - (j' : ℤ) * wa| * q := by
      rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))]
    nlinarith [habs, hfact.symm.le, hfact.le]
  -- 14·|k|·q < 27·q ⟹ |k| ≤ 1
  by_contra hcon
  push Not at hcon
  have h2 : 2 ≤ |(i' : ℤ) * wb - (j' : ℤ) * wa| := by
    have : |(j' : ℤ) * wa - (i' : ℤ) * wb| = |(i' : ℤ) * wb - (j' : ℤ) * wa| :=
      abs_sub_comm _ _
    omega
  nlinarith [hscale, hsumZ, hqZ, h2]

/-- **The witness–mod bridge**: for `x` with residue `r = (x mod q)`, the
witness form of the `14M`-narrow band equals the two-arc residue form.  Only
`s ∈ {⌊x/q⌋, ⌊x/q⌋+1}` can win the strict bound. -/
theorem witness_mod_bridge (x : ℤ) (M q : ℕ) (hM : 1 ≤ M) (hq : 0 < q) :
    (∃ s : ℤ, 14 * M * |x - s * q| < q)
      ↔ (14 * M * ((x % (q : ℤ)).toNat) < q ∨
         14 * M * (q - ((x % (q : ℤ)).toNat)) < q) := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hr_nonneg : (0 : ℤ) ≤ x % (q : ℤ) := Int.emod_nonneg _ (by omega)
  have hr_lt : x % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
  have hrval : (((x % (q : ℤ)).toNat : ℕ) : ℤ) = x % (q : ℤ) :=
    Int.toNat_of_nonneg hr_nonneg
  have hdm := Int.ediv_add_emod x (q : ℤ)
  set r : ℤ := x % (q : ℤ) with hr
  set d : ℤ := x / (q : ℤ) with hd
  have hMZ : (1 : ℤ) ≤ (M : ℤ) := by exact_mod_cast hM
  constructor
  · rintro ⟨s, hs⟩
    -- x − s·q = r + (d − s)·q; the strict bound forces d − s ∈ {0, −1}
    have hx : x - s * q = r + (d - s) * q := by
      have : (q : ℤ) * d + r = x := hdm
      linarith [this]
    rw [hx] at hs
    have habs : |r + (d - s) * q| < q := by
      have h14 : (14 : ℤ) * M ≥ 14 := by linarith
      nlinarith [abs_nonneg (r + (d - s) * q), hs]
    have hds : d - s = 0 ∨ d - s = -1 := by
      rcases le_or_gt 1 (d - s) with h | h
      · exfalso
        have : r + (d - s) * q ≥ q := by nlinarith [hr_nonneg]
        rcases abs_cases (r + (d - s) * q) with ⟨he, _⟩ | ⟨he, _⟩ <;> omega
      · rcases le_or_gt (d - s) (-2) with h2 | h2
        · exfalso
          have : r + (d - s) * q ≤ -q := by nlinarith [hr_lt]
          rcases abs_cases (r + (d - s) * q) with ⟨he, _⟩ | ⟨he, _⟩ <;> omega
        · omega
    rcases hds with h0 | h1
    · left
      rw [h0] at hs
      simp only [zero_mul, add_zero] at hs
      rw [abs_of_nonneg hr_nonneg] at hs
      have : 14 * (M : ℤ) * ((((x % (q : ℤ)).toNat : ℕ)) : ℤ) < q := by
        rw [hrval]; exact hs
      exact_mod_cast this
    · right
      rw [h1] at hs
      have hval : r + (-1) * q = -(q - r) := by ring
      rw [hval, abs_neg, abs_of_nonneg (by linarith : (0:ℤ) ≤ (q : ℤ) - r)] at hs
      have hcast : ((q - (x % (q : ℤ)).toNat : ℕ) : ℤ) = (q : ℤ) - r := by
        have hrq : (x % (q : ℤ)).toNat ≤ q := by omega
        push_cast [hrq]
        rw [hrval]
      have : 14 * (M : ℤ) * (((q - (x % (q : ℤ)).toNat : ℕ)) : ℤ) < q := by
        rw [hcast]; linarith [hs]
      exact_mod_cast this
  · rintro (h | h)
    · refine ⟨d, ?_⟩
      have hx : x - d * q = r := by linarith [hdm]
      rw [hx, abs_of_nonneg hr_nonneg]
      have hZ : 14 * (M : ℤ) * r < q := by
        have := h
        have hcast : (14 * M * ((x % (q : ℤ)).toNat) : ℕ) < (q : ℕ) := h
        have : ((14 * M * ((x % (q : ℤ)).toNat) : ℕ) : ℤ) < ((q : ℕ) : ℤ) := by
          exact_mod_cast hcast
        push_cast at this
        rw [hrval] at this
        linarith [this]
      exact hZ
    · refine ⟨d + 1, ?_⟩
      have hx : x - (d + 1) * q = r - q := by linarith [hdm]
      rw [hx, abs_of_nonpos (by linarith : r - (q : ℤ) ≤ 0)]
      have hcast : ((q - (x % (q : ℤ)).toNat : ℕ) : ℤ) = (q : ℤ) - r := by
        have hrq : (x % (q : ℤ)).toNat ≤ q := by omega
        push_cast [hrq]
        rw [hrval]
      have hZ : ((14 * M * (q - (x % (q : ℤ)).toNat) : ℕ) : ℤ) < ((q : ℕ) : ℤ) := by
        exact_mod_cast h
      have hZ' : 14 * (M : ℤ) * ((q - (x % (q : ℤ)).toNat : ℕ) : ℤ) < (q : ℤ) := by
        push_cast at hZ ⊢
        linarith [hZ]
      rw [hcast] at hZ'
      have hneg : -(r - (q : ℤ)) = (q : ℤ) - r := by ring
      rw [hneg]
      linarith [hZ']

/-- **The locked rational pair, narrow form**: with `gcd(i', j') = 1` and
`i' + j' ≤ 13`, joint failure of `g·i'` and `g·j'` is EXACTLY the GCD-speed
narrow band at `M = max i' j'`. -/
theorem rational_pair_fail_iff (g : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 13)
    (hcop : Nat.Coprime i' j') :
    ((∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
     (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q))
      ↔ (∃ s : ℤ, 14 * (max i' j' : ℕ) * |g * (p : ℤ) - s * q| < q) := by
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (0 : ℤ) < (j' : ℤ) := by exact_mod_cast hj
  constructor
  · rintro ⟨⟨wa, ha⟩, ⟨wb, hb⟩⟩
    have hrel : (j' : ℤ) * (g * i') = (i' : ℤ) * (g * j') := by ring
    have hlock : (j' : ℤ) * wa = (i' : ℤ) * wb :=
      rational_lock (g * i') (g * j') wa wb i' j' q p hi hj hsum hrel ha hb
    -- coprimality splits the Bézout ray: wa = i'·s
    have hdvd : (i' : ℤ) ∣ (j' : ℤ) * wa := ⟨wb, hlock⟩
    have hcopZ : IsCoprime (i' : ℤ) (j' : ℤ) := by
      rw [Int.isCoprime_iff_gcd_eq_one]
      simpa using hcop
    obtain ⟨s, hs⟩ := hcopZ.dvd_of_dvd_mul_left hdvd
    refine ⟨s, ?_⟩
    -- |g·p − s·q| = |X_a| / i' via X_a = i'·(g p − s q)
    have hXa : (g * i') * (p : ℤ) - wa * q = (i' : ℤ) * (g * p - s * q) := by
      rw [hs]; ring
    rw [hXa, abs_mul, abs_of_pos hiZ] at ha
    -- and the b-side: wb = j'·s from the lock
    have hwb : wb = (j' : ℤ) * s := by
      have h1 : (i' : ℤ) * wb = (i' : ℤ) * ((j' : ℤ) * s) := by
        calc (i' : ℤ) * wb = (j' : ℤ) * wa := hlock.symm
          _ = (j' : ℤ) * ((i' : ℤ) * s) := by rw [hs]
          _ = (i' : ℤ) * ((j' : ℤ) * s) := by ring
      exact mul_left_cancel₀ (by linarith) h1
    have hXb : (g * j') * (p : ℤ) - wb * q = (j' : ℤ) * (g * p - s * q) := by
      rw [hwb]; ring
    rw [hXb, abs_mul, abs_of_pos hjZ] at hb
    -- the max side follows from whichever of a/b carries it
    rcases max_cases i' j' with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm]
    · calc 14 * (i' : ℤ) * |g * (p : ℤ) - s * q|
          = 14 * ((i' : ℤ) * |g * (p : ℤ) - s * q|) := by ring
        _ < q := ha
    · calc 14 * (j' : ℤ) * |g * (p : ℤ) - s * q|
          = 14 * ((j' : ℤ) * |g * (p : ℤ) - s * q|) := by ring
        _ < q := hb
  · rintro ⟨s, hs⟩
    have habs : (0 : ℤ) ≤ |g * (p : ℤ) - s * q| := abs_nonneg _
    have hmax_i : (i' : ℤ) ≤ (max i' j' : ℕ) := by
      have := le_max_left i' j'
      exact_mod_cast this
    have hmax_j : (j' : ℤ) ≤ (max i' j' : ℕ) := by
      have := le_max_right i' j'
      exact_mod_cast this
    constructor
    · refine ⟨(i' : ℤ) * s, ?_⟩
      have hXa : (g * i') * (p : ℤ) - ((i' : ℤ) * s) * q
          = (i' : ℤ) * (g * p - s * q) := by ring
      rw [hXa, abs_mul, abs_of_pos hiZ]
      nlinarith [hs, habs]
    · refine ⟨(j' : ℤ) * s, ?_⟩
      have hXb : (g * j') * (p : ℤ) - ((j' : ℤ) * s) * q
          = (j' : ℤ) * (g * p - s * q) := by ring
      rw [hXb, abs_mul, abs_of_pos hjZ]
      nlinarith [hs, habs]

open Classical in
/-- **THE RATIONAL PAIR COUNT**: at `gcd(g, q) = 1`, the joint-failure count
of the reduced-ratio pair `(g·i', g·j')` with `i' + j' ≤ 13` is EXACTLY
`2·⌊(q−1)/(14·max(i',j'))⌋` — the discrete face of LEM-044's μ-law. -/
theorem rational_pair_count (g : ℤ) (i' j' : ℕ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 13)
    (hcop : Nat.Coprime i' j') (hgcd : Int.gcd g (q : ℤ) = 1) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
        (∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
        (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q)).card
      = 2 * ((q - 1) / (14 * max i' j')) := by
  classical
  have hM : 1 ≤ max i' j' := le_trans hi (le_max_left i' j')
  -- rewrite the filter through the iff and the bridge into residue form
  have hfilter : ((Finset.Ioo 0 q).filter fun p : ℕ =>
      (∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
      (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q))
      = ((Finset.Ioo 0 q).filter fun p : ℕ =>
          14 * (max i' j') * ((g * (p : ℤ) % (q : ℤ)).toNat) < q ∨
          14 * (max i' j') * (q - ((g * (p : ℤ) % (q : ℤ)).toNat)) < q) := by
    apply Finset.filter_congr
    intro p _
    rw [rational_pair_fail_iff g i' j' q p hi hj hsum hcop,
      witness_mod_bridge (g * (p : ℤ)) (max i' j') q hM hq]
  rw [hfilter]
  have htrans := card_mod_filter_eq g q hq hgcd
    (fun r => 14 * (max i' j') * r < q ∨ 14 * (max i' j') * (q - r) < q)
  rw [htrans]
  exact narrowFailN_count (max i' j') q hM hq

/-! ## Axiom audit -/
#print axioms rational_lock
#print axioms rational_branch_bound
#print axioms witness_mod_bridge
#print axioms rational_pair_fail_iff
#print axioms rational_pair_count

end LRC14Concrete
end LonelyRunner
