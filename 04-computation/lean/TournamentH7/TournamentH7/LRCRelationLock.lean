/-
  TournamentH7.LRCRelationLock — THE RELATION LOCK BY COEFFICIENT WEIGHT AND
  THE MEDIANT TRIPLE COUNT (death-star-2026-07-17-S51, HYP-7260; the owner's
  projective hint decoded).

  THE PRINCIPLE: witnesses inherit every light integer relation of the speeds.
  If `Σ αᵢ·vᵢ = 0` with COEFFICIENT weight `Σ|αᵢ| ≤ 14`, and every involved
  runner fails the safe band at `p/q`, then `Σ αᵢ·wᵢ = 0` EXACTLY — the
  identity `Σ αᵢ·Xᵢ = −(Σ αᵢwᵢ)·q` plus the per-term integer gap
  `14|Xᵢ| ≤ q−1` (no strictness gymnastics: `14·|Σαw|·q ≤ (Σ|α|)(q−1) ≤
  14(q−1) < 14q` kills any nonzero value).

  * `relation_lock` — the general Finset form (any index set, any ℤ weights).
  * `relation_lock3` — the three-term workhorse.
  * `sum_triple_lock` — sum-triples `{a, b, a+b}` have weight 3: ALWAYS
    locked, `w_c = w_a + w_b` — including triples whose PAIRS are sparse
    (e.g. `(5,6,11)`: pairs at weight 16, 17 — yet the triple locks).
  * `rational_lock_weight14` — the corrected pair boundary: the ratio
    relation `j'v_a − i'v_b = 0` has weight `i'+j'`, so pairs lock through
    weight 14 (extends THM-967's ≤ 13; the weight-14 pairs of {1..13} —
    (1,13), (3,11), (5,9) — have EMPTY Bézout branches, recon (b)).
  * `mediant_triple_fail_iff` / `mediant_triple_count` — the projective
    object: the mediant chain `(g·i', g·j', g·(i'+j'))` (a Farey/Stern–Brocot
    line; the sum-lines of {1..7} form the Fano incidence whose Levi 6-cycles
    are line-triangles).  All three fail ⟺ the gcd-speed fails the
    `14·(i'+j')`-narrow band; count `2·⌊(q−1)/(14·(i'+j'))⌋` at coprime
    moduli.  The TRIPLE layer of the {1..13} ledger opens: recon (d) exact on
    all 9 test triples (`N = 2⌊(q−1)/(14·c/g)⌋`).

  Recon (`relationlock_recon_deathstar_S51.out`): relation lock 4471/4471;
  weight-14 pairs branch-empty 3/3; sum-triple lock 4135/4135 over all 36
  sum-triples of {1..13}; mediant counts exact 9/9.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCBranchInterval

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **THE RELATION LOCK (general form)**: witnesses inherit every vanishing
integer relation of the speeds whose coefficient weight is ≤ 14. -/
theorem relation_lock {ι : Type*} (s : Finset ι) (α v w : ι → ℤ) (q p : ℕ)
    (hrel : ∑ i ∈ s, α i * v i = 0)
    (hwt : ∑ i ∈ s, |α i| ≤ 14)
    (hband : ∀ i ∈ s, 14 * |v i * (p : ℤ) - w i * q| < q) :
    ∑ i ∈ s, α i * w i = 0 := by
  rcases Nat.eq_zero_or_pos q with rfl | hq
  · -- q = 0: each band bound reads 14|…| < 0, impossible unless s empty-ish
    rcases Finset.eq_empty_or_nonempty s with rfl | ⟨i, hi⟩
    · simp
    · have h := hband i hi
      have habs : (0 : ℤ) ≤ |v i * (p : ℤ) - w i * 0| := abs_nonneg _
      push_cast at h
      linarith
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- the exact identity: Σ αX = −(Σ αw)·q
  have hkey : (∑ i ∈ s, α i * w i) * (q : ℤ)
      = -(∑ i ∈ s, α i * (v i * (p : ℤ) - w i * q)) := by
    have hexp : ∑ i ∈ s, α i * (v i * (p : ℤ) - w i * q)
        = (∑ i ∈ s, α i * v i) * (p : ℤ) - (∑ i ∈ s, α i * w i) * q := by
      rw [Finset.sum_mul, Finset.sum_mul]
      rw [← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro i _
      ring
    rw [hexp, hrel]
    ring
  -- the ≤-chain with the integer gap
  have hterm : ∀ i ∈ s, 14 * |v i * (p : ℤ) - w i * q| ≤ (q : ℤ) - 1 := by
    intro i hi
    have h := hband i hi
    omega
  have hsum : 14 * |∑ i ∈ s, α i * (v i * (p : ℤ) - w i * q)|
      ≤ (∑ i ∈ s, |α i|) * ((q : ℤ) - 1) := by
    calc 14 * |∑ i ∈ s, α i * (v i * (p : ℤ) - w i * q)|
        ≤ 14 * ∑ i ∈ s, |α i * (v i * (p : ℤ) - w i * q)| := by
          have := Finset.abs_sum_le_sum_abs
            (fun i => α i * (v i * (p : ℤ) - w i * q)) s
          linarith
      _ = ∑ i ∈ s, |α i| * (14 * |v i * (p : ℤ) - w i * q|) := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro i _
          rw [abs_mul]
          ring
      _ ≤ ∑ i ∈ s, |α i| * ((q : ℤ) - 1) := by
          apply Finset.sum_le_sum
          intro i hi
          exact mul_le_mul_of_nonneg_left (hterm i hi) (abs_nonneg _)
      _ = (∑ i ∈ s, |α i|) * ((q : ℤ) - 1) := by
          rw [Finset.sum_mul]
  by_contra hne
  have h1 : 1 ≤ |∑ i ∈ s, α i * w i| := Int.one_le_abs hne
  have h2 : (q : ℤ) ≤ |(∑ i ∈ s, α i * w i) * (q : ℤ)| := by
    rw [abs_mul, abs_of_pos hqZ]
    nlinarith [h1]
  rw [hkey, abs_neg] at h2
  have hwtZ : (∑ i ∈ s, |α i|) ≤ (14 : ℤ) := hwt
  have hq1 : (0 : ℤ) ≤ (q : ℤ) - 1 := by omega
  nlinarith [hsum, h2, hwtZ, hq1]

/-- **The three-term relation lock** — the workhorse form (direct proof). -/
theorem relation_lock3 (va vb vc wa wb wc α β γ : ℤ) (q p : ℕ)
    (hrel : α * va + β * vb + γ * vc = 0)
    (hwt : |α| + |β| + |γ| ≤ 14)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q)
    (hc : 14 * |vc * (p : ℤ) - wc * q| < q) :
    α * wa + β * wb + γ * wc = 0 := by
  rcases Nat.eq_zero_or_pos q with rfl | hq
  · exfalso
    push_cast at ha
    linarith [abs_nonneg (va * (p : ℤ) - wa * 0)]
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hkey : (α * wa + β * wb + γ * wc) * (q : ℤ)
      = -(α * (va * (p : ℤ) - wa * q) + β * (vb * (p : ℤ) - wb * q)
          + γ * (vc * (p : ℤ) - wc * q)) := by
    linear_combination (p : ℤ) * hrel
  have habs_add : ∀ x y : ℤ, |x + y| ≤ |x| + |y| := by
    intro x y
    have h := abs_sub x (-y)
    rw [sub_neg_eq_add, abs_neg] at h
    exact h
  have htri : |α * (va * (p : ℤ) - wa * q) + β * (vb * (p : ℤ) - wb * q)
      + γ * (vc * (p : ℤ) - wc * q)|
      ≤ |α| * |va * (p : ℤ) - wa * q| + |β| * |vb * (p : ℤ) - wb * q|
        + |γ| * |vc * (p : ℤ) - wc * q| := by
    calc |α * (va * (p : ℤ) - wa * q) + β * (vb * (p : ℤ) - wb * q)
        + γ * (vc * (p : ℤ) - wc * q)|
        ≤ |α * (va * (p : ℤ) - wa * q) + β * (vb * (p : ℤ) - wb * q)|
          + |γ * (vc * (p : ℤ) - wc * q)| := habs_add _ _
      _ ≤ |α * (va * (p : ℤ) - wa * q)| + |β * (vb * (p : ℤ) - wb * q)|
          + |γ * (vc * (p : ℤ) - wc * q)| := by
          have := habs_add (α * (va * (p : ℤ) - wa * q))
            (β * (vb * (p : ℤ) - wb * q))
          linarith
      _ = |α| * |va * (p : ℤ) - wa * q| + |β| * |vb * (p : ℤ) - wb * q|
          + |γ| * |vc * (p : ℤ) - wc * q| := by
          rw [abs_mul, abs_mul, abs_mul]
  have hga : 14 * |va * (p : ℤ) - wa * q| ≤ (q : ℤ) - 1 := by omega
  have hgb : 14 * |vb * (p : ℤ) - wb * q| ≤ (q : ℤ) - 1 := by omega
  have hgc : 14 * |vc * (p : ℤ) - wc * q| ≤ (q : ℤ) - 1 := by omega
  by_contra hne
  have h1 : 1 ≤ |α * wa + β * wb + γ * wc| := Int.one_le_abs hne
  have h2 : (q : ℤ) ≤ |(α * wa + β * wb + γ * wc) * (q : ℤ)| := by
    rw [abs_mul, abs_of_pos hqZ]
    nlinarith [h1]
  rw [hkey, abs_neg] at h2
  have hq1 : (0 : ℤ) ≤ (q : ℤ) - 1 := by omega
  nlinarith [htri, h2, hga, hgb, hgc, hwt, hq1,
    abs_nonneg α, abs_nonneg β, abs_nonneg γ,
    abs_nonneg (va * (p : ℤ) - wa * q), abs_nonneg (vb * (p : ℤ) - wb * q),
    abs_nonneg (vc * (p : ℤ) - wc * q)]

/-- **The sum-triple lock**: `{a, b, a+b}` has coefficient weight 3 — the
witnesses ALWAYS satisfy `w_c = w_a + w_b`, even when both extreme pairs are
sparse (e.g. `(5, 6, 11)`). -/
theorem sum_triple_lock (va vb wa wb wc : ℤ) (q p : ℕ)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q)
    (hc : 14 * |(va + vb) * (p : ℤ) - wc * q| < q) :
    wc = wa + wb := by
  have h := relation_lock3 va vb (va + vb) wa wb wc 1 1 (-1) q p
    (by ring) (by norm_num) ha hb hc
  linarith [h]

/-- **The corrected pair-lock boundary**: the ratio relation has weight
`i'+j'`, so pairs lock through weight 14 — extending THM-967's ≤ 13. -/
theorem rational_lock_weight14 (va vb wa wb : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 14)
    (hrel : (j' : ℤ) * va = (i' : ℤ) * vb)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q) :
    (j' : ℤ) * wa = (i' : ℤ) * wb := by
  rcases Nat.eq_zero_or_pos q with rfl | hq
  · exfalso
    push_cast at ha
    linarith [abs_nonneg (va * (p : ℤ) - wa * 0)]
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hzero : 14 * |(0 : ℤ) * (p : ℤ) - 0 * q| < q := by
    simp only [zero_mul, sub_zero, abs_zero, mul_zero]
    exact hqZ
  have h := relation_lock3 va vb 0 wa wb 0 (j' : ℤ) (-(i' : ℤ)) 0 q p
    (by linarith [hrel])
    (by
      rw [abs_of_nonneg (by positivity : (0 : ℤ) ≤ (j' : ℤ)), abs_neg,
        abs_of_nonneg (by positivity : (0 : ℤ) ≤ (i' : ℤ)), abs_zero]
      have hZ : (j' : ℤ) + (i' : ℤ) ≤ 14 := by exact_mod_cast (by omega : j' + i' ≤ 14)
      linarith)
    ha hb hzero
  linarith [h]

/-- **The mediant triple, narrow form**: the projective chain
`(g·i', g·j', g·(i'+j'))` with `gcd(i',j') = 1`, `i'+j' ≤ 13`: all three fail
⟺ the gcd-speed fails the `14·(i'+j')`-narrow band. -/
theorem mediant_triple_fail_iff (g : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 13)
    (hcop : Nat.Coprime i' j') :
    ((∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
     (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q) ∧
     (∃ wc : ℤ, 14 * |(g * (i' + j' : ℕ)) * (p : ℤ) - wc * q| < q))
      ↔ (∃ s : ℤ, 14 * ((i' + j' : ℕ) : ℤ) * |g * (p : ℤ) - s * q| < q) := by
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (0 : ℤ) < (j' : ℤ) := by exact_mod_cast hj
  constructor
  · rintro ⟨⟨wa, ha⟩, ⟨wb, hb⟩, ⟨wc, hc⟩⟩
    -- the pair (a, b) is locked (i'+j' ≤ 13), giving s with wa = i's, wb = j's
    obtain ⟨s, hs⟩ := (rational_pair_fail_iff g i' j' q p hi hj hsum hcop).mp
      ⟨⟨wa, ha⟩, ⟨wb, hb⟩⟩
    -- the sum-triple lock pins wc = wa + wb = (i'+j')·s
    have hlock : wc = wa + wb := by
      have hc' : 14 * |(g * i' + g * j') * (p : ℤ) - wc * q| < q := by
        have : (g * (i' + j' : ℕ) : ℤ) = g * i' + g * j' := by push_cast; ring
        rwa [this] at hc
      exact sum_triple_lock (g * i') (g * j') wa wb wc q p ha hb hc'
    -- witnesses on the locked ray: wa = i'·s (from the pair-iff proof shape)
    -- recover the narrow bound directly from the c-band
    have hwa : wa = (i' : ℤ) * s := by
      have hlockp : (j' : ℤ) * wa = (i' : ℤ) * wb :=
        rational_lock (g * i') (g * j') wa wb i' j' q p hi hj hsum
          (by ring) ha hb
      have hdvd : (i' : ℤ) ∣ (j' : ℤ) * wa := ⟨wb, hlockp⟩
      have hcopZ : IsCoprime (i' : ℤ) (j' : ℤ) := by
        rw [Int.isCoprime_iff_gcd_eq_one]
        simpa using hcop
      obtain ⟨t, ht⟩ := hcopZ.dvd_of_dvd_mul_left hdvd
      -- both s and t witness the narrow band; uniqueness pins t = s
      have hXa : (g * i') * (p : ℤ) - wa * q = (i' : ℤ) * (g * p - t * q) := by
        rw [ht]; ring
      rw [hXa, abs_mul, abs_of_pos hiZ] at ha
      have hnarrow_t : 14 * |g * (p : ℤ) - t * q| < q := by
        nlinarith [ha, abs_nonneg (g * (p : ℤ) - t * q), hiZ]
      have hnarrow_s : 14 * |g * (p : ℤ) - s * q| < q := by
        have hM : (0 : ℤ) < ((max i' j' : ℕ) : ℤ) := by
          exact_mod_cast lt_of_lt_of_le hi (le_max_left i' j')
        nlinarith [hs, abs_nonneg (g * (p : ℤ) - s * q), hM]
      have := witness_unique g t s q p (by
        rcases Nat.eq_zero_or_pos q with rfl | hq
        · exfalso
          have habs : (0 : ℤ) ≤ |g * (p : ℤ) - t * 0| := abs_nonneg _
          push_cast at hnarrow_t
          linarith
        · exact hq) hnarrow_t hnarrow_s
      rw [ht, this]
  -- from wc = wa + wb = i's + j's, the c-band gives the (i'+j') narrow bound
    have hwb : wb = (j' : ℤ) * s := by
      have hlockp : (j' : ℤ) * wa = (i' : ℤ) * wb :=
        rational_lock (g * i') (g * j') wa wb i' j' q p hi hj hsum
          (by ring) ha hb
      have h1 : (i' : ℤ) * wb = (i' : ℤ) * ((j' : ℤ) * s) := by
        calc (i' : ℤ) * wb = (j' : ℤ) * wa := hlockp.symm
          _ = (j' : ℤ) * ((i' : ℤ) * s) := by rw [hwa]
          _ = (i' : ℤ) * ((j' : ℤ) * s) := by ring
      exact mul_left_cancel₀ (by linarith) h1
    refine ⟨s, ?_⟩
    have hXc : (g * (i' + j' : ℕ)) * (p : ℤ) - wc * q
        = ((i' : ℤ) + j') * (g * p - s * q) := by
      rw [hlock, hwa, hwb]
      push_cast
      ring
    rw [hXc, abs_mul, abs_of_pos (by linarith : (0 : ℤ) < (i' : ℤ) + j')] at hc
    calc 14 * ((i' + j' : ℕ) : ℤ) * |g * (p : ℤ) - s * q|
        = 14 * (((i' : ℤ) + j') * |g * (p : ℤ) - s * q|) := by push_cast; ring
      _ < q := hc
  · rintro ⟨s, hs⟩
    have habs : (0 : ℤ) ≤ |g * (p : ℤ) - s * q| := abs_nonneg _
    have hsumZ : ((i' + j' : ℕ) : ℤ) = (i' : ℤ) + j' := by push_cast; ring
    have hib : (i' : ℤ) ≤ (i' : ℤ) + j' := by linarith
    have hjb : (j' : ℤ) ≤ (i' : ℤ) + j' := by linarith
    rw [hsumZ] at hs
    refine ⟨⟨(i' : ℤ) * s, ?_⟩, ⟨(j' : ℤ) * s, ?_⟩, ⟨((i' : ℤ) + j') * s, ?_⟩⟩
    · have hXa : (g * i') * (p : ℤ) - ((i' : ℤ) * s) * q
          = (i' : ℤ) * (g * p - s * q) := by ring
      rw [hXa, abs_mul, abs_of_pos hiZ]
      nlinarith [hs, habs]
    · have hXb : (g * j') * (p : ℤ) - ((j' : ℤ) * s) * q
          = (j' : ℤ) * (g * p - s * q) := by ring
      rw [hXb, abs_mul, abs_of_pos hjZ]
      nlinarith [hs, habs]
    · have hXc : (g * (i' + j' : ℕ)) * (p : ℤ) - (((i' : ℤ) + j') * s) * q
          = ((i' : ℤ) + j') * (g * p - s * q) := by push_cast; ring
      rw [hXc, abs_mul, abs_of_pos (by linarith : (0 : ℤ) < (i' : ℤ) + j')]
      nlinarith [hs, habs]

open Classical in
/-- **THE MEDIANT TRIPLE COUNT**: at coprime moduli the sum-triple
`(g·i', g·j', g·(i'+j'))` jointly fails on EXACTLY
`2·⌊(q−1)/(14·(i'+j'))⌋` multipliers — the triple layer's first exact rung. -/
theorem mediant_triple_count (g : ℤ) (i' j' : ℕ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 13)
    (hcop : Nat.Coprime i' j') (hgcd : Int.gcd g (q : ℤ) = 1) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
        (∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
        (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q) ∧
        (∃ wc : ℤ, 14 * |(g * (i' + j' : ℕ)) * (p : ℤ) - wc * q| < q)).card
      = 2 * ((q - 1) / (14 * (i' + j'))) := by
  have hM : 1 ≤ i' + j' := by omega
  have hfilter : ((Finset.Ioo 0 q).filter fun p : ℕ =>
      (∃ wa : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q) ∧
      (∃ wb : ℤ, 14 * |(g * j') * (p : ℤ) - wb * q| < q) ∧
      (∃ wc : ℤ, 14 * |(g * (i' + j' : ℕ)) * (p : ℤ) - wc * q| < q))
      = ((Finset.Ioo 0 q).filter fun p : ℕ =>
          14 * (i' + j') * ((g * (p : ℤ) % (q : ℤ)).toNat) < q ∨
          14 * (i' + j') * (q - ((g * (p : ℤ) % (q : ℤ)).toNat)) < q) := by
    apply Finset.filter_congr
    intro p _
    rw [mediant_triple_fail_iff g i' j' q p hi hj hsum hcop,
      witness_mod_bridge (g * (p : ℤ)) (i' + j') q hM hq]
  rw [hfilter]
  have htrans := card_mod_filter_eq g q hq hgcd
    (fun r => 14 * (i' + j') * r < q ∨ 14 * (i' + j') * (q - r) < q)
  rw [htrans]
  exact narrowFailN_count (i' + j') q hM hq

/-! ## Axiom audit -/
#print axioms relation_lock
#print axioms relation_lock3
#print axioms sum_triple_lock
#print axioms rational_lock_weight14
#print axioms mediant_triple_fail_iff
#print axioms mediant_triple_count

end LRC14Concrete
end LonelyRunner
