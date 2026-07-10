/-
  TournamentH7.LRCHyperbolaBox — LEM-022's combinatorial heart in Lean
  (death-star-2026-07-09-S10, HYP-5870).

  LEM-022 (death-star-S9, klein-S226's handoff (a)) bounds the interval multiplicative
  pair-correlation `|C_w − b²/q| ≤ 6q(1+log₂q)²/P(w)` with `P(w) = min_{h≠0} |h|·|wh|` the
  minimal ratio-lattice product.  Its combinatorial heart — the piece that does all the work —
  is THE SEPARATION COUNT:

  > solutions of `|k| ≤ K`, `|wk| ≤ M` are pairwise `P/(2M)`-separated on the circle, hence
  > `(N − 1)·P ≤ 4·K·M`.

  This file proves that count kernel-pure and DIVISION-FREE:
    · the circle-metric API on `ZMod q` (`cdist`, subadditivity, negation-invariance, the
      signed representative and its two bridge lemmas);
    · `card_mulsep_in_Icc` — the multiplicative-separation pigeonhole (division appears only
      inside the fiber map `t ↦ ((t−a)·D)/S`; all bounds assemble multiplicatively) — the
      third instance of the clamped-fiber pattern (`exists_free_piece`, `dvd_Ioc_card_le`);
    · `hyperbola_box_count` — the separation count itself, `(N − 1)·P ≤ 4·K·M`.

  The remaining LEM-022 Lean surface (honestly flagged): the Fourier completion
  `|C_w − b²/q| ≤ (q/4)·Σ_{k≠0} 1/(|k|·|wk|)` (finite characters over `ℤ/q` + the sine bound)
  and the dyadic assembly — bounded, but ℂ-infrastructure-heavy; this file's count is the
  input both consume.

  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib

namespace LonelyRunner
namespace HyperbolaBox

variable {q : ℕ} [NeZero q]

/-- The circle metric to `0` on `ZMod q`: `cdist z = min(val z, q − val z)`. -/
def cdist (z : ZMod q) : ℕ := min z.val (q - z.val)

theorem cdist_zero : cdist (0 : ZMod q) = 0 := by
  simp [cdist, ZMod.val_zero]

theorem cdist_le_val (z : ZMod q) : cdist z ≤ z.val := min_le_left _ _

/-- Negation-invariance of the circle metric. -/
theorem cdist_neg (z : ZMod q) : cdist (-z) = cdist z := by
  rcases eq_or_ne z 0 with rfl | hz
  · simp
  · have : NeZero z := ⟨hz⟩
    have hval : (-z).val = q - z.val := ZMod.val_neg_of_ne_zero z
    have hlt : z.val < q := ZMod.val_lt z
    have hpos : 0 < z.val := ZMod.val_pos.mpr hz
    unfold cdist
    rw [hval]
    omega

/-- **Subadditivity of the circle metric**: `cdist (x + y) ≤ cdist x + cdist y`. -/
theorem cdist_add_le (x y : ZMod q) : cdist (x + y) ≤ cdist x + cdist y := by
  have hq : 0 < q := NeZero.pos q
  have hx : x.val < q := ZMod.val_lt x
  have hy : y.val < q := ZMod.val_lt y
  have hval : (x + y).val = (x.val + y.val) % q := ZMod.val_add x y
  unfold cdist
  rw [hval]
  rcases Nat.lt_or_ge (x.val + y.val) q with hlt | hge
  · rw [Nat.mod_eq_of_lt hlt]
    omega
  · have h2q : x.val + y.val < 2 * q := by omega
    have : (x.val + y.val) % q = x.val + y.val - q := by
      rw [Nat.mod_eq_sub_mod hge, Nat.mod_eq_of_lt (by omega)]
    rw [this]
    omega

/-- The signed representative: the unique lift of `z` to `(−q/2, q/2]` (as an integer). -/
def sgnRep (z : ZMod q) : ℤ :=
  if 2 * z.val ≤ q then (z.val : ℤ) else (z.val : ℤ) - q

/-- `|sgnRep z| = cdist z`. -/
theorem natAbs_sgnRep (z : ZMod q) : (sgnRep z).natAbs = cdist z := by
  have hlt : z.val < q := ZMod.val_lt z
  unfold sgnRep cdist
  split_ifs with h
  · omega
  · omega

/-- The signed representative is a genuine representative: `(sgnRep z : ZMod q) = z`. -/
theorem sgnRep_cast (z : ZMod q) : ((sgnRep z : ℤ) : ZMod q) = z := by
  unfold sgnRep
  split_ifs with h
  · push_cast
    simp [ZMod.natCast_val, ZMod.cast_id]
  · push_cast
    simp [ZMod.natCast_val, ZMod.cast_id, ZMod.natCast_self]

/-- The signed representative is injective (it inverts through the cast). -/
theorem sgnRep_injective : Function.Injective (sgnRep (q := q)) := by
  intro x y hxy
  have hx := sgnRep_cast x
  have hy := sgnRep_cast y
  rw [← hx, ← hy, hxy]

/-- **The circle metric is bounded by any integer representative**:
if `(n : ZMod q) = z` then `cdist z ≤ n.natAbs`. -/
theorem cdist_le_natAbs {n : ℤ} {z : ZMod q} (h : (n : ZMod q) = z) :
    cdist z ≤ n.natAbs := by
  rcases Int.natAbs_eq n with hn | hn
  · -- n = natAbs n
    have h2 : ((n.natAbs : ℤ) : ZMod q) = z := by rw [← hn]; exact h
    have hcast : ((n.natAbs : ℕ) : ZMod q) = z := by
      rw [← Int.cast_natCast]; exact h2
    have hval : z.val = n.natAbs % q := by
      rw [← hcast, ZMod.val_natCast]
    calc cdist z ≤ z.val := cdist_le_val z
      _ = n.natAbs % q := hval
      _ ≤ n.natAbs := Nat.mod_le _ _
  · -- n = −natAbs n : use negation invariance
    have hneg : -n = (n.natAbs : ℤ) := by omega
    have h2 : ((n.natAbs : ℤ) : ZMod q) = -z := by
      rw [← hneg]
      push_cast
      rw [h]
    have hcast : ((n.natAbs : ℕ) : ZMod q) = -z := by
      rw [← Int.cast_natCast]; exact h2
    have hval : (-z).val = n.natAbs % q := by
      rw [← hcast, ZMod.val_natCast]
    calc cdist z = cdist (-z) := (cdist_neg z).symm
      _ ≤ (-z).val := cdist_le_val _
      _ = n.natAbs % q := hval
      _ ≤ n.natAbs := Nat.mod_le _ _

/-! ## The multiplicative-separation pigeonhole -/

/-- **The multiplicative-separation pigeonhole.**  A finite set of integers in `[a, b]` whose
pairwise differences satisfy `S ≤ |x − y|·D` has at most `(b−a)·D/S + 1` members — stated
division-free: `(|T| − 1)·S ≤ (b − a)·D`.  The fiber map `t ↦ ((t−a)·D)/S` is the only place
division appears; same-fiber points would violate the separation (their `ediv/emod`
decompositions force `|x−y|·D < S`). -/
theorem card_mulsep_in_Icc (T : Finset ℤ) (a b S D : ℤ) (hS : 0 < S) (hD : 0 < D)
    (hab : a ≤ b)
    (hT : ∀ t ∈ T, a ≤ t ∧ t ≤ b)
    (hsep : ∀ x ∈ T, ∀ y ∈ T, x ≠ y → S ≤ |x - y| * D) :
    ((T.card : ℤ) - 1) * S ≤ (b - a) * D := by
  classical
  set X : ℤ := (b - a) * D / S with hX
  have hX0 : 0 ≤ X := Int.ediv_nonneg (by positivity) (le_of_lt hS)
  -- the fiber injection
  have hcard : T.card ≤ (Finset.Icc (0 : ℤ) X).card := by
    apply Finset.card_le_card_of_injOn (fun t => (t - a) * D / S)
    · intro t ht
      obtain ⟨hta, htb⟩ := hT t ht
      refine Finset.mem_Icc.mpr ⟨Int.ediv_nonneg (by positivity) (le_of_lt hS), ?_⟩
      exact Int.ediv_le_ediv hS (by nlinarith)
    · intro x hx y hy hxy
      by_contra hne
      have hsepxy := hsep x hx y hy hne
      have hxy' : (x - a) * D / S = (y - a) * D / S := hxy
      -- decompose both fibers
      have hxd := Int.ediv_add_emod ((x - a) * D) S
      have hyd := Int.ediv_add_emod ((y - a) * D) S
      have hxr0 : 0 ≤ (x - a) * D % S := Int.emod_nonneg _ (ne_of_gt hS)
      have hxr1 : (x - a) * D % S < S := Int.emod_lt_of_pos _ hS
      have hyr0 : 0 ≤ (y - a) * D % S := Int.emod_nonneg _ (ne_of_gt hS)
      have hyr1 : (y - a) * D % S < S := Int.emod_lt_of_pos _ hS
      -- same fiber ⟹ |x − y|·D = |r_x − r_y| < S
      have hdiff : (x - y) * D = (x - a) * D % S - (y - a) * D % S := by
        have h1 : (x - a) * D = S * ((x - a) * D / S) + (x - a) * D % S := hxd.symm
        have h2 : (y - a) * D = S * ((y - a) * D / S) + (y - a) * D % S := hyd.symm
        have h3 : (x - a) * D - (y - a) * D = (x - y) * D := by ring
        rw [hxy'] at h1
        omega
      have habs : |x - y| * D = |(x - y) * D| := by
        rw [abs_mul, abs_of_pos hD]
      have : |(x - y) * D| < S := by
        rw [hdiff]
        rw [abs_lt]
        omega
      rw [← habs] at this
      omega
  have hXcard : ((Finset.Icc (0 : ℤ) X).card : ℤ) = X + 1 := by
    rw [Int.card_Icc]
    omega
  have hle : (T.card : ℤ) - 1 ≤ X := by
    have : (T.card : ℤ) ≤ ((Finset.Icc (0 : ℤ) X).card : ℤ) := by exact_mod_cast hcard
    omega
  calc ((T.card : ℤ) - 1) * S ≤ X * S := by
        apply mul_le_mul_of_nonneg_right hle (le_of_lt hS)
    _ ≤ (b - a) * D := Int.ediv_mul_le _ (ne_of_gt hS)

/-! ## The hyperbola box count (LEM-022 Step 2) -/

/-- **The separation count / hyperbola box count (LEM-022's combinatorial heart).**
If every nonzero `h` has ratio-lattice product `≥ P` (`P ≤ cdist h · cdist (w·h)`), then the
box `{k ≠ 0 : cdist k ≤ K, cdist (w·k) ≤ M}` holds at most `1 + 4KM/P` points — division-free:
`(N − 1)·P ≤ 4·K·M`.  Mechanism: signed representatives live in `[−K, K]` and are pairwise
multiplicatively separated (`P ≤ |Δ|·2M`, via subadditivity of the circle metric applied to
the difference), so the pigeonhole above applies with `D = 2M`. -/
theorem hyperbola_box_count (w : ZMod q) (K M P : ℕ)
    (hPmin : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist (w * h)) :
    ((((Finset.univ : Finset (ZMod q)).filter
        fun k => k ≠ 0 ∧ cdist k ≤ K ∧ cdist (w * k) ≤ M).card : ℤ) - 1) * (P : ℤ)
      ≤ 4 * K * M := by
  classical
  set Sf := (Finset.univ : Finset (ZMod q)).filter
      fun k => k ≠ 0 ∧ cdist k ≤ K ∧ cdist (w * k) ≤ M with hSf
  set T : Finset ℤ := Sf.image sgnRep with hT
  rcases Nat.eq_zero_or_pos P with rfl | hP
  · -- P = 0 : the bound is trivial
    have hz : ((Sf.card : ℤ) - 1) * ((0 : ℕ) : ℤ) = 0 := by
      simp
    rw [hz]
    exact mul_nonneg (mul_nonneg (by norm_num) (Int.natCast_nonneg K)) (Int.natCast_nonneg M)
  rcases Nat.eq_zero_or_pos M with rfl | hM
  · -- M = 0 : the box is empty (P > 0 forbids cdist (w·k) = 0 for k ≠ 0)
    have hempty : Sf = ∅ := by
      rw [hSf, Finset.filter_eq_empty_iff]
      rintro k - ⟨hk0, -, hkM⟩
      have hfloor := hPmin k hk0
      have hz : cdist (w * k) = 0 := Nat.le_zero.mp hkM
      rw [hz, Nat.mul_zero] at hfloor
      omega
    rw [hempty]
    simp only [Finset.card_empty, Nat.cast_zero, Nat.cast_ofNat, mul_zero]
    have hPnn : (0 : ℤ) ≤ (P : ℤ) := Int.natCast_nonneg P
    nlinarith
  have hTcard : T.card = Sf.card :=
    Finset.card_image_of_injective _ sgnRep_injective
  have hbound := card_mulsep_in_Icc T (-(K : ℤ)) (K : ℤ) (P : ℤ) (2 * M)
    (by exact_mod_cast hP)
    (by positivity)
    (by omega)
    ?_ ?_
  · calc ((Sf.card : ℤ) - 1) * (P : ℤ) = ((T.card : ℤ) - 1) * (P : ℤ) := by rw [hTcard]
      _ ≤ ((K : ℤ) - (-(K : ℤ))) * (2 * M) := hbound
      _ = 4 * K * M := by ring
  · -- membership: signed reps live in [−K, K]
    intro t ht
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.mp ht
    obtain ⟨_, hkK, _⟩ := (Finset.mem_filter.mp hk).2
    have habs : (sgnRep k).natAbs ≤ K := by
      rw [natAbs_sgnRep]; exact hkK
    constructor
    · omega
    · omega
  · -- separation: P ≤ |Δ|·2M
    intro x hx y hy hxy
    obtain ⟨k₁, hk₁, rfl⟩ := Finset.mem_image.mp hx
    obtain ⟨k₂, hk₂, rfl⟩ := Finset.mem_image.mp hy
    have hk12 : k₁ ≠ k₂ := fun h => hxy (by rw [h])
    obtain ⟨hk₁0, _, hk₁M⟩ := (Finset.mem_filter.mp hk₁).2
    obtain ⟨hk₂0, _, hk₂M⟩ := (Finset.mem_filter.mp hk₂).2
    set δ : ZMod q := k₁ - k₂ with hδ
    have hδ0 : δ ≠ 0 := sub_ne_zero_of_ne hk12
    -- the ratio-lattice floor at δ
    have hPδ := hPmin δ hδ0
    -- cdist (w δ) ≤ 2M by subadditivity
    have hwδ : cdist (w * δ) ≤ 2 * M := by
      have halg : w * δ = w * k₁ + (-(w * k₂)) := by rw [hδ]; ring
      calc cdist (w * δ) = cdist (w * k₁ + (-(w * k₂))) := by rw [halg]
        _ ≤ cdist (w * k₁) + cdist (-(w * k₂)) := cdist_add_le _ _
        _ = cdist (w * k₁) + cdist (w * k₂) := by rw [cdist_neg]
        _ ≤ M + M := Nat.add_le_add hk₁M hk₂M
        _ = 2 * M := by ring
    -- cdist δ ≤ |sgnRep k₁ − sgnRep k₂|
    have hrep : ((sgnRep k₁ - sgnRep k₂ : ℤ) : ZMod q) = δ := by
      push_cast
      rw [sgnRep_cast, sgnRep_cast, hδ]
    have hδle : cdist δ ≤ (sgnRep k₁ - sgnRep k₂).natAbs := cdist_le_natAbs hrep
    -- assemble: P ≤ cdist δ · cdist(wδ) ≤ |Δ| · 2M
    have hchain : P ≤ (sgnRep k₁ - sgnRep k₂).natAbs * (2 * M) := by
      calc P ≤ cdist δ * cdist (w * δ) := hPδ
        _ ≤ (sgnRep k₁ - sgnRep k₂).natAbs * (2 * M) :=
            Nat.mul_le_mul hδle hwδ
    have : (P : ℤ) ≤ ((sgnRep k₁ - sgnRep k₂).natAbs : ℤ) * (2 * M) := by
      exact_mod_cast hchain
    rwa [Int.natCast_natAbs] at this

/-! ## The dyadic assembly (LEM-022 Step 3, ℂ-free)

The harmonic ratio sum `S = Σ_{k≠0} 1/(cdist k · cdist (wk))` is what the Fourier completion
multiplies by `q/4`.  The per-fiber dichotomy makes the assembly SIMPLER than the paper proof:
an empty dyadic fiber contributes nothing, and a nonempty one contains a witness forcing
`P ≤ cdist·cdist < 2^{i+j+2}`, so its `2^{−(i+j)}` head is `< 4/P` — no geometric series
needed.  Final form, division-free: `S · P ≤ 20 · L²` with `L = log₂ q + 1`. -/

/-- Nonzero elements have `cdist ≥ 1`. -/
theorem one_le_cdist {z : ZMod q} (hz : z ≠ 0) : 1 ≤ cdist z := by
  have hval : z.val < q := ZMod.val_lt z
  have hpos : 0 < z.val := ZMod.val_pos.mpr hz
  unfold cdist
  omega

/-- **The dyadic assembly of the harmonic ratio sum** (LEM-022 Step 3):
`(Σ_{k≠0} 1/(cdist k · cdist(wk))) · P ≤ 20·(log₂ q + 1)²`, given the ratio-lattice floor
`P ≤ cdist h · cdist (wh)` for all `h ≠ 0`.  Together with the Fourier completion
(`error ≤ (q/4)·S`, the remaining ℂ step) this yields LEM-022's
`|C_w − b²/q| ≤ 5·q·(log₂q + 1)²/P`. -/
theorem harmonic_ratio_sum_mul_le (w : ZMod q) (P : ℕ)
    (hPmin : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist (w * h)) :
    (∑ k ∈ (Finset.univ : Finset (ZMod q)).filter (fun k => k ≠ 0),
        (1 : ℚ) / (cdist k * cdist (w * k))) * (P : ℚ)
      ≤ 20 * (Nat.log 2 q + 1) ^ 2 := by
  classical
  set L : ℕ := Nat.log 2 q + 1 with hL
  set S0 := (Finset.univ : Finset (ZMod q)).filter (fun k => k ≠ 0) with hS0
  set f : ZMod q → ℚ := fun k => (1 : ℚ) / (cdist k * cdist (w * k)) with hf
  set idx : ZMod q → ℕ × ℕ :=
    fun k => (Nat.log 2 (cdist k), Nat.log 2 (cdist (w * k))) with hidx
  -- every nonzero k lands in range L ×ˢ range L
  have hmaps : ∀ k ∈ S0, idx k ∈ Finset.range L ×ˢ Finset.range L := by
    intro k _
    have h1 : cdist k ≤ q := by
      unfold cdist; omega
    have h2 : cdist (w * k) ≤ q := by
      unfold cdist; omega
    have hlog1 : Nat.log 2 (cdist k) < L := by
      rw [hL]; exact Nat.lt_succ_of_le (Nat.log_mono_right h1)
    have hlog2 : Nat.log 2 (cdist (w * k)) < L := by
      rw [hL]; exact Nat.lt_succ_of_le (Nat.log_mono_right h2)
    simp only [hidx]
    exact Finset.mem_product.mpr ⟨Finset.mem_range.mpr hlog1, Finset.mem_range.mpr hlog2⟩
  -- fiberwise decomposition
  have hsplit : (∑ k ∈ S0, f k)
      = ∑ p ∈ Finset.range L ×ˢ Finset.range L, ∑ k ∈ S0.filter (fun k => idx k = p), f k :=
    (Finset.sum_fiberwise_of_maps_to hmaps f).symm
  -- per-fiber bound: fiber sum · P ≤ 20
  have hfiber : ∀ p ∈ Finset.range L ×ˢ Finset.range L,
      (∑ k ∈ S0.filter (fun k => idx k = p), f k) * (P : ℚ) ≤ 20 := by
    intro p _
    obtain ⟨i, j⟩ := p
    set F := S0.filter (fun k => idx k = (i, j)) with hF
    rcases Finset.eq_empty_or_nonempty F with hemp | hne
    · rw [hemp]
      simp
    · -- the witness forces P < 2^{i+j+2}
      obtain ⟨k₀, hk₀⟩ := hne
      have hk₀0 : k₀ ≠ 0 := (Finset.mem_filter.mp ((Finset.mem_filter.mp hk₀).1)).2
      have hk₀idx : idx k₀ = (i, j) := (Finset.mem_filter.mp hk₀).2
      have hi₀ : Nat.log 2 (cdist k₀) = i := by
        have := congrArg Prod.fst hk₀idx; simpa [hidx] using this
      have hj₀ : Nat.log 2 (cdist (w * k₀)) = j := by
        have := congrArg Prod.snd hk₀idx; simpa [hidx] using this
      have hc₀1 : 1 ≤ cdist k₀ := one_le_cdist hk₀0
      have hcw₀pos : P ≤ cdist k₀ * cdist (w * k₀) := hPmin k₀ hk₀0
      have hup1 : cdist k₀ < 2 ^ (i + 1) := by
        rw [← hi₀]
        exact Nat.lt_pow_succ_log_self (by norm_num) _
      have hup2 : cdist (w * k₀) < 2 ^ (j + 1) := by
        rw [← hj₀]
        exact Nat.lt_pow_succ_log_self (by norm_num) _
      have hPsmall : P < 2 ^ (i + j + 2) := by
        calc P ≤ cdist k₀ * cdist (w * k₀) := hcw₀pos
          _ < 2 ^ (i + 1) * 2 ^ (j + 1) :=
              Nat.mul_lt_mul_of_lt_of_le hup1 (le_of_lt hup2) (by positivity)
          _ = 2 ^ (i + j + 2) := by ring
      -- fiber cardinality via the box count at K = 2^{i+1}, M = 2^{j+1}
      have hsub : F ⊆ (Finset.univ : Finset (ZMod q)).filter
          (fun k => k ≠ 0 ∧ cdist k ≤ 2 ^ (i + 1) ∧ cdist (w * k) ≤ 2 ^ (j + 1)) := by
        intro k hk
        have hk0 : k ≠ 0 := (Finset.mem_filter.mp ((Finset.mem_filter.mp hk).1)).2
        have hkidx : idx k = (i, j) := (Finset.mem_filter.mp hk).2
        have hi' : Nat.log 2 (cdist k) = i := by
          have := congrArg Prod.fst hkidx; simpa [hidx] using this
        have hj' : Nat.log 2 (cdist (w * k)) = j := by
          have := congrArg Prod.snd hkidx; simpa [hidx] using this
        refine Finset.mem_filter.mpr ⟨Finset.mem_univ _, hk0, ?_, ?_⟩
        · rw [← hi']
          exact le_of_lt (Nat.lt_pow_succ_log_self (by norm_num) _)
        · rw [← hj']
          exact le_of_lt (Nat.lt_pow_succ_log_self (by norm_num) _)
      have hboxZ := hyperbola_box_count w (2 ^ (i + 1)) (2 ^ (j + 1)) P hPmin
      have hcardle : F.card ≤ ((Finset.univ : Finset (ZMod q)).filter
          (fun k => k ≠ 0 ∧ cdist k ≤ 2 ^ (i + 1) ∧ cdist (w * k) ≤ 2 ^ (j + 1))).card :=
        Finset.card_le_card hsub
      -- per-term bound: f k ≤ 1/2^{i+j} on the fiber
      have hterm : ∀ k ∈ F, f k ≤ (1 : ℚ) / 2 ^ (i + j) := by
        intro k hk
        have hk0 : k ≠ 0 := (Finset.mem_filter.mp ((Finset.mem_filter.mp hk).1)).2
        have hkidx : idx k = (i, j) := (Finset.mem_filter.mp hk).2
        have hi' : Nat.log 2 (cdist k) = i := by
          have := congrArg Prod.fst hkidx; simpa [hidx] using this
        have hj' : Nat.log 2 (cdist (w * k)) = j := by
          have := congrArg Prod.snd hkidx; simpa [hidx] using this
        rcases Nat.eq_zero_or_pos (cdist (w * k)) with hz | hwpos
        · -- zero denominator: f k = 1/0 = 0 in ℚ
          rw [hf]
          push_cast
          rw [hz]
          norm_num
        · have hlow1 : 2 ^ i ≤ cdist k := by
            rw [← hi']
            exact Nat.pow_log_le_self 2 (by have := one_le_cdist hk0; omega)
          have hlow2 : 2 ^ j ≤ cdist (w * k) := by
            rw [← hj']
            exact Nat.pow_log_le_self 2 (by omega)
          have hprodN : 2 ^ (i + j) ≤ cdist k * cdist (w * k) := by
            rw [pow_add]
            exact Nat.mul_le_mul hlow1 hlow2
          rw [hf]
          apply one_div_le_one_div_of_le (by positivity)
          exact_mod_cast hprodN
      calc (∑ k ∈ F, f k) * (P : ℚ)
          ≤ (F.card • ((1 : ℚ) / 2 ^ (i + j))) * (P : ℚ) := by
            apply mul_le_mul_of_nonneg_right _ (by positivity)
            exact Finset.sum_le_card_nsmul F f _ hterm
        _ = (F.card : ℚ) * (P : ℚ) / 2 ^ (i + j) := by
            rw [nsmul_eq_mul]; ring
        _ ≤ 20 := by
            -- (card)·P ≤ (1 + 16·2^{i+j}/P-form)·... : use the box count + the witness
            have hcardZ : ((F.card : ℤ) - 1) * (P : ℤ) ≤ 4 * 2 ^ (i + 1) * 2 ^ (j + 1) := by
              calc ((F.card : ℤ) - 1) * (P : ℤ)
                  ≤ ((((Finset.univ : Finset (ZMod q)).filter
                      (fun k => k ≠ 0 ∧ cdist k ≤ 2 ^ (i + 1) ∧ cdist (w * k) ≤ 2 ^ (j + 1))).card : ℤ)
                      - 1) * (P : ℤ) := by
                    apply mul_le_mul_of_nonneg_right _ (Int.natCast_nonneg P)
                    have := hcardle
                    omega
                _ ≤ 4 * 2 ^ (i + 1) * 2 ^ (j + 1) := by
                    have := hboxZ
                    push_cast at this ⊢
                    linarith
            -- convert: card·P ≤ P + 16·2^{i+j} < 2^{i+j+2} + 16·2^{i+j} = 20·2^{i+j}
            have h4 : (4 : ℤ) * 2 ^ (i + 1) * 2 ^ (j + 1) = 16 * 2 ^ (i + j) := by
              rw [pow_succ, pow_succ, pow_add]
              ring
            have hcardP : (F.card : ℤ) * (P : ℤ) ≤ (P : ℤ) + 16 * 2 ^ (i + j) := by
              have := hcardZ
              rw [h4] at this
              nlinarith [Int.natCast_nonneg P]
            have hPZ : (P : ℤ) < 2 ^ (i + j + 2) := by exact_mod_cast hPsmall
            have hfinal : (F.card : ℤ) * (P : ℤ) ≤ 20 * 2 ^ (i + j) := by
              have h2 : (2 : ℤ) ^ (i + j + 2) = 4 * 2 ^ (i + j) := by
                rw [pow_add]; ring
              omega
            have hQ : (F.card : ℚ) * (P : ℚ) ≤ 20 * 2 ^ (i + j) := by
              exact_mod_cast hfinal
            rw [div_le_iff₀ (by positivity : (0:ℚ) < 2 ^ (i + j))]
            linarith
  -- assemble
  rw [hsplit, Finset.sum_mul]
  calc ∑ p ∈ Finset.range L ×ˢ Finset.range L,
        (∑ k ∈ S0.filter (fun k => idx k = p), f k) * (P : ℚ)
      ≤ ∑ _p ∈ Finset.range L ×ˢ Finset.range L, (20 : ℚ) :=
        Finset.sum_le_sum hfiber
    _ = ((Finset.range L ×ˢ Finset.range L).card : ℚ) * 20 := by
        rw [Finset.sum_const, nsmul_eq_mul]
    _ = 20 * (L : ℚ) ^ 2 := by
        rw [Finset.card_product, Finset.card_range]
        push_cast
        ring
    _ = 20 * (Nat.log 2 q + 1) ^ 2 := by
        rw [hL]
        push_cast
        ring

/-! ## Axiom audit -/
#print axioms cdist_add_le
#print axioms cdist_le_natAbs
#print axioms card_mulsep_in_Icc
#print axioms hyperbola_box_count
#print axioms harmonic_ratio_sum_mul_le

end HyperbolaBox
end LonelyRunner
