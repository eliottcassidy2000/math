/-
  TournamentH7.LRCClusterGcd — THE CLUSTER-GCD LADDER, |S| ≤ 3 RUNGS
  (kind-pasteur-2026-07-05-S18, HYP-4227; executing the S17 draft handoff,
  drafts/cluster-gcd-ladder-kps-S17.md).

  THE THEOREM (this file, lossy-count form): a 12-family with no `2/25`-margin
  point satisfies, for every `S` with `1 ≤ |S| ≤ 3` and every common divisor
  `d` of the complement,

      (25 − 8·|S|) · d  ≤  25 · (Σ_{i∈S} |vᵢ| + |S|).

  Rungs: |S| = 1: `17d ≤ 25(G+1)` — THE HEADLINE: the gcd of ANY 11-subfamily
  is bounded by ~1.5× the remaining runner; |S| = 2: `9d ≤ 25(Σ+2)`;
  |S| = 3: `d ≤ 25(Σ+3)`.

  MECHANISM (the absolute-height engine the residue filters cannot supply —
  mac-mini-S55's ray-periodicity obstruction): the citation gives the big side
  `T = W∖S` (≤ 11 runners) margin `1/12 > 2/25` at some `t₀`; that margin is
  `1/d`-PERIODIC (every `T`-runner is a multiple of `d`), so all `d` copies
  `t₀ + j/d` are `T`-good; each must sit in an `S`-tooth; and one runner `w`
  can host at most `(8/25)d + w + 1` of the `d` equally-spaced copies
  (`≤ w+1` teeth in the window, `≤ (4/25)d/w + 1` copies per tooth).
  Pigeonhole closes.

  The SHARP count `(4/25)d + 2w` (pole |S| ≤ 6) is proved in the draft via
  exact equidistribution; upgrading this file's count lemma to it extends the
  rungs to |S| ≤ 6 with constant 50 — flagged as the follow-up.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCTeethR
import TournamentH7.LRCTripleWalk

namespace LonelyRunner
namespace ClusterGcd

open TeethR TripleWalk
open scoped Classical

/-- Margin transport along the `1/d`-period: a runner divisible by `d` keeps its
whole `∀ m` margin at every shifted point `t + j/d`. -/
lemma periodic_margin (v d : ℤ) (hd : 0 < d) (hdvd : d ∣ v) (t : ℝ) (j : ℤ)
    (β : ℝ) (h : ∀ m : ℤ, β ≤ |(v : ℝ) * t - m|) :
    ∀ m : ℤ, β ≤ |(v : ℝ) * (t + (j : ℝ) / d) - m| := by
  intro m
  obtain ⟨c, rfl⟩ := hdvd
  have hdR : ((d : ℝ)) ≠ 0 := by
    have : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd
    exact ne_of_gt this
  have heq : ((d * c : ℤ) : ℝ) * (t + (j : ℝ) / d)
      = ((d * c : ℤ) : ℝ) * t + ((c * j : ℤ) : ℝ) := by
    push_cast
    field_simp
  rw [heq, show ((d * c : ℤ) : ℝ) * t + ((c * j : ℤ) : ℝ) - (m : ℝ)
      = ((d * c : ℤ) : ℝ) * t - ((m - c * j : ℤ) : ℝ) by push_cast; ring]
  exact h (m - c * j)

/-- Integers of `Finset.range D` in an open real interval `(A, B)` (with `A ≤ B`)
number at most `B − A + 1`. -/
lemma range_filter_interval_card (D : ℕ) (A B : ℝ) (hAB : A ≤ B)
    (P : ℕ → Prop) [DecidablePred P]
    (hP : ∀ j ∈ Finset.range D, P j → A < (j : ℝ) ∧ (j : ℝ) < B) :
    (((Finset.range D).filter P).card : ℝ) ≤ B - A + 1 := by
  have hsub : ((Finset.range D).filter P).image (fun j : ℕ => (j : ℤ))
      ⊆ Finset.Icc ⌈A⌉ ⌊B⌋ := by
    intro z hz
    rw [Finset.mem_image] at hz
    obtain ⟨j, hj, rfl⟩ := hz
    rw [Finset.mem_filter] at hj
    obtain ⟨hjr, hjP⟩ := hj
    obtain ⟨h1, h2⟩ := hP j hjr hjP
    rw [Finset.mem_Icc]
    refine ⟨Int.ceil_le.mpr ?_, Int.le_floor.mpr ?_⟩
    · exact_mod_cast le_of_lt h1
    · exact_mod_cast le_of_lt h2
  have hinj : (((Finset.range D).filter P).image (fun j : ℕ => (j : ℤ))).card
      = ((Finset.range D).filter P).card := by
    apply Finset.card_image_of_injective
    exact fun a b hab => by exact_mod_cast hab
  have hcard : ((Finset.range D).filter P).card ≤ (Finset.Icc ⌈A⌉ ⌊B⌋).card := by
    rw [← hinj]
    exact Finset.card_le_card hsub
  have hIcc : ((Finset.Icc ⌈A⌉ ⌊B⌋).card : ℝ) ≤ B - A + 1 := by
    rw [Int.card_Icc]
    rcases le_or_gt (⌊B⌋ + 1 - ⌈A⌉) 0 with hneg | hpos
    · rw [Int.toNat_of_nonpos hneg]
      simp only [Nat.cast_zero]
      linarith
    · have hnat : (((⌊B⌋ + 1 - ⌈A⌉).toNat : ℕ) : ℝ) = ((⌊B⌋ + 1 - ⌈A⌉ : ℤ) : ℝ) := by
        exact_mod_cast Int.toNat_of_nonneg (le_of_lt hpos)
      rw [hnat]
      push_cast
      have hf : (⌊B⌋ : ℝ) ≤ B := Int.floor_le B
      have hc : A ≤ (⌈A⌉ : ℝ) := Int.le_ceil A
      linarith
  calc (((Finset.range D).filter P).card : ℝ)
      ≤ ((Finset.Icc ⌈A⌉ ⌊B⌋).card : ℝ) := by exact_mod_cast hcard
    _ ≤ B - A + 1 := hIcc

/-- **The tooth-visit count (lossy form)**: among `d` equally-spaced points
`c + j/d`, `j < d`, at most `(8/25)d + w + 1` are within `2/25` (scaled by `w`)
of an integer. -/
lemma tooth_visit_count (w : ℤ) (hw : 0 < w) (D : ℕ) (hD : 0 < D) (c : ℝ) :
    (((Finset.range D).filter
        (fun j : ℕ => ∃ m : ℤ, |(w : ℝ) * (c + (j : ℝ) / D) - m| < 2/25)).card : ℝ)
      ≤ (8/25) * D + w + 1 := by
  classical
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hw1 : (1 : ℝ) ≤ (w : ℝ) := by exact_mod_cast hw
  have hDR : (0 : ℝ) < (D : ℝ) := by exact_mod_cast hD
  set V : Finset ℕ := (Finset.range D).filter
    (fun j : ℕ => ∃ m : ℤ, |(w : ℝ) * (c + (j : ℝ) / D) - m| < 2/25) with hV
  -- classify by the (unique) witnessing integer = the round
  set F : ℕ → ℤ := fun j => round ((w : ℝ) * (c + (j : ℝ) / D)) with hF
  set M : Finset ℤ := Finset.Icc ⌈(w : ℝ) * c - 2/25⌉ ⌊(w : ℝ) * c + w + 2/25⌋ with hM
  have hFmem : ∀ j ∈ V, F j ∈ M := by
    intro j hj
    rw [hV, Finset.mem_filter] at hj
    obtain ⟨hjr, m, hm⟩ := hj
    have hround : F j = m := by
      rw [hF]
      exact GridAttainment.round_eq_of_abs_sub_lt_half _ m (by linarith [hm] )
    have hj0 : (0 : ℝ) ≤ (j : ℝ) / D := by positivity
    have hjlt : (j : ℝ) / D < 1 := by
      rw [div_lt_one hDR]
      exact_mod_cast Finset.mem_range.mp hjr
    have hx1 : (w : ℝ) * c ≤ (w : ℝ) * (c + (j : ℝ) / D) := by nlinarith
    have hx2 : (w : ℝ) * (c + (j : ℝ) / D) < (w : ℝ) * c + w := by nlinarith
    have habs := abs_lt.mp hm
    rw [hround, hM, Finset.mem_Icc]
    constructor
    · apply Int.ceil_le.mpr
      push_cast
      linarith
    · apply Int.le_floor.mpr
      push_cast
      linarith
  -- fiberwise count
  have hsplit : V.card = ∑ m ∈ M, (V.filter (fun j => F j = m)).card :=
    Finset.card_eq_sum_card_fiberwise hFmem
  -- per-fiber: the j's live in an open interval of length (4/25)·D/w
  have hfiber : ∀ m ∈ M, ((V.filter (fun j => F j = m)).card : ℝ)
      ≤ (4/25) * D / w + 1 := by
    intro m _
    have hsubset : V.filter (fun j => F j = m)
        ⊆ (Finset.range D).filter (fun j : ℕ =>
            ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w < (j : ℝ) ∧
            (j : ℝ) < ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w) := by
      intro j hj
      rw [Finset.mem_filter] at hj ⊢
      obtain ⟨hjV, hjF⟩ := hj
      rw [hV, Finset.mem_filter] at hjV
      obtain ⟨hjr, m', hm'⟩ := hjV
      have hround : F j = m' := by
        rw [hF]
        exact GridAttainment.round_eq_of_abs_sub_lt_half _ m' (by linarith)
      have hmm : m' = m := by rw [← hround, hjF]
      subst hmm
      have hexp : (w : ℝ) * (c + (j : ℝ) / D) = (w : ℝ) * c + (w : ℝ) * (j : ℝ) / D := by
        ring
      rw [hexp] at hm'
      have habs := abs_lt.mp hm'
      refine ⟨hjr, ?_, ?_⟩
      · rw [div_lt_iff₀ hwR]
        have h1 : ((m' : ℝ) - (w : ℝ) * c - 2/25) < (w : ℝ) * (j : ℝ) / D := by
          linarith [habs.1]
        have h2 := (lt_div_iff₀ hDR).mp h1
        linarith [h2]
      · rw [lt_div_iff₀ hwR]
        have h1 : (w : ℝ) * (j : ℝ) / D < ((m' : ℝ) - (w : ℝ) * c + 2/25) := by
          linarith [habs.2]
        have h2 := (div_lt_iff₀ hDR).mp h1
        linarith [h2]
    have hABm : ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w
        ≤ ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w := by
      apply div_le_div_of_nonneg_right ?_ hwR.le
      nlinarith [hDR.le]
    have hbound := range_filter_interval_card D
      (((m : ℝ) - (w : ℝ) * c - 2/25) * D / w)
      (((m : ℝ) - (w : ℝ) * c + 2/25) * D / w) hABm
      (fun j : ℕ => ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w < (j : ℝ) ∧
            (j : ℝ) < ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w)
      (fun j hj hP => hP)
    have hlen : ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w
        - ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w = (4/25) * D / w := by
      field_simp
      ring
    calc ((V.filter (fun j => F j = m)).card : ℝ)
        ≤ (((Finset.range D).filter (fun j : ℕ =>
            ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w < (j : ℝ) ∧
            (j : ℝ) < ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w)).card : ℝ) := by
          exact_mod_cast Finset.card_le_card hsubset
      _ ≤ ((m : ℝ) - (w : ℝ) * c + 2/25) * D / w
            - ((m : ℝ) - (w : ℝ) * c - 2/25) * D / w + 1 := hbound
      _ = (4/25) * D / w + 1 := by rw [hlen]
  -- the m-range is small: |M| ≤ w + 1
  have hMZ : (M.card : ℤ) ≤ w + 1 := by
    have hfc : (⌊(w : ℝ) * c + w + 2/25⌋ : ℝ) - (⌈(w : ℝ) * c - 2/25⌉ : ℝ) < (w : ℝ) + 1 := by
      have hf : (⌊(w : ℝ) * c + w + 2/25⌋ : ℝ) ≤ (w : ℝ) * c + w + 2/25 := Int.floor_le _
      have hc : (w : ℝ) * c - 2/25 ≤ (⌈(w : ℝ) * c - 2/25⌉ : ℝ) := Int.le_ceil _
      linarith
    have hZ : ⌊(w : ℝ) * c + w + 2/25⌋ - ⌈(w : ℝ) * c - 2/25⌉ < w + 1 := by
      exact_mod_cast hfc
    rw [hM, Int.card_Icc]
    omega
  have hMcard : ((M.card : ℕ) : ℝ) ≤ (w : ℝ) + 1 := by
    have : ((M.card : ℤ) : ℝ) ≤ ((w + 1 : ℤ) : ℝ) := by exact_mod_cast hMZ
    push_cast at this ⊢
    linarith
  -- assemble
  have htotal : (V.card : ℝ) ≤ ((M.card : ℕ) : ℝ) * ((4/25) * D / w + 1) := by
    rw [hsplit]
    push_cast
    calc (∑ m ∈ M, ((V.filter (fun j => F j = m)).card : ℝ))
        ≤ ∑ m ∈ M, ((4/25) * (D : ℝ) / w + 1) := by
          apply Finset.sum_le_sum
          intro m hm
          exact hfiber m hm
      _ = (M.card : ℝ) * ((4/25) * D / w + 1) := by
          rw [Finset.sum_const, nsmul_eq_mul]
  have hDw : (4/25 : ℝ) * D / w ≤ (4/25) * D := by
    rw [div_le_iff₀ hwR]
    nlinarith
  have hterm : (0 : ℝ) ≤ (4/25) * D / w + 1 := by positivity
  have hfinal : ((M.card : ℕ) : ℝ) * ((4/25) * D / w + 1) ≤ (8/25) * D + w + 1 := by
    have hmul : ((M.card : ℕ) : ℝ) * ((4/25) * D / w + 1)
        ≤ ((w : ℝ) + 1) * ((4/25) * D / w + 1) :=
      mul_le_mul_of_nonneg_right hMcard hterm
    have hexpand : ((w : ℝ) + 1) * ((4/25) * D / w + 1)
        = (4/25) * D + (4/25) * D / w + w + 1 := by
      field_simp
      ring
    linarith [hmul, hexpand.le, hexpand.ge, hDw]
  linarith [htotal, hfinal]

/-- **THE CLUSTER-GCD RUNG (|S| ≤ 3)**: a 12-family with no `2/25`-margin point
satisfies `(25 − 8|S|)·d ≤ 25·(Σ_S |vᵢ| + |S|)` for every `S` with `1 ≤ |S| ≤ 3`
and every positive common divisor `d` of the complement. -/
theorem gap_gcd_rung (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (S : Finset (Fin 12)) (hS1 : S.Nonempty) (hS3 : S.card ≤ 3)
    (d : ℤ) (hd0 : 0 < d) (hdvd : ∀ i, i ∉ S → d ∣ v i) :
    (25 - 8 * (S.card : ℤ)) * d ≤ 25 * (∑ i ∈ S, |v i| + S.card) := by
  classical
  -- the big side T, cited
  set T : Finset (Fin 12) := Finset.univ \ S with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ S := by
    intro i
    rw [hT, Finset.mem_sdiff]
    simp only [Finset.mem_univ, true_and]
  have hTcard : T.card ≤ 12 := by
    calc T.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ T)
      _ = 12 := by rw [Finset.card_univ, Fintype.card_fin]
  obtain ⟨t₀, hmargin⟩ := cite_margin_gen cite v T hTcard (fun i _ => hv i)
  -- the citation margin exceeds 2/25 (T.card ≤ 11 since S nonempty)
  have hT11 : T.card ≤ 11 := by
    have hScard1 : 1 ≤ S.card := Finset.card_pos.mpr hS1
    have : T.card = 12 - S.card := by
      rw [hT, Finset.card_univ_diff, Fintype.card_fin]
    omega
  have hmargin' : ∀ i ∈ T, ∀ m : ℤ, (2:ℝ)/25 < |(v i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i hi m
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (2:ℝ)/25 < 1/12 := by norm_num
      _ ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(v i : ℝ) * t₀ - m| := h
  -- every shifted copy is T-good; the violator at each copy is in S
  set D : ℕ := d.toNat with hD
  have hDpos : 0 < D := by rw [hD]; omega
  have hDd : (D : ℤ) = d := by rw [hD]; omega
  have hDdR : ((D : ℕ) : ℝ) = ((d : ℤ) : ℝ) := by exact_mod_cast hDd
  have hpick : ∀ j : ℕ, j ∈ Finset.range D →
      ∃ i ∈ S, ∃ m : ℤ, |(v i : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25 := by
    intro j _
    obtain ⟨i, m, hm⟩ := hnl (t₀ + (j : ℝ) / D)
    by_cases hiS : i ∈ S
    · exact ⟨i, hiS, m, hm⟩
    · exfalso
      have hiT : i ∈ T := (hTmem i).mpr hiS
      have hper := periodic_margin (v i) d hd0 (hdvd i hiS) t₀ (j : ℤ)
        ((2:ℝ)/25) (fun m' => le_of_lt (hmargin' i hiT m'))
      have := hper m
      rw [show ((j : ℤ) : ℝ) = (j : ℝ) by push_cast; ring, ← hDdR] at this
      linarith
  -- the classifier and the fiberwise count
  have hchoice : ∀ j : ℕ, ∃ i : Fin 12, j ∈ Finset.range D →
      (i ∈ S ∧ ∃ m : ℤ, |(v i : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25) := by
    intro j
    by_cases hj : j ∈ Finset.range D
    · obtain ⟨i, hiS, hm⟩ := hpick j hj
      exact ⟨i, fun _ => ⟨hiS, hm⟩⟩
    · exact ⟨0, fun hc => absurd hc hj⟩
  choose F hFspec using hchoice
  have hFmem : ∀ j ∈ Finset.range D, F j ∈ S := fun j hj => (hFspec j hj).1
  have hsplit : (Finset.range D).card
      = ∑ i ∈ S, ((Finset.range D).filter (fun j => F j = i)).card :=
    Finset.card_eq_sum_card_fiberwise hFmem
  -- per-fiber: inside the visit filter of |v i|
  have hfiber : ∀ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
      ≤ (8/25) * D + |v i| + 1 := by
    intro i _
    have hsub : (Finset.range D).filter (fun j => F j = i)
        ⊆ (Finset.range D).filter
          (fun j => ∃ m : ℤ, |((|v i| : ℤ) : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25) := by
      intro j hj
      rw [Finset.mem_filter] at hj ⊢
      obtain ⟨hjr, hjF⟩ := hj
      have hm := (hFspec j hjr).2
      rw [hjF] at hm
      exact ⟨hjr, (cov_abs_iff (v i) (t₀ + (j : ℝ) / D)).mp hm⟩
    have hcount := tooth_visit_count (|v i|) (abs_pos.mpr (hv i)) D hDpos t₀
    calc (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
        ≤ (((Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |((|v i| : ℤ) : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25)).card : ℝ) := by
          exact_mod_cast Finset.card_le_card hsub
      _ ≤ (8/25) * D + ((|v i| : ℤ) : ℝ) + 1 := hcount
  -- assemble over S
  have hsum : ((D : ℕ) : ℝ)
      ≤ (8/25) * D * S.card + (∑ i ∈ S, ((|v i| : ℤ) : ℝ)) + S.card := by
    have hcards : ((D : ℕ) : ℝ)
        = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
      have h0 : ((Finset.range D).card : ℝ)
          = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
        exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) hsplit
      rw [← h0, Finset.card_range]
    conv_lhs => rw [hcards]
    calc (∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ))
        ≤ ∑ i ∈ S, ((8/25) * (D : ℝ) + ((|v i| : ℤ) : ℝ) + 1) :=
          Finset.sum_le_sum hfiber
      _ = (8/25) * D * S.card + (∑ i ∈ S, ((|v i| : ℤ) : ℝ)) + S.card := by
          rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_const,
            Finset.sum_const]
          push_cast
          ring
  -- cast down to the integer conclusion
  have hcast : (((25 - 8 * (S.card : ℤ)) * d : ℤ) : ℝ)
      ≤ ((25 * (∑ i ∈ S, |v i| + S.card) : ℤ) : ℝ) := by
    have hdD : ((d : ℤ) : ℝ) = ((D : ℕ) : ℝ) := hDdR.symm
    push_cast
    push_cast at hsum
    rw [hdD]
    nlinarith [hsum]
  exact_mod_cast hcast

#print axioms tooth_visit_count
#print axioms gap_gcd_rung

end ClusterGcd
end LonelyRunner
