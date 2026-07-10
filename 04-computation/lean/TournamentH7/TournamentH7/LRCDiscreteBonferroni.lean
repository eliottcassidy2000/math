/-
  TournamentH7.LRCDiscreteBonferroni — THM-671's HISTOGRAM BOUND in Lean
  (death-star-2026-07-09-S5, HYP-5780; klein-S210's named handoff "Lean next:
  LRCDiscreteBonferroni.lean (decide-shaped, no analysis)").

  THE OBJECT (klein-S210, THM-671).  Fix a modulus `q > 0` and count, for every nonzero
  multiplier `p ∈ (0, q)`, how many runners' residues `(v i · p) mod q` FALL OUTSIDE the safe
  band `[q/14, 13q/14]` — the coverage count `C(p) = bandCount`.  A LIVE multiplier (`C(p) = 0`)
  is an explicit loneliness witness `t = p/q` (kps-S114 `mreach_ge_of_pairsum_band`).  The
  quintic Bonferroni functional

      B5(v, q) = Σ_{d ≤ 5} (−1)^d · S_d,     S_d = Σ_{p ≠ 0} C(C(p), d)

  is computable in ONE pass over the coverage histogram, and `B5 ≤ liveCount` — so `B5 > 0` is a
  DECIDABLE loneliness certificate.  The depth ladder explains the certificate history: at the
  iid heuristic, `B1(13) = −6/7`, `B3(13) = −34/343` (why every union-bound/depth-2 ledger
  broke), `B5(13) = 2052/7⁵ = +0.1221` — depth 5 is the first truncation that clears 13 runners
  at threshold `1/7`.

  THE PROOF SHAPE.  Pointwise engine: mac-mini-S101's `BonferroniTruncation`
  (`odd_truncation_le_uncovered`: for odd `D`, `Σ_{d≤D} (−1)^d C(c,d) ≤ 1{c=0}` — THM-599's
  core, NOT re-proved here).  Sum it over the multipliers, swap the two finite sums, and the
  moment form `B5` appears; positivity then EXHIBITS a live `p`, and kps-S114's pair-sum
  dispatch turns it into `Mreach ≥ 1/14`.

  END-TO-END DEMO.  The AP extremal `{1,…,13}` at its unique live ruler `q = 14` (mac-mini
  THM-668: the AP kills every other ruler): coverage histogram `{0: 6, 1: 6, 6: 1}` — the six
  units of `(ℤ/14)ˣ` are live, the six even multipliers each cover only `v = 7`, and `p = 7`
  covers the six even speeds — so `B5 = 6·1 + 6·0 + 1·(−C(5,5)) = 5 > 0`, and the machine
  certifies `Mreach ≥ 1/14` for the equality extremal by `decide` + the dispatch.

  Kernel-pure target (plain `decide` for the demo; no `native_decide`, no `sorry`).
-/
import Mathlib
import TournamentH7.BonferroniTruncation
import TournamentH7.LRCPairSumDispatch

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Runner `i` passes the safe band at multiplier `p`, modulus `q`:
`(v i · p) mod q ∈ [q/14, 13q/14]`, as exact integer inequalities (the exact
shape `mreach_ge_of_pairsum_band` consumes). -/
def inBand (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13) : Prop :=
  (q : ℤ) ≤ 14 * ((v i * (p : ℤ)) % (q : ℤ)) ∧
    14 * ((v i * (p : ℤ)) % (q : ℤ)) ≤ 13 * (q : ℤ)

instance (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13) : Decidable (inBand v q p i) :=
  inferInstanceAs (Decidable (_ ∧ _))

/-- **The coverage count `C(p)`** (klein-S210): how many runners fail the safe band at
multiplier `p`.  `C(p) = 0` ⟺ `p` is a live multiplier ⟺ `t = p/q` is a loneliness witness. -/
def bandCount (v : Fin 13 → ℤ) (q p : ℕ) : ℕ :=
  (univ.filter fun i : Fin 13 => ¬ inBand v q p i).card

/-- The live-multiplier count `LM(q)`: nonzero multipliers whose coverage vanishes. -/
def liveCount (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 0).card

/-- The `d`-th factorial moment `S_d` of the coverage histogram (one pass). -/
def momentS (v : Fin 13 → ℤ) (q : ℕ) (d : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d

/-- **THM-671's quintic Bonferroni functional** `B5 = S₀ − S₁ + S₂ − S₃ + S₄ − S₅`. -/
def B5 (v : Fin 13 → ℤ) (q : ℕ) : ℤ :=
  ∑ d ∈ range 6, (-1 : ℤ) ^ d * (momentS v q d : ℤ)

/-- **The histogram bound (THM-671, klein-S210): `B5 ≤ LM(q)`.**  Sum mac-mini-S101's
pointwise odd-depth truncation over the multipliers and swap the finite sums. -/
theorem B5_le_liveCount (v : Fin 13 → ℤ) (q : ℕ) : B5 v q ≤ (liveCount v q : ℤ) := by
  have hswap : B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold B5 momentS
    calc ∑ d ∈ range 6, (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range 6, ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = ∑ p ∈ Finset.Ioo 0 q,
            ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) :=
          Finset.sum_comm
  rw [hswap]
  have hpoint : ∀ p ∈ Finset.Ioo 0 q,
      ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
        ≤ (if bandCount v q p = 0 then 1 else 0) := fun p _ =>
    TournamentH7.BonferroniTruncation.odd_truncation_le_uncovered (bandCount v q p) 5 rfl
  calc ∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
      ≤ ∑ p ∈ Finset.Ioo 0 q, (if bandCount v q p = 0 then (1 : ℤ) else 0) :=
        Finset.sum_le_sum hpoint
    _ = (liveCount v q : ℤ) := by
        unfold liveCount
        rw [Finset.sum_boole]

/-- **Positivity exhibits a live multiplier.** -/
theorem exists_live_of_B5_pos (v : Fin 13 → ℤ) (q : ℕ) (h : 0 < B5 v q) :
    ∃ p ∈ Finset.Ioo 0 q, bandCount v q p = 0 := by
  by_contra hno
  push_neg at hno
  have hzero : liveCount v q = 0 := by
    unfold liveCount
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    exact hno
  have hle := B5_le_liveCount v q
  rw [hzero] at hle
  omega

/-- **The decidable loneliness certificate (THM-671 → THM-668 dispatch): `B5(v, q) > 0` ⟹
`Mreach v ≥ 1/14`.**  A positive quintic Bonferroni functional at ANY modulus `q > 0` exhibits a
live multiplier `p`, and the rational `t = p/q` clears every runner (kps-S114).  Everything to
the left of `Mreach` is integer arithmetic — `decide`-shaped, no analysis. -/
theorem mreach_ge_of_B5_pos (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q) (h : 0 < B5 v q) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  obtain ⟨p, _, hlive⟩ := exists_live_of_B5_pos v q h
  apply mreach_ge_of_pairsum_band v (p : ℤ) (q : ℤ) (by exact_mod_cast hq)
  intro i
  have hempty : (univ.filter fun i : Fin 13 => ¬ inBand v q p i) = ∅ :=
    Finset.card_eq_zero.mp hlive
  have hi : ¬ ¬ inBand v q p i := by
    intro hni
    have hmem : i ∈ (univ.filter fun i : Fin 13 => ¬ inBand v q p i) :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ i, hni⟩
    rw [hempty] at hmem
    exact absurd hmem (Finset.notMem_empty i)
  exact not_not.mp hi

/-! ## The depth-general ladder (death-star-S6, HYP-5800; THM-675's escalation)

klein-S212/THM-675 found the B5 boundary: the 40→41 near-dilation (V = 260, primitive,
covering) has **B5 = 0/260 certified** while 244/260 moduli are exactly live — and **depth
escalation rescues it** (B7 = 100 ≈ LM = 102 at q = 280).  The odd-depth ladder is the same
machine at every depth, and it ENDS AT EXACTNESS: the coverage count is at most 13 (thirteen
runners), so the depth-13 truncation is the indicator ITSELF — `bonf 13 = liveCount` exactly,
a complete decidable live-ruler test with no Bonferroni loss. -/

/-- The odd-depth Bonferroni functional at general depth `D`:
`bonf D = Σ_{d ≤ D} (−1)^d S_d`.  `bonf 5 = B5` definitionally. -/
def bonf (D : ℕ) (v : Fin 13 → ℤ) (q : ℕ) : ℤ :=
  ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * (momentS v q d : ℤ)

theorem B5_eq_bonf5 (v : Fin 13 → ℤ) (q : ℕ) : B5 v q = bonf 5 v q := rfl

/-- **The histogram bound at every odd depth**: `bonf D ≤ LM(q)` for odd `D`. -/
theorem bonf_le_liveCount (D : ℕ) (hD : D % 2 = 1) (v : Fin 13 → ℤ) (q : ℕ) :
    bonf D v q ≤ (liveCount v q : ℤ) := by
  have hswap : bonf D v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold bonf momentS
    calc ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range (D + 1), ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = ∑ p ∈ Finset.Ioo 0 q,
            ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) :=
          Finset.sum_comm
  rw [hswap]
  calc ∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
      ≤ ∑ p ∈ Finset.Ioo 0 q, (if bandCount v q p = 0 then (1 : ℤ) else 0) :=
        Finset.sum_le_sum fun p _ =>
          TournamentH7.BonferroniTruncation.odd_truncation_le_uncovered (bandCount v q p) D hD
    _ = (liveCount v q : ℤ) := by
        unfold liveCount
        rw [Finset.sum_boole]

/-- The coverage count never exceeds the thirteen runners. -/
theorem bandCount_le_thirteen (v : Fin 13 → ℤ) (q p : ℕ) : bandCount v q p ≤ 13 := by
  unfold bandCount
  calc (univ.filter fun i : Fin 13 => ¬ inBand v q p i).card
      ≤ (univ : Finset (Fin 13)).card := Finset.card_filter_le _ _
    _ = 13 := by simp

/-- **The ladder ends at exactness (pointwise)**: for `c ≤ 13` the depth-13 truncation IS the
uncovered indicator — `Σ_{d ≤ 13} (−1)^d C(c,d) = 1{c = 0}` (the tail `C(c−1, 13)` vanishes
below the coverage cap). -/
theorem truncation13_eq_uncovered (c : ℕ) (hc : c ≤ 13) :
    ∑ d ∈ range 14, (-1 : ℤ) ^ d * (c.choose d) = (if c = 0 then 1 else 0) := by
  rcases Nat.eq_zero_or_pos c with rfl | hcpos
  · simp only [if_pos rfl]
    have hterms : ∀ d ∈ range 14, (-1 : ℤ) ^ d * ((0 : ℕ).choose d)
        = (if d = 0 then (1 : ℤ) else 0) := by
      intro d _
      match d with
      | 0 => simp
      | (e + 1) =>
        rw [Nat.choose_eq_zero_of_lt (Nat.succ_pos e)]
        simp
    rw [Finset.sum_congr rfl hterms, Finset.sum_ite_eq' (range 14) 0 (fun _ => (1 : ℤ))]
    simp
  · rw [if_neg (by omega : c ≠ 0),
      TournamentH7.BonferroniTruncation.partial_alternating_choose c 13 hcpos,
      Nat.choose_eq_zero_of_lt (by omega : c - 1 < 13)]
    simp

/-- **`bonf 13 = LM(q)` exactly** — the complete decidable live-ruler test.  At depth 13 the
Bonferroni truncation loses nothing, because no multiplier can be covered more than 13 times. -/
theorem bonf13_eq_liveCount (v : Fin 13 → ℤ) (q : ℕ) :
    bonf 13 v q = (liveCount v q : ℤ) := by
  have hswap : bonf 13 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 14, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold bonf momentS
    calc ∑ d ∈ range 14, (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range 14, ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = _ := Finset.sum_comm
  rw [hswap]
  calc ∑ p ∈ Finset.Ioo 0 q,
        ∑ d ∈ range 14, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ)
      = ∑ p ∈ Finset.Ioo 0 q, (if bandCount v q p = 0 then (1 : ℤ) else 0) :=
        Finset.sum_congr rfl fun p _ =>
          truncation13_eq_uncovered (bandCount v q p) (bandCount_le_thirteen v q p)
    _ = (liveCount v q : ℤ) := by
        unfold liveCount
        rw [Finset.sum_boole]

/-- **The certificate at every odd depth**: `bonf D (v, q) > 0` for odd `D` at any modulus
`q > 0` forces `Mreach ≥ 1/14`.  Depth is an analysis dial (5 generic, 7+ near-coherent,
13 exact); the conclusion is depth-free. -/
theorem mreach_ge_of_bonf_pos (v : Fin 13 → ℤ) (D : ℕ) (hD : D % 2 = 1) (q : ℕ)
    (hq : 0 < q) (h : 0 < bonf D v q) : (1 : ℝ) / 14 ≤ Mreach v := by
  have hlive : 0 < (liveCount v q : ℤ) := lt_of_lt_of_le h (bonf_le_liveCount D hD v q)
  have hlive' : 0 < liveCount v q := by exact_mod_cast hlive
  unfold liveCount at hlive'
  obtain ⟨p, hp⟩ := Finset.card_pos.mp hlive'
  have hzero : bandCount v q p = 0 := (Finset.mem_filter.mp hp).2
  apply mreach_ge_of_pairsum_band v (p : ℤ) (q : ℤ) (by exact_mod_cast hq)
  intro i
  have hempty : (univ.filter fun i : Fin 13 => ¬ inBand v q p i) = ∅ :=
    Finset.card_eq_zero.mp hzero
  have hi : ¬ ¬ inBand v q p i := by
    intro hni
    have hmem : i ∈ (univ.filter fun i : Fin 13 => ¬ inBand v q p i) :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ i, hni⟩
    rw [hempty] at hmem
    exact absurd hmem (Finset.notMem_empty i)
  exact not_not.mp hi

/-! ## The part-6 sockets (THM-673's target shape + dispersal counting core) -/

/-- **The aggregation pigeonhole (THM-673's target, consumable form)**: a positive
AGGREGATED Bonferroni mass `Σ_{q ∈ Q} bonf D (v, q) > 0` over any finite family of positive
moduli exhibits a certified modulus, hence `Mreach ≥ 1/14`.  This is the exact statement shape
the a-priori supply (part 6) proves: the mean over `q ∈ (V, 2V]` is positive, so SOME `q`
certifies — existence from an average of the CERTIFICATE, not of the event (no MISTAKE-129
max/mean fallacy: `bonf` itself is the certified lower bound). -/
theorem mreach_ge_of_bonf_sum_pos (v : Fin 13 → ℤ) (D : ℕ) (hD : D % 2 = 1)
    (Q : Finset ℕ) (hQ : ∀ q ∈ Q, 0 < q)
    (h : 0 < ∑ q ∈ Q, bonf D v q) : (1 : ℝ) / 14 ≤ Mreach v := by
  have hex : ∃ q ∈ Q, 0 < bonf D v q := by
    by_contra hno
    push_neg at hno
    exact absurd h (not_lt.mpr (Finset.sum_nonpos hno))
  obtain ⟨q, hqQ, hpos⟩ := hex
  exact mreach_ge_of_bonf_pos v D hD q (hQ q hqQ) hpos

/-- **The dispersal counting core (THM-673 (A), kernel-pure)**: a nonzero relation value
`N ≤ K·V` has at most `K` divisors in the modulus window `(V, 2V]` — each divisor `q` is
pinned by its cofactor `N/q ∈ [1, K]`, injectively.  This is the one-line counting that makes
the `m ≠ 0` resonances vanish from the aggregated mean as `V` grows. -/
theorem dvd_Ioc_card_le (N K V : ℕ) (hN : 0 < N) (hNK : N ≤ K * V) :
    ((Finset.Ioc V (2 * V)).filter fun q => q ∣ N).card ≤ K := by
  classical
  have hcard : ((Finset.Ioc V (2 * V)).filter fun q => q ∣ N).card
      ≤ (Finset.Icc 1 K).card := by
    apply Finset.card_le_card_of_injOn (fun q => N / q)
    · -- maps into [1, K]
      intro q hq
      obtain ⟨hqIoc, hdvd⟩ := Finset.mem_filter.mp hq
      obtain ⟨hVq, _⟩ := Finset.mem_Ioc.mp hqIoc
      have hqpos : 0 < q := lt_of_le_of_lt (Nat.zero_le V) hVq
      have hqleN : q ≤ N := Nat.le_of_dvd hN hdvd
      have hc1 : 1 ≤ N / q := (Nat.one_le_div_iff hqpos).mpr hqleN
      have hcq : N / q * q = N := Nat.div_mul_cancel hdvd
      have hcK : N / q ≤ K := by
        -- (N/q)·V < (N/q)·q = N ≤ K·V, then cancel V (weakly)
        have hlt : N / q * V < N / q * q :=
          mul_lt_mul_of_pos_left hVq (by omega : 0 < N / q)
        rw [hcq] at hlt
        have hltN : N / q * V < K * V := lt_of_lt_of_le hlt hNK
        have := lt_of_mul_lt_mul_right hltN (Nat.zero_le V)
        omega
      exact Finset.mem_Icc.mpr ⟨hc1, hcK⟩
    · -- injective: the cofactor determines the divisor
      intro q₁ hq₁ q₂ hq₂ hceq
      obtain ⟨hqIoc₁, hdvd₁⟩ := Finset.mem_filter.mp hq₁
      obtain ⟨_, hdvd₂⟩ := Finset.mem_filter.mp hq₂
      have h₁ : N / q₁ * q₁ = N := Nat.div_mul_cancel hdvd₁
      have h₂ : N / q₂ * q₂ = N := Nat.div_mul_cancel hdvd₂
      have hcpos : 0 < N / q₁ := by
        obtain ⟨hVq₁, _⟩ := Finset.mem_Ioc.mp hqIoc₁
        have hqpos₁ : 0 < q₁ := lt_of_le_of_lt (Nat.zero_le V) hVq₁
        exact (Nat.one_le_div_iff hqpos₁).mpr (Nat.le_of_dvd hN hdvd₁)
      have hceq' : N / q₁ = N / q₂ := hceq
      have hmul : N / q₁ * q₁ = N / q₁ * q₂ := by
        rw [h₁]
        calc N = N / q₂ * q₂ := h₂.symm
          _ = N / q₁ * q₂ := by rw [← hceq']
      exact Nat.eq_of_mul_eq_mul_left hcpos hmul
  simpa [Nat.card_Icc] using hcard

/-! ## End-to-end demo: the AP extremal at its unique live ruler -/

/-- The LRC(14) equality extremal `{1, …, 13}`. -/
def apSpeeds : Fin 13 → ℤ := fun i => (i : ℕ) + 1

/-- `B5(AP, 14) = 5 > 0`: the coverage histogram at the AP's unique live ruler `q = 14` is
`{0: 6, 1: 6, 6: 1}` (units live; even multipliers cover only `v = 7`; `p = 7` covers the six
even speeds), so `B5 = 6 − C(5,5) = 5`.  Kernel `decide` — 13 multipliers × 13 runners. -/
theorem b5_ap_fourteen_pos : 0 < B5 apSpeeds 14 := by decide

/-- **The histogram machine certifies the equality extremal**: `Mreach(AP) ≥ 1/14`, recovering
kps-S110's `mreach_AP_ge` through the THM-671 → THM-668 pipe with a `decide` certificate. -/
theorem mreach_AP_via_histogram : (1 : ℝ) / 14 ≤ Mreach apSpeeds :=
  mreach_ge_of_B5_pos apSpeeds 14 (by norm_num) b5_ap_fourteen_pos

/-! ## Axiom audit -/
#print axioms B5_le_liveCount
#print axioms exists_live_of_B5_pos
#print axioms mreach_ge_of_B5_pos
#print axioms b5_ap_fourteen_pos
#print axioms mreach_AP_via_histogram
#print axioms bonf_le_liveCount
#print axioms bonf13_eq_liveCount
#print axioms mreach_ge_of_bonf_pos
#print axioms mreach_ge_of_bonf_sum_pos
#print axioms dvd_Ioc_card_le

end LRC14Concrete
end LonelyRunner
