/-
  TournamentH7.LRCDenseCoreRelationTrap — THE RELATION TRAP on the chain-dense core:
  connecting kps's generic certification law (THM-934) to codex's
  `DenseCoreDissociatedB5Supply` (death-star-2026-07-16-S34, HYP-7171).

  kps's THM-934 mechanism (THM-930 leverage identity): a B5 certificate dies only under
  SYSTEMATIC coordinated relations — dilate chains, AP cores — never incidental ones.
  codex's `DenseCoreDissociatedB5Supply` asks for `∃ q, 0 < B5 v q` on primitive,
  dissociated, CHAIN-DENSE families.  This module proves the chain-dense certificate
  itself FORBIDS the killing structure above its dense pair, in kernel-checkable form,
  and hands that to the supplier:

  * `ladder_top_step` — one ladder step: above the dense pair, w(t) ≥ 3·w(t−1).
  * `no_low_mass_relation_above_pair` (A1) — NO relation `∑ coeff·w = 0` whose top
    support element t sits ≥ 2 above the dense pair and whose below-mass satisfies
    `∑_{i<t} |coeff i| ≤ 2` — this kills the DOUBLING-DILATE shapes (w_t = 2·w_s,
    coeff = (1, −2)) and every unit-pair shape (w_t = w_s + w_r) with top above the
    pair: the ladder's factor 3 beats mass 2.
  * `no_unit_relation_high` (A2) — NO unit-coefficient relation (∀ i, |coeff i| ≤ 1)
    whose top support element t is ≥ 4 above the dense pair: the geometric ladder sum
    `2·∑_{j+1≤i<t} w(i) ≤ w(t) − w(j+1)` plus the bottom-block bound `≤ 12·w(j+1)`
    give `∑_{i<t} w(i) < w(t)` once `w(t) ≥ 27·w(j+1)`.
  * `TrappedDenseCoreB5Supply` / `denseCoreDissociatedB5Supply_of_trapped` — the
    adapter: the B5 supplier may additionally assume the trap facts on the sorted
    speeds.  The analytic singular-series bridge (relation masses → B5) remains the
    named separate obligation (codex-S28's honest note in LRCB5RelationBudget).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDenseCoreEndgame

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- One ladder step above the dense pair: if `ChainDenseCore`'s upper-ladder clause
holds at witness `j`, then any element index `t` with `j + 2 ≤ t ≤ 12` satisfies
`3 · w(t−1) ≤ w(t)`. -/
theorem ladder_top_step (w : Fin 13 → ℤ) (j : Fin 12)
    (hladder : ∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ)
    (t : ℕ) (ht2 : (j : ℕ) + 2 ≤ t) (ht13 : t < 13) :
    3 * w ⟨t - 1, by omega⟩ ≤ w ⟨t, ht13⟩ := by
  have htm : t - 1 < 12 := by omega
  have hjk : j < (⟨t - 1, htm⟩ : Fin 12) := by
    have : (j : ℕ) < t - 1 := by omega
    exact this
  have h := hladder ⟨t - 1, htm⟩ hjk
  have hcs : ((⟨t - 1, htm⟩ : Fin 12)).castSucc = (⟨t - 1, by omega⟩ : Fin 13) := rfl
  have hsc : ((⟨t - 1, htm⟩ : Fin 12)).succ = (⟨t, ht13⟩ : Fin 13) := by
    apply Fin.ext
    show (t - 1) + 1 = t
    omega
  rwa [hcs, hsc] at h

/-- The sum over all thirteen indices splits as the top term plus the strictly-below
sum, when every strictly-above coefficient vanishes. -/
theorem sum_split_at_top (w coeff : Fin 13 → ℤ) (t : Fin 13)
    (hhigh : ∀ i : Fin 13, (t : ℕ) < (i : ℕ) → coeff i = 0) :
    (∑ i, coeff i * w i)
      = coeff t * w t
        + ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            coeff i * w i := by
  have hpart := Finset.sum_filter_add_sum_filter_not (Finset.univ : Finset (Fin 13))
    (fun i : Fin 13 => (i : ℕ) < (t : ℕ)) (fun i => coeff i * w i)
  have hnotset : Finset.univ.filter (fun i : Fin 13 => ¬ ((i : ℕ) < (t : ℕ)))
      = insert t (Finset.univ.filter (fun i : Fin 13 => (t : ℕ) < (i : ℕ))) := by
    ext i
    simp only [Finset.mem_insert, Finset.mem_filter, Finset.mem_univ, true_and, not_lt]
    constructor
    · intro h
      rcases Nat.eq_or_lt_of_le h with h1 | h1
      · exact Or.inl (Fin.ext h1.symm)
      · exact Or.inr h1
    · intro h
      rcases h with rfl | h
      · omega
      · omega
  have hzero' : ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (t : ℕ) < (i : ℕ)),
      coeff i * w i = 0 :=
    Finset.sum_eq_zero (fun i hi => by
      rw [hhigh i (Finset.mem_filter.mp hi).2]
      ring)
  have hins : ∑ i ∈ insert t (Finset.univ.filter
      (fun i : Fin 13 => (t : ℕ) < (i : ℕ))), coeff i * w i = coeff t * w t := by
    rw [Finset.sum_insert (by
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      omega), hzero']
    ring
  rw [hnotset, hins] at hpart
  linarith [hpart]

/-- **(A1) The low-mass relation trap.**  On a chain-dense family, no vanishing
integer relation has its top support element `t ≥ j+2` (strictly above the dense
pair) with total below-mass `∑_{i<t} |coeff i| ≤ 2`: the ladder's factor `3` beats
mass `2`.  Kills the doubling-dilate and unit-pair shapes above the pair. -/
theorem no_low_mass_relation_above_pair (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) (j : Fin 12)
    (hladder : ∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ)
    (coeff : Fin 13 → ℤ) (t : Fin 13)
    (ht2 : (j : ℕ) + 2 ≤ (t : ℕ))
    (htop : coeff t ≠ 0)
    (hhigh : ∀ i : Fin 13, (t : ℕ) < (i : ℕ) → coeff i = 0)
    (hmass : (∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
      |coeff i|) ≤ 2) :
    (∑ i, coeff i * w i) ≠ 0 := by
  have ht13 : (t : ℕ) < 13 := t.isLt
  have htm13 : (t : ℕ) - 1 < 13 := by omega
  set s : Fin 13 := ⟨(t : ℕ) - 1, htm13⟩ with hs
  have hstep : 3 * w s ≤ w t := by
    have := ladder_top_step w j hladder (t : ℕ) ht2 ht13
    have hteq : (⟨(t : ℕ), ht13⟩ : Fin 13) = t := Fin.ext rfl
    rwa [hteq] at this
  have hws : 0 < w s := hpos s
  have hwt : 0 < w t := hpos t
  have hbelow : ∀ i : Fin 13, (i : ℕ) < (t : ℕ) → w i ≤ w s := by
    intro i hi
    apply hmono
    show (i : ℕ) ≤ (t : ℕ) - 1
    omega
  intro hzero
  have hsplit := sum_split_at_top w coeff t hhigh
  have hbound : |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
      coeff i * w i| ≤ 2 * w s := by
    calc |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
        coeff i * w i|
        ≤ ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            |coeff i * w i| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            |coeff i| * w s := by
          apply Finset.sum_le_sum
          intro i hi
          rw [abs_mul, abs_of_pos (hpos i)]
          have hile := hbelow i (Finset.mem_filter.mp hi).2
          have habs : (0 : ℤ) ≤ |coeff i| := abs_nonneg _
          nlinarith [hpos i]
      _ = (∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            |coeff i|) * w s := by rw [Finset.sum_mul]
      _ ≤ 2 * w s := by nlinarith [hws]
  have htopabs : (1 : ℤ) ≤ |coeff t| := Int.one_le_abs htop
  have htopterm : w t ≤ |coeff t * w t| := by
    rw [abs_mul, abs_of_pos hwt]
    nlinarith [hwt]
  rw [hsplit] at hzero
  have heq : |coeff t * w t|
      = |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
          coeff i * w i| := by
    have h1 : coeff t * w t
        = -(∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            coeff i * w i) := by linarith
    rw [h1, abs_neg]
  nlinarith [hstep, hws]

/-- **(A2) The unit-coefficient trap.**  On a chain-dense family, no vanishing
unit-coefficient relation has its top support element `t ≥ j+4`: the geometric
ladder sum plus the bottom block stay strictly under `w(t)`. -/
theorem no_unit_relation_high (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) (j : Fin 12)
    (hladder : ∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ)
    (coeff : Fin 13 → ℤ) (t : Fin 13)
    (ht4 : (j : ℕ) + 4 ≤ (t : ℕ))
    (htop : coeff t ≠ 0)
    (hhigh : ∀ i : Fin 13, (t : ℕ) < (i : ℕ) → coeff i = 0)
    (hunit : ∀ i, |coeff i| ≤ 1) :
    (∑ i, coeff i * w i) ≠ 0 := by
  have ht13 : (t : ℕ) < 13 := t.isLt
  have hj1 : (j : ℕ) + 1 < 13 := by omega
  set p : Fin 13 := ⟨(j : ℕ) + 1, hj1⟩ with hp
  have hwp : 0 < w p := hpos p
  -- the geometric invariant, by induction up the ladder
  have hgeom : ∀ s : ℕ, (hs : s < 13) → (j : ℕ) + 1 ≤ s → s ≤ (t : ℕ) →
      2 * (∑ i ∈ Finset.univ.filter
        (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < s), w i)
      ≤ w ⟨s, hs⟩ - w p := by
    intro s
    induction s with
    | zero => intro hs h1 _; omega
    | succ n ih =>
      intro hs h1 hle
      rcases Nat.lt_or_ge n ((j : ℕ) + 1) with hn | hn
      · have hempty : Finset.univ.filter
            (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < n + 1) = ∅ := by
          apply Finset.filter_eq_empty_iff.mpr
          intro i _
          omega
        rw [hempty]
        have hpe : (⟨n + 1, hs⟩ : Fin 13) = p := Fin.ext (by
          show n + 1 = (j : ℕ) + 1
          omega)
        rw [hpe]
        simp
      · have hn13 : n < 13 := by omega
        have hstepn : 3 * w ⟨n, hn13⟩ ≤ w ⟨n + 1, hs⟩ := by
          have := ladder_top_step w j hladder (n + 1) (by omega) hs
          have he : (⟨n + 1 - 1, by omega⟩ : Fin 13) = ⟨n, hn13⟩ := Fin.ext (by
            show n + 1 - 1 = n
            omega)
          rwa [he] at this
        have hsplit : Finset.univ.filter
            (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < n + 1)
            = insert (⟨n, hn13⟩ : Fin 13) (Finset.univ.filter
                (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < n)) := by
          ext i
          simp only [Finset.mem_insert, Finset.mem_filter, Finset.mem_univ, true_and]
          constructor
          · intro ⟨hia, hib⟩
            rcases Nat.lt_or_ge ((i : ℕ)) n with h | h
            · exact Or.inr ⟨hia, h⟩
            · left
              apply Fin.ext
              show (i : ℕ) = n
              omega
          · intro h
            rcases h with rfl | ⟨hia, hib⟩
            · exact ⟨hn, Nat.lt_succ_self n⟩
            · exact ⟨hia, by omega⟩
        rw [hsplit, Finset.sum_insert (by
          simp only [Finset.mem_filter, Finset.mem_univ, true_and]
          omega)]
        have hihn := ih hn13 hn (by omega)
        nlinarith [hihn, hstepn]
  -- bottom block: at most twelve elements, each ≤ w(p)
  have hcard12 : (Finset.univ.filter
      (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)).card ≤ 12 := by
    have hsub : Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)
        ⊆ Finset.univ.erase t := by
      intro i hi
      have hi' := (Finset.mem_filter.mp hi).2
      rw [Finset.mem_erase]
      refine ⟨?_, Finset.mem_univ i⟩
      intro h
      rw [h] at hi'
      omega
    have hle := Finset.card_le_card hsub
    have herase : (Finset.univ.erase t).card = 12 := by
      rw [Finset.card_erase_of_mem (Finset.mem_univ t)]
      simp
    omega
  have hbottom : (∑ i ∈ Finset.univ.filter
      (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1), w i) ≤ 12 * w p := by
    have hle : (∑ i ∈ Finset.univ.filter
        (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1), w i)
        ≤ (Finset.univ.filter
            (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)).card • w p := by
      apply Finset.sum_le_card_nsmul
      intro i hi
      apply hmono
      show (i : ℕ) ≤ (j : ℕ) + 1
      have := (Finset.mem_filter.mp hi).2
      omega
    have hcast : ((Finset.univ.filter
        (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)).card : ℤ) ≤ 12 := by
      exact_mod_cast hcard12
    have hsmul : (Finset.univ.filter
        (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)).card • w p
        = ((Finset.univ.filter
            (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1)).card : ℤ) * w p :=
      nsmul_eq_mul _ _
    rw [hsmul] at hle
    nlinarith [hle, hcast, hwp]
  -- the 3-step magnification: w(t) ≥ 27·w(p)
  have h27 : 27 * w p ≤ w t := by
    have h1 : 3 * w p ≤ w ⟨(j : ℕ) + 2, by omega⟩ := by
      have := ladder_top_step w j hladder ((j : ℕ) + 2) (by omega) (by omega)
      have he : (⟨(j : ℕ) + 2 - 1, by omega⟩ : Fin 13) = p := Fin.ext (by
        show (j : ℕ) + 2 - 1 = (j : ℕ) + 1
        omega)
      rwa [he] at this
    have h2 : 3 * w ⟨(j : ℕ) + 2, by omega⟩ ≤ w ⟨(j : ℕ) + 3, by omega⟩ := by
      have := ladder_top_step w j hladder ((j : ℕ) + 3) (by omega) (by omega)
      have he : (⟨(j : ℕ) + 3 - 1, by omega⟩ : Fin 13)
          = ⟨(j : ℕ) + 2, by omega⟩ := Fin.ext (by
        show (j : ℕ) + 3 - 1 = (j : ℕ) + 2
        omega)
      rwa [he] at this
    have h3 : 3 * w ⟨(j : ℕ) + 3, by omega⟩ ≤ w ⟨(j : ℕ) + 4, by omega⟩ := by
      have := ladder_top_step w j hladder ((j : ℕ) + 4) (by omega) (by omega)
      have he : (⟨(j : ℕ) + 4 - 1, by omega⟩ : Fin 13)
          = ⟨(j : ℕ) + 3, by omega⟩ := Fin.ext (by
        show (j : ℕ) + 4 - 1 = (j : ℕ) + 3
        omega)
      rwa [he] at this
    have h4 : w ⟨(j : ℕ) + 4, by omega⟩ ≤ w t := by
      apply hmono
      show (j : ℕ) + 4 ≤ (t : ℕ)
      exact ht4
    nlinarith [h1, h2, h3, h4]
  -- assemble: ∑_{i<t} w(i) < w(t)
  have hwt : 0 < w t := hpos t
  have hsum_lt : (∑ i ∈ Finset.univ.filter
      (fun i : Fin 13 => (i : ℕ) < (t : ℕ)), w i) < w t := by
    have hsplit2 : Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ))
        = Finset.univ.filter
            (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1) ∪
          Finset.univ.filter
            (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < (t : ℕ)) := by
      ext i
      simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and]
      omega
    have hdisj : Disjoint
        (Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (j : ℕ) + 1))
        (Finset.univ.filter
          (fun i : Fin 13 => (j : ℕ) + 1 ≤ (i : ℕ) ∧ (i : ℕ) < (t : ℕ))) := by
      rw [Finset.disjoint_left]
      intro i hi hi'
      have h1 := (Finset.mem_filter.mp hi).2
      have h2 := (Finset.mem_filter.mp hi').2
      omega
    rw [hsplit2, Finset.sum_union hdisj]
    have hg := hgeom (t : ℕ) ht13 (by omega) (le_refl _)
    have hteq : (⟨(t : ℕ), ht13⟩ : Fin 13) = t := Fin.ext rfl
    rw [hteq] at hg
    nlinarith [hbottom, hg, h27, hwp]
  -- conclude
  intro hzero
  have hsplit := sum_split_at_top w coeff t hhigh
  have hbound : |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
      coeff i * w i|
      ≤ ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)), w i := by
    calc |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
        coeff i * w i|
        ≤ ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            |coeff i * w i| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)), w i := by
          apply Finset.sum_le_sum
          intro i _
          rw [abs_mul, abs_of_pos (hpos i)]
          have h1 := hunit i
          nlinarith [hpos i, abs_nonneg (coeff i)]
  have htopabs : (1 : ℤ) ≤ |coeff t| := Int.one_le_abs htop
  have htopterm : w t ≤ |coeff t * w t| := by
    rw [abs_mul, abs_of_pos hwt]
    nlinarith [hwt]
  rw [hsplit] at hzero
  have heq : |coeff t * w t|
      = |∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
          coeff i * w i| := by
    have h1 : coeff t * w t
        = -(∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
            coeff i * w i) := by linarith
    rw [h1, abs_neg]
  nlinarith [hsum_lt]

/-- **The trapped-supply obligation**: `DenseCoreDissociatedB5Supply`'s clauses PLUS
the two relation-trap facts on the sorted speeds — the supplier may assume no
low-mass relation tops out above the dense pair and no unit relation tops out `≥ 4`
above it.  This is the formal interface through which kps's THM-934 generic law
("only systematic top structure kills B5") meets the chain-dense core. -/
def TrappedDenseCoreB5Supply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g → nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), ∃ j : Fin 12,
      Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|) ∧
      (∀ coeff : Fin 13 → ℤ, ∀ t : Fin 13, (j : ℕ) + 2 ≤ (t : ℕ) →
        coeff t ≠ 0 → (∀ i : Fin 13, (t : ℕ) < (i : ℕ) → coeff i = 0) →
        (∑ i ∈ Finset.univ.filter (fun i : Fin 13 => (i : ℕ) < (t : ℕ)),
          |coeff i|) ≤ 2 →
        (∑ i, coeff i * |v (σ i)|) ≠ 0) ∧
      (∀ coeff : Fin 13 → ℤ, ∀ t : Fin 13, (j : ℕ) + 4 ≤ (t : ℕ) →
        coeff t ≠ 0 → (∀ i : Fin 13, (t : ℕ) < (i : ℕ) → coeff i = 0) →
        (∀ i, |coeff i| ≤ 1) →
        (∑ i, coeff i * |v (σ i)|) ≠ 0)) →
    ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q

/-- **The adapter**: a trapped supplier discharges codex's
`DenseCoreDissociatedB5Supply` — the trap facts are THEOREMS on the chain-dense
core, so handing them over costs nothing. -/
theorem denseCoreDissociatedB5Supply_of_trapped
    (h : TrappedDenseCoreB5Supply) : DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcert
  obtain ⟨σ, hmono, j, hjbad, hladder, hfee⟩ := hcert
  set w : Fin 13 → ℤ := fun i => |v (σ i)| with hw
  have hpos : ∀ i, 0 < w i := fun i => abs_pos.mpr (hv (σ i))
  apply h v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc
  refine ⟨σ, j, hmono, ⟨j, hjbad, hladder, hfee⟩, ?_, ?_⟩
  · intro coeff t ht2 htop hhigh hmass
    exact no_low_mass_relation_above_pair w hpos hmono j hladder coeff t
      ht2 htop hhigh hmass
  · intro coeff t ht4 htop hhigh hunit
    exact no_unit_relation_high w hpos hmono j hladder coeff t
      ht4 htop hhigh hunit

/-- **LRC(14) through the trapped supply** — the composed headline: citation +
four-detuned dispatch + B5 on the trapped chain-dense dissociated core. -/
theorem lrc14_from_four_detuned_and_trapped_B5
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchFour)
    (hB5 : TrappedDenseCoreB5Supply) : LRC14.LRC14Statement :=
  lrc14_from_four_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_trapped hB5)

/-! ## Axiom audit -/
#print axioms no_low_mass_relation_above_pair
#print axioms no_unit_relation_high
#print axioms denseCoreDissociatedB5Supply_of_trapped
#print axioms lrc14_from_four_detuned_and_trapped_B5

end LRC14Grand
end LonelyRunner
