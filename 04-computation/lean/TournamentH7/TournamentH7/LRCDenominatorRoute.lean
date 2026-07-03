/-
  TournamentH7.LRCDenominatorRoute — LRC(14) AS A BOUNDED-DENOMINATOR SEARCH
  (kind-pasteur-2026-07-03-S29).

  The cleanest single-obligation framing of LRC(14).  Composing
    * `lrc14_of_covering_lonely` (the covering reduction — non-covering families are
      lonely at `t = 1/q` for the missed `q ≤ 14`), and
    * `lonely14_of_ratio` (the rational-witness sieve — `t = p/q` is lonely iff every
      `vᵢ·p` keeps integer-distance `≥ q/14` from every multiple of `q`),

  LRC(14) reduces to ONE statement:

      every covering family has a witness `p/q` with `q ≤ Q`.

  This is EXACTLY the route the machine-checked `n ≤ 13` proofs take (a finite search
  over bounded denominators), and kps-S28's `lrc14_denominator_bound/stress` gives the
  empirical `Q ≈ 35`: over 407 hard covering-compressed instances, all lonely at some
  `q ≤ 35`, concentrated at `q = 17, 19`, INDEPENDENT of speed magnitude.

  The obligation `hden` below is the whole remaining content of LRC(14).  It is a
  FINITE check — `p/q`-loneliness depends only on the speeds' residues mod `q` — but a
  large one, whose structured reduction (symmetry + covering + case trees) is the open
  computational frontier.  Kernel-pure; no `sorry`.
-/
import TournamentH7.LRCSpread13
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace LRC14

/-- **LRC(14) AS A BOUNDED-DENOMINATOR SEARCH.**  If every nonzero covering family
admits a rational witness `p/q` with `q ≤ Q` at which every runner keeps
integer-distance `≥ q/14` from every lattice point, then LRC(14) holds.  The
hypothesis `hden` is the entire remaining content of the conjecture, cast as a finite
denominator search. -/
theorem lrc14_of_covering_bounded_denom (Q : ℤ)
    (hden : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      ∃ p q : ℤ, 0 < q ∧ q ≤ Q ∧ ∀ i, ∀ k : ℤ, q ≤ 14 * |v i * p - k * q|) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  obtain ⟨p, q, hq, _hqQ, hres⟩ := hden v hv hcov
  exact ⟨(p : ℝ) / q, lonely14_of_ratio v p q hq hres⟩

/-- **The non-covering half is FREE.**  Without the covering hypothesis the witness is
immediate: a `q ∈ {2,…,14}` dividing no speed gives `t = 1/q` (denominator `≤ 14`).
So the bounded-denominator obligation lives ENTIRELY on covering families — the
`q ≥ 15` regime — matching the empirical concentration at `q = 17, 19`. -/
theorem bounded_denom_of_not_covering (v : Fin 13 → ℤ)
    (hnc : ¬ CoveringFamily v) :
    ∃ p q : ℤ, 0 < q ∧ q ≤ 14 ∧ ∀ i, ∀ k : ℤ, q ≤ 14 * |v i * p - k * q| := by
  -- ¬ covering: some q ∈ [2,14] divides no speed
  unfold CoveringFamily at hnc
  push_neg at hnc
  obtain ⟨q, hq2, hq14, hqnd⟩ := hnc
  refine ⟨1, (q : ℤ), by exact_mod_cast (by omega : 0 < q), by exact_mod_cast hq14, ?_⟩
  intro i k
  -- q ∤ v i, so |v i - k q| ≥ 1, hence 14·|v i·1 − k·q| ≥ 14 ≥ q
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast (by omega : 0 < q)
  have hndvd : ¬ (q : ℤ) ∣ v i := hqnd i
  have hne : v i * 1 - k * q ≠ 0 := by
    rw [mul_one, sub_ne_zero]
    intro hz
    exact hndvd ⟨k, by rw [hz]; ring⟩
  have h1 : (1 : ℤ) ≤ |v i * 1 - k * q| := Int.one_le_abs hne
  have hq14Z : (q : ℤ) ≤ 14 := by exact_mod_cast hq14
  calc (q : ℤ) ≤ 14 := hq14Z
    _ = 14 * 1 := by ring
    _ ≤ 14 * |v i * 1 - k * q| := by
        apply mul_le_mul_of_nonneg_left h1 (by norm_num)

end LRC14
end LonelyRunner
