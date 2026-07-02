/-
  TournamentH7.ForcedOverlap — MODULE 2 REMAINDER: the forced direction of the dangerous-pattern
  dichotomy (THM-605(ii)).
  kind-pasteur-2026-07-02-S8 (HYP-3965), transcribing mac-mini HYP-3868's recorded case-free
  construction (their Lean attempt yielded per the anti-stuck discipline; this is the queued
  fresh-context transcription).  Single-writer discipline: satellite file.

  THE STATEMENT.  For coprime P, Q with `1 < 2r(P+Q)` (at r = 1/14: exactly P+Q ≥ 8), EVERY
  phase θ admits a common near-integer point: ∃ x a b, |Px − a| < r ∧ |Qx − θ − b| < r.
  This is the anti-covering FLOOR direction — a bound on the adversary's power valid for all
  shifts (audit-safe ∀-IMP form): no phase choice can decorrelate a non-dangerous pattern to
  empty overlap.  Together with CombPatterns' `pattern_overlap_zero` (the avoidance direction,
  P+Q ≤ 7) it completes the dichotomy of THM-605(i).

  THE CONSTRUCTION (fully explicit, no case analysis):
   * the open window `(Pθ − r(P+Q), Pθ + r(P+Q))` has length `2r(P+Q) > 1`, so it contains an
     integer `z` STRICTLY inside; set `s := Pθ − z`, so `|s| < r(P+Q)`;
   * `u := s/(P+Q)` has `|u| < r`;
   * Bézout `cP + dQ = 1` gives `a := z·d`, `b := −z·c` with `Qa − Pb = z`;
   * `x := (a + u)/P` then satisfies `Px − a = u` and `Qx − θ − b = −u` EXACTLY.
-/
import TournamentH7.RatIntervals

namespace LonelyRunner
namespace ForcedOverlap

/-- An open interval of length `> 1` contains an integer strictly inside. -/
theorem exists_int_strictly_inside {α β : ℚ} (h : β - α > 1) :
    ∃ z : ℤ, α < (z : ℚ) ∧ (z : ℚ) < β := by
  refine ⟨⌊α⌋ + 1, ?_, ?_⟩
  · push_cast
    have := Int.floor_le α
    have := Int.lt_floor_add_one α
    linarith [Int.lt_floor_add_one α]
  · push_cast
    have h1 : (⌊α⌋ : ℚ) ≤ α := Int.floor_le α
    linarith

/-- **The forced direction (THM-605(ii), the anti-covering floor)**: a non-dangerous coprime
pattern has nonempty overlap at EVERY phase — with explicit rational witness and exact residues
`±u`. -/
theorem forced_overlap_exists (P Q : ℤ) (hP : 0 < P) (hQ : 0 < Q)
    (c d : ℤ) (hcd : c * P + d * Q = 1) (r θ : ℚ) (hr : 0 < r)
    (hbig : 1 < 2 * r * ((P : ℚ) + (Q : ℚ))) :
    ∃ (x : ℚ) (a b : ℤ), |(P : ℚ) * x - (a : ℚ)| < r ∧ |(Q : ℚ) * x - θ - (b : ℚ)| < r := by
  have hPQ : (0 : ℚ) < (P : ℚ) + (Q : ℚ) := by
    have h1 : (0 : ℚ) < (P : ℚ) := by exact_mod_cast hP
    have h2 : (0 : ℚ) < (Q : ℚ) := by exact_mod_cast hQ
    linarith
  have hPQne : ((P : ℚ) + (Q : ℚ)) ≠ 0 := ne_of_gt hPQ
  have hPne : ((P : ℚ)) ≠ 0 := by
    have : (0 : ℚ) < (P : ℚ) := by exact_mod_cast hP
    exact ne_of_gt this
  -- the window (Pθ − r(P+Q), Pθ + r(P+Q)) has length 2r(P+Q) > 1
  obtain ⟨z, hz1, hz2⟩ := exists_int_strictly_inside
    (α := (P : ℚ) * θ - r * ((P : ℚ) + (Q : ℚ)))
    (β := (P : ℚ) * θ + r * ((P : ℚ) + (Q : ℚ)))
    (by linarith)
  set s : ℚ := (P : ℚ) * θ - (z : ℚ) with hs
  have hsabs : |s| < r * ((P : ℚ) + (Q : ℚ)) := by
    rw [abs_lt]
    exact ⟨by linarith, by linarith⟩
  set u : ℚ := s / ((P : ℚ) + (Q : ℚ)) with hu
  have huabs : |u| < r := by
    rw [hu, abs_div, abs_of_pos hPQ, div_lt_iff₀ hPQ]
    linarith [hsabs]
  -- Bézout: c·P + d·Q = 1 (supplied as certificate data)
  set a : ℤ := z * d with ha
  set b : ℤ := -(z * c) with hb
  have hQaPb : (Q : ℚ) * (a : ℚ) - (P : ℚ) * (b : ℚ) = (z : ℚ) := by
    have hcdQ : (c : ℚ) * (P : ℚ) + (d : ℚ) * (Q : ℚ) = 1 := by exact_mod_cast hcd
    have hz1 : (z : ℚ) * ((c : ℚ) * (P : ℚ) + (d : ℚ) * (Q : ℚ)) = (z : ℚ) := by
      rw [hcdQ]; ring
    rw [ha, hb]
    push_cast
    ring_nf
    ring_nf at hz1
    linarith
  refine ⟨((a : ℚ) + u) / (P : ℚ), a, b, ?_, ?_⟩
  · -- Px − a = u exactly
    have hPx : (P : ℚ) * (((a : ℚ) + u) / (P : ℚ)) - (a : ℚ) = u := by
      field_simp
      ring
    rw [hPx]
    exact huabs
  · -- Qx − θ − b = −u exactly
    have hQx : (Q : ℚ) * (((a : ℚ) + u) / (P : ℚ)) - θ - (b : ℚ) = -u := by
      have hsu : ((P : ℚ) + (Q : ℚ)) * u = s := by
        rw [hu]; field_simp
      have hzs : (z : ℚ) = (P : ℚ) * θ - s := by rw [hs]; ring
      field_simp
      ring_nf
      ring_nf at hQaPb hsu hzs
      linarith
    rw [hQx, abs_neg]
    exact huabs

/-- **THM-605(ii) at the LRC(14) working band**: at `r = 1/14`, every coprime pattern with
`P + Q ≥ 8` has nonempty overlap at every phase.  (The complement of CombPatterns'
`pattern_overlap_zero`, which handles `P + Q ≤ 7`; together: the exact dichotomy.) -/
theorem forced_overlap_lrc14 (P Q : ℤ) (hP : 0 < P) (hQ : 0 < Q)
    (c d : ℤ) (hcd : c * P + d * Q = 1) (hPQ : 8 ≤ P + Q) (θ : ℚ) :
    ∃ (x : ℚ) (a b : ℤ),
      |(P : ℚ) * x - (a : ℚ)| < 1/14 ∧ |(Q : ℚ) * x - θ - (b : ℚ)| < 1/14 := by
  apply forced_overlap_exists P Q hP hQ c d hcd (1/14) θ (by norm_num)
  have h8 : (8 : ℚ) ≤ (P : ℚ) + (Q : ℚ) := by exact_mod_cast hPQ
  nlinarith

end ForcedOverlap
end LonelyRunner
