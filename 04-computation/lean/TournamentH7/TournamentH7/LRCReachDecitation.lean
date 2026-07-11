/-
  TournamentH7.LRCReachDecitation — de-citing the ruler embedding through klein's
  transfer machinery (mac-mini-2026-07-09-S65 cont.27).

  cont.26 left `THM527ARulerEmbedding` as a citation carrying ALL of the reach
  content.  This file PROVES that content and shrinks the citation to pure
  arithmetic certificate existence:

  PROVED:
  * `minReach_ge_of_strictWitness` — a strict-margin time (kps's `StrictWitness`)
    HAS `minReach ≥ 1/14`: instantiate the ∀-integer margin at `⌊x⌋` and `⌊x⌋+1`
    to pin both `fract` and `1 − fract`.
  * `strictWitness_abs` — the |v|-normalization: a strict witness for the absolute
    speeds is one for the signed speeds (`|−x − m| = |x − (−m)|`; the ∀m margin is
    reflection-invariant).
  * `rulerEmbedding_of_certificateSupply` — the DE-CITATION REDUCTION: certificate
    existence ⟹ the ruler embedding, through klein's `strictWitness_of_cert`
    (S242: SafeIvStrict × 13 → grid point → strictly-live ruler → StrictWitness).

  NARROWED CITATION (`THM527ACertificateSupply`): every positive-witness shape
  realization admits a rational interval `[x/D, y/D]` with thirteen strict
  band certificates.  This is the classical program's OUTPUT SHAPE (the engine's
  intervals; klein's THM-693 constructive witnesses at large V; the census below):
  no measure theory, no reach geometry, no compactness remains in the citation —
  only integer inequalities.

  `lrc14_from_certificate_citations`: LRC(14) from THM-661 + the ≤7-arcs
  pigeonhole + certificate existence.
-/
import Mathlib
import TournamentH7.LRCReachCitation
import TournamentH7.LRCMeasureTransfer

set_option maxHeartbeats 800000

namespace LonelyRunner
namespace LRC14
namespace ReachDecitation

open MomentCitation ReachCitation LRC14Grand

/-- **A strict-margin time has `minReach ≥ 1/14`.**  The ∀-integer margin at
`⌊x⌋` pins `fract x ≥ 1/14 + ε`, at `⌊x⌋ + 1` pins `1 − fract x ≥ 1/14 + ε`. -/
theorem minReach_ge_of_strictWitness {v : Fin 13 → ℤ} (h : StrictWitness v) :
    ∃ τ : ℝ, (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v τ := by
  obtain ⟨t₀, _ht, ε, hε, hmargin⟩ := h
  refine ⟨t₀, ?_⟩
  unfold LRC14Concrete.minReach
  apply le_ciInf
  intro i
  unfold LRC14Concrete.nearInt
  have hfl : (⌊(v i : ℝ) * t₀⌋ : ℝ) ≤ (v i : ℝ) * t₀ := Int.floor_le _
  have hfu : (v i : ℝ) * t₀ < ⌊(v i : ℝ) * t₀⌋ + 1 := Int.lt_floor_add_one _
  have h1 := hmargin i ⌊(v i : ℝ) * t₀⌋
  have h2 := hmargin i (⌊(v i : ℝ) * t₀⌋ + 1)
  have e1 : |(v i : ℝ) * t₀ - (⌊(v i : ℝ) * t₀⌋ : ℝ)|
      = (v i : ℝ) * t₀ - ⌊(v i : ℝ) * t₀⌋ :=
    abs_of_nonneg (by linarith)
  have e2 : |(v i : ℝ) * t₀ - ((⌊(v i : ℝ) * t₀⌋ + 1 : ℤ) : ℝ)|
      = 1 - ((v i : ℝ) * t₀ - ⌊(v i : ℝ) * t₀⌋) := by
    push_cast
    rw [abs_of_nonpos (by linarith)]
    ring
  rw [e1] at h1
  rw [e2] at h2
  refine le_min ?_ ?_
  · show (1 : ℝ) / 14 ≤ (v i : ℝ) * t₀ - ⌊(v i : ℝ) * t₀⌋
    linarith
  · show (1 : ℝ) / 14 ≤ 1 - ((v i : ℝ) * t₀ - ⌊(v i : ℝ) * t₀⌋)
    linarith

/-- **The |v|-normalization**: a strict witness for the absolute speeds is a strict
witness for the signed speeds (the ∀-integer margin is reflection-invariant). -/
theorem strictWitness_abs {v : Fin 13 → ℤ}
    (h : StrictWitness (fun i => |v i|)) : StrictWitness v := by
  obtain ⟨t₀, ht, ε, hε, hmargin⟩ := h
  have hm' : ∀ (j : Fin 13) (m : ℤ), (1 : ℝ) / 14 + ε ≤ |(|v j| : ℝ) * t₀ - m| := by
    intro j m
    have := hmargin j m
    simpa using this
  refine ⟨t₀, ht, ε, hε, ?_⟩
  intro i m
  rcases abs_choice (v i) with habs | habs
  · have hR : |((v i : ℤ) : ℝ)| = ((v i : ℤ) : ℝ) := by
      rw [← Int.cast_abs, habs]
    have h := hm' i m
    rw [hR] at h
    exact h
  · have hR : |((v i : ℤ) : ℝ)| = -((v i : ℤ) : ℝ) := by
      rw [← Int.cast_abs, habs]
      push_cast
      ring
    have h := hm' i (-m)
    rw [hR] at h
    have heq : (-((v i : ℤ) : ℝ)) * t₀ - ((-m : ℤ) : ℝ)
        = -(((v i : ℤ) : ℝ) * t₀ - (m : ℝ)) := by
      push_cast
      ring
    rw [heq, abs_neg] at h
    exact h

/-- **THE NARROWED CITATION** — pure arithmetic certificate existence: every
positive-witness shape realization admits a rational interval `[x/D, y/D] ⊆ [0,1)`
of length `> 1/q` with thirteen strict band certificates for the absolute speeds.
This is the classical program's output shape (canon THM-527-A's good ruler periods
have strict rational sub-intervals; klein's THM-693 writes the large-`V` two-scale
witnesses down constructively; the windows/banks census covers bounded `V`). -/
def THM527ACertificateSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < witnessG2 (shapeOf v) →
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, MeasureTransfer.SafeIvStrict |v i| D x y

/-- **THE DE-CITATION REDUCTION**: certificate existence implies the ruler
embedding.  Klein's transfer (S242) turns the thirteen certificates into a
strictly-live ruler and a `StrictWitness` for `|v|`; the normalization and the
margin lemma finish. -/
theorem rulerEmbedding_of_certificateSupply (h : THM527ACertificateSupply) :
    THM527ARulerEmbedding := by
  intro v hv hpos
  obtain ⟨D, x, y, q, hD, hq, hx, hy, hlen, hcert⟩ := h v hv hpos
  have habs : ∀ i, 0 < |v i| := fun i => abs_pos.mpr (hv i)
  have hsw : StrictWitness (fun i => |v i|) :=
    MeasureTransfer.strictWitness_of_cert (fun i => |v i|) D x y q
      habs hD hq hx hy hlen hcert
  exact minReach_ge_of_strictWitness (strictWitness_abs hsw)

/-- **LRC(14) with the reach citation de-cited to arithmetic**: THM-661 + the
≤7-arcs pigeonhole + certificate existence.  All reach geometry, strictness
transfer, measure theory, and compactness are THEOREMS. -/
theorem lrc14_from_certificate_citations
    (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull)
    (hcerts : THM527ACertificateSupply) :
    LRC14Statement :=
  ReachCitation.lrc14_from_citations_only h661 hsmall7
    (rulerEmbedding_of_certificateSupply hcerts)

end ReachDecitation
end LRC14
end LonelyRunner
