/-
# LRCStrictRuler — the strict rounding identity: a strictly-live ruler is a strict witness

klein's THM-685 (A1) is the *rounding identity*: for integer speeds, `c·v mod q` lies in the rounded band
iff `frac(v·c/q)` lies in the real band.  This file formalizes its **strict** form and its consequence:

  `StrictlyLive v q p`  :  `0 < p < q`  and  `∀ i,  q < 14·((vᵢ·p) % q) < 13·q`
  ⟹  `StrictWitness v`  at  `t₀ = p/q`,  with the *uniform* margin  `ε = 1/(14q)`.

**Integrality supplies the margin for free.**  The residues are integers, so `q < 14r` already forces
`q + 1 ≤ 14r`; no minimum over `i` is needed.  The whole content is the integer inequality

  `q + 1 ≤ 14·|r − j·q|`   for every `j : ℤ`,  given `0 ≤ r < q`, `q < 14r`, `14r < 13q`

(case `j ≤ 0`: `|r − jq| ≥ r` and `14r ≥ q+1`;  case `j ≥ 1`: `|r − jq| ≥ q − r` and
`14(q−r) = 14q − 14r ≥ 14q − (13q−1) = q+1`).  Dividing by `q` turns this into
`1/14 + 1/(14q) ≤ |vᵢ·(p/q) − m|` for every integer `m`.

**Consequence.**  Composing with `LRCStrictWitnessFloor`:

  `LRC14Statement  ⟸  LRCUpTo13  +  (every residual family has a strictly-live ruler)`,  kernel-pure.

The remaining obligation is now a **purely integer, Diophantine** statement: for each residual family
exhibit one modulus `q` and multiplier `p` with all thirteen residues strictly inside `(q/14, 13q/14)`.
No measure theory, no continuum, no Fourier.  It is NOT proved here — it is the open case of LRC(14).

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRCStrictWitnessFloor

namespace LonelyRunner
namespace LRC14Grand

/-- `v` has a **strictly-live ruler** at modulus `q` with multiplier `p`: `0 < p < q`, and every residue
`(vᵢ·p) % q` lies STRICTLY inside the safe band `(q/14, 13q/14)`. -/
def StrictlyLive (v : Fin 13 → ℤ) (q p : ℤ) : Prop :=
  0 < p ∧ p < q ∧ ∀ i, q < 14 * ((v i * p) % q) ∧ 14 * ((v i * p) % q) < 13 * q

/-- **The integer core of the strict rounding identity.**  A residue strictly inside the band stays
`(q+1)/14` away from every multiple of `q`.  Integrality is what produces the `+1`. -/
theorem int_band_bound {q r : ℤ} (hq : 0 < q) (hr0 : 0 ≤ r) (hrq : r < q)
    (hb1 : q < 14 * r) (hb2 : 14 * r < 13 * q) (j : ℤ) :
    q + 1 ≤ 14 * |r - j * q| := by
  rcases le_or_gt j 0 with hj | hj
  · have hs : j * q ≤ 0 := mul_nonpos_iff.mpr (Or.inr ⟨hj, le_of_lt hq⟩)
    rw [abs_of_nonneg (by omega)]
    omega
  · have hj1 : (1 : ℤ) ≤ j := hj
    have hs : q ≤ j * q := le_mul_of_one_le_left (le_of_lt hq) hj1
    rw [abs_of_nonpos (by omega)]
    omega

/-- **The strict rounding identity.**  A strictly-live ruler `(q, p)` gives a strict witness at
`t₀ = p/q`, with the uniform margin `ε = 1/(14q)`. -/
theorem strictWitness_of_strictlyLive {v : Fin 13 → ℤ} {q p : ℤ} (h : StrictlyLive v q p) :
    StrictWitness v := by
  obtain ⟨hp0, hpq, hband⟩ := h
  have hq : (0 : ℤ) < q := lt_trans hp0 hpq
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp0
  refine ⟨(p : ℝ) / (q : ℝ), ⟨by positivity, ?_⟩, 1 / (14 * (q : ℝ)), by positivity, ?_⟩
  · rw [div_lt_one hqR]; exact_mod_cast hpq
  · intro i m
    set r : ℤ := (v i * p) % q with hrdef
    set k : ℤ := (v i * p) / q with hkdef
    have hr0 : 0 ≤ r := Int.emod_nonneg _ (ne_of_gt hq)
    have hrq : r < q := Int.emod_lt_of_pos _ hq
    obtain ⟨hb1, hb2⟩ := hband i
    have hdec : (v i) * p = q * k + r := by
      rw [hrdef, hkdef]; exact (Int.ediv_add_emod (v i * p) q).symm
    have hc : (v i : ℝ) * (p : ℝ) = (q : ℝ) * (k : ℝ) + (r : ℝ) := by exact_mod_cast hdec
    -- phase decomposition:  vᵢ·(p/q) − m  =  r/q − (m − k)
    set j : ℤ := m - k with hjdef
    have hval : (v i : ℝ) * ((p : ℝ) / (q : ℝ)) - (m : ℝ)
        = ((r : ℝ) - (j : ℝ) * (q : ℝ)) / (q : ℝ) := by
      rw [hjdef]
      push_cast
      field_simp
      linarith [hc]
    rw [hval, abs_div, abs_of_pos hqR]
    -- integer core, cast to ℝ
    have hint : q + 1 ≤ 14 * |r - j * q| := int_band_bound hq hr0 hrq hb1 hb2 j
    have hintR : ((q : ℝ) + 1) ≤ 14 * |(r : ℝ) - (j : ℝ) * (q : ℝ)| := by
      exact_mod_cast hint
    -- 1/14 + 1/(14q) = (q+1)/(14q)
    have hrw : (1 : ℝ) / 14 + 1 / (14 * (q : ℝ)) = ((q : ℝ) + 1) / (14 * (q : ℝ)) := by
      field_simp
    rw [hrw, div_le_div_iff₀ (by positivity) hqR]
    nlinarith [hintR, hqR, abs_nonneg ((r : ℝ) - (j : ℝ) * (q : ℝ))]

/-- **THE REMAINING OBLIGATION, as an integer statement.**  Every residual family admits a strictly-live
ruler: one modulus `q` and multiplier `p` with all thirteen residues strictly inside `(q/14, 13q/14)`. -/
def StrictlyLiveSupply : Prop :=
  ∀ v : Fin 13 → ℤ, IsResidual v → ∃ q p : ℤ, StrictlyLive v q p

/-- A strictly-live-ruler supply gives the strict-witness supply. -/
theorem strictWitnessSupply_of_strictlyLiveSupply (h : StrictlyLiveSupply) : StrictWitnessSupply := by
  intro v hres
  obtain ⟨q, p, hl⟩ := h v hres
  exact strictWitness_of_strictlyLive hl

/-- **LRC(14) from the citation and one strictly-live ruler per residual family.**  Kernel-pure; the
entire remaining mathematical content is `StrictlyLiveSupply`, a purely integer statement. -/
theorem lrc14_of_strictlyLiveSupply (cite : LRCUpTo13) (h : StrictlyLiveSupply) :
    LRC14.LRC14Statement :=
  lrc14_of_strictWitnessSupply cite (strictWitnessSupply_of_strictlyLiveSupply h)

end LRC14Grand
end LonelyRunner
