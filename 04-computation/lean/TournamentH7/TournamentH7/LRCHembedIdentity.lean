/-
  TournamentH7.LRCHembedIdentity — the exact reduction at the core of THM-527 Part A / hembed
  (opus-2026-07-09-S176).

  `hembed` (klein-S203, the SHARED blocker of the good-period AND density routes) realizes a
  fast-phase gap as a real lonely time.  Its EXACT analytic core, grounded numerically this session:
  for a speed `v = Vmax − e` (co-offset `e`, ruler `Vmax`), the loneliness margin of runner `v` at
  real time `τ` equals the clearance of the SLOW tooth `frac(e·τ)` by the FAST phase `frac(Vmax·τ)`:

      nearInt(v·τ) = nearInt(frac(Vmax·τ) − frac(e·τ)),   v = Vmax − e.

  So "runner `Vmax−e` is `1/14`-lonely at `τ`" ⟺ "the fast phase `frac(Vmax·τ)` is `>1/14` from the
  slow tooth `frac(e·τ)`" — the two-scale reduction, EXACT (no error; the finite-`Vmax` "coupling" is
  only the tooth WOBBLE `≤ spread/Vmax` as `τ` sweeps a ruler cell, opus-S176).  `nearInt` is the
  `LRCMreachConcrete` loneliness kernel `min (fract x) (1 − fract x)`.  This file proves the identity
  (kernel-pure); the remaining hembed content is the realization (a `τ` in the ruler cell putting the
  fast phase in the gap — mac-mini-S64's widest-arc, klein-S204).

  CONVERGENCE (same day): kps-S105 `LRCSlowFast.nearInt_speed_eq_phase_sub` + `drift_eq` (the wobble
  `e·φ/Vmax`) and klein-S204 `LRCCriterionC` (co-offset identity ⟹ `Mreach_ge_of_fastphase_clears`)
  independently formalized this SAME identity.  This file is an independent third verification (kept
  unimported to avoid aggregate redundancy); its unique companion is the opus-S176 NUMERICAL grounding
  (`lrc14_hembed_ruler_embedding_opus_S176`): hembed is TRUE (`M(S) ≥ 1/14` realized at some good
  period), the coupling is only the tooth wobble `≤ spread/Vmax`, and the naive embed is clean when
  `Vmax ≫ spread` — confirming kps-S105's "hembed is a FORMALIZATION gap, not open analysis."
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner.LRC14Concrete

/-- `nearInt` is invariant under integer shifts (`Int.fract` is). -/
theorem nearInt_add_int (x : ℝ) (n : ℤ) : nearInt (x + n) = nearInt x := by
  unfold nearInt
  rw [Int.fract_add_intCast]

/-- **The hembed core identity.**  For a speed `Vmax − e` (real `Vmax, e`), its loneliness margin at
`τ` is the fast-phase clearance of the slow tooth: `nearInt((Vmax − e)·τ) = nearInt(frac(Vmax·τ) −
frac(e·τ))`.  The exact two-scale reduction underlying THM-527 Part A. -/
theorem nearInt_speed_eq_fastPhase_clear (Vmax e τ : ℝ) :
    nearInt ((Vmax - e) * τ) =
      nearInt (Int.fract (Vmax * τ) - Int.fract (e * τ)) := by
  have h : (Vmax - e) * τ
      = (Int.fract (Vmax * τ) - Int.fract (e * τ)) + ((⌊Vmax * τ⌋ - ⌊e * τ⌋ : ℤ) : ℝ) := by
    simp only [Int.fract]; push_cast; ring
  rw [h, nearInt_add_int]

/-- Consequence, in loneliness form: runner `Vmax − e` is `≥1/14`-lonely at `τ` **iff** the fast phase
`frac(Vmax·τ)` clears the slow tooth `frac(e·τ)` by `≥1/14`. -/
theorem lonely_iff_fastPhase_clears (Vmax e τ : ℝ) :
    (1 : ℝ) / 14 ≤ nearInt ((Vmax - e) * τ) ↔
      (1 : ℝ) / 14 ≤ nearInt (Int.fract (Vmax * τ) - Int.fract (e * τ)) := by
  rw [nearInt_speed_eq_fastPhase_clear]

end LonelyRunner.LRC14Concrete
