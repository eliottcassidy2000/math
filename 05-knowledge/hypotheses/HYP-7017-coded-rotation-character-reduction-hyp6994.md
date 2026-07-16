# HYP-7017 — Uniform HYP-6994 via the coded-rotation character expansion

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S18)
**Verify-first:** stub until the script + results land this session.

## The idea (the owner-directed uniform residual, my angle)

Per THM-882-klein (assault), uniform HYP-6994 = flatness of the per-owner sign words
`u_e ∈ {−1,0,+1}^{Z_{7e}}` (positions j/(7e); the joint S(n) descends by klein's D3).
Structure: `u_e(j) = F(σ₁(j),…,σ₆(j); j mod 7)` where `σᵢ(j) = ⌊7{eᵢ j/(7e)}⌋` is runner
eᵢ's section at the j-th boundary of owner e, and F is the FIXED enter/leave sign rule of
the miss-pattern logic (alphabet Z₇^7, independent of the cluster sizes).

Expand F in Z₇^7 characters and each `ω^{c·σᵢ(j)}` in frequencies of the rotation
`j ↦ jeᵢ/(7e)`:
- **(Step 1, one line):** the section-mean of every NONTRIVIAL Z₇ character vanishes
  (Σ_{σ∈Z₇} ω^{cσ} = 0), so `g_c(0) = δ_{c=0}` — zero-frequency legs force trivial
  character legs.
- **(Step 2, compute):** the FIRST-ORDER (single-coordinate) Fourier coefficients of F.
  If they vanish (conjectured from klein's measured flatness C ≤ 14 with no e-drift),
  then every frequency peak of `û_e` needs ≥ 2 nonzero rotation legs, weight
  `≲ 1/(k₁k₂)`, and the peak height at m is controlled by the RESONANCE COUNT of
  `Σ kᵢeᵢ ≡ −m (mod 7e)` over the hyperbola lattice — the uniform O(M·polylog) form of
  HYP-6994 modulo an explicit relation-count bound.
- **(Step 3):** assemble + verify against klein's five scanned clusters; state the
  residual as the exact relation-count lemma (resonant clusters = small integer relations
  among speeds — meeting the repo's resonance taxonomy where it should).

This fuses klein's routes (i) (automaton) and (ii) (silent-boundary density) into one
Fourier-analytic mechanism; silent boundaries = F-zeros; the automaton = the coded
rotation orbit.

-> THM-882-klein (assault), THM-881, HYP-6994, THM-729; death-star-S18.
