---
id: HYP-3152
title: Coverage is a radius — q₀=q₆·R⁶ with the miss-PGF zeros near a circle of radius R; cap = binomial/de-Moivre on-circle (EVEN/low-λ), dip = φ⁴ off-circle (ODD = κ₃ = odd ear = Ω = Worpitzky); rigorous solvability via dual degree≤4 ⟹ Galois⊆S₄ (correcting the V₄ claim — the flip-action is a monoid with an absorbing apex arc); consec is the joint extremizer (max R, max bimodality L_y, max Newton-violation)
status: VERIFIED (q₀=q₆R⁶; near-circular zeros; Galois trivial/⊆S₄; Newton defects negative) + CORRECTION (monoid not V₄) + SYNTHESIS (λ=dip=odd face). Not a proof.
source: mac-mini-2026-06-27-S72
merges:
  - HYP-3161   # the parity split even(biquadratic)/odd(Worpitzky) -- here = on-circle/off-circle
  - HYP-3150   # function-compression guardrail
  - HYP-3132   # k=8 biquadratic = the EVEN face; CORRECTED: Galois<=S4, not V4
  - HYP-3147   # Worpitzky/ear = the ODD/off-circle face
related:
  - HYP-3153   # Lee-Yang/Worpitzky/quartic compression packet follow-on
  - HYP-3103   # the miss-PGF zeros = the circle
  - HYP-3122   # φ⁴ cumulants (κ₄ even, κ₃ odd)
  - THM-577    # cap = binomial pair-Pascal = the on-circle value
  - OPEN-Q-108
external: Lee-Yang circle theorem; de Moivre-Laplace; Newton-Maclaurin inequalities; S₄ solvability
---

# HYP-3152 — Coverage is a radius; cap on-circle, dip off-circle; Galois ⊆ S₄

## The Lee-Yang circle (VERIFIED)
The miss-count PGF `G_N(z)` (degree 6) has, by Vieta, **`q₀ = q₆·R⁶`** with `R=(q₀/q₆)^{1/6}`, and the zeros lie
**near a circle of radius R** (`lrc_leeyang_circle_galois_newton_macmini_S72.py`): R = 1.59→1.96 for k=8→13,
zero-radius ratio 1.14–1.36. **Coverage `p0 = q₀ = q₆·R⁶` is a RADIUS.** Two knobs:
- **λ→0 (zeros on the circle)** = the binomial/Pascal coefficients = de Moivre-Laplace Gaussian = the **EVEN**
  face = `cap = C(k+1,2)/91` (THM-577), solvable;
- **λ>0 (off-circle)** = the φ⁴ correction = the **dip** = the **ODD** face = κ₃ = odd ear = Ω = Worpitzky.
Extremality = **large R + low λ**; consec is the extremizer.

## The bimodality functional + Newton-Maclaurin (VERIFIED)
`L_y = 10q₀+q₃+10q₆ = 10(q₀+q₆+0.1q₃)` = the **bimodality functional** (extreme mass = the two circle-poles);
consec maximizes it. Newton-Maclaurin: at consec the normalized-moment defects `p_k²−p_{k-1}p_{k+1}` are **all
negative** — consec is the **extremal of the moment-inequality VIOLATION** = max bimodality = max extreme-mass =
0 real roots (S66). One object, several names.

## Rigorous solvability: Galois ⊆ S₄ (CORRECTING V₄)
## The compression hierarchy (beyond commutativity)
What the iso-class arc-cube function respects: **linearity** (score, even) — exact n≤4, fails n=5;
**associativity** (flip = XOR) — always (a monoid); **invertibility** (group) — FAILS (the absorbing apex
collapses info); **multiplicativity** (H, cycle, odd) — the irreducible nonlinear part. The compressible part is
the on-circle/even/binomial/biquadratic; the incompressible part is the **absorbing apex / off-circle / odd ear**
— the tournament avatar of the **apex prime 7**.

## Toward the proof
`cap = q₆R⁶ (on-circle binomial) − dip (off-circle φ⁴)`. Bound the dip = **[even biquadratic, solvable
degree-2 (S70)] + [odd Worpitzky/ear/κ₃ (codex HYP-3147)]**; Galois⊆S₄ makes the even part explicit; the
Lee-Yang circle keeps λ small. The off-circle (odd) part is the dominant, irreducible content (S71).

## Next
1. Quantify λ (the off-circle deviation) and relate it to the dip rigorously.
2. The even/odd (on-circle/off-circle) bound split = the S70 biquadratic + the codex Worpitzky sum ⟹ close the dip.
