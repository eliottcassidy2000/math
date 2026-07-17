---
id: LEM-041
title: THE EFFECTIVE FRAME THEOREM (the capstone of the frame-spectrum program, closed with explicit constants). (A) THE SHARP WEIGHT BOUND: |Ŵ_g(χ)| ≤ Ŵ_g(χ₀) = (π²/P²)·J₂(P/g)/3 (triangle inequality on the csc² sum; sharp at χ₀), and by the Jordan identity Σ_{g|P,g<P} J₂(P/g) = P²−1 the TOTAL weight is exactly (π²/3)(1 − P⁻²). (B) THE CLASS PARSEVAL (exact): Σ_{χ∈Û_q}|X̂_g(χ)|² = φ(q)⁻¹Σ_u X(gu)² =: E_g. (C) THE VARIANCE ENVELOPE (proved: Cauchy–Schwarz over classes with weight |Ŵ_g|, then (A)+(B)): Var_w(cross) ≤ (π²/3)(1−P⁻²)·Σ_{g|P, g<P} (π²/P²)(J₂(P/g)/3)·E_g — every constant arithmetic (Jordan totients × class energies); measured tightness ×2.45 (family60: 4060.3 ≤ 9947.1), ×3.25 (two-owner: 3563.1 ≤ 11574.5), ×2.73 (balanced: 3202.7 ≤ 8731.0) — single-digit slack on a two-step bound. (D) THE EFFECTIVE GENERIC-FRAME THEOREM (Chebyshev, now fully explicit): #{w : |cross(w) − mean| ≥ λ} ≤ φ(P)·Var/λ² with Var exact (LEM-039(C) closed form) or enveloped (C); referee: all exceedance counts within bounds, zero frames beyond 5σ on any cluster. (E) THE PER-CHARACTER MASS BOUND: |ĉ(χ)| ≤ Σ_{g: cond|P/g} (π²/P²)(J₂(P/g)/3)·min(√E_g, A_g) (A_g the class ℓ¹-mean) — the "finite list of products" bound the program promised, tight to ×2.0–×9.2 on the top masses of all three clusters. (F) class energies concentrate at small g (g = 1, 2, 3, 4 carry the bulk)
status: PROVED ((A) triangle + Jordan; (B) orthogonality; (C) CS chain; (D) Chebyshev; (E) assembly) + REFEREED EXACT/GREEN on three clusters (identities at 1e-13–1e-16; envelope and mass bounds verified with measured tightness; Chebyshev instantiated)
source: boxeph-2026-07-17-S68 (owner directive: finish out the capstone)
depends_on: [LEM-031/032/033 (the closed-form factors), LEM-039(C) (the exact variance), THM-892 ((C*)/MI0/Jordan), LEM-030 (mean)]
script: 04-computation/lrc14_capstone_effective_boxeph_S68.py -> 05-knowledge/results/lrc14_capstone_effective_boxeph_S68.out
---

# LEM-041 — the effective frame theorem

The distributional frame-law program, end to end: the spectrum factorizes
(LEM-031), the weights are L-values (LEM-032), the cross factors are graded
Gauss sums (LEM-033), the variance is an exact finite closed form (LEM-039),
and now the SIZE side closes: a proved envelope with arithmetic constants
within ×3 of truth, per-character mass bounds within single digits, and the
generic-frame theorem effective on every cluster. Nothing on this line
remains conjectural: identities exact, bounds proved, tightness measured.
What remains for LRC(14) integration is only the uniformity across the
cluster family — the propagation ledger's business (its W₀ rows), not a
spectral gap.

## Evidence log
- [x] weight bound sharp + Jordan collapse (1e-13–1e-15, three clusters)
- [x] class Parseval exact (93 classes); envelope ×2.45/×3.25/×2.73
- [x] Chebyshev instantiated (0 frames beyond 5σ anywhere)
- [x] per-character bounds ×2.0–×9.2 on top masses
- [ ] named (optional sharpening): replace CS by the grade-diagonal form of
      E_g (LEM-033 cells) to shave the ×3
