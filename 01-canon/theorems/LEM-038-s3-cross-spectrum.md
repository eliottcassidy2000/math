---
id: LEM-038
title: THE s = 3 CROSS-SPECTRUM (LEM-037's named open RESOLVED). (A) THE MASTER REFLECTION IDENTITY: S_e^{(6−s)}(n) = −conj(S_e^{(s)}(n)) for every owner and frequency (the reflection maps endpoint systems with positions p ↦ P−p, signs flipped, owners preserved) — subsumes LEM-037(D). (B) THE IMAGINARY LAW: at the central section, Re S_e^{(3)}(n) = 0 for ALL owners and ALL n — the s = 3 system is a PURE SINE SYSTEM, S_e(n) = 2i Σ_{pairs} ε sin(2πnp/P); referee: max |Re| ≈ 10⁻¹⁵ over 15 owners × up to 17640 frequencies, four clusters. (C) THE QUADRUPLE LAW: R_3 has no endpoint at 0 or 1/2 (at x = 1/2 every odd runner sits in the forbidden section 3; an all-even cluster collapses occupancy to {0}; x = 0 likewise), so the reflection pairing is fixed-point-free and intervals pair I ↔ 1−I with no self-symmetric interval: M(3) ≡ 0 (mod 4) — verified M = 20, 8, 4, 32, 40 across five geometries. (D) THE PHASE-FREE CROSS: S_e^{(3)} = i·a_e (a real) makes every owner cross-term real: X(m) = (Σa_e)² − Σa_e² — a GRAM FORM in the imaginary parts; interference without phases (exact, 10⁻¹⁵). (E) THE BASELINE COLLAPSE: σ_e(3) = 0 (LEM-037) kills the −(π²/3)Σσ² term: measured means at s = 0..3 (balanced) −3.01/−11.93/−2.38/−0.23 with σ²-terms −210.6/−59.2/−32.9/0 — the s = 3 cross-mean is pure (K−1/6)-torsion and nearly zero. (F) THE VERDICT (census, s = 0 vs 3): the s = 3 CROSS is systematically tamer — variance strictly smaller on every multi-owner cluster (3202.7→2520.4; 3563.1→0 EXACTLY, R_3 being single-owner; 940.1→81.0) with mean ≈ 0 — while ⟨Q⟩/M stays at the universal ~π²/3 (2.83→3.13; 3.00→3.13; 1.95→3.52) and frame-extremes are not uniformly better: s = 3 SUPPRESSES THE INTERFERENCE, NOT THE DIAGONAL. Parity law holds verbatim at s = 3 (proof never used s). (G) measured: N(h) = M for most h at s = 3 — off-diagonal signed coincidences mostly cancel under the pairing; N ≡ M (mod 2) always
status: PROVED ((A)/(B)/(D) one reflection computation; (C) 3 lines) + REFEREED EXACT (identities at 10⁻¹²–10⁻¹⁵; quadruple law ×5; censuses on three clusters both sections)
source: boxeph-2026-07-17-S66 (owner directive: the s = 3 cross-spectrum; find other things to prove)
depends_on: [LEM-037 (reflection + σ(3) = 0), LEM-030 (the collapsing baseline), S61 census machinery]
script: 04-computation/lrc14_s3_cross_spectrum_boxeph_S66.py -> 05-knowledge/results/lrc14_s3_cross_spectrum_boxeph_S66.out
---

# LEM-038 — the s = 3 cross-spectrum

One reflection identity, five laws. The central section is the PURE SINE
SECTION: spectra imaginary, cross a real Gram form, endpoint count ≡ 0 mod 4,
owner imbalances gone, and the measured cross-interference collapses (variance
down 21%–100%, mean ≈ 0). The frame program's cross-term analysis is easiest
at s = 3 — and the two-owner cluster's R_3 is literally interference-free
(single owner). What s = 3 does NOT improve: the diagonal ⟨Q⟩/M ≈ π²/3
(THM-892's universal law is section-blind).

## Evidence log
- [x] master identity + imaginary law (15 owners, 4 clusters, all n)
- [x] quadruple law ×5; phase-free Gram ×3; baseline table ×3
- [x] census s = 0 vs 3 ×3 (variance ↓ everywhere; parity law verbatim)
- [x] both named opens CLOSED (LEM-039, S67): the odd class law T_{-r} = -T_r (65,107 checks) + the conductor-resolved variance-drop tables
