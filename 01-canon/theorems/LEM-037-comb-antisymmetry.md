---
id: LEM-037
title: THE COMB ANTISYMMETRY LAW — the reflection x ↦ 1−x (sections σ ↦ 6−σ, R_s ↦ R_{6−s}, classes c ↦ −c, signs swap, owner attribution PRESERVED since the co-owner lattice is reflection-invariant) gives: (A) N_c^{(e)}(R_s) = −N_{7−c mod 7}^{(e)}(R_{6−s}) for every owner, class, section (1421 exact checks, 6 clusters); (B) ν̂_e^{(s)}(a) = −conj(ν̂_e^{(6−s)}(a)); (C) σ_e(s) = −σ_e(6−s), hence σ_e(3) = 0 for EVERY owner — THE OWNER-IMBALANCE BASELINE OF LEM-030 VANISHES IDENTICALLY AT THE CENTRAL SECTION s = 3 (verified: all 21 owners across the five clusters with R_3 ≠ ∅); (D) |S^{(s)}(n)| = |S^{(6−s)}(n)| pointwise, hence Q_s(w) = Q_{6−s}(w) EXACTLY for every dilation frame (referee 10⁻¹⁵) — the entire frame-spectrum theory (THM-880/892, LEM-030..033) needs only s ∈ {0,1,2,3}
status: PROVED (the reflection argument of LEM-034(A), 3 lines per part) + REFEREED EXACT (1421 + 1421 instances; 21 owners at s = 3; Q-pairs at 1e-15)
source: boxeph-2026-07-17-S65 (self-directed; resolves LEM-034's σ_e-decomposition named item)
depends_on: [LEM-034(A) (the reflection), LEM-030 (the baseline it refines), THM-880 (bilinear referee)]
script: 04-computation/lrc14_comb_antisymmetry_boxeph_S65.py -> 05-knowledge/results/lrc14_comb_antisymmetry_boxeph_S65.out
---

# LEM-037 — comb antisymmetry

One reflection, four laws. The deepest consequence: the LEM-030 owner-
imbalance term −(π²/3)Σσ_e² is a function of s that vanishes at the center
— cross-term analysis at s = 3 starts from a ZERO baseline, and every
s-battery in the program can be halved (s ≤ 3 suffices, with s = 3
self-reflective).

## Evidence log
- [x] (A)/(B) 1421 + 1421 exact; (C) 21 owners; (D) 1e-15 on four clusters
- [ ] named: the s = 3 cross-spectrum — with σ ≡ 0, LEM-030's baseline is pure
      (K−1/6)-torsion; measure whether s = 3 frames are systematically better
