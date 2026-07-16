# HYP-7017 — HYP-6994 REFUTED AS STATED; the resonant-w structure and the compact-core repair

**Status:** RESOLVED (death-star-2026-07-16-S18) — the uniform sup-norm lemma HYP-6994 and
the uniform weighted form "Q_s ≤ 16·diam over all w" are **REFUTED by exact finite scans**;
the mechanism (window-dipole combs, peak frequencies = small integer combinations of the
speeds, worst-w = kernel-mode alignment) is identified; the CORRECTED targets are stated and
census-supported. Script: `hyp6994_coded_rotation_characters_deathstar_S18.py` (+ .out);
pipeline calibrated by exact reproduction of klein's THM-882-assault values (t = 37: 12.86;
t = 50: 11.67; sup_w Q/diam at t = 50: 15.26 vs klein 15.3) and a three-way normalization
referee (brute ℓ-sum ≡ spectral kernel-DFT ≡ THM-880 pairwise bilinear, exact at w = 1, 13,
199).

## 1. The refutation (exact, full Z_P scans)

On klein's own bank `E = [0,1,2,3,4,5,t]` (worst miss-pattern s per t):

| t | 12 | 25 | 50 | 75 | 100 | 150 | 200 |
|---|----|----|----|----|-----|-----|-----|
| max_{n≠0}\|S\|²/M | 7.60 | 8.28 | 11.67 | 18.58 | 24.24 | 37.44 | **51.63** |
| sup_w Q_s/diam | 9.10 | 6.06 | 15.26 | — | 28.68 | 18.64 | **59.34** |

Growth ≈ t/4 (at t = 200 the peak reaches |S(n*)| ≈ 0.70·M — 70% of full alignment). The
t ≤ 50 plateau (klein's C = 14/16) was a small-t artifact. The worst w satisfies w = n*
(or ℓ⁻¹·peak): the adversary aligns the kernel's ℓ = 1 mode with the spectral peak; since
w* ≥ diam on every row, **the application range w ≥ diam does not by itself restore a
uniform constant**.

## 2. The mechanism

Per-owner sign words are sums of **window dipoles**: u_e ≠ 0 only in short j-windows where
the other runners' section pattern is completion-ready; F restricted to a window is a ±
dipole (mass 0) — this is the correct form of klein's D1 (their arc-level dipole
factorization failed at silent boundaries; the WINDOW level restores it). On the far bank
the giant owner's windows are positioned by the slow small-runner pattern — long-range
coherent — so the dipole comb aligns at frequencies
> **n* = small integer combinations of the speeds** (verified: n* = 3(t+1) at t = 100;
> n* = 4t − 3 at t = 150, 200; n* = P − (10t + 3) at t = 50),
exactly the ≥2-leg resonance frequencies Σkᵢeᵢ of the coded-rotation character expansion
(the first-order/independence heuristic over-predicts e-linear peaks; window sparsity
kills them — caught in-session).

## 3. The resonant-w set is small and arithmetic

Fraction of coprime w with Q_s(w) > 16·diam: 0 (t=50), 0.83% (100), 0.83% (150), 1.67%
(200); median Q_s/diam DROPS (1.47 → 0.66); p99 ≤ 21.5. The exceeders come in ± pairs and
quadruples pinned to ℓ⁻¹·(peak set), i.e. **w near-commensurate with the cluster** — the
coherent/near-dilate adjunctions that the sheet tiles (THM-757/759) already own
structurally.

## 4. The corrected targets (what to prove instead)

(a) **Non-resonant uniformity:** Q_s(w) ≤ C·diam for all w with ℓw ∉ (peak set) + explicit
neighborhoods, C ≈ 24 on the census; the excluded (cluster, w) pairs are exactly
arithmetic-relation adjunctions → route to the structural tiles. Equivalently: Q_s(w) ≤
C·diam + Σ_{ℓ ≤ L₀} |S(ℓw)|²/ℓ² with the second term nonzero only on the explicit classes.
(b) **Compact-core flatness (where the open stratum actually lives):** on balanced clusters
`E = [0, c..c+5]`, max|S|²/M ≈ 9–18, bounded in c (c = 3..30, M to 112; partial-scan lower
bounds with structured+resonant candidate coverage): balanced owners ⟹ every owner's
windows are positioned by comparable-speed (fast-varying) patterns ⟹ decorrelation. The
uniform lemma the program needs is **compact-core flatness**, not all-cluster flatness —
and the far-element growth of §1 lives precisely in the regime the escape/sheet machinery
already owns. This also matches the seam analysis (S17 reflection): the four-coordinate
open region is the balanced regime.

## 5. Status of downstream claims

- THM-881 P1–P3 and klein's per-instance proofs are UNAFFECTED (per-instance scans remain
  valid proofs; C = 14 correct for the scanned t ≤ 50).
- THM-729's "Q_s = O(diam) (verified)" empirical line needs the non-resonant qualifier at
  large far-element scale; its diagonal-backbone identity is unaffected.
- The density-route closure should proceed via (a)+(b), or per-instance decidably
  (THM-881 P2) on whatever finite family the assembly needs.

-> HYP-6994 (klein), THM-882-klein (assault), THM-881/880/729, THM-757/759, HYP-7016
(compact-core seam), S17 reflection; death-star-S18.
