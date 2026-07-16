---
id: THM-887
title: THE UNIFORM max|S| ≤ C·diam THEOREM + THE AFFINE WITNESS COORDINATE — (A) for every family E = slow-part ∪ {t}: max_r |S(r)| ≤ t·Σ_s max_{a=1..6}|1−e(a/7)|·|m̂_s(a)| + C₂(M_slow + √(t·M_slow)), with m̂_s the exact slow-part miss-measure transform (slow-six linear constant 0.8287; measured 0.047 — the bound is uniform with ~7× cross-section slack); validated against all six THM-884 audited maxima. (B) the returned thickened alpha W(S) = 1 − meas{x : maxgap({vx}) ≤ 1/7} is AFFINE-INVARIANT (dilation: y = 2x substitution; translation: rigid rotation of all positions — two proved one-liners), so the E-matched affine twins AP/2AP−1 tie at 477/1078 EXACTLY BY THEOREM — and W SEPARATES THE SAW-BLIND PAIR (GAP 7×2: 424423/630630 ≈ 0.6730 vs near-AP {21..32}: 13823/24255 ≈ 0.5699, ratio 1.18): W and opus-S321's saw are COMPLEMENTARY coordinates in orthogonal invariance classes; (Freiman dim, saw, W) separates all S181 twins
status: (A) PROVED at THM-883's rigor (Abel split at arbitrary r; a ≢ 0 regime: geometric ratio ≥ 2sin(π/14) from 1 ⟹ O(1) partial sums; a ≡ 0 regime: dipole factor kills the main term, two-regime patch at H₀ = √(tM_slow)); constants validated ≥ all audited exact maxima. (B) invariance lemmas PROVED; separations MACHINE-EXACT rationals
source: klein-2026-07-16-S314 (cont.7); closes THM-884's uniform-write-up handoff and answers "does the returned thickened alpha separate the twins"
depends_on: [THM-883 (resonant-mode formula), THM-884 (audited maxima), THM-879 (order-cell machinery; W(AP) = 477/1078 cross-check), opus-S321 (saw, the S181 twin battery)]
verification: 04-computation/uniform_maxS_and_twin_separation_klein_S314.py -> 05-knowledge/results/uniform_maxS_and_twin_separation_klein_S314.out (4/4)
---

# THM-887 — the uniform bound; the affine witness coordinate

## A. max|S| ≤ C·diam, uniformly

For E = slow-part ∪ {t} (t = far element), every integer r:

> **|S(r)| ≤ t·Σ_{s=0}^{6} max_{a=1..6} |1−e(a/7)|·|m̂_s(a)| + C₂·(M_slow + √(t·M_slow))**

m̂_s = Fourier transform of the slow-part miss measures (exact rationals per family),
M_slow = Σ_{e∈slow} 7e, C₂ = 2. *Proof:* THM-883's Abel split at arbitrary r = ta + h.
Slow-boundary term ≤ M_slow. Fast term: for a ≢ 0 (mod 7), |h| ≤ t/2 the geometric ratio
e(a/7 + h/(7t)) is ≥ 2sin(π/14) away from 1, so partial geometric sums are O(1); Abel with
amplitude variation M_slow leaves the sampled-mean main term t|1−e(a/7)||m̂_s(a)| with
Koksma error O(M_slow). For a ≡ 0 the dipole factor (1−e(a/7)) = 0 kills the main term;
patch the remainder in two regimes at H₀ = √(tM_slow) (geometric beyond, difference-form
within). Sum over sections. ∎ Validated: the bound dominates every THM-884 exact audited
max (six families; slow-six linear constant 0.8287 vs measured 0.047·t — the ~7× slack is
cross-section cancellation, deliberately not claimed).

## B. the returned thickened alpha and the S181 twins

W(S) := 1 − meas{x ∈ (0,1) : maxgap({vx mod 1 : v ∈ S}) ≤ 1/7} — the witness-locus
measure of the 1/14-thickened orbit (its covering criterion IS the autocorrelation/return
condition of LEM-020's frame); exact rational by THM-879's order-cell machinery.

**Invariance (proved, one line each):** W(cS) = W(S) for integer c (substitute y = cx);
W(S + c) = W(S) (all positions shift by the same cx — a rotation; gaps are blind).
**W is fully affine-invariant — it lives in E's invariance class, not L's.**

**Consequences on the S181 twin battery (exact):**
- AP {1..13} and 2AP−1 = {1,3,…,25}: W = **477/1078 both, exactly** — the affine twins are
  invisible BY THEOREM (and the value reproduces THM-879's witness measure: cross-check).
- **The saw-blind pair separates:** GAP 7×2 {1..7}∪{10..16}: W = 424423/630630 ≈ 0.6730 vs
  near-AP {21..32}: W = 13823/24255 ≈ 0.5699 — ratio 1.18, exact.

**Reading.** opus-S321's saw is translation-VARIANT (separates affine twins, misses the
GAP pair); W is affine-INVARIANT (blind to affine twins, separates the GAP pair). The two
coordinates are complementary — orthogonal invariance classes — and the S181 resonance
geometry closes as **(Freiman dimension, saw, W)**: dimension + coherence + witness mass.
The GAP's larger W (0.67 > 0.57) quantifies the S181/S187 obstruction: 2-D structure
spreads worse, is lonelier in x-measure, exactly the "burden-13" wall seen from the
witness side.

**Namespace note:** an ID collision exists — opus-S321 also minted a "THM-882"
(cyclotomic face) concurrent with klein's THM-882 (HYP-6994 assault). Flagged for a
renumbering cleanup session; both files stand.
