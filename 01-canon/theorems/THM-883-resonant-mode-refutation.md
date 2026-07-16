---
id: THM-883
title: THE RESONANT-MODE THEOREM — the miss-pattern automaton induction, executed: conditioning on the fastest runner t and Abel-summing its boundaries yields S(t·a) = t·(1−e(a/7))·m̂_s(a) + O(M_slow), where m̂_s(a) is the Fourier transform of the EXACT 6-runner miss-measures A_c(s) + B(s) (computed: e.g. s=6: B = 2/35, A = (0, 1/28, 11/210, 5/84, 1/35, 1/35), max|m̂| = 0.1556); the formula is VALIDATED to O(M_slow) across t = 25..400 (|S| predicted 62.81 vs actual 65.12 at t = 400) — and it REFUTES the uniform lemma: max|S|²/M ~ t·|m̂|² → ∞ (HYP-6994 sup-version FALSE), and sup_w Q_s ~ diam² on the resonant w-classes w ≡ t·a·ℓ̄ (mod P) (15.3 → 28.7 ·diam from t = 50 → 100): THM-729's empirical O(diam) held because tested w avoided those classes. THE REPAIRED LAW: Q_s = O(diam) OFF the explicitly-listed resonant classes; ON them Q_s ≤ C·diam²·|m̂|²/ℓ² + C·diam
status: the resonant-mode formula PROVED (Abel/Koksma; the fast-boundary Weyl sum degenerates at n = ta since e(nj/(7t)) = e(aj/7) is constant per residue class; Koksma error V(F_c)·t·D*_t = O(M_slow)) + VALIDATED numerically to the predicted precision; the two refutations MACHINE-EXACT; the repaired off-resonance law CONJECTURED with the resonant classes explicit
source: klein-2026-07-16-S314 (cont.4); owner directive "prove the uniform lemma via the miss-pattern automaton induction" — the induction WORKS and the honest outcome is a refutation-with-replacement
depends_on: [THM-882 (w-freeness, per-instance scans), THM-881, THM-880, THM-729]
verification: 04-computation/automaton_induction_resonant_mode_klein_S314.py -> 05-knowledge/results/automaton_induction_resonant_mode_klein_S314.out (4/4)
---

# THM-883 — the resonant mode: the induction proves the formula, the formula kills the lemma

## 1. The induction step (proved)

Write f = Σ_c 1[c_t = c]·F_c(σ_slow) (conditioning on the fast runner). d(uv)-splitting gives
S(n) = Σ_c [∫F_c e(nx) d1[c_t=c] + ∫1[c_t=c] e(nx) dF_c]. The second term is ≤ M_slow
always. In the first, the fast boundaries sit at j/(7t) and at the t-RESONANT frequencies
n = t·a (a = 1..6): e(nj/(7t)) = e(aj/7) is CONSTANT on each residue class j mod 7 — the
Weyl sum degenerates into sampled means of F_c on the t-grid, and Koksma (V(F_c) = O(M_slow),
grid discrepancy 1/t) gives

> **S(t·a) = t·(1 − e(a/7))·m̂_s(a) + O(M_slow),  m̂_s(a) = Σ_{c≠s}(A_c(s) + B(s))·e(ac/7),**

with A_c(s) = meas{slow six miss exactly {s,c}}, B(s) = meas{slow six miss exactly {s}} —
fixed exact rationals of the slow cluster. Validated: predicted vs actual |S(ta)| =
(3.07, 8.17), (5.96, 9.84), (14.11, 16.31), (28.13, 31.75), (62.81, 65.12) at
t = 25..400 — errors O(M_slow) as proved.

## 2. The two refutations

- **HYP-6994 (uniform sup-version) is FALSE:** since some m̂_s(a) ≠ 0 for the slow six
  {0..5} (max 0.1556 at s = 6, a = 1, 6), max_{n≠0}|S(n)|² ≥ |S(ta)|² ~ t²|m̂|²·|1−e(a/7)|²,
  so max|S|²/M ~ t → ∞. The per-instance bounds (THM-882, C = 14 at t ≤ 50) were the
  pre-asymptotic regime.
- **Uniform Q_s = O(diam) over ALL w is ALSO FALSE:** choosing w in a resonant class
  (ℓw ≡ ta mod P for small ℓ) pushes the t²-mode through the 1/ℓ² weight:
  sup_w Q_s/diam = 15.3 (t=50) → 28.7 (t=100) — diam-linear growth of the ratio, i.e.
  sup_w Q_s ~ diam². THM-729's empirical O(diam) (w = 997-type samples) stands as an
  OFF-RESONANCE law.

## 3. The repaired law (the replacement conjecture)

> **Q_s(w) ≤ C·diam for all w OUTSIDE the resonant classes {w ≡ t·a·ℓ̄ (mod P/gcd) :
> a = 1..6, ℓ ≤ ℓ₀}; on a resonant class at depth ℓ, Q_s ≤ C·diam²·|m̂|²/ℓ² + C·diam.**

The resonant classes are finitely many, explicit, and computable per cluster family (the
m̂ vector is a one-page exact computation). The density route must check that its actual
peels w avoid the resonant classes — or route them through the existing resonant-peel
handling (S275); this is the named follow-up. The automaton induction did exactly what an
induction should: it identified the invariant that survives (off-resonance cancellation)
and the exact obstruction that does not (the miss-measure mode m̂).
