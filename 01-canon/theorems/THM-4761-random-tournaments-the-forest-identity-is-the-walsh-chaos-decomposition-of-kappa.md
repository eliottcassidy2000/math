---
id: THM-4761
title: "Random tournaments: the forest identity is the Walsh chaos decomposition of kappa. Under independent fair arcs, kappa = H 2^(n-1)/n! = sum_L N_L/(n)_L splits into orthogonal homogeneous chaoses kappa_L = N_L/(n)_L, with E N_L^2 = sum_p 2^p (n)_(L+p) C(L/2 - 1, p - 1)/p!. Consequently Var kappa = 2/n - 14/(3 n^3) - 58/(3 n^4) + ... (no n^-2 term), and the score chaos kappa_2 = 1 - W/(n(n-1)) = 8 c_3/(n(n-1)) - (n-2)/3 (W = sum_v (2 s_v - n + 1)^2, c_3 = number of cyclic triangles) carries the fraction 1 - 1/n + (4/3) n^-2 + O(n^-3) of the variance; with Moran's central limit theorem for c_3, sqrt(n)(kappa - 1) => N(0, 2)."
status: >
  PROVED (opus-2026-10-22-S54; elementary; independently audited 2026-10-22,
  04-computation/experiments/loop_nonregular_20261022/audit_S54/REPORT.md, which corrected the sign in (1) and supplied the
  classical link in (5)). (1)-(2) are exact identities. Brute force over all tournaments gives E kappa^2 = 4/3, 4/3, 79/60,
  58/45 for n = 3..6 (random_chaos.out) and 635/504 at n = 7 (audit). (3) is a finite expansion of (2). Sample check
  (FINITE, uniform random tournaments with exact H): corr(kappa, 1 + kappa_2) = 0.9628 at n = 12 (1500 samples) and 0.9705
  at n = 14 (800), against sqrt(share) = 0.9649 and 0.9687; the audit's 20,000 samples at n = 12 give 0.9654.
  (5) is CONDITIONAL on the CITED central limit theorem for c_3 (Kendall-Babington Smith 1940, Moran 1947; Moon, Topics on
  Tournaments, section 6, Theorem 6).
  PRIOR ART:
  - The variance of the score chaos is classical. It is Kendall-Babington Smith's (3/16) C(n,3) for c_3 (Moon, section 6).
  - Concentration (Var H = o((E H)^2) for Hamilton cycles, "routine calculations") is stated by Wormald, Tournaments with
    many Hamilton cycles (arXiv:2610.16245), concluding remarks.
  - Janson, Hamilton cycles in a random tournament, Random Structures Algorithms 6 (1995) 213-219: the abstract (search
    snippet only, full text not read) states asymptotic normality of the Hamilton-cycle count.
  - The orthogonality is the one Session30 used for E kappa = 1 and Var W
    (05-knowledge/results/session30_nonregular_loops_20261022/session30-root/SCORE_SCALE_AND_RANDOM_CONTROLS.md, section 2).
  - The exact second-moment formula (2), the n^-3 coefficient and the score share were NOT FOUND (targeted web search and
    the audit's reading of Moon's book); not claimed new beyond that.
session: opus-2026-10-22-S54
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/loop_nonregular_20261022/random_chaos.py (+ random_chaos.out): formula, brute force n <= 6, samples
  - 04-computation/experiments/loop_nonregular_20261022/random_tournaments.py (+ random_tournaments.out): loop-expansion reading
  - 04-computation/experiments/loop_nonregular_20261022/ar_tool.c: uniform random tournaments with exact H (mode random)
  - 04-computation/experiments/loop_nonregular_20261022/audit_S54/ (chaos_check.out, ekappa2_brute.c, rand_reading.out): independent audit
related:
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 05-knowledge/hypotheses/HYP-9369-almost-regular-loop-predictor-at-scale-complete-n10-census-and-the-joint-remainder.md
---

# THM-4761 — the forest identity as a chaos decomposition

Notation as in THM-4757: `S = A − Aᵀ`, `κ(T) = H(T)·2^(n−1)/n!`, and `N_L(S)` is the signed count of directed linear forests of `K_n` with `L` arcs. For a uniformly random tournament, every arc of `K_n` is oriented by an independent fair coin. Fix a reference orientation `i → j` for every pair `i < j`, and let `s_ij = S_ij = ±1` be the corresponding independent signs.

## 1. Statement

**(1) Chaos decomposition.** A directed path and its reversal carry the products `∏ S_e` and `(−1)^ℓ ∏ S_e`, where `ℓ` is the length. Odd paths therefore cancel against their reversals, and
```
N_L = Σ_F 2^(p(F)) ε_F χ_F,     χ_F = ∏_(ij ∈ F, i < j) s_ij,
ε_F = (−1)^(number of arcs traversed from a larger to a smaller label).
```
* The sum runs over the undirected linear forests `F` of `K_n` with `L` edges whose `p(F)` paths all have even length.
* Each path is traversed in either direction. `ε_F` is well defined because every path has even length.
* Distinct forests have distinct edge sets, so the Walsh characters `χ_F` are orthonormal: `E χ_F χ_F' = [F = F']`.

So `κ = Σ_L κ_L`, with `κ_L = N_L/(n)_L`, is the decomposition of `κ` into homogeneous chaoses of degree `L`. In particular `E κ = 1` and `E κ_L κ_L' = 0` for `L ≠ L'`.

(The audit corrected an earlier display without `ε_F`. For transitive `T_3` it gives `N_2 = 6` instead of `−2`.)

**(2) Exact second moment.** Since `ε_F² = 1`,
```
E N_L² = Σ_F 4^(p(F)) = Σ_(λ) 2^p (n)_(L+p) / ∏_ℓ m_ℓ! = Σ_p 2^p (n)_(L+p) C(L/2 − 1, p − 1) / p!,
E κ² = Σ_L E N_L² / ((n)_L)².
```
* Here `λ` runs over the multisets of even positive path lengths with `|λ| = L`, `p` parts, multiplicities `m_ℓ`, and `L + p ≤ n`.
* The count of forests of shape `λ` is `(n)_(L+p)/(2^p ∏ m_ℓ!)`: order the `L + p` vertices, then forget the path directions and the order of equal paths.
* The last form counts compositions of `L/2` into `p` parts.

**(3) Asymptotics.** A shape of `p` parts and size `L` has order `n^(p−L) ≤ n^(−p)`. So only `λ = (2)`, `(2,2)`, `(4)` and `(2,2,2)` contribute through `n^(−3)`:
* `(2)`: `2(n−2)/(n(n−1)) = 2/n − 2/n² − 2/n³ + O(n^(−4))`;
* `(2,2)`: `2(n)_6/((n)_4)² = 2/n² − 6/n³ + O(n^(−4))`;
* `(4)`: `2/n³ + O(n^(−4))`;
* `(2,2,2)`: `(4/3)/n³ + O(n^(−4))`.

Summing,
```
Var κ = 2/n − (14/3) n^(−3) − (58/3) n^(−4) + O(n^(−5)).
```
The `n^(−2)` terms cancel between the one-cherry and two-cherry chaoses. The `n^(−4)` coefficient is the audit's. Exact values of `n·Var κ` are 1.8714 (n = 8), 1.9529 (12), 1.9855 (20), 1.9940 (30), 1.9980 (50).

**(4) The score chaos.** `N_2 = Σ_(u ≠ w) (S²)_(uw) = n(n−1) − W`, with `W = ‖S1‖² = Σ_v (2s_v − n + 1)²` (Session30). So
```
κ_2 = 1 − W/(n(n−1)),    Var κ_2 = 2(n−2)/(n(n−1)),
Var κ_2 / Var κ = 1 − 1/n + (4/3) n^(−2) + O(n^(−3)),
E (κ − 2 + W/(n(n−1)))² = 2/n² + O(n^(−3)).
```

**(5) The classical link and a central limit theorem.** Kendall and Babington Smith's identity for the number `c_3` of cyclic triangles gives `W = n(n² − 1)/3 − 8c_3`. Hence
```
κ_2 = 8c_3/(n(n−1)) − (n−2)/3.
```
The score chaos is an affine function of the cyclic-triangle count. Its variance is the classical `(3/16) C(n,3)` for `c_3`, rescaled; this was noted by the audit from Moon §6, Theorem 6.

Moran's theorem (1947) says `c_3` is asymptotically normal. By (4), `√n(κ − 1 − κ_2) → 0` in `L²`. Therefore
```
√n (κ − 1) ⇒ N(0, 2).
```
This is CONDITIONAL on the cited theorem.

## 2. Reading in the loop expansion

For random tournaments a typical score deviation is `|(S1)_v| ≍ √n`; the maximum is `≍ √(n log n)`. So the path-tree class `T(1) = −δᵀ(I + M²)^(−1)δ` (`δ = S1/n`) is of order one. This is Session30's power counting (section 1 of SCORE_SCALE_AND_RANDOM_CONTROLS.md).

At leading order:
* `T(1) = −W/n² + O(n^(−1))`;
* `log ρ = (n−1)/n + O_(L¹)(n^(−1))`.

So `log ρ + T(1) = 1 − W/n² + …` is the logarithm of `1 + κ_2` to first order. The spectral gain `ρ → e` is paid back by the score energy `W/n² → 1`, tournament by tournament, in `L²` up to `O(n^(−1))`. The average identity `E κ = 1` is the mean of this statement.

Measured on uniform samples with exact `H` (`random_tournaments.out`; the audit's larger samples in brackets, `audit_S54/rand_reading.out`):

| quantity | n = 12 (1500 samples [20,000]) | n = 14 (800 samples [3,000]) |
|---|---|---|
| corr(`log κ`, `log ρ + T(1)`) | 0.957 [0.964] | 0.972 [0.973] |
| corr(`log κ`, `−W`) | 0.955 [0.962] | 0.970 [0.971] |
| corr(`log κ`, `log ρ`) | −0.02 [−0.04] | −0.05 [−0.06] |
| mean of `log(κ/ρ) − T(1)` | +0.034 [+0.049] | +0.045 [+0.048] |

Adding the regular two-loop terms `X(1,0) + E_0[a]/n` changes little. Adding also the lollipop, `X(0,2)` and the star worsens the mean residual to −0.70 and −0.53. The almost-regular grading `w = e + leaves` is the wrong truncation at score scale `√n`.

**Not claimed.** No statement is made for random *regular* tournaments. That is Wormald's Lemma 3 setting, where `κ → e` (THM-4757 prior-art note).
