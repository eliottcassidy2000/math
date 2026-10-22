---
id: HYP-9365
title: "Two-loop law for Hamiltonian paths in regular tournaments: H(T) 2^(n-1)/n! = rho(T) (1 - y(T)/n + O(n^-2)) with y(T) = -n X_1(T) - E_0[C(L,2)](T) explicit (THM-4757). Consequences: (a) every doubly regular tournament (Paley included) has y -> 1/2, so H(T) 2^(n-1)/n! = e (1 - 3/(2n) + O(n^-2)); (b) rotational R_n: y -> tan(1) tanh(1) + tan(1)^2/2 - 1 = 1.398872...; (c) ternary towers: y -> 1.345726...; (d) all doubly regular tournaments of one order agree through three loops, H(D)/H(D') = 1 + O(n^-3); (e) the n mod 4 wobble of H(R_n) is non-perturbative and tracks the alternating count F(R_n,-1)/n!."
status: >
  OPEN, pre-audit (mac-mini-2026-10-22-sierp). The expansion is derived (THM-4757: exact forest identity, proved one- and
  two-loop cluster sums); its asymptotic validity (that excess-e clusters are O(n^-e) uniformly and the expansion is
  asymptotic) is the open part. Evidence (FINITE-EXACT): all 208 circulant classes at n = 25, corr(y, y_pred) = 0.9984,
  scatter 0.009 against a range 0.46-1.16, offset -0.037 = O(1/n); all 1223 regular 11-tournaments, corr 0.978, and
  the parameter-free corrected predictor ranks kappa with Spearman 0.896 (rho alone -0.743) and finds P_11; towers to
  729 vertices, y = 1.338620 against 1.339738; R_n to n = 33; P_q to q = 31; doubly regular 19-tournaments
  H(D_19)/H(P_19) = 1 + 3.2e-5. Supersedes the empirical second-order fit a + b M4 of HYP-9364 (f).
session: mac-mini-2026-10-22-sierp
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/two_loop_predictor.py (+ .out): circulants, R_n, P_q
  - 04-computation/experiments/sierpinski_ocf_20261022/two_loop_regular.py (+ .out): all regular n = 9, 11
  - 04-computation/experiments/sierpinski_ocf_20261022/two_loop_towers.py, two_loop_fft.py (+ .out): towers, limits
  - 04-computation/experiments/sierpinski_ocf_20261022/drt_test.py, drt_anneal.c, drt_switch.py, hk_big.c, hk64.c: doubly regular tournaments
  - 04-computation/experiments/sierpinski_ocf_20261022/rot_limit_fit.py (+ .out), fneg1.c (+ .out), paley_y.py (+ .out)
related:
  - 05-knowledge/hypotheses/HYP-9364-self-blow-up-principle-tower-limit-rotational-cosh-over-cos.md
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
  - 05-knowledge/hypotheses/HYP-9274-paley-rate-three-halves-wheel-w4-apery-form.md
---

# HYP-9365 — the two-loop law

## Statement

Let `T` be regular on `n` vertices, `Ŝ = S/n` with eigenvalues `i x_k`, and `G, H, O, d` as in THM-4757. Put

    X_1 = Σ_v [2 G_vv H_vv − (3/2) G_vv²] + Σ_(u,v) G_uv (H_uv² − O_uv² − (2/3) G_uv²) − dᵀ G d,
    D_1 = Σ_k 2x_k²/(1 − x_k⁴),   D_2 = Σ_k [2x_k²/(1 − x_k²)² + 2x_k²/(1 + x_k²)²],
    y(T) = −n X_1 − (D_2 + D_1² − D_1)/2.

**Conjecture.** `κ(T) = H(T) 2^(n−1)/n! = ρ(T) (1 − y(T)/n + O(n^(−2)))`, uniformly over regular tournaments.

**Derivation.** THM-4757 (1) gives `κ = Σ_L φ_L n^L/(n)_L`, where `Φ_T(z) = Σ_L φ_L z^L = exp(X_0 + X_1 + X_2 + …)`.
* `n^L/(n)_L = 1 + C(L,2)/n + O(L⁴/n²)`.
* `X_0(1) = log ρ`, and `X_1 = O(1/n)`.
* At order `1/n`, the mean of `C(L,2)` may be taken under the one-loop measure. There `log Φ_0(z) = ½ Σ log((1+z²x²)/(1−z²x²))`, so its `z∂_z`-cumulants are `D_1, D_2`.

## Consequences

**(a) Doubly regular tournaments, Paley included.** `S² = J − nI`, so every `x_k² = 1/n` for `k ≠ 0`, and `G, H, O` are explicit elements of `span{I, J, Ŝ}`.
* Then `n·F8 → −7/2`, the theta and dumbbell terms are `O(1/n)`, and `E_0[C(L,2)] → 3`, so `y → 1/2`.
* Since `Σ_k x_k² = (n−1)/n` for every regular tournament, `ρ = e(1 − 1/n + O(n^(−2)))`. Hence

      H(T) 2^(n−1)/n! = e (1 − 3/(2n) + O(n^(−2))).

* This sharpens HYP-9274 (A) with an explicit rate. The `1/n` coefficient is checked by an exact brute-force cluster sum: figure-eight `−7/2`, theta `+1·O(1/q²)` (`cluster_shapes.py`).
* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`; the two-point extrapolation gives 0.502.

**(b) Rotational tournaments.** The large `x_k` tend to `±2/(πj)`, `j` odd.
* So `tr G → tan 1` and `tr H → −tanh 1`, giving `n·F8 → −2 tan 1 tanh 1 − (3/2)tan² 1`.
* `E_0[C(L,2)] → 1 + tan² 1 + tan 1 tanh 1`, and `n·Θ → 0` (like `1.55/n`).
* Hence

      y_∞(R) = tan 1 · tanh 1 + ½ tan² 1 − 1 = 1.398872…

* Exact data: `y(R_33) = 1.2106`, with extrapolations 1.3975 (3 points) and 1.406 (4 points).
* Together with `ρ(R_n) → cosh 1/cos 1`, this settles the shape of the approach to Wormald's constant: `κ(R_n) = (cosh 1/cos 1)(1 − c/n + …)`.

**(c) Ternary towers `S_k = C_3^[k]`.** Here the theta term survives, `n·Θ → −0.23703`, because the outlier frequencies `3^j·u` are additively closed.
* The two-loop limit is `y_∞ = 1.345726` (`two_loop_fft.out`, `k ≤ 9`).
* Richardson on the exact data at `n = 243, 729` gives 1.345606.

**(d) Universality of doubly regular tournaments.** For every doubly regular `T` of order `n`, all propagators lie in the 3-dimensional algebra `span{I, J, Ŝ}`.
* Every three-loop core has at most four vertices.
* So every two- and three-loop term is a signed count of subtournaments on at most 4 vertices. For doubly regular tournaments these counts are fixed by `n`, because out-neighbourhoods are regular of degree `(n−3)/4`.
* Hence non-isomorphic doubly regular tournaments of one order should agree to relative `O(n^(−3))`.
* **Test (order 19).** `H(P_19) = 1172695746915`. A non-Paley, self-converse doubly regular 19-tournament `D_19` (simulated annealing on `SSᵀ = nI − J`, nauty-certified non-isomorphic, `drt_anneal.c`) has `H(D_19) = 1172732961771`, a ratio of `1 + 3.17·10^(−5)`. Regular 19-tournaments spread over several percent.
* **Test (order 23).** Three further non-Paley classes (nauty-distinct) give `(H − H(P_23))/H(P_23) = −1.61·10^(−5), −1.94·10^(−5), −1.93·10^(−5)`, with `H(P_23) = 15760206976379349` (`drt_results.out`).
* The scaled gaps `n³·(H/H(P) − 1)` are `+0.218` at order 19 and `−0.196, −0.236, −0.235` at order 23. They have the same size at both orders, as an `O(n^(−3))` (four-loop) difference requires.
* At order 15 the annealer found one doubly regular tournament up to isomorphism and converse (the two nauty classes are converses of each other, so their equal `H` is trivial).
* `D_19` has more Hamiltonian paths than `P_19`, though fewer than `R_19`.
* The census agrees with the mechanism (`drt_subcounts.out`). `P_19` and `D_19` have identical 4-vertex subtournament counts (1026, 570, 570, 1710), and their 5-vertex counts differ (for example 342 against 360 transitive 5-sets). Four-loop cores have up to 6 vertices, so that is where the difference first enters.

**(e) A non-perturbative term for `R_n`.** The values `κ(R_n)` split by `n mod 4`.
* Fitted separately, `n ≡ 1 (mod 4)` (7 points) gives `κ → 2.855930` and `n ≡ 3 (mod 4)` gives 2.85611–2.85631. The target is `cosh 1/cos 1 = 2.855958`. A joint polynomial fit in `1/n` is unstable.
* The alternating component `w_n` has sign `(−1)^((n+1)/2)` and decays geometrically, at ratio ≈ 0.6 per vertex (`rot_limit_fit.out`).
* It tracks the top term of the forest series, `N_(n−1)/n! = F(R_n, −1)/n!`, the alternating forward-step count (THM-4756). That count has the same sign pattern (+ at `n = 7, 11, 15, 19`; − at 9, 13, 17, 21), and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.72` for `n = 13–21` (`fneg1.out`).
* The `1/n` loop expansion cannot produce such a term. It comes from forests of length `L ≈ n`, and is exponentially small.

## Evidence summary

| family | exact `y` at largest `n` | two-loop `y_pred` (same `n`) | two-loop limit | exact extrapolation |
|---|---|---|---|---|
| `P_q`, `q ≤ 31` | 0.4400 (31) | 0.4882 | 1/2 | 0.43–0.50 |
| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3975 / 1.406 |
| `S_k`, `n ≤ 729` | 1.338620 (729) | 1.339738 | 1.345726 | 1.345606 |

| population | classes | corr(`y`, `y_pred`) | Spearman | scatter of `y − y_pred` | mean offset |
|---|---|---|---|---|---|
| circulants `n = 25` | 208 | 0.9984 | 0.9974 | 0.0090 | −0.037 |
| circulants `n = 21` | 88 | 0.9969 | 0.9920 | 0.0142 | −0.041 |
| all regular `n = 11` | 1223 | 0.9777 | 0.9731 | 0.0227 | −0.061 |
| all regular `n = 9` | 15 | 0.9700 | 0.9643 | 0.0718 | −0.053 |

**Refutation targets.**
* A family of regular tournaments in which `n·(y − y_pred)` grows without bound would refute the conjecture.
* So would two doubly regular tournaments of order `n ≥ 23` whose `H` differ by `≫ n^(−3)` relative.
* For the rotational and circulant families, `n·(y − y_pred)` (the three-loop coefficient) drifts slowly: −0.86 to −1.02 for `R_25..R_33`, and −0.82 for the tower at 729.
