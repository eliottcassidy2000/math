---
id: HYP-9365
title: "Two-loop law for Hamiltonian paths in regular tournaments: H(T) 2^(n-1)/n! = rho(T) (1 - y(T)/n + O(n^-2)) with y(T) = -n X_1(T) - E_0[C(L,2)](T) explicit (THM-4757). Consequences: (a) every doubly regular tournament (Paley included) has y -> 1/2, so H(T) 2^(n-1)/n! = e (1 - 3/(2n) + O(n^-2)); (b) rotational R_n: y -> tan(1) tanh(1) + tan(1)^2/2 - 1 = 1.398872...; (c) ternary towers: y -> 1.345726...; (d) all doubly regular tournaments of one order agree through three loops, H(D)/H(D') = 1 + O(n^-3); (e) the n mod 4 wobble of H(R_n) is non-perturbative and tracks the alternating count F(R_n,-1)/n!."
status: >
  OPEN (mac-mini-2026-10-22-sierp; independently audited 2026-10-22, 05-knowledge/results/sierpinski_audit2_20261022.md:
  every exact H, y and y_pred recomputed and confirmed; wording and credits corrected). The expansion is derived (THM-4757: exact forest identity, proved one- and
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

*Notation below:* `y` is the exact residual `n(1 − κ/ρ)` and `y_pred` the two-loop value `y(T)` defined above; the conjecture says `y = y_pred + O(1/n)`.

**Derivation.** THM-4757 (1) gives `κ = Σ_L φ_L n^L/(n)_L`, where `Φ_T(z) = Σ_L φ_L z^L = exp(X_0 + X_1 + X_2 + …)`.
* `n^L/(n)_L = 1 + C(L,2)/n + O(L⁴/n²)`.
* `X_0(1) = log ρ`, and `X_1 = O(1/n)`.
* At order `1/n`, the mean of `C(L,2)` may be taken under the one-loop measure. There `log Φ_0(z) = ½ Σ log((1+z²x²)/(1−z²x²))`, so its `z∂_z`-cumulants are `D_1, D_2`.

## Consequences

**(a) Doubly regular tournaments, Paley included.** `S² = J − nI`, so every `x_k² = 1/n` for `k ≠ 0`, and `G, H, O` are explicit elements of `span{I, J, Ŝ}`.
* Then `n·F8 → −7/2`, the theta and dumbbell terms are `O(1/n)`, and `E_0[C(L,2)] → 3`, so `y → 1/2`.
* Since `Σ_k x_k² = (n−1)/n` for every tournament, and here `Σ_k x_k^(2m) = (n−1)/n^m`, `ρ = ((n+1)/(n−1))^((n−1)/2) = e(1 − 1/n + O(n^(−2)))`. Hence

      H(T) 2^(n−1)/n! = e (1 − 3/(2n) + O(n^(−2))).

* This sharpens HYP-9274 (A) with an explicit rate. The constant `e` itself is the classical random-regular value: a uniformly random regular tournament has `E[H] ∼ e·n!/2^(n−1)` (Wormald, arXiv:2610.16245: the case `d = 1` of Theorems 1 and 5, via Lemma 3 and McKay's method; Alon's observation that the Adler–Alon–Ross bound holds for regular tournaments, recorded by Yuster, arXiv:1511.07784). New here are the single-tournament statement and the rate `3/(2n)`. The `1/n` coefficient is checked by an exact brute-force cluster sum: figure-eight `−7/2`, theta `+1·O(1/q²)` (`cluster_shapes.py`).
* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`. Polynomial extrapolation in `1/q` gives 0.502, 0.479, 0.425 from the last 2, 3, 4 points, so it fixes the limit only loosely. For every doubly regular `T` the two-loop value is exactly `y_pred = (3n² − 8n + 1)/(6(n − 1)²) = 1/2 − 1/(3n) − 1/n² + O(n^(−3))`, and the residual `q(y − y_pred) = −1.55, −1.53, −1.50, −1.49` is stable, as an `O(1/q)` three-loop term should be; `κ(P_q)/(e(1 − 3/(2q))) − 1 = 0.0094, 0.0063, 0.0045, 0.0034`, about `3.3/q²` (`audit2/f_families.out`; `H(P_27)` recomputed in `audit2/d_paley27.out`).

**(b) Rotational tournaments.** The large `x_k` tend to `±2/(πj)`, `j` odd.
* So `tr G → tan 1` and `tr H → −tanh 1`, giving `n·F8 → −2 tan 1 tanh 1 − (3/2)tan² 1`.
* `E_0[C(L,2)] → 1 + tan² 1 + tan 1 tanh 1`, and `n·Θ → 0` (like `1.55/n`).
* Hence

      y_∞(R) = tan 1 · tanh 1 + ½ tan² 1 − 1 = 1.398872…

* Exact data: `y(R_33) = 1.2106`. Mixed polynomial fits in `1/n` are unstable (1.381, 1.3975, 1.406, 1.475, 1.968 from the last 2–6 points) because of the `n mod 4` term of (e); fitted within each residue class mod 4, the 4- to 7-point fits give 1.3971–1.4002 (`audit2/f_families.out`).
* Since `ρ(R_n) = (cosh 1/cos 1)(1 − 1/n + O(n^(−2)))` (Bhattacharjee's expansion of the regular-block gain `e^(1/n) ρ(R_n)` has no `1/n` term), the conjecture predicts `κ(R_n) = (cosh 1/cos 1)(1 − (1 + y_∞)/n + O(n^(−2)))`, with `1 + y_∞ = 2.398872`. This is conditional on HYP-9365; it does not yet settle the approach to Wormald's constant.

**(c) Ternary towers `S_k = C_3^[k]`.** Here the theta term survives, `n·Θ → −0.237054`, because the outlier frequencies `3^j·u` are additively closed.
* The two-loop limit is `y_∞ = 1.345726` (`two_loop_fft.out`, `k ≤ 9`).
* Richardson on the exact data at `n = 243, 729` gives 1.345606.

**(d) Universality of doubly regular tournaments.** For every doubly regular `T` of order `n`, all propagators lie in the 3-dimensional algebra `span{I, J, Ŝ}`.
* Every three-loop core has at most four vertices.
* So every two- and three-loop term is a signed count of subtournaments on at most 4 vertices. For doubly regular tournaments these counts are fixed by `n`, because out-neighbourhoods are regular of degree `(n−3)/4`.
* Hence `X_0`, `X_1`, `X_2` are the same formal power series for all doubly regular tournaments of one order (a theorem, given THM-4757 (2)). That non-isomorphic ones then agree to relative `O(n^(−3))` in `H` is the conjectural part.
* **Test (order 19).** `H(P_19) = 1172695746915`. A non-Paley, self-converse doubly regular 19-tournament `D_19` (simulated annealing on `SSᵀ = nI − J`, nauty-certified non-isomorphic, `drt_anneal.c`) has `H(D_19) = 1172732961771`, a ratio of `1 + 3.17·10^(−5)`. Regular 19-tournaments spread over several percent.
* **Test (order 23).** Three further non-Paley classes (nauty-distinct) give `(H − H(P_23))/H(P_23) = −1.61·10^(−5), −1.94·10^(−5), −1.93·10^(−5)`, with `H(P_23) = 15760206976379349` (`drt_results.out`).
* The scaled gaps `n³·(H/H(P) − 1)` are `+0.218` at order 19 and `−0.196, −0.236, −0.235` at order 23. They have the same size at both orders, which is consistent with an `O(n^(−3))` (four-loop) difference; two orders, with a sign change between them, cannot separate `n^(−3)` from `n^(−2)` or `n^(−4)`. At order 23, `D_23.2` and `D_23.3` have identical 5-vertex censuses (their 6-vertex censuses differ), and their `H` differ by only `2.8·10^(−8)` relative (`audit2/e_drt_noIE.out`, `audit2/e_drt_census6.out`).
* At order 15 the annealer found one doubly regular tournament up to isomorphism and converse (the two nauty classes are converses of each other, so their equal `H` is trivial).
* `D_19` has more Hamiltonian paths than `P_19`, though fewer than `R_19`.
* The census agrees with the mechanism (`drt_subcounts.out`). `P_19` and `D_19` have identical 4-vertex subtournament counts (1026, 570, 570, 1710), and their 5-vertex counts differ (for example 342 against 360 transitive 5-sets). Four-loop cores have up to 6 vertices, so that is where the difference first enters.

**(e) A non-perturbative term for `R_n`.** The values `κ(R_n)` split by `n mod 4`.
* Fitted separately, `n ≡ 1 (mod 4)` (7 points) gives `κ → 2.855930` and `n ≡ 3 (mod 4)` gives 2.85611–2.85631. The target is `cosh 1/cos 1 = 2.855958`. A joint polynomial fit in `1/n` is unstable.
* The alternating component `w_n` has sign `(−1)^((n+1)/2)` and decays geometrically, at ratio ≈ 0.6 per vertex (`rot_limit_fit.out`).
* It tracks the top term of the forest series, `N_(n−1)/n! = F(R_n, −1)/n!`, the alternating forward-step count (THM-4756). That count has the same sign pattern (+ at `n = 7, 11, 15, 19`; − at 9, 13, 17, 21), and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.69` for `n = 13–21` with the cubic smooth fit of `rot_limit_fit.py` (`fneg1.out`). The ratio depends on the fit: 1.60–1.70 with a quartic, 0.40–2.76 with a quadratic (`audit2/g_wobble.out`). The values and signs of `F(R_n, −1)` are FINITE-EXACT (`n ≤ 21`); the ratio is NUMEROLOGY.
* The `1/n` loop expansion cannot produce such a term. It comes from forests of length `L ≈ n`, and is exponentially small.

## Evidence summary

| family | exact `y` at largest `n` | two-loop `y_pred` (same `n`) | two-loop limit | exact extrapolation |
|---|---|---|---|---|
| `P_q`, `q ≤ 31` | 0.4400 (31) | 0.4882 | 1/2 | 0.43–0.50 |
| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3971–1.4002 (within residue classes mod 4) |
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

## Addendum (same session, after audit 2)

**More populations, no fitted parameters.**
* Random regular tournaments from a 3-cycle-reversal chain: `n = 15` (200 samples) corr 0.9926, Spearman 0.9898; `n = 21` (150) corr 0.9975, scatter 0.0028 (`two_loop_random_regular_*.out`).
* **Circulant maximisers.** Ranking circulant classes by `ρ(1 − y_pred/n)` picks the true `H`-maximiser at every `n = 7, …, 25`: `P_7`, `{2,3,4,8}`, `P_11`, then `R_n` (`two_loop_argmax.out`). This explains THM-4751's flat-to-rotational crossover. From `n = 19` the predicted winning margins match the true ones to `10^(−4)`.
* **All regular 13-tournaments** (1,495,297, exact `H`). The maximum `H = 3719893` (a converse pair, as in HYP-9357) ranks 7th under the two-loop predictor and 1,476,889th under `ρ` alone (`reg13_argmax.out`).

**The next order is computed** (HYP-9367).
* A general loop engine (15 three-loop cores) reproduces the brute-force excess-2 cluster sums exactly.
* With it, the second-order prediction leaves residuals `O(n^(−3))`: `n³·resid ≈ 6` on Paley, 3.2–3.7 on `R_11..R_35`. This resolves the drift of `n(y − y_pred)` listed under the refutation targets.
* For doubly regular tournaments it gives the universal curve `κ = e(1 − 3/(2n) + 25/(8n²) + O(n^(−3)))`, with `q² X_2 → 29`. Exact data: 3.09–3.13 for the `n^(−2)` coefficient.
* At order 23 the non-universal remainder is visibly a four-loop effect. `H` is linear in the transitive-5-set count across 12 non-Paley classes, to `3.7·10^(−7)`.

**Rotational tournaments.** `H(R_35) = 1610822294345031078388661430315` (`R35.out`) gives:
* `n(y − y_pred) = −1.046` (the second order accounts for this);
* residue-class fits for `n ≡ 3 (mod 4)` within `2·10^(−3)` of `1.398872` (audit 2, `f_R35_check.out`);
* `κ(R_n) → 2.856045` from the `n ≡ 3` class (`rot_limit_fit.out`).

**The cosine twin** (HYP-9366). `R_n` and the transitive tournament have the same skew spectrum. The alternating end `F(R_n,−1)/n!` decays at rate `2/π`.
