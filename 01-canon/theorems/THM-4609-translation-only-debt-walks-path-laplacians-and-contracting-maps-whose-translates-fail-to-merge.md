---
id: THM-4609
title: "Translation-only debt walks: if every multiplier is 1 mod d, the coupling of two Haar d-adic orbits is a translation by e mod d and the debt step is uniform over the cyclic roots v_(j+b) - v_j, with covariance (1/d)(2I - A_b) for independent multipliers (A_b = adjacency of the lag-b cycle with the unit positions deleted, a union of paths for prime d); a one-step Lamperti form makes the debt walk transient against ANY sequence of lags; for prime d >= 5 every independent configuration of rank >= 3 has one (the standard form for rank >= 4, two explicit integer forms for the two rank-3 types), so contracting maps such as x/5, (x+4)/5, (6x+3)/5, (11x+2)/5, (16x+1)/5 have translates that fail to merge with positive probability; rank zero merges almost surely iff a finite offset chain reaches 0"
status: >
  PROVED: statements 1-5 (elementary: d-adic Terras property, Lamperti's Lyapunov function |x|^(-alpha), optional stopping,
  path-graph spectra 2 - 2cos(pi m/(n+1)) < 4, exact Sylvester checks of the explicit 3x3 forms). FINITE-EXACT: the two
  rank-3 forms (rank3_universal.py); every AGL(1,d)-type of independent configuration with rank >= 3 for d = 5, 7, 11 has an
  exactly verified form (balanced_types.py: 2 + 6 + 24 types). NUMERICAL: for the Z_5 example, 12.8-13.5% of 2000-digit integers
  n merge with n+1 within 100-4000 steps (flat), while every orbit from n <= 20000 reaches one of four cycles; an adversarial
  lag chooser does not stop the escape (median |x|_A = 24 after 4000 steps). Not yet independently audited.
session: mac-mini-2026-10-08-rank
source: 05-knowledge/results/debt_rank_least_requirements_20261008.md
scripts:
  - 04-computation/experiments/rank_20261008/rank3_universal.py, balanced_types.py, balanced_lamperti.py, balanced_whitened.py (+ .out)
  - 04-computation/experiments/rank_20261008/z5_types_exact.py, adversarial_lag.py, z5_example_integers.py (+ .out)
related:
  - HYP-9244 (the trichotomy; this proves its rank >= 3 row for the translation-only class, without any mixing hypothesis)
  - THM-4607 (expanding maps), THM-4608 (Z_2), THM-4581 (rank one, d = 2)
  - Lamperti (1960); Menshikov, Popov and Wade, Non-homogeneous Random Walks (2016); Peres, Popov and Sousi (2013) (adaptive laws)
---

# THM-4609 — translation-only debt walks

## Setting

* A Matthews–Watts map `T(x) = (m_i x + r_i)/d` on `x ≡ i mod d` is **translation-only** if every `m_i ≡ 1 (mod d)`.
* Pair chain `u = M v + e` (HYP-9244). Let `v_j` be the log-vector of `m_j` in `Γ ⊗ R`, where `Γ = ⟨m_a/m_b⟩` has rank `ρ`.
* **Lag covariance:** `C_b = (1/d) Σ_j (v_(j+b) − v_j)(v_(j+b) − v_j)^T` for `b ∈ Z/d`.
* **Independent configuration:** the multipliers equal 1 at a set `Z` of positions, and the remaining `ρ = d − |Z|` multipliers are multiplicatively independent. Only differences `v_(j+b) − v_j` matter, so any single repeated value can play the role of 1.
* **Balanced:** some positive definite `Q` has

      S_b := ½ tr(C_b Q^(−1)) Q − C_b  ≻ 0   for every b with C_b ≠ 0.

## Statements

1. **Translation coupling (PROVED).**
   * For a translation-only map, `M ≡ 1 (mod d)` always, so `u`'s digit is `i = j + ē`, where `ē = e mod d` is past-measurable and `j` is fresh and uniform.
   * So given the past, the debt step `v_(j+ē) − v_j` is uniform over the `d` cyclic roots of lag `ē`. Its covariance is `C_ē`; at `ē = 0` the step is 0.
   * For independent configurations, in the basis of the independent vectors, `C_b = (1/d)(2I − A_b)`. Here `A_b` is the adjacency matrix of the graph on the non-unit positions, with `p ~ p'` iff `p' − p ≡ ±b`.
   * For prime `d` and `b ≠ 0`, this graph is the `d`-cycle `j → j + b` with the unit positions deleted, a disjoint union of paths. So `C_b` has spectrum `(1/d)(2 − 2cos(πm/(n+1)))` over paths of `n` vertices, all eigenvalues are `< 4/d`, and `tr C_b = 2ρ/d`.
2. **Rank zero (PROVED; any d).**
   * If all `m_i` equal `m`, then `M ≡ 1`. The offset is an integer chain `e' = (m e + r_(j+e) − r_j)/d` driven by the fresh digit `j`.
   * If `m < d` (contracting), `|e|` enters `{|e| ≤ 2R/(d − m)}` and the chain is a finite Markov chain.
   * `y` and `y + e_0` merge almost surely iff `0` is reachable from every state reachable from `e_0`, with a geometric tail. Otherwise they merge with probability `< 1`.
   * The congruence obstruction of HYP-9244 produces closed classes avoiding 0.
3. **Lamperti transience (PROVED).**
   * If a translation-only map is balanced, choose `0 < α < min_b tr(C_b Q^(−1))/λ_max(C_b Q^(−1)) − 2` and put `V(x) = (x^T Q^(−1) x)^(−α/2)` on the debt lattice.
   * There is `r_1` such that, while `|x| ≥ r_1`, `E[V(x_(t+1)) − V(x_t) | F_t] ≤ 0` whatever the lag `ē_t` is.
   * Hence, from debt `x` with `|x|_(Q^(−1)) = r > r_1 + B` (B the step bound), the debt walk ever comes within `r_1` of the origin with probability at most `(r_1/(r − B))^α`.
   * Since absorption needs debt 0, every non-absorbed start, in particular `(1, e_0)` with `e_0 ≠ 0`, is absorbed with probability `< 1`. No hypothesis on how the lag `ē_t` evolves is used: an adversary choosing the lag at every step cannot prevent the escape.
4. **The balance criterion for prime `d ≥ 5` (PROVED + FINITE-EXACT).** Every independent configuration of rank `ρ ≥ 3` is balanced.
   * **`ρ ≥ 4`.** The standard form `Q = I` works, since `λ_max(C_b) < 4/d ≤ ρ/d = tr C_b/2`.
   * **`ρ = 3`.** By statement 1, the lag covariances are scale-free members of one of two universal families:
     * AP type (the three positions form an arithmetic progression mod d; for `d = 5` every 3-set does): `{2I − A(P_3), 2I − E_13, 2I}`. Here `Q = I` gives equality. The form
       `Q_AP = [[6, −2, 0], [−2, 6, −1], [0, −1, 5]]` (basis: end, middle, end) is balanced.
     * Non-AP type: `{2I − E_12, 2I − E_23, 2I − E_13, 2I}`. Again `Q = I` gives equality. The form
       `Q_nonAP = [[4, −1, −1], [−1, 5, 0], [−1, 0, 5]]` is balanced.
     * Both forms are verified exactly by Sylvester's criterion (`rank3_universal.py`). A form valid for the larger family covers the cases where `2I` does not occur (`d = 5`, `d = 7`).
5. **Contracting maps whose translates fail to merge (PROVED).**
   * `Z_5`, multipliers `(1, 1, 6, 11, 16)`: `T(x) = x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5`. `Λ = −0.217`, rank 3 (6, 11, 16 independent; AP type).
   * `Z_7`, multipliers `(1, 1, 1, 1, 8, 15, 22)` (rank 3, `Λ = −0.820`) and `(1, 1, 1, 8, 15, 22, 29)` (rank 4, `Λ = −0.339`). For both, `r_i = (0, 6, 5, 4, 3, 2, 1)`.
   * In each, Haar `y` and `y + e` merge at equal time with probability `< 1`, although the map contracts. The same holds for every contracting translation-only map of independent rank `≥ 3` on a prime base `d ≥ 5`.
   * NUMERICAL, `Z_5` example:
     * only 12.8–13.5% of random 2000-digit integers `n` merge with `n + 1` within 100–4000 steps (flat);
     * every orbit from `n ≤ 20000` reaches one of four cycles (minimal elements 1, 8, 9, 33);
     * 28.6% of consecutive pairs `n ≤ 20000` merge at equal time.

## Proofs

**1.**
* `M` is a product of ratios of multipliers `≡ 1`, so `M ≡ 1 (mod d)`. The coupling `i ≡ M j + e` is then a translation.
* `j` is fresh and uniform (d-adic Terras), and `ē` is `F_t`-measurable.
* **Covariance.** Each non-unit position `p` lies in exactly two lag-`b` roots, `v_p − v_(p−b)` and `v_(p+b) − v_p`, with coefficients `±1` on its basis vector. Hence the diagonal is `2/d`.
* An off-diagonal entry is `−1/d` for each root `±(e_p' − e_p)`, i.e. each `p' − p ≡ ±b`.
* For prime `d` and `b ≠ 0`, `j ↦ j + b` is one `d`-cycle. Deleting `|Z| ≥ 1` vertices leaves paths. If `Z` is empty, shift all vectors by `−v_0` to get `|Z| = 1`.
* Path adjacency spectra are `2cos(πm/(n+1)) > −2`. ∎

**2.**
* `e' = u' − v'` is an integer: `m e + r_i − r_j ≡ m(e − i + j) ≡ 0 (mod d)`, using `r_i ≡ −m i`.
* `|e'| ≤ (m|e| + 2R)/d` with `m < d`.
* Finite-state Markov chain theory gives the rest. ∎

**3.**
* In the coordinates `y = Q^(−1/2) x`, the step `η` has conditional mean 0, covariance `Σ = Q^(−1/2) C_ē Q^(−1/2)` and is bounded.
* Taylor expansion gives

      E[ΔV] = −(α/2)|y|^(−α−2)(tr Σ − (α + 2) ŷ^T Σ ŷ) + O(|y|^(−α−3))
            ≤ −(α/2)|y|^(−α−2)(tr Σ − (α + 2)λ_max(Σ)) + O(|y|^(−α−3)).

* Balance is equivalent to `tr Σ > 2λ_max Σ` for each `Σ ≠ 0`: `λ_max(Σ) = max_w w^T C w/w^T Q w`, and `S_b ≻ 0` says this is `< ½ tr`. So the bracket is positive for small `α`, and `E[ΔV] ≤ 0` for `|y| ≥ r_1`. When `ē = 0` the step is 0.
* Optional stopping for the bounded supermartingale `V(x_(t∧τ))` gives the return bound.
* **Escape with positive probability, from any non-absorbed state.**
  * At debt 0 with `ē = 0`, `e' = m_j e/d` lowers `v_d(e)`, so the walk moves within `v_d(e)` steps.
  * Elsewhere, staying unmoved forever along every continuation would make a positive-measure cylinder merge with `M ≠ 1`, which is null.
  * At each moving step some root of the current lag has `⟨ξ, x⟩_(Q^(−1)) ≥ 0`. The roots of a lag span the lattice and have mean 0, and for `x = 0` any nonzero root will do.
  * Each such root raises `|x|²` by at least `min |ξ|² > 0`, and the debt never revisits 0 along the way.
  * So a finite digit word reaches `|x| ≥ r` without absorption. Then apply the return bound. ∎

**4.**
* `ρ ≥ 4`: `λ_max(2I − A_b) < 4 ≤ ρ = ½ tr(2I − A_b)`.
* `ρ = 3`: the lag graphs on three positions are as listed. A 2-edge path occurs exactly at the lags `±c` of an arithmetic progression with difference `c`, and its end-to-end edge at `±2c`. For a non-AP set the three pair differences are distinct up to sign.
* The exact checks are in `rank3_universal.py`. ∎

**5.**
* `6 = 2·3`, `11` and `16 = 2^4` are independent. So are `8 = 2^3`, `15 = 3·5`, `22 = 2·11` and `29`.
* `Λ < 0` is direct: `Π m_i = 1056 < 5^5`, `2640 < 7^7` and `76560 < 7^7`.
* Apply 3 and 4. ∎

## Reading

* **The least-requirement principle.** For translation-only maps every possible one-step law of the debt is a cyclic-root law. So one fixed quadratic form controls all of them at once: no mixing, equidistribution or decorrelation of the coupling state is needed. This is why rank `≥ 3` transience is provable here, while for general maps (couplings `j ↦ aj + b`, including involutions with rank-2 covariances) it is still a conjecture.
* **Path Laplacians.** The debt covariance is the Laplacian-like matrix `2I − A` of a path system cut out of a cycle. A path's spectrum never reaches 4, so from rank 4 on the walk escapes by pure dimension counting. Rank 3 needs a tilted form because the path `P_3` and the single edge saturate the standard inequality.
* **What this does to the conjecture.** The rank `≥ 3` row of HYP-9244 holds for the whole translation-only class, with no hypothesis on the coupling dynamics. Contraction does not imply coalescence once the digit base allows three independent multiplier directions. On `Z_2` it does (THM-4608).
