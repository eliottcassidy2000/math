---
id: THM-4712
title: "Mean cluster size on the Mersenne line: the expected number of deletion children merged by Terras time T is the linearity sum sum_D P(tau_D <= T), bounded above by sqrt(pi T/2) (Azuma on THM-4581's level martingale) and below by c sqrt(T/log T) (THM-4581 (4)); the SRW skeleton predicts sqrt(T/pi) = 0.564 sqrt(T); the exact horizon to K = 50000 (reproducing K <= 12800 byte for byte) has 231 clusters = level sets of (sigma_T - K, odd count) = distinct odd counts, orphan fraction x sqrt(K) in [0.76, 0.99] for K >= 5000, members below K / sqrt(sigma_T(x_K)) in [0.41, 0.62], least shift 1 in 95.5% of cases with the D^-2 law for D >= 3, interleaving depth <= 3, exactly three terminal channels 5/85/341, and payment (pay - K)/K = 2.8213 +- 0.07 against the drift value 2.8188"
status: >
  PROVED: statements 2 and 3 (linearity; the upper bound, from the martingale property of the level process proved in
  THM-4581 3(i) and the Azuma-Hoeffding maximal inequality). PROVED at the sketch level of THM-4581 (4): statement 4 (the
  lower bound; it inherits the |k_0| first-passage term of that proof). INHERITED/PROVED: statement 1 (clusters are the
  level sets; mersenne_line_barriers_20261007 Prop. 4.1 and HYP-9242's partner criterion). HEURISTIC + NUMERICAL:
  statement 5 (the SRW-skeleton constant 1/sqrt(pi) against the data). FINITE-EXACT: statement 6 (exact Terras orbits of
  every 2^K - 1, 2 <= K <= 50000, C + GMP, reproducing the previous K <= 12800 table exactly). Not claimed: any grounding
  of a Mersenne exponent (deletion certificates preserve the orbit length: barriers Prop. 4.1, codex's finite-anchor
  closure bound), the orphan exponent beyond its numerical range, or any statement about integers outside the Haar model
  for statements 2-5. Pre-audit.
session: mac-mini-2026-10-19-mersenne
source: 05-knowledge/results/mersenne_horizon_20261019_clusters_orphans_reanchoring.md
scripts:
  - 04-computation/experiments/mersenne_horizon_20261019/horizon.c (exact orbits, 6 threads; data mersenne_horizon_50000.txt.gz: K odd terras last_j peak pay)
  - 04-computation/experiments/mersenne_horizon_20261019/analyze_horizon.py (+ .out, summary .json), analyze_horizon_extra.py (+ .out)
related:
  - THM-4581 (Haar coalescence; the pair chain, its level martingale, the rate (4), the one-big-jump heuristic (7), the Mersenne chain 6(d))
  - THM-4556 (the Mersenne line: repunit transit, reset pairs, odd least shift, 2-adic periodicity, census to 6000)
  - THM-4555 (uniform switches are collisions at -1), THM-4640 (the repair algebra at -1), THM-4601 (two anchors, the +1 barrier)
  - 05-knowledge/results/mersenne_line_barriers_20261007.md (certificate engine, stopping-time invariant Prop. 4.1, horizon clusters, HYP-9243)
  - HYP-9242 (orphan law), HYP-9243 (diffusive coalescence), HYP-9332 (this session's horizon laws), HYP-9333 (terminal channels vs THM-4517)
  - THM-4517 (trunk-entry basin densities e_2 = 0.93796, e_4 = 0.0236, e_5 = 0.0378)
  - 05-knowledge/results/mersenne_joint_guards_20261014c.md, mersenne_clocked_completion_20261014c.md (codex: short paired repairs, fixed-head classifier, finite-anchor closure bound)
---

# THM-4712 — mean cluster size on the Mersenne line, and the exact horizon to K = 50000

## Setting

* `M_K = 2^K − 1`, `T` is the Terras map, `U` the odd map, `x_K = U^(K−1)(M_K) = 2·3^(K−1) − 1` the post-run state.
* `σ_T(n)` = Terras steps from `n` to 1; `odd(n)` = odd steps from `n` to 1. The **deletion invariant** of `K` is `(σ_T(M_K) − K, odd(M_K))`.
* A **deletion certificate** `M_K ⇝ M_(K−D)` is an equal-`U`-time merge of the two orbits at a value `≥ 3` (THM-4555 (iv), THM-4556 (iv)). The **cluster** of `K` is its class under the transitive closure of certificates among exponents whose orbits reach 1.
* **Haar model** (THM-4581 6(d)): `x` is Haar-random in the post-run cell `{x ≡ 1, 5 mod 16}` and the children are `y_D = (x+1)/3^D − 1`, `D ≥ 1`; the `D`-chain is THM-4581's pair chain started at level `k_0 = D` with translation `3^D − 1`; `τ_D` is its absorption time (Terras steps), and

      N_T = #{ D ≥ 1 : τ_D ≤ T }.

  On the integer Mersenne line, `N_T` with `T = σ_T(x_K)` counts the cluster members below `K` ("members below `K`").

## Statements

1. **Clusters are level sets (INHERITED, PROVED in the barriers note, Prop. 4.1; HYP-9242's partner criterion).** For exponents whose orbits reach 1, `M_K` and `M_(K')`, `K' < K`, have a deletion certificate iff their deletion invariants agree. A certificate removes exactly `D` rising Terras steps and no odd step, so the invariant is preserved; conversely equal invariants force a merge at a value `≥ 8` (HYP-9242, audit C). Hence the clusters are exactly the level sets of the invariant, and the **orphans** (odd `K` with no deletion child merging before 1) are the minima of the level sets.

2. **Linearity (PROVED).** In the Haar model, for every `T`,

       E[N_T] = Σ_(D ≥ 1) P(τ_D ≤ T).

   No independence between the chains is used or available: all chains are driven by the one source `x`.

3. **Upper bound (PROVED).** `P(τ_D ≤ T) ≤ exp(−D²/(2T))`, hence

       E[N_T] ≤ Σ_(D≥1) exp(−D²/(2T)) < ∫_0^∞ exp(−x²/(2T)) dx = sqrt(π T / 2).

   So the expected number of merged deletion children is `O(sqrt T)`, with the explicit constant `sqrt(π/2) = 1.2533`.

4. **Lower bound (PROVED at the sketch level of THM-4581 (4)).** There are constants `c, T_0 > 0` with

       E[N_T] ≥ c · sqrt(T / log T)      for T ≥ T_0.

   Hence `E[N_T] = sqrt(T)^(1 + o(1))`: the exponent `1/2` of HYP-9243 holds for the mean number of merged children in the Haar model.

5. **The SRW skeleton and the constant (HEURISTIC + NUMERICAL).** THM-4581's proofs describe the level process at flip times as a simple random walk with fresh coins, flips being about half of the Terras steps (THM-4581 (7)). The probability that the level started at `D` reaches 0 within `T` Terras steps is then `≈ erfc(D / sqrt T)`, and

       Σ_D erfc(D / sqrt T) ≈ sqrt(T / π) = 0.5642 sqrt T

   is an upper envelope for `E[N_T]` (absorption also needs the translation to vanish, which costs `E[J] ≈ 10` further excursions for a fixed `D`, THM-4581 (7)). On the integer Mersenne line the ratio (mean members below `K`) / `sqrt(σ_T(x_K))` over the windows of statement 6 is `0.44, 0.26, 0.44, 0.47, 0.41, 0.51, 0.47, 0.56, 0.53, 0.62` (first window `3–400`, then dyadic windows to `50000`): the envelope `0.564` is matched to within the window noise from `K ≈ 800` on, with no sign of the `E[J]` deficit at the top of the range.

6. **The exact horizon to `K = 50000` (FINITE-EXACT).** The Terras orbits of every `M_K`, `2 ≤ K ≤ 50000`, were computed exactly (`horizon.c`, GMP, 6 threads, about 15 minutes); the columns `(K, odd, σ_T)` reproduce `mersenne_sigma_12800.txt` byte for byte for `K ≤ 12800`.
   * (a) **231 clusters.** The odd count alone separates them: the number of distinct `odd(M_K)` for `K ≤ 100, 400, 1200, 6000, 12800, 25000, 50000` is `22, 36, 57, 103, 136, 177, 231` (THM-4556 (vi) counted `K = 1` as well: `23, 37, 58, 104`), and it equals the number of clusters at every cutoff. Every even `K ≥ 4` shares its cluster with `K − 1` (THM-4556 (ii)); the least shift to a smaller cluster-mate is odd in all 24,769 cases (THM-4556 (iii)).
   * (b) **Orphans.** Among odd `K`, the orphan fraction times `sqrt(K_mid)` is `2.50, 1.59, 1.13, 1.04, 1.17, 0.92, 0.99, 0.84, 0.87, 0.76` on the windows `3–400, 401–800, 801–1600, 1601–3200, 3201–6400, 6401–12800, 12801–19200, 19201–25600, 25601–38400, 38401–50000`, and `0.97, 0.98, 0.95, 0.78, 0.86, 0.87, 0.93, 0.66, 0.79` on consecutive windows of 5000 from 5001. The least-squares exponent over the dyadic windows `400–50000` is `0.63`; the ratio is slowly decreasing, so the exponent is above `1/2` on this range (HYP-9242 fitted `0.59` with CI `[0.42, 0.75]` to 12800). There are 95 orphans in `(12800, 50000]`, listed in `analyze_horizon_extra.out`.
   * (c) **Cluster sizes.** For clusters whose minimum lies in the window, the mean size is `12.4, 32.0, 65.7, 90.4, 118.7, 220.2, 247.2, 368.2, 416.5, 534.7` on the ten windows above, i.e. `2/(mean size)·sqrt(K_mid)` is `2.29, 1.53, 1.06, 1.09, 1.17, 0.89, 1.02, 0.81, 0.86, 0.79` (the barriers identity: orphan fraction `= 2/(mean cluster size)` per odd `K`). Over clusters with minimum `≥ 1000`, `size / sqrt(K_min)` has median `2.05`, mean `2.17`, maximum `6.32` (the largest cluster: 1276 exponents, `40713 ≤ K ≤ 42002`, 14 holes).
   * (d) **Least-shift law.** Of the 24,769 odd `K ≤ 50000` with a smaller cluster-mate, `95.52%` have least shift `D = 1`; the counts for odd `D = 3, 5, …, 41` are `492, 204, 129, 75, 50, 35, 16, 13, 9, 11, 10, 5, 7, 4, 7, 5, 5, 2, 3, 2`, fitting `D^(−2.01)`; the largest least shift is `243`, and only 6 exceed 100.
   * (e) **Interleaving.** 23 of the 231 clusters are intervals of exponents; the others have holes (hole fraction `0.147` of the total span); the number of clusters whose span covers a given exponent is at most `3` (mean `1.1–1.4`). So the clusters are not intervals, but the interleaving is shallow.
   * (f) **Three terminal channels.** Every orbit with `K ≥ 3` enters 1 through `5` (46,186 exponents, 215 clusters), `85` (1,912 exponents, 8 clusters) or `341` (1,900 exponents, 8 clusters); no cluster mixes channels; `21, 1365, 87381` are unreachable (`3 | (4^m − 1)/3` for `3 | m`) and `5461, 21845` never occur. THM-4556 (ii) recorded `{5, 85, 341}` for odd `K ≤ 6000`.
   * (g) **Payment and peak.** The first Terras time at which the orbit drops below `M_K` satisfies `(pay − K)/K = 2.8213` with standard deviation `0.0714` over `100 ≤ K ≤ 50000` (minimum `1.887`, maximum `4.336`); the drift value for a random post-run orbit is `(log_2 3 − 1)/(1 − (log_2 3)/2) = 2.8188`. The peak bit length is `1.5851 K` on average (the run end `log_2 3 = 1.5850`), at most `1.6154 K`.
   * (h) **Orphans and the two-run.** By two-run length `J = 1 + ⌊v_2(K−1)/2⌋` (THM-4601), the orphan rate over odd `K ∈ [1001, 50000]` is `0.0056, 0.0085, 0.0113, 0.0078` for `J = 1, 2, 3, ≥ 4` (68, 78, 26, 6 orphans): a weak increase, no strong concentration of the uncovered branch in the deep reset-2 shells.

## Proofs

**2.** `N_T = Σ_D 1{τ_D ≤ T}`; take expectations. ∎

**3.** Absorption of the `D`-chain at `(0, 0)` requires its level `k` to reach 0. THM-4581, proof of 3(i), proves that `k` is a martingale with increments in `{−1, 0, 1}` (the coins at flip times are fresh; Durrett 4.3.1 is applied there). Started at `k_0 = D`, the event `{τ_D ≤ T}` is contained in `{min_(t≤T) (k_t − k_0) ≤ −D}`. The Azuma–Hoeffding maximal inequality for a martingale with increments bounded by 1 gives `P(min_(t≤T)(k_t − k_0) ≤ −D) ≤ exp(−D²/(2T))`. Summing over `D ≥ 1` and comparing with the integral proves the bound; `∫_0^∞ e^(−x²/(2T)) dx = sqrt(πT/2)`. ∎

**4 (sketch level).** THM-4581 (4), proof, upper bound: for a start at level `k_0 = D`, `P(no merge by T) ≤ C (log T)² T^(−1/2) + C' D (log T / T)^(1/2)` (the displayed union bound plus the first-passage tail from `|k_0|`, whose killed-walk Green function is `2 min(|k_0|, h)`). Put `D_T = ⌊ sqrt T / (4 C' sqrt(log T)) ⌋`. For `T ≥ T_0` and `D ≤ D_T` the right side is `≤ 1/2`, so `P(τ_D ≤ T) ≥ 1/2` and `E[N_T] ≥ D_T / 2`. ∎

**5.** The skeleton description is THM-4581, proof of 2 ("the excursion"): the flip skeleton `|k|` is a simple random walk with fresh coins. With `n ≈ T/2` flips, `P(SRW from D hits 0 within n steps) ≈ 2 P(S_n ≥ D) ≈ erfc(D / sqrt(2n)) = erfc(D / sqrt T)`, and `∫_0^∞ erfc(x / sqrt T) dx = sqrt(T/π)`. This is a heuristic (the walk and the flip count are not exactly those of a fair SRW on `T/2` steps), and absorption needs more than reaching level 0. The data column is `mean_below_over_sqrtT` in `analyze_horizon.out`. ∎

**6.** `horizon.c` runs the Terras map with GMP on `2^K − 1`, counting odd and total steps, recording the last odd value before 1 (as `j` with `3x + 1 = 2^j`), the peak bit length and the first time the orbit is below the source; `analyze_horizon.py` and `analyze_horizon_extra.py` form the level sets of `(σ_T − K, odd)` and compute every number above. The agreement with `mersenne_sigma_12800.txt` is an exact `diff`. ∎

## Remarks

* **Coverage versus grounding.** On this range, `99.5%` of odd exponents above 5000 have a deletion child, and the least shift is 1 in `95.5%` of cases: the single-step repair `M_K ⇝ M_(K−1)` almost always succeeds, longer blocks (`D ≥ 3`, law `D^(−2)`) repair most of the rest, and orphans are the residue (`0.4%` near 50000). None of this grounds an exponent: a certificate preserves `σ_T(x_K)` (barriers Prop. 4.1) and zero-defect receipts cannot ground exponents above the largest supplied terminal rank (codex's finite-anchor closure bound). Universal coverage of the Mersenne line by Mersenne children is therefore a statement about density (→ 1 by THM-4581 6(d) and statement 4), never about ROOT.
* **Blocks that succeed together.** A cluster is the set of exponents that share one future; when the `D = 1` chain of `M_K` is absorbed, every exponent already in the cluster of `K − 1` joins at that moment (the barriers note's absorbed groups of 1916 deletions). The least-shift law measures how often the single step fails and a longer block succeeds.
