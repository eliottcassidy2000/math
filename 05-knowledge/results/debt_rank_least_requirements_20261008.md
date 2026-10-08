# The debt-lattice rank with the least proof requirements: every expanding map escapes, Z_2 is classified, and translation-only maps of rank three or more never coalesce (THM-4607, THM-4608, THM-4609)

Session `mac-mini-2026-10-08-rank`, 2026-10-08.

Brief from the owner: "keep investigating the structure of the debt lattice rank and coming up with creative new angles to satisfy the least proof requirements".

Scripts are in `04-computation/experiments/rank_20261008/`.

## 0. The answer in one screen

**The angle.** HYP-9244 says coalescence of affinely related Haar orbits = accessibility × contraction × recurrence of the debt walk, with the rank of the debt lattice deciding recurrence. Its hard parts all seemed to need control of the coupling state `(M mod d, e mod d)`, an adaptive and non-Markov process. The least-requirement route avoids that control:
- **Expansion** is felt through a normalization (`F = |e|/(1 + M)`) whose conditional log-drift is exactly `Λ`, whatever the coupling. Only a convexity correction near `M = 1` remains, and a martingale with increments bounded away from 0 spends `O(√n)` moves there.
- **Transience** can be forced one step at a time, provided every possible one-step law of the debt is controlled by one quadratic form. For translation-only maps (all `m_i ≡ 1 mod d`) the possible laws are a short explicit list, the cyclic-root laws, so an adversary choosing the coupling at every step still cannot stop the escape.

**Proved.**
- **THM-4607.** Every expanding Matthews–Watts map, of any digit base and any rank, fails to coalesce: `q(e) → 0` as `|e| → ∞`.
- **THM-4608.** The complete classification on `Z_2`. Translates coalesce almost surely iff the map contracts (`m_0 m_1 ≤ 3`) and the offset is accessible (explicit divisibility).
  - The tail is geometric for `x + s` and diffusive for `3x + s`.
  - Up to affine conjugacy `3x + 1` is the only 2-adic map with diffusive coalescence.
- **THM-4609.** Translation-only debt walks.
  - The covariance under lag `b` is `(1/d)(2I − A_b)`, with `A_b` the lag-`b` cycle with the unit positions cut out: a union of paths for prime `d`.
  - Path spectra stay below 4, so from rank 4 on the standard form is a Lamperti form. Rank 3 has two universal families (AP / non-AP) with explicit integer forms.
  - Hence every translation-only map with independent multipliers of rank `≥ 3` on a prime base `d ≥ 5` has a transient debt walk. That includes contracting maps: `x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5` contracts (`Λ = −0.217`), yet `y` and `y + e` merge with probability `< 1`.
  - Rank zero coalesces almost surely iff a finite offset chain reaches 0.

**The structural picture.**
- The debt walk is a random walk on (a projection of) the root lattice `A_(d−1)`. Its step under coupling `π` is a root `v_π(j) − v_j`, and its covariance is the transported Laplacian of the cycle graph of `π`.
- For `d = 2` the rank is at most 1, which is why Collatz-type coalescence is recurrent and diffusive.
- From `d = 5` on, contraction and coalescence come apart.

## 1. THM-4607: expansion beats any debt walk

**The normalization.**
- With `F = |e|/(1 + M)` and `x = ln M`:
  - `F' ≥ μF − R/d`;
  - `ln μ = (1 − g'(x)) ln(m_i/d) + g'(x) ln(m_j/d) − ρ`, with `g(y) = ln(1 + e^y)`.
- Given the past, both digits `i` and `j` are uniform. So the first two terms have conditional mean exactly `Λ`, whatever the coupling.
- The convexity remainder satisfies `0 ≤ ρ ≤ (Δ²/2) min(1/4, e^(Δ − |x|))`, and `ρ = 0` unless the debt moves.

**The skeleton.**
- The debt at its successive moves is a martingale with increments of absolute value in `[δ_min, Δ]`.
- For the even convex `H` with `H'' = κ min(1/4, e^(2Δ − |y|))`, the bounded slope gives `E Σ_(k<N) min(1/4, e^(Δ − |X_k|)) ≤ C√N`.
- So the total convexity loss is `O(√n)` in expectation, and `≤ εn + C` uniformly with high probability (dyadic Markov bounds).
- Azuma–Hoeffding handles the martingale part.

**What it uses.** Only the uniformity of both digits and bounded nonzero debt moves. There is no rank condition, no mixing and no accessibility. This is THM-4606's argument with the simple-random-walk skeleton replaced by any bounded martingale skeleton.

## 2. THM-4608: the 2-adic world is closed

| multipliers `(m_0, m_1)` | class | conjugate to | coalescence of `y`, `y + e` |
|---|---|---|---|
| (1, 1) | rank 0, contracting | `x/2, (x + s_0)/2`, `s_0 = r_1 − r_0` | a.s. iff `s_0 \| e`; geometric tail |
| (1, 3), (3, 1) | rank 1, contracting | `3x + s`, `s = r_0 + r_1` | a.s. iff `s/3^(v_3 s) \| e`; diffusive tail |
| `m_0 m_1 ≥ 5` | expanding | — | `q(e) → 0` (THM-4607) |

**Tools.**
- Translations, including a parity-swapping one for `(3, 1)`, preserve offsets and Haar measure.
- Scalings `x ↦ s x` rescale offsets.
- An `ℓ`-adic valuation argument shows that offsets with an odd prime `ℓ ∤ 3` in the denominator of `e/s` never reach 0.
- The departure rule needed for the expanding maps was checked exactly (`z2_step_table.py`): at a departure, both coins give `min(m_0, m_1)/2`.
  - A first draft said "`m_0/2`". It is the smaller multiplier.

## 3. THM-4609: translation-only debt walks

### 3.1 The cyclic-root walk and path Laplacians

- If `m_i ≡ 1 (mod d)` for all `i`, then `M ≡ 1 (mod d)` and `u`'s digit is `j + ē`. The debt step is uniform over the lag-`ē` roots.
- For independent multipliers placed at positions off a set `Z` of unit positions:
  - `C_b = (1/d)(2I − A_b)`;
  - the diagonal is `2/d` (each position lies in two roots);
  - off-diagonal entries are `−1/d` for each pair at distance `±b`.
- For prime `d`, deleting `Z` from the `d`-cycle `j ↦ j + b` leaves paths. Their adjacency spectra `2cos(πm/(n+1))` are `> −2`.
- So `λ_max(C_b) < 4/d` and `tr C_b = 2ρ/d`, which settles rank `≥ 4` with `Q = I`.
- Only differences of the multiplier vectors enter, so an empty `Z` is the case `|Z| = 1` shifted.

### 3.2 Rank 3: two universal families

| type | lag graphs on the three positions | `Q = I` | explicit form (verified exactly) |
|---|---|---|---|
| AP (`p, p+c, p+2c`; every 3-set of `Z/5`) | path `P_3` at `±c`, end-to-end edge at `±2c`, edgeless otherwise | equality (fails) | `[[6,−2,0],[−2,6,−1],[0,−1,5]]` (end, middle, end) |
| non-AP (distinct pair differences) | one edge at each of `±δ_12, ±δ_23, ±δ_13`, edgeless otherwise | equality (fails) | `[[4,−1,−1],[−1,5,0],[−1,0,5]]` |

- The families are scale-free (the factor `1/d` cancels), so the same forms serve every prime `d ≥ 5`.
- **Census.** Every AGL(1, d)-orbit of unit-position sets with rank `≥ 3` has an exactly verified balanced form:
  - `d = 5`: 2 orbits;
  - `d = 7`: 6 orbits;
  - `d = 11`: 24 orbits, ranks 3–7.
- **Map-level checks on `Z_5`.** Whitened search verified 49 of the 54 contracting rank-3 maps directly. The other 5 are search misses of the same configuration type, covered by the type's exact form: balance is invariant under change of basis.
- **Correction to an earlier attempt.** A first, coordinate-dependent random search reported margin "0" for many maps (`balanced_lamperti.out`, partial). Equivalent configurations got different answers, which exposed the search, not the mathematics.

### 3.3 Lamperti and the adversary

- With `V(x) = (x^T Q^(−1) x)^(−α/2)` and `0 < α < min_b tr/λ_max − 2` (for the `Z_5` type, `α_max ≈ 0.21`), `V` is a supermartingale far out under every lag, including lag 0, which does not move.
- So `P(return within r_1 | debt at distance r) ≤ (r_1/(r − B))^α`.
- A positive-probability finite digit word carries any non-absorbed start to large debt without absorption.
- **Adversarial test** (`adversarial_lag.py`, `Z_5` AP type, 300 runs of 4000 steps):

  | lag chooser | mean returns to the origin | median `|x|_A` at the end |
  |---|---|---|
  | worst lag each step | 1.82 | 24.2 |
  | uniform lags | 0.52 | 24.5 |

  The adversary costs a few early returns but not the escape.

### 3.4 A contracting map whose translates fail to merge

`T(x) = x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5`, with `Λ = −0.217` and rank 3 (6, 11, 16 independent):

| quantity | value |
|---|---|
| random 2000-digit `n` merging with `n + 1` at equal time, within 100 / 1000 / 4000 steps | 12.8% / 13.0% / 13.5% (flat) |
| `n ≤ 20000`: cycles reached | all orbits reach one of 4 cycles (minimal elements 1, 8, 9, 33) |
| `n ≤ 20000`: consecutive pairs merging at equal time | 28.6% |

**The picture.** The map contracts and every small orbit enters a cycle. But equal-time merging, the coalescence that drives Collatz's tree structure, fails for most pairs: their debts drift apart in three dimensions.

## 4. What remains (the minimal open requirements)

1. **Rank 1, `d ≥ 3`.**
   - Two-valued translation-only maps have an exactly fair simple-random-walk skeleton.
   - The offset multiplier's conditional `θ`-moment is exactly `κ(θ) = (1/d)Σ(m_i/d)^θ < 1` at every step away from zero debt, whatever the lag; at zero debt it is at most `κ(θ)`.
   - So THM-4581's architecture (a level weight `s^|k|`, runs controlled by `v_d(μ^h − 1)`, the SRW Green function) should give a.s. coalescence for accessible starts. Not yet written.
2. **Rank 2.** The critical two-dimensional case. The two lag families `{2I − E, 2I}` cannot share an isotropizing form. Recurrence needs decorrelation between the lag process and the debt direction.
3. **General couplings, rank `≥ 3`.** Involutions `j ↦ −j + b` give low-rank covariances, so a block Lamperti form is needed.
4. **Accessibility.** Does the absence of congruence obstructions suffice?

## 5. Reproduction

All runs were in `04-computation/experiments/rank_20261008/`, with Python 3.10, standard library only.

| command | runtime | what it does |
|---|---|---|
| `python3 z2_step_table.py` | seconds | exact 2-adic step table; departures give `min(m_0, m_1)/2` |
| `python3 rank3_universal.py` | about 5 min | explicit rank-3 forms and the rank `≥ 4` spectral bound |
| `python3 balanced_types.py` | about 10 min | AGL-orbit census for `d = 5, 7, 11` |
| `python3 z5_types_exact.py` | about 1 min | |
| `python3 adversarial_lag.py` | about 5 min | |
| `python3 z5_example_integers.py` | about 10 min | |

- `balanced_lamperti.py` and `balanced_whitened.py` hold the per-map searches. Their `.out` files are partial or superseded by `balanced_types.py`, and are kept as records.
