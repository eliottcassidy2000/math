# The debt-lattice rank with the least proof requirements: every expanding map escapes, Z_2 is classified, and one-step Lamperti forms make debt walks of rank three or more transient (THM-4607, THM-4608, THM-4609)

Session `mac-mini-2026-10-08-rank`, 2026-10-08. Independently audited by audit G (CORRECT WITH FIXES; MISTAKE-590). Sections 3.4, 4 and 5 were added after the audit and are not independently audited.

Brief from the owner: "keep investigating the structure of the debt lattice rank and coming up with creative new angles to satisfy the least proof requirements".

Scripts are in `04-computation/experiments/rank_20261008/`.

## 0. The answer in one screen

**The angle.** HYP-9244 says coalescence of affinely related Haar orbits = accessibility × contraction × recurrence of the debt walk, with the rank of the debt lattice deciding recurrence. Its hard parts all seemed to need control of the coupling state `(M mod d, e mod d)`, an adaptive and non-Markov process. The least-requirement route avoids that control:
- **Expansion** is felt through a normalization (`F = |e|/(1 + M)`) whose conditional log-drift is exactly `Λ`, whatever the coupling. Only a convexity correction near `M = 1` remains, and a martingale with increments bounded away from 0 spends `O(√n)` moves there.
- **Transience** can be forced one step at a time, provided every possible one-step law of the debt is controlled by one quadratic form. Every possible law is the root law of an affine coupling `j ↦ aj + b`, a short explicit list. So an adversary choosing the coupling at every step still cannot stop the escape. This is the trace criterion of Peres, Popov and Sousi (2013, Theorem 1.3), applied with the lag-0 freeze.

**Proved.**
- **THM-4607.** Every expanding Matthews–Watts map, of any digit base and any rank, fails to coalesce: `q(e) → 0` as `|e| → ∞`. Also `q(e) < 1` at every offset for `px + 1` and for rank zero.
- **THM-4608.** The classification on `Z_2` for positive multipliers. Contracting maps (`m_0 m_1 ≤ 3`) coalesce almost surely iff the offset is accessible (explicit divisibility).
  - The tail is geometric for `x + s` and diffusive for `3x + s`.
  - Among positive-multiplier 2-adic maps, up to affine conjugacy, only the `3x + s` class coalesces at the diffusive rate.
- **THM-4609.** Root-law debt walks.
  - Under coupling `π` the covariance is `(1/d)(2I_nf − A_π)`: the Laplacian of the cycle graph of `π` with the unit positions cut out. Its spectrum lies in `[0, 4]`, with 4 only for even cycles.
  - For prime `d` and independent multipliers, the standard form `Q = I` balances every coupling from rank 4 (translation-only maps), rank 5 (odd-order coupling groups) and rank 6 (every map). Rank-3 translation-only maps have two universal families with explicit integer forms.
  - Hence contracting maps whose translates fail to merge: `x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5` (`Λ = −0.217`), and on `Z_7` the generic map with multipliers `1, 2, 3, 5, 11, 13, 17` (`Λ = −0.346`).
  - Rank zero coalesces almost surely iff a finite offset chain reaches 0.

**The structural picture.**
- The debt walk is a random walk on (a projection of) the root lattice `A_(d−1)`. Its step under coupling `π` is a root `v_π(j) − v_j`, and its covariance is the transported Laplacian of the cycle graph of `π`.
- What costs dimension: a fixed point of the coupling removes one unit of trace, and an even cycle lets the top eigenvalue reach 4.
- For `d = 2` the rank is at most 1, which is why Collatz-type coalescence is recurrent and diffusive. From `d = 5` on, contraction and coalescence come apart.

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

**Small offsets.** `q(e) < 1` at every offset needs a positive-probability word to a large normalized offset. For rank zero a greedy coin works for every `d`: the additive terms `r_(j+ē) − r_j` sum to 0, so some digit gives `|e'| ≥ (m/d)|e|`. For `px + 1` it is THM-4606. For other expanding maps it is numerical (audit G: `q ≤ 0.28` in every tested case).

**What it uses.** Only the uniformity of both digits and bounded nonzero debt moves. There is no rank condition, no mixing and no accessibility. This is THM-4606's argument with the simple-random-walk skeleton replaced by any bounded martingale skeleton.

## 2. THM-4608: the 2-adic world (positive multipliers)

| multipliers `(m_0, m_1)` | class | conjugate to | coalescence of `y`, `y + e` |
|---|---|---|---|
| (1, 1) | rank 0, contracting | `x/2, (x + s_0)/2`, `s_0 = r_1 − r_0` | a.s. iff `s_0 \| e`; geometric tail |
| (1, 3), (3, 1) | rank 1, contracting | `3x + s`, `s = r_0 + r_1` | a.s. iff `s/3^(v_3 s) \| e`; diffusive tail |
| `m_0 m_1 ≥ 5` | expanding | — | `q(e) → 0` (THM-4607); `q(e) < 1` at every `e` for `px + 1` and rank 0 |

**Tools.**
- Translations, including a parity-swapping one for `(3, 1)`, preserve offsets and Haar measure.
- Scalings `x ↦ s x` rescale offsets.
- An `ℓ`-adic valuation argument shows that offsets with an odd prime `ℓ ∤ 3` in the denominator of `e/s` never reach 0.
- THM-4581 3′ (offsets in `Z[1/3]`) gives the accessible case of the Collatz class; `3x ± 3^a` make every integer offset accessible.
- The departure rule needed for the expanding maps was checked exactly (`z2_step_table.py`): at a departure, both coins give `min(m_0, m_1)/2`. A first draft said "`m_0/2`".

**Negative multipliers are outside.** Matthews and Watts allow them. `(1, −3)` and `(−1, 3)` contract with rank one, are not conjugate to `3x + s` (branch slopes are conjugacy invariants), and are numerically diffusive (audit G).

## 3. THM-4609: root-law debt walks

### 3.1 Root laws and cycle-graph Laplacians

- At each step the coupling is affine: `u`'s digit is `π(j) = M̄ j + ē`, with `M̄` in the coupling group `G = ⟨m_i/m_j mod d⟩`. The debt step is uniform over the roots `v_π(j) − v_j`.
- Translation-only maps (`G = {1}`, i.e. all `m_i` congruent mod `d`) only ever use translations.
- For independent multipliers placed off a set `Z` of unit positions, `C_π = (1/d) D_π` with `D_π = 2I_nf − A_π`:
  - the diagonal is 2 at non-unit positions moved by `π` (each lies in two roots) and 0 at fixed ones;
  - `A_π` is the adjacency, with multiplicity, of the cycle graph of `π` on the non-unit positions (a transposition gives a double edge).
- Spectra: `2 − 2cos(2πk/ℓ)` on an `ℓ`-cycle, `2 − 2cos(πk/(n+1))` on a path. So `λ_max ≤ 4`, with 4 only on even cycles; `tr D_π = 2·#(non-unit positions moved)`.
- For prime `d`, a translation is one odd `d`-cycle; deleting `Z` leaves paths. Only differences of the multiplier vectors enter, so an empty `Z` is the case `|Z| = 1` shifted.

### 3.2 Rank 3, translation-only: two universal families

| type | lag graphs on the three positions | `Q = I` | explicit form (verified exactly) | `α_max` |
|---|---|---|---|---|
| AP (`p, p+c, p+2c`; every 3-set of `Z/5`) | path `P_3` at `±c`, end-to-end edge at `±2c`, edgeless otherwise | fails strictly at the `P_3` lags (`λ_max − ½ tr = √2 − 1`), equality at the single-edge lags | `[[6,−2,0],[−2,6,−1],[0,−1,5]]` (end, middle, end) | 0.039 |
| non-AP (distinct pair differences) | one edge at each of `±δ_12, ±δ_23, ±δ_13`, edgeless otherwise | equality (fails) | `[[4,−1,−1],[−1,5,0],[−1,0,5]]` | 0.078 |

- The families are scale-free (the factor `1/d` cancels), so the same forms serve every prime `d ≥ 5`. Audit G checked all 12,065 three-point configurations with `5 ≤ d ≤ 31` prime.
- **Census.** AGL(1, d)-orbits of unit-position sets with rank `≥ 3`, each with an exactly verified balanced form: `d = 5`: 2 of 2; `d = 7`: 6 of 6; `d = 11`: 26 orbits, 22 verified in `balanced_types.out` (ranks 8–10 skipped there), all 26 verified by audit G. A first version said "24 orbits".
- **Map-level checks on `Z_5`.** Whitened search verified 49 of the 54 contracting rank-3 maps directly. The other 5 are search misses of the same configuration type, covered by the type's exact form: balance is invariant under change of basis.
- **Correction to an earlier attempt.** A first, coordinate-dependent random search reported margin "0" for many maps (`balanced_lamperti.out`, partial). Equivalent configurations got different answers, which exposed the search, not the mathematics.

### 3.3 Lamperti and the adversary

- With `V(x) = (x^T Q^(−1) x)^(−α/2)` and `0 < α < min_π tr/λ_max − 2`, `V` is a supermartingale far out under every coupling, including the identity, which does not move. The values of `α_max` are per form: 0.039 (`Q_AP`), 0.078 (`Q_nonAP`), about 0.21 for a form valid only on `Z_5`.
- So `P(return within r_1 | debt at distance r) ≤ (r_1/(r − B))^α`.
- A positive-probability finite digit word carries any non-absorbed start to large debt without absorption, so absorption has probability `< 1`.
- **Adversarial test** (`adversarial_lag.py`, `Z_5` AP type, 300 runs of 4000 steps):

  | lag chooser | mean returns to the origin | median `|x|_A` at the end |
  |---|---|---|
  | worst lag each step | 1.82 | 24.2 |
  | uniform lags | 0.52 | 24.5 |

  The adversary costs a few early returns but not the escape. Audit G ran six adversaries (no late returns in rank 3).

### 3.4 Beyond translations: the standard-form hierarchy (added after audit G)

For prime `d`, independent multipliers of rank `ρ`, and a coupling `π ≠ id`:
- an affine map with `a ≠ 1` has exactly one fixed point, so `½ tr D_π ≥ ρ − 1`; a translation has none, so `½ tr D_π = ρ`;
- if the coupling group has odd order (`−1 ∉ G`), every coupling has only odd cycles (length `ord(a)`), so `λ_max < 4`.

| coupling class | standard form works from | reason |
|---|---|---|
| translations | rank 4 | `λ_max < 4 ≤ ρ = ½ tr` |
| odd-order coupling group | rank 5 | `λ_max < 4 ≤ ρ − 1 ≤ ½ tr` |
| every map | rank 6 | `λ_max ≤ 4 < 5 ≤ ρ − 1 ≤ ½ tr` |

- `rank_hierarchy.py` checks every unit set of each rank for `d = 7, 11, 13` exactly. The standard form fails somewhere at ranks 3, 4 and 5 respectively, so the thresholds are sharp for it.
- Below the thresholds explicit forms can exist (`squares_group_examples.py`, exact): for `G` the squares, rank 4, `Q = 8I − J` on `Z_7` (units `{0, 1, 2}`), and integer forms on `Z_11` (units `{0..6}` and `{0..5, 7}`).
- One-step balance needs every coupling covariance to have rank ≥ 3 (`tr Σ ≤ rank·λ_max`). Involutions `j ↦ −j + b` have rank ≤ `(d − 1)/2`, so on `Z_5` every map with `−1 ∈ G` (every non-translation map there) is out of reach of one-step forms. `coupling_groups.out` lists searches by subgroup; its "no" entries are proofs only when the minimal covariance rank is ≤ 2. Its `|G| = 1` "no" at `d = 7`, units `(0, 1, 2, 3)`, is a search miss: `Q_AP` balances that configuration exactly.

### 3.5 Contracting maps whose translates fail to merge

| `d` | multipliers by residue | `Λ` | rank | coupling class | form |
|---|---|---|---|---|---|
| 5 | `(1, 1, 6, 11, 16)` | −0.217 | 3 | translation-only, AP | `Q_AP` |
| 7 | `(1, 1, 1, 1, 8, 15, 22)` | −0.820 | 3 | translation-only, AP | `Q_AP` |
| 7 | `(1, 1, 1, 8, 15, 22, 29)` | −0.339 | 4 | translation-only | `Q = I` |
| 7 | `(1, 1, 1, 2, 11, 23, 29)` | −0.575 | 4 | `G = {1, 2, 4}` | `8I − J` |
| 7 | `(1, 1, 2, 11, 23, 29, 37)` | −0.060 | 5 | `G = {1, 2, 4}` | `Q = I` |
| 7 | `(1, 2, 3, 5, 11, 13, 17)` | −0.346 | 6 | all units | `Q = I` |

In every row `r_i ≡ −m_i i (mod d)`. Direct integer checks for the `Z_5` map `x/5, (x + 4)/5, (6x + 3)/5, (11x + 2)/5, (16x + 1)/5` (`z5_example_integers.py`):

| quantity | value |
|---|---|
| random `n` with 2000 base-5 digits merging with `n + 1` at equal time, within 100 / 1000 / 4000 steps (400 samples each) | 12.8% / 13.0% / 13.5% (flat) |
| `n ≤ 20000`: cycles reached | all orbits reach one of 4 cycles (minimal elements 1, 8, 9, 33) |
| `n ≤ 20000`: consecutive pairs merging at equal time | 28.6% |

Audit G: `0.140 ± 0.008` at `T = 16384` from offset 1, and ten other offsets flat at 0.045–0.274.

**The picture.** The map contracts and every small orbit enters a cycle. But equal-time merging, the coalescence that drives Collatz's tree structure, fails for most pairs: their debts drift apart in three dimensions.

## 4. Accessibility: the twisted obstruction (audit G)

- **Lemma (PROVED).** Let `ℓ ∤ d` be prime, all `m_i ≡ d (mod ℓ)`, and `r_i ≡ d·ψ(m_i) (mod ℓ)` for a homomorphism `ψ` from the multiplicative group generated by the `m_i` to `Z/ℓ`. Then `e_n − ψ(M_n) ≡ e_0 (mod ℓ)`. So if `ℓ ∤ e_0`, the chain never reaches `(1, 0)`.
- Both this lemma and audit F's constant-`c` lemma are cocycle invariants: `e − f(M)` is multiplied by `m_i/d` mod `ℓ`, with `f(M) = (1 − M)c` or `f = ψ`. The twisted `f` reads the exponents of `M`, not its residue.
- **Examples.**
  - `Z_3`, m = (1, 5, 7), r = (0, 1, 1): the "1.000" expanding row of HYP-9244 was this obstruction at `ℓ = 2`. Even offsets merge: `q(2) = 0.032`.
  - `Z_3`, m = (1, 1, 5), r = (0, 2, 5): contracting (`Λ = −0.562`), rank one, no constant-`c` obstruction for any admissible `ℓ`, and odd offsets never merge (`twisted_check.out`: 1,971,984 chain steps without an invariant violation; offsets 1 and 3 merge 0 of 300 times by `T = 3000`, offsets 2 and 4 merge 286 and 278 times).
- So HYP-9244's question "does the absence of congruence obstructions suffice?" has a negative answer for the constant-`c` kind. The refined question is whether every inaccessible start is detected by some cocycle invariant mod `ℓ^k`.

## 5. Rank one on `Z_3` and `Z_5`: numerics for the open row

`rank1_two_valued.py` (exact pair chain from `(1, 1)`, 1500 chains, `T = 16 … 16384`):

| map | `μ` vs `d²` | `Λ` | `√T·P(no merge by T)` |
|---|---|---|---|
| `Z_3` (1, 1, 4) | 4 < 9 | −0.637 | 1.37–1.64 |
| `Z_5` (1, 1, 1, 1, 26) | 26 > 25 | −0.958 | 1.20–1.96 |
| `Z_5` (1, 1, 1, 1, 6) | 6 < 25 | −1.251 | 0.79–1.02 |
| `Z_3` (1, 4, 4) | — | −0.174 | 2.5 rising to 7.9, 9.1, 8.1 at `T = 1024, 4096, 16384` |

All four are consistent with the `T^(−1/2)` tail. The sketch route (THM-4581's architecture with a level weight `s^|k|`) works for every `μ`: at `s = 1` the per-step weighted moment equals `κ(θ) < 1` for every lag, so `s` just below 1 leaves room (audit G). An earlier version of this note restricted the route to `μ ≤ d²`; that was wrong.

## 6. What remains (the minimal open requirements)

1. **Rank 1, `d ≥ 3`.** Write the THM-4581 transfer for two-valued translation-only maps (section 5), then general rank-one maps.
2. **Rank 2.** The critical two-dimensional case. The two lag families `{2I − E, 2I}` cannot share an isotropizing form. Recurrence needs decorrelation between the lag process and the debt direction.
3. **Ranks 3–5 with an even-order coupling group** (all non-translation maps on `Z_5`, such as `Z_5` (1, 2, 3, 7, 1)). One-step forms are impossible when an involution's covariance has rank ≤ 2; a multi-step form must use how the coupling state moves.
4. **Dependent multipliers and composite `d`** (e.g. `d = 4` with multipliers 1, 3, 5, 7, contracting of rank 3).
5. **Accessibility.** Is every inaccessible start detected by a cocycle invariant mod some `ℓ^k`?

## 7. Audit G and the corrections (MISTAKE-590)

Audit G (independent, 2026-10-08; report in `audit_G/REPORT.md`) found the mathematics of THM-4607–4609 correct and required these fixes, all applied:
- THM-4608 first claimed the per-offset "iff" for expanding maps; only `q(e) → 0`, and `q(e) < 1` for `px + 1` and rank zero, are proved. It also omitted "positive multipliers".
- HYP-9244's update missed the twisted obstruction (section 4), carried a wrong "`μ ≤ d²`" caveat for rank one, and said one-step balance fails for every involution of rank ≤ `⌊d/2⌋` (true only when that rank is ≤ 2).
- THM-4609: the Lamperti criterion is Peres–Popov–Sousi's Theorem 1.3; the census at `d = 11` is 26 orbits (22 author-verified); `Q = I` fails strictly at the AP `P_3` lags; `α_max` is per form; "2000-digit" meant 2000 base-5 digits with 400 samples.
- Process: audit G downloaded a survey PDF (about 450 KB) without asking, against the download-permission rule, then deleted it.

## 8. Reproduction

All runs were in `04-computation/experiments/rank_20261008/`, with Python 3.10, standard library only.

| command | runtime | what it does |
|---|---|---|
| `python3 z2_step_table.py` | seconds | exact 2-adic step table; departures give `min(m_0, m_1)/2` |
| `python3 rank3_universal.py` | about 5 min | explicit rank-3 forms and the rank `≥ 4` spectral bound |
| `python3 balanced_types.py` | about 10 min | AGL-orbit census for `d = 5, 7, 11` |
| `python3 rank_hierarchy.py` | not recorded | standard-form hierarchy, exact, `d = 7, 11, 13` |
| `python3 squares_group_examples.py` | not recorded | exact rank-4 forms for the squares mod 7 and 11 |
| `python3 coupling_groups.py` | not recorded | form searches by coupling subgroup (re-runs `balanced_types`' census on import) |
| `python3 twisted_check.py` | not recorded | the twisted invariant along the chain; merge counts by offset |
| `python3 rank1_two_valued.py 1500 16384` | not recorded | rank-one tails |
| `python3 z5_types_exact.py` | about 1 min | |
| `python3 adversarial_lag.py` | about 5 min | |
| `python3 z5_example_integers.py` | about 10 min | |

- `balanced_lamperti.py` and `balanced_whitened.py` hold the per-map searches. Their `.out` files are partial or superseded by `balanced_types.py`, and are kept as records.
- `audit_G/` holds the auditor's scripts and outputs (`forms_exact.py`, `run_z5.py`, `run_z7_hunt.py`, `adversary.py`).
