# Coalescence needs accessibility, contraction and recurrence: the debt-lattice rank decides when affinely related Collatz-type orbits merge (THM-4606, HYP-9244), with a battery of quick reframe tests on the Mersenne line, the +1 barrier and periodic anchors

Session `mac-mini-2026-10-08-reframes`, 2026-10-07/08.

Brief from the owner: creative hypothesis generation and testing; look through past results for inspiration; find reframes that unlock Collatz-adjacent proofs.

Scripts are in `04-computation/experiments/reframes_20261007/`.

## 0. The answer in one screen

**The reframe.**
* THM-4581 (two 2-adic orbits related by `u = 3^k v + e` merge almost surely, with tail `T^(-1/2)`) has two separate inputs:
  1. **Contraction of the offset.** At zero debt, the difference `e = u − v` is multiplied by 1/2 or 3/2 whenever the parities agree. Its geometric mean is `√3/2 < 1`.
  2. **Recurrence of the debt walk.** The debt `k` is a one-dimensional simple random walk in flip time.
* Replace 3 by another multiplier, or Collatz by a generalized map on `Z_d`, and the two inputs come apart. Each can fail on its own.
* A third input, implicit for `3x + 1` with integer offsets, is **accessibility**. Audit F found congruence obstructions: if a prime `ℓ ∤ d·Πm_i` and some `c` satisfy `(m_i − d)c + r_i ≡ 0 (mod ℓ)` for every branch, then the merge state is unreachable from offsets prime to `ℓ`. For example, under `3x + 5`, `y` and `y + 1` never meet.

**Proved (THM-4606).**
* For `px + 1` with `p ≥ 5`, the orbits of Haar `y` and `y + e` meet with probability `q_p(e) < 1`, and `q_p(e) → 0` as `|e| → ∞`.
* For `p = 1, 3` they meet almost surely.
* So Haar coalescence among the `px+1` maps holds exactly at the contracting multipliers. Through the conjugacy `x ↦ rx`, the same holds for `px + r` with translations by multiples of `r`. Other offsets can be obstructed: `3x + 5` with offset 1 never merges.
* Numerically:
  * `q_5 = 0.00833`, `q_7 = 0.1086`, `q_9 = 0.00225`, `q_11 ≈ 2.43e-5`, `q_13 ≈ 1.6e-6` (audit E: exact enumeration plus Monte Carlo tail; rigorous floors 0.004387, 0.10727, 0.002180, 2.356e-5, 1.059e-6);
  * the spikes at `p = 2^a − 1` come from a reset path of length `a + 1`.

**Conjectured, with all five regimes observed (HYP-9244).**
* For generalized (Matthews–Watts) maps and accessible starts (no congruence obstruction; the obstruction lemma is PROVED), coalescence holds iff the map contracts and the debt walk on the multiplicative group `Γ = ⟨m_i/m_j⟩` is recurrent.
* When the coupling states equidistribute, the rank of `Γ` decides it, as in Pólya's theorem:

  | rank | coalescence | tail of P(no merge by T) |
  |---|---|---|
  | 0 | almost sure | exponential |
  | 1 | almost sure | `T^(-1/2)` |
  | 2 | almost sure | `1/log T` |
  | ≥ 3 | fails with positive probability, even for contracting maps | — |

* Expanding maps fail at every rank.
* Caveats: the debt walk's step law is chosen adaptively, and adaptive choice can change the type (Peres–Popov–Sousi 2013, for laws with invertible covariance; zero-drift walks in `d = 2` can be transient, Georgiou–Menshikov–Mijatović–Wade 2016). The conjecture is that this arithmetic rule does not conspire. Audit F's anisotropy test (effective dimension = rank, uniform coupling states) supports this on five maps.

**Readings.**
* **Collatz is rank one in its natural clocks.**
  * Terras clock: debt `3^Z`, with the 2-powers spent on time.
  * Odd-step clock: debt `2^Z`.
  * Grand orbits reduce to the same (THM-4581 6(b)).
  * This holds for offsets in `Z[1/3]`: `y` and `y + 1/5` never meet under `3x + 1`.
* The `T^(−1/2)` of THM-4581 is the rank-one first-return exponent (upper bound at sketch level). HEURISTIC: so is the orphan law's (HYP-9242 measures 0.42–0.63 at accessible sizes).
* The atlas's P12 "rank two" (2 and 3 independent) splits as clock + debt = 1 + 1 (ANALOGY). For maps on `Z_d` the debt rank is at most `d − 1`, so a larger debt rank needs a larger digit base, not merely more primes.
* **P1, corrected.** Coalescence laws are invariant under affine conjugacy, which rescales offsets. They depend on the sheet through which offsets are accessible. The pair `3x ± 1` (conjugate by `x ↦ −x`) has identical laws; `3x + 5` does not coalesce from offset 1.

**Battery (section 1).**
* Exact facts, with the classical identity named:
  * Mersenne "rivers" are level sets of the odd count;
  * a river's tuning offset is `{−o log_2 3}`;
  * jump tuning is automatic.
* Numerical laws:
  * river spans satisfy span/√K ≤ 5.9 for `K ≤ 12800` (the mean ratio drifts from 1.85 to 2.68, so `O(√K)` is not established);
  * staircase jump rate `~K^(−0.42)`;
  * river count `≈ 2.8–3.1 K^0.41` (local exponents 0.39, 0.39, 0.36);
  * the +1-barrier children satisfy the Mersenne lag-1 relation and coalesce in rivers (195 for `D < 4500`); per river the landing integer has the Kesten–Goldie tail `≈ 2.87/n` (`κ = 1`).
* Anchor census: of 60 periodic anchors whose first letter is 2, 43 are barriers to depth 40. Every one of the 94 anchors with first letter ≥ 3 merges at depth 1 in three steps.

## 1. The battery (quick tests of reframes on the Mersenne line and the barrier)

Data: `runcompress_20261007/mersenne_sigma_12800.txt`, which gives the odd count `o` and Terras length `σ_T` of `M_K = 2^K − 1` for `K ≤ 12800`.

### 1.1 Rivers are level sets of the odd count (FINITE-EXACT for K ≤ 12800; overlaps THM-4605 statement 8)

**The identity.** For `n` reaching 1, let the orbit's odd values be `x`. Then

    σ_T(n) = o(n) log_2 3 + log_2 n + ε(n),     ε(n) = Σ_(odd x > 1 in the orbit) log_2(1 + 1/(3x)) ∈ (0, 0.33)  (empirically: max 0.3256 for n ≤ 10^6, at n = 993; 0.3217 over M_K).

* This is the classical product formula `2^(σ_T) = 3^o n Π(1 + 1/(3x))` (Terras, Everett, Lagarias).
* It was checked to 6 decimals at `K = 101, 1001, 5001`.
* For `n = M_K`: `σ_T − K = o log_2 3 + ε + log_2(1 − 2^(−K))`, and the right side lies in an interval of length < 1.

**Consequences.**
* The partner key `(o, σ_T − K)` is a function of `o` alone, so rivers (partner classes) are exactly the level sets of `o(M_K)`. Checked for every `K ≥ 2` (audit F reproduced the data file byte for byte).
* A river's offset is `ε = {−o log_2 3}` up to `2^(−K)/ln 2`.
  * battery1.py reported a "spread 0.0995 inside a river". That is the `log_2(1 − 2^(−K))` term at `K = 3, 4`, not a real spread.
* Jumps between rivers satisfy `‖Δo log_2 3‖ ≤ max ε − min ε ≈ 0.32`; that bound is automatic. H1's median 0.063 (against 0.253 for random integers) reflects in addition the concentration of `ε` (interquartile range 0.16–0.27).
* The concurrent session opus-S21 proved the exact form of this level-set statement (THM-4605 statement 8: certificates with merge value ≥ 3 preserve (σ_T(x_K), o(M_K)), and conversely). Its HYP-9243 studies the clusters.

### 1.2 River spans relative to √K (NUMERICAL)

* `K` and `K − D` are partners only if the lag-D debt walk, started at `D`, returns to 0 within the `~7.6K` Terras steps the child has left. That takes `~D^2` steps.

  | max K in river | rivers | mean span/√K | max span/√K |
  |---|---|---|---|
  | [100, 400) | 13 | 1.85 | 5.00 |
  | [400, 1600) | 27 | 2.11 | 4.84 |
  | [1600, 6400) | 44 | 2.38 | 5.29 |
  | [6400, 12800] | 31 | 2.68 | 5.87 |

* The means with errors are 1.85 ± 0.45, 2.11 ± 0.24, 2.38 ± 0.21, 2.68 ± 0.25 (audit F). The rise is about 2σ, and the last window is right-censored at 12800. So the data support `span/√K ≤ 5.9` for `K ≤ 12800`, not `O(√K)`.
* Widest by span: `[10367, 10972]` (605 exponents, 360 members). Widest by span/√K and largest by membership: `[7763, 8298]` (532 members, 5.87√K).
* This is consistent with the diffusive picture of HYP-9243.

### 1.3 Staircase jumps and river births (NUMERICAL)

* **Jumps** (`o(M_K) ≠ o(M_(K−1))`) per unit `K`: 0.110, 0.098, 0.054, 0.038 on [100,400), …, [6400,12800]. That is `rate·√K = 1.7, 3.1, 3.4, 3.7`, a local exponent of about 0.42, not 1/2 yet.
* **River births** (orphans): `R(X)/X^0.41 = 3.09, 3.01, 2.92, 2.82` at `X = 400, 1600, 6400, 12800`; local exponents 0.39, 0.39, 0.36. All 106 births at `K ≥ 200` occur at odd `K` (audit F).
* **Orphan rate by the two-run length** `m = v_2(K − 1)` (`K ≥ 200`): 0.0137, 0.0235, 0.0127, 0.0203, 0.0204 for `m = 1, 2, 3, 4, ≥5`. Flat but borderline (χ² = 7.4 on 4 degrees of freedom, p ≈ 0.11; audit F) (H2).
* **Braid width**: the number of distinct rivers in windows of 50, 200, 1000 exponents averages 1.6, 2.7, 8.4, with maxima 4, 7, 18 (H6).

### 1.4 The +1 barrier: the children form their own Mersenne line (PROVED structure + FINITE-EXACT + NUMERICAL)

* **Level recursion.** The child `y* = 2/3^D − 1` of HYP-9240 sheds one factor 3 per odd step. At 3-adic level m its numerator obeys
  `a_(m−1) = (3^(m−1) + odd(a_m))/2`, with `a_D = 2 − 3^D` and landing `N_D = a_0`, reached at time `s_0`.
  * The landing debt is `⌈s_0/2⌉ ≥ D/2`.
  * `N_D ≥ 0` and `N_D + 1 < (3/2)^D` were already in HYP-9240.
* **Rivers.**
  * `y*_(D−1) = 3y*_D + 2` holds at equal times; it is the Mersenne lag-1 pair-chain state.
  * So the children coalesce in rivers, and `(N_D, s_0)` and the entry debt are constant on each river.
  * PROVED from THM-4581 3 (transfer by equidistribution of `3^(−D)` in `Z_2`): consecutive depths share their landing for a density-one set of D.
  * For `3 ≤ D < 4500`: 195 rivers. 3929 of 4496 consecutive depths share a landing. The largest rivers have 144, 126, 124, 114, 114 depths.
  * Records: `N = 880` on D ∈ [253, 268]; `732` on [815, 852]; `412` on {43, 44}.
* **Tail.**
  * The real recursion `b_(m−1) = 1/2 + (3/2)2^(−v)b_m` is a perpetuity with `E[A] = 1`, so `κ = 1`. Its Kesten–Goldie constant is `C = 0.5/0.1744 = 2.87`.
  * Per river, `P(N > 16, 64, 256) = 0.22, 0.067, 0.026`, against `2.87/x = 0.18, 0.045, 0.011`. That is generous at 64–256: 5 rivers above 256 against 2.2 expected (audit F). There are only 40 distinct `N_D` values but 195 distinct landings `(N_D, s_0)`.
  * Per depth the tail looked truncated (no `N_D > 1024` for `D < 4500`). That is a sample-size effect: one river counts once.
  * From a uniform start `c = 3^m(v+1) ∈ [1, 2·3^m]`, the exhaustive landing tail is close to the i.i.d. model, and its maximum is `≈ 2(3/2)^m` (landing_map.out).
* **Reading.**
  * The barrier is decided river by river.
  * Absorption needs the orbit of `N_D` to carry an odd-step excess of about `s_0/2 ≈ D`. For `N ≤ 880` the excess is at most 9.
  * An all-D proof would have to control the breaks of the children's rivers. That is the Mersenne coalescence problem transported to the rational line `2/3^D − 1`.

### 1.5 Census of periodic anchors (FINITE-EXACT, anchor_census.py)

* The run end `x` may shadow any periodic point `c` that is 1 mod 4: 154 such points on 82 cycles (U-words of sum ≤ 9, length ≤ 4).
* The depth-D child is then the exact rational pair `(c, (c+1)/3^D − 1)`.
* Over `D ≤ 40`, the outcomes are: absorption 296, re-anchoring 192, shift 4222, drift 1450.
* **First letter ≥ 3** (`c ≡ 5 mod 8`), 94 anchors: all merge at `D = 1` at time 3. This is the classical three-step merge of `x ≡ 5 (mod 8)` with `(x − 2)/3`. Its proof is three steps: (odd, odd), (even, even), then `u` even and `v` odd, landing on `(3c+1)/8` (audit F checked all 2048 test residues).
* **First letter 2** (`c ≡ 1 mod 8`), 60 anchors:
  * 43 never absorb for `D ≤ 40`. Among them are `+1`, `−7`, `−29/11`, `37/5`, `53/37`, … (The 192 ANCHOR outcomes are numeric meetings with unpaid debt, not absorptions.)
  * 17 merge at some larger `D`. For example `7/23` (word (2,3)) at `D = 5, 6`, and `7/503` (word (2,7)) at `D = 15–18, 25, 26, 35, 36`.
* So the first letter 2 is necessary for a barrier but not sufficient. The +1 barrier (HYP-9240, `D ≤ 8000`) is the extreme case.

### 1.6 The 3x − 1 diagnostic (atlas P1; NUMERICAL, sheet_diagnostic.py)

* Under `x ↦ −x`, the Mersenne line's conjugate for 3x − 1 is `2^K + 1`. Its residual class is the even `K` (run end `≡ 9 mod 16` after conjugation).
* Orphans in the residual class, `K ≤ 3000`: 75 for 3x − 1 against 77 for 3x + 1. Window fractions: 0.107/0.093, 0.040/0.043, 0.020/0.023.
* The 3x − 1 orbits end about evenly in its three cycles (1037/1096/866).
* So the orphan law agrees for the conjugate pair `3x ± 1`, as it must (`x ↦ −x` conjugates them, with offsets ±1). This does not establish general sheet-blindness, which fails (audit F; the obstruction lemma of HYP-9244).
* A first count over odd `K` gave "2 orphans". That was the wrong parity class, since odd `K` are the three-step class for `2^K + 1`.

## 2. The coalescence phase diagram (NUMERICAL; coalescence_phase.py, coalescence_visits.py)

**The simulated chain.** The exact pair chain `u = M v + e` of a map `T(x) = (m_i x + r_i)/d` (`x ≡ i mod d`), driven by fresh uniform digits:

    j uniform,   i = (M j + e) mod d,   M' = M m_i/m_j,   e' = (m_i e + r_i − (m_i/m_j) M r_j)/d,

with start `(1, 1)` (`y` and `y + 1`) and exact rationals.

**Diagnostics.**
* `q(T)`, the probability of no merge by `T`.
* The visit window: the mean number of debt visits to `M = 1` during `(T/4, T]`, per chain still unmerged at `T/4`. Pólya predicts a window growing like `√T`, flat (log-recurrent), or decaying to 0 for rank 1, 2, ≥ 3.
  * The first version quoted the mean over all samples, frozen at absorption. That saturates even in rank 1 (3x+1: 9.0, 11.9, 14.1, 15.2 at T = 1024 … 65536), so it was replaced (audit F).
* Values marked [F] are audit F's independent integer-orbit runs (`y` uniform mod `d^(T+64)`, 1500–5000 samples). They reproduce the author's chain runs within about 2σ.

| map | ρ | contracting | q(T) | visit window |
|---|---|---|---|---|
| `x+1` | 0 | yes | 0 by T = 16 | — |
| `3x+1` | 1 | yes | `√T q` = 2.48, 4.13, 6.44, 8.60, 10.35, 10.69, 10.75 (T = 16 … 65536) [F] | about 8–13 [F] |
| `Z_3` (1,2,4) | 1 | yes | `√T q` 1.59 → 2.50 (T = 4096), 2.11 ± 0.36 (T = 16384) [F] | — |
| `Z_3` (1,2,5) | 2 | yes | 0.643, 0.582, 0.530, 0.487, 0.451, 0.421, 0.391 (T = 16 … 65536, ±0.011) [F]; author's three runs 0.60–0.66 → 0.34–0.43 | 0.21–0.26, flat [F] |
| `Z_5` (1,2,3,1,1) | 2 | yes | 0.396 → 0.210 (T = 16384) [F] | log-recurrent |
| `Z_5` (1,2,4,7,1) | 2 | yes | 0.699 → 0.487 (T = 16384) [F] | log-recurrent |
| `Z_5` (1,2,3,7,1) | 3 | yes | 0.689, 0.6875, 0.686, 0.6855 (T = 1024 … 65536) [F]; author 0.654 | 0.117 → 0.001 [F] |
| `5x+1` | 1 | no | 0.9930 ± 0.0012 [F]; `1 − q_5 = 0.99167` | `≈ 1.08 √T` |
| `Z_3` (1,4,16) | 1 | no | 0.803 ± 0.010 [F] (author 0.77) | `0.544 √T` |
| `Z_3` (1,5,7) | 2 | no | 1.000 to T = 4096 | `~log T` (2.86 at T = 4096) |

* **What the last three rows show.** The debt walk can be recurrent and still not coalesce, because the offset expands.
* **What the rank-3 row shows.** A contracting map can still fail to coalesce, because the debt walk is transient.
* **What audit F's obstructions show.** A contracting, recurrent map can still fail to coalesce from a given offset, because the merge state is inaccessible: `3x + 5` from offset 1; `Z_3` (1,2,5) with `r = (0,7,14)`.
* So coalescence needs all three factors: accessibility, contraction and recurrence.
* **Hunt for counterexamples** (audit F): 12 more rank-2 maps (none transient), 6 rank-3 maps and a rank-4 map (`q` flat at 0.45–0.98, windows to 0), and tied multipliers behaving by their rank.
  * Weakly contracting rank-2 maps can look frozen at accessible `T` (`Z_3` (1,2,11): 0.757 → 0.735), but their windows stay flat. So the window, not `q`, is the discriminator.

### 2b. Far starts: a refuted prediction and its repair (NUMERICAL + HEURISTIC; far_starts.py)

* **Prediction (refuted).** From offset `E` (3x+1, start `(0, E)`), the offset would need `~c log E` debt excursions, each with tail `t^(−1/2)`, so the merge time would scale like `(log E)^2`.
* **Data** (2000 chain paths each):

  | `E` | `2^1+1` | `2^2+1` | `2^4+1` | `2^8+1` | `2^16+1` | `2^32+1` | `2^64+1` | `2^128+1` |
  |---|---|---|---|---|---|---|---|---|
  | median merge time | 93 | 340 | 580 | 820 | 1053 | 1417 | 2223 | 3399 |
  | median/(ln E)² | 93 | 131 | 72 | 27 | 8.6 | 2.9 | 1.1 | 0.43 |
  | `√T·q(T)` at T = 16000 | 10.9 | 14.7 | 17.0 | 18.9 | 21.1 | 22.3 | 27.8 | 31.4 |

* **Repair.**
  * Every step, not only those at zero debt, multiplies the normalized offset by 1/2 or 3/2 on a fair coin (THM-4581 1). So `log|f|` drifts at `½ ln(3/4) = −0.144` per step, and the offset dies in about `7 ln E` steps.
  * By then the debt has diffused to `~√(3.5 ln E)`, and the return to zero debt costs another time of order `ln E` (heavy-tailed).
  * So the merge time should be linear in `ln E`. The data are consistent with this (slope 26–36 per unit of `ln E` between `2^32` and `2^128`), and also with `~230 (ln E)^0.6`; they do not distinguish the two.
  * The tail constant grows only slowly, roughly like `a + b√(ln E)`.

## 3. THM-4606 — coalescence needs contraction (PROVED)

See the theorem file. The proof has four steps:

1. **Step table.** Off departures, each step multiplies `f = e/p^(max(k,0))` by 1/2 or `p/2` on a fresh fair coin, with additive error < 1. Departures (flips at `k = 0`) give `×1/2` for both coins. The `k`-skeleton is a simple random walk.
2. **Escape lemma.**
   * The coin walk `Σξ` has drift `½ ln(p/4) > 0` (Hoeffding).
   * Departures are visits of the SRW to 0. They exceed `εn + n_0` with probability `≤ Σ_(n>n_0) e^(−ε√n/2)`, by the return tail `C(2m,m)/4^m ≥ 1/(2√m)`.
   * Induction then gives `|f_n| ≥ 4e^(Λn/3)` forever from `|f_0| ≥ A_p(η)`, with probability `≥ 1 − η`.
3. **Growth lemma.** From `(0, e)`: run at `k = 0` with `×p/2`, depart to `k = −1`, run with `×p/2`, return. This gives `|e*| ≥ (p/4)|e| − 1`. It grows geometrically above `x* = 4/(p−4)`, which is `< 1` for `p ≥ 9`. For `p = 5, 7`, explicit paths lift `|e| ≤ 4` above `x*`.
4. **Integers.** The density of `n` merging with `n + e` by time `K` is `≤ q_p(e)`. The numerical bound 0.0087 is for `e = ±1` only: `q_5(3) ≥ 1/32` (audit E).

Checks (escape_check.out; independently re-derived by audit E, see audit_E/REPORT.md, corrections in MISTAKE-589):
* the growth lemma for odd `p ≤ 31`, `|e| ≤ 2000`;
* escape paths for `p = 5, 7`;
* the step table on 20,000 random states.

Shortest merges from `(0, 1)` (shortest_merge.out; its pruning at `|f| > 10^4` is not rigorous beyond depth about 14, but audit E's rigorously pruned search confirms every value):

| p | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| shortest merge time | 3 | 11 | 4 | 11 | 18 | 23 | 5 | 13 | none to depth 30 | 8 | 18 |

## 4. HYP-9244 — the Pólya trichotomy (OPEN)

* See the hypothesis file for the statement, the evidence table, the caveat and the refutation tests.
* The hypothesis now carries an accessibility clause. Its obstruction lemma (found by audit F) is PROVED: if a prime `ℓ ∤ d·Πm_i` and some `c` satisfy `(m_i − d)c + r_i ≡ 0 (mod ℓ)` for all `i`, then `e_n − (1 − M_n)c ≡ (Π m_(i_k)/d) e_0 (mod ℓ)`.
* The parts most worth proving next:
  1. **The expanding case for every `d`.** Normalize `F = |e|/(1 + M)`. Then `E[log` of the step multiplier `| past] = Λ − (`a convexity correction concentrated near `M = 1)`. What remains is to show that the debt walk's occupation of a strip around `log M = 0` is `o(n)`, with quantitative control.
  2. **Transience in one explicit rank-3 example**, such as `Z_5` (1,2,3,7,1). Audit F's equidistribution of the coupling state `(M mod d, e mod d)`, independent of the debt, is the natural route: a Lamperti-type Lyapunov function averaged over that state.
  3. **Whether the absence of local obstructions implies accessibility.**

## 5. Dictionary (DICTIONARY / ANALOGY)

* **P12 (rank two).**
  * Collatz has two independent units at the adic places; one is consumed by the clock, the other is the debt.
  * So Collatz coalescence is the borderline-recurrent case: rank-one debt and a heavy `T^(−1/2)` tail.
  * For maps on `Z_d` the debt rank is at most `d − 1`, so larger debt ranks need larger digit bases.
* **P1, corrected.**
  * Coalescence laws are invariant under affine conjugacy, which rescales offsets. So `3x + 1` and `3x − 1` (conjugate by `x ↦ −x`) cannot be distinguished by coalescence: the orphan check in §1.6 agrees.
  * Other sheets differ through accessibility: `3x + 5` from offset 1 never coalesces (HYP-9244 obstruction lemma).
  * The atlas's P1 concerns divergence arguments, and those concern all of `3x + k` at once. The obstruction shows that "sheet-blind" must be read modulo the rescaling of offsets.
* **THM-4581's Lyapunov weight** (HEURISTIC). `ρ(θ) = 1 − √(1 − (3/4)^θ)` is the 3/4 contraction made explicit. At `p ≥ 5` the same algebra with `p/4 > 1` suggests `|f|^(−θ)` as the escaping quantity; THM-4606 proves escape by the coin walk instead.
* **Orphans** (HEURISTIC). The orphan law's measured exponents (0.42–0.63, HYP-9242) sit near the rank-one first-return exponent 1/2. The upward drift with many children may be a several-pursuers effect inside rank one.

## 6. Reproduction

All runs were in `04-computation/experiments/reframes_20261007/`, with Python 3.10, standard library only.

| command | runtime | output |
|---|---|---|
| `python3 battery1.py`, `python3 battery2.py`, `python3 staircase.py` | seconds | — |
| `python3 anchor_census.py` | about 15 min | anchor_census.out; `python3 anchor_letters.py` gives the first-letter split |
| `python3 coalescence_phase.py NAME NSAMP TMAX SEED` | — | phase_*.out |
| `python3 coalescence_visits.py NAME NSAMP TMAX SEED` | up to about 10 min | visits_*.out (the 3x+1 and Z3_124 runs quoted earlier had no saved output; the table now uses audit F's runs for those rows) |
| `python3 far_starts.py 2000 16000 31`, `python3 sheet_diagnostic.py 3000` | about 15 min each | far_starts.out, sheet_diagnostic.out |
| audits | — | audit_E/REPORT.md, audit_F/REPORT.md |
| `python3 escape_check.py` | about 30 s | escape_check.out |
| `python3 qp_estimate.py 200000 12` | about 40 min | qp_estimate.out |
| `python3 shortest_merge.py` | — | shortest_merge.out |
| `python3 landing_map.py`, `landing_windows.py`, `landing_rivers.py`, `landing_growth.py` | 1–15 min | landing_*.out (section 1.4) |
