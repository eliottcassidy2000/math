# Coalescence is contraction times recurrence: the debt-lattice rank decides when affinely related Collatz-type orbits merge (THM-4606, HYP-9244), with a battery of quick reframe tests on the Mersenne line, the +1 barrier and periodic anchors

Session `mac-mini-2026-10-08-reframes`, 2026-10-07/08.

Brief from the owner: creative hypothesis generation and testing; look through past results for inspiration; find reframes that unlock Collatz-adjacent proofs.

Scripts are in `04-computation/experiments/reframes_20261007/`.

## 0. The answer in one screen

**The reframe.**
* THM-4581 (two 2-adic orbits related by `u = 3^k v + e` merge almost surely, with tail `T^(-1/2)`) has two separate inputs:
  1. **Contraction of the offset.** At zero debt, the difference `e = u − v` is multiplied by 1/2 or 3/2 whenever the parities agree. Its geometric mean is `√3/2 < 1`.
  2. **Recurrence of the debt walk.** The debt `k` is a one-dimensional simple random walk in flip time.
* Replace 3 by another multiplier, or Collatz by a generalized map on `Z_d`, and the two inputs come apart. Each can fail on its own.

**Proved (THM-4606).**
* For `px + 1` with `p ≥ 5`, the orbits of Haar `y` and `y + e` meet with probability `q_p(e) < 1`, and `q_p(e) → 0` as `|e| → ∞`.
* For `p = 1, 3` they meet almost surely.
* So Haar coalescence among the `px+1` maps holds exactly at the contracting multipliers. It is sheet-blind (`px + r` is conjugate to `px + 1` on `Z_2`).
* Numerically:
  * `q_5 = 0.0083`, `q_7 = 0.108`, `q_9 = 0.0023`, `q_11 ≈ 2.5e-5`, `q_13 ≈ 5e-6`;
  * the spikes at `p = 2^a − 1` come from a reset path of length `a + 1`.

**Conjectured, with all five regimes observed (HYP-9244).**
* For generalized (Matthews–Watts) maps, coalescence holds iff the map contracts and the debt walk on the multiplicative group `Γ = ⟨m_i/m_j⟩` is recurrent.
* For non-degenerate maps the rank of `Γ` decides it, as in Pólya's theorem:

  | rank | coalescence | tail of P(no merge by T) |
  |---|---|---|
  | 0 | almost sure | exponential |
  | 1 | almost sure | `T^(-1/2)` |
  | 2 | almost sure | `1/log T` |
  | ≥ 3 | fails with positive probability, even for contracting maps | — |

* Expanding maps fail at every rank.
* Caveat (Peres–Popov–Sousi): the debt walk's step law is chosen adaptively, and adaptive choice can change the type. The conjecture is that this arithmetic rule does not.

**Readings.**
* **Collatz is rank one in every clock.**
  * Terras clock: debt `3^Z`, with the 2-powers spent on time.
  * Odd-step clock: debt `2^Z`.
  * Grand orbits reduce to the same (THM-4581 6(b)).
* The `T^(−1/2)` of THM-4581 and the `1/2` of the orphan law (HYP-9242) are the rank-one first-return exponent.
* The atlas's P12 "rank two" (2 and 3 independent) splits as clock + debt = 1 + 1. One more independent unit would slow coalescence to `1/log T`; two more would break it.
* Coalescence-based arguments are sheet-blind (P1).

**Battery (section 1).**
* Exact facts, with the classical identity named:
  * Mersenne "rivers" are level sets of the odd count;
  * a river's tuning offset is `{−o log_2 3}`;
  * jump tuning is automatic.
* Numerical laws:
  * river spans are `O(√K)`, with span/√K at most 5.9;
  * staircase jump rate `~K^(−0.42)`;
  * river count `≈ 2.9 K^0.41`;
  * the +1-barrier landing integers have a Kesten–Goldie tail `P(N_D > n) ≈ 3/n`, with exponent `κ = 1`.
* Anchor census: of 60 periodic anchors whose first letter is 2, 43 are barriers to depth 40. Every one of the 94 anchors with first letter ≥ 3 merges at depth 1 in three steps.

## 1. The battery (quick tests of reframes on the Mersenne line and the barrier)

Data: `runcompress_20261007/mersenne_sigma_12800.txt`, which gives the odd count `o` and Terras length `σ_T` of `M_K = 2^K − 1` for `K ≤ 12800`.

### 1.1 Rivers are level sets of the odd count (exact; overlaps THM-4605 statement 8)

**The identity.** For `n` reaching 1, let the orbit's odd values be `x`. Then

    σ_T(n) = o(n) log_2 3 + log_2 n + ε(n),     ε(n) = Σ_(odd x > 1 in the orbit) log_2(1 + 1/(3x)) ∈ (0, 0.33).

* This is the classical product formula `2^(σ_T) = 3^o n Π(1 + 1/(3x))` (Terras, Everett, Lagarias).
* It was checked to 6 decimals at `K = 101, 1001, 5001`.
* For `n = M_K`: `σ_T − K = o log_2 3 + ε + log_2(1 − 2^(−K))`, and the right side lies in an interval of length < 1.

**Consequences.**
* The partner key `(o, σ_T − K)` is a function of `o` alone, so rivers (partner classes) are exactly the level sets of `o(M_K)`. Checked for every `K ≥ 5`.
* A river's offset is `ε = {−o log_2 3}` up to `2^(−K)`.
  * battery1.py reported a "spread 0.0995 inside a river". That is the `log_2(1 − 2^(−K))` term at `K = 3, 4`, not a real spread.
* Jumps between rivers satisfy `‖Δo log_2 3‖ ≤ max ε − min ε`. H1's "tuning of jumps" (median 0.063 against 0.253 for random integers) is therefore automatic.
* The concurrent session opus-S21 proved the exact form of this level-set statement (THM-4605 statement 8: certificates with merge value ≥ 3 preserve (σ_T(x_K), o(M_K)), and conversely). Its HYP-9243 studies the clusters.

### 1.2 River spans are O(√K) (NUMERICAL)

* `K` and `K − D` are partners only if the lag-D debt walk, started at `D`, returns to 0 within the `~7.6K` Terras steps the child has left. That takes `~D^2` steps.

  | max K in river | rivers | mean span/√K | max span/√K |
  |---|---|---|---|
  | [100, 400) | 13 | 1.85 | 5.00 |
  | [400, 1600) | 27 | 2.11 | 4.84 |
  | [1600, 6400) | 44 | 2.38 | 5.29 |
  | [6400, 12800] | 31 | 2.68 | 5.87 |

* Widest: `[7763, 8298]`, 532 members, 5.87√K.
* This is the diffusive picture of HYP-9243 seen from the span.

### 1.3 Staircase jumps and river births (NUMERICAL)

* **Jumps** (`o(M_K) ≠ o(M_(K−1))`) per unit `K`: 0.110, 0.098, 0.054, 0.038 on [100,400), …, [6400,12800]. That is `rate·√K = 1.7, 3.1, 3.4, 3.7`, a local exponent of about 0.42, not 1/2 yet.
* **River births** (orphans): `R(X)/X^0.41 = 3.09, 3.01, 2.92, 2.82` at `X = 400, 1600, 6400, 12800`.
* **Orphan rate by the two-run length** `m = v_2(K − 1)` (`K ≥ 200`): 0.0137, 0.0235, 0.0127, 0.0203, 0.0204 for `m = 1, 2, 3, 4, ≥5`. Flat within noise (H2).
* **Braid width**: the number of distinct rivers in windows of 50, 200, 1000 exponents averages 1.6, 2.7, 8.4, with maxima 4, 7, 18 (H6).

### 1.4 The +1 barrier: the child lands on a small integer (Kesten–Goldie tail; HEURISTIC + NUMERICAL)

* The child `y* = 2/3^D − 1` of HYP-9240 sheds one factor 3 per odd step:
  * `a/3^m` is odd ↦ `(a + 3^(m−1))/(2·3^(m−1))`.
  * So it becomes an integer `N_D` after exactly `D` odd steps, at time `s_0`.
  * The landing debt is `⌈s_0/2⌉ ≥ D/2`.
* **Elementary bound:** `N_D < (3/2)^D`, from the maximal carry `2^(s_0−D)(3^D − 2^D)`.
* **Data, `3 ≤ D ≤ 4000`:**
  * every `N_D` is positive; max 880; median 4;
  * `s_0/D ∈ [1.5, 2.88]`;
  * `n·P(N_D > n) = 1.5, 1.7, 2.8, 3.5, 3.6, 3.6, 2.5, 2.8, 3.8` for `n = 2, …, 512`.
* **Reading.**
  * Along the landing, the real recursion `v ↦ v/2` or `(3v+1)/2` with fair parities is a perpetuity. `E[A^κ] = ½((1/2)^κ + (3/2)^κ) = 1` holds at `κ = 1`, so `P(N > n) ~ C/n` (Kesten 1973, Goldie 1991).
  * Absorption needs the orbit of `N_D` to have odd excess about `D`. For `N ≤ 880` the excess is at most 9.
  * Paying a debt `≥ D/2` needs `N_D ≳ e^(cD)`, which has probability `~e^(−cD)`: summable.
  * This is a quantitative heuristic for HYP-9240. It is not a proof: a proof would need the orbits of all `N_D`.

### 1.5 Census of periodic anchors (FINITE-EXACT, anchor_census.py)

* The run end `x` may shadow any periodic point `c` that is 1 mod 4: 154 such points on 82 cycles (U-words of sum ≤ 9, length ≤ 4).
* The depth-D child is then the exact rational pair `(c, (c+1)/3^D − 1)`.
* Over `D ≤ 40`, the outcomes are: absorption 296, re-anchoring 192, shift 4222, drift 1450.
* **First letter ≥ 3** (`c ≡ 5 mod 8`), 94 anchors: all merge at `D = 1` at time 3. This is the classical three-step merge of `x ≡ 5 (mod 8)` with `(x − 2)/3`. Its proof: two common odd steps, a common even step, then `u` even and `v` odd land on `(3c+1)/8`.
* **First letter 2** (`c ≡ 1 mod 8`), 60 anchors:
  * 43 never merge for `D ≤ 40`. Among them are `+1`, `−7`, `−29/11`, `37/5`, `53/37`, …
  * 17 merge at some larger `D`. For example `7/23` (word (2,3)) at `D = 5, 6`, and `7/503` (word (2,7)) at `D = 15–18, 25, 26, 35, 36`.
* So the first letter 2 is necessary for a barrier but not sufficient. The +1 barrier (HYP-9240, `D ≤ 8000`) is the extreme case.

## 2. The coalescence phase diagram (NUMERICAL; coalescence_phase.py, coalescence_visits.py)

**The simulated chain.** The exact pair chain `u = M v + e` of a map `T(x) = (m_i x + r_i)/d` (`x ≡ i mod d`), driven by fresh uniform digits:

    j uniform,   i = (M j + e) mod d,   M' = M m_i/m_j,   e' = (m_i e + r_i − (m_i/m_j) M r_j)/d,

with start `(1, 1)` (`y` and `y + 1`) and exact rationals.

**Diagnostics.**
* `q(T)`, the probability of no merge by `T`.
* The mean number of visits of the debt to `M = 1`. Pólya predicts `√T`, `log T` or bounded for rank 1, 2, ≥ 3.

| map | ρ | contracting | q(T) | visits |
|---|---|---|---|---|
| `x+1` | 0 | yes | 0 by T = 16 | — |
| `3x+1` | 1 | yes | `√T q`: 1.6, 2.4, 3.9, 6.0, 7.8, 9.6 (T = 4..4096) | `√T` |
| `Z_3` (1,2,4) | 1 | yes | `√T q`: 1.5, 2.0, 2.4, 2.6, 2.75, 2.69 (T = 16..16384) | `√T` |
| `Z_3` (1,2,5) | 2 | yes | 0.598, 0.524, 0.480, 0.442, 0.394, 0.366, 0.338 (T = 16..65536) | `+0.22..0.34` per ×4 among unmerged chains |
| `Z_5` (1,2,3,1,1) | 2 | yes | 0.423 → 0.243 (T = 16384) | log |
| `Z_5` (1,2,4,7,1) | 2 | yes | 0.692 → 0.478 (T = 16384) | log |
| `Z_5` (1,2,3,7,1) | 3 | yes | 0.674, 0.664, 0.660, 0.658, 0.656, 0.654, 0.654 (T = 16..65536) | 0.47 → 0.54, saturated |
| `5x+1` | 1 | no | → 0.9917 by T = 64 | — |
| `Z_3` (1,4,16) | 1 | no | 0.77 from T = 16 on | `0.54 √T` |
| `Z_3` (1,5,7) | 2 | no | 1.000 to T = 4096 | 2.70 at T = 4096, `~log T` |

* **What the last two rows show.** The debt walk can be recurrent and still not coalesce, because the offset expands.
* **What the rank-3 row shows.** A contracting map can still fail to coalesce, because the debt walk is transient.
* So coalescence needs both factors.

## 3. THM-4606 — coalescence needs contraction (PROVED)

See the theorem file. The proof has four steps:

1. **Step table.** Off departures, each step multiplies `f = e/p^(max(k,0))` by 1/2 or `p/2` on a fresh fair coin, with additive error < 1. Departures (flips at `k = 0`) give `×1/2` for both coins. The `k`-skeleton is a simple random walk.
2. **Escape lemma.**
   * The coin walk `Σξ` has drift `½ ln(p/4) > 0` (Hoeffding).
   * Departures are visits of the SRW to 0. They exceed `εn + n_0` with probability `≤ Σ_(n>n_0) e^(−ε√n/2)`, by the return tail `C(2m,m)/4^m ≥ 1/(2√m)`.
   * Induction then gives `|f_n| ≥ 4e^(Λn/3)` forever from `|f_0| ≥ A_p(η)`, with probability `≥ 1 − η`.
3. **Growth lemma.** From `(0, e)`: run at `k = 0` with `×p/2`, depart to `k = −1`, run with `×p/2`, return. This gives `|e*| ≥ (p/4)|e| − 1`. It grows geometrically above `x* = 4/(p−4)`, which is `< 1` for `p ≥ 9`. For `p = 5, 7`, explicit paths lift `|e| ≤ 4` above `x*`.
4. **Integers.** The density of `n` merging with `n + e` by time `K` is `≤ q_p(e)`.

Checks (escape_check.out):
* the growth lemma for odd `p ≤ 31`, `|e| ≤ 2000`;
* escape paths for `p = 5, 7`;
* the step table on 20,000 random states.

Shortest merges from `(0, 1)` (shortest_merge.out):

| p | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| shortest merge time | 3 | 11 | 4 | 11 | 18 | 23 | 5 | 13 | none to depth 26 | 8 | 18 |

## 4. HYP-9244 — the Pólya trichotomy (OPEN)

* See the hypothesis file for the statement, the evidence table, the caveat and the refutation tests.
* The two parts most worth proving next:
  1. The expanding case for every `d`. The obstacle is the number of times the debt crosses `M = 1`. For `d = 2` it is the SRW skeleton; for general `d`, a mean-zero walk chosen adaptively.
  2. Transience in one explicit rank-3 example, such as `Z_5` (1,2,3,7,1). This needs a Lamperti-type Lyapunov function averaged over the modulating state `(M, e) mod 5`. Peres–Popov–Sousi show that some averaging is unavoidable.

## 5. Dictionary (DICTIONARY / ANALOGY)

* **P12 (rank two).**
  * Collatz has two independent units at the adic places; one is consumed by the clock, the other is the debt.
  * So Collatz coalescence is the borderline-recurrent case: rank-one debt, heavy `T^(−1/2)` tail, orphans with density `~L^(−1/2)`.
  * A "Collatz with a third unit" would have rank-two debt. A fourth would make orphans a positive proportion.
* **P1 (sheet-blindness).** Coalescence depends only on the multipliers, so it can never prove a statement that distinguishes `3x + 1` from `3x − 1`.
* **THM-4581's Lyapunov weight.** `ρ(θ) = 1 − √(1 − (3/4)^θ)` is the 3/4 contraction made explicit. At `p ≥ 5` the same algebra with `p/4 > 1` produces the escaping quantity `|f|^(−θ)`.
* **Orphans.** The orphan law's exponent 1/2 is the rank-one first-return exponent. Its upward drift with many children (0.54–0.63, HYP-9242) is a several-pursuers effect inside rank one, not a change of rank.

## 6. Reproduction

All runs were in `04-computation/experiments/reframes_20261007/`, with Python 3.10, standard library only.

| command | runtime | output |
|---|---|---|
| `python3 battery1.py`, `python3 battery2.py`, `python3 staircase.py` | seconds | — |
| `python3 anchor_census.py` | about 15 min | anchor_census.out; `python3 anchor_letters.py` gives the first-letter split |
| `python3 coalescence_phase.py NAME NSAMP TMAX SEED` | — | phase_*.out |
| `python3 coalescence_visits.py NAME NSAMP TMAX SEED` | up to about 10 min | visits_*.out |
| `python3 escape_check.py` | about 30 s | escape_check.out |
| `python3 qp_estimate.py 200000 12` | about 40 min | qp_estimate.out |
| `python3 shortest_merge.py` | — | shortest_merge.out |
