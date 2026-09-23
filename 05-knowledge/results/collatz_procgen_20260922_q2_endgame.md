# The Q2 endgame: one expensive link, a sharp budget, and where Diophantine input can and cannot enter

**Status.**
**PROVED** (hand proofs below; each is also machine-checked where marked):
(1) every 3-adic unit lies in the class of exactly one of the two primary hostile points `1`, `1/2`, and
the canonical escape `Psi` (that thread's optimal `k`-move escape) has exactly one expanding branch, the
`1/2`-transfer with `k >= 3`, of price `2c(k-1)/3 in (1,2)`; every `Psi` step consumes exactly `k`
ternary digits. The escapes with `k >= 3` use a minimal loop of length `k-1`, i.e. HYP-9122(`k-1`),
FINITE-EXACT for `k <= 6001`; every other statement below is conditional on this only where marked;
(2) the lift lemma and the exact recursion for repeated hostile landings, for every dyadic thread. By
exact computation, the 7 known threads and the 98-point backward census are all `Psi`-preimages of
`1` or `1/2`;
(3) the frontier-item-1 recursion is correct for its route but not optimal (it overstates the landing by
exactly a factor 2 whenever `floor(k log2 3) - floor((k-1) log2 3) = 1`);
(4) the budget theorem: along canonical chains the cumulative price is at most `exp(c* D)` with
`c* = ln(128/81)/4 = 0.114395`, where `D` counts the digits consumed by expensive links. The bound holds for
all threads, known or not, and is attained;
(5) the value criterion, the sharp integer-safety Lemma `S_b*`, the backward perturbation Theorem `P_b`,
and a finite hostility criterion for dyadic points below `3/2`;
(6) Theorem `F_b`: for every `i >= 3` the points `1/2 + 3^i/2^(K0(i-1)+2)` and `1/2 + 3^i/2^(K0(i-1)+1)`
lie in `Bad_inf`, so the backward exceptional set is infinite and accumulates at `1/2` (mirror of the
forward Theorem F);
(7) conditional theorems `Q2 <= X_min`, `Q2 <= HYP-9122 & X_T` and `Q2 <= HYP-9122 & X_Psi`, with all
thresholds handled;
(8) the exact 3-adic form of LTE; `dim_H Bad_Psi <= 0.6415` (upper bound by a Chernoff estimate).
**REFUTED:**
(i) "a chain is limited by the digits of `m`": `m = 829812238225934675`, a 37.6-digit integer, runs 15
consecutive expensive links consuming 60 digits, and `m = 962019445183081187` consumes 73 digits in
expensive links before descending;
(ii) "chains raise the value by at most about `m^0.104`": `6082250` reaches `6.57 > m^0.104 = 5.08`, and
the random-model exponent is `1/theta* = 0.1157`;
(iii) the budget "along any certified descent": it fails for arbitrary descending paths and holds for
canonical links.
**FINITE-EXACT:**
* every `2 <= m <= 10^18` prime to 3 has a multiplicative `Psi`-descent, by an exhaustive DFS over
  `Psi`-alive classes mod `3^38`. Hence **Q2 holds below `10^18`** (previous bound `7.87*10^17`);
* the complete scan of `m = 14 (mod 27)` up to `10^11`;
* samples near all 98 census threads up to `10^18`;
* the digit-exhausted family `(1+3^k w)/2`, `w <= 49`, `k <= 400`, with `m` up to `10^192`;
* 93 of the 97 non-integer census points PROVED hostile by complete exact certificates;
* exact counts of `Psi`-alive classes to `n = 200`.

**CITED:** Yu (2007), Stewart (1980), Senge–Straus (1973), Lagarias (2009); statements verified as
recorded in section 5.
**OPEN:** Q2; `X_min`; `X_Psi`; HYP-9122; whether `Bad_inf` is countable.
**UPDATE (orchestrator, 2026-09-23).** The four census points above `3/2` that §4 leaves undecided all
**descend**. The [Q1-mirror lane](collatz_procgen_20260922_q1_mirror.md) §6.2 gives certificates: a
first move `k_1`, then the greedy map. The orchestrator replayed them independently with exact
arithmetic: `793585/2^19` (`k_1 = 33`, `2^84 < 3^53`), `1675672562339/2^40` (`k_1 = 45`,
`2^103 < 3^65`), `419868489953/2^38` (`k_1 = 101`, `2^336 < 3^212`) and `1690578594467/2^40`
(`k_1 = 57`, `2^168 < 3^106`). Each path ends at the node `1`. Each uses a root child with
`k >= 9`, exactly as §4's reduction requires. So candidate **C5's conjecture that these points are
hostile is REFUTED**, and no hostile point above `3/2` is known on either side (Q1 mirror M3).

Session `collatz-procgen-20260922` (mac-mini), endgame lane, 2026-09-22/23. Scripts:
`04-computation/experiments/collatz_procgen_20260922_endgame_{chain.py, hostile.py, diophantine.py, psi.c, run.sh}`.
Output: [collatz_procgen_20260922_endgame.out](collatz_procgen_20260922_endgame.out). Reproduction: §9.

## 0. Inheritance and setting

* **Objects.** The graph `E` has arrows `n -> n/2` (`n` even) and `n -> 3n+1` (every `n`). The reverse
  move is `x -> (2^k x - 1)/3`; it is legal iff the result is a 3-adic unit, which for a positive integer
  means a positive integer prime to 3. A path of `s` moves with `K` halvings has multiplier `R = 2^K/3^s`
  and ends at `R(x - beta)`, where `beta = sum_t 3^(t-1)/2^(K_t)` (Lemma C_b). Q2 asks whether `1`
  reaches every `m` with `3` not dividing `m`. `Bad_inf` is the set of 3-adic units from which no legal
  path ever has `R < 1`.
* **Inherited and used (not re-proved).**
  * The Lean mod-27 lemma: every `m > 1` outside `1, 14 (mod 27)` descends in at most two moves
    ([choice ladder](collatz_procgen_20260922_choice_ladder.md) §5).
  * The loops lane
    ([loops and escapes](collatz_procgen_20260922_loops_and_escapes.md)):
    * Prop. 1.2: every loop through 1 of length `s >= 2` has `K >= K0(s) = floor((s+1) log2 3)`.
      This is unconditional and is used in Theorem F_b.
    * The 1-escape; Lemma 7.1 (1/2-transfer) and Lemma 7.2 (its sharpness); the prices `2c(k-1)/3`
      and `64c(k-4)/81`; the level-12 see-saw data; HYP-9122, FINITE-EXACT for `s <= 6000`.
  * The dimension lane's backward threads to `r = 41`, including the dump behind the 98-point
    census, and Lemmas S_b, C_b
    ([exceptional dimension](collatz_procgen_20260922_exceptional_dimension.md) §2, §4.7).
  * The greedy map `G` ([three-adic G map](collatz_mod6_20260917_three_adic_g_map.md), audited):
    * stationary drift `log(2/3)` per step and sharp large-deviation rate `0.758751`;
    * exceptional dimension `0.748`, and the sharp peak bound `K_n <= 2n + [k_1 = 3]`, i.e. growth at
      most `ln(4/3) = 0.2877` nats per digit;
    * hostile family `3^j + 1`; worst peak `133.03 m` at `m = 4847486`;
    * "every `m` reaches 1 under `G`" implies Q2.

    The wave-one note ([extended E-graph](collatz_mod6_20260917_extended_collatz_scc.md)) supplies the
    word-determination law and the density-1 stopping theorem.
* **Notation.**
  * `L = log2 3`; `K0(s) = floor((s+1)L)`.
  * `K*(0)=0`, `K*(1)=2`, and `K*(s) = K0(s)` for `s >= 2` (the least halving count of a loop of
    length `s`, assuming HYP-9122).
  * `c(s) = 2^K0(s)/3^s in (3/2, 3)`.
  * "Precision" of `x` at `h`: `v_3(x - h)`.

## 1. The canonical escape `Psi` and the single expensive link (PROVED)

**Lemma 1.1 (two primary hostile points).** Every 3-adic unit `x` satisfies exactly one of
`x = 1 (mod 3)` and `x = 2 = 1/2 (mod 3)`. Put `h(x) = 1` or `1/2` accordingly, and
`k(x) = v_3(x - h(x)) >= 1`. (Trivial. In base 3, the 1-class at precision `k` reads `...d 0^(k-1) 1` and
the 1/2-class reads `...d 1^(k-1) 2`, with `d` not continuing the run.)

**Definition (canonical escape).** `Psi(x) = rho_h(k) (x - h)`, where
`rho_1(k) = 2^(K*(k-1))/3^k` and `rho_(1/2)(k) = 2^(K*(k-1)+1)/3^k`.

**Proposition 1.2.** `Psi(x)` is the endpoint of an explicit legal path of exactly `k` moves.
* For `k = 1, 2` these are the Lean routes:
  * `k = 1`: move 0 for the 1-class, move 1 for the 1/2-class;
  * `k = 2`: moves `(2,0)` and `(3,0)`.
* For `k >= 3` the path is a minimal loop of length `k-1` through 1, with its first move increased by
  one for `h = 1/2`, followed by the move 0. This is the 1-escape, respectively the loops lane's
  Lemma 7.1; it needs HYP-9122(`k-1`).

Among all legal paths of length `<= k` it has the least endpoint (inherited Prop. 2.1 and Lemma 7.2;
re-checked by brute force on 160 random shells, `k = 3..6`). The multipliers are:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 12 | general `k >= 3` |
|---|---|---|---|---|---|---|---|---|---|
| `rho_1(k)` | 1/3 | 4/9 | 16/27 | 64/81 | 128/243 | 512/729 | 2048/2187 | 0.9865 | `c(k-1)/3 in (1/2,1)` |
| `rho_(1/2)(k)` | 2/3 | 8/9 | 32/27 | **128/81** | 256/243 | 1024/729 | 4096/2187 | 1.9731 | `2c(k-1)/3 in (1,2)` |
| `ln rho_(1/2)(k) / k` | -0.405 | -0.059 | 0.0566 | **0.1144** | 0.0104 | 0.0566 | 0.0896 | 0.0566 | `< ln2/k` |

**Proposition 1.3 (one expensive link).** `rho < 1` on every branch except `(h, k) = (1/2, k >= 3)`.
There `1 < rho < 2`. So **the 1/2-transfer at precision `>= 3` (an "L1 link") is the only step of
`Psi` that raises the value.** *Proof.* `c(s) in (3/2,3)` for `s >= 2`. ∎

**Lemma 1.4 (digit consumption).** On each shell `S_(h,k) = {x : v_3(x-h) = k}`, the map `Psi` is a
bijection onto the 3-adic units, and it multiplies distances by `3^k`. So `Psi(x) mod 3^n` is
determined by `x mod 3^(n+k)`, and it is not determined by fewer digits. *Proof.* `x - h` runs over
`3^k Z_3^x` and `rho = 2^a/3^k`. ∎ (Checked on all residues, `k <= 5`, `n <= 4`.) Under Haar measure the
branches of successive `Psi` steps are i.i.d., with `P(h,k) = 3^-k` (§3.4).

## 2. Exact chain calculus (task 1; PROVED, machine-checked)

**Lemma 2.1 (Psi-lift).** Let `x` be a 3-adic unit with `Psi`-itinerary
`(h_1,k_1), ..., (h_T,k_T)`, `D_T = k_1 + ... + k_T` and `P_T = rho_1 ... rho_T = 2^(A_T)/3^(D_T)`. If
`v_3(m - x) >= D_T + 1`, then `m` has the same first `T` steps and
`Psi^T(m) - Psi^T(x) = P_T (m - x)`. In particular, if `Psi^T(x) = tau in {1, 1/2}` and
`v_3(m - x) = j`, then `Psi^T(m) = tau + P_T(m - x)` lies on the thread of `tau` at precision exactly
`j - D_T`.

*Proof.* The branch `(h, k)` of a unit `y` is determined by `y mod 3^(k+1)`. If `v_3(m_t - x_t) >= k_t + 1`,
the two share the branch, and `m_(t+1) - x_(t+1) = rho_t (m_t - x_t)` loses exactly `k_t` in valuation.
Induct. ∎

**Theorem 2.2 (universal landing recursion).** Write `m_t = h_t + 3^(k_t) u_t`, where `u_t` is a 3-adic
unit in `Z[1/2]` (`u = w/2` on the 1/2 thread). Then `m_(t+1) = Psi(m_t) = 2^(a_t) u_t` with
`a_t = K*(k_t-1) + [h_t = 1/2]`, and

    u_(t+1) = (2^(a_t) u_t - h_(t+1)) / 3^(k_(t+1)),   k_(t+1) = v_3(2^(a_t) u_t - h_(t+1)).

On the 1/2 thread this reads `w_(t+1) = (2^(K0(k_t - 1)+1) w_t - 1)/3^(k_(t+1))`. It was checked on
3300 consecutive steps of random big-integer orbits.

**Corollary 2.3 (every dyadic thread).** Let `h = p/2^e` be a dyadic hostile point with
`Psi^T(h) = tau`, digits `D_h`, multiplier `P_h = 2^(A_h)/3^(D_h)`. Let `m = (p + 3^j W)/2^e` with
`3` not dividing `W` and `j > D_h`. Then:
* **the link:** `m` follows `h`'s itinerary for `D_h` digits. It then makes the `tau`-escape at
  precision `j - D_h`, landing at `y' = 2^(A_h + a_tau(j - D_h)) W / 2^e`.
* **price:** `P_h rho_tau(j - D_h) = 2^(A_h + a_tau(j-D_h))/3^j`, as a 3-adic multiplier. The actual
  ratio `y'/m` is smaller by the factor `1 - beta/m`.
* **consumption:** exactly `j` ternary digits.
* **repeated landing:** if `y'` is at thread `h'' = p''/2^(e'')` with precision `j''`, then

      W'' = (2^(A_h + a_tau(j - D_h) + e'' - e) W - p'') / 3^(j'').

Verified on 6984 random integers near all 98 census threads (`chain.py` §F).

**Thread table (PROVED membership: §4; itineraries exact).**

| thread `h` | value | itinerary | `D_h` | `P_h` | `tau` | link price at precision `j` |
|---|---|---|---|---|---|---|
| `1/2` | 0.500 | base | 0 | 1 | 1/2 | `2c(j-1)/3` |
| `43/32` | 1.344 | `1/2@3` | 3 | 32/27 | 1 | `(32/27) c(j-4)/3` (descends iff `c(j-4) < 81/32`) |
| `59/64` | 0.922 | `1/2@3` | 3 | 32/27 | 1/2 | `64 c(j-4)/81` |
| `145/128` | 1.133 | `1/2@4` | 4 | 128/81 | 1 | `(128/81) c(j-5)/3` |
| `209/256` | 0.816 | `1/2@4` | 4 | 128/81 | 1/2 | `(256/243) c(j-5)` |
| `371/256` | 1.449 | `1/2@5` | 5 | 256/243 | 1 | `(256/729) c(j-6)` |
| `499/512` | 0.975 | `1/2@5` | 5 | 256/243 | 1/2 | `(512/729) c(j-6)` |

* **Census.** All 98 census points reach `1` (63 points) or `1/2` (35 points) exactly, in at most 3 `Psi`
  steps and at most 26 digits. Every prefix multiplier is `>= 1.039`. **Every step of every census
  itinerary is a 1/2-class step.** The shapes are:

  | steps | ending at 1 | ending at 1/2 | example |
  |---|---|---|---|
  | 1 | 23 | 22 | `1055729/2^20`, `1580017/2^21` |
  | 2 | 37 | 12 | `9958534033/2^33`, `4754357/2^22` |
  | 3 | 2 | 0 | `1675672562339/2^40` |

  Later steps may be the descending branches `1/2@1`, `1/2@2`, e.g. `661/512 -> 5/4 -> 1/2` and
  `T_12^(-1)(2) -> 2 -> 1`.
* **Per-digit price.** The maximum of `ln(link price)/j` over all threads and precisions is `0.114395`.
* **Relation to the loops lane.** Its prices `2c(k-1)/3` and `64c(k-4)/81` are the first two rows. Its
  routes `(3,2,0)` for `43/32` and `59/64` are the single `Psi` step `1/2@3`.

**Proposition 2.4 (the frontier recursion).** Frontier item 1 of the synthesis states
`w_(t+1) = (2^(K_t+3) w_t - 1)/3^(j_(t+1))`, with `K_t` the halving count of a 1-loop of length `j_t - 2`.
This recursion is **correct for the route** `1/2 --(move 3)--> 1 --(loop of length k-2)--> 1 --(0)--> 2^(K_t+2) w`.
Its landing is `2^(K*(k-2)+2) w`, whereas the optimal transfer lands at `2^(K*(k-1)) w`. The two agree
iff `K0(k-1) = K0(k-2) + 2`. Otherwise the frontier route lands **twice as high**; this happens for 32 of
the 78 values `k = 4..81`, a fraction tending to `2 - log2 3 = 0.415`. The optimal recursion replaces
`K_t + 3` by `K0(k_t - 1) + 1`. (Executed move by move on random integers, `k <= 81`.)

## 3. The budget theorem and the exact endgame condition (task 2)

**Theorem 3.1 (budget; PROVED, sharp).** Let `m_0` be a positive integer, and assume HYP-9122 for the
precisions `k_t - 1` used (FINITE-EXACT for `k_t <= 6001`). Then for every `T`,

    m_T / m_0  <=  prod_(t<T) rho_t  <=  exp(c* D_L1(T)),   c* = ln(128/81)/4 = 0.1143953,

where `D_L1(T)` is the number of digits consumed by the L1 links among the first `T` steps. The constant
is optimal, and it is attained exactly (up to the factors `1 - h_t/m_t`) by runs of links with `k = 4`,
i.e. `v_3(2 m_t - 1) = 4`.

*Proof.* `m_(t+1) = rho_t (m_t - h_t) < rho_t m_t`. For `k >= 8`, `ln rho_(1/2)(k)/k < ln 2/8 < c*`. The
table of §1 covers `k <= 7`. Non-L1 steps have `rho < 1`. ∎

The price of a single link has no floor above 1: `rho_(1/2)(k) = 2c(k-1)/3` equals `1.0535` at `k = 5`,
`1.0115` at `k = 41` and `1.0010` at `k = 306`. Its infimum is 1, approached at the records of
`log2 3`. The synthesis's floor `32/27` came from the route through `1`, which the loops lane already
superseded. What is bounded is the price per digit.

Because every hostile thread's link factors through `Psi` (Cor. 2.3), **the constant covers all threads,
known or unknown.** The loops lane's `0.114` was computed on the 1/2 thread alone; it is the universal
constant. The example `m = 60091390742` attains `15.5722 = (128/81)^6` with six consecutive `k = 4`
links.

**What the budget does not give (REFUTED statements).**

1. *Digits do not run out.* A run of L1 links does not consume the digits of `m_0`. Each link replaces
   `m` by `2^a u` with `u = (2m-1)/3^k`, and the value grows by the price. The next precision is bounded
   only by `log_3(2 m_t)`, not by what is left of `m_0`. Examples of L1 digits exceeding the ternary
   length of `m_0`:
   * `m_0 = 19296859109`: 21.5 digits, 38 digits in L1 links (scan to `10^11`);
   * `m_0 = 829812238225934675`: 37.6 digits, a single run of 15 consecutive L1 links consuming 60
     digits;
   * `m_0 = 962019445183081187`: 37.7 digits, 73 digits in L1 links (DFS to `10^18`; orbits in the
     `.out`).
2. *"At most about `m^0.104`" is false.* `m = 6082250` reaches excursion `6.5695 > m^0.104 = 5.08`.
   `m = 909072835151420597` reaches `90.93 > m^0.104 = 73.7`. Under Haar measure the log-price is a
   random walk with drift `-0.4114` per digit and Lundberg exponent `theta* = 8.6434`
   (`E[rho^theta*] = 1`). So the random-model size of the worst excursion over `m <= X` is
   `X^(1/theta*) = X^0.1157`. It predicts `7.3, 18.7, 121` at `X = 3*10^7, 10^11, 10^18`; the observed
   values are `6.57, 15.57, 90.93`. Mixed chains, which carry entropy, beat pure `k = 4` chains, which
   is why `0.1157 > 0.104`.
3. *Not along arbitrary descents.* The bound is a property of the canonical links. A single legal move
   `2i+1` from the 1/2 class costs `2^(2i+1)/3` for one digit. The half-chain lane's shortest descending
   paths reach `959 m` within 16 moves, against `exp(16 c*) = 6.2`.

**Proposition 3.2 (what a chain leaves at its end; PROVED).** Take a maximal run of L1 links
`m_0 -> ... -> m_T`. Then:
* `m_T = 2^(K0(k-1)) w_final`, where `k = k_(T-1)` and `w_final = (2 m_(T-1) - 1)/3^k` is odd and prime
  to 3;
* either `m_T = 1 (mod 3)`, or `m_T = 2 (mod 3)` with `v_3(2m_T - 1) <= 2`. In both cases the next
  `Psi` step descends;
* `w_final = 1` (a pure power of 2) iff `m_(T-1) = (3^k + 1)/2`, i.e. iff the last link consumed all
  remaining digits. More generally, `w_final` is small iff `k ~ log_3(2 m_(T-1))`.

**The exact endgame condition (3-adic).** For `theta > 0` put

    Bad(theta) = { xi in Z_3^x : every legal reverse path from xi has all prefix multipliers >= theta },

a closed set, with `Bad(1) = Bad_inf`. Then the following hold.

Let `P_T = prod_(t<T) rho_t = 2^(A_T)/3^(D_T)` be the 3-adic multiplier of the run (note
`m_T = P_T(m_0 - beta_T) < P_T m_0`).

* (E1) `m_0 in Bad_inf` implies `m_T in Bad(1/P_T)` for the end point of every canonical run (exact,
  3-adic).
* (E2) If `m_T` is not in `Bad(1/P_T)`, i.e. some legal path from `m_T` has multiplier `< 1/P_T`, then
  `m_0` has a multiplicative descent through its run. This is the certificate form; E1 is its
  contrapositive for hostile points.
* (E3) In the digit-exhausted case the end point is `2^(floor(k log2 3)) w` with `w` small. The
  condition then says that **the 3-adic numbers `2^(floor(k log2 3)) w` avoid
  `Bad(3^k/2^(floor(k log2 3)+1))`**. This is the exact sense in which "the endgame is the ternary
  digits of `2^K`". It is a transversality statement: the Beatty exponents `floor(k log2 3)` against a
  thin 3-adic set.

**§3.4 The renewal picture and the canonical-strategy endgame (PROVED upper bound, FINITE-EXACT counts).**
Under Haar measure, branch `(h,k)` has mass `3^-k` and consumes `k` digits. This gives:
* `E[digits] = 3/2`, `E[ln rho] = -0.617109`, so the drift is **`-0.411406` nats per digit**.
* The `Psi`-exceptional set `Bad_Psi = {xi : prod_(t<T) rho_t(xi) >= 1 for all T}` satisfies
  `dim_H Bad_Psi <= min_theta s(theta) = 0.64150`, where `sum_b 3^(-s k_b) rho_b^theta = 1`.
  *Proof of the upper bound.* An alive class mod `3^(n+1)` is one of 2 residue classes attached to an
  alive branch word of total length `D <= n`. For `theta >= 0`,
  `#{words: D <= n, P >= 1} <= 3^(sn) sum_w 3^(-s D(w)) P(w)^theta = 3^(sn)/(1 - Z(s,theta))`.
  This is finite when `Z < 1`. ∎

  The matching lower bound is the standard tilted-measure argument, not written out. The exact counts
  of alive classes (renewal DP, confirmed at `n <= 19` by the independent class DFS in C) give
  `log_3(count)/n = 0.509, 0.538, 0.552, 0.571, 0.595, 0.614` at `n = 20, 30, 40, 60, 100, 200`, still
  rising toward `0.6415`.
* The Haar probability of surviving `D` digits without a `Psi`-certificate is `3^(-0.3585 D) = 0.6745^D`.

| strategy | typical drift (nats/digit) | survival rate per digit | exceptional dimension | adversarial growth (nats/digit) | worst excursion found |
|---|---|---|---|---|---|
| greedy `G` (inherited) | `ln(2/3) = -0.4055` | `0.758751` | `0.748` | `ln(4/3) = 0.2877` (sharp) | `133.03` at `m <= 10^7` |
| canonical `Psi` (this note) | `-0.4114` | `0.6745` | `<= 0.6415` (counts: `0.614` at `n=200`) | `c* = 0.1144` (sharp) | `6.57` (`<= 3*10^7`), `15.57` (`<= 10^11`), `90.93` (`<= 10^18`) |
| full choice `E` (dimension lane) | – | – | evidence: 0 (countable, `~r^1.7` classes) | – | – |

`Psi` resolves the inherited `G`-hard cases at once:
* `m = 4847486` (`G`-peak `133.03 m`): 2 steps, excursion `1.87`;
* the family `3^j + 1` (`G` climbs by `(4/3)^(j-1)`): one descending step.

## 4. The hostile set near integers (task 4)

**Lemma V (value criterion; PROVED).** Let `x > 0`. A node `z > 0` of `x`'s reverse tree with `z >= x`
has multiplier `R > 1`. Hence a prefix with `R <= 1` ends at a node of real value `< x`. *Proof.*
`z = R(x - beta)` with `beta > 0`, so `x - beta = z/R > 0` and `R = z/(x-beta) > z/x >= 1`. ∎

**Lemma S_b\* (sharp safety at integer nodes; PROVED).** Let `y >= 1` be an integer node with multiplier
`R`.
* **(S1)** If `R > max(1, 3y/4)`, every continuation keeps the total multiplier `> 1`.
* **(S2)** If `3R >= 2y + 1`, every continuation keeps the total multiplier `> 1`.

*Proof.* A continuation of `s' >= 1` moves ends at an integer `z >= 1` with `2^(K')y = 3^(s')z + B'` and
`B' >= (3^(s')-1)/2 >= 3^(s'-1)`.
* (S1): `R' >= (z + 1/3)/y >= 4/(3y)`.
* (S2), case `z >= 2`: `RR' >= 7R/(3y) > 1`.
* (S2), case `z = 1`: every move satisfies `(2^k v - 1)/3 >= (v-1)/3`, so `3^(s') >= (2y+1)/3`. Then
  `R' >= (3 - 3^(-s'))/(2y) >= 3/(2y+1)` and `RR' >= 1`. Equality is impossible because
  `2^K != 3^s`. ∎

At a node of start `x`, `y = R(x - beta)`, so (S1) holds when `x - beta < 4/3` and (S2) when
`R(3 - 2(x - beta)) >= 1`. In particular, **for `x <= 4/3` every positive integer node is either a
descent or SAFE.**

**Theorem P_b (backward perturbation; PROVED).** Let `h = p/2^e > 0` be in `Bad_inf`, with a tree of
positive nodes. Let `x = h + 3^i c/2^j > 0`, where `i >= 1`, `j >= 0` and `c != 0`. Assume:
* **(B_b)** every legal path of length `i-1` from `h` makes at least `max(e, j)` halvings;
* **(A_b)** for every path from `h`, at the first time `tau` with `K_tau >= max(e, j)`, the integer
  `y = x_tau(h) + 2^(K_tau - j) 3^(i - tau) c` is `>= 1` and satisfies (S1) or (S2).

Then `x in Bad_inf`.

*Proof.* Paths of length `<= i-1` from `x` are paths from `h`, since the residues agree mod
`3^(i-t) >= 9`, and they carry the same multipliers, `> 1`. By (B_b) each path reaches, at some
`tau <= i-1`, the integer node `y`. That node is SAFE by (A_b) and Lemma S_b\*. ∎

*What it generates.* The strongest choice is `h = 1/2`, `c = 1`, `j = K0(i-1) + 1`, which is exactly
`T_i^(-1)(1)` (Theorem F_b). Condition (B_b) caps `j` at `g_(i-1)(h)`, so every perturbation has real
size `>= 3/M_h(i-1)`, with `M_h` the least multiplier over paths of length `i-1`. For `h = 1/2` this size
exceeds `1/2`. Unlike the forward Theorem P, which has the slack `1/2 - delta`, **the backward
perturbations have no slack below the value `3/2`**. Most census points are not perturbations: they
are `Psi`-preimages, see below.

**Theorem F_b (an infinite family; PROVED for all `i >= 3`, machine-checked for `i <= 22`).** Let
`a = K0(i-1)`, so `2^a < 3^i < 2^(a+1)`. Then

    x_i^+ = 1/2 + 3^i/2^(a+2) = T_i^(-1)(1/2)  (value 1/2 + 3/(4c(i-1)) in (3/4, 1))
    x_i^- = 1/2 + 3^i/2^(a+1) = T_i^(-1)(1)    (value 1/2 + 3/(2c(i-1)) in (1, 3/2))

both lie in `Bad_inf`. The pairs are `(59/64, 43/32)`, `(209/256, 145/128)` and `(499/512, 371/256)` for
`i = 3, 4, 5`. Since `v_3(x_i^+- - 1/2) = i`, **`Bad_inf` is infinite and `1/2` is an accumulation
point**: the backward mirror of the forward Theorem F.

*Proof.*
* **Depth `s <= i-1`.** A node is `N + 3^(i-s) 2^(K_s - j)`, with `N >= 1` the node of `1/2` along the
  same moves (the paths are shared). Its multiplier is `2^(K_s)/3^s = 2(N + B_s/3^s) >= 8/3`.
* **Depth `i`.** A node is `z = M + 2^(K_i - j)` with `M = (2^(k_i) N - 1)/3 >= 0` an integer. From
  `2^(K_(i-1)) >= 3^(i-1)(2N+1) - 1` one gets:
  * `z` is non-integral only when `k_i = 0` and `N = 1`, i.e. the path from `1/2` reached `1` at depth
    `i-1`. Then `K_(i-1) >= K0(i-1)+1` by the unconditional loop bound;
  * the cases `N = 4` with `k_i = 0`, `N = 2` with `k_i = 1`, and `k_i >= 2` give
    `2^(K_i) >= 3^(i+1) - 1`, `(10/3)3^i - 2` and `4 * 3^i - 4` respectively. Each exceeds
    `2 * 3^i > 2^(a+1)`, so `K_i >= a+2`: the node is an integer (larger `N` only increase `K_i`).

  So for `x^+` the only non-integral depth-`i` node is `1/2`, reached with multiplier `2c(i-1)/3 > 1`,
  whose subtree is that of the hostile point `1/2`. For `x^-` there is none.
* **Integer nodes.** For `x^+ < 1` every positive integer node is SAFE (Lemma C_b with S1). For `x^-`,
  every path of length `i-1` from `1/2` has `2^K >= 3^i - 1 > 2^a`, since `3^i - 1` is not a power of 2
  for `i >= 3`. So every path from `x^-` is integral by depth `i-1`: `x^-` is exactly Theorem P_b with
  `h = 1/2`, `c = 1`, `j = a+1`. (S2) holds at every integer node at depth `<= i`:
  * at the first integer node of a path, at depth `s <= i-1`, it reduces to
    `N(4c-6) + 6(c-1)B_s/3^s - c >= 0`. This is true when a halving follows the first move
    (`B_s/3^s > 1/2`, so the left side is `>= 3(2c-3) > 0`), and otherwise it reduces to
    `2(2^(a+1) - 3^i) >= 1`;
  * the same computation at depth `i` gives `2^(a+1) - 3^i >= c - 1`, which is true because
    `2^(a+1) - 3^i` is odd and `>= 3`. It is only needed as a check, since for `x^-` every path is
    already integral by depth `i-1`.
* No node is negative. ∎

**Criterion H (a finite hostility test; PROVED).** For a dyadic `x < 3/2` every node satisfies
`z/R = x - beta < 3/2`, i.e. `3R > 2z`. So the integer children `y_k` of any node satisfy (S2) for all but
finitely many `k`. On unsafe integer nodes `R < 1/(2 beta_1)` holds and `beta` grows by `> 2 beta_1/3`
per step, so unsafe paths have bounded length. The search in `endgame_hostile.py` is therefore a
terminating decision procedure for dyadic `x < 3/2`; negative nodes, which never occurred, would need a
separate descent. Results:
* **93 of the 97 non-integer census points are PROVED hostile.** Their values lie in `[0.5, 1.4622]`,
  and they are all the census points below `3/2`.
* **The four census points above `3/2` are undecided.** They are `793585/2^19`,
  `419868489953/2^38`, `1675672562339/2^40` and `1690578594467/2^40`, with values `1.514–1.538`.
  * All four are `Psi`-preimages of the **integer 2**: itinerary `1/2@k -> 2 -> 1/2@1 -> 1`, e.g.
    `793585/2^19 = 1/2 + 3^12/2^19 = T_12^(-1)(2)` and `419868489953/2^38 = T_24^(-1)(2)`.
  * They sit at the shells `i = 12` and `i = 24`, where `2^19/3^12 = 0.9865` and `2^38/3^24 = 0.9731`
    are close to 1 (the lower convergent `19/12` of `log2 3` and its double). There `c(i-1) = 2.96, 2.92`
    is close to 3, and the value `1/2 + 3/c(i-1)` is just above `3/2`. This explains the dimension
    lane's "near-coincidence" of numerators: `793585 = 3^12 + 2^18` is also the numerator of the forward
    generation-1 point `-1 - 2^18/3^12`.
  * Their Bad_41 membership excludes certificates of length `<= 41`. Capped searches (depth `<= 20`,
    `4*10^5` nodes each) found none.
  * The obstruction was located exactly. For `793585/2^19`, `419868489953/2^38` and
    `1690578594467/2^40`, an exhaustive search of the rest of the tree proves every node safe. The only
    exceptions are the root's children `y_k = (2^k x - 1)/3` with `2^k(2x - 3) >= 2`: moves `k >= 9`
    (`k >= 5` for the last point), integral or not; the listing was checked up to `k = 60`. For each of
    them `3R <= 2z` (the **3/2 wall**). For the fourth point the search stopped at `2*10^6` nodes.
  * Hence `x in Bad_inf` iff no such child has a continuation of multiplier `< 3/2^k`. For an integer
    child (`k >= e`) this means a path from `y_k` to 1 whose carry ratio `B'/3^(s')` is
    `< x - 1 - 2^(-k)`, i.e. within `0.014–0.038` of the absolute minimum `1/2`, which only the climb
    `(3^J - 1)/2 -> ... -> 1` attains.
* **The candidates of `Psi`-preimage type.** Consider all `Psi`-preimages of `{1, 1/2}` with every prefix
  multiplier `> 1`, value in `[0.3, 2]` and denominator `<= 2^24`: there are 842.
  * All 35 census points with `e <= 24` are among them.
  * Of the other 807, 805 carry explicit descents (choice beyond `Psi`); 2 are undecided (values
    `1.61`, `1.62`).
  * Descending candidates lie mostly above `3/2` (607) or in `(1, 3/2]` (195).
  * The census points lie in `[0.5, 1.5376]`.

**The integer question: test and obstruction.**
* *Test (FINITE-EXACT).* No integer `2 <= m <= 10^18` lies in `Bad_inf`, because each has a
  multiplicative `Psi`-descent (§7).
* *Why the suggested route cannot close.*
  * (a) "Points of `Bad_inf` are dyadic rationals and their limits" is vacuous. `Z[1/2]` is dense in
    `Z_3`, and the positive integers are themselves dyadic (`e = 0`).
  * (b) The usable structure statement is `Bad_inf subset Z[1/2] cap (0,2)` (**S_D**): every point,
    including every accumulation point, is a dyadic rational of real value `< 2`. S_D implies
    `Bad_inf cap Z_>=2 = {}` trivially. But its `e = 0` instance *is* that statement, and nothing in
    Lemmas V, S_b\*, C_b or Theorem P_b reduces the `e = 0` instance to `e >= 1`.
  * (c) The obstruction is precise. At an integer `m >= 2` the root is an unsafe integer node:
    (S1) needs `1 > 3m/4` and (S2) needs `3 >= 2m+1`. Its children `y_k = (2^k m - 1)/3` have safety
    ratio `R_k/y_k -> 1/m <= 1/2 < 2/3`. So **infinitely many unsafe branches** carry no finite
    certificate either way, and the budget method can neither certify hostility nor exclude it.
  * Integers sit beyond the **3/2 wall**. This is the same wall that leaves the four census points
    undecided, where their fate is decided by exceptionally efficient paths to 1, i.e. near-solutions
    of `2^A p - 2^B ~ 3^C` (an archimedean linear form in `log 2, log 3, log p`).

## 5. Diophantine inputs (task 3)

The four inputs, with the citations checked this session:

**(a) Lifting the exponent (PROVED; elementary).** For `3` not dividing `w`, put `n0 = [w = 2 mod 3]`,
`w' = 2^(n0) w = 1 (mod 3)` and `alpha_w = -log(w')/log 4 in Z_3` (3-adic logarithm; `v_3(log 4) = 1`).
Then

    v_3(2^n w - 1) = [n = n0 mod 2] * (1 + v_3((n - n0)/2 - alpha_w)).

(Checked in 79,800 cases, `w < 400`, `n < 600`.)
* If `w` is a power of 2, `alpha_w` is a non-positive integer and the depth is `1 + v_3(M + a)`, which is
  `<= 1 + log_3 n`. This is the pure-power endgame: the first landing of `2^(K0(k-1))` is predicted
  exactly (398 of 398 cases in the family `(1+3^k)/2`).
* Otherwise `alpha_w` is a 3-adic irrational. The depth is the length of a coincidence between `n/2`
  and the digits of `alpha_w`. Over `w < 400` the excess over the trivial `log_3 n` is at most 4 at
  `n < 3^20` (maximum at `w = 341 = (4^5-1)/3`, a trunk number).
* **Role:** it controls the pure-power case completely, and for general `w` it reduces the depth to a
  3-adic approximation property of `alpha_w`, which needs (b).

**(b) `p`-adic linear forms (Kunrui Yu, "p-adic logarithmic forms and group varieties III", Forum Math.
19 (2007) 187–280, doi:10.1515/FORUM.2007.009; CITED).** Two rationals, `p = 3`, as quoted in arXiv:2107.00971, Thm A.2:
`v_3((x1/y1)^b1 (x2/y2)^b2 - 1) < (16e)^6 2^(3/2) (log 4)^2 (3/(log 3)^2) (log A1)(log A2) max(log T, delta B/B_m)`.
The constant is `9.1*10^10`. Applied to `2^n w - p = p(2^(n-r)(2^r w/p) - 1)`, with `r` chosen so that
`3` does not divide `n - r`:

    v_3(2^n w - p) <= 9.2*10^10 * log(4 max(w,p,e)) * (log n + O(1)).

* **Role:** after a digit-exhausting link, whose remainder `w` has bounded height, the next hostile
  landing near any point `p/2^e` of bounded height has depth `O(log n log H)`. Its price is therefore at
  most `exp(c* O(log n log H))` (Theorem 3.1).
* It does **not** control the landings after that. Their remainders `(2^n w - 1)/3^k` are integers of
  size `2^n`, of unbounded height. It also gives no uniform bound near hostile points of unbounded
  height, which accumulate (Theorem F_b).
* Numerically it beats the trivial bound `n log_3 2` only for `n > ~10^14`.

**(c) Senge–Straus ("PV-numbers and sets of multiplicity", Period. Math. Hungar. 3 (1973) 93–100,
doi:10.1007/BF02018464) and Stewart ("On the representation of an integer in two different bases",
J. reine angew. Math. 319 (1980) 63–72, doi:10.1515/crll.1980.319.63); CITED.** Senge–Straus: for multiplicatively independent `a, b` only finitely many
integers have both digit sums bounded. Stewart's effective form (his lecture notes, Thm 12) states: for
`0 <= alpha < a` and `0 <= beta < b`, the number `L` of base-`a` digits `!= alpha` plus base-`b` digits
`!= beta` satisfies `L > log log n/(log log log n + C) - 1` for `n > 25`.

*Proposition 5.1 (no double exhaustion; PROVED from Stewart).* Suppose `n = 2^a u = (1 + 3^k w')/2`: a
digit-exhausted landing which is again a digit-exhausted 1/2-point. Then the ternary digits of `n` are
`1^(k-1) 2` preceded by those of `(w'-1)/2`, so `L_(1,3)(n) <= 2 + log_3 w'`, and
`L_(0,2)(n) <= 1 + log_2 u`. Stewart with `(a, alpha, b, beta) = (3, 1, 2, 0)` gives

    log_2 u + log_3 w' >= log log n / (log log log n + C) - 4.

The 1-thread version (`n = 1 + 3^k w'`) is the same with `alpha = 0`. **Consecutive digit-exhausting
links with bounded remainders happen only finitely often.** Numerically, with `u, w' <= 6000` and
`a <= 300`, there are 719 (1/2 thread) and 830 (1 thread) solutions, all with `n <= 3.7*10^7` and
`a <= 13`.
* **Role:** it caps the "adversary uses the last digits twice" scenario. Like (b), it says nothing about
  generic remainders.

**(d) Erdős's ternary problem and Lagarias, "Ternary expansions of powers of 2" (J. London Math. Soc.
79 (2009) 562–588, doi:10.1112/jlms/jdn080; arXiv math/0512006, Thm 1.5 and Conj. B, read in the
text).**
`E(Z_3) = {lambda : infinitely many lambda 2^n omit the digit 2}`, with `dim E^(1) = log_3 2`,
`dim E^(2) <= 1/2`, and Conjecture B: `dim E(Z_3) = 0`.
* **Role:** it is the same *type* of statement as the canonical-strategy endgame `X_Psi` (§6). Both ask
  whether an integer orbit under `x2`-type maps avoids a Cantor set of positive dimension (`log_3 2` for
  `Sigma_3`, `<= 0.6415` for `Bad_Psi`). But neither statement implies the other: the sets and the maps
  differ.
* With full choice the endgame is not of Erdős type. If `Bad_inf` is a countable set of dyadics, as the
  census and Theorem F_b suggest, integers avoid it by the structure statement S_D. No ternary-digit
  transversality is needed.

**Verdict.**
* The controlling inputs of the *digit-exhausted sub-case* are (b), with (a) as its exact special case,
  and (c): they bound the first post-exhaustion landing and forbid repeated exhaustion.
* **No Diophantine input controls the endgame as a whole.** After one non-exhausting link the remainders
  are generic. The remaining obstruction is structural (the shape of `Bad_inf`) or, for a fixed
  strategy, of Collatz/Erdős transversality type (`dim Bad_Psi > 0`).
* One place where Diophantine approximation genuinely enters the structure: the points just above `3/2`,
  via the convergents of `log2 3` (§4).

## 6. Conditional theorems (tasks 3 and 4; PROVED implications)

**Theorem 6.1.** `Q2` follows from
`X_min`: *no integer `m >= 2` lies in `Bad_inf`* (equivalently, every `m >= 2` prime to 3 has a legal
reverse path with `2^K < 3^s`).

*Proof.* For `m >= 2`, the path ends at `y = (2^K m - B)/3^s < (2^K/3^s) m < m`. By Lemma S_b, `y` is a
positive integer prime to 3. Induct on `m` from the base `m = 1`.

*Thresholds.* None. The backward carry `B > 0` only subtracts, so a multiplicative certificate gives
descent for every positive member. No HYP-9122 is needed.

*Status of `X_min`.* OPEN; FINITE-EXACT for `m <= 10^18` (this note, §7).

**Theorem 6.2.** Assume HYP-9122. Then `X_min` holds iff **S_1/2**: *no integer `m >= 14` with
`m = 14 (mod 27)` lies in `Bad_inf`*. `S_1/2` in turn follows from **`X_T`**: *for every `k >= 3` and
every odd `w` prime to 3, the integer `2^(floor(k log2 3)) w` has a legal reverse path with multiplier
`< 3^k/2^(floor(k log2 3)+1)`.* So `Q2 <= HYP-9122 & X_T`.

*Proof.* Take `m >= 2` prime to 3.
* **(i)** `m` not `1, 14 (mod 27)`: the Lean lemma gives descent by a factor `<= 8/9` within two moves
  (`m > 1` is needed).
* **(ii)** `m = 1 (mod 27)`, `k = v_3(m-1) >= 3`: the 1-escape uses a loop of length `k-1` with `K0(k-1)`
  halvings (HYP-9122(`k-1`); FINITE-EXACT when `k <= 6001`, i.e. for all `m < 3^6002`). It gives
  `y = 2^(K0(k-1))(m-1)/3^k < m`, since `2^(K0(k-1)) < 3^k`.
* **(iii)** `m = 14 (mod 27)`, `k = v_3(2m-1) >= 3`: the transfer gives
  `y = 2^(K0(k-1)) w = (2^(K0(k-1)+1)/3^k)(m - 1/2)`. By `X_T` some path from `y` has multiplier
  `rho' < 3^k/2^(K0(k-1)+1)`. Its endpoint `z < rho' y < m - 1/2` is a positive integer prime to 3, and
  the total multiplier is `< 1`.

Induction as in 6.1. Every certificate is multiplicative, so there are no thresholds. ∎

Under HYP-9122 the chain of implications is `X_T => S_1/2 <=> X_min => Q2`.

**Theorem 6.3 (canonical strategy).** Assume HYP-9122. Then `Q2` follows from
`X_Psi`: *no integer `m >= 2` lies in `Bad_Psi`* (every `Psi`-orbit reaches a prefix multiplier `< 1`).
Moreover `X_Psi => X_min`.

*Proof.* Every `Psi`-prefix is a legal path; apply Theorem 6.1. ∎

* **Status.** OPEN; FINITE-EXACT for `m <= 10^18`.
* **Nature.** `X_Psi` asks the integers to avoid a closed set of dimension `<= 0.6415` defined by a
  deterministic digit-shift map. It is a statement of Collatz/Erdős type (§5(d)).

**Relation to "no divergent positive `G`-orbit"** (the inherited open problem).
* `X_G`: *every positive `G`-orbit reaches 1*. It implies Q2 by reversing the orbit (inherited
  Theorem 6.1 of the G-note).
* `X_G` and `X_Psi` are both "no non-descending orbit" statements for deterministic inverse maps. They
  differ only on the two hostile threads:
  * `Psi` replaces the greedy climb (factor `4^(k-1)/3^k`) by the optimal escape (factors `c(k-1)/3` and
    `2c(k-1)/3`), which uses choice through the loops;
  * on all other classes the two maps take the same word, the greedy Lean routes.
* Neither statement implies the other.
* `X_Psi` is the more plausible of the two, though it is not logically weaker. Its exceptional set is
  thinner (`<= 0.6415` against `0.748`), its adversarial growth rate is smaller (`0.1144` against
  `0.2877` per digit), and its survival rate is smaller (`0.6745` against `0.7588` per digit).
* `X_min` is weaker than both. It uses all choice, and its exceptional set is conjecturally countable.
* The weakest statement is Q2 itself (`X_min => Q2`, not conversely). Q2 also allows
  non-multiplicative (threshold) descents. Such a descent from `m` ends at `y = R(m - beta) < m` with
  every prefix multiplier `> 1`. Hence `beta < s/3`, so `R - 1 < s/(3m - s)`: the multiplier
  `2^K/3^s` must lie within about `s/(3m)` of 1. Only very long paths can do this for large `m`.

## 7. Numerics (task 5; FINITE-EXACT, two independent code paths)

C (`endgame_psi.c`, 128-bit) and Python (`chain.py`) give identical histograms on
`m = 14 (mod 27) <= 3*10^6`. The C class-DFS and the Python renewal DP give identical alive-class counts
(24, 274, 3050, 31438 at `n = 7, 11, 15, 19`). In every run `ln(m_t/m_0) - c* D_L1(t) <= 0`: the budget
held with maximal excess `-1*10^-6`.

**(a) The whole 1/2 neighbourhood to `10^11`** (every `m = 14 (mod 27)`, 3,703,703,704 orbits; 3.5 min):

* All orbits descend. The mean number of `Psi` steps is 2.2848; the maximum is 19, at
  `m = 26426067953`.
* L1 links before descent: `1: 3.46e9, 2: 2.27e8, 3: 1.65e7, 4: 1.38e6, 5: 135448, 6: 13630, 7: 1468,
  8: 166, 9: 10, 10: 1`.
* Longest consecutive L1 run: **8 links, 28 digits**, price `12.30`, at `m = 23853568379`.
* Most L1 digits: 38, at `m = 19296859109` (21.5 digits).
* **Largest excursion `15.5722` at `m = 60091390742`**, which equals the budget `exp(24 c*)` exactly.

**(b) Every `m <= 10^18`** (exhaustive DFS over the `Psi`-alive classes mod `3^38`; `1.28*10^10` nodes;
29 min; 1 MB):

* `3,498,679,527` alive classes have a least representative `m <= 10^18`. Every one of them descends,
  and in every case the first descent is multiplicative (0 exceptions).
* All other classes die by a multiplicative `Psi`-certificate within their first 38 digits.
* Together with strong induction: **Q2 holds for every `m <= 10^18`.**
* Largest number of digits consumed before descent: **87**, at `m = 5728828425613736` (33 ternary
  digits; 30 `Psi` steps, printed in the `.out`).
* **Largest excursion `90.93`** at `m = 909072835151420597`, after 11 steps of which 10 are L1 links.
* **Longest chain: 15 consecutive L1 links consuming 60 digits** (price `59.81`) at
  `m = 829812238225934675`, a 37.6-digit integer.
* Most L1 links in one orbit: 17, at `m = 466685813006909240`. Most L1 digits: 73, at
  `m = 962019445183081187`.
* Digits consumed until descent: the histogram falls from `1.50*10^9` orbits at `D = 38` to 1 orbit at
  `D = 87`.
* An earlier run of the same DFS flagged 69 "non-multiplicative" leaves. They were exactly the leaves
  with `D >= 80`, where a 128-bit power table was silently skipped; the check now uses
  `2^A < 3^D <=> A <= floor(D log2 3)` and finds none.

**(c) Near the other hostile threads, up to `10^18`** (for all 98 census threads, 2000 random `W` per
thread and per precision `j <= 37`, `m = (p + 3^j W)/2^e <= 10^18`; 7,216,965 orbits; 1.5 s):
* all descend;
* largest excursion `10.94`, near the thread `6283/2^13`;
* at most 9 L1 links; the longest run has 7 links, with cumulative price `7.29` near `4754357/2^22`;
* the budget excess is `<= 0` in every orbit.

Near every thread the first steps follow the thread's itinerary, as Lemma 2.1 requires. The excursions
stay far below the worst ones of the whole range (b): the hostile threads are not where the long chains
live, because a chain needs fresh 1/2-landings, which the itinerary does not supply.

**(d) The digit-exhausted family** `m = (1 + 3^k w)/2`, odd `w <= 49` prime to 3, `3 <= k <= 400` (6766
members, `m` up to `10^192`):
* all descend under `Psi`;
* worst excursion `5.09` (`w = 49`, `k = 184`); at most 7 steps;
* the precision of the post-exhaustion landing `2^(K0(k-1)) w` is `<= 12` throughout (LTE and Yu predict
  `O(log k)`).

**(e) Comparison with the budget.**
* Observed worst excursions: `6.57`, `15.57`, `90.93` for `m <= 3*10^7, 10^11, 10^18`.
* Budget along the same orbits: `exp(c* D_L1)` = `8.79`, `15.57`, and `>= 171` for 45 L1 digits. The
  bound is attained at `10^11`.
* In terms of `m` there is no budget. The random model predicts `X^0.1157` (7.3, 18.7, 121), consistent
  with the data.
* With (b), the `Psi`-exceptional integers are rare but they exist at every scale. The number of classes
  that survive `n` digits grows like `3^(0.6 n)`.

## 8. Hypothesis candidates (not filed as HYP files)

* **C1 (`Psi`-preimage structure).** `Bad_inf cap Z[1/2]` consists exactly of those `Psi`-preimages of
  `{1, 1/2}` whose prefix multipliers all exceed `1` and whose trees pass Criterion H. The real values of
  its points lie in `[1/2, 3/2 + 0.04]`.
  * Support: census `subset` candidates for `e <= 40`; 805 of 807 non-census candidates descend
    (`e <= 24`); 93/97 census points proved.
* **C2 (value wall).** No point of `Bad_inf` has real value `>= 1.6`. The census maximum is `1.5376`,
  attained at `Psi`-preimages of 2 at the convergent shells of `log2 3`.
* **C3 (`X_Psi`).** Every integer `m >= 2` has a multiplicative `Psi`-descent. FINITE-EXACT to `10^18`.
  Heuristic support: the survival probability decays as `0.6745^D` while the integer's digits are
  exhausted after `log_3 m`, but `dim Bad_Psi` is positive.
* **C4 (worst excursion).** `max_(m <= X)` of the `Psi`-excursion is `X^(0.1157 + o(1))`. The observed
  values are within a factor 1.3 of the prediction at `X = 3*10^7, 10^11, 10^18`.
* **C5 (the points above 3/2).** `T_12^(-1)(2) = 793585/2^19` lies in `Bad_inf` iff no root child
  `y_k = (2^k x - 1)/3` with `k >= 9` has a continuation of multiplier `< 3/2^k`. For `k >= 19` this
  means a path from `y_k` to 1 whose carry ratio is `< x - 1 - 2^(-k)`. The equivalence is PROVED by
  exhaustive search of the rest of the tree (the root's children with `k < 9` and all their
  descendants); the same holds for `T_24^(-1)(2)`. We conjecture that these points are hostile (Bad_41; no descent found). If so, hostility
  above `3/2` is decided by near-optimal climbs, i.e. by near-solutions of `2^A p - 2^B ~ 3^C`.
  **REFUTED (2026-09-23):** all four descend. The first moves are `k_1 = 33, 45, 101, 57`, and the
  greedy continuations reach multipliers `2^84/3^53`, `2^103/3^65`, `2^336/3^212` and `2^168/3^106`,
  all `< 1` (see the status update). The near-optimal climbs the conjecture anticipated do exist.

## 9. Reproduction

```bash
bash 04-computation/experiments/collatz_procgen_20260922_endgame_run.sh > 05-knowledge/results/collatz_procgen_20260922_endgame.out
FULL=1 bash 04-computation/experiments/collatz_procgen_20260922_endgame_run.sh   # + scan to 1e11, DFS to 1e18 (~35 min)
```

The quick mode took 90 s on one core (mac-mini, other lanes running), with peak memory `< 300 MB`. The
slow part is the exact prover on the census (about 55 s). Set `DUMP=<bwd_bad_r41.txt>` to recompute the 98-point census from
the dimension lane's dump (the census is otherwise stored as data in `chain.py`). The `.out` file
concatenates the quick run and the FULL-mode logs (scan to `10^11`, DFS to `10^18`), with their timings.

| program | content |
|---|---|
| `endgame_chain.py` | `Psi`, `c*`, explicit escapes with loop words (`s <= 80`), digit consumption, optimality inside `k` moves, frontier check, universal recursion, thread table, lift lemma on 6984 integers |
| `endgame_hostile.py` | exact prover (Lemmas V, S_b\*, criterion H), census, Theorem F_b families and `T^(-1)(2)` points (`i <= 22`), `Psi`-candidates, capped searches above `3/2` |
| `endgame_diophantine.py` | LTE in 3-adic form, `alpha_w` digits, Yu's constant, double exhaustion, the digit-exhausted family to `10^192`, renewal drift, `dim Bad_Psi`, Lundberg exponent, exact alive counts to `n = 200` |
| `endgame_psi.c` | 128-bit `Psi` orbits: scan, samples near threads, exhaustive alive-class DFS, single orbits |
