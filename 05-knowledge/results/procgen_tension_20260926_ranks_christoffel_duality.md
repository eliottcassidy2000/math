# Tensions are certificates: periodic ranks, rank defects and flows on the strategy cube; the Christoffel claims refuted, exact maximal densities, and the negation census

**Status.**
- **PROVED** (hand proofs below; every lemma is also exercised by code):
  - **Theorem R (rank theorem).** For a level-`k` strategy `sigma` the following are equivalent: (a) class (i); (b) a periodic rank `R(n) = a log n + h(n mod 2^k)`, `a > 0`, decreases by a fixed `eps > 0` at every step `n -> T_sigma(n)`, `n >= n_1`; (b') some rank `a log n + h(n)` with `h` merely **bounded** (any period, or none) never increases at any step `n >= n_1`; (c) an integer potential certificate (Lemma P) for some `F < log_3 2`. The integer correction `log(1 + sigma/(3n)) <= 1/(3n)` is absorbed by half of the certificate margin.
  - **Theorem D (rank defect = maximum cycle mean).** For every strategy, `min over periodic h of limsup_n [R(Tn) - R(n)] = lambda*(sigma) = rho_max log 3 - log 2`, attained by `h = (log 3 / r) psi` with `psi` the least potential at `F = rho_max`; bounded `h` cannot do better. Flow side: `lambda*` is the maximum of the mean weight over normalized circulations (LP duality), and exactly one of "tension with `w + dh <= lambda`" and "circulation of positive `(w - lambda)`-weight" exists. For **Collatz** the least defect is **`log(3/2)`, the mean weight of the self-loop at `-1`, attained by `h = 0`**: no periodic correction improves the plain logarithm.
  - **Corollary B.** The Bernoulli-boundary obstruction (no `a log n + h(n)`, `h` bounded, descends at every step) is Theorem R's (b') => (a) for the self-loop at `-1`: its witnesses `2^H - 1` are the integer shadow of that loop.
  - **Corollary K (finite banks).** No rank `a log n + h(n mod 2^K) + sum_i c_i v2(n - beta_i)` (finitely many centers `beta_i`, arbitrary real `c_i`, any periodic `h`) is non-increasing at every large Collatz step: an expanding periodic orbit `1^(L-1)0` avoiding the bank keeps every counter periodic. The least asymptotic defect over all such ranks is still `log(3/2)`. This contains the Kuratowski lane's finite-center theorem (§5 there), with a one-period proof.
  - **Proposition H.** The Haar (uniform) edge flow is a circulation of `G_sigma` iff `sigma` is constant (Collatz or 3n-1); for Collatz it has odd density exactly `1/2`.
  - **Theorem M (maximizers of `G_(sigma_k)`).** A cycle of `G_(sigma_k)` has the maximal density `F_k = a/d` iff its orbit factors into first-descent blocks of density `F_k` (plus flip blocks `10` when `F_k = 1/2`). These blocks are exactly the words of length `md <= k` with `ma` ones lying strictly above the line of slope `F_k`; those of length `d` are one per necklace (`C(d,a)/d` of them). The upper Christoffel word is one of them, the unique balanced one and the lowest lattice path.
  - **Lemma S (sibling splitting) and the numerator bound.** A maximum-density cycle can be chosen with its nodes distinct mod `2^(k-1)`, so `rho_max = a/p` with `a <= 2^(k-2)`, `p <= 2^(k-1)`. Hence `F_k <= M_k <= G_(2^(k-2))`, where `M_k` is the largest `rho_max` over class (i) and `G_A` the best lower approximation of `log_3 2` with numerator `<= A`.
  - **Exact maxima for `k <= 6`.** `M_k = G_(2^(k-2)) = 1/2, 1/2, 3/5, 5/8, 5/8` for `k = 2..6`: the numerator bound is attained by explicit strategies (`chi_(-4)`, a level-4 and a level-5 optimum, and its lift), each certified by exact Karp.
  - **Proposition N (negation in u/d coordinates).** With `U(sigma)` = the residues where `sigma != chi_(-4)`, negation acts by `U -> -U`; `sigma` is self-dual iff `U = -U`. Every self-dual strategy is at Haar distance exactly `1/2` from Collatz and from 3n-1, and the nearest self-dual class-(i) strategy to `sigma_k` is `chi_(-4)`, at Haar distance `1/2 - |Bad_k|/2^(k-1) -> 1/2`.
- **REFUTED.**
  - **Q2 as posed** ("every maximum-density cycle of `G_(sigma_k)` is an upper Christoffel orbit") is **false for every `k >= 4`**. The counterexamples are `1100` at `k = 4`, `11100` at `k = 5..7`, six non-Christoffel blocks of length 8 for `8 <= k <= 26` (plus longer ones from `k = 16`), and 312,454 of length 27 for `27 <= k <= 45`. It is true for `k = 2, 3`. The maximizing set even has **positive entropy** for `k >= 4`.
  - **Q3's rigidity** ("`rho_max(sigma) <= F_k` for every class-(i) level-`k` strategy") is **false for every `4 <= k <= 26`**: certified witnesses have `rho_max = 3/5` at level 4, `5/8` at level 5 and `17/27` at level 7 (lifted to higher levels). The class-(i) strategies of largest `rho_max` at levels 4 and 5 have **no** Christoffel maximal cycle at all.
- **FINITE-EXACT.**
  - `M_2..M_5` re-derived exhaustively (all 65,812 strategies of levels 2–5), `M_6 = 5/8` re-confirmed by CP-SAT; `M_7` is `17/27` or `29/46`: the lower bound is an explicit certified strategy, the upper bound is the numerator bound.
  - The complete value sets `V_k` of `rho_max` over class (i): `V_4 = {1/2, 3/5}`, `V_5 = {1/2, 5/9, 4/7, 3/5, 5/8}`, and `V_6` = 16 fractions with denominators up to 22 (every realized value re-verified by exact Karp, every absent one CP-SAT-certified).
  - Negation census, levels 2–5: class (i) = `1, 1, 2 + 2·7, 12 + 2·520` (self-dual + pairs); class (i) is **connected** under single sign flips at every level `<= 5`; at level 6, 927 of the 65,536 self-dual strategies are class (i), with largest `rho_max = 5/8`.
- **EMPIRICAL / OPEN.** The law of `M_k` (data: the numerator bound is attained for `k <= 6`; at `k = 7` either the bound `29/46` is attained or `M_7 = 17/27`, and CP-SAT did not decide which); whether the rigidity fails for `k >= 27`; the observed "self-dual maximum at level `k` equals `M_(k-1)`" (`k = 3..6`).
- **ANALOGY / NUMEROLOGY.** The Tutte-type dictionary of §6: only the flow/tension (Farkas/Gallai) and max-plus entries are genuine maps.
- No HYP or THM file was created. Nothing here bears on Collatz itself (§7).

Session `collatz-procgen-20260922`, tension lane, 2026-09-26. Scripts `04-computation/experiments/procgen_tension_20260926_{lib,q1,q2,q3,q4,mip,run}.py`; output [procgen_tension_20260926.out](procgen_tension_20260926.out). Parents: [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) and its [note](procgen_cube_20260925_strategy_cube.md); the [cube-distance note](procgen_cubedist_20260925_distance_to_provability.md) (Lemma P, `sigma_k`, now [THM-4479](../../01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md)); the [Bernoulli-boundary note](bernoulli_boundary_20260925.md) §4; the [Kuratowski reframe](kuratowski_reframe_20260925.md) §5; the [Kuratowski/Tait note](procgen_kuratowski_20260925_tait_kempe_triples.md) §3.

## 0. The answers in brief

| question | answer | status |
|---|---|---|
| Q1 ranks = certificates | class (i) ⟺ periodic rank with margin ⟺ bounded-correction rank ⟺ integer potential below `c`. The least rank defect is the maximum cycle mean; for Collatz it is `log(3/2)`, at `-1`, attained by `h = 0`. Expanding cycles = positive circulations. The Haar circulation exists only for the two constant strategies, with density `1/2`. | PROVED |
| Q1 Bernoulli / Kuratowski | The bounded-correction obstruction is the `-1` loop; finite banks of dyadic counters fail on an orbit `1^(L-1)0` outside the bank. | PROVED |
| Q2 Christoffel maximizers | False for `k >= 4`. Maximizers = periodic concatenations of strictly-above blocks of slope `F_k`; Christoffel = the unique balanced, lowest block. Positive entropy. | REFUTED + PROVED characterization |
| Q2 Sturmian optimization | Bousch's unique Sturmian maximizer has no analogue here; the potential is a finite sub-action and the maximizers are its tight cycles. | PROVED (finite) + CITED |
| Q3 rigidity `rho_max <= F_k` | False for `4 <= k <= 26`. `F_k <= M_k <= G_(2^(k-2))`; `M_k` = 1/2, 1/2, 3/5, 5/8, 5/8 for `k = 2..6`; `M_7` is `17/27` or `29/46` (open). | REFUTED + PROVED bounds + FINITE-EXACT |
| Q4 duality | `nu` = negation of the u-set. Census above; class (i) connected; self-duals are Haar-`1/2` from Collatz, so no self-dual approximant of Collatz exists. | PROVED + FINITE-EXACT |
| Q5 dictionary | cycles ↔ circulations, potentials ↔ tensions, `lambda*` ↔ max-plus eigenvalue: PROVED. Tutte flow–colouring duality, planar duality, {K5, K3,3, Petersen}, {U2,4, F7, F7*}: ANALOGY or NUMEROLOGY. | typed table §6 |

## 1. Setting and notation

Notation follows THM-4474 §1 and the cube-distance note §1.

* **Strategies.** A level-`k` strategy is `sigma : {1, 3, ..., 2^k - 1} -> {+1, -1}`; `T = T_sigma` is `n/2` on even `n` and `(3n + sigma(n mod 2^k))/2` on odd `n`. It maps positive integers to positive integers, and it acts on `Z_2`. Collatz is `sigma = +`, 3n−1 is `sigma = -`.
* **The graph.** `G_sigma` has nodes `Z/2^k`; from `s` there are two edges, to the two lifts mod `2^k` of `T(s) mod 2^(k-1)`. For every integer `n`, `(n mod 2^k, Tn mod 2^k)` is an edge: if `n = s (mod 2^k)` then `Tn = T(s) (mod 2^(k-1))`, because the sign used is `sigma(s)` and `3(n - s)/2 = 0 (mod 2^(k-1))`.
* **Weights.** `w(s) = log(3/2)` for odd `s`, `-log 2` for even `s`. A closed walk with `a` odd nodes and length `p` has weight `a log 3 - p log 2`, never `0`; it is **expanding** iff `3^a > 2^p`, i.e. its odd density exceeds `c = log_3 2`.
* `rho_max(sigma)` is the largest odd density of a cycle, and `lambda*(sigma) = rho_max log 3 - log 2` the maximum cycle mean of `w`. **Class (i)** means `rho_max < c` (THM-4474 Theorem A).
* **Terras (THM-4474 Lemma 1).** The classes mod `2^(k+L)` correspond bijectively to the paths with `L` edges of `G_sigma`. On the class of a path, `T^j(n) = (3^(a_j) n + c_j)/2^j` for `j <= L`, with `a_j, c_j` constant.
* **Periodic points (THM-4474 Lemma 2).** A cycle `gamma` with signed word `(s_0, ..., s_(p-1))` (letters `0`, or the sign at an odd node) carries exactly one periodic point `x_gamma`, the fixed point of the composite affine map; it is rational with odd denominator.
* **Best approximations.** `F_D` is the largest `a/d < c` with `d <= D`; `G_A` the largest `a/p < c` with `a <= A`. The `F_k` are `1/2 (k <= 4), 3/5 (5..7), 5/8 (8..26), 17/27 (27..45), 29/46, 41/65, ...`
* **Lemma P** (cube-distance note §1). For `F = q/r` (`0 < q <= r`) put `w_F = r - q` on odd and `-q` on even nodes. Every cycle has density `<= F` iff there is an integer `psi >= 0` with `psi(t) <= psi(s) - w_F(s)` on every edge. The least such `psi(s)` is the supremum over walks from `s` of the partial sums of `w_F`. Moreover `(log 3 / r) w_F = w - mu` with `mu = F log 3 - log 2`.
* **Siblings.** `t` and `t' = t + 2^(k-1) (mod 2^k)` are the two lifts of one class mod `2^(k-1)`. They have the **same in-neighbours**.
* **u/d.** At an odd residue `r` the step is *d* (at least two halvings follow) iff `sigma(r) = chi_(-4)(r)`, where `chi_(-4)(r) = +1` for `r = 1` and `-1` for `r = 3 (mod 4)`; otherwise it is *u*. The **u-set** is `U(sigma) = {r : sigma(r) != chi_(-4)(r)}`. Collatz has `U = {3 mod 4}`, 3n−1 has `U = {1 mod 4}`, all-d `chi_(-4)` has `U = {}`, all-u has every odd residue.
* **Negation.** `(nu sigma)(m) = -sigma(-m)`. By THM-4474 Theorem E, `s -> -s` is a parity-preserving isomorphism `G_sigma = G_(nu sigma)`.

## 2. Q1: ranks are potentials; the rank defect is the maximum cycle mean; flows

### 2.1 The rank theorem

**Theorem R.** For a level-`k` strategy `sigma` the following are equivalent.
* **(a)** `sigma` is in class (i).
* **(b)** There are `a > 0`, `h : Z/2^k -> R`, `eps > 0` and `n_1` such that `R(n) = a log n + h(n mod 2^k)` satisfies `R(Tn) <= R(n) - eps` for every integer `n >= n_1`.
* **(b')** There are `a > 0`, a **bounded** function `h : Z_(>0) -> R` and `n_1` such that `R(n) = a log n + h(n)` satisfies `R(Tn) <= R(n)` for every integer `n >= n_1`.
* **(c)** There are a rational `F = q/r < c` and an integer potential `psi` as in Lemma P.

Under (c), (b) holds with `a = 1`, `h = (log 3 / r) psi`, `eps = |mu|/2` and `n_1 = ceil(2/(3|mu|))`, where `mu = F log 3 - log 2 < 0`.

*Proof.*
* **(a) ⇒ (c).** Take `F = rho_max(sigma)`. It is rational (the density of a simple cycle) and `F < c`, so Lemma P gives `psi`.
* **(c) ⇒ (b).** By Lemma P's rescaling, `h = (log 3/r) psi` satisfies `h(t) <= h(s) - w(s) + mu` on every edge. Let `n >= 1`, `s = n mod 2^k`, `t = Tn mod 2^k`; `(s, t)` is an edge (§1).
  * The **integer correction.** `log Tn - log n = w(s) + eta(n)`, where `eta(n) = 0` for even `n` and `eta(n) = log(1 + sigma(s)/(3n))` for odd `n`. Since `log(1 + x) <= x`, `eta(n) <= 1/(3n)` in both signs (for `sigma(s) = -1`, `eta < 0`).
  * Hence `R(Tn) - R(n) = w(s) + eta(n) + h(t) - h(s) <= mu + 1/(3n) <= mu/2` as soon as `1/(3n) <= |mu|/2`, i.e. `n >= 2/(3|mu|)`.
* **(b) ⇒ (b').** A periodic `h` is bounded, and `eps > 0` gives non-increase.
* **(b') ⇒ (a).** Suppose `G_sigma` has a cycle `gamma` with `a` odd nodes, length `p` and `3^a > 2^p`. Fix `m >= 1`.
  * The integers whose itinerary follows `gamma` `m` times form one class mod `2^(k+mp)` (Terras). On it, `T^j(n) = (3^(a_j) n + c_j)/2^j`, so for `n` large every `T^j(n)`, `j <= mp`, is `>= n_1`, and `T^(mp)(n)/n -> 3^(ma)/2^(mp)` as `n -> infinity` in the class.
  * Chaining the `mp` inequalities, `R(T^(mp) n) <= R(n)`, i.e. `a log(T^(mp)(n)/n) <= h(n) - h(T^(mp) n) <= 2 ||h||`.
  * Letting `n -> infinity` gives `a m log(3^a/2^p) <= 2 ||h||` for every `m`, impossible since `log(3^a/2^p) > 0`. ∎

**Remarks.**
* **Any period.** (b') contains every periodic correction, whatever its period (not only powers of 2). So *a periodic rank of any period exists iff class (i)*. For a periodic `h` the proof needs one period only: `h(T^p n) = h(n)` when `h` has period dividing `2^k`.
* **Explicit lookahead from a rank.** Under (b) with `a = 1`, `log(T^L n / n) <= osc(h) - eps L < 0` once `L > osc(h)/eps`, as long as the orbit stays above `n_1`. An orbit that drops below `n_1 <= n` has descended already. So every `n >= n_1` descends below itself within `floor(osc(h)/eps) + 1` steps. With the certificate of (c) this is `L <= 2 (log 3/r) max psi / |mu| + 1`, an alternative to THM-4474's bound `2^k (1 + log(3/2)/|mu|) + 1`.
* **Lookahead ranks.** Requiring only `R(T^j n) <= R(n)` for some `j <= L` gives the same class: chain the jumps along `gamma^m` and compare at a time in `[mp - L, mp]`.

### 2.2 The rank defect and the flow side

For `h : Z/2^k -> R` let `delta(h) = max over edges s -> t of [w(s) + h(t) - h(s)]`.

**Theorem D.**
* **(i)** For every `h`, `limsup over n -> infinity of [log Tn - log n + h(Tn mod 2^k) - h(n mod 2^k)] = delta(h)`.
* **(ii)** `min over h of delta(h) = lambda*(sigma) = rho_max log 3 - log 2`, attained by `h* = (log 3/r) psi`, where `psi` is the least potential at `F = rho_max = q/r`.
* **(iii)** Bounded arbitrary corrections do no better: for every bounded `h : Z_(>0) -> R`, `limsup_n [log Tn - log n + h(Tn) - h(n)] >= lambda*(sigma)`.
* **(iv) Alternative.** For real `lambda` exactly one holds: (T) a tension `dh` (`dh(s -> t) = h(t) - h(s)`) with `w + dh <= lambda` on every edge; (F) a nonzero circulation `f >= 0` with `sum_e f_e (w(tail e) - lambda) > 0`. For `lambda = 0`: *a descent tension exists iff there is no positive circulation, iff there is no expanding cycle.*
* **(v) LP duality.** `lambda*(sigma) = max { sum_e f_e w(tail e) : f >= 0 a circulation, sum_e f_e = 1 }`.

*Proof.*
* **(i).** Every edge is the 1-path of a class mod `2^(k+1)`, which contains arbitrarily large integers, and `eta(n) -> 0`.
* **(ii).** Lower bound: along any cycle the terms `h(t) - h(s)` telescope, so the largest edge value is at least the cycle's mean weight, and at least `lambda*` for a cycle of maximal density. Upper bound: Lemma P's rescaling at `F = rho_max` gives `w(s) + h*(t) - h*(s) <= mu = lambda*` on every edge.
* **(iii).** Follow a maximal cycle `gamma` `m` times, as in Theorem R. The average of the `mp` step increments is `(1/mp)[log(T^(mp) n / n) + h(T^(mp) n) - h(n)] >= lambda* - 2||h||/(mp) - o(1)` as `n -> infinity`. Some step is at least the average; let `m -> infinity`.
* **(iv).** If both held, `0 < sum f (w - lambda) <= sum_e f_e (h(s) - h(t)) = 0`, because a circulation is orthogonal to every tension. If (T) fails, some cycle has mean `> lambda`: otherwise `h(s) = sup` over walks from `s` of the partial sums of `(w - lambda)` (finite, attained on paths) would satisfy `h(s) >= w(s) - lambda + h(t)`, a valid (T) (Lemma P's argument with real weights). Its indicator is a circulation as in (F).
* **(v).** A circulation `f >= 0` decomposes into cycles, `f = sum lambda_gamma 1_gamma` (THM-4474 Theorem C). So the normalized mean weight is a mediant of cycle means, hence `<= lambda*`; a maximal cycle, normalized, attains it. ∎

*Remark (finite versus asymptotic defect).* Theorem D uses the `limsup` over `n`. The maximum over all `n >= 1` adds the integer correction `eta(n) <= 1/(3n)`, which only matters for small `n`. For Collatz with `h = 0` the maximum is `log 2`, at `n = 1 -> 2`. Along the loop at `-1` it is `log(3/2) + log(1 + 1/(3n))` for `n = -1 (mod 2^(k+1))`. The `limsup` is the level-independent quantity, and it is the one that decides the existence of a rank above some `n_1`.

This is the directed case of the classical potential / negative-cycle alternative (Gallai; the flow–tension pair of Hoffman's circulation theorem and Minty's painting lemma). The proof above is self-contained. In ergodic-optimization language, `h*` is a **sub-action** and (v) is the finite Mañé/Conze–Guivarc'h principle (§3.3).

### 2.3 Collatz: the defect is `log(3/2)`, at `-1`, and the Bernoulli obstruction is that loop

For Collatz at every level `k >= 2`, node `2^k - 1` has the self-loop `-1 -> -1`: `T(-1) = -1`, and `2^k - 1` is a lift of `-1 mod 2^(k-1)`. Hence `rho_max = 1` and `lambda* = log(3/2)`. At `F = 1`, `w_F = 0` on odd and `-1` on even nodes, so the least potential is `psi = 0`.

**Corollary (Collatz rank defect).** At every level, the least achievable per-step defect of a rank `log n + h(n mod 2^k)` is exactly `log(3/2)`, the mean weight of the self-loop at `-1`. It is attained by `h = 0`: **no periodic correction improves on the plain logarithm**, and a bounded correction cannot either (Theorem D(iii)).

For comparison, the least defect `lambda* = rho_max log 3 - log 2` of some strategies (the best margin of a periodic rank is `-lambda*` when it is negative):

| strategy | `rho_max` | `lambda*` | maximal cycles |
|---|---|---|---|
| Collatz (3n+1) | 1 | +0.40547 | the loop at `-1` |
| 3n−1 | 1 | +0.40547 | the loop at `1` |
| all-u (`-chi_(-4)`) | 1 | +0.40547 | all-odd cycles (loops at `1` and `-1`) |
| `chi_(-4)` (all-d) | 1/2 | −0.14384 | `{1,2}`, `{-1,-2}`, ... |
| `sigma_5`, `sigma_6`, `sigma_7` | 3/5 | −0.03398 | blocks `11010`, `11100` |
| `sigma_8`, ..., `sigma_26` | 5/8 | −0.00651 | 7 blocks of length 8 (and longer from `k = 16`) |
| the level-7 witness (§4.3); `sigma_27`, ..., `sigma_45` | 17/27 | −0.00143 | a 27-cycle; 312,455 blocks of length 27 |

**Corollary B (the Bernoulli-boundary obstruction).** The Bernoulli note (§4) proves that no `V(n) = a log n + h(n)` with bounded `h` decreases at every positive odd Collatz step `n > 1`, with witnesses `n_H = 2^H - 1`. In the present language:
* The witnesses are the **integers following the self-loop at `-1`**. `2^H - 1 = -1 (mod 2^k)` for `H >= k`, and `T(-1 + 2^H m) = -1 + 3·2^(H-1) m`, so the orbit stays on node `-1` for `H - k + 1` steps (checked, `k <= 11`, `H <= k + 39`).
* Their proof is Theorem R's (b') ⇒ (a) applied to `gamma = {-1}`. Theorem R shows this is the **only** kind of obstruction: a bounded-correction rank exists iff every cycle of `G_sigma` is contracting. For Collatz the obstruction is the densest cycle available (density 1), and Theorem D says it costs exactly `log(3/2)` per step.
* It is also the cube's reading of THM-4474 Proposition F (`sigma(-1) = +` forbids class (i)) and of the pair-0 obstruction of THM-4470(5).

### 2.4 Finite banks of dyadic counters

**Corollary K.** Fix `K >= 1`, `h : Z/2^K -> R`, finitely many centers `beta_1, ..., beta_s` in `Q_2` that are not positive integers, real `c_i` and `a > 0`. Put `R(n) = a log n + h(n mod 2^K) + sum_i c_i v2(n - beta_i)`. Then `R(Tn) <= R(n)` fails for some arbitrarily large `n`, for Collatz `T`. Moreover, the least asymptotic defect `limsup_n [R(Tn) - R(n)]/a` over all such ranks is `log(3/2)`.

*Proof.*
* **An orbit outside the bank.** For `L >= 3` the Collatz periodic orbit `O_L` of the word `1^(L-1) 0` is expanding (`3^(L-1) > 2^L`) and has exactly `L` points. Distinct `L` give disjoint orbits, so all but finitely many `O_L` avoid the bank. Fix one, with points `x_0, ..., x_(L-1)`, and let `V = max v2(x_j - beta_i)`, which is finite.
* **Its integer shadow.** Take `N >= K + L + V + 1` and a large integer `n = x_0 (mod 2^N)`. For `j <= L`, `T^j(n) = T^j(x_0) (mod 2^(N-j))`, because `T` maps a class mod `2^m` into a class mod `2^(m-1)`. At `j = L`, `T^L(n) = x_0 = n (mod 2^(K+V+1))`.
* **Everything periodic cancels.** Therefore `h(T^L n) = h(n)` and `v2(T^L n - beta_i) = v2(x_0 - beta_i) = v2(n - beta_i)` for every `i`. So `R(T^L n) - R(n) = a log(T^L n / n)`, which is `> 0` for large `n`, since `T^L n / n -> 3^(L-1)/2^L > 1`. The `L` step inequalities would force it to be `<= 0`.
* **Defect.** The same computation over the period gives a step with increment at least `(a/L) log(T^L n/n) -> a((L-1) log 3 - L log 2)/L`. This tends to `a log(3/2)` as `L -> infinity` through the orbits that avoid the bank; `h = 0` without a bank attains `log(3/2)`. ∎

**Scope.**
* This contains the Kuratowski reframe's §5 theorem (negative rational centers, `h = 0`, descent at every odd step `n > 1`) with a one-period proof. Its version for Syracuse steps is the same argument: over one period the odd-to-odd steps of `O_L` compose to `T^L`.
* It is still only a statement about **finitely many** counters. Infinitely many expanding periodic orbits (THM-4479 Theorem 2 counts `>= N_k` of them at level `k`) would each need one; unbounded banks, adaptive centers and nonlinear couplings are not excluded.

### 2.5 The Haar circulation

**Proposition H.** The uniform edge flow (`1` on every edge) is a circulation of `G_sigma` iff `sigma` is constant. For Collatz it is the uniform flow on the de Bruijn graph `B(2,k)`, and its odd density is exactly `1/2`: the Haar drift `(1/2) log 3 - log 2 < 0`.

*Proof.* The in-neighbours of a node `t` are:
* the even node `2t mod 2^k`, always;
* the odd node `(2t - 1)/3`, iff its sign is `+`;
* the odd node `(2t + 1)/3` (inverses mod `2^k`), iff its sign is `-`.

So in-degrees lie in `{1, 2, 3}`, and the uniform flow is a circulation (every in-degree 2) iff `sigma((2t-1)/3) = sigma((2t+1)/3)` for all `t`, i.e. `sigma(s) = sigma(s + 2/3)` for every odd `s`. Since `2/3 = 2·(odd) (mod 2^k)`, the translation by `2/3` is a single cycle on the odd residues, which forces `sigma` constant. For constant `sigma`, every in-degree is 2, and half of the nodes are odd. ∎

So the Haar measure sits exactly at the drift `1/2` for Collatz, inside its window `[rho_min, rho_max] = [0, 1]`. For every non-constant strategy the invariant flows are the non-uniform stationary flows of THM-4474 Theorem C.

**Checks (output §A–§C).**
* **(a) ⇒ (c)** on all 1,070 class-(i) strategies of levels 2–5 (potential at `F = rho_max`, every edge checked).
* **not (a) ⇒ not (c)** on all 258 non-(i) strategies of levels 2–4 and 3,000 sampled at level 5 (no potential at `F_(2^k)`, the largest possible cycle density below `c`).
* **(c) ⇒ (b)**, numerically: `R(Tn) - R(n) <= -eps` on 2,000 consecutive integers above `n_1`, for every class-(i) strategy of levels 2–4, 200 at level 5, and `sigma_8`, `sigma_10`, `sigma_12`.
* **(b) ⇒ (a)** witnesses: for every non-(i) strategy of levels 2–4, 3,000 at level 5 and 180 random ones at levels 6–8, an integer `n > 10^12` with `T^p(n) = n (mod 2^k)` and `T^p(n) > n`, checked exactly.
* Theorem D(ii) on all 276 strategies of levels 2–4.
* The Collatz loop and `psi = 0` at `k = 2..16`.
* The bank witnesses for `L = 5, 7, 9, 13, 23` against a 20-center bank.
* Proposition H on all 276 strategies of levels 2–4.
* LP duality (v), numerically, on 276 strategies.

## 3. Q2: the maximum-density cycles of `G_(sigma_k)`

`sigma_k` is Collatz with the signs flipped exactly on `Bad_k`, the residues whose Collatz word of length `k` has every prefix multiplier `> 1` (cube-distance Theorem 1, THM-4479). Its orbits split into blocks:
* an **F-block** `x -> (3x-1)/2 -> (3x-1)/4`, word `10`, at `x` in `Bad_k`;
* a **C-block**, elsewhere: the Collatz segment up to the first descent, a *first-descent word* of length `d(x) <= k`.

Every cycle's density is a mediant of its block densities, and `rho_max(sigma_k) = F_k`.

### 3.1 The theorem

Write `F_k = a/d` in lowest terms. A binary word of length `L` is **strictly above** slope `F` if its prefix counts satisfy `a_j > F j` for `0 < j < L`.

**Theorem M.**
* **(i)** A cycle of `G_(sigma_k)` has density `F_k` iff every block of its periodic orbit is a C-block of density `F_k`, or an F-block (the latter possible only when `F_k = 1/2`, i.e. `k <= 4`).
* **(ii)** The first-descent words of density `F_k` and length `<= k` are exactly the words of length `md <= k` with `ma` ones that are strictly above slope `F_k`.
* **(iii)** Every periodic concatenation of such words is the word of a periodic `T_(sigma_k)`-orbit (with Collatz signs throughout). Hence it gives a closed walk of density `F_k` in `G_(sigma_k)`, whose simple cycles are all maximal.
* **(iv)** The strictly-above words of length `d` are exactly one per necklace of `(d, a)`-words, `C(d,a)/d` of them. The upper Christoffel word `C_F` (prefix counts `ceil(F j)`) is one of them. It is the unique balanced one, and its lattice path is pointwise lowest: every strictly-above word `w` of length `d` has `a_j(w) >= ceil(F j)`, hence also the smallest maximal excursion `max_j 3^(a_j)/2^j`.
* **(v)** Consequently every maximal cycle of `G_(sigma_k)` carries a power of the Christoffel word `10` when `k = 2, 3`. For every `k >= 4` there are non-Christoffel maximizers: `1100` at `k = 4`, and `C(d,a)/d - 1 >= 1` non-Christoffel blocks of length `d >= 5` for `k >= 5`. For `k >= 4` the set of maximal closed walks has **positive entropy**.

*Proof.*
* **(i).** The blocks have densities `1/2` (F) and `a'/d' <= F_k` (C). A mediant equals `F_k` iff every term does.
* **(ii).** For `j <= k` the interval `(F_k j, c j]` contains no integer: an integer there would give a lower approximation of `c` with denominator `j <= k` exceeding `F_k` (and `c j` is irrational). So for `j <= k`, `a_j > c j` iff `a_j > F_k j`. A first-descent word of length `L <= k` and density `F_k` has `3^(a_j) > 2^j` for `j < L` and `a_L = F_k L < c L`; since `a_L` is an integer, `d | L`.
* **(iii).** Let `x` be the Collatz periodic point of the concatenation. At a block start the next `<= k` letters contain a descent, so `x` is not in `Bad_k`. By heredity (cube-distance Lemma 1.1) no point inside a block is in `Bad_k` either. So `sigma_k = +` along the whole orbit, and it is a `T_(sigma_k)`-orbit.
* **(iv).** Give a letter `b` the weight `d b - a`. The total over a length-`d` word with `a` ones is `0`, and the partial sums `P_j = d a_j - a j` are distinct mod `d` for `0 <= j < d` (as `gcd(a,d) = 1`). So the minimum over a period is attained exactly once, and the rotation starting there is the only one with all interior partial sums `> 0`: this is the cycle lemma. All `(d,a)`-words are primitive, so there are `C(d,a)/d` necklaces.
  * `C_F` has `ceil(aj/d) > aj/d` for `0 < j < d`, so it is strictly above.
  * Any strictly-above `w` has `a_j(w) >= floor(aj/d) + 1 = ceil(aj/d)`.
  * The balanced words of slope `a/d` form exactly the conjugacy class of the Christoffel word (CITED: Berstel, Lauve, Reutenauer, Saliola, *Combinatorics on Words: Christoffel Words and Repetitions in Words*, CRM Monograph Series 27, AMS 2008). So exactly one block is balanced, namely `C_F`. This is also checked directly for `d = 2, 5, 8, 27`.
* **(v).** For `k = 2, 3` the only blocks are `10` (C and F). At `k = 4`, `1100` is strictly above slope `1/2` and has length `4 <= k`. For `k >= 5`, `d >= 5` and `2 <= a <= d - 2`, so `C(d,a)/d >= 2`. With `B >= 2` blocks, the free concatenations of total length `L` number at least `2^(floor(L/d_max))`, so the entropy is positive. ∎

### 3.2 The catalogue (FINITE-EXACT)

| `k` | `F_k` | blocks by length | Christoffel block | simple maximal cycles (words) | growth rate of maximal walks |
|---|---|---|---|---|---|
| 2 | 1/2 | {2: 1} | `10` | 2 (`10`, `10`) | 1 (Christoffel only) |
| 3 | 1/2 | {2: 1} | `10` | 3 (`10` ×2, `1010`) | 1 |
| 4 | 1/2 | {2: 1, 4: 1} | `10` | 8 (`10` ×2, `1010`, **`1100`**, `110010` ×2, `11001010` ×2) | 1.5538 (flip blocks tight too) |
| 5–7 | 3/5 | {5: 2} | `11010` | 3 (`11010`, **`11100`**, `1110011010`) | `2^(1/5)` = 1.14870 |
| 8–15 | 5/8 | {8: 7} | `11011010` | — | `7^(1/8)` = 1.27537 |
| 16–23 | 5/8 | {8: 7, 16: 476} | `11011010` | — | 1.49976 (`7 lambda^-8 + 476 lambda^-16 = 1`) |
| 24–26 | 5/8 | {8: 7, 16: 476, 24: 51,033} | `11011010` | — | 1.60517 |
| 27–45 | 17/27 | {27: 312,455} | `110110110101101101011011010` | — | `312455^(1/27)` = 1.598 |

* The last column is the root `lambda > 1` of `sum over blocks of lambda^(-|w|) = 1` (the growth rate of block concatenations). It equals the Perron root of the tight subgraph for `k = 5..11`, to `1e-6` (numerical).
* All counts, the uniqueness of the balanced block and the lowest-path property are checked for `k = 2..30`. The realization (iii) is checked for every block with `k <= 14`, and the exhaustive cycle list for `k <= 7`.

**The global maximizers are anti-Christoffel.** Over all class-(i) strategies (§4), the largest `rho_max` at level 4 is `3/5`, attained by 8 strategies. Every maximal cycle of every one of them has the word `11100`. At level 5 (`5/8`, 220 strategies) the maximal words are `11100110`, `11111000`, `11101100`, `11110010`. The Christoffel words `11010` and `11011010` never occur. The strategies furthest from Collatz-type provability put their ones in runs.

### 3.3 Relation to Sturmian ergodic optimization

* **The finite sub-action principle (PROVED).** Let `psi` be the least potential at `F = rho_max` and call an edge *tight* if `psi(t) = psi(s) - w_F(s)`.
  * Around any closed walk the slacks `psi(s) - w_F(s) - psi(t) >= 0` sum to `r p (F - a/p)`.
  * So a closed walk has maximal density iff all its edges are tight. The maximizing invariant measures are exactly those carried by the tight subgraph (its "Mather set").
  * This is the finite-graph form of the sub-action / revelation principle of ergodic optimization, and `(log 3 / r) psi` is a sub-action for `w - lambda*`.
* **Bousch (CITED).** T. Bousch, *Le poisson n'a pas d'arêtes*, Ann. Inst. H. Poincaré Probab. Statist. 36 (2000), no. 4, 489–508 (zbMATH 0971.37001). Among the probability measures invariant under `x -> 2x` on the circle, for each `omega` the function `cos(2 pi (t - omega))` has exactly one maximizing measure. It is Sturmian (supported in a semicircle), and periodic except for `omega` in a set of measure and Hausdorff dimension zero.
* **The contrast (PROVED here).**
  * In Collatz parity coordinates, the C-block orbits of `sigma_k` are shift orbits, and `T` is the doubling map on `Z_2`.
  * The functional is the frequency of the letter `1`. The constraint is the subshift of finite type that forbids the `Bad_k` windows (the ballot words of length `k`).
  * For this pair the maximizing set is **not** unique and has positive entropy for `k >= 4` (Theorem M(v)), so there is no Bousch-type rigidity.
  * The Sturmian (balanced) maximizer exists and is unique among single blocks: it is the upper-Christoffel orbit, whose lattice path hugs the critical line.
* **Summary.** In ergodic-optimization terms, `sigma_k`'s maximizing set is a positive-entropy Mather set containing exactly one Sturmian orbit. The Christoffel orbit is distinguished by balance and minimal excursion, not by maximality.

## 4. Q3: the rigidity `rho_max <= F_k` is false; exact maxima

Let `M_k` be the largest `rho_max` over class-(i) level-`k` strategies, and `V_k` the set of values of `rho_max` on class (i).

### 4.1 Bounds

**Lemma S (sibling splitting).** Let a simple cycle `gamma = (v_0 -> ... -> v_(p-1) -> v_0)` contain siblings `v_i`, `v_j = v_i + 2^(k-1)` (`i < j`). Then `v_(i-1) -> v_j` and `v_(j-1) -> v_i` are edges. So `gamma` splits into the two simple cycles `(v_i, ..., v_(j-1))` and `(v_j, ..., v_(i-1))`, which partition its nodes; odd counts and lengths add.

*Proof.* The out-edges of `v_(i-1)` go to both lifts of `v_i mod 2^(k-1)`, one of which is `v_j`; similarly for `v_(j-1)`. ∎

**Corollary S.** Every strategy has a maximal cycle whose nodes are pairwise distinct mod `2^(k-1)`: split a shortest maximal cycle, and note that both parts must be maximal, since their mediant is. Such a cycle has at most `2^(k-2)` odd nodes (the odd residues mod `2^(k-1)`). Hence:
* `rho_max(sigma) = a/p` with `a <= 2^(k-2)` and `p <= 2^(k-1)`, for **every** strategy;
* for class (i), `rho_max <= G_(2^(k-2))`, the best lower approximation of `log_3 2` with numerator `<= 2^(k-2)`.

Also `rho_max >= 1/2` on class (i): Proposition F gives `sigma(1) = +`, and then `1 -> 2 -> 1` is a cycle.

**Bounds.** `F_k <= M_k <= G_(2^(k-2)) <= F_(2^(k-1))`, and `M_k` is nondecreasing in `k`, since a lift keeps the map. The lower bound is `sigma_k`.

| `k` | `F_k` | `M_k` | `G_(2^(k-2))` | `F_(2^(k-1))` | status of `M_k` |
|---|---|---|---|---|---|
| 2 | 1/2 | 1/2 | 1/2 | 1/2 | PROVED (`chi_(-4)`); exhaustive |
| 3 | 1/2 | 1/2 | 1/2 | 1/2 | PROVED (`chi_(-4)`); exhaustive |
| 4 | 1/2 | **3/5** | 3/5 | 5/8 | PROVED (bound + witness); exhaustive: 8 strategies attain |
| 5 | 3/5 | **5/8** | 5/8 | 5/8 | PROVED; exhaustive: 220 attain |
| 6 | 3/5 | **5/8** | 5/8 | 17/27 | PROVED (bound + lifted witness); CP-SAT confirms; the flow-formulation MIP (HiGHS) reproduces the answers at levels 4–5 and finds a level-6 class-(i) strategy with a cycle above `3/5` |
| 7 | 3/5 | `17/27` or `29/46` | 29/46 | 29/46 | lower bound: certified witness (§4.3); `29/46` undecided (refined CP-SAT, 2 h, see §8) |
| 8 | 5/8 | `>= 17/27` (lift) | 41/65 | 41/65 | no witness above `17/27` found (refined search for `41/65`, complete, 30 min; unrefined search for density `> 17/27`, 15 min; both undecided) |
| 9 | 5/8 | `>= 17/27` | 94/149 | 147/233 | — |

So the **numerator bound is attained for every `k <= 6`**, and `M_k` is known exactly there (PROVED: Corollary S is a hand proof, the witnesses are checked by exact Karp).

### 4.2 The value sets

* `V_2 = V_3 = {1/2}`.
* `V_4 = {1/2, 3/5}`.
* `V_5 = {1/2, 5/9, 4/7, 3/5, 5/8}` (exhaustive: 128, 32, 104, 568, 220 strategies).
* `V_6` (CP-SAT, one query per fraction in `[1/2, 5/8]` with denominator `<= 32`; 42 queries): exactly
  `1/2, 8/15, 7/13, 6/11, 5/9, 9/16, 4/7, 11/19, 7/12, 10/17, 13/22, 3/5, 11/18, 8/13, 13/21, 5/8`.
  * Each realized value is re-verified by exact Karp on the returned strategy.
  * The 26 absent candidates are `9/17`, `10/19`, `11/20`, `11/21` and all 22 candidates with denominator `>= 23`. Nine of them (numerator `>= 17`) are excluded by the numerator bound alone; CP-SAT excludes the others.
* Denominators relative to `2^k`: the largest denominators in `V_k` are `2, 5, 9, 22` for `k = 3..6`, against the proved bound `2^(k-1) = 4, 8, 16, 32`. The largest numerators are `1, 3, 5, 13`, against `2^(k-2) = 2, 4, 8, 16`. The maximizing cycle of `M_7 >= 17/27` has 17 odd nodes out of the 32 allowed.

### 4.3 The level-7 witness

The strategy with mask `12481914054424834826` (bit `i` = sign minus at residue `2i+1`) is class (i) with `rho_max = 17/27`, by exact Karp and an edge-checked integer potential at `17/27`.
* Its maximal cycle is `2, 65, 97, 18, 73, 45, 3, 68, 98, 113, 41, 61, 27, 40, 84, 106, 117, 111, 39, 59, 25, 37, 120, 60, 30, 15, 87`, with canonical word `111111000110110111001111000`: very unbalanced, not Christoffel.
* Its periodic point is `-326096134/5077565`, where `5077565 = 2^27 - 3^17`.
* Lifted, it refutes `rho_max <= F_k` at every level `7 <= k <= 26`. The level-4 and level-5 optima cover `k = 4..6`.

### 4.4 Answer to Q3

* **REFUTED** for every `4 <= k <= 26`, **true** for `k <= 3`, **OPEN** for `k >= 27`, where `F_k >= 17/27` and no witness above `F_k` is known.
* **The true maximum:** `M_k = G_(2^(k-2)) = 1/2, 1/2, 3/5, 5/8, 5/8` for `k = 2..6` (PROVED); `M_7` is `17/27` or `29/46`, open (§4.3). The values of `rho_max` are fractions `a/p` with `a <= 2^(k-2)` and `p <= 2^(k-1)` (PROVED); which of them occur is irregular (`V_6` is not an interval of any Farey sequence).
* Two of my own conjectures failed on `V_6`: "all fractions in `[1/2, M_k]` with denominator `<= 2^(k-2) + 1`" and "at most `2^(k-3)` even nodes" (§7).

### 4.5 How the solver queries are encoded

All queries are exact integer models; every positive answer is re-verified by exact Karp.
* **Variables.** One Boolean per odd residue (its sign), and an integer potential `psi` per node.
* **The certificate side.** Each of the (two or four) candidate edges out of a node carries Lemma P's inequality `psi(t) <= psi(s) - w_F(s)`, enforced when the edge is present. `F = F_(2^k)` encodes class (i), and `F = f` encodes "all cycles `<= f`".
* **The obstruction side.** It is encoded in one of two ways.
  * **CP-SAT:** an `AddCircuit` cycle on present edges, with the linear density constraint `p_0 · #odd - a_0 · #nodes >= 1` (for "a cycle of density `> a_0/p_0`") or `>= 0` (for "exactly `f`").
  * **MIP (HiGHS, `..._mip.py`):** a nonzero 0/1 **circulation** of positive `(p_0 [odd] - a_0)`-weight. It is an independent formulation: the obstruction is written as a flow, and Theorem D(iv) makes the two equivalent.
* **The refined CP-SAT query** (`cpsat_refined`) adds Lemma S: the cycle has no two siblings, length exactly `p` and exactly `a` odd nodes. It is complete when `2p > 2^(k-1)` (then every sibling-free maximal cycle has length exactly `p`), and otherwise only a restricted search.
* CP-SAT runs with 2 workers.

## 5. Q4: the negation duality

**Proposition N.**
* **(i)** `U(nu sigma) = -U(sigma)`, so `sigma` is self-dual iff `U(sigma) = -U(sigma)`. There are `2^(2^(k-2))` self-dual strategies: the `2^(k-1)` odd residues form `2^(k-2)` pairs `{m, -m}`, and each pair has one member `= 1` and one `= 3 (mod 4)`.
* **(ii)** The Hamming distance from `sigma` to the nearest self-dual strategy is the number of pairs on which `U(sigma)` is asymmetric, `|U Δ (-U)|/2`. For Collatz and 3n−1 it is `2^(k-2)`, and **every** self-dual strategy is at exactly this distance (Haar `1/2`) from both.
* **(iii)** The nearest self-dual class-(i) strategy to `sigma_k` is `chi_(-4)` (all-d), at distance `2^(k-2) - |Bad_k|`, i.e. Haar `1/2 - |Bad_k|/2^(k-1) -> 1/2`.
* **(iv)** Every class-(i) strategy puts a minus sign on an odd point of each of the integer 3n+1 cycles `{-1}`, `{-5,-7,-10}`, `{-17,...,-34}`, and a plus sign on an odd point of each 3n−1 cycle `{1}`, `{5,7,10}`, `{17,...,34}`.

*Proof.*
* **(i).** `r` is in `U(nu sigma)` iff `-sigma(-r) != chi(r)` iff `sigma(-r) != chi(-r)` (as `chi` is odd) iff `-r` is in `U(sigma)`. The pairs statement: `-m = 3 (mod 4)` iff `m = 1 (mod 4)`.
* **(ii).** A self-dual `tau` has `U(tau)` a union of whole pairs, and `d(sigma, tau) = |U(sigma) Δ U(tau)|`; each asymmetric pair costs exactly 1. For Collatz, `U = {3 mod 4}` meets every pair in exactly one element.
* **(iii).** `U(sigma_k) = {3 mod 4} \ Bad_k` (Bad_k lies in `3 mod 4`), so `U ∩ (-U)` is empty and every pair meeting `U` is asymmetric. So the distance to any self-dual is at least `|U(sigma_k)| = 2^(k-2) - |Bad_k|`, attained by `U = {}`, i.e. `chi_(-4)`, which is class (i).
* **(iv).** If `sigma = +` at every odd residue of a negative integer cycle, its itinerary is a closed walk of `G_sigma` with density `1`, `2/3` or `7/11`, all `> c`. The positive case is its `nu`-image. ∎

**Census (FINITE-EXACT, all strategies of levels 2–5; exact Karp).**

| level | class (i) | self-dual | `nu`-pairs | `rho_max` of the self-dual ones | `rho_max` of the pairs (number of pairs) | flip-components | nearest self-dual class (i): distance histogram | symmetric core `U ∩ -U` class (i) |
|---|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 0 | 1/2 | — | 1 | {0: 1} | 1/1 |
| 3 | 1 | 1 | 0 | 1/2 | — | 1 | {0: 1} | 1/1 |
| 4 | 16 | 2 | 7 | 1/2 ×2 | 1/2 ×3, 3/5 ×4 | 1 | {0: 2, 1: 6, 2: 6, 3: 2} | 16/16 |
| 5 | 1,052 | 12 | 520 | 1/2 ×8, 3/5 ×4 | 1/2 ×60, 5/9 ×16, 4/7 ×52, 3/5 ×282, 5/8 ×110 | 1 | {0: 12, 1: 82, 2: 222, 3: 314, 4: 258, 5: 126, 6: 34, 7: 4} | 972/1052 |

* `rho_max` is `nu`-invariant on all 65,812 strategies of levels 2–5 (THM-4474 Theorem E, re-checked).
* **Level 6, all `2^16` self-dual strategies (exact Karp):** 927 are class (i), with `rho_max` distribution `1/2: 256, 5/9: 34, 4/7: 214, 3/5: 276, 5/8: 147`. The self-dual maximum is `5/8 = M_6 = M_5`.
* **Observed (EMPIRICAL):** the largest `rho_max` over self-dual class-(i) strategies at level `k` equals `M_(k-1)` for `k = 3..6` (`1/2, 1/2, 3/5, 5/8`).

**Answer to "does every class-(i) class contain, or lie close to, a self-dual member?"**
* **As a set, yes, at every level checked.** Class (i) is a single connected component of the flip graph at each level `k <= 5`, and it contains `chi_(-4)` and all self-dual class-(i) strategies. At level 5 every class-(i) strategy is within 7 flips (Haar `7/16`) of a self-dual class-(i) strategy. Keeping only the symmetric part of the u-set works for 972 of the 1,052.
* **Metrically near Collatz, no (PROVED).** Every self-dual strategy is at Haar distance exactly `1/2` from Collatz, while class (i) comes within Haar `2^(-0.05k)` of it (THM-4479). By the triangle inequality, a class-(i) strategy at distance `delta` from Collatz is at distance `>= 2^(k-2) - delta` from every self-dual one. So the approximants of Collatz are necessarily far from self-duality, and their `nu`-partners approach 3n−1.
* **The level-2 square** (Collatz ↔ 3n−1 a `nu`-pair, `chi_(-4)` and `-chi_(-4)` self-dual) is, in u/d coordinates, the Boolean square on the u-set with `nu` swapping the two atoms `{1}` and `{3}`. Its fixed points are the bottom (all-d, class (i)) and the top (all-u, class (ii)). From level 3 on, the self-dual class-(i) strategies are a thin, Collatz-distant part of class (i): 12 of 1,052 at level 5.

## 6. Q5: a typed dictionary

| correspondence | what is exactly true | type |
|---|---|---|
| cycles of `G_sigma` ↔ flows (circulations) | Every nonnegative circulation is a positive combination of simple cycles (THM-4474 Thm C). `lambda*` = max mean weight over normalized circulations (Thm D(v)). Expanding cycle ⟺ circulation of positive weight (Thm D(iv)). | **PROVED** |
| potentials ↔ tensions | A certificate `psi` is a potential, and `d psi` a tension with the one-sided bound `d psi <= -w_F`. Lemma P and Thm D(iv) are the directed feasible-tension (Gallai / Farkas) alternative: tension certificate ⟺ no positive circulation. The rank theorem R identifies ranks with these tensions. | **PROVED** |
| `rho_max`, `lambda*` ↔ max-plus spectral data | `lambda*` is the max-plus eigenvalue of the weighted adjacency matrix. `h*` is a max-plus sub-eigenvector, and the tight subgraph is its critical graph (Karp's theorem = the min-max of Thm D(ii)). | **PROVED** |
| Tutte's flow–colouring duality (planar `G`: nowhere-zero `k`-flows ↔ proper `k`-colourings of `G*`) | The pair (potential certificate, expanding-cycle obstruction) is Farkas duality **inside one digraph**; there is no colouring and no dual graph. The common ancestor is genuine: Lemma P and Minty's cycle-ratio criterion for colourability are both instances of the feasible-tension theorem (Minty 1962, Amer. Math. Monthly 69; recalled, not re-read). No map from strategies to colouring problems is claimed. | **ANALOGY** |
| planar duality ↔ ? | Cycle space ⊥ cut space inside `G_sigma` (circulations ⊥ tensions), which planar duality would realize geometrically. `G_sigma` is not planar in general and no dual graph is used. The involution `nu` is a graph isomorphism `G_sigma = G_(nu sigma)`, not a planar duality. | orthogonality **PROVED**; "planar duality" **ANALOGY** |
| {K5, K3,3, Petersen}: circular flow numbers 2, 3, 5; girths 3, 4, 5; Kuratowski's two planarity obstructions + Tutte's 4-flow obstruction | Standard facts about these graphs (CITED). No map to the cube is known or claimed. | reference only |
| {K5, K3,3, Petersen} ↔ the integer expanding cycles `{-1}`, `{-5,-7,-10}`, `{-17,...,-34}` (periods 1, 3, 11; densities 1, 2/3, 7/11) | Every class-(i) strategy must break all three (Prop N(iv), PROVED). But they are **not** a finite obstruction list: class (i) must break every expanding cycle, at least `N_k -> infinity` of them at level `k` (THM-4479 Thm 2). No minor order, no ordering among them, and none among the graphs that matches. `-1` and `-5` are integral because `3 - 2 = 9 - 8 = 1` (Gersonides); `-17` is sporadic (`139 | c_w`; wave-13 check). Their densities are best *upper* approximations of `c`, while `sigma_k`'s maximal cycles are best *lower* approximations (5/8 and 2/3 have mediant 7/11). | **NUMEROLOGY** (the must-break statement is PROVED) |
| Tutte's regular-matroid excluded minors {U2,4 self-dual; F7 ↔ F7*} ↔ the level-2 cube {Collatz ↔ 3n−1; `chi_(-4)`, `-chi_(-4)` fixed} | Both are `Z/2`-actions with one free orbit. The cube has **two** fixed points (bottom and top of the u-set square), the matroid list one. In Tutte's list all three are obstructions; in the cube the fixed points are the two provable corners and the free orbit the open/obstructed pair. | **ANALOGY** (orbit shape only) |
| Kuratowski-type characterization ↔ class (i) | Class (i) = the strategies containing no *expanding-cycle pattern* (the signs forced on the odd nodes of an expanding cycle of the union graph). This is THM-4474 Theorem A restated, so it is exact but tautological. Unlike Kuratowski's list, the minimal patterns grow with `k`: at least `N_k` node-disjoint ones at level `k` (THM-4479 Thm 2). | **PROVED** (restatement); "finite obstruction list" **false** |
| "dual pair + self-dual" (KT-b of the Kuratowski note) ↔ class (i) under `nu` | Class (i) = self-dual + pairs (Prop N, census §5), a literal instance of a `Z/2`-set. Nothing Kuratowski-like is attached. | **PROVED** (trivial) / KT reading **ANALOGY** |

## 7. What this does and does not say; failures

**What it does.**
* It makes the "tension" reading of the strategy cube exact. A periodic (or merely bounded) Lyapunov correction to `log n` exists iff the Farkas dual, a positive circulation, does not.
* The cost of the best correction is the maximum cycle mean. Collatz pays `log(3/2)` per step, because of the loop at `-1`, and no periodic correction and no finite bank of dyadic counters reduces this.
* It replaces the Christoffel picture by the correct one: positive-entropy maximizing sets, with the Christoffel orbit as the unique balanced member.
* It shows that provable strategies can be far more near-critical than `sigma_k` at the same level. `M_k` equals the numerator bound `G_(2^(k-2))` for `k <= 6`, and at level 7 a class-(i) strategy reaches `17/27` (denominator 27, margin `0.0014` per step), where `sigma_7` has `3/5` (margin `0.034`). Whether the denominators of `M_k` grow exponentially, as the numerator bound allows, is open.

**What it does not do.** Collatz is class (iv) at every level (THM-4474), and every statement above about Collatz is an obstruction to a *kind* of proof, not a statement about its orbits. No HYP or THM file is proposed. The Bernoulli and Kuratowski obstructions are recovered and slightly extended, not overturned.

**Failures and refuted intermediate guesses** (all checked):
* Q2 as posed (§3) and Q3's rigidity (§4) are refuted.
* **"`V_k` is the set of fractions in `[1/2, M_k]` with denominator `<= 2^(k-2) + 1`."** This fits `k <= 5`, but at `k = 6` the values `11/18`, `11/19`, `13/21`, `13/22` occur and `9/17` does not.
* **"A maximal cycle has at most `2^(k-3)` even nodes."** True for every strategy of levels 3–5 (65,808 checked), false at level 6: `13/22` has 9 even nodes.
* **"Class (i) is closed under turning u into d"** (a down-set in the u-lattice). False: 4 of 30 single u→d moves at level 4 and 380 of 4,794 at level 5 leave class (i).
* **"The symmetric core `U ∩ -U` of a class-(i) strategy is class (i)."** False for 80 of the 1,052 at level 5.
* **Unfinished solver runs.** The unrefined CP-SAT model ("a class-(i) cycle of density `> 17/27`" at level 7) ran about 30 minutes without a decision and was stopped. The refined model (sibling-free, length exactly 46, exactly 29 odd nodes; complete since `2·46 > 64`) stayed undecided after its 2-hour limit. At level 8, neither the refined query for `41/65` (30 min) nor the unrefined query for density `> 17/27` (15 min) decided. The flow-formulation MIP decides levels 4–5 in seconds, but did not decide "density `> 5/8`" at level 6 within 5 CPU-minutes.

## 8. Reproduction

```bash
python3 04-computation/experiments/procgen_tension_20260926_run.py > 05-knowledge/results/procgen_tension_20260926.out
```

* **Requirements.** `numpy`, `scipy` (HiGHS through `scipy.optimize.linprog` / `milp`), `ortools` (CP-SAT, 2 workers). No compiled engine, no caches.
* **Options.**
  * `--no-level6` skips the level-6 CP-SAT scan.
  * `--level7-bound` runs the refined level-7 query (time limit 2 h).
  * `--level8-search` runs `cpsat_refined` at level 8 for `41/65` (a complete query) and then `29/46` (a restricted search), 30 min each.
  * `--mip6` adds the level-6 query "density > 5/8" to the MIP cross-check.
* **The default run** re-derives everything else, including the level-7 witness, which is stored as a mask and re-verified.
* **Cost.** One process. Wall time 245 s (Q1 7 s, Q2 27 s, Q3 205 s including the level-6 CP-SAT scan, MIP 3 s, Q4 4 s); peak RSS 327 MB by the runner's own `getrusage` (343 MB by `/usr/bin/time -l`) (checked in the output to be below 700 MB).
* **Solver queries outside the default run** (same code, `procgen_tension_20260926_q3.py`): the level-7 query `cpsat_refined(7, 29/46)` (option `--level7-bound`) was run with the identical model from a scratch driver for its 2-hour limit; see §4.1 for its status. At level 8, `cpsat_refined(8, 41/65)` (30 min) and the unrefined `cpsat_query(8, 17/27)` (15 min) were run through drivers calling the committed functions, and both ended undecided (no witness above `17/27`). None of these runs is repeated by the default runner, and none of their outcomes is used as a proof.
* **Checks.** Every `ok:` line of the output is a `check(...)` that aborts the run on failure. The unlabelled indented lines print the computed tables quoted above.
* **Scratch.** Exploration scripts live in `scratch/procgen_tension/` and are not committed. The level-7 witness was first found there by the unrefined CP-SAT model (29 s); the runner only re-verifies it.
* **SHA-256** (raw bytes):
  * `procgen_tension_20260926_lib.py` `4f77ca9906e019273a37f5d1a06034b8805b248d2013258ae57928cee009fc3d`
  * `procgen_tension_20260926_q1.py` `3afb8dd826245eeecf86f7b78871b30bf04f92ed8a714cee7ec23d2f74ba7e24`
  * `procgen_tension_20260926_q2.py` `4c40f2bd9206f8f6644b4c2bb678302ec5ed2937be9ccf0d3ac6def4f8c8c518`
  * `procgen_tension_20260926_q3.py` `89eb4b5249d8025b2ffee40d15a3a56f45c7b3f81bc76ce65b86a164e7241dbd`
  * `procgen_tension_20260926_q4.py` `d2cb8a58f2d6d417c326486db5c03d9101defc5862612a2e222ff7f32bd15233`
  * `procgen_tension_20260926_mip.py` `b717e18ea089fcba87243c10ed6902fdf2fac975c26ea36f84cfe5c59f535f7d`
  * `procgen_tension_20260926_run.py` `1e87325490781d3bfbd805a98edaea2377247053afdf1a98f4b618221ad2767a`
  * `procgen_tension_20260926.out` `05de1d98418511a7736fde97ec8fa4f0632a62de173191c44afbd2472df030b5`

