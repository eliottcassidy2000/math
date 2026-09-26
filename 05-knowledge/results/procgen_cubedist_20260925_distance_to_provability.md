# Collatz's distance to provability in the strategy cube tends to zero at the exact exponential rate of the undecided density: flip the undecided residues (HYP-9138 proved)

**Status.**
- **PROVED** (hand proofs below; every lemma is also checked by code):
  - **Theorem 1 (HYP-9138).** For every `k >= 2`, let `sigma_k` agree with Collatz except `sigma_k = -` on `Bad_k`, the residues mod `2^k` on which Collatz has no descent within `k` steps. Then `sigma_k` is in class (i). Every cycle of its parity graph has odd density at most `F_k`, the best lower rational approximation of `log_3 2` with denominator `<= k`, and this is attained. Every integer `n > n_0(k)` descends below itself within `k` steps. Hence `delta_k <= |Bad_k| <= 2^(hk)`, and the Haar distance satisfies `delta_k / 2^(k-1) <= 2^(1-(1-h)k) -> 0`, with `h = h(log_3 2) = 0.9499555`. **HYP-9138 holds.**
  - **Theorem 2 (matching lower bound).** `delta_k >= N_k`, the number of binary necklaces of length `k` with more than `k log_3 2` ones, and `N_k >= 2^(hk)/(3k^2)`. So `log_2 delta_k = hk + O(log k)`: the exponent `1 - h = 0.0500445` of the Haar distance is **sharp**. This is the cube analogue of HYP-9137. When this lane finished, only the upper half of HYP-9137 was proved; the crossroads session has since proved it (THM-4478, 2026-09-26) by a different mechanism, an integer-capacity first-hit cut. The necklace packing here needs no capacity argument, because the cube is periodic.
  - **Potential certificates (Lemma P).** For rational `F = q/r < log_3 2`, "every cycle has odd density `<= F`" is equivalent to an integer potential `psi >= 0` with `psi(t) <= psi(s) - w_F(s)` on every edge, where `w_F = r - q` on odd and `-q` on even nodes. Then `phi = (log 3 / r) psi` satisfies the min-mean-cycle inequality `phi(t) <= phi(s) - (w(s) - mu)` with `mu = F log 3 - log 2 < 0`.
  - **Proposition 3 (SHEET).** Negation carries Collatz's neighbourhood onto 3n−1's, `sigma_k` onto the 3n−1 undecided construction, and preserves `delta_k`.
  - **Proposition 4 (DRIFT: stationary flip mass).** In the 5n±1 cube, every closed class of a class-(i) strategy has stationary mass at least `(1/2 - log_5 2)/k = 0.0693/k` on the residues where the strategy differs from 5n+1. There is no analogue for 3n+1.
- **FINITE-EXACT.**
  - Integer potential certificates for `sigma_k`, `k <= 20`, each checked edge by edge; exact Karp for `k <= 15`; `rho_max(sigma_k) = F_k`.
  - `sigma_k` is **transitive**, i.e. every positive orbit reaches 1, for every `k <= 300`. The threshold `n_0(k)` changes only at the denominators `2, 5, 8, 27, 46, 65, 149, 233` of the `F_k`.
  - `delta_k` for `k <= 9` re-derived by an independent solver: `1, 2, 2, 4, 5, 9, 14, 23`.
  - `40 <= delta_10 <= 44`, with the lower bound certified by a stored, re-derivable no-good set. The exact value is open.
  - `52 <= delta_11 <= 72`.
  - DRIFT: the 5n±1 provable class is **empty** at levels 2–6 (UNSAT no-good certificates, re-checked by a SAT solver; exhaustive for `k <= 5`) and **nonempty** from level 7 on. Collatz's analogue there, 5n+1, is at exact distance `29/64 = 0.453` at `k = 7`.
- **EMPIRICAL.**
  - Greedy pruning inside `Bad_k` gives class-(i) sets of size `44, 72, 131, 233, 398, ...` at `k = 10, 11, 12, 13, 14, ...` (§5).
  - For 5n+1, pruned class-(i) sets have Haar fraction about `3.2/k` for `7 <= k <= 14`, against a proved necklace lower bound of about `2/k`.
- **OPEN.**
  - The constant: `N_k <= delta_k <= |Bad_k|`, and the ratio `|Bad_k|/N_k` grows slowly (`2.4` at `k = 9`, `7.5` at `k = 200`). Where `delta_k` sits between them is open.
  - For 5n+1: whether its distance to class (i) tends to 0. The data suggest `≍ 1/k`, i.e. the drift changes the rate from exponential to polynomial.
- **REFUTED.** Nothing.
- No HYP or THM file was created. Collatz itself is untouched (§7).

Session `collatz-procgen-20260922`, cube-distance lane, 2026-09-25. Scripts `04-computation/experiments/procgen_cubedist_20260925_{lib,bad,exact,prune,controls,run}.py` and the C engine `procgen_cubedist_20260925_engine.c`. Output [procgen_cubedist_20260925.out](procgen_cubedist_20260925.out). Parents: [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) and its [note](procgen_cube_20260925_strategy_cube.md); [HYP-9138](../hypotheses/HYP-9138-cube-distance-to-provability.md); [THM-4475](../../01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md) and its [note](procgen_price_20260925_provability_price.md).

## 0. The answer in brief

**Setting** (THM-4474). A level-`k` strategy `sigma` puts a sign on each odd residue mod `2^k`:
* `T_sigma(n) = n/2` for even `n`, and `(3n + sigma(n mod 2^k))/2` for odd `n`.
* The parity graph `G_sigma` has nodes `Z/2^k` and edges from `s` to the two lifts of `T(s) mod 2^(k-1)`. Node weights are `w = log(3/2)` (odd) and `-log 2` (even).
* **Class (i)**, bounded-lookahead provability, holds iff every cycle of `G_sigma` has odd density `< c = log_3 2` (Theorem A).
* `delta_k` is the least number of residues at which a class-(i) strategy differs from Collatz (`sigma = +`). HYP-9138 asks whether `delta_k / 2^(k-1) -> 0`.

**The construction.** Flip exactly the residues on which Collatz is undecided for `k` steps:

`sigma_k(r) = -`  iff  `r mod 2^k` lies in `Bad_k = { r : 3^(a_j(r)) > 2^j for j = 1..k }`.

This is THM-4475's A-rescue ("flip the undecided number itself") frozen into residues. It works for two reasons.
* **The flip is itself a descent.** An undecided residue is `3 (mod 4)`, so `(3r-1)/2` is even and `T^2 = (3r-1)/4`: two steps with multiplier `3/4`.
* **Decided residues never meet a flip before deciding** (heredity, Lemma 1.1). If `x` is decided, its Collatz orbit descends at some `d(x) <= k`. If it met an undecided point first, the orbit would stay up for `k` more steps, which is impossible.

So every orbit of `T_(sigma_k)` is a concatenation of blocks of multiplier `< 1` and length `<= k`. The rewiring caused by the flips, the obstacle named in HYP-9138, is never seen inside a block.

**Sharpness.** In parity-word coordinates Collatz's parity graph is the de Bruijn graph `B(2,k)`. The rotation classes of `k`-words (necklaces) are node-disjoint closed walks, and each necklace with more than `ck` ones is an expanding cycle. Any class-(i) strategy must change a sign on each of them, so `delta_k >= N_k`.

| `k` | `N_k` (Thm 2, lower) | `delta_k` (exact, or bounds) | pruned (upper, §5) | `\|Bad_k\|` (Thm 1, upper) | Haar `delta_k/2^(k-1)` | Haar of `sigma_k` |
|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 1 | 0.500 | 0.500 |
| 3 | 2 | 2 | 2 | 2 | 0.500 | 0.500 |
| 4 | 2 | 2 | 2 | 3 | 0.250 | 0.375 |
| 5 | 2 | 4 | 4 | 4 | 0.250 | 0.250 |
| 6 | 5 | 5 | 5 | 8 | 0.156 | 0.250 |
| 7 | 5 | 9 | 9 | 13 | 0.141 | 0.203 |
| 8 | 6 | 14 | 15 | 19 | 0.109 | 0.148 |
| 9 | 16 | 23 | 25 | 38 | 0.090 | 0.148 |
| 10 | 19 | 40–44 | 44 | 64 | 0.078–0.086 | 0.125 |
| 11 | 52 | 52–72 | 72 | 128 | 0.102–0.141 | 0.125 |
| 12 | 70 | | 131 | 226 | | 0.110 |
| 13 | 85 | | 233 | 367 | | 0.090 |
| 14 | 251 | | 398 | 734 | | 0.090 |
| 16 | 434 | | | 2,114 | | 0.065 |
| 20 | 6,910 | | | 27,328 | | 0.052 |
| 40 | 1.11e9 | | | 6.40e9 | | 0.0117 |
| 100 | 4.21e25 | | | 3.03e26 | | 4.8e-4 |
| 200 | 6.60e53 | | | 4.92e54 | | 6.1e-6 |

* Both bounds are `2^(hk) k^(-3/2)` up to slowly varying factors: `N_k k^1.5 / 2^(hk)` oscillates in `[0.70, 1.65]`, and `|Bad_k| k^1.5 / 2^(hk)` creeps from 2.2 to 8.9 over `8 <= k <= 200`.
* The apparent `2^(-k/3)` decay reported by the cube lane for `k <= 9` is mostly the `k^(-3/2)` factor, whose local slope `-1.5/(k ln 2)` is about `-0.3` per level at `k = 7`. The exponential rate is only `2^(-0.05k)`.

## 1. Setting, and potentials as certificates

Notation follows THM-4474 §1. For a closed walk with `a` odd nodes and length `p`, the weight is `a log 3 - p log 2`. It is expanding iff `3^a > 2^p`, and equality never occurs.

**Lemma P (potential certificates).** Let `F = q/r` with `0 < q < r` integers, and put `w_F(s) = r - q` for odd `s` and `-q` for even `s`. For any strategy `sigma` the following are equivalent:
* (a) every cycle of `G_sigma` has odd density `<= F`;
* (b) there is an integer `psi : Z/2^k -> Z_(>=0)` with `psi(t) <= psi(s) - w_F(s)` on every edge `s -> t`.

The least such `psi` is `psi(s) = sup` over walks from `s` of the partial sums of `w_F`, computed by value iteration.
* If `F < log_3 2`, then (b) certifies class (i).
* Moreover `phi = (log 3 / r) psi` satisfies the form `phi(t) <= phi(s) - (w(s) - mu)` with `mu = F log 3 - log 2 < 0`.
* Conversely, a class-(i) strategy satisfies (a) with `F = rho_max(G_sigma) < log_3 2`, so such a certificate always exists. This is the min-mean-cycle (Karp) duality in integer form.

*Proof.*
* (b) ⇒ (a). Sum the inequality around a cycle with `a` odd nodes and length `p`: `0 <= -(a(r-q) - (p-a)q) = pq - ar`.
* (a) ⇒ (b). Every cycle has `w_F`-weight `ar - pq <= 0`, so the supremum of partial sums over walks from `s` is attained on a simple path and is finite. Prepending an edge gives `psi(s) >= w_F(s) + psi(t)`, and the empty walk gives `psi >= 0`.
* The rescaling. `(log 3/r)(r - q) = log 3 - F log 3 = log(3/2) - mu`, and `-(log 3/r) q = -log 2 - mu`. So `(log 3/r) w_F = w - mu` exactly. ∎

* The engine computes `psi` by value iteration on the reversed graph. Every certificate is re-checked in Python on all `2^(k+1)` edges.
* For `F_k = 5/8` the margin is `mu = -0.00652` per step. It shrinks to `-0.00143` (`17/27`), `-0.00054` (`29/46`), `-0.00018` (`41/65`), and so on: the provable strategies near Collatz are near-critical.

## 2. Theorem 1: flipping the undecided residues

**Notation.**
* For `x in Z_2` let `a_j(x)` be the number of odd terms among `x, Tx, ..., T^(j-1)x` (Collatz `T`), and `M_j(x) = 3^(a_j(x))/2^j`.
* The first `k` parities of `x` depend only on `x mod 2^k` (Terras). So `Bad_k = { r mod 2^k : M_j(r) > 1, 1 <= j <= k }` is a set of residues, the classes on which Collatz has no `k`-step descent (choice ladder, price note).
* For `x` outside `Bad_k`, let `d(x) = min{ j <= k : M_j(x) < 1 }` (never `M_j = 1`).
* For `k >= 2`, `Bad_k ⊂ {3 mod 4}`: `M_2 > 1` forces `x` and `Tx` odd.
* `sigma_k = -` on `Bad_k` and `+` elsewhere; write `T_k = T_(sigma_k)`.

**Lemma 1.1 (heredity).** If `x` is outside `Bad_k`, then `T^i x` is outside `Bad_k` for `0 <= i < d(x)`.

*Proof.* Suppose `T^i x` lies in `Bad_k` with `1 <= i < d(x)`.
* Then `M_i(x) > 1`, since `i < d(x)`.
* Also `M_(i+t)(x) = M_i(x) M_t(T^i x) > 1` for `1 <= t <= k`.
* So `M_j(x) > 1` for all `1 <= j <= i + k`, contradicting `M_(d(x))(x) < 1` with `d(x) <= k`. ∎

**Lemma 1.2 (blocks).** Every `x in Z_2` falls in exactly one case.
* **(F)** `x` in `Bad_k`. Then `T_k x = (3x-1)/2` is even and `T_k^2 x = (3x-1)/4`: two steps of multipliers `3/2` and `1/2`, product `3/4`.
* **(C)** `x` outside `Bad_k`. Then `T_k^j x = T^j x` for `0 <= j <= d(x)`, with `M_j(x) > 1` for `1 <= j < d(x)` and `M_(d(x))(x) < 1`. The parity word of the block is a *first-descent word* of length `d(x) <= k`.

*Proof.*
* (F): `x = 3 (mod 4)`, so `3x - 1 = 0 (mod 4)`.
* (C): by Lemma 1.1, `sigma_k = +` at every odd `T^i x` with `i < d(x)`. ∎

**Theorem 1.** Let `k >= 2`.
* **(a)** Every `T_k`-orbit in `Z_2` is a concatenation of (F)- and (C)-blocks, each of multiplier `< 1` and length `<= k`. Hence `3^(a_j)/2^j <= (3/2)^(k-1)` along every orbit and at every time, and `sigma_k` is in class (i).
* **(b)** `rho_max(G_(sigma_k)) = F_k := max{ a/d : 1 <= d <= k, 3^a < 2^d }`, the best lower approximation of `log_3 2` with denominator `<= k`.
* **(c)** Every integer `n > n_0(k)` satisfies `T_k^j(n) < n` for some `j <= k`, with `j = 2` on `Bad_k`. Here `n_0(k) = max floor(c_w/(2^d - 3^a))` over first-descent words `w` of length `d <= k`, and `T^d(n) = (3^a n + c_w)/2^d` on the class of `w`. In Theorem A's terms, `Bad_k(sigma_k)` is empty and the lookahead is `L_min <= k`.
* **(d)** `delta_k <= |Bad_k| <= 2^(hk)`, so `delta_k/2^(k-1) <= 2^(1-(1-h)k) -> 0`. **HYP-9138 holds.**

*Proof.*
* **(a)** By Lemma 1.2 applied at each block start. A block that is still in progress at time `j` has multiplier at most `(3/2)^(k-1)`, and the completed blocks multiply to less than 1. Bounded walk weights exclude expanding cycles.
* **(b), upper bound.**
  * An (F)-block has density `1/2 <= F_k`; a (C)-block has density `a/d` with `3^a < 2^d`, `d <= k`.
  * A cycle `gamma` of `G_(sigma_k)` carries a periodic point `x_gamma` of period `p = |gamma|` (THM-4474 Lemma 2). Block starts are determined by the current point, so two of them coincide within the orbit.
  * The blocks between them cover a whole number of periods, and `a/p` is a mediant of block densities.
* **(b), equality.** Let `F_k = a/d`.
  * The word with prefix counts `a_j = floor(c j) + 1` (`1 <= j < d`) and `a_d = a` is a first-descent word. Its increments are 0 or 1 because `a/(d-1) > c` (maximality of `F_k`) and `(a-1)/(d-1) < c`.
  * Its Collatz periodic orbit descends within `d <= k` steps from every point, so it avoids `Bad_k`. It is therefore a `T_k`-orbit, and it yields a cycle of density exactly `F_k`.
* **(c)** On `Bad_k`, `(3n-1)/4 < n`. Off `Bad_k`, `T^d(n) < n` iff `n(2^d - 3^a) > c_w`. Even `n` halve.
* **(d)** `sigma_k` differs from Collatz exactly on `Bad_k`. And `Bad_k ⊂ {a_k > ck}`, with `sum_(j > ck) C(k,j) <= 2^(h(c)k)` (Chernoff). ∎

**Remarks.**
* **Why the rewiring is harmless.** Flipping `r` moves its target from `T(r)` to `T(r) - 1`, a 2-adic jump that rewires `G_sigma` globally. But an orbit meets a flipped residue only at the start of a block: blocks of type (C) avoid `Bad_k` by Lemma 1.1. At a flipped residue the orbit takes a `3/4` descent before anything else. No block "sees" the rewiring.
* **Comparison with THM-4475.**
  * `sigma_k` is the A-rescue alone.
  * In the pairing family the A-flip of `n` pushes its partner up, which forces the F/B rescues and the `-5 (mod 2^9)` analysis. A sign flip has no partner.
  * The pair-0 obstruction (no periodic provable pairing) has no analogue. The two fixed-point conditions of Proposition F are automatic: `-1` lies in `Bad_k` and `1` does not.
  * `sigma_k` also breaks the `-5` cycle: its word `(110)^inf` is ballot at every length, so `-5` lies in `Bad_k`.
* **Near-criticality.** `rho_max(sigma_k) = F_k` tends to `log_3 2`, and the certificate margin `mu_k = F_k log 3 - log 2` tends to 0. The exact optima of §4 also sit at `rho_max = 5/8` for `5 <= k <= 9`.

**FINITE-EXACT checks** (output §A).
* `|Bad_k|` by enumeration equals the ballot DP (`k <= 20`); every element is `3 (mod 4)`.
* The block lemma holds on the parity graph for all lifts (`k <= 13`, 0 violations).
* Class (i) of `sigma_k`:
  * an integer certificate at threshold `F_k`, checked on every edge, for `2 <= k <= 20` (`max psi` grows from 1 to 36; the bound is `(k-1)(r-q)`);
  * exact Karp confirms `rho_max = F_k` for `k <= 15`. Here `F_k = 1/2` for `k <= 4`, `3/5` for `5 <= k <= 7`, `5/8` for `8 <= k <= 26`, `17/27` for `27 <= k <= 45`, and so on.
* **Transitivity.** `n_0(k)` is computed exactly by a DP over first-descent words and agrees with enumeration for `k <= 20`. It takes the values
  * `1` (`k >= 2`), `4` (`k >= 5`), `24` (`k >= 8`), `108` (`k >= 27`), `281` (`k >= 46`), `867` (`k >= 65`), `2419` (`k >= 149`), `4862` (`k >= 233`),
  
  changing exactly at the denominators of the `F_k`. For every `k <= 300`, every `n <= n_0(k)` reaches 1 under `T_k`. With (c) and strong induction, **every positive orbit of `T_(sigma_k)` reaches 1, for every `k <= 300`**.
  * The descent within `k` steps is also spot-checked on `(n_0, n_0 + 20000]` for `k <= 20`.

## 3. Theorem 2: the matching lower bound

**De Bruijn coordinates.** By the Terras bijection, residues mod `2^k` correspond to parity words of length `k`. Collatz's `G_0` is the de Bruijn graph `B(2,k)`: `v_0...v_(k-1) -> v_1...v_(k-1) b`, and a node is odd iff its first letter is 1.

**Theorem 2.** `delta_k >= N_k`, the number of binary necklaces of length `k` with more than `k log_3 2` ones. Moreover `N_k >= (1/k) C(k, ceil(ck))`, and `N_k >= 2^(hk)/(3k^2)` for every `k >= 2`.

*Proof.*
* For a `k`-word `v`, the rotations `v -> rot(v) -> ... -> v` form a closed walk of length `k` in `G_0`, since `rot(v)` is a successor of `v`. Its odd density is the density of `v`, and it is expanding iff `3^(|v|_1) > 2^k`.
* Distinct necklaces give node-disjoint closed walks.
* If `sigma` agrees with `+` at every odd node of such a walk, the walk survives in `G_sigma`: edges out of even nodes never change. So a class-(i) `sigma` changes a sign at an odd node of every expanding necklace, and these residues are distinct.
* The bound: each necklace has at most `k` words, and `C(k,m) >= 2^(k h(m/k))/(k+1)`. Moreover `h(ceil(ck)/k) >= h(c) - 1.442/k` for `k >= 10`, since `|h'| <= log_2(0.7309/0.2691)` on `[c, c + 0.1]`. This gives `N_k >= 2^(hk)/(2.717 k(k+1)) >= 2^(hk)/(3k^2)` for `k >= 10`; `k <= 9` is checked directly. ∎

**Corollary (sharp exponent).** `(1/k) log_2(delta_k / 2^(k-1)) -> -(1-h) = -0.0500445`. Theorems 1 and 2 give `2^(1-(1-h)k)/(3k^2) <= delta_k/2^(k-1) <= 2^(1-(1-h)k)`.

**Remarks.**
* The pairing family has the same exponent in its upper bound (THM-4475). Its lower bound was `0.774` when this lane ran (then `0.1445` by THM-4477); HYP-9137 is now proved by THM-4478 (crossroads, 2026-09-26), with the sharp exponent `1-h`, via growth bands and actual-integer capacity. In the cube, periodicity supplies a finite graph with disjoint necklace cycles, so the lower bound is elementary.
* Theorem 2 uses nothing about the multiplier; it is reused for 5n+1 in §6.
* **Orders of magnitude.**
  * `N_k ≍ k^(-3/2) 2^(hk)` is elementary (binomial tail), with a factor oscillating with `frac(ck)` in `[0.70, 1.65]` (§A, B4).
  * For `|Bad_k| = 2^k P(tau > k)`, the first-passage time of a non-lattice walk with steps `log(3/2), -log 2`, the classical asymptotics `P(tau > n) ~ C n^(-3/2) gamma^n` (Iglehart 1974) would give `|Bad_k| ≍ k^(-3/2) 2^(hk)`, hence `delta_k ≍ k^(-3/2) 2^(hk)`. This is not re-proved here; the data (`2.2 -> 8.9`) are consistent with slow convergence of the constant.

## 4. Exact values and the structure of the optima

**Method** (`..._exact.py`): an implicit hitting set, the cube lane's lazy MaxSAT with a MIP solver.
* A cycle of `G_sigma` survives whenever `sigma` is unchanged on its odd nodes, so each expanding cycle met is a valid no-good.
* The hitting-set optimum over the no-goods found so far (HiGHS, integer data) is a lower bound. The first optimum without expanding cycles is optimal.
* No-goods are harvested from the optimum by greedy repair phases (node-disjoint cycles from the engine, plus all expanding cycles of length `<= 14` through changed nodes). The seeds are all expanding de Bruijn cycles of length `<= k+3`.
* Every no-good is stored as (cycle, signs on its odd nodes) and is re-derivable on its own. Every optimum is re-verified by exact Karp.

**Results.**
* `delta_k = 1, 2, 2, 4, 5, 9, 14, 23` for `k = 2..9`, equal to the cube lane's MaxSAT values (independent code and solver). `k = 9` takes about a minute.
* **`40 <= delta_10 <= 44`** (Haar `0.078–0.086`).
  * The lower bound is certified by a stored no-good set, [procgen_cubedist_20260925_k10_certificate.json.gz](procgen_cubedist_20260925_k10_certificate.json.gz): 1949 no-goods, plus the 241 regenerated seeds.
  * The runner re-derives every no-good as a genuine expanding cycle and re-proves that no hitting set of size `<= 39` exists (105 s).
  * The upper bound is the pruned set of §5, re-verified class (i).
  * Closing the gap needs more IHS rounds: each MIP round took 2–6 minutes at the end, and the run was stopped at LB 40. A 2-swap local search and a potential-based MIP (threshold `5/8`, 25 min) found nothing below 44.
* **`52 = N_11 <= delta_11 <= 72`** (Theorem 2 and the pruned set). No IHS run was attempted at `k = 11`.

**Structure** (output §B; "necklace profile" = how many expanding necklaces carry `m` flipped rotations).
* **What every optimum contains.** Every optimum found (`k <= 9`) and the best `k = 10` set flip `-1` (forced, Proposition F). For `k >= 3` they also flip `-5` (not forced). Each carries at least one flip on every expanding necklace (forced, Theorem 2).
* **Where the flips sit.**
  * Most expanding necklaces carry one or two flips. The profiles are `{1: 5}` at `k = 6`, `{1: 9, 2: 7}` at `k = 9`, and `{1: 3, 2: 9, 3: 5, 4: 2}` for the `k = 10` set.
  * Some flips sit on non-expanding necklaces: 5 of the 14 at `k = 8`.
  * Some lie outside `Bad_k`: 5 of 14 at `k = 8`, 5 of 23 at `k = 9`.
  * Almost all flips are `u -> d` (residues `3 mod 4`). The exception is one `d -> u` flip (`9`) in the `k = 5` optimum.
* **So `Bad_k` is a proof device, not the shape of the optimum.** The optimum uses about `1.4–2.3` flips per expanding necklace (`delta_k/N_k = 1, 1, 1, 2, 1, 1.8, 2.3, 1.4` for `k = 2..9`), against `|Bad_k|/N_k = 1–3.4` over the same range.

## 5. Smaller explicit sets, and what fails

* **Greedy pruning inside `Bad_k`** (restore `sigma = +` wherever class (i) survives; exact test by certificate; final sets re-verified by an edge-checked certificate and by Karp). The resulting sizes are
  * `2, 2, 4, 5, 9, 15, 25, 44, 72, 131, 233, 398` for `k = 3..14` (six pruning orders each).
  * These are optimal for `k <= 7`, and within `1, 2` of `delta_k` at `k = 8, 9`. All have `rho_max = 5/8` for `k >= 6`.
  * A pure exponential fit over `8 <= k <= 14` gives `2^(-0.207 k)`. A fit `a + b k + c log2 k` gives `b = -0.060`, `c = -1.10`, close to the theoretical `-(1-h) = -0.050` and `-3/2`. For `sigma_k` itself, over `8 <= k <= 200`, the same fit gives `b = -0.054`, `c = -0.98`.
* **One flip per expanding necklace fails.** Flip only the lowest ballot rotation, or only the highest, of each expanding necklace, so that `|R| = N_k`.
  * Both rules leave an expanding cycle for `k = 5` and for every `7 <= k <= 14`.
  * The lowest rotation works only at `k = 3, 6`; the highest only at `k = 3, 4`.
  * So the lower bound is not attained by the obvious choice, and indeed `delta_k > N_k` at `k = 5, 7, 8, 9`.
* **Flipping only the rising runs near `-1` and `-5` fails.** These are the 2-adic neighbourhoods of the two smallest expanding cycles, the natural first construction. Checked at `k = 8, 10, 12`:
  * `R_j = {r = -1 mod 2^j}` is class (i) only for `j = 2`, where `R_2` is all residues `3 mod 4`, i.e. the all-d strategy (Haar `1/2`). The words `(1^(j-1) 0)^inf` avoid `R_j` and have density `(j-1)/j > c`.
  * Adding the `-5` neighbourhood, `R'_j = R_j ∪ {r = -5 mod 2^j}`, works for `j <= 4`: `R'_3 = R_2`, and `R'_4 = {11, 15 mod 16}` is the lifted level-4 optimum (Haar `1/4`). It fails from `j = 5` on.
  * So the neighbourhoods of `-1` and `-5` alone stop working at depth 5. The undecided set, which follows every expanding orbit and not only the two smallest ones, is what works.
* **The optima are not subsets of `Bad_k`** (§4).
  * 5 of the 14 flips lie outside `Bad_8`, and 5 of 23 outside `Bad_9`.
  * `d -> u` flips (`1 mod 4`) also occur: `9` in the `k = 5` optimum here, and `121` and `249` in the cube lane's optima at `k = 7, 8`.

## 6. Controls

**SHEET** (Proposition 3; Theorem E of THM-4474).
* `nu(x) = -x` gives `T_(nu sigma) = nu T_sigma nu` with `(nu sigma)(m) = -sigma(-m)`. It maps Collatz with changes on `R` to 3n−1 with changes on `-R`, and preserves class (i).
* So the distance of 3n−1 to class (i) is `delta_k`. The 3n−1 undecided set is `-Bad_k`, and the image of `sigma_k` is "3n−1 with `+` on its undecided residues", a class-(i) strategy.
* Verified for `k <= 9`:
  * independent IHS runs started from the all-minus strategy, without seeds, return the same `delta_k`;
  * the `nu`-image of each Collatz optimum is class (i) with the same `rho_max`;
  * `-Bad_k(3n+1) = Bad_k(3n-1)`, and its flip is class (i).

**DRIFT: the 5n±1 cube.** `T(n) = n/2` or `(5n + sigma)/2`. A cycle is expanding iff `5^a > 2^p`, i.e. its odd density exceeds `log_5 2 = 0.430677`. The Collatz analogue is 5n+1 (`sigma = +`), whose drift is positive.
* Theorem A of THM-4474 holds verbatim with threshold `log_5 2`: its proof uses only that the multiplier is odd (Terras bijection, Banach periodic points, cycle lemma).
* **The provable class is empty at levels 2–6.** At each such level the expanding-cycle no-goods admit no hitting set.
  * The MIP is infeasible, and the UNSAT is re-confirmed by the SAT solver Glucose4 on the re-derived no-goods.
  * For `k <= 5` it is also confirmed by exhaustive Karp over all `4, 16, 256, 65536` strategies.
* **It is nonempty at level 7, hence at every level `>= 7`** (a lift has the same map). The first provable strategies are near-critical: `rho_max = 3/7 = 0.4286`, just below `log_5 2`.
* **5n+1's exact distance at `k = 7` is `29 of 64` residues (Haar `0.453`).** At `k = 8` the IHS did not finish in 20 minutes; the pruned upper bound is `50 of 128` (`0.391`).
  * The optimum changes 14 of the 16 residues `5 (mod 8)`, 9 of `1 (mod 8)`, 6 of `7 (mod 8)` and none of `3 (mod 8)`.
* **Why Theorem 1 fails for 5n+1.**
  * The undecided set `Bad_k(5n+1)` has positive density, since 5n+1 has positive drift: `|Bad_k(5)|/2^k = 0.26, 0.22, 0.196, 0.178` at `k = 10, 20, 40, 160`, tending to about `0.176` as in the price note.
  * The flip block is not a descent: at `r = 1 (mod 8)`, `(5r-1)/4` has multiplier `5/4 > 1`.
* **Proposition 4 (stationary flip mass).** Let `sigma` be a class-(i) strategy of the 5n±1 cube and `R` the residues where it differs from 5n+1. Then every stationary law `pi` of the uniform-lift chain on a closed class satisfies `pi(R) >= (1/2 - log_5 2)/k = 0.0693/k`.

  *Proof.*
  * Let `P_0` be the 5n+1 chain. `P_0^k(s,·)` is uniform (`T^k` maps `s + 2^k Z_2` affinely onto `Z_2`), and `P_sigma - P_0 = Delta` is supported on the rows of `R`.
  * From `pi P_sigma = pi`: `pi - U = sum_(j<k) (pi Delta) P_0^j`. Hence `||pi - U||_1 <= 2k pi(R)` and `|pi(odd) - 1/2| <= k pi(R)`.
  * The sandwich (THM-4474 Theorem C, whose proof is multiplier-free) gives `pi(odd) <= rho_max < log_5 2`. ∎
  
  For 3n+1 the same computation gives nothing, since `1/2 < log_3 2`. Positive drift is exactly what forces the orbits to meet changed residues at frequency `≳ 1/k`.
* **Data** (output §D2; lower bound = Theorem 2's necklaces for multiplier 5; upper = greedy pruning of the lifted level-7 optimum, every set certificate-checked; Proposition 4 checked on every closed class for `k <= 11`):

| `k` | `N_k(5)` (lower) | Haar | pruned (upper) | Haar | `k` × Haar (upper) | min `pi(R)` over closed classes (bound `0.0693/k`) |
|---|---|---|---|---|---|---|
| 7 | 10 | 0.156 | 29 (exact) | 0.453 | 3.17 | 0.147 (≥ 0.0099) |
| 8 | 23 | 0.180 | 50 | 0.391 | 3.13 | 0.143 (≥ 0.0087) |
| 9 | 44 | 0.172 | 94 | 0.367 | 3.31 | 0.136 (≥ 0.0077) |
| 10 | 67 | 0.131 | 158 | 0.309 | 3.09 | 0.125 (≥ 0.0069) |
| 11 | 136 | 0.133 | 297 | 0.290 | 3.19 | 0.119 (≥ 0.0063) |
| 12 | 216 | 0.105 | 534 | 0.261 | 3.13 | — |
| 13 | 448 | 0.109 | 1005 | 0.245 | 3.19 | — |
| 14 | 714 | 0.087 | 1901 | 0.232 | 3.25 | — |

* **Reading.**
  * For 5n+1 the necklace lower bound is `N_k(5)/2^(k-1) ~ 2/k`, since the fraction of words with density `> 0.43` tends to 1; it is `1.1/k–1.5/k` at `k <= 14`. This is only polynomially small, so **no positive floor can come from cycle packing**.
  * On the data, the distance of 5n+1 to provability lies between `N_k(5)/2^(k-1)`, i.e. `1.1/k–1.5/k`, and about `3.2/k` for `7 <= k <= 14`.
  * Whether it tends to 0 is OPEN. If it does, it does so polynomially, against the exponential `2^(-0.05k)` of 3n+1.
  * The drift changes the **rate**. It does not visibly change the limit: both neighbourhoods may accumulate on the provable class, 3n+1's exponentially fast.

## 7. What this does and does not say

* **What it does.**
  * HYP-9138 is proved with an explicit, one-line construction and a two-line lemma. The exponent `1 - h(log_3 2)` of the Haar distance is proved sharp.
  * In the strategy cube, Collatz lies in the Haar closure of the provable class, i.e. of the strategies with a bounded-lookahead descent proof.
  * The approximants `sigma_k` are transitive (every positive orbit reaches 1) for every `k <= 300`, and they certify it by `k`-step descent above a threshold that grows only at the convergent denominators.
* **What it does not do.**
  * Nothing here bears on Collatz itself. `sigma_k` differs from Collatz exactly on the undecided classes, i.e. exactly where Collatz's behaviour is unknown (DEFECT, as in THM-4475 §5).
  * Collatz stays in (iv) at every level (THM-4474).
* **Place in the session.** The same exponent `1 - h` is now proved from both sides in all three settings:
  * the pairing family: THM-4475 from above, THM-4478 from below;
  * arbitrary fixed-horizon edits: THM-4478;
  * the strategy cube: this note.
  The three lower bounds use three different mechanisms: integer capacity (THM-4478), second moments (THM-4477, which stops at `2(1-h)`), and necklace packing (here).

## 8. Reproduction

```bash
python3 04-computation/experiments/procgen_cubedist_20260925_run.py > 05-knowledge/results/procgen_cubedist_20260925.out
```

* **Engine.** The runner compiles `procgen_cubedist_20260925_engine.c` into `scratch/procgen_cubedist/` (caches there are not to be committed). It needs `numpy`, `scipy`, `highspy` (the HiGHS MIP) and `pysat` (Glucose4).
* **Long runs.** The `k = 10` lower bound comes from a long IHS run. The runner re-verifies its stored certificate from scratch: every no-good is re-derived as a genuine expanding cycle, the hitting set below the bound is re-proved infeasible, and the upper-bound set is re-checked by Karp.
  * The `k = 10` IHS run is checkpointed and resumable:
    ```bash
    python3 04-computation/experiments/procgen_cubedist_20260925_exact.py 10 ckpt short=14 ub=prune time=5400
    python3 04-computation/experiments/procgen_cubedist_20260925_exact.py export 10
    ```
  * The runner does not repeat it. It re-verifies the stored certificate `procgen_cubedist_20260925_k10_certificate.json.gz` (27 KB).
* **Cost.** One process at a time. Wall time 610 s. Peak RSS of the runner 424 MB; the HiGHS MIP uses 2 threads; engine subprocesses are small. The `k = 10` IHS run took about 1.5 hours in total and was stopped at lower bound 40.
* **Checks.** Every claim in the output is a `check(...)` that raises on failure.
* **SHA-256** (raw bytes):
  * `procgen_cubedist_20260925_lib.py` `3ea102495a2265a2e4e17c4cd0d58b700b4d444d01bc427ca788f72c00e050ec`
  * `procgen_cubedist_20260925_bad.py` `f088c44a23e31a27aa07725eba9f312cad465ad8d901985ecdb64834844bd74e`
  * `procgen_cubedist_20260925_exact.py` `0885737724aaa98b01a98153d5420cd6581762b264d367b506c76bad9f036dd1`
  * `procgen_cubedist_20260925_prune.py` `945bbc0ab6608757ab63bd8d50f2a88ee0e9973150d021e849374168a890b788`
  * `procgen_cubedist_20260925_controls.py` `d88bd92564df7d37234226a80d740943d24cc430e80413b2c5330de3a54a5e64`
  * `procgen_cubedist_20260925_run.py` `5d556f152569e5f3fc33d0d89157626a7a674674f62da92087f613b65f3335d0`
  * `procgen_cubedist_20260925_engine.c` `9b1f2d4d6b927a10143e696b171ed8a33c21722e5178a0ce256f2f0d1a9db054`
  * `procgen_cubedist_20260925.out` `4596290838f30a09589b1c1d0072004793bec4ff10f9c8a13494ff74ac2a2dd5`
  * `procgen_cubedist_20260925_k10_certificate.json.gz` `b351d6dccf135c585f9a3fdd00a5c8fdb7a4e7082e868e08094a93de3e3ce198`
