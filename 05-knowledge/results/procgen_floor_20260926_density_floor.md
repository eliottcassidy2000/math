# The floor of the strategy cube: the min-max cycle density is the value of a mean-payoff game; exact values to level 22, q = 5 settled at every level, a negative-integer adversary proving rho_max >= log_(q+1) 2, and the limits of entropy-type floors

**Status.**
- **PROVED** (hand proofs in §§2-7; the runner checks their finite content):
  - **Theorem G (game form of rho\*).** For odd `q` and level `k`, `rho*(q,k)` (the least `rho_max` over level-`k` sign strategies) is the value of a finite mean-payoff game: Min fixes the signs, Max fixes a lift bit at every residue pair. Two elementary certificate lemmas turn every computed value into a checkable fact:
    - an **upper certificate** is a strategy with an integer potential (THM-4474's Lemma P);
    - a **lower certificate** is a Max lift strategy `tau` with a potential on a closed node set. It proves that **every** sign strategy has a cycle of density `>= F`. It is a compact, `2^(k-1)`-bit alternative to the drift lane's UNSAT sets of cycle no-goods.
    - (CITED: Ehrenfeucht–Mycielski positional determinacy, as used by energy-game algorithms.) Matching certificates always exist, so `rho*(q,k) = max_tau min_sigma`. The values below do not depend on this citation.
  - **Proposition S (the q-dependence).**
    - (a) The level-`k` arena depends only on `q mod 2^(k-1)`: swapping the two lifts of every odd pair is an isomorphism `arena(q) ≅ arena(q + 2^(k-1))`.
    - (b) `rho*(-q,k) = rho*(q,k)`. The least fixed point of Min's operator is negation-invariant, and `(-sigma*, f*)` certifies the same value for `-q`. This is not an arena isomorphism: the full `rho_max` distributions of `q` and `-q` differ.
    - (c) For `k >= 3` there is no arena isomorphism between `q` and `q'` unless `q' ≡ q (mod 2^(k-1))`.
  - **Theorem N (the negative-integer adversary).** For every odd `q >= 3`, every level `k >= 2` and every sign strategy, `G_sigma` has a cycle of odd density `>= log_(q+1) 2`.
    - Max always takes the lift with top bit 1, i.e. plays the negative integers `-u`, `1 <= u <= 2^(k-1)`. The potential `u` grows by a factor of at most `(q+1)/2` per odd step and halves on even steps.
    - The bound is exact for this adversary when `q = 2^j - 1`: the free cycle `1 -> 2^(j-1) -> ... -> 1` has density `1/j`.
    - `log_(q+1) 2 > p0 = 0.2270922` exactly for `q <= 19`. So this raises THM-4481's universal floor for `q = 3, 5, ..., 19`. It **settles no new q**, since `log_(q+1) 2 < log_q 2`. It confines every class-(i) strategy of `q n ± 1` to `rho_max ∈ [log_(q+1) 2, log_q 2)`.
  - **Corollary 5 (q = 5 is settled at every level).** `rho*(5,k) = 1/2` for `k <= 6`, `3/7` for `7 <= k <= 10`, `5/12` for `11 <= k <= 14`, and **`2/5` for every `k >= 15`**.
    - The lower bound for all `k` comes from a corrected potential (`u^2`, with value 8 at `u = 3`) for the same adversary.
    - The obstruction is the sporadic 5x+1 cycle `(1, 3, 8, 4, 2)` of THM-4484, of density `2/5`. The adversary meets it as the 5x−1 cycle `(−1, −3, −8, −4, −2)`.
  - **Upper bound 1/2 (elementary, for completeness).** The level-2 max-halving strategy makes every odd step followed by an even one, so `rho*(q,k) <= 1/2` for all `q` and `k >= 2`.
  - **Proposition F (limits of stationary-law floors).**
    - For every odd `q`, the level-2 max-halving strategy has a unique uniform-lift stationary law, with `pi(odd) = 1/3` exactly (and `rho_max = 1/2`). Hence no floor derived from uniform-lift stationary laws (THM-4481's entropy law, the gain identity at `theta = 1/2`, block entropies of that chain) can exceed `1/3`, and none can settle `q <= 7`.
    - Certified strategies at levels 10 and 12 push this cap below `log_q 2` also for `q = 9` and `q = 11`.
  - **Lemma T (exact traps).** A finite set of rationals closed under `x -> x/2` and `x -> (qx ± 1)/2` gives a level-independent floor. It is weak: `1/2, 1/3, 1/4, 1/4` for `q = 3, 5, 7, 9`, and nothing found for `q = 11, 13`.
- **FINITE-EXACT** (both certificates re-checked with exact integers for every entry; runner §B, §C):
  - `rho*(q,k)` for `q = 3` (`k <= 16`), `q = 5` (`k <= 21`), `q = 7` (`k <= 22`), `q = 9, 11, ..., 21` (`k <= 18`).
  - All residue classes `±q mod 2^(k-1)` for `k <= 11`. Hence **for every odd q**:
    - `rho*(q,k) = 1/2` for `k <= 6`;
    - `rho*(q,k) >= 3/7` for `k <= 9`;
    - `rho*(q,k) >= 2/5` for `k <= 11`.
  - In particular **no q n ± 1 with odd q >= 7 has a provable sign strategy at any level k <= 11**. For `q = 7` the same holds for all `k <= 22`, and for `q = 9..21` for all `k <= 18`.
- **VERIFIED** (runner §A):
  - the drift lane's values are reproduced;
  - exhaustive enumeration agrees for `k <= 5` (all `q`);
  - an independent pure-Python Karp and the drift lane's own engine confirm `rho_max` of the extracted strategies;
  - corrupted certificates are rejected.
- **REFUTED**:
  - "`rho*(q,k) >= 3/7` for every odd `q >= 5` and every `k`" (Q1): `rho*(7,10) = 2/5`, `rho*(9,10) = rho*(11,10) = 5/12`, `rho*(5,11) = 5/12`.
  - The q-independence of `rho*` for `q = 5..11`, which holds only for `k <= 9`. It breaks at `k = 10`: `3/7, 2/5, 5/12, 5/12` for `q = 5, 7, 9, 11`.
- **EMPIRICAL**:
  - `rho*(q,k)` keeps decreasing slowly for every `q` in `7..21`: `rho*(7,22) = 14/37 = 0.3784`, still `0.022` above `log_7 2`.
  - The best stationary odd frequency `m(q,k)` is about `0.28–0.30` at `k = 16`.
  - The `tau = 1` values are stable for `14 <= k <= 16` when `q <= 19`, and are attained there by exact integer cycles.
- **OPEN**:
  - `lim_k rho*(q,k)` for `q >= 7`. For `q = 3, 5` the limit equals the value of the negative-integer adversary (the sparsest sign-choice cycle on the positive integers). If that held for `q = 7`, the limit would be `1/3 < log_7 2`, and 7n±1 would have a provable strategy at some finite level. The data to `k = 22` neither confirm nor exclude this (§8).
    - This naive guess is **REFUTED for `q = 31`**: the free cycle `(1, 16, 8, 4, 2)` gives the adversary only `1/5`, while `rho*(31,k) >= p0 = 0.227` by THM-4481.
  - A universal floor above `p0` valid for all `q` (§7).
- No HYP or THM file was created. Collatz itself is untouched (`q = 3`: `rho* = 1/2`, reproved by Theorem N).

Session `collatz-procgen-20260922`, floor lane, 2026-09-26.
- Scripts: `04-computation/experiments/procgen_floor_20260926_{lib,run}.py` and the C engine `procgen_floor_20260926_game.c`.
- Output: [procgen_floor_20260926.out](procgen_floor_20260926.out).
- Parents: [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) (Theorems A, C, E), [THM-4481](../../01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md) and its note [procgen_drift_20260926_positive_drift_provability.md](procgen_drift_20260926_positive_drift_provability.md) (§§1–3, 7), and [THM-4484](../../01-canon/theorems/THM-4484-free-and-sporadic-cycles.md) (free and sporadic cycles).

## 0. The answer in brief

| question | answer | status |
|---|---|---|
| Q1: is `rho*(q,k) >= 3/7` for every odd `q >= 5` and every `k`? | **No.** `rho*(7,10) = 2/5`; `rho*(5,k) = 2/5` for all `k >= 15`; `rho*(7,22) = 14/37`. | FINITE-EXACT, PROVED |
| Q1: why is `rho*` the same for `q = 5, 7, 9, 11` at `k <= 9`? | The level-`k` arena depends only on `±q mod 2^(k-1)` (Prop. S). At `k = 7` the floor `3/7` is attained exactly by the classes `±5, ±7, ±9, ±11 mod 64`, and all other classes have `1/2`. It is a coincidence of minima, not an isomorphism, and it breaks at `k = 10`. | PROVED + FINITE-EXACT |
| Q2: raise the floor `p0 = 0.2271` | **Theorem N**: `rho_max >= log_(q+1) 2` for every strategy of every `q`. This beats `p0` for `q <= 19` but settles no new `q`. | PROVED |
| Q2: can the entropy route reach `log_7 2`? | **No.** Any floor through uniform-lift stationary laws is `<= 1/3` (level-2 greedy); for `q = 9, 11` it is also below `log_q 2`. Settling `q = 7..11` needs an adversary that exploits cycles, not averages. | PROVED |
| Q3: what is the true floor? | `q = 3`: `1/2`; **`q = 5`: `2/5`** (from `k = 15` on); `q >= 7`: unknown. `rho*(7,k)` decreases to `0.378` at `k = 22` and is bounded below by `1/3`. | PROVED / OPEN |

The whole picture in one table (runner §B, §C, §D; `tau = 1` is the negative-integer adversary of Theorem N):

| `q` | `log_q 2` | `log_(q+1) 2` (Thm N) | `tau = 1` value (`k = 16`) | `rho*(q,k)`, last computed | first `k` with `rho* < 3/7` |
|---|---|---|---|---|---|
| 3 | 0.6309 | 1/2 | 1/2 | 1/2 (every `k`) | never |
| 5 | 0.4307 | 0.3869 | **2/5** | **2/5 (every `k >= 15`)** | 11 |
| 7 | 0.3562 | 1/3 | 1/3 | 14/37 = 0.3784 (`k = 22`) | 10 |
| 9 | 0.3155 | 0.3010 | 5/16 | 5/13 = 0.3846 (`k = 18`) | 10 |
| 11 | 0.2891 | 0.2789 | 2/7 | 41/110 = 0.3727 (`k = 18`) | 10 |
| 13 | 0.2702 | 0.2627 | 17/63 | 8/21 = 0.3810 (`k = 18`) | 11 |
| 15 | 0.2560 | 1/4 | 1/4 | 39/100 = 0.3900 (`k = 18`) | 13 |
| 17 | 0.2447 | 0.2398 | 1/4 | 24/61 = 0.3934 (`k = 18`) | 13 |
| 19 | 0.2354 | 0.2314 | 4/17 | 43/112 = 0.3839 (`k = 18`) | 11 |
| 21 | 0.2277 | 0.2242 | 9/35 | 25/66 = 0.3788 (`k = 18`) | 12 |

## 1. Setting

Fix an odd `q` and a level `k >= 2`; put `N = 2^k` and `H = 2^(k-1)`.
- A **sign strategy** `sigma` assigns `±1` to the odd residues mod `N`. It gives `T(n) = n/2` for even `n` and `(qn + sigma(n mod N))/2` for odd `n`. `R` is the set of residues with `sigma = -1`.
- The **parity graph** `G_sigma` has nodes `Z/N` and edges `s ->` the two lifts `P`, `P + H` of the target pair `P = T(s) mod H`.
- `rho_max(sigma)` is the largest odd density of a cycle of `G_sigma`, and `rho*(q,k) = min_sigma rho_max(sigma)`.
- Class (i) (bounded-lookahead provability, THM-4474 A) at level `k` is nonempty iff `rho*(q,k) < log_q 2`, iff `q^a < 2^p` for `rho* = a/p`. The check is exact.
- `rho*(q,k)` is non-increasing in `k`, because lifting a strategy keeps the map (drift lane §7).
- Threshold weights for `F = fn/fd`: `e(s) = fd - fn` for odd `s` and `-fn` for even `s`. A cycle has `e`-weight `> 0`, `= 0`, `< 0` iff its density is `> F`, `= F`, `< F`.
- `nu` is negation `s -> -s`, an automorphism of `G_sigma ≅ G_(nu sigma)` (THM-4474 E).

## 2. Theorem G: rho* is the value of a mean-payoff game

**The arena.**
- Nodes are `Z/N`. The pairs `P ∈ Z/H` have lifts `P` and `P + H`.
- **Options.** An even node `s` has one option, the pair `s/2`. An odd node `s` has two: `(qs+1)/2 mod H` (sign `+`) and `(qs-1)/2 mod H` (sign `-`).
- **Min** chooses an option at every odd node once and for all. That is exactly a sign strategy `sigma`, a positional strategy.
- **Max** chooses a lift at every pair once and for all: a lift strategy `tau : Z/H -> {0,1}`, also positional.
- The play from a node is eventually periodic, and its cycle is a cycle of `G_sigma`. Conversely, a densest cycle of `G_sigma` can be chosen to visit each pair at most once (a closed walk in the node–pair graph is a mediant of simple ones), and it is then such a play for a suitable `tau` and start.

**Lemma G1 (upper certificate).** Suppose `sigma` and integers `psi` satisfy `psi(t) + e(s) <= psi(s)` on every edge `s -> t` of `G_sigma`. Then `rho_max(sigma) <= F`.
*Proof.* Summing around a cycle gives `e`-weight `<= 0`. ∎

**Lemma G2 (lower certificate).** Let `W ⊆ Z/N` be nonempty, `tau` a lift strategy and `f : W -> Z`. Suppose that for every `s ∈ W` and every option `P` of `s` (both options if `s` is odd), `t = P + tau(P)H` lies in `W` and `f(t) <= f(s) + e(s)`. Then **every** sign strategy has a cycle of density `>= F`, i.e. `rho*(q,k) >= F`.

*Proof.* Fix any `sigma` and any `s_0 ∈ W`.
- Put `s_(i+1) = P_i + tau(P_i) H`, where `P_i` is the `sigma`-option of `s_i`.
- By hypothesis every `s_i` lies in `W` and `f(s_(i+1)) - f(s_i) <= e(s_i)`. Also `s_i -> s_(i+1)` is an edge of `G_sigma`.
- The sequence is eventually periodic. Summing over its cycle gives `0 <= Σ e(s_i) = fd·a - fn·p`, so the cycle has density `a/p >= F`. ∎

**Theorem G.** `rho*(q,k) = max_tau max_(s_0) min_sigma (density of the (sigma,tau)-cycle from s_0)`. At `F = rho*(q,k)` both certificates exist.

*Proof.*
- `>=` and `<=` for any certificate pair are Lemmas G1 and G2.
- Existence of matching certificates is positional determinacy of finite mean-payoff games with uniform optimal strategies (Ehrenfeucht–Mycielski 1979; used in this form by energy-game algorithms, Brim–Chaloupka–Doyen–Gentilini–Raskin 2011): CITED.
- The upper certificate at the value is elementary: any optimal `sigma` has a potential. The least-fixed-point Lemma below shows it can be chosen canonical. ∎

**Lemma L (least fixed point).** Let `(Φf)(s) = max(0, e(s) + min_(P option of s) max(f(P), f(P+H)))` on nonnegative integer functions.
- If some `sigma` has `rho_max(sigma) <= F`, then Kleene iteration from `f = 0` stabilizes at a finite `f*`.
- For every `argmin` choice `sigma*`, the pair `(sigma*, f*)` is an upper certificate.
- `f*` is invariant under every arena automorphism, in particular `nu`.

*Proof.*
- A nonnegative potential `psi` of `sigma` (longest-path potential, shifted) satisfies `Φ psi <= psi`. `Φ` is monotone, so the iterates stay `<= psi`, increase, and stabilize.
- At the fixed point `f*(s) >= e(s) + max over the lifts of the chosen pair`, which is the potential condition.
- `Φ` commutes with arena automorphisms (they preserve parity, options and lifts), so every iterate from 0 is invariant. ∎

**The algorithm** (`procgen_floor_20260926_game.c`, driven by `..._lib.py`).
- It runs a Stern–Brocot search on `F`. At each mediant it computes two least fixed points by FIFO chaotic iteration with a cap:
  - Min's operator `Φ` (upper certificate if no node exceeds the cap);
  - Max's operator `f(s) = max(0, -e(s) + max_P min(f(P), f(P+H)))` (lower certificate on `W = {f < cap}`, with `tau = argmin`).
- Caps only shrink the certified sets, so a too-small cap can only make the search fail, never certify a wrong value. The failure is detected and the search restarts with 8× caps; this never happened in the runs.
- Every certificate is re-checked in Python with exact integers.
- Cost: `k = 21` takes about 10 s and `k = 22` about 230 s in the runner. The SAT descent of the drift lane stopped at `k = 8–9`.
- In every computed case `W` is all of `Z/N`: the game value is the same from every node (EMPIRICAL).

## 3. Proposition S: what the value depends on

**(a) Only `q mod 2^(k-1)` matters.** Write `m_q(s) = (qs+1)/2 mod H`. Then `m_(q+H)(s) = m_q(s) + H/2 = m_q(s+H)` for odd `s`, and likewise for the `-` option.
- So in `arena(q+H)`, node `s` has exactly the options of node `s+H` in `arena(q)`. Even nodes are unchanged.
- Swapping the two lifts of every odd pair is therefore an isomorphism. It maps `G_sigma(q)` to `G_(sigma')(q+H)` with `sigma'(s) = sigma(s+H)`.
- Hence `rho*(q,k) = rho*(q + 2^(k-1), k)`, and the whole level-`k` classification except the threshold `log_q 2` depends only on `q mod 2^(k-1)`. The runner checks the option identity for `k <= 12` and all `q` (§C1).

**(b) `rho*(-q,k) = rho*(q,k)`.**
- *Step 1.* Let `F = rho*(q,k)`, and let `(sigma*, f*)` be the Lemma L certificate for `q`; `f*` is `nu`-invariant.
- *Step 2.* In `arena(-q)`, the target of odd `s` under the sign `-sigma*(s)` is `(-qs - sigma*(s))/2 = -((qs + sigma*(s))/2)`, minus its target in `arena(q)`.
- *Step 3.* So every edge `s -> t'` of `G_(-sigma*)(-q)` has `t' = -t` for an edge `s -> t` of `G_(sigma*)(q)`. Even edges are identical.
- *Step 4.* Then `f*(t') + e(s) = f*(t) + e(s) <= f*(s)`, so `(-sigma*, f*)` certifies `F` for `-q`, and `rho*(-q,k) <= rho*(q,k)`.
- *Step 5.* Exchange `q` and `-q`. ∎
- *Checks.* Runner §C3 checks the transported certificates for ten `(q,k)`; the `nu`-symmetric Max certificates transport too. §C2 checks both symmetries on all residues for `k <= 10`.
- *Not an isomorphism.* The full `rho_max` distributions differ, e.g. at `k = 4` for `q = 3` versus `5 = -3 mod 8`. Only the `nu`-invariant strategies (`sigma(-s) = -sigma(s)`) correspond, via `sigma -> -sigma` (§C4).

**(c) No other structural coincidence.** Let `k >= 3` and suppose `phi` is an arena isomorphism `arena(q) -> arena(q')`. It is a bijection of nodes preserving parity, pairs and option sets; let `pi` be the induced map on pairs.
- *Step 1.* The odd option sets are exactly the edges `{m-1, m}` of the cycle `C_H`, each once, because `s -> m_q(s)` is a bijection from odd residues onto `Z/H`. So `pi` is a dihedral map `x -> εx + c`.
- *Step 2.* The even node `2x` has the option `{x}` and lies in pair `2x mod H`. This forces `2 pi(x) ≡ pi(2x) (mod H)`, hence `c ≡ 0` and `pi = ±id`.
- *Step 3, `pi = id`.* The images of the two nodes of an odd pair `P` are the two nodes of `P`. An edge `{m-1, m}` of `C_H` determines `m`, and `m_(q')(P+H) = m_(q')(P) + H/2`.
  - So `m_(q')(P) ≡ m_q(P) (mod H/2)` for every odd `P`, i.e. `(q'P + 1)/2 ≡ (qP + 1)/2 (mod H/2)`.
  - Since `P` is odd, this is `q' ≡ q (mod H)`.
- *Step 3, `pi = -id`.* The option `{m-1, m}` goes to `{-m, 1-m}`, so `m_(q')(φ(s)) = 1 - m_q(s)` with `φ(s) ≡ -s (mod H)`.
  - This gives `(-q's + 1)/2 ≡ (1 - qs)/2 (mod H/2)`, again `q' ≡ q (mod H)`.
  - The map is then the negation automorphism `nu` composed with an isomorphism of type (a). ∎

**The q-independence observed by the drift lane.**
- The classes `±q mod 2^(k-1)` number `2^(k-3)`. Runner §C5 computes `rho*` on all of them:

| `k` | classes | level floor `min_q rho*(q,k)` | classes attaining it |
|---|---|---|---|
| 3–6 | 1, 2, 4, 8 | **1/2 for every class** | all |
| 7 | 16 | 3/7 | exactly `±5, ±7, ±9, ±11` mod 64 (all others 1/2) |
| 8 | 32 | 3/7 | `±5, ±7, ±9, ±11, ±53, ±55, ±57, ±59` mod 128 |
| 9 | 64 | 3/7 | 20 classes |
| 10 | 128 | 2/5 | `±7, ±59, ±117, ±245` mod 512 |
| 11 | 256 | 2/5 | 12 classes: `±7, ±59, ±69, ±117, ±245, ±267, ±373, ±395, ±443, ±453, ...` mod 1024 |

- So `q = 5, 7, 9, 11` agree at `k <= 9` because at `k = 7` they are exactly the four classes where the floor first drops. They stay at the floor through `k = 9`, then split at `k = 10`.
- The floors themselves are **universal statements**, since every odd integer `q` lies in one class:
  - every sign strategy of every `q n ± 1` has `rho_max >= 1/2` at `k <= 6`, `>= 3/7` at `k <= 9` and `>= 2/5` at `k <= 11`;
  - hence **class (i) is empty at every level `k <= 11` for every odd `q >= 7`** (`log_q 2 <= 0.3562 < 2/5`).

## 4. Exact values (runner §B; every entry carries both certificates)

`rho*(q,k)`, listed at the levels where it changes (the value holds up to the next listed level):

| `q` | values | last `k` |
|---|---|---|
| 3 | 1/2 | 16 (and every `k`, Theorem N plus greedy) |
| 5 | 1/2 (2), 3/7 (7), 5/12 (11), **2/5 (15)** | 21 (and every `k`, Corollary 5) |
| 7 | 1/2 (2), 3/7 (7), 2/5 (10), 15/38 (14), 7/18 (16), 19/49 (18), 13/34 (19), 34/89 (20), 8/21 (21), **14/37 (22)** | 22 |
| 9 | 1/2, 3/7 (7), 5/12 (10), 2/5 (14), 12/31 (17), 5/13 (18) | 18 |
| 11 | 1/2, 3/7 (7), 5/12 (10), 9/22 (12), 2/5 (13), 13/33 (14), 45/116 (15), 21/55 (16), 29/77 (17), 41/110 (18) | 18 |
| 13 | 1/2, 4/9 (9), 11/26 (11), 5/12 (12), 20/49 (13), 2/5 (14), 35/88 (15), 16/41 (16), 5/13 (17), 8/21 (18) | 18 |
| 15 | 1/2, 5/11 (8), 9/20 (10), 7/16 (11), 3/7 (12), 21/50 (13), 7/17 (14), 2/5 (15), 32/81 (17), 39/100 (18) | 18 |
| 17 | 1/2, 5/11 (8), 4/9 (9), 7/16 (11), 10/23 (12), 29/69 (13), 7/17 (14), 15/37 (15), 27/67 (16), 53/133 (17), 24/61 (18) | 18 |
| 19 | 1/2, 4/9 (9), 3/7 (10), 8/19 (11), 7/17 (12), 15/37 (14), 2/5 (15), 13/33 (16), 23/59 (17), 43/112 (18) | 18 |
| 21 | 1/2, 4/9 (9), 10/23 (10), 3/7 (11), 7/17 (12), 11/27 (13), 2/5 (14), 24/61 (15), 33/85 (16), 38/99 (17), 25/66 (18) | 18 |

- For every `q` in `7..21` and every computed `k`, `q^a > 2^p` for `rho* = a/p`: **class (i) is empty there**. So 7n±1 has no provable sign strategy at `k <= 22`, and 9n±1, ..., 21n±1 have none at `k <= 18`.
- For `q = 5`, class (i) is nonempty exactly from `k = 7`.
- The certificate hashes of every entry with `k >= 14` are printed in the .out.
- The critical cycles move with `k` (exploratory, not in the runner). For `q = 7` at `k = 17`, the tight part of the Max certificate has three small strongly connected components, of sizes 18, 18 and 28. The two 18-cycles, with 7 odd steps each, are negatives of each other.
- For `q = 5` at `k = 16` the tight part contains the cycle `{1, 3, 8, 4, 2}`.

## 5. Theorem N: the negative-integer adversary

**Theorem N.** For every odd `q >= 3`, every `k >= 2` and every sign strategy `sigma`, `G_sigma` contains a cycle of odd density `>= log_(q+1) 2 = 1/log2(q+1)`.

*Proof.*
- **The adversary.** Max always takes the lift with top bit 1 (`tau ≡ 1`). The set `W = {x : H <= x < N}` is closed under every option followed by `tau`.
- **Coordinates.** Write `x = N - u` with `1 <= u <= H`; these are the negative integers `-u`. For even `x` the next node is `N - u/2`. For odd `x` with sign `sigma(x)`, the next node is `N - rho(m)`, where `m = (qu - sigma(x))/2 >= (q-1)/2 >= 1` and `rho(m) ∈ [1, H]` is the representative of `m mod H`.
- So the Min graph of this adversary is the graph `U_k`: `u -> u/2`, and `u -> rho((qu ± 1)/2)` for odd `u`. Runner §D1 checks this identification for `k <= 14`, `q <= 65`.
- **The potential.** Since `rho(m) <= m`, every odd step has `v <= (qu+1)/2 <= (q+1)u/2`, and every even step has `v = u/2`.
- **A cycle.** For every `sigma`, the play from any `u` closes a cycle of `G_sigma`. Along it `Π v/u = 1`, so `1 <= ((q+1)/2)^a (1/2)^(p-a)`, i.e. `2^p <= (q+1)^a`. ∎

**Remarks.**
1. **Exactness for the adversary.**
   - The inequality is attained at `u = 1` by `(q·1 + 1)/2`.
   - For `q = 2^j - 1`, the free cycle `1 -> 2^(j-1) -> ... -> 2 -> 1` (THM-4484, shape `(j,1)`) has density exactly `1/j = log_(q+1) 2`. So the value of `tau ≡ 1` is exactly `log_(q+1) 2` (checked for `q = 3, 7, 15, 31, 63`).
   - For `q = 3` this is `rho*(3,k) = 1/2`: the cycle `{1, 2}` that pins `q = 3`.
2. **Comparison with the entropy floor.**
   - `log_(q+1) 2 > p0 = 0.2270922` exactly for `q <= 19` (runner §D). So for `q = 5..19` Theorem N is the better universal-in-`k` floor. For example `q = 7` rises from `0.227` to `1/3`, and `q = 9` to `0.301`.
   - For `q >= 21` the entropy floor is larger.
   - `log_(q+1) 2 < log_q 2` always, so Theorem N settles no new `q`. It narrows the class-(i) window to `rho_max ∈ [log_(q+1) 2, log_q 2)`; for `q = 19` that is `[0.2314, 0.2354)`.
3. **What the adversary sees: the sparsest sign-choice cycles of `q u ∓ 1` on the positive integers** (runner §D5, EMPIRICAL at the listed levels). The value of `tau ≡ 1` is stable in `k` and is attained by an exact integer cycle without wrap-around:

| `q` | `tau = 1` value | critical cycle (positive integers) | type |
|---|---|---|---|
| 3 | 1/2 | `(1, 2)` | free (`4 - 3 = 1`) |
| 5 | 2/5 | `(1, 3, 8, 4, 2)` | the sporadic 5x+1 cycle of THM-4484 (gap `2^5 - 5^2 = 7`) |
| 7 | 1/3 | `(1, 4, 2)` | free (`8 - 7 = 1`) |
| 9 | 5/16 | `(1, 5, 22, 11, 50, 25, 112, 56, 28, 14, 7, 32, 16, 8, 4, 2)` | mixed signs, contracting |
| 11 | 2/7 | `(1, 6, 3, 16, 8, 4, 2)` | mixed signs, contracting |
| 13 | 17/63 | a 63-cycle through `3, ..., 1840` | contracting (`13^17 < 2^63`) |
| 15, 17 | 1/4 | `(1, 8, 4, 2)` | free (`16 ∓ 1`); expanding for 17 |
| 19 | 4/17 | `(9, 86, 43, 408, ..., 485, 4608, ..., 18)` | contracting (`19^4 = 130321 < 2^17`) |
| 21 | 9/35 (`k = 16`) | a 70-cycle with 3 wrap-arounds | not stable in `k` |

   - In `u`-coordinates the adversary sees the orbits `u -> (qu ∓ 1)/2` of the negatives `-u`. So the 5x+1 cycle `(1, 3, 8, 4, 2)` is met as the 5x−1 cycle `(−1, −3, −8, −4, −2)`, and `(1, 4, 2)` as `(−1, −4, −2)` of 7x−1.
   - **Heuristic.** Sign-choice cycles with `a` odd steps and density near `c = log_q 2` should exist at large heights in number growing like `2^(a (h(c) + c - 1)/c)`: `C(p,a) 2^a` words and signs, each integral with probability about `1/|2^p - q^a|`. The exponent is positive iff `c > p0`, i.e. iff `q <= 21`. This is the entropy law's counting again.
   - So for `q <= 21` the adversary's value should tend to at most about `log_q 2` as `k` grows, and it cannot settle any `q`. For `q = 17` its value `1/4 > log_17 2` at `k <= 16` is then a small-height effect. EMPIRICAL/heuristic.

## 6. Corollary 5: q = 5 at every level

**Claim.** Every cycle of `U_k` (`q = 5`, any `k >= 2`) has density `>= 2/5`.

*Proof.* Take `Phi(u) = u^2` for `u != 3` and `Phi(3) = 8`. We show `Phi(v) <= 8 Phi(u)` on odd steps and `4 Phi(v) <= Phi(u)` on even steps.
- **Even `u`, `v = u/2`.** If `v != 3`, then `4 Phi(v) = u^2 = Phi(u)`. If `v = 3`, then `u = 6` and `4·8 = 32 <= 36`.
- **Odd `u`.** Here `v = rho(m)` with `2 <= m = (5u ± 1)/2`, and `v <= m`.
  - `u = 1`: `v ∈ {1, 2, 3}`, and `Phi(v) <= 8 = 8 Phi(1)`.
  - `u = 3`: `v <= 8`, and `Phi(v) <= 64 = 8 Phi(3)`.
  - `u >= 5`: `Phi(v) <= max(8, ((5u+1)/2)^2) <= 8u^2`, since `7u^2 - 10u - 1 >= 0`.
- **A cycle.** Along a cycle `Π Phi(v)/Phi(u) = 1`, so `1 <= 8^a 4^(-(p-a))`, i.e. `5a >= 2p`. ∎

Equivalently, `f = log2 Phi` is a real-valued lower certificate at `F = 2/5` (`e = +3` odd, `-2` even); Lemma G2 holds verbatim for real `f`.

The runner checks the edge inequalities for `k <= 22` (§E1). The cycle `(1, 3, 8, 4, 2)` attains `2/5`.

**Corollary 5.**
- Lemma G2 with this adversary gives `rho*(5,k) >= 2/5` for every `k`.
- Lifting the certified level-15 optimum gives `rho*(5,k) <= 2/5` for `k >= 15`.
- With the certified values for `k <= 14`, the whole sequence is `1/2, 3/7, 5/12, 2/5` on `k <= 6`, `7..10`, `11..14`, `>= 15`.

It is the first `q >= 5` whose min-max density is known at every level. The limit `2/5` lies inside `(log_6 2, log_5 2) = (0.387, 0.431)`.

## 7. Q2: raising the universal floor — what works and what cannot

**What is proved.** Theorem N (§5) is an "unavoidable dense cycle" lemma of exactly the kind asked for: a single adversary that forces, in every strategy, a cycle through a family of integer orbits. It raises the floor from `p0` to `log_(q+1) 2` for `q <= 19`. It settles no new `q`.

**Proposition F: entropy-type arguments are capped at 1/3.** A *stationary-law floor* is a number `x` with `pi(odd) >= x` for every stationary law of the uniform-lift chain of every strategy at every level. THM-4481's `p0` is one.
- **The cap.** The level-2 max-halving strategy (flip `r` iff `4 ∤ qr + 1`) has chain `0, 1, 3 -> {0, 2}` and `2 -> {1, 3}`, independent of `q`.
  - Its unique stationary law is `(1/3, 1/6, 1/3, 1/6)`, so `pi(odd) = 1/3` exactly, while `rho_max = 1/2`.
  - Hence every stationary-law floor is `<= 1/3`: runner §F1, exact fractions, `q = 3..23`.
- **Consequence.** Block entropies, merge structure and the gain identity of the uniform chain cannot settle `q <= 7`.
- **q = 9 and 11.** Certified strategies close the route too.
  - `q = 9`, `k = 10`: every closed class has `pi(odd) <= 0.31 < log_9 2`.
  - `q = 11`, `k = 12`: every closed class has `pi(odd) <= 0.2885 < log_11 2`.
  - *Certificate:* integer functions `h` with `2^24[s odd] + (h(t1) + h(t2))/2 <= 2^24 c + h(s)` on each closed class, checked exactly (runner §F2).
- **The empirical cap.** `m(q,k)`, the least stationary odd frequency (float value iteration, runner §F3):
  - `0.303, 0.296, 0.291, 0.277` at `k = 16` for `q = 5, 7, 9, 11`;
  - `0.279–0.287` for `q = 13..23`.
  - So for `q = 13..21` the entropy route is not yet excluded: its cap is still above `log_q 2`, and decreasing slowly.
- **The gap to the true values.** The min-max values `rho* ≈ 0.37–0.39` are far above these stationary densities. The optimal strategies themselves have `pi(odd) ≈ 0.30` (exploratory: `q = 7`, `k = 12` gives `0.3015` against `rho_max = 2/5`). The floor is a worst-cycle phenomenon; averages cannot see it.

**What an improvement must use.** By Theorem G the best floor at level `k` is attained by a positional adversary `tau`. A proof valid at all levels needs a family of adversaries with a level-independent potential. Three candidate families were tried:
- **Exact traps** (Lemma T: `max_S min-cycle(S)` over finite closed rational sets `S`). *Proof.* For any `sigma`, put `ε(x) = sigma(x mod 2^k)` on `S`. The `T_ε`-orbits stay in `S`, and each of their cycles projects to a closed walk of `G_sigma` with the same parities. No separation of residues is needed: a `sigma` that is constant on residue classes is just a special `ε`. ∎
  - Best found (denominators `<= 31` in the runner, `<= 61` in an exploratory run with the same results): `q = 3`: 1/2 (`{1, 2}`); `q = 5`: 1/3 (`{1/3, 4/3, 2/3}`); `q = 7`: 1/4 (`D = 9`); `q = 9`: 1/4 (`{1/7, ...}`); `q = 11, 13`: none except `{0}`.
  - All are weaker than Theorem N, because the adversary may not use lift freedom.
- **The negative-integer adversary** (Theorem N, corrected potentials as in §6). At the computed levels (`k <= 16`) its value is the sparsest sign-choice cycle on the positive integers, which is below `log_q 2` for `q = 5..15` and for 19.
- **Rational windows** (Max plays `a/D` with `a` in a window of length `H`, i.e. the dynamics of `qx ± D` with wrap-around). The best value for `q = 7` is `1/3`, over `D ∈ {1, 3, 5, 9, 11, 13, 15}` and windows `A/H ∈ {-1, -1/2, 0}` (plus exploratory windows): the free cycle `{D, 4D, 2D}` of `7x ± D` is always available to Min (runner §G2).

The optimal adversaries at `k = 17–22` (value `≈ 0.38` for `q = 7`) mix both lifts irregularly near small integers. For example, the extracted optimal `tau` at `q = 7`, `k = 17` takes `H - 1` instead of `-1` and `1 - H` instead of `1` (exploratory).
- No level-independent description was found.
- A height-type potential (a function of the absolute value of the representative) cannot certify such an adversary: a jump from height 1 to height `≈ H` must be paid by the potential within one step.
- So an all-level proof beyond `log_(q+1) 2` would need a different kind of potential (heuristic).

## 8. Q3: the true floor, the evidence, and the obstruction

- **q = 3:** `1/2` at every level, pinned by the free cycle `{1, 2}` (Theorem N with `q + 1 = 4`).
- **q = 5:** `2/5` from `k = 15` on, pinned by the 5x+1 cycle `(1, 3, 8, 4, 2)`. It is reached after the plateaus `3/7` (`k = 7–10`) and `5/12` (`11–14`).
- **q = 7:** `rho*(7,k)` takes the values `3/7, 2/5, 15/38, 7/18, 19/49, 13/34, 34/89, 8/21, 14/37` (`k = 9..22`).
  - Proven lower bound: `1/3`, the free cycle `(1, 4, 2)` of 7x+1 (Theorem N). The adversary value is exactly `1/3`.
  - The gap to `log_7 2 = 0.3562` at `k = 22` is `0.0222`. The last four changes of value (`k = 19..22`) were `0.0054, 0.0003, 0.0011, 0.0026`.
  - Some of the values are Fibonacci ratios (`13/34, 34/89, 8/21`), but `8/21 < 1/phi^2 = 0.38197`, so the limit is not `1/phi^2` (REFUTED as a guess).
- **q = 9..21:** the values decrease slowly (`0.37–0.39` at `k = 18`), far above `log_q 2` (`0.23–0.32`) and above the `tau = 1` values.

**The dichotomy.** Either `lim_k rho*(7,k) > log_7 2`, and 7n±1 is never provable, or some finite level has a provable 7n±1 strategy.
- A proof of the first needs, at every level, an adversary of value `> log_7 2`. The only level-independent adversaries found (§7) are capped at `1/3` by the free cycle `(1, 4, 2)`.
- For `q = 3, 5` the limit equals the negative-integer value.

**Conjecture candidate (OPEN, weak evidence).** Is `lim_k rho*(q,k)` the least density of a sign-choice cycle of `u -> u/2`, `(qu ± 1)/2` on the positive integers?
- It holds for `q = 3, 5`.
- For `q = 7` it predicts `1/3 < log_7 2`, i.e. a provable 7n±1 strategy at a finite level.
- It is **false for `q = 31`**, and for every `q = 2^j - 1` with `1/j < p0`, i.e. `j >= 5`. The free cycle `1 -> 2^(j-1) -> ... -> 1` gives the adversary exactly `1/j`, while `rho*(q,k) >= p0 = 0.2271` (THM-4481). So at best it is a statement for small `q`.
- The data do not support extrapolation: at `k = 22` the `q = 7` value is still `0.045` above `1/3`, and `q = 5` needed 15 levels to reach its limit.

**The obstruction to proving any floor above `log_(q+1) 2`.**
- Every level-independent argument found so far is a height potential for an adversary that plays small integers.
- Min can always use the small contracting integer cycles, e.g. `(1, 4, 2)` for 7x+1, against such an adversary.
- The adversaries that do better at finite levels escape to heights near `2^k`, where height potentials fail (heuristic, from the extracted optimal `tau`). Their advantage decays slowly with `k`, and its limit is exactly the open question.

## 9. Failures and caught mistakes

- **A cap issue during development.** The first search used the full energy bound `H·fd` as the cap. It climbed `O(N^2)` when losing, e.g. 13 s at `k = 14`. Moderate caps with detected failure and retry fixed it; caps only shrink certified sets.
- **The Fibonacci guess.** From `13/34, 34/89` I guessed a golden-ratio limit `1/phi^2` for `q = 7`. `rho*(7,21) = 8/21 < 1/phi^2` refutes it. It was never recorded as a claim.
- **A premature symmetry proof.** I first tried to prove `rho*(q) = rho*(-q)` by an arena isomorphism. The runner shows the full `rho_max` distributions differ (§C4), and Prop. S(c) excludes an isomorphism. The correct proof goes through `nu`-invariant least fixed points (§3b).
- **Equivariance at one pair.** `tau` cannot be `nu`-equivariant at the pair `H/2`, whose two lifts `nu` swaps. There `f` ties automatically, so the transport still works (§C3). The equivariance check exempts that pair.
- **Memory.** The `k = 23` run for `q = 7` was stopped before completion because its projected RSS (about 800 MB) exceeded the lane limit. `k = 22` peaked near 400 MB standalone.
- **Nothing refuted from the canon.**
  - THM-4474 and THM-4481 are used as stated. THM-4481's table entries (`rho* = 3/7` at `k = 7, 8`; 7n+1 at `k = 9`) are reproduced.
  - Its sentence "for q = 7, 9, 11 no provable strategy exists at the computed levels" is now extended to `k <= 22` (q = 7) and `k <= 18`.
  - Its open question whether any `q ∈ 7..21` is ever provable remains OPEN (§8).

## 10. Reproduction

```bash
/usr/bin/time -l python3 -u 04-computation/experiments/procgen_floor_20260926_run.py > 05-knowledge/results/procgen_floor_20260926.out
```

- **Requirements.** `numpy`, `scipy`, `mpmath` and a C compiler.
  - The engine is compiled into `scratch/procgen_floor/game_<source-hash>.so`, which is not committed.
  - Section A3 imports the drift lane's `procgen_drift_20260926_lib.py` read-only; its engine build already exists in `scratch/procgen_drift/`.
- **Environment.** `FLOOR_KMAX5` (21), `FLOOR_KMAX7` (22), `FLOOR_KMAXQ` (18, for `q = 9..21`), `FLOOR_KALLQ` (11, all residue classes).
- **Checks.** Every claim in the output is a `check(...)` that raises on failure. The output ends with `ALL CHECKS PASSED`.
- **Cost.** One process; wall time `710 s` on the shared 8-core machine (load average about 5 from other lanes).
  - Peak RSS is `482 MB` self-reported, and `/usr/bin/time -l` reports a maximum resident set of `505,200,640` bytes. The peak comes from `q = 7`, `k = 22` (`231 s`).
  - The `k = 23` computation was not attempted in the runner: its projected RSS of about 800 MB exceeds the lane limit.
- **Output.** 41 checks, ending with `ALL CHECKS PASSED`.
- **SHA-256** (raw bytes; the runner prints the first three itself):
  - `procgen_floor_20260926_lib.py` `be7024fa7dc147c21d3eacabfc59d09c558943dcfa3cb30b421abda79e3e4f2f`
  - `procgen_floor_20260926_game.c` `dc2b9e7eec12971ff5a7c54e55cd40144867856d9a33053cfa49f8f088d8cf71`
  - `procgen_floor_20260926_run.py` `311862d4c19e37378e35c8a64f8f023ec29da875156338ae30576a40ae6aa8cf`
  - `procgen_floor_20260926.out` `34374e14665bb374c0c27c22e91d867515abb780c15a3c466084d607d8de12d5`
- **Timing lines** in the .out (seconds, RSS) vary between runs; everything else is deterministic.
