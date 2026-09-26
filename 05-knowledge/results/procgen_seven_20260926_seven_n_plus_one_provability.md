# Is 7n±1 ever bounded-lookahead provable? Exact min-max densities to level 29 (a 3/8 plateau at k = 24–26, then 35/94 and 13/35), a lean exact solver, the rational relaxation, and all-level floors 5/16 (q = 9) and 2/7 (q = 11)

**Status.**
- **Answer to the lane question: UNDECIDED.** Neither an all-level lower certificate above `log_7 2` (route A) nor a level with `rho*(7,k) < log_7 2` (route B) was found. What is new:
- **FINITE-EXACT** (runner §B). Every exact value is backed by a lower and an upper certificate, re-checked on all `2^k` nodes by an independent exact checker; at `k = 29` the upper bound comes from `k = 28` by monotonicity:
  - `rho*(7,k)` for `k = 19..29`: `13/34, 34/89, 8/21, 14/37, 14/37, 3/8, 3/8, 3/8, 35/94, 13/35, 13/35`. The values for `k <= 22` reproduce THM-4486; **`k = 23..29` are new**.
    - `k = 28`: a lower and an upper certificate at the same threshold.
    - `k = 29`: a lower certificate at `13/35`; the upper bound is monotonicity.
  - `rho*(7,30) >= 7/19 = 0.36842`.
  - Hence **7n±1 has no bounded-lookahead-provable sign strategy at any level `k <= 30`** (`7^13 > 2^35`, `7^7 > 2^19`). At `k = 29` the value `13/35` is still `0.0152` above `log_7 2 = 0.356207`.
  - `rho*(9,k)`, `k = 17..24`: `12/31, 5/13, 8/21, 17/45, 3/8, 19/51, 27/73, 18/49`; `rho*(11,k)`, `k = 17..24`: `29/77, 41/110, 7/19, 4/11, 71/197, 29/81, 11/31, 31/88`. **No provable 9n±1 or 11n±1 strategy at `k <= 24`.**
- **PROVED** (hand proofs below; the runner checks their finite content):
  - **Theorem Q (q = 9 and q = 11 at every level).** `rho*(9,k) >= 5/16` and `rho*(11,k) >= 2/7` for **every** `k`. The proof is a corrected potential `u^5` (resp. `u^2`) for THM-4486's negative-integer adversary, exact on `u < 16384` (resp. `u < 1024`). These are exactly that adversary's values: its critical cycles are `(1, 5, 22, 11, 50, 25, 112, 56, 28, 14, 7, 32, 16, 8, 4, 2)` for 9u∓1 and `(1, 6, 3, 16, 8, 4, 2)` for 11u∓1. This improves Theorem N (`log_10 2 = 0.30103`, `log_12 2 = 0.27894`) at every level. It settles nothing: `5/16 < log_9 2` and `2/7 < log_11 2`.
  - **Proposition S (the rational relaxation).** The real-sign rule `sigma(x) = sgn(x)` on rationals with odd denominator keeps **no** expanding cycle: every cycle of `x -> x/2, (qx + sgn x)/2` has `2^p = prod(q + 1/|x_i|) > q^a`. Hence **no finite set of rational cycles can obstruct provability at every level**: once their points have distinct residues mod `2^k`, some level-`k` strategy breaks them all. An all-level proof that 7n±1 is never provable must use residue collisions (rationals congruent mod `2^k`) at every level. Consistently, the critical cycles found in the certificates have huge denominators, such as `7^35 - 2^94` at `k = 27` (exploratory).
  - **Proposition P (pair-form certificates).** The compact certificate format checked by the independent checker (a lift bit and a potential per pair `P mod 2^(k-1)`; a sign per odd node and a potential per pair) proves `rho* >= F` and `rho* <= F`. It is Lemma G1/G2 of THM-4486 restricted to potentials that depend only on the pair.
- **VERIFIED / EMPIRICAL** (runner §A, §C, §E):
  - The lean solver reproduces 20 certified values of THM-4486 and 6 fresh values of the floor lane's engine.
  - Karp's algorithm re-evaluates 5 certified strategies at `k = 8..10` and agrees.
  - The floor lane's Python verifiers accept the same certificates as the checker, and the checker rejects corrupted ones.
  - **The 3/8 plateau** (`k = 24, 25, 26`). At `k = 24` (runner §C) the certified optimal sign strategy keeps the rational 7x+1 cycle `(-169, -548, -274, -137, -436, -218, -109, -338)/87`: shape `(8,3)`, and `2^8 - 7^3 = -87`, so it is an integer cycle of 7x+87 of the free shape `(8,3)` (THM-4484). All its edges are tight for the certified potential. It keeps 4 of the 54 simple `(8,3)`-cycles. Min's best responses to the certified Max strategy are density-3/8 cycles of shapes `(8,3), (24,9), (32,12), ..., (64,24)`.
  - At `k = 27` all 54 `(8,3)`-cycles are broken (runner §C). Exploratory tight walks find only critical cycles of shapes `(94,35)` and `(188,70)`: real-bounded rationals (`|x| < 120` along the orbits) with denominators dividing `7^35 - 2^94 ~ 3.6 * 10^29` (resp. `7^70 - 2^188`).
  - Every level-independent adversary tried for q = 7 has value at most `1/3` at `k = 8, 10, 12`: the top lift, five shifted windows, random lifts, and the top lift switched at the free cycle's classes. The last one, the suggested fix, is strictly worse at `k = 10, 12` (`7/24, 5/19` and `7/24, 11/38`).
- **REFUTED:** "`lim rho*(7,k) = 3/8`". It is the natural guess from the plateau, as `2/5` was final for q = 5, but `rho*(7,27) = 35/94 < 3/8`.
- **OPEN:** whether `lim_k rho*(7,k) < log_7 2`. The proved floor is still `1/3` (Theorem N). The certified values decrease irregularly, by about `0.0012` per level over `k = 20..29`. A linear extrapolation would cross `log_7 2` near `k ~ 42`, but the decrements appear to shrink; the data do not decide it (§8).
- No HYP or THM file was created. Collatz is untouched.

Session `collatz-procgen-20260922`, lane **seven**, 2026-09-26.
- Scripts: `04-computation/experiments/procgen_seven_20260926_{game,verify,tight}.c`, `procgen_seven_20260926_{run,potential,cycles}.py`. The floor lane's `procgen_floor_20260926_lib.py` is reused read-only (its engine, its certificate verifiers, `value_of_tau`).
- Output: [procgen_seven_20260926.out](procgen_seven_20260926.out).
- Parents:
  - [THM-4486](../../01-canon/theorems/THM-4486-min-max-cycle-density-game.md) and its note [procgen_floor_20260926_density_floor.md](procgen_floor_20260926_density_floor.md): the game, Theorem N, Corollary 5;
  - [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) (class (i) iff `rho_max < log_q 2`);
  - [THM-4481](../../01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md), [THM-4484](../../01-canon/theorems/THM-4484-free-and-sporadic-cycles.md) (free shapes).

## 0. The answer in brief

| question | answer | status |
|---|---|---|
| (A) a lower certificate family above `log_7 2` at every level? | Not found. Every level-independent adversary tried is capped at `1/3`. Proposition S shows why an all-level certificate cannot be a finite list of rational cycles: it must use residue collisions at every level. | OPEN (+ PROVED obstruction to one approach) |
| (B) does `rho*(7,k)` go below `log_7 2`? | Not up to `k = 30`: `rho*(7,27) = 35/94`, `rho*(7,28) = rho*(7,29) = 13/35`, `rho*(7,30) >= 7/19`. The sequence plateaus at exactly `3/8` for `k = 24..26`, then drops. | FINITE-EXACT; limit OPEN |
| the Fibonacci-like denominators `34, 89` and the break at `k = 22` | Not structural. The later values `14/37, 3/8, 35/94, 13/35` are not Fibonacci ratios. The values for `k <= 23` lie in `(3/8, 2/5)` and approach `3/8`, which is attained for `k = 24..26`; then the sequence moves below it (§3). | EMPIRICAL |
| (C) q = 9, 11 | Exact values to `k = 24` (still `> 0.35`, far above `log_q 2`). **All-level floors `5/16` and `2/7`** (Theorem Q), just below `log_9 2 = 0.3155` and `log_11 2 = 0.2891`. | FINITE-EXACT; PROVED |

## 1. Setting; the pair form of the game (Proposition P)

Notation as in THM-4486.
- Level `k`, odd `q`, `N = 2^k`, `H = 2^(k-1)`. A sign strategy `sigma` fixes a sign at each odd residue mod `N`.
- `G_sigma` has nodes `Z/N`, and an edge from `s` to both lifts of the target pair `T_sigma(s) mod H`.
- `rho*(q,k) = min_sigma rho_max(sigma)`. Class (i) (bounded-lookahead provability, THM-4474 A) is nonempty at level `k` iff `rho*(q,k) < log_q 2`, iff `q^a < 2^p` for `rho* = a/p`.
- For a threshold `F = fn/fd` put `e(x) = fd - fn` (odd `x`), `-fn` (even `x`).

**The pair game.**
- A pair `P in Z/H` has the lifts `x_b = P + bH` (`b = 0, 1`), both of the parity of `P`.
- The options of `x_b` are the pairs `x_b/2 mod H` (even), or `(q x_b + s)/2 mod H` with `s = ±1` (odd).
- Since `q H/2 ≡ H/2 (mod H)`, the options of `x_1` are those of `x_0` shifted by `H/2`.
- Max chooses `b` at each pair; Min chooses `s` at each odd node.

**Proposition P (pair-form certificates).**
- **(lower)** Let `W ⊆ Z/H` be nonempty, `tau : W -> {0,1}` and `g : W -> Z`. Suppose that for every `P in W` and every option `Q` of `x = P + tau(P) H`, we have `Q in W` and `g(Q) <= g(P) + e(x)`. Then every sign strategy has a cycle of density `>= F`, so `rho*(q,k) >= F`.
- **(upper)** Let `sigma` be a sign strategy and `Psi : Z/H -> Z`. Suppose that `Psi(Q) + e(x) <= Psi(x mod H)` for every node `x`, where `Q` is the target pair of `x` under `sigma`. Then `rho_max(sigma) <= F`.

*Proof.*
- *(lower)* Fix `sigma` and `P_0 in W`, and let `P_(i+1)` be the `sigma`-option of `x_(P_i) = P_i + tau(P_i) H`.
  - Every `P_i` lies in `W`, and `x_(P_i) -> x_(P_(i+1))` is an edge of `G_sigma`, since `x_(P_(i+1))` is a lift of the target pair of `x_(P_i)`.
  - The sequence is eventually periodic. Summing `g(P_(i+1)) - g(P_i) <= e(x_(P_i))` around its cycle gives `e`-weight `>= 0`, i.e. density `>= F`.
- *(upper)* Any cycle `x_0 -> x_1 -> ...` of `G_sigma` has `x_(i+1) mod H` equal to the target pair of `x_i`, so `Psi(x_(i+1) mod H) + e(x_i) <= Psi(x_i mod H)`. Summing gives `e`-weight `<= 0`. ∎
- *Remark.* These are THM-4486's Lemmas G2/G1 with the potential of a node taken to depend only on its pair. The least fixed points of the two energy operators give such potentials: `g(P) = min_b f(x_b)` and `Psi(P) = max_b f(x_b)`. So nothing is lost at the value.

## 2. The lean solver and the independent checker

- **Solver** (`procgen_seven_20260926_game.c`, mode `lsearch2`). It runs a Stern–Brocot search on `F` between Farey neighbours, skipping mediants above a known upper bound. At each mediant it computes the least fixed points of the two pair-form energy operators:
  - MAX: `g(P) = max(0, -e(P) + min_b max_(Q option of x_b) g(Q))`;
  - MIN: `G(P) = max(0, e(P) + max_b min_(Q option of x_b) G(Q))`.
- **Engine details.**
  - Values are `uint16`, with a cap; exceeding the cap means TOP.
  - The fixed points are computed on the representatives `min(P, H - P)` of the negation orbits. The negation `x -> -x` maps the two lifts of `P` onto the two lifts of `-P` and the option set of `x` onto minus the option set of `-x`. So both operators commute with it, and their least fixed points from 0 are negation-invariant.
  - Iteration is by in-order sweeps over a dirty bitset, so memory is about `2.1` bytes per representative and state: 142 MB for the two-state search at `k = 27`, 285 MB for the one-state lower bound at `k = 29`.
- **Search.** The two tests are dovetailed with doubling work budgets; the first to finish decides the side (`rho* >= m`, `< m`, `<= m` or `> m`). The other test then gets a bounded confirmation budget. If both certificates exist, `m` is the value. An endpoint left open by a one-sided verdict is re-tested to completion when the search has moved towards it three times.
- **Safety.** Nothing the solver reports is used as a proof step.
  - Caps only turn finite values into TOP, so a small cap, a wrong hint or a wrong one-sided verdict can only make the search fail. It cannot produce a value lacking two valid certificates.
  - The value's certificates are written in a compact format: lift bits and `uint8/uint16` potentials, `tau` and `g` at the argmin, `sigma` at the argmin with ties to `+`.
- **Checker** (`procgen_seven_20260926_verify.c`). It is written separately, uses no symmetry, walks all `2^k` nodes, and checks Proposition P's inequalities exactly. It reports `RESULT CERTIFIED fn/fd` iff both certificates hold.
- **Cross-checks** (runner §A).
  - 20 values of THM-4486's table (q = 5, 7, 9, 11; `k = 8..18`) are reproduced. 6 fresh instances (q = 7, 13, 15, 17, 21; `k = 11..13`) agree with the floor lane's engine.
  - Karp's algorithm (the floor lane's `rho_max_exact`, and its C Karp on `G^tau`) re-evaluates 5 certified strategy pairs at `k = 8..10`: `rho_max(sigma)` equals the value, and Min's best density against `tau` is `>=` the value.
  - The floor lane's Python verifiers accept the same certificates converted to node form.
  - Corrupting one potential entry at a tight edge makes the checker fail, and the checker rejects the neighbouring value `5/13`.
- **One-sided modes.** `lower` and `upper` compute a single energy fixed point at a given `F` and write only that half of the certificate. They are used where the exact search is too slow.
  - At `k = 28` the mediants near the value took over 20 minutes each. An exploratory search (42 min) found `13/35`, and the runner certifies it by a lower and an upper certificate at `13/35`.
  - At `k = 29` a lower certificate at `13/35` pins the value, since `rho*` is non-increasing.
  - At `k = 30` only a lower certificate with `uint8` potentials fits the memory cap: solver 546 MiB, checker 577 MiB (605 MB).
- **Cost** (final run).
  - At `k = 22` the full search takes about 8–18 s and under 10 MB (the floor lane: 230 s, 480 MB).
  - `k = 27`: 426 s and 140 MB; its checker takes 1 s and 145 MB (the target pairs of consecutive pairs lie in a few sequential streams, so the check is memory-bandwidth bound).
  - The one-sided certificates at `k = 28, 29, 30` take 11–72 s each. Their checkers take up to 577 MiB (`k = 30`).
  - Timings and RSS for every entry are in the .out.

## 3. Exact values (runner §B)

`rho*(q,k)` at the levels where it changes (it holds until the next listed level); bold entries are new.

| `q` | `log_q 2` | values | last `k` |
|---|---|---|---|
| 7 | 0.35621 | 3/7 (7), 2/5 (10), 15/38 (14), 7/18 (16), 19/49 (18), 13/34 (19), 34/89 (20), 8/21 (21), 14/37 (22, **23**), **3/8 (24, 25, 26)**, **35/94 (27)**, **13/35 (28, 29)**; **`rho*(7,30) >= 7/19`** | 29 |
| 9 | 0.31546 | ..., 12/31 (17), 5/13 (18), **8/21 (19)**, **17/45 (20)**, **3/8 (21)**, **19/51 (22)**, **27/73 (23)**, **18/49 (24)** | 24 |
| 11 | 0.28906 | ..., 29/77 (17), 41/110 (18), **7/19 (19)**, **4/11 (20)**, **71/197 (21)**, **29/81 (22)**, **11/31 (23)**, **31/88 (24)** | 24 |

- Every listed value or lower bound `a/p` has `q^a > 2^p`, so **class (i) is empty at all these levels**. For `q = 7` this covers every `k <= 30`, the lower levels by THM-4486 and monotonicity.
- **Distance to the threshold** at the last exact level: `q = 7`: `0.0152` (`k = 29`, `13/35 - log_7 2`); `q = 9`: `0.0519`; `q = 11`: `0.0632`. Among `q = 7, 9, 11`, 7n±1 is the closest to provability in absolute terms, although its proved floor `1/3` is the farthest below its threshold.
- **The Fibonacci-like denominators** (`13/34`, `34/89`, `8/21` at `k = 19..21`) are not structural.
  - The values between `k = 14` and `k = 23` all lie in `(3/8, 2/5)`, between `2/5` and its Stern–Brocot left child `3/8`. The Fibonacci ratios are simply points of small denominator in this interval.
  - `14/37` (the "break at `k = 22`") is followed by `14/37` at `k = 23` and the plateau `3/8` at `k = 24..26`. Then the value leaves the interval: `35/94 < 3/8` at `k = 27`, and `13/35` at `k = 28, 29`.
  - Continued fractions: `3/8 = [0;2,1,2]`, `35/94 = [0;2,1,2,5,2]`, `13/35 = [0;2,1,2,4]`.
  - `log_7 2 = [0;2,1,4,5,4,...]`. Values below it need a third partial quotient `>= 4`. The sequence has moved from `[0;2,1,1,...]` (`k <= 23`) to `[0;2,1,2,...]` (`k >= 24`).

## 4. The 3/8 plateau and its end (runner §C; EMPIRICAL facts about certified objects)

- **C1. The optimal strategy keeps a rational cycle of denominator 87.** At `k = 24` the certified optimal `sigma` contains the 7x+1 cycle through `x0 = -169/87`:
  - the cycle is `(-169, -548, -274, -137, -436, -218, -109, -338)/87` with signs `+, +, +` at its three odd points;
  - its shape is `(8,3)`, density `3/8`, and `2^8 - 7^3 = -87`. In THM-4484's terms it is a cycle of 7x+87, for which the shape `(8,3)` is free;
  - every one of its 8 edges is tight for the certified potential `psi`.
- **C2. Min's best responses to the certified Max strategy.**
  - The tight core of the lower certificate has 38044 pairs. Its cycles, which are exactly the Min responses of density `3/8`, have shapes `(8,3), (24,9), (32,12), (40,15), (48,18), (56,21), (64,24)`.
  - The `(8,3)` ones have denominators dividing 87. All are real-bounded rationals (`|x| <= 170`).
- **C3. Census of small-shape rational cycles kept by the certified optimal strategies.** A cycle is kept iff `sigma` agrees with its sign at the residue of each odd point.

| shape `(p,a)` | density | expanding? | simple cycles | kept, `k = 24` | kept, `k = 27` |
|---|---|---|---|---|---|
| (2,1) | 1/2 | yes | 2 | 0 | 0 |
| (5,2) | 2/5 | yes | 8 | 0 | 0 |
| (7,3) | 3/7 | yes | 38 | 0 | 0 |
| (13,5) | 5/13 | yes | 3164 | 0 | 0 |
| **(8,3)** | **3/8** | yes | **54** | **4** | **0** |
| (11,4) | 4/11 | yes | 480 | 61 | 33 |
| (14,5) | 5/14 = 0.35714 | yes (barely) | 4568 | 374 | 323 |

- **Reading.** Every cycle denser than the value is broken, as it must be. On the plateau the optimal strategy can afford some `(8,3)`-cycles; one level of lookahead later it breaks all 54.
  - A class-(i) strategy would have to break all 4568 `(14,5)`-cycles (`7^5 = 16807 > 2^14 = 16384`), all 480 `(11,4)`-cycles, and infinitely many longer expanding ones.
  - The census shows what "one more bit of lookahead" buys. It says nothing about whether this can be completed.
- **At `k = 27`** (random tight walks in both certificates, exploratory, not in the runner): all critical cycles found (19 against `tau`, 5 in `G_sigma`) have shapes `(94,35)` or `(188,70)`.
  - Their periodic points are rationals like `-53911781070525928505242544720/359011651637098697284331638359 = -0.150...`, with denominators dividing `7^35 - 2^94 = 359011651637098697284331638359` (resp. `7^70 - 2^188`).
  - Along the orbits `|x| < 120`. The residues of the orbit points reproduce the cycle nodes exactly.

## 5. Proposition S: the rational relaxation, and what an all-level proof would need

**Proposition S.** Let `q >= 3` be odd. For a nonzero rational `x` with odd denominator put `sigma_inf(x) = sgn(x)`. Then:
- **(a)** every cycle `x_0, ..., x_(p-1)` (`p >= 1`, `a >= 1` odd steps) of `T(x) = x/2` (`x` even), `(q x + sgn x)/2` (`x` odd) satisfies `2^p = prod_(odd i) (q + 1/|x_i|) > q^a`, so its density is `< log_q 2`;
- **(b)** for every finite family `C` of rational sign-choice cycles of density `>= log_q 2`, and every level `k` at which the points of the members of `C` have pairwise distinct residues mod `2^k`, some level-`k` sign strategy contains no member of `C`.

*Proof.*
- (a) For odd `x`, `q x + sgn(x) = sgn(x)(q|x| + 1)`, so `T(x)` has the sign of `x` and `|T(x)| = (q|x| + 1)/2`. For even `x`, `|T(x)| = |x|/2`.
  - Around a cycle, `1 = prod |x_(i+1)|/|x_i| = prod_(odd) (q + 1/|x_i|) / 2^p`. An all-even cycle is `{0}`.
  - Hence `2^p > q^a` when `a >= 1`.
- (b) Each `c in C` has density `>= log_q 2`, so by (a) it is not a cycle of `T` with the signs `sgn`. So `c` requires the sign `-sgn(x_c)` at some odd point `x_c`. Give the residue of `x_c` the sign `sgn(x_c)`, which is consistent because the residues are distinct. Choose the other signs arbitrarily. ∎

**Consequences.**
- **What an all-level certificate cannot be.** A proof that 7n±1 is never provable cannot consist of a finite list of expanding rational cycles that every strategy must keep, the analogue of q = 5's sporadic cycle `(1,3,8,4,2)`. That cycle is contracting, and indeed `sgn` keeps it.
  - For every finite list, `sigma = sgn` on its points breaks it at all large levels.
  - The obstruction, if it exists, lives in residue collisions: rationals `x != y` with `x ≡ y (mod 2^k)` and opposite real signs, with growing denominators as `k` grows.
- **The data fit this.** Long critical cycles with huge denominators appear at every computed level.
  - `k = 20`: shape `(89,34)`, denominator dividing `7^34 - 2^89`.
  - `k = 27`: shapes `(94,35)` and `(188,70)`, `7^35 - 2^94 ~ 3.6 * 10^29`.
  - On the plateau (`k = 24`) the short `(8,3)`-cycles (denominator 87) occur together with long ones.
- **Min's point of view.** A level-`k` strategy is continuous on `Z_2` and cannot know the real sign. `rho*(7,k) - log_7 2` measures how well a `k`-bit rule can imitate `sgn` on the expanding cycles.
- **Checked** (runner C4). Among all 24653 expanding simple cycles of 7x±1 with `p <= 11`, `sgn` keeps none. It keeps 132 of the 395 contracting ones, e.g. the free cycle `(1,4,2)`.

## 6. Theorem Q: q = 9 and q = 11 at every level

- **The graph.** `U_k` (THM-4486 §5) has nodes `u = 1..H` and edges `u -> u/2` (even `u`) and `u -> rho((qu ± 1)/2)` (odd `u`), where `rho(m) ∈ [1, H]` and `rho(m) ≡ m (mod H)`.
- **The reduction.** `U_k` is the Min graph of the top-lift adversary. So a bound "every cycle of `U_k` has density `>= F`" gives `rho*(q,k) >= F` (Lemma G2 with the real potential `log Phi`, THM-4486 §6).

**Theorem Q.** Let `(q, F, beta, U0) = (9, 5/16, 5, 16384)` or `(11, 2/7, 2, 1024)`. Write `F = a0/p0`, `B = 2^-beta`, `A = 2^(beta(p0 - a0)/a0)`; this is `2^11` for `q = 9` and `2^5` for `q = 11`, so `A^a0 B^(p0-a0) = 1`.
- Let `Phi(u) = u^beta` for `u >= U0`. On `1 <= u < U0` let `Phi` be the least positive solution of `(E) Phi(u/2) <= B Phi(u)` (u even) and `(O) Phi((qu ± 1)/2) <= A Phi(u)` (u odd). This is an exact rational fixed point, computed by the runner.
- Then for every `k` with `H = 2^(k-1) >= (q U0 + 1)/2` (i.e. `k >= 18`, resp. `k >= 14`), every cycle of `U_k` has density `>= F`.
- Consequently `rho*(q,k) >= F` for **every** `k`.

*Proof.*
- **The finite checks.** The runner checks exactly (runner §D):
  - C1: (E), (O) on `u < U0`, with the unreduced successors and `Phi(m) = m^beta` for `m >= U0`;
  - C2: `Phi(w) <= w^beta` on `U0/2 <= w < U0`;
  - C3: `(q U0 + 1)^beta <= A 2^beta U0^beta`;
  - C4: `max_(u<U0) Phi(u) <= A U0^beta`;
  - C5: `Phi > 0`.
- **Every edge `u -> v` of `U_k` satisfies `Phi(v) <= lambda Phi(u)`** with `lambda = B` (even `u`) and `A` (odd `u`):
  - *`u < U0` even:* C1.
  - *`u < U0` odd:* `m = (qu ± 1)/2 <= (q(U0-1)+1)/2 < H` and `m >= 1`, so `v = m`. Then C1 applies.
  - *`u >= U0` even:* if `u/2 >= U0`, then `Phi(u/2) = B u^beta`. Otherwise `U0/2 <= u/2 < U0`, and C2 gives `Phi(u/2) <= (u/2)^beta = B Phi(u)`.
  - *`u >= U0` odd:* `v = rho(m) <= m <= (qu+1)/2`.
    - If `v >= U0`, then `Phi(v) = v^beta <= ((qu+1)/2)^beta`. This is `<= A u^beta` by C3, because `(qu + 1)/u` decreases in `u`.
    - If `v < U0`, then `Phi(v) <= A U0^beta <= A u^beta` by C4.
- **A cycle.** For a cycle with `a` odd steps among `p`: `1 = prod Phi(next)/Phi(cur) <= A^a B^(p-a) = 2^(beta((p0-a0)a/a0 - (p-a)))`. Hence `p0 a >= a0 p`, i.e. `a/p >= F`.
- **Every level.** Lemma G2 gives `rho*(q,k) >= F` for `k >= 18` (resp. 14). For smaller `k`, `rho*(q,k) >= rho*(q,18)` (resp. `rho*(q,14)`) because `rho*` is non-increasing in `k`. ∎

**Remarks.**
- **Exactness.** The cycles `(1, 5, 22, 11, 50, 25, 112, 56, 28, 14, 7, 32, 16, 8, 4, 2)` (9u∓1; `5/16`, `9^5 = 59049 < 2^16`) and `(1, 6, 3, 16, 8, 4, 2)` (11u∓1; `2/7`, `121 < 128`) lie in `U_k` without wrap-around for `k >= 8` (resp. 5). So Theorem Q is the exact value of the negative-integer adversary at every `k >= 18` (resp. 14).
  - This parallels THM-4486 Corollary 5 for q = 5. The difference: for q = 5 the adversary's value `2/5` is also `lim rho*`, while for q = 9, 11 the certified `rho*` are still far above it (`18/49` and `31/88` at `k = 24`).
- **Why these U0.** `max Phi(u)/u^beta` over the corrected range is `0.9312` for q = 9 (`0.9683` for q = 11).
  - Smaller `U0` fail C2. For q = 9 with `U0 = 4096`, `Phi(3072)/3072^5 = 1.033`; the step `1 -> 5` has ratio `5^5/2^11 = 1.53 > 1`, and it propagates to `1, 3, 7`.
  - For q = 13 the same scheme (`17/63`, `beta = 17`) still fails C2 at `U0 = 2^20` (exploratory), so q = 13 is not claimed.
- **Neither floor reaches `log_q 2`** (`9^5 < 2^16`, `11^2 < 2^7`). The gaps are `0.0030` and `0.0034`. So Theorem Q narrows the class-(i) window of 9n±1 to `[5/16, log_9 2) = [0.3125, 0.3155)`, and of 11n±1 to `[0.2857, 0.2891)`.

## 7. Route (A) for q = 7: what was tried

The goal: at every level, an adversary with value `> log_7 2`, together with a checkable potential.
- **The top lift (Theorem N) is exactly `1/3`.** The free cycle `(-1, -4, -2)` of 7x−1 is always available to Min.
- **Switching the top lift at the free cycle** (the orchestrator's suggestion; runner §E). Take the top lift everywhere except at the class of `-1` (resp. `-1, -2, -4`), where the other lift is taken. This avoids `(-1,-4,-2)`, but the values at `k = 8, 10, 12` are `7/24, 7/24, 5/19` (resp. `1/3, 7/24, 11/38`), below `1/3` from `k = 10` on: the adversary gets worse.
  - Exploratory values at `k = 14`: `2/7` and `9/31`.
  - Adding the shadow classes `-1 - H/2^j` (`j <= 5`) does not help either: `9/29, 9/29, 5/19, 2/7` at `k = 8..14`.
- **Shifted windows** (runner §E). The adversary plays `-x` in a window `[aH, (a+1)H)` for `a = 7/8, 3/4, 1/2, 1/4, 1/8`. Values at `k = 8, 10, 12` are at most `1/3`:
  - `a = 1/2`: exactly `1/3`;
  - `a = 1/4`: `1/8, 1/10, 1/12`;
  - `a = 7/8` and `a = 1/8`: `3/10, 1/4, 1/4`;
  - random lifts: `2/7, 2/7, 1/4`.
  - The floor lane's rational windows and exact traps also stayed at or below `1/3`.
- **Why height potentials cannot work** (heuristic, supplementing THM-4486 §7). Any adversary with value `> 1/3` must break both free cycles `±(1,4,2)`. It must lift one of the classes of `±1, ±2, ±4 mod H` to the other lift, at height `~H`. Such a node is 2-adically within `2^-(k-1)` of the free cycle.
  - Whatever Max does afterwards, Min can shadow the free cycle for about `(k-1)/3` rounds. The congruence depth drops by one per step: Max's lifts act only on the top bit, so no lift can shorten the shadowing.
  - In height terms this is a free ascent of about `k` bits, after which Min may descend along halving-rich integer paths.
  - So a potential that is a function of height alone cannot certify such an adversary above `1/3`. A branching-random-walk estimate (descents at the entropy density `p0 = 0.227`) even suggests values near `1/4`, consistent with the switched top lifts and the windows above (`1/4` to `0.29` at `k = 10, 12`). This is heuristic.
  - The certified optimal adversaries survive only because their potentials are 2-adic, not height-based. At `F = 3/8` the maximal Max potential is exactly `3(k-1)`: `69, 72, 75` at `k = 24, 25, 26`. It sits at the pair 0, whose lifts force `k - 1` halvings.
  - An all-level certificate would need a uniform description of such 2-adic potentials. Proposition S shows it cannot come from finitely many rational cycles.
- **Inter-level structure of the certificates** (exploratory). At `F = 3/8`, `k = 24 -> 25 -> 26`:
  - `g_(k+1)(P) - g_k(P mod 2^(k-2))` is `0` on 52% (resp. 35%) of the pairs and `1..3` elsewhere; no projection to high bits matches.
  - 62% of the pairs have a forced lift, split exactly half/half, with no dependence on the top bits.
  - The Min potential on the plateau is bounded (`max psi = 21`), consistent with a lifted strategy; the Max potential is not.

## 8. Route (B): the trend (EMPIRICAL)

| `k` | 14 | 16 | 18 | 19 | 20 | 21 | 22 | 24 | 27 | 28, 29 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `rho*(7,k)` | 0.39474 | 0.38889 | 0.38776 | 0.38235 | 0.38202 | 0.38095 | 0.37838 | 0.37500 | 0.37234 | 0.37143 | `>= 0.36842` |
| minus `log_7 2` | 0.0385 | 0.0327 | 0.0316 | 0.0261 | 0.0258 | 0.0247 | 0.0222 | 0.0188 | 0.0161 | 0.0152 | `>= 0.0122` |

- The decrease is irregular: plateaus of 1–3 levels, drops of `0.0003–0.005`. Over `k = 20..29` the mean drop is about `0.0012` per level.
- A linear extrapolation reaches `log_7 2` near `k ~ 42`. For comparison, q = 5 converged to its adversary value by `k = 15`.
- The three odd `q` behave alike: values near `0.35–0.37` at `k = 24`, decreasing. Their proved floors (`1/3, 5/16, 2/7`) are far below.
- **Limits of computation.**
  - The exact search at `k = 28` is slow: the mediants near the value took over 20 minutes each.
  - The solver needs about `2.1 * 2^(k-2)` bytes per state (546 MiB at `k = 30`). The checker needs `1.1 * 2^(k-1)` bytes with `uint8` potentials (577 MiB at `k = 30`), or `2.1 * 2^(k-1)` with `uint16`.
  - So `k = 30` admits only lower certificates with small potentials, and `k = 31` exceeds the 700 MB cap.
- **Verdict.** The data are compatible with a limit anywhere in `[1/3, 13/35]`. Whether 7n±1 becomes provable at some level is OPEN.
- **Suggested next steps** (suggestions, not claims).
  - Route A: look for Max certificates of a restricted form. The lift bit would be computed by a finite automaton reading the pair's binary expansion, and the potential would be a multiple of a 2-adic valuation plus an automaton-computed bounded correction. At `F = 3/8` the maximum `3(k-1)` sits at the pair 0 and comes from its forced halvings, which suggests this shape.
    - For such certificates the edge inequalities at all levels at once form a sentence of Büchi arithmetic `(N, +, V_2)`, which is decidable.
    - The inter-level data of §7 (partial agreement of `g_(k+1)` with `g_k` on low bits) are a starting point.
  - Route B: a disk-backed streaming version of the sweep engine could reach `k ~ 32` and show whether the decrements keep shrinking.

## 9. Failures and caught mistakes

- **A wrong structural claim, caught by the runner.** From an early exploratory run (SCC shortest cycles `32, 40`) I asserted that Min's best responses at `k = 24` have no `(8,3)` cycle. The runner's BFS from 3000 core nodes found `(8,3)` cycles with denominator 87. The check was corrected to state the shapes found (§4, C2).
- **A misused known upper bound.** A test run passed `3/8` as the known upper bound at `k = 22`, where the value is `14/37 > 3/8`. The search correctly failed without certifying anything, which illustrates the safety argument of §2.
- **The plateau guess.** From `3/8` at `k = 24, 25, 26` (a compact Max potential of maximum `3(k-1)`, a Min potential of maximum `21`) I expected `3/8` to be final, as `2/5` is for q = 5. A lean test at `F = 3/8` for `k = 27` failed on the MAX side, and the search then certified `35/94`. REFUTED, and recorded.
- **Engine.**
  - The first search (full MIN completion after every MAX win) spent 150–400 s per mediant at `k = 27` on losing climbs. It was replaced by one-sided verdicts with confirmation budgets and the endpoint guard (§2).
  - A `tauval` mode with large caps was too slow and was abandoned for the floor lane's `value_of_tau` (§E).
- **Process hygiene.** An exploratory `pgrep`-style kill matched its own shell (the pattern occurred in the heredoc text) and killed it. No data were lost; later kills used explicit PIDs.
- **q = 13:** the corrected-potential scheme failed (C2 not satisfied up to `U0 = 2^20`); no claim is made.
- **Nothing in the canon is contradicted.**
  - THM-4486's values for `k <= 22` are reproduced exactly. Its FINITE-EXACT range for 7n±1 (`k <= 22`) is extended: exact values to `k = 29`, `rho*(7,30) >= 7/19`, and no provable strategy at any `k <= 30`.
  - Its OPEN question (`lim rho*(7,k)` vs `log_7 2`) remains OPEN.

## 10. Reproduction

```bash
/usr/bin/time -l python3 -u 04-computation/experiments/procgen_seven_20260926_run.py > 05-knowledge/results/procgen_seven_20260926.out
```

- **Requirements.** `numpy`, `scipy`, a C compiler.
  - The C programs are compiled into `scratch/procgen_seven/run/` with the source hash in their names.
  - The floor lane's `procgen_floor_20260926_lib.py` builds its engine into `scratch/procgen_floor/`.
  - Certificates are written to `scratch/procgen_seven/run/` and deleted after checking. The `k = 24` and `k = 27` certificates are kept until §C has used them.
- **Environment.** `SEVEN_KMAX7` (27, the last level of the exact search chain), `SEVEN_KBR7` (1: the one-sided certificates at `k = 28, 29, 30`), `SEVEN_KMAXQ` (24), `SEVEN_KSTRUCT` (24).
- **Hints.** The runner passes the Farey bracket `[16/43, 3/8]` for `k = 27` (from an exploratory run) to shorten the search. The thresholds `13/35` (`k = 28, 29`) and `7/19` (`k = 30`) come from exploratory runs. A wrong hint can only make a computation fail.
- **Output.** 118 checks, ending with `ALL CHECKS PASSED`.
- **Cost.** One runner process plus one child at a time. Wall time 1209 s (about 20 min) on the shared 8-core machine, with load about 3–5 from other lanes.
  - Peak child RSS 577 MiB: the checker of the `k = 30` lower certificate. Runner RSS 278 MiB.
  - `/usr/bin/time -l` of the whole run reports a maximum resident set of 605,110,272 bytes. That is the same checker, below the 700 MB cap.
- **SHA-256** (raw bytes):
  - `procgen_seven_20260926_game.c` `4d2a193f96126ca5c13e4377c9d9018da7c25034cf0feda1152198750e5d7b34`
  - `procgen_seven_20260926_verify.c` `53576fbcf5472a989055623af881defafe9c16b54df28fb2b836989c54cdef39`
  - `procgen_seven_20260926_tight.c` `1ff6968c4d7509d079259f4f7eb771309dd9078df0f93ed07d439f1fda9efb14`
  - `procgen_seven_20260926_potential.py` `f322b3bb001d6fa05c0eb9f041b2d6e98290baa473f116d1e795cc9e559b9682`
  - `procgen_seven_20260926_cycles.py` `a7a66ac72b86c9ee1083aca3dfa6d1e3d5cbf8e1245fa3ab18127ee5d814daae`
  - `procgen_seven_20260926_run.py` `8325e6d49c5de577b7dde1a4f5cd7f202511cd6fada187281850d109a4dfe0fc`
  - `procgen_seven_20260926.out` `f48efef142e902aa6bd88cd84c270828d3b8275b8744a9e01c0c11369d6fae96`
  - reused read-only: `procgen_floor_20260926_lib.py` `be7024fa7dc147c21d3eacabfc59d09c558943dcfa3cb30b421abda79e3e4f2f` (unchanged from THM-4486)
- **Timing and RSS lines** in the .out vary between runs; values, digests of certificates and all checks are deterministic.
