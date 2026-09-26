# The structure of the optimal 5n±1 and 7n±1 sign strategies: the max-halving skeleton, a flip calculus, explicit optimal rules, and the levels 30–31

**Status.**
- **Answer to the lane question (is 7n±1 ever bounded-lookahead provable, i.e. is `lim rho*(7,k) < log_7 2`?): UNDECIDED.** Neither direction was proved. What is new:
- **PROVED** (hand proofs in §§1–3; the runner checks their finite content):
  - **Lemma C (cycles are periodic points).** The closed walks of `G_sigma` correspond one-to-one to the periodic points of the 2-adic map `T_sigma : Z_2 -> Z_2`, with the same odd density. So `rho_max(sigma)` depends only on the map `T_sigma`: a rule that reads `d` digits has the same `rho_max` at every level `k >= d`. This makes an explicit finite rule a level-independent object.
  - **Lemma MH (the max-halving skeleton).** Max-halving (MH: the sign with `q x + s = 0 mod 4`) has accelerated steps `x -> (q x + s)/2^v` with `v >= 2`. Its itinerary `((s_i, v_i))` is a bijection onto all sequences in `{±1} x {2,3,...}`, each symbol using `v_i` bits. `rho_max(MH) = 1/2` at every level, attained exactly on the periodic points of the set `S_inf` where every valuation is 2. `S_inf` is a full 2-shift (sign sequences). It has exactly `2^n` periodic points of period dividing `n`, rationals in `[-1/(q-4), 1/(q-4)]` (q >= 5), with fixed points `±1/(q-4)`: `±1` for q = 5, `±1/3` for q = 7. Every strategy with `rho_max < 1/2` must change the MH sign on every periodic orbit of `S_inf`.
  - **Lemma F and Corollary F (the flip calculus, q = 7; the q = 5 analogue).** Flipping the MH sign at `x` gives `y = 2^(v_1 - 1) x_1 - s_1`. Its MH valuation is determined by the first three itinerary symbols of `x` (one sub-case, `>= 7`, depends on more). Consequence: a flip followed by max-halving gains halvings over two max-halving steps **exactly** when the itinerary starts `(s,2),(s,2)`. It is density-neutral exactly on `(s,4),(-s,.)` (the orbit rejoins), and it loses otherwise. Modulo 32 the gaining pattern is exactly the classes 11 and 21, i.e. those of `1/3` and `-1/3`. At the fixed points the flip gives:
    - for q = 7, the valuation pair `(1,5)`: the cycle `(-1,-5,-16,-8,-4,-2)/3` of density `1/3`;
    - for q = 5, the pair `(1,4)`: the sporadic cycle `(1,3,8,4,2)` of density `2/5`, which is exactly the q = 5 limit.
- **FINITE-EXACT** (every entry carries certificates re-checked by an independent exact checker; runner §C, §D, §G):
  - **Levels 30 and 31 (new):** **`rho*(7,30) = 37/100` exactly** (mean valuation `100/37 = 2.703`), and `7/19 <= rho*(7,31) <= 37/100`. Hence **7n±1 has no provable sign strategy at any level `k <= 31`.** The value drops again at `k = 30`: `37/100 < 13/35 = rho*(7,29)`.
  - **An explicit optimal 5n±1 rule for every level `k >= 15`.** Max-halving with the sign flipped on 46 residue classes of depth at most 15 (a 122-leaf decision tree on the low bits; §4 lists the classes) has `rho_max = 2/5`. By Lemma C this one finite rule realizes `rho*(5,k) = 2/5` at every `k >= 15`. It flips MH at the `S_inf` fixed points `±1`, turning `(1,2)` into the sporadic cycle, and at the period-2 orbit `±1/9`.
  - **Why 15** (for this rule). Its 8 depth-15 leaves form 4 classes mod `2^14`. Each is split into two halves with opposite signs; the halves contain `19/99 | -163/13` and `-889/1129 | 9/7` (and negatives). Giving a whole class mod `2^14` one sign, either one, creates a cycle denser than 2/5. Among them is the all-outward `(12,5)`-cycle of 5x+971, of density `5/12 = rho*(5,14)`. The two conflict cycles of each class pass through rationals that are congruent mod `2^14` but not mod `2^15`, e.g. `4613/971` and `37501/52947`: a residue collision resolved at level 15.
  - **Explicit 7n±1 rules:** 46 leaves (depth 10) for `2/5`, 300 leaves (depth 14) for `15/38`, 508 leaves (depth 16) for `7/18`, each certified optimal at its level. The 2/5 rule is max-halving flipped on 18 classes (§5). The 46-leaf rule's depth-5 flips are exactly the classes of `±1/3` (the gaining pattern of Corollary F). It also flips `±1/11` (the alternating `S_inf` orbit, where every flip loses). Its densest cycles are the `(5,2)`-cycles of 7x±1 with denominator 17, which are flip + valuation 4.
- **EMPIRICAL** (runner §D; items marked "exploratory" in §§5–6 are not in the runner):
  - **Min.** The certified optimal Min potentials admit the MH sign at 85–88% of the odd residues (k = 19..27). The exceptions sit where MH has valuation 2 (x = 3, 5 mod 8: 24–27%, against 2–3% at x = 1, 7 mod 8) and on the classes of `±1/3` (55–61%).
  - **Nesting.** Across consecutive levels with the same value, the classification forced MH / forced flip / free agrees with the lift on 94–96% of the residues, with at most 0.24% direct conflicts. The optima are refined, not rebuilt.
  - **Max.** The optimal Max lift strategies show no dependence on the real height: the forced-lift shares vary by less than 0.01 across height bins.
  - **Critical cycles** (the tight cycles of both certificates, k = 20–24, exploratory) almost all contain flips (valuation 1), and use valuations 1..10.
  - **Mean valuation.** `1/rho*` rises from 2.5 (k = 10) through `phi^2 = 2.618` (k = 20–21, the Fibonacci values) to `100/37 = 2.703` at k = 30. The threshold is `log_2 7 = 2.807`.
- **REFUTED / negative (FINITE-EXACT at the listed levels):**
  - *Real-sign imitation.* The rule "sign of the simplest rational in the residue class" has `rho_max = 1` (all-odd cycles) for q = 5, 7 at k = 10..20. Correcting MH by it on classes of height `<= h` is never better than MH.
  - *Small automata.* Flip sets recognised by automata with 2 or 3 transient states (all), or 4 (20000 random), never beat MH (`1/2`) at k = 12.
  - *Low-bit adversaries.* `tau(P) = f(P mod 2^m)` never beats `1/3` at k = 10 (m <= 3) and k = 12 (m <= 2). Hill-climbing over m = 6 at k = 12 did not beat `1/3` either (exploratory).
- **OPEN:** `lim_k rho*(7,k)` versus `log_7 2`. Proved floor `1/3` (Theorem N); certified `rho*(7,30) = 37/100` and `rho*(7,31) >= 7/19`.
- No HYP or THM file was created. Collatz is untouched.

Session `collatz-procgen-20260922`, lane **seven2**, 2026-09-26.
- **Scripts** (`04-computation/experiments/`): `procgen_seven2_20260926_{rhomax,restrict,lean8,verify8}.c` and `procgen_seven2_20260926_{lib,run}.py`.
- **Reused read-only:** `procgen_seven_20260926_{game,verify,tight}.c` and `procgen_floor_20260926_lib.py`.
- **Output:** [procgen_seven2_20260926.out](procgen_seven2_20260926.out).
- **Parents:**
  - [THM-4486](../../01-canon/theorems/THM-4486-min-max-cycle-density-game.md) (the game, Theorem N, Corollary 5, Proposition S);
  - the floor note [procgen_floor_20260926_density_floor.md](procgen_floor_20260926_density_floor.md);
  - the seven note [procgen_seven_20260926_seven_n_plus_one_provability.md](procgen_seven_20260926_seven_n_plus_one_provability.md).

## 0. The answer in brief

| question (lane brief) | answer | status |
|---|---|---|
| decide `lim rho*(7,k)` vs `log_7 2` | not decided; no provable strategy at `k <= 31`; `rho*(7,30) = 37/100`, `rho*(7,31) >= 7/19` | OPEN; FINITE-EXACT |
| (1) structure of the optima, q = 5 | **one explicit 122-leaf rule is optimal at every `k >= 15`**. It is max-halving plus flips: first at the `S_inf` fixed points `±1`, where the flip creates the sporadic cycle, then in shrinking 2-adic neighbourhoods. Depth 15 is needed to split 4 classes mod `2^14` (2 negation pairs) that contain pairs of rationals needing opposite signs. | PROVED (Lemma C) + FINITE-EXACT |
| (1) structure, q = 7 | the same skeleton. Max-halving on 85–88% of residues; flips concentrated where max-halving has valuation 2; the first flips exactly on `±1/3` (the unique gaining pattern) and on the alternating orbit `±1/11`. Explicit rules of 46 / 300 / 508 leaves for `2/5`, `15/38`, `7/18`. | PROVED lemmas + FINITE-EXACT + EMPIRICAL |
| (1) self-similarity, renormalization, Fibonacci | stage 0 = MH (`1/2`, obstruction `S_inf`); stage 1 = flip the gaining pattern and the alternating orbit (`2/5`, obstruction: flip-then-valuation-4 cycles of denominator 17); later stages are refinements (94–96% nested) whose obstructions are long flip-containing cycles with huge denominators. The Fibonacci values are `1/rho*` crossing `phi^2`; no Sturmian or Fibonacci renormalization was found. | EMPIRICAL |
| (2) automatic strategies | real-sign imitation is catastrophic (`rho_max = 1`); small automata never beat max-halving. A rule reading `d` digits *is* its level-`d` graph (Lemma C). The honest finite-state content is the flip calculus: a single flip's effect is a finite function of three itinerary symbols. | PROVED + FINITE-EXACT |
| (3) better adversaries | none above `1/3` among low-bit adversaries (m ≤ 3 at k = 10; hill-climbed m = 6); the certified optimal adversaries are height-independent; their critical cycles contain flips | EMPIRICAL |
| (4) computation | **`rho*(7,30) = 37/100`** and `rho*(7,31) >= 7/19`, via the seven engine (16-bit), a one-byte engine (with an offset `g - fn·v_2`) and a mirror-streaming checker, all under 610 MB. Exact values at `k >= 31` need 16-bit potentials (1 GiB), which do not fit in 700 MB. | FINITE-EXACT |

## 1. Lemma C: cycles of `G_sigma` are periodic points of `T_sigma`

Setting as in THM-4486.
- Level `k`, odd `q`, `N = 2^k`, `H = 2^(k-1)`. A sign strategy `sigma` is a sign on each odd residue mod `N`.
- `T_sigma(x) = x/2` for even `x` and `(q x + sigma(x mod N))/2` for odd `x`, a map of `Z_2`.
- `G_sigma` has nodes `Z/N` and edges `w -> ` both lifts of `T_sigma(w) mod H`.

**Lemma C.**
- (a) If `x_0` is a periodic point of `T_sigma` with orbit `x_0, ..., x_(p-1)`, then `(x_i mod N)` is a closed walk of `G_sigma` with the same parities.
- (b) Conversely, every closed walk `w_0 -> ... -> w_(p-1) -> w_0` of `G_sigma` comes from exactly one periodic point, the one with `x_i ≡ w_i (mod N)`, and `x_0 = c/(2^p - q^a)` is rational.
- (c) Hence `rho_max(sigma)` is the largest odd density of a periodic orbit of `T_sigma`, and it does not change when `sigma` is read at a higher level (lifting).

*Proof.*
- (a) `T_sigma(x) mod H` depends only on `x mod N`. So `x_(i+1) mod N` is a lift of the pair `T_sigma(x_i mod N) mod H`.
- (b) Put `eps_i = w_i mod 2` and `s_i = sigma(w_i)` for odd `w_i`. Let `A_i(y) = y/2` (`eps_i = 0`) or `(q y + s_i)/2` (`eps_i = 1`).
  - The inverse branches `A_i^(-1)(z) = 2z` or `(2z - s_i)/q` are 2-adic contractions by `1/2`, with image in the class `eps_i mod 2`.
  - So `A_0^(-1) ∘ ... ∘ A_(p-1)^(-1)` has a unique fixed point `x_0`. Its orbit `x_(i+1) = A_i(x_i)` is periodic and has the parities `eps_i`.
  - *Induction on `m <= k`: `x_i ≡ w_i (mod 2^m)` for all `i`.* The case `m = 1` holds by construction.
  - *Inductive step.* If it holds for `m <= k - 1`, then `A_i(x_i) = x_(i+1) ≡ w_(i+1) ≡ A_i(w_i) (mod 2^m)`. The last congruence is the edge condition mod `H`, and `m <= k-1`.
  - `A_i(y) - A_i(y')` is `(y - y')/2` or `q(y - y')/2`, so `v_2(x_i - w_i) >= m + 1`.
  - Hence `x_i ≡ w_i (mod N)`, `sigma(x_i) = s_i` and `x_(i+1) = T_sigma(x_i)`. Uniqueness is the uniqueness of the fixed point.
- (c) This follows from (a) and (b), since the densest cycle of the finite graph `G_sigma` is attained. ∎

*Checked:* every witness cycle printed by the runner (hundreds) is re-walked on `G_sigma`, and its rational periodic point is shown to have the walk's residues (`procgen_seven2_20260926_lib.cycle_rational`).

**Consequence.** A sign rule that depends on the lowest `d` bits is one map `T_sigma`, and its `rho_max` is a single number. A certificate at level `d` proves the same value at every level `k >= d`. This is how the explicit rules below are "level-independent optimal families".

## 2. The max-halving skeleton (Lemma MH)

Fix odd `q >= 5`.
- For odd `x`, let `s(x)` be the unique sign with `q x + s ≡ 0 (mod 4)` (MH). For `q ≡ 3 (mod 4)`, `s(x) ≡ x (mod 4)`; for `q ≡ 1 (mod 4)`, `s(x) ≡ -x`.
- Let `v(x) = v_2(q x + s(x)) >= 2`, and `A(x) = (q x + s(x))/2^(v(x))`, the accelerated MH map: one odd step and `v - 1` halvings.

**Lemma MH.**
- (i) *Density.* Every cycle of `G_MH` has density `a / sum v_i <= 1/2`, with equality iff all `v_i = 2`.
- (ii) *Itineraries.* `x -> ((s_i, v_i))_(i >= 0)` (the MH itinerary) is a bijection from the odd 2-adic integers onto `({±1} x {2,3,...})^N`. A prefix `(s_1,v_1),...,(s_n,v_n)` is exactly one residue class mod `2^(1 + sum v_i)`.
- (iii) *The set `S_inf`.* Let `S_inf = {x : all v_i = 2}`. The map `x -> (s_i)` is a bijection `S_inf -> {±1}^N` conjugating `A` to the shift.
  - For each sign word `(s_1..s_n)` there is exactly one periodic point with that period word: the fixed point of the composition of the 2-adic contractions `g_s(y) = (4y - s)/q`.
  - So there are exactly `2^n` points of period dividing `n`. They are rationals `c/(4^n - q^n)` and lie in `[-1/(q-4), 1/(q-4)]`, the attractor of the real contractions `g_s` (ratio `4/q`) with fixed points `-s/(q-4)`.
- (iv) `rho_max(MH) = 1/2` at every level. Every `sigma` with `rho_max(sigma) < 1/2` differs from MH at some point of every periodic orbit of `S_inf`.

*Proof.*
- (i) A cycle of `G_MH` is a periodic MH orbit (Lemma C), a concatenation of blocks of `v_i` steps with one odd step each.
- (ii) Given `x mod 4`, `s_1` is determined. Given `s_1`, the condition `q x + s_1 ≡ 2^(v_1) (mod 2^(v_1 + 1))` fixes `x mod 2^(v_1 + 1)`, and `x_1 = (q x + s_1)/2^(v_1)` mod `2^m` corresponds to `x mod 2^(m + v_1)`. Induct.
- (iii) *Periodic points.* `q g_s(y) + s = 4y`. For odd `y`, `g_s(y)` is odd, its MH sign is `s` (since `4 | q g_s(y) + s`) and its valuation is exactly 2. So every sign word is realised, and the fixed point of a contraction is unique.
  - *Real values.* The real maps `g_s` are contractions (ratio `4/q < 1` for `q >= 5`) with fixed points `-s/(q-4)`. Every periodic point of the IFS lies in its attractor `[-1/(q-4), 1/(q-4)]`.
  - The rational periodic points of the 2-adic and the real dynamics coincide: they are the fixed points of the same affine composition, with the same value `c/(4^n - q^n)`.
- (iv) This follows from (i), (iii) and Lemma C. ∎

*Checked* (runner §B):
- `rho_max(MH) = 1/2` with `S_inf` witnesses, for q = 5, 7, 9, 11 and k = 8, 12, 16;
- exactly `2^n` periodic points for n = 1..8, all in the interval, for q = 5, 7.

## 3. The flip calculus (Lemma F, Corollary F)

A **flip** at odd `x` means using `-s(x)` instead of the MH sign. Let `x` have MH itinerary `(s_1,v_1),(s_2,v_2),(s_3,v_3),...` with MH points `x_1, x_2, x_3`, and let `y = (q x - s_1)/2` be the flipped successor. It is odd, so a flip is an odd step of valuation 1.

**Lemma F (q = 7).** `y = 2^(v_1 - 1) x_1 - s_1`, and the MH valuation `v(y)` is:

| first symbols of `x` | `v(y)` | MH successor of `y` |
|---|---|---|
| `v_1 = 2, s_2 ≠ s_1` | 2 | `2^(v_2-1) x_2 - s_1` |
| `v_1 = 2, s_2 = s_1, v_2 >= 3` | 3 | `2^(v_2-2) x_2 - s_1` |
| `v_1 = 2, s_2 = s_1, v_2 = 2, s_3 ≠ s_1` | 4 | `(x_2 - s_1)/2` |
| `v_1 = 2, s_1 = s_2 = s_3, v_2 = 2` | 5 (`v_3 = 2`), 6 (`v_3 >= 4`), `>= 7` (`v_3 = 3`) | |
| `v_1 = 3` | 2 | `7 x_1 - 2 s_1` |
| `v_1 = 4, s_2 ≠ s_1` | `3 + v_2` | `x_2` (the orbit rejoins) |
| `v_1 = 4, s_2 = s_1` | 4 | `2^(v_2-1) x_2 - s_1` |
| `v_1 >= 5` | 3 | `2^(v_1-4) 7 x_1 - s_1` |

*Proof.* For q = 7, `s(z) ≡ z (mod 4)`, and `7 x + s_1 = 2^(v_1) x_1`, `7 x_1 + s_2 = 2^(v_2) x_2`, `7 x_2 + s_3 = 2^(v_3) x_3`.
- `y = (7x - s_1)/2 = (2^(v_1) x_1 - 2 s_1)/2`.
- **Case `v_1 >= 3`.** Then `y ≡ -s_1 (mod 4)`, so `s(y) = -s_1`, and
  `7y + s(y) = 2^(v_1 - 1) · 7x_1 - 8 s_1 = 2^(v_1 - 1 + v_2) x_2 - 2^(v_1 - 1) s_2 - 8 s_1`.
  - `v_1 = 3`: this is `4(2^(v_2) x_2 - s_2 - 2 s_1)`, with the bracket odd.
  - `v_1 = 4`: `8(2^(v_2) x_2 - s_2 - s_1)`. This equals `2^(3+v_2) x_2` when `s_2 = -s_1`, and `16(2^(v_2 - 1) x_2 - s_1)`, with the bracket odd, when `s_2 = s_1`.
  - `v_1 >= 5`: `8(2^(v_1 - 4) 7 x_1 - s_1)`, with the bracket odd.
- **Case `v_1 = 2`.** Then `y = 2x_1 - s_1 ≡ 2 - s_1 (mod 4)`, so `s(y) = s_1`, and `7y + s_1 = 2(2^(v_2) x_2 - s_2 - 3 s_1)`.
  - `s_2 = -s_1`: this is `4(2^(v_2 - 1) x_2 - s_1)`, with the bracket odd.
  - `s_2 = s_1`: `8(2^(v_2 - 2) x_2 - s_1)`. The bracket is odd if `v_2 >= 3`. If `v_2 = 2`, it is `x_2 - s_1 ≡ s_3 - s_1 (mod 4)`.
  - If moreover `s_3 = s_1`: `7(x_2 - s_1) = 2^(v_3) x_3 - 8 s_1`, whose valuation is `min(v_3, 3)` for `v_3 ≠ 3` and `>= 4` for `v_3 = 3`. ∎

**Corollary F (where flips pay).** Compare the flipped pair of odd steps `x -> y -> (successor of y)` with the MH pair `x -> x_1 -> x_2`. Both have two odd steps; the number of even steps differs by the **gain** `(1 + v(y)) - (v_1 + v_2)`.
- The gain is `v(y) - 3 >= 1` exactly on the pattern `(s,2),(s,2)`.
- It is `0` exactly on `(s,4),(-s,·)`, where the flipped orbit rejoins `x_2`.
- It is `<= -1` in every other case.
- *Proof:* read off the table.
  - `1 - v_2` for `(s,2),(-s,·)`;
  - `2 - v_2` for `(s,2),(s,>=3)`;
  - `-v_2` for `v_1 = 3`;
  - `1 - v_2` for `(s,4),(s,·)`;
  - `4 - v_1 - v_2` for `v_1 >= 5`. ∎

*Remarks.*
- (a) The gaining pattern is a union of classes mod 32 (it fixes bits 0–4). Those classes are exactly 11 and 21, the classes of `1/3` and `-1/3` (runner D1').
- (b) **At the fixed points.** At `-1/3` (all signs `+`, all `v = 2`) the flip gives `v(y) = 5`: the cycle `-1/3 -> -5/3 -> (v = 5) -> -1/3`, i.e. `(-1,-5,-16,-8,-4,-2)/3`, of density `2/6 = 1/3`.
- (c) **q = 5** (runner B3'). The same computation, with `s(z) ≡ -z (mod 4)`, gives `v(y)`:
  - 2 for `(s,2),(-s,·)`; 3 for `(s,2),(s,>=3)`; 4 for `(s,2),(s,2),(s,·)`; `>= 5` for `(s,2),(s,2),(-s,·)`;
  - `2 + v_2` with rejoining for `(s,3),(-s,·)`; 3 for `(s,3),(s,·)`; 2 for `v_1 >= 4`.
  - At the fixed point `1` (all signs `-`) the flip gives `(1,4)`: `1 -> 3` (flip), then `(5·3 + 1)/2^4 = 1`. This is the sporadic cycle `(1,3,8,4,2)`, of density exactly `2/5`.
  - **So for q = 5 the first renormalization step already lands on the final value `2/5`** (THM-4486 Corollary 5 shows nothing lower is possible). For q = 7 it lands at `1/3` at the fixed point, but the neighbouring patterns `(s,2),(s,2),(-s,·)` gain only 1 (the `(1,4)` pattern, density 2/5). This is where the next stage begins (§5).
- (d) The alternating orbit `±1/11` (signs `+-+-...`) never shows the gaining pattern. Every flip on it loses locally, yet `sigma` must flip somewhere on it (Lemma MH iv). Flipping exactly at `-1/11` gives `(1,2,4,3,5)`, density 1/3; its 2-adic neighbours are what costs.
- *Checked* (runner B3): the table and the gain classification on 100000 random 2-adic integers, exact to `2^-192`, for q = 7; and 100000 for q = 5.

## 4. q = 5: an explicit optimal rule for every level, and why 15

**Construction (runner C1).** A greedy top-down decision tree on the low bits.
- Classes `c mod 2^d` are processed breadth first (MH sign tried first).
- A class becomes a leaf if Min's energy game at `F = 2/5`, with that class forced to one sign, still has a finite least fixed point.
- The engine is `procgen_seven2_20260926_restrict.c`, warm-started. The result is an exploratory object, re-established by the certificates.

**Result (FINITE-EXACT + Lemma C).**
- The tree has **122 leaves** (61 negation pairs) with depths `{3: 2, 5: 2, 6: 6, 8: 10, 9: 12, 10: 14, 11: 24, 12: 12, 13: 12, 14: 20, 15: 8}`.
- 46 leaves flip the max-halving sign.
- Its strategy has `rho_max = 2/5` at k = 15, 16, 17, 18. The witness cycle is re-walked, and the potential certificate is accepted by the seven lane's independent checker.
- By Lemma C, **this single rule has `rho_max = 2/5` at every level `k >= 15`**. It is a level-independent optimal family, and reproves the upper half of Corollary 5 with an explicit object. The greedy finds the same tree at k = 15, 16, 17.

**The rule, explicitly.** `sigma = MH`, except that the sign is flipped on the 46 classes

| modulus | flipped classes |
|---|---|
| `2^6` | 31, 33 |
| `2^8` | 57, 199 |
| `2^9` | 1, 65, 111, 401, 447, 511 |
| `2^10` | 7, 71, 953, 1017 |
| `2^11` | 63, 191, 327, 337, 761, 1287, 1711, 1721, 1857, 1985 |
| `2^12` | 575, 879, 1345, 2751, 3217, 3521 |
| `2^13` | 2735, 3921, 4271, 5457 |
| `2^14` | 1873, 2927, 3409, 5999, 10385, 12975, 13457, 14511 |
| `2^15` | 7023, 10095, 22673, 25745 |

Here MH for q = 5 is `-` on `x ≡ 1` and `+` on `x ≡ 3 (mod 4)`. The flipped classes are pairwise disjoint and closed under negation; the full leaf list is printed in the .out.

**What it does.**
- It is max-halving except on these 46 classes.
- The flips include the `S_inf` fixed points `±1`: class `1 mod 2^9` gets `+` and `511 mod 2^9` gets `-`. This turns `(1,2)` into the sporadic cycle (Corollary F, remark c).
- They also include the period-2 orbit `±1/9` (`57, 199 mod 2^8`) and nested classes around `±1`.
- The sign at `±9/7`, on the density-2/5 cycle `(9,26,13,36,18)/7` of 5x+7, is the MH sign.

**Why 15 (runner C3).**
- The 8 depth-15 leaves form 4 classes mod `2^14` (2 negation pairs). Each class is split into two halves with opposite signs.
- Giving a whole class mod `2^14` one sign creates a cycle denser than 2/5, whichever sign is used (exact witnesses):

| class mod `2^14` | halves (simplest rational, sign) | all `+` | all `-` |
|---|---|---|---|
| 6289 | 6289 (`19/99`, -), 22673 (`-163/13`, +) | 7/17, denominator 17649 | **5/12**, denominator 971 |
| 7023 | 7023 (`-889/1129`, -), 23407 (`9/7`, +) | 3/7, denominator 253 | 7/17, denominator 52947 |
| 9361 | 9361 (`-9/7`, -), 25745 (`889/1129`, +) | 7/17, denominator 52947 | 3/7, denominator 253 |
| 10095 | 10095 (`163/13`, -), 26479 (`-19/99`, +) | **5/12**, denominator 971 | 7/17, denominator 52947 |

- The 5/12 witnesses are the all-outward `(12,5)`-cycles of 5x±971 through `±1651/971`. Their density is `5/12 = rho*(5,14)`, and cycles of this family lie in the tight core of THM-4486's level-14 Max certificate (recomputed here, exploratory).
- **The collision, exactly** (runner C3). In each merged class the two conflict cycles pass through two rationals that are congruent mod `2^14` but not mod `2^15`:
  - class 10095: `4613/971` on the 5/12-cycle, which needs `-` there, and `37501/52947` on a `(17,7)`-cycle of 5x±52947 (`52947 = 2^17 - 5^7`), which needs `+`. Their difference is `2^14 · 12685/(971 · 52947)`.
  - class 7023: `7347/253` on a `(14,6)`-cycle of density 3/7 (`253 | 2^14 - 5^6`), and `44669/52947` on a 7/17-cycle. Their difference is `2^14 · 23053/(253 · 52947)`.
  - The simplest rationals of the two halves are separated at the same bit: `9/7 + 889/1129 = 2^14/7903` and `163/13 + 19/99 = 2^14/1287`.
- **Reading.** With the rest of this rule fixed, a level-14 strategy must give such a merged class one sign, and so keeps one of the two conflict cycles. Level 15 separates them.
- This is a residue collision of exactly the kind Proposition S says every all-level obstruction needs, here resolved at a finite level.
- *Scope.* This explains the depth of this rule. That **every** strategy needs level 15 is `rho*(5,14) = 5/12`, THM-4486's certificate.

## 5. q = 7: the same skeleton, explicit rules, and the structure of both certificates

**The 2/5 rule, explicitly.** `sigma = MH` (for q = 7: `+` on `x ≡ 1`, `-` on `x ≡ 3 (mod 4)`), except that the sign is flipped on the 18 classes
- `11, 21 (mod 2^5)` (the classes of `1/3, -1/3`: the gaining pattern);
- `35, 67, 91, 165, 189, 221 (mod 2^8)`;
- `93, 157, 349, 381, 413, 611, 643, 675, 867, 931 (mod 2^10)` (among them `±1/11 ≡ 931, 93`).

By Lemma C this 18-class rule has `rho_max = 2/5` at every level `k >= 10`.

**Explicit rules (runner D1–D3, FINITE-EXACT + Lemma C).**

| value | level of the certificate | leaves | maximal depth | remarks |
|---|---|---|---|---|
| 2/5 | 10 (re-checked at 10, 12, 14, 16) | 46 | 10 | flips `±1/3` (depth 5 = the gaining pattern, D1'), `±1/11` (depth 10), and `±11/7`, `±7/19`, `±17/3`, `±7/13`, `±23/3`, `±23/35`, `±17/5` |
| 15/38 | 14 | 300 | 14 | |
| 7/18 | 16 | 508 | 16 | 154 leaves at the maximal depth |

- **The residual obstruction of the 2/5 rule (D2).** It consists of the `(5,2)`-cycles of 7x±1 with denominator 17 (`2^5 - 7^2 = -17`), e.g. `(-9,-40,-20,-10,-5)/17`.
  - Such a cycle is a flip at a point of the class of `1/3` with `s_3 ≠ s_1`, followed by valuation 4: the gain-1 case of Corollary F.
  - Going below 2/5 means giving up these gain-1 flips. The alternative flips cost more, and the rule grows from 46 to 300 leaves.
  - The residual obstruction of the 15/38 rule (exploratory, k = 14) is already a family of 68 `(38,15)`-cycles with denominator `7^15 - 2^38 = 4472683602999`. Each contains 4 flips, e.g. with valuations `(4,2,2,2,4,2,2,3,3,1,4,1,1,1,6)`: no longer a single local pattern.
- **Growth.** The description length grows fast as the value decreases: 46 → 300 → 508 leaves for `2/5 → 15/38 → 7/18`. At the maximal depth, 154 of the 508 leaves are single residues.
  - The greedy trees are not canonical: they depend on the order, and the shallow leaves of the 15/38 tree differ from those of the 2/5 tree.
  - So the growth rate is only indicative.

**The certified optima (runner D4–D6; the seven lane's engine re-solved at the certified values; EMPIRICAL statistics of certified objects).**
- *Min.*
  - The least Min potential admits the MH sign at 84.8–87.1% of the odd residues at k = 19..26 (runner D4), and at 87.9% at k = 27 (exploratory). The exception rate decreases slowly: 0.152, 0.144, 0.138, 0.145, 0.142, 0.130, 0.130, 0.129, 0.121 for k = 19..27.
  - The exceptions:

| x mod 8 | 1 | 3 | 5 | 7 |
|---|---|---|---|---|
| exception rate within the class (k = 19) | 0.030 | 0.274 | 0.274 | 0.030 |
| (k = 24) | 0.022 | 0.239 | 0.239 | 0.022 |

  - MH has valuation 2 exactly on `x ≡ 3, 5 (mod 8)`. On the classes of `±1/3` mod 32 the exception rate is 0.55–0.61.
  - By the length `r` of the leading run of valuation-2 steps (k = 24, exploratory): `P(exception | r) = 0.022, 0.078, 0.31, 0.47, 0.48, 0.51, 0.57, 0.62, 0.63` for r = 0..8.
  - On `S_inf` cylinders, the forced-flip share by the first three signs (k = 24) is `+++`: 0.67, `++-`: 0.60, `+-+`: 0.49, `+--`: 0.22 (and negatives). This is the bias Corollary F predicts (equal leading signs gain), but it is not a function of three signs.
- *Nesting* (D5). Between consecutive levels with the same value (k = 22→23, 24→25, 25→26), the forced-MH / forced-flip / free classification agrees with the lift of the previous level on 94.0%, 95.4%, 95.6% of the residues. Direct conflicts are 0.24%, 0.09%, 0.07%.
- *Max* (D6). The share of pairs whose lift is forced to the "near" representative (the one of smaller absolute value) is the same, to within 0.01, in each of 8 height bins, at every tested level.
  - The optimal adversaries carry no information about real height. This is why height potentials (THM-4486 §7) cannot certify them.
  - The Max potential grows linearly in `v_2(P)`, with slope about `0.78 fn` (k = 22–26, exploratory): the forced halvings.
- *Critical cycles* (exploratory; the tight cores of both certificates at k = 20, 22, 24).
  - Almost all contain flips (valuation 1) and valuations up to 8–10.
  - Only the `(8,3)` cycle of denominator 87 at 3/8 and a handful of `(24,9)`, `(32,12)`, `(40,15)` cycles are pure max-halving.
  - All are real-bounded rationals (`|x| < 1100`) with large denominators.

**Renormalization and Fibonacci.** The stages are:
- stage 0: MH, value `1/2`, obstruction `S_inf`;
- stage 1: flip the gaining pattern and the alternating orbit, value `2/5`, obstruction: the denominator-17 cycles of type (flip, 4);
- later stages: refinements, nested at 94–96%, whose obstructions are long, flip-containing, real-bounded cycles with huge denominators.

In mean-valuation terms, `1/rho*(7,k)` is:

| k | 10 | 14 | 16 | 18 | 19 | 20 | 21 | 22 | 24 | 27 | 28 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `1/rho*` | 2.500 | 2.533 | 2.571 | 2.579 | 2.615 | 2.618 | 2.625 | 2.643 | 2.667 | 2.686 | 2.692 | 2.703 |

- The Fibonacci values `13/34, 34/89, 8/21` (k = 19–21) are the passage through `phi^2 = 2.618`. Nothing stops there: `1/rho*(7,30) = 100/37 = 2.7027`.
- No Sturmian or Fibonacci structure was found in the critical cycles: they use valuations 1..10, not two consecutive ones.
- The threshold for provability is `log_2 7 = 2.807`. Least-squares fits over k = 10..30 (runner H) give limits on both sides of the threshold: `A - B/k` gives `2.792` and `A - B/sqrt(k)` gives `2.994`. They do not decide anything.
- The same holds for the decrements. From `3/8` (k = 24) to `37/100` (k = 30) the value falls by 0.00083 per level on average, and the gap to `log_7 2` at k = 30 is 0.01379.
  - Constant decrements would cross the threshold near k ≈ 47.
  - Decrements decaying like `1/k` would cross near k ≈ 55.
  - Geometrically decaying decrements with ratio below about 0.94 per level would never cross.
  - The data cannot distinguish these (EMPIRICAL, no claim).

## 6. Automatic strategies (route 2) and adversaries (route 3)

**Real-sign imitation fails completely (E1, E2; FINITE-EXACT at each listed level).**
- Take `sigma(x)` = the real sign of the simplest rational in the class of `x`: the shortest vector `(a, D)`, `D` odd, of the lattice `a ≡ D x (mod 2^k)`.
- Then `rho_max = 1` for q = 5 and 7 at k = 10, 12, ..., 20, with all-odd cycles of real value near `±1/(q-2)`.
- Using that sign only on classes whose simplest rational has height `<= h` (h = 3, 9, 33, 129), and MH elsewhere, never beats MH.
- *Why.* Corollary F: a flip away from the gaining pattern loses halvings. The real-sign rule flips wholesale.
- Proposition S (THM-4486 update) concerns exact rationals, and a residue class does not know which rational it represents.

**Small automata (E3).** Consider flip sets recognised by a DFA that reads bits 1, 2, ... of `x`, with absorbing FLIP/KEEP states and a default for unresolved transient states. None beats MH (`1/2`) at k = 12:
- all DFAs with 2 transient states (98 distinct strategies);
- all with 3 (1970);
- 20000 random with 4 (5514 distinct).

The 2/5 rule already needs depth 10 and 23 negation pairs of leaves.

**The "finite product graph" (discussion).**
- By Lemma C, a rule that reads `d` digits is exactly its level-`d` graph. Composing it with the multiplication transducer gives nothing smaller in general.
- A rule reading unboundedly many digits (an infinite regular flip tree) defines `T_sigma` off a closed null set. Its truncations at level `k` are ordinary level-`k` strategies, and their `rho_max` is what the game computes.
- The finite-state content that does exist is Lemma F. The effect of one flip is a finite function of three MH itinerary symbols, and a flip rule that is a union of MH-itinerary cylinders acts on the full shift of itineraries by finitely many local rewrites.

**Adversaries (F1; exploratory hill-climb).**
- `tau(P) = f(P mod 2^m)` for all `f`:
  - best value `1/3` for m = 1, 2, 3 at k = 10 and m = 1, 2 at k = 12;
  - at k = 8 some m = 3 rules reach 5/13, a small-level effect: the best one falls to 2/7, 17/65, 1/4 at k = 10, 12, 14.
- A hill-climb over m = 6 at k = 12 (exploratory) did not exceed 1/3. Its local optimum from the top lift had 5/18 and 19/69 at k = 14, 16.
- *Reading.* The certified optimal adversaries at k = 19..26 are height-independent (D6) and not low-bit (their forced lifts are balanced in every class mod 8). They are 2-adically global objects, as Proposition S requires of any all-level obstruction.

## 7. Levels 30 and 31 (route 4)

**Engines.**
- `procgen_seven2_20260926_lean8.c` is the seven lane's one-sided fixed-point computation with one byte per negation representative (`2^(k-2)` bytes plus a dirty bitset). It has an optional offset for the Max potential: `h = g - fn·v(P)` is stored, where `v(P)` is the number of forced halvings (`v_2(P)`, and `k-1` for `P = 0`).
  - The offset is sound because `g >= fn·v(P)` at the least fixed point: each forced halving costs `fn`.
  - The iteration from `g_0 = fn·v` still converges to the least fixed point, since `g_0 <= Phi(g_0)` and `g_0 <= g*`.
  - Values of `h` above the cap become TOP, which can only shrink `W`.
- `procgen_seven2_20260926_verify8.c` is a checker written separately from the engines.
  - It first verifies that the potential file is invariant under `P -> H - P`, streaming the second half against the first.
  - It then checks every inequality for all `H` pairs (lower) or all `N` nodes (upper), reading potentials from the kept half and bits in streaming order.
  - It needs about `2^(k-2)` bytes per potential byte.
  - It accepts the seven lane's `lo_g16/up_psi16` files and the new `lo_h8`.
  - It rejects mirror-asymmetric files and a raised entry at a tight edge (runner A3).

**Results** (FINITE-EXACT; runner §G; times and RSS of the final run, under the load of other lanes):

| level | bound | certificate | engine time / RSS | checker |
|---|---|---|---|---|
| 30 | `rho*(7,30) <= 10/27 = 0.37037` | upper, one byte, max `psi` 76 | 72 s, 293 MiB | verify8 (261 MiB) and the seven lane's checker (577 MiB) |
| 30 | `rho*(7,30) >= 37/100` | lower, the seven engine (16-bit), max `g` 1073 = 37·29 | 64 s, 548 MiB | verify8 (16-bit half, 521 MiB) |
| 30 | `rho*(7,30) <= 37/100` | upper, the seven engine (16-bit), max `psi` 372; 1511 sweeps | 690 s, 548 MiB | verify8 (521 MiB) |
| 31 | `rho*(7,31) >= 7/19 = 0.36842` | lower, one byte with the offset, max `g` 210 = 7·30 | 71 s, 581 MiB | verify8 (517 MiB) |

- Hence **`rho*(7,30) = 37/100`** (`= [0;2,1,2,2,1,3]`, mean valuation `100/37`), and `7/19 <= rho*(7,31) <= 37/100`.
- **No provable 7n±1 strategy at any level `k <= 31`** (`7^37 > 2^100`, `7^7 > 2^19`).
- The value drops again at `k = 30` (`37/100 < 13/35 = rho*(7,29)`), by `0.00143`.
- The search. Farey bracketing with both one-sided tests run concurrently at each mediant:
  - `upper 10/27` won;
  - `lower 27/73` won (exploratory, subsumed);
  - `lower 37/100` won;
  - `upper 37/100` won.
  Several losing tests were stopped without a verdict (§8).

**Limits.**
- At k = 30 the offset `h` for `F = 27/73` still reaches 710: 6.5% of the pairs exceed 254. So the lower certificates near the value need 16-bit potentials, i.e. `2^(k-1)` bytes on the negation quotient: 1 GiB at k = 31.
- An exact value at `k >= 31` therefore does not fit in 700 MB with these methods. Only one-sided bounds with small potentials do, for example `7/19` above.
- The Min side of `F = 37/100` also needs 16 bits (`psi` up to 372).

## 8. Failures and caught mistakes

- **A too-coarse search objective.** The first attempt to find simple strategies was a greedy ball-flipping on `rho_max`. It stalled: `rho_max` stays at 1/2 while any `S_inf` orbit survives. The "number of pairs reaching a dense cycle" objective was also useless, since all pairs reach one. The restricted-game tree builder replaced both.
- **Non-canonical trees.** At k = 13 the greedy returned a 59-leaf tree for 2/5, although the 46-leaf tree of k = 10 is valid there (Lemma C). The greedy order, not a mathematical effect, explains the difference. So tree sizes are reported as upper bounds on the description length, not as intrinsic quantities.
- **An over-strong adversary claim, caught by the runner.** The first version of F1 claimed "best value 1/3 at k = 8, 10". At k = 8 some low-bit rules reach 5/13. The check now states the small-level effect and follows the rule to k = 14.
- **The k = 30 and k = 31 searches.** Losing one-sided tests are slow, so several were stopped without a verdict:
  - at k = 30: `lower 13/35` (30 min of CPU, niced); `lower 10/27` and `upper 17/46` (about 7 min each); `upper 27/73` (15 min); `lower 47/127` (5 min);
  - at k = 31: `upper 7/19` (32 min) and `lower 17/46` with one-byte potentials (10 min). The latter may also have failed only because of the one-byte cap.
  - The protocol that worked: run both sides at the same Farey mediant concurrently and take the first to finish. A winning upper test can itself be slow: `upper 37/100` at k = 30 needed 1511 sweeps (14 min).
  - So `rho*(7,31)` is only bracketed: `[7/19, 37/100]`.
- **Memory.**
  - The seven checker needs `2^(k-1)` bytes of potential and does not fit at k = 30 with 16-bit potentials. This motivated the mirror-streaming checker; the seven checker is still used where it fits (the k = 30 upper certificate at 10/27, 577 MiB).
  - The first complete run of the runner passed every check but its own process peaked at 702 MiB, above the lane cap. The cause was numpy temporaries in the potential statistics (k = 26) and in the lattice reduction (k = 20). Both were chunked, and the run was repeated; the .out is from the repeated run.
- **Numerology recorded, not claimed.** `1/rho*(7,30) = 100/37 = 2.7027` is close to `e = 2.718`, and `phi^2` was passed at k = 20. No claim is made.
- **Nothing in the canon is contradicted.** THM-4486's values and the seven lane's values are reproduced wherever recomputed (k = 19..27 in D4, 24 in A3). Its OPEN question stays OPEN.

## 9. Reproduction

```bash
/usr/bin/time -l python3 -u 04-computation/experiments/procgen_seven2_20260926_run.py > 05-knowledge/results/procgen_seven2_20260926.out
```

- **Requirements.** `numpy` and a C compiler. The floor lane's `procgen_floor_20260926_lib.py` builds its engine into `scratch/procgen_floor/`. All other binaries go to `scratch/procgen_seven2/run/`, with source hashes in their names; certificates are deleted after checking.
- **Environment.**
  - `SEVEN2_TREE7_KMAX` (16): the q = 7 tree series;
  - `SEVEN2_K30` (1): levels 30 and 31;
  - `SEVEN2_STRUCT_KMAX` (26): the potential statistics.
- **Output.** 28 checks, ending with `ALL CHECKS PASSED`.
- **Cost.**
  - Wall time 1411 s (23.5 min) on the shared 8-core machine, niced, with other lanes running. About 11 min of it is the level-30 upper certificate at 37/100 (1511 sweeps).
  - Peak child RSS is 581 MiB: `lean8 lower 7 31 7 19`. `/usr/bin/time -l` of the whole run reports a maximum resident set of 609,140,736 bytes.
  - The runner process itself peaks at 463 MiB, from the potential statistics at k = 26.
  - No process exceeds 700 MB.
- **SHA-256** (raw bytes; the .out prints the script hashes itself):
  - `procgen_seven2_20260926_rhomax.c` `ea89e7a8aefe13d0f0ed9824cb3d3297e6c2b2367710d94d0ec94b17876a673b`
  - `procgen_seven2_20260926_restrict.c` `305aa42892d0758249c7bcb3228ed83ffaff76bae4289604a23bfd0000e672aa`
  - `procgen_seven2_20260926_lean8.c` `a4dc2b4b1bb3ab4afbd6a9a0798962f2fc2d5532651cb692c69e4b168d5d1801`
  - `procgen_seven2_20260926_verify8.c` `51c97489cadc4c98c7eae7b533e690d1c4b102e8e518fe781d29d727618c5da4`
  - `procgen_seven2_20260926_lib.py` `e59943587822f54d0d201788016a46d1b7d0ee89b443a83c7173103780e96c71`
  - `procgen_seven2_20260926_run.py` `fa2a28d46cf6a3c2b8060b8849180f5ed584b847beabc1d6db1715ea0093a653`
  - `procgen_seven2_20260926.out` `a60b79a44a4e2189625d1532eb43b432f8916efc36a7199a388a005a39f349db`
  - reused read-only, unchanged:
    - `procgen_seven_20260926_game.c` `4d2a193f96126ca5c13e4377c9d9018da7c25034cf0feda1152198750e5d7b34`
    - `procgen_seven_20260926_verify.c` `53576fbcf5472a989055623af881defafe9c16b54df28fb2b836989c54cdef39`
    - `procgen_seven_20260926_tight.c` `1ff6968c4d7509d079259f4f7eb771309dd9078df0f93ed07d439f1fda9efb14`
    - `procgen_floor_20260926_lib.py` `be7024fa7dc147c21d3eacabfc59d09c558943dcfa3cb30b421abda79e3e4f2f`
- **Timing and RSS lines** in the .out vary between runs. Values, witnesses, tree shapes and all checks are deterministic.
