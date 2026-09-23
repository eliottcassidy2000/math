# The E-game exceptional set looks countable: exact threads to `2^64` and `3^41`, a perturbation lemma, and a certified census of hostile rationals

**Status: `dim Bad_inf(E)` is still OPEN (neither direction is proved), but
the evidence now favours dimension `0`: a countable set of finite
Cantor–Bendixson rank, with polynomial covering numbers.** The rest is typed
as follows.

* **PROVED:**
  * the exact integer reformulation and the Beatty doubling law (§1.1);
  * the safety lemma S and the carry identity C (§4.1);
  * the perturbation lemma, Theorem P (§4.2);
  * Theorem F (§4.3): `x_i = -1 - 2^i/3^alpha(i)` is hostile for every
    `i >= 3`, so `Bad_inf(E)` is infinite and `-1` is an accumulation
    point;
  * the backward lemmas S_b and C_b (§4.7).
* **FINITE-EXACT** (where marked, independent code paths agree):
  * forward counts to `m=64` and backward counts to `r=41` (§2);
  * `alpha(i) = ceil((i+1) log_3 2)` for `9 <= i <= 56`;
  * the closure of `{-1}` under Theorem P has exactly four generations
    `0..3` up to 2-adic scale `56`;
  * a complete census of the hostile points `-p/3^j` with
    `1 <= |x| < 3/2`: 343 points for `j <= 25` and 382 for `j <= 27`,
    each with an exact certificate. They form a tree of canonical
    perturbations whose depth never exceeds 7;
  * every class of the projected exceptional set at levels `20` and `30`
    contains a `Z[1/3]` point;
  * Q2 verified below `7.87*10^17`.
* **HEURISTIC / CONDITIONAL:** countability (§4.6), and Conjecture G:
  a finite seed generates every hostile rational through Theorem P
  (verified for `15 <= j <= 27`). The "credit construction" fails at an
  identified step (§4.6).

Session `collatz-procgen-20260922` (mac-mini), dimension lane. Scripts:
`04-computation/experiments/collatz_procgen_20260922_dim_*`. Output:
`05-knowledge/results/collatz_procgen_20260922_dim.out`. Parent note:
[choice ladder](collatz_procgen_20260922_choice_ladder.md), §§1, 2, 4, 5.

## 0. The answer in brief

1. **Algorithm.** Cost decreases only at halvings, so a class `x mod 2^m`
   is good iff `f_b(x) <= floor(b log_3 2)` for some `b <= m`. Here `f_b`
   is the least number of `x3`-moves over E-paths making exactly `b`
   halvings. `f_b` is integer valued, exact, and nondecreasing in `b`. So
   the thread computation is a small integer DFS:
   * the DFS stops at an exact oracle table `f_30` (512 MB);
   * it prunes with the bound `f_p(y) >= f_30(y mod 2^30)`;
   * it skips every level where `floor(m log_3 2)` does not grow. There
     the count provably doubles.

   The old DFS needed about `7*10^10` nodes to reach `m=40`. The new one
   needs `7.7*10^6` nodes at `m = 40` itself (`P0 = 30`); a from-scratch
   run to `m = 40` with `P0 = 26` takes about one second. Forward counts
   now reach `m=64` (`1.34*10^12` nodes, 93 min for that level alone),
   and backward counts reach `r=41`. Every overlapping count matches
   the old programs exactly.
2. **Growth.** Raw counts are misleading: they double at every
   non-Beatty level and contain slow descenders.
   * The covering number `K_m` of `Bad_inf` lies between two curves:
     * **below,** a certified bound (classes containing a proved hostile
       point) that is essentially linear, `-88 + 10.9 m` on
       `m = 20..40`;
     * **above,** the converged projected count `P_m = |Bad_61 mod 2^m|`,
       about `13.5` per level. Its best fit is `1.18 m^1.60` (log-RMS
       `0.044`, against `0.100` exponential).
   * Over this window polynomial and `2^(0.06 m)` growth differ by only
     about 20%, so the fits lean polynomial without deciding.
   * The Horton–Strahler number of the splitting tree is constant, `5`
     for `L = 16..58`. The backward side behaves the same way:
     `~r^1.7`, Strahler number `3` for `L = 6..38`.
3. **Theory.** Every hostile point we can identify is a negative element
   of `Z[1/3]` with real value in `[-3/2, -1]`.
   * Lemma C turns the "credit" into a **real-valued budget**: at an
     integer state the branch is safe iff the carry ratio `beta` exceeds
     `delta = |x| - 1`, and `beta <= 1/2 + O(3^-a)` on the cheapest
     branches.
   * Theorem P builds hostile points by perturbation
     `h -> h - 2^i c/3^j`. Each perturbation spends part of the budget
     `1/2`. Its closure dies out after three generations.
   * The certified census shows about 17 hostile points per 3-adic
     denominator exponent `j`, with no growth.
   * The census is a tree rooted at `-1` of **bounded depth**: at most 7
     overall, and at most 6 for every `j = 15..25`. Each step is a
     perturbation `h -> h ± 2^i/3^(alpha_h(i)+e)` with `e in {-1,0,1}`,
     and 91% are the canonical step `e = 0`, downward.
   * So the credit construction does not concatenate: the credit banked
     near `-1` is the slack `1/2 - delta`, and every block spends a
     definite part of it.

## 1. Exact reformulation and the fast algorithm

### 1.1 Value functions (PROVED)

**Forward.** For `y in Z_2` and `b >= 0`, let `f_b(y)` be the least
number of `x3`-moves over E-paths from `y` that make exactly `b` halvings
(the path stops right after its `b`-th halving). The legality of such a
path is decided by `y mod 2^b`. Every prefix of a path making `b`
halvings is a path making `b-1` halvings, so
`f_b(x mod 2^b) >= f_(b-1)(x mod 2^(b-1))` (**monotonicity**).

The running multiplier drops below `1` only at a halving, and
`3^a < 2^b` iff `a <= A_b := floor(b log_3 2)`. Hence:

* `x mod 2^m` is exceptional iff `f_b(x mod 2^b) >= A_b + 1` for all
  `b <= m`;
* `Bad_m = { x : x mod 2^(m-1) in Bad_(m-1)` and `f_m(x) >= A_m + 1 }`.

**Beatty doubling law (PROVED).** If `A_m = A_(m-1)`, every lift of a
class in `Bad_(m-1)` has `f_m >= f_(m-1) >= A_m + 1`. So
`|Bad_m| = 2|Bad_(m-1)|` exactly, which is visible in every row of the
§2 table. Only at the `63%` of levels where `A_m` grows can classes die,
and then only the **tight** lifts, those with `f_m = A_m`.

The recursion is exact and integer:

* `f_0 = 0`;
* for odd `y`: `f_p(y) = 1 + f_p(3y+1)`;
* for even `y`: `f_p(y) = min_j [2j + f_(p-1)(y_j/2)]`, where `y_0 = y`
  and `y_(j+1) = 9y_j + 4 (mod 2^p)`.

**Backward.** `g_s(x)` is the least total `K` over legal `s`-move reverse
paths, decided by `x mod 3^(s+1)`. It is monotone in `s` for the same
reason, and `x` is exceptional at level `r` iff `g_s >= floor(s log_2 3)
+ 1` for all `s <= r`. The recursion is
`g_s(x) = min_(legal k) [k + g_(s-1)(y_k)]`, with
`y_(k+2) = 4y_k + 1`. Since `floor(s log_2 3)` grows at every level,
there is no doubling law on this side.

### 1.2 The thread algorithm

* **Phase 1.** Build the full DP table of `f_p` for `p <= P0`, storing
  even entries only as `uint8` (`P0 = 30`: 512 MB). Keep `Bad_p` as
  threads.
* **Phase 2.** At a Beatty level `m`, test each lift by a DFS over its
  top `m - P0` bits:
  * the tail is exact: a state at precision `P0` contributes
    `f_P0(y)`;
  * the pruning bound is `f_p(y) >= f_P0(y mod 2^P0)`, by the prefix
    argument.
* **Both lifts at once.** The two lifts `x`, `x + 2^(m-1)` of a class
  have identical parities and identical pruning lookups at every
  internal node, since their states differ only in the top bit. At a
  leaf they differ exactly in bit `P0 - 1`. One search therefore
  decides both lifts. This is the default (`dfs2`); `SINGLE_LIFT=1`
  tests them separately.
* **Backward.** Same structure, with a table `g_17` modulo `3^18`
  (258 MB) and `__int128` arithmetic for `r >= 40`.
* **Transposition-table variant.** `..._dim_forward_tt.c` keeps exact
  IDA*-style lower and upper bounds on `f_p(y)` in a hash table that
  persists across levels. It cuts nodes 7- to 13-fold (13-fold at
  `m = 48` with the table at every precision, 7-fold at `m = 50` with it
  only at `p >= P0 + 4`). Memory latency eats the gain, so it serves as
  an independent check only. At `m=43`
  the plain DFS visits `5.3*10^8` nodes but only `2.8*10^7` distinct
  states: hostile points map to hostile points.

### 1.3 Validation (FINITE-EXACT, independent paths)

* The new integer DP equals the old float DP (`e_forward_dp_full`) for
  all `m <= 24`, and the old reverse DP for all `r <= 14`.
* Combined-lift and single-lift searches give identical exceptional
  lists at `m = 27, 32, 35, 40, 43, 44, 45, 48, 50`.
* The TT variant, with separate search code, gives identical counts to
  `m = 50`.
* Every count quoted in the choice-ladder note is reproduced exactly:
  * forward: `124, 391, 369, 255, 561, 454, 908, 1030` at
    `m = 16, 21, 24, 27, 32, 35, 36, 40`;
  * backward: `14, 18, 32, 67, 94, 134, 190, 157` at
    `r = 6, 8, 10, 19, 24, 27, 30, 31`, with minimal representatives
    `7.17*10^12` and `2.02*10^13` at `r = 30, 31`.

### 1.4 Cost

| run | memory | time |
|---|---|---|
| forward, `P0 = 26`, `m <= 50` | 0.1 GB | 30 s (old DFS: about 20 min for `m <= 40`) |
| forward, `P0 = 30`, `m <= 64` | 0.8 GB | levels `59`, `61`, `62`, `64`: `1.3*10^11`, `4.4*10^11`, `4.5*10^11`, `1.34*10^12` nodes (11, 31, 33, 93 min); 3.0 h in total |
| backward, `P0 = 17`, `r <= 41` | 0.4 GB | `9.5*10^10` and `2.1*10^11` nodes at `r = 40, 41` (19 and 42 min) |

Node counts grow about `1.5`-fold per forward Beatty level and `2.1`-fold
per backward level, so `m ~ 80` stays out of reach. The theory below
makes it unnecessary.

## 2. Counts (FINITE-EXACT)

**Forward, `|Bad_m|` (classes mod `2^m`).** `A_m = floor(m log_3 2)`.
A starred level is non-Beatty, where the count doubles exactly.
`minrep` is the smallest nonnegative representative of an exceptional
class other than `-1`. `minabs` is the least `|n|` over integers `n`, of
either sign, in such a class.

| m | 30 | 31 | 32 | 33* | 34 | 35 | 36* | 37 | 38* | 39 | 40 | 41* | 42 | 43 | 44* | 45 | 46 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| count | 1052 | 1197 | 561 | 1122 | 1195 | 454 | 908 | 936 | 1872 | 2152 | 1030 | 2060 | 2194 | 809 | 1618 | 1667 | 518 |

| m | 47* | 48 | 49* | 50 | 51 | 52* | 53 | 54 | 55* | 56 | 57* | 58 | 59 | 60* | 61 | 62 | 63* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| count | 1036 | 1054 | 2108 | 2344 | 988 | 1976 | 2073 | 758 | 1516 | 1548 | 3096 | 3523 | 1601 | 3202 | 3379 | 1198 | 2396 |

`m = 64` (Beatty): 2454. `m <= 29` as in the parent note.

* **Local minima:** `255, 454, 1030, 809, 518, 988, 758, 1601, 1198` at
  `m = 27, 35, 40, 43, 46, 51, 54, 59, 62`. They are not monotone,
  because the raw count at a level depends on `frac(m log_3 2)`.
* **`minrep`:** `1.50*10^8` (`m = 35`), `3.70*10^9` (`40`),
  `1.43*10^11` (`46`), `3.91*10^12` (`51`), `6.12*10^13` (`56`),
  `3.19*10^14` (`61`–`63`), `2.20*10^16` (`64`). So every positive
  integer below `2.20*10^16` lies in a class with a multiplicative
  certificate of precision 64. `minabs` at `m = 64` is
  `2.93*10^14`, from a negative integer.
* **Correction to the parent note.** Its sentence "the only exceptional
  class mod `2^36` containing an integer of absolute value below
  `1.5*10^8` is `-1`" is false for negative integers. `minabs = 17797505`
  at `m = 35..38`: the class of `-17797505` first gets a certificate at
  precision `37`. Exactly, `-17797505` descends, with `3^46 < 2^73`. The
  positive-integer statement is correct.

**Backward, `|Bad_r|` (classes mod `3^(r+1)`)** and the smallest positive
representative of an exceptional class other than `1`.

| r | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| count | 190 | 157 | 181 | 183 | 233 | 235 | 198 | 212 | 214 | 258 | 260 | 338 |
| minrep | 7.17e12 | 2.02e13 | 2.02e13 | 2.02e13 | 2.02e13 | 6.32e14 | 6.32e14 | 1.95e15 | 5.96e16 | 5.96e16 | 9.53e16 | 7.87e17 |

(`r = 41`: `2.1*10^11` nodes, 42 min. The exact `minrep` at `r = 41` is
`786759327538998251`.) The parent note's argument applies
unchanged: backward certificates give actual descent to a positive
integer prime to `3`, and strong induction finishes. **Hence `1` reaches
every `m < 7.87*10^17` with `3` not dividing `m` in `E`**, up from
`2.02*10^13`. By the same step, Le–Smith's Conjecture 1 holds below
this bound.

## 3. Growth exponent: fits, and what the data support

**Why raw counts are the wrong statistic.** First, they oscillate with
the Beatty word of `log_3 2`: exact doubling at non-Beatty levels,
deaths of tight classes at Beatty levels. Second, they contain **slow
descenders** that stay exceptional for many levels without being
hostile (§4.5). For example, `-9905/6561` first descends after 233
halvings.

**Projected counts.** `P_m(M) = |Bad_M mod 2^m|` counts the classes that
still have an exceptional descendant at level `M`. It is an upper bound
for the covering number `K_m` of `Bad_inf` at scale `2^-m`, and is exact
as `M -> inf`. It agrees for `M = 56` and `M = 61`, and for `M = 59`
and `M = 63`, at every `m <= 54`, so it has converged there.

| m | 16 | 20 | 24 | 28 | 32 | 36 | 40 | 44 | 48 | 52 | 54 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `P_m` (upper bound for `K_m`) | 91 | 147 | 190 | 250 | 300 | 354 | 416 | 468 | 571 | 659 | 758 |
| certified lower bound: classes containing a census point (§4.5) | 83 | 132 | 170 | 215 | 257 | 300 | 343 | 343* | – | – | – |
| weaker proved lower bound: Theorem-P closure (§4.4) | 55 | 77 | 96 | 119 | 139 | 163 | 191 | 213 | 249 | 269 | – |

\* The census stops at `j = 25`, so the census bound saturates beyond
`m ≈ 1.585*25 - 3 ≈ 37`.

* **Where the gap sits.** Between the two bounds lie the classes whose
  only `Z[1/3]` point is a slow-descender candidate: 15, 37 and 54
  classes at `m = 20, 30, 36`.
* **Why the upper curve is inflated.** These classes persist beyond
  `M = 61`, so `P_m` overstates `K_m` by a growing amount.

**Fits, log residuals.**

| data | window | power law | exponential | linear |
|---|---|---|---|---|
| upper bound `P_m` | `m = 16..54` | `1.18 m^1.600`, RMS `0.044` | `56 * 2^(0.071 m)`, RMS `0.100` | RMS `0.123` |
| certified lower bound | `m = 20..40` | `m^1.41`, RMS `0.024` | `2^(0.070 m)`, RMS `0.046` | `-88 + 10.9 m`, RMS `0.022` |
| Theorem-P closure | `m = 16..56` | `m^1.36`, RMS `0.021` | RMS `0.083` | – |
| backward `Q_r` (converged) | `r = 8..37` | `r^1.68`, RMS `0.059` | `2^(0.121 r)` (`0.076` trit/level), RMS `0.128` | – |

**Local exponential rates** (windows of 10 levels).

* **Upper curve `P_m`,** at `m = 16, 18, ..., 44`: `0.136, 0.098, 0.088,
  0.088, 0.084, 0.060, 0.062, 0.062, 0.060, 0.046, 0.055, 0.058, 0.061,
  0.054, 0.070` bit/level. It falls until `m ~ 34` and then levels off
  near `0.055–0.07`. Slow descenders feed that plateau.
* **Certified lower bound,** at `m = 16, ..., 30`: `0.129, 0.089, 0.082,
  0.082, 0.076, 0.056, 0.061, 0.056`, still falling.
* **Backward** (windows of 8): from `0.18` to about `0.09` by `r = 17`, then fluctuating in `0.09–0.13` for `r = 18..23`.

Over `m = 20..40` a linear law and `2^(0.06 m)` differ by only about
20% at the ends. **So the fits favour polynomial growth, with
linear-plus-intercept best, but they are not decisive.**

**Tree shape.** The Horton–Strahler number of the compressed splitting
tree of the projected set is `4, 5, 5, 5, 5, 5, 5, 5, 5` at
`L = 10, 16, ..., 58`. The number of branch nodes of orders `4` and `5`
has been frozen at 7 and 7 since `L = 18`, while the leaves grow from
127 to 637. The backward Strahler number is `3` at every
`L = 6..38`.

For comparison, a branching Cantor set of dimension `d` gains between
about `dL/2` orders (random-like branching, Strahler about `log_4` of
the number of leaves) and `dL` orders (regular branching). With
`d = 0.06` that is 1 to 2 orders over `L = 16..52`, where the leaves grow
7-fold (`log_4 7 ~ 1.4`). Instead we see no gain at all. A countable
comb-of-combs, by contrast, keeps a fixed Strahler number.

**What the data support.**

* **Support:**
  * covering numbers squeezed between a certified, essentially linear
    lower bound (about `11` new classes per level) and an upper bound
    about `13.5` per level on `m = 20..40`;
  * a splitting tree of fixed Strahler number;
  * a census of hostile points that is bounded per 3-adic exponent and
    of bounded tree depth (§4.5).

  All of these point to a countable set of finite Cantor–Bendixson rank,
  so dimension `0`. The backward side agrees.
* **Do not support:** a dimension near `0.1`. That value came from local
  minima of raw counts, which are Beatty artefacts.
* **Do not decide:** counts alone cannot exclude a dimension below about
  `0.07` over this window. Nor do they exclude a Cantor set that starts
  branching only beyond `m ~ 60`. The structural facts of §4.4–4.5 are
  the stronger evidence.

## 4. Theory

Throughout, for `x in Q` with odd denominator (so `x in Z_2`) and an
E-path with `a` multiplications, `b` halvings and carry `B`, the state is
`v = (3^a x + B)/2^b`. Write `R = 3^a/2^b`, `beta = B/3^a`, and
`B = sum_(s<a) 3^(a-1-s) 2^(b(s))`, where `b(s)` is the number of
halvings before the `s`-th multiplication. `beta` is nondecreasing along
a path: a halving leaves it fixed, and a `x3`-move adds `2^b/3^(a+1)`.
Also `beta >= (1 - 3^-a)/2`, with equality iff all multiplications
precede the first halving.

### 4.1 Two lemmas (PROVED)

**Lemma S (safety at negative integers).**

1. The negative integers are closed under both E-moves: `3v+1 <= -2` for
   `v <= -1`, and `v/2 <= -1` for even `v <= -2`.
2. Suppose a path reaches an integer `v <= -1` with multiplier `R`. Then
   every continuation has running multiplier at least `R/|v|` at each of
   its halvings, and at least `R/(|v| - 1/3)` when `v` is odd. So if
   `R > |v|` (or `v` is odd and `R > |v| - 1/3`), **no continuation ever
   descends**.

*Proof.* A continuation value `(3^a' v + B')/2^b'` is an integer
`<= -1`, so `3^a'|v| >= 2^b' + B'`. Hence `R' >= (1 + B'/2^b')/|v|
>= 1/|v|`. If `v` is odd its first move is `x3`, so `B' >= 3^(a'-1)`,
which gives `R' >= (1 + R'/3)/|v|`. ∎

**Lemma C (carry identity).** If `x < 0` and the state `v` is negative,
then `|v| = R(|x| - beta)`. So at an integer state `v <= -1`,
`R > |v|` iff `beta > delta := |x| - 1`.

*Proof.* `|v| = (3^a|x| - B)/2^b`. ∎

In words: for a negative rational start, the whole "credit" question at
integer states is whether the carry ratio `beta` has passed the fixed
real number `delta = |x| - 1`. Lemma S with `v = -1` and `R = 1` gives
`R' >= 3/2`. That is the Applegate–Lagarias bound, and it proves that
`-1` is hostile.

### 4.2 The perturbation lemma (PROVED)

**Theorem P.** Let `h = -p/3^(j_h) < 0` satisfy

> **(H1)** every E-path from `h` has all values `< 0`, and running
> multiplier `> 1` at every halving.

Let `i >= 1`, `j >= max(j_h, 1)` and `c >= 1` satisfy

> **(B)** every E-path from `h` has made at least `j` multiplications at
> the moment of its `i`-th halving; that is, `j <= alpha_h(i) :=
> f_i(h mod 2^i)`;
>
> **(A)** `(1 - 3^-j)/2 > delta_h + 2^i c/3^j`, where
> `delta_h = |h| - 1`.

Then `x = h - 2^i c/3^j` satisfies (H1). In particular
`x in Bad_inf(E)`.

*Proof.* `x ≡ h (mod 2^i Z_2)`. Take an E-path from `x` with states
`v_t`. Let `u_t` be the states of the same move sequence applied to `h`.
Then `v_t - u_t = -3^(a_t) 2^(i - b_t) c/3^j`.

1. While `b_t <= i - 1` the difference lies in `2Z_2`, so the parities
   agree. Hence the moves up to and including the `i`-th halving form an
   E-path from `h`. By (H1), `u_t < 0` and `R_t > 1` at those halvings.
   Also `v_t < u_t < 0`: no escape and no descent up to the `i`-th
   halving.
2. By (B) the path has made `>= j` multiplications by its `i`-th
   halving. Let `tau` be the first time with `a = j`, so `b_tau <= i`.
   Then `u_tau = (-3^(j - j_h) p + B)/2^(b_tau)` is a 2-adic integer
   with a 2-power denominator, hence an integer, and `u_tau <= -1`. So
   `v_tau = u_tau - 2^(i - b_tau) c <= -2` is an integer.
3. `beta_tau >= (1 - 3^-j)/2 > delta_x`, by (A). By Lemma C,
   `R_tau > |v_tau|`. By Lemma S, nothing after `tau` descends, and all
   later values are integers `<= -1`. ∎

### 4.3 The first generation: `-1` is not isolated (PROVED)

**Lemma A.** For `i >= 2`, every E-path from `-1` has `3^a >= 2^(i+1) + 1`
at its `i`-th halving. Hence `3^(alpha(i)) >= 2^(i+1) + 1`, where
`alpha(i) := alpha_(-1)(i)`.

*Proof.* Let `a >= 1` count the multiplications made before the `i`-th
halving; the first move from `-1` is `x3`. The value right after the
`i`-th halving, `(-3^a + B)/2^i`, is an integer `<= -1` by Lemma S, so
`3^a - B >= 2^i`.

Suppose all `a` multiplications came before the first halving, so the
path begins `M^a H^i`. The value before halving would be
`-(3^a + 1)/2`, whose 2-adic valuation is at most `1`, so the second
halving would be illegal. Hence some multiplication has `b(s) >= 1`,
and `B >= (3^a - 1)/2 + 1`. So `(3^a - 1)/2 >= 2^i`. ∎

**Theorem F.** For every `i >= 3`, the point
`x_i := -1 - 2^i/3^(alpha(i))` lies in `Bad_inf(E)`. The `x_i` are
pairwise distinct, since `v_2(x_i + 1) = i`, and `x_i -> -1` in `Z_2`.
So `Bad_inf(E)` is infinite and `-1` is an accumulation point. Their real
values lie in `(1, 3/2)`.

*Proof.* Apply Theorem P with `h = -1` (which satisfies (H1) by
Lemma S), `j = alpha(i)` and `c = 1`. (B) holds by definition. (A) is
`3^alpha > 2^(i+1) + 1`. By Lemma A the inequality holds with `>=`, and
equality would give `3^a - 2^n = 1` with `n = i + 1 >= 4`, which is
impossible: odd `a` gives `n = 1`; `a = 2k` gives
`(3^k - 1)(3^k + 1) = 2^n`, so `n = 3`. ∎

**FINITE-EXACT (`..._dim_alpha.c`).**

* `alpha(i)` for `i = 1..56` is
  `1,2,3,4,5,6,7,7,7,7,8,9,9,10,11,11,12,12,13,14,14,15,16,16,17,18,18,
  19,19,20,21,21,22,23,23,24,24,25,26,26,27,28,28,29,30,30,31,31,32,33,
  33,34,35,35,36,36`.
* For `9 <= i <= 56`, `alpha(i) = ceil((i+1) log_3 2)`: Lemma A is sharp.
  So `R_min(i) = 3^alpha/2^i` lies in `(2, 6)` and
  `delta_i = 2^i/3^alpha` in `(1/6, 1/2)`.
* Examples: `-3211/2187 = -1 - 2^10/3^7` and
  `-793585/531441 = -1 - 2^18/3^12`. These are exactly the lone
  exceptional classes found at distance `2^-i` from `-1`.

### 4.4 The closure of `{-1}` under Theorem P (FINITE-EXACT; every point PROVED)

`..._dim_generations.py` iterates Theorem P over all `(i, j, c)`, with
`alpha_h` computed exactly, up to 2-adic scale `i <= 56`:

| generation | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| new points | 1 | 61 | 231 | 25 | **0** |
| real values `|x|` | 1 | `[1.0585, 1.4933]` | `[1.1951, 1.4996]` | `[1.3609, 1.4989]` | – |

Scale `<= 44` gives the same picture: `49/164/20/0`.

* The closure's covering counts `55, 77, 96, 119, 139, 163, 191, 213,
  249, 269, 317` at `L = 16, 20, ..., 56` are PROVED lower bounds for
  `K_L`. They grow about linearly, 6 per level.
* **Why it stops.** A generation-`(g+1)` child of `h` needs
  `R_min^h(i) = 3^(alpha_h(i))/2^i > 1/(1/2 - delta_h - 3^-j/2)`. The
  third generation has `delta >= 0.36`, so it would need
  `R_min > 7.2`, which never happens up to scale 56.
* **Examples.**
  * generation 1: `-97/81 = -1 - 2^4/3^4`;
  * generation 2: `-355/243 = -97/81 - 2^6/3^5`;
  * generation 3: the points with `delta` in `[0.36, 0.5)`.

**The closure is a proper subset of the hostile points.** For example,
`-13/9`, `-87209/59049`, `-781553/531441` and `-2374355/1594323` are
PROVED hostile by the exact prover but lie outside the closure.

* Against its nearest closure point `h`, `-87209/59049` reads
  `h - 2^11*7/3^10` with `h = -2699/2187`, and `j = 10 > alpha_h(11) = 8`.
  So (B) fails for that `h`: a path passes the 11th halving at a
  non-integer state, and its later values still never reach
  `(-1/3, 0)`.
* `-13/9 = -1 - 2^2/3^2` meets (B) but misses (A) by exact equality,
  since `9 = 8 + 1`.
* §4.5 accounts for all of these with a single tree.

### 4.5 A certified census of hostile rationals, and slow descenders (FINITE-EXACT)

**Exact prover (`..._dim_hostile_prover.py`).** It runs a DFS over
E-paths with exact values `N/3^e`. A branch is:

* **SAFE** by Lemma S: an integer `v <= -1` with `R > |v|`, or `v` odd
  with `R > |v| - 1/3`;
* an **ESCAPE**: a value `>= 0`, which is followed explicitly on
  positives to a descent;
* a **DESCENT**: a prefix with `3^a < 2^b`.

`hostile` is a complete certificate. For `|x| < 3/2`, i.e. `delta < 1/2`,
the unsafe tree is finite: after `T > log_3(1/(1 - 2 delta))`
multiplications, `beta >= 1/2 - 3^-T/2 > delta`. So the search always
terminates.

**Census (`..._dim_z13scan.py`).**

* **Method.** Take all `x = -p/3^j` with `j <= 25` and
  `1 <= p/3^j <= 1.6` that lie in a class of `Bad_61`. For each class and
  each `j` there is exactly one candidate, and random hits are below
  `8*10^-4`.
* **Hostile points.** Every candidate with `|x| < 3/2` is **PROVED
  hostile**: 343 points, the largest being `|x| = 1.49839`. None
  descends. Any hostile `-p/3^j` in this window lies in `Bad_61`, so
  **these 343 are exactly the hostile points `-p/3^j` with `j <= 25`
  and `1 <= p/3^j < 3/2`.** Per `j`:

  | j | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | hostile | 1 | 0 | 1 | 1 | 2 | 3 | 6 | 10 | 15 | 18 | 22 | 12 | 26 |

  | j | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | hostile | 14 | 15 | 20 | 11 | 26 | 12 | 18 | 17 | 12 | 26 | 12 | 29 | 14 |

  Block sums over `j = 9–13, 14–18, 19–23` are `92, 84, 85`: **bounded,
  no growth.** 186 of the 343 lie in the Theorem-P closure; 157 lie
  outside it, 6–13 per `j`.
* **Slow descenders.** The other 62 candidates have
  `3/2 < |x| <= 1.525`. The exact prover cannot certify hostility when
  `|x| > 3/2`: excursion-heavy branches never become SAFE. For four of
  them the excursion route (`M^a`, then Collatz on negatives with at most
  two extra `x3`-moves, down to `-1`) gives **PROVED explicit
  descents**:
  * `-9905/6561` at `(a,b) = (147, 233)`;
  * `-90025/59049`, `-89641/59049`, `-9953/6561` at `(41, 65)`.

  The first case is confirmed independently: the exact thread function
  gives `f_65(-90025/59049) <= A_65` and `f_i > A_i` for all `i <= 64`.
  All 62 are still exceptional at `m = 64`. So they are consistent with
  slow descenders whose first certificate needs precision `>= 65`; the
  other 58 are undetermined.

  Heuristic for the route: `M^a` reaches `w ~ -3^a(delta + 1/2)`, and a
  continuation back to `-1` ending in a long run of halvings has total
  multiplier about `(1 + B'/2^b')/(delta + 1/2) < 1` when
  `delta > 1/2`. The descent precision diverges as `|x| -> 3/2+`, which
  is why such points persist in `Bad_m`. Making this a proof needs a
  reachability lemma for `-2^k` (`k` large) in the negative E-graph.
* **Cross-check.** All 186 points of the Theorem-P closure (§4.4) with
  `j <= 25` belong to the census and are certified by the exact prover,
  whose code path is independent of Theorem P.
* **Extension to `j = 26, 27`** (on `Bad_63`, `..._dim_z13scan.py` with
  `JMIN=26`).
  * `j = 26`: 14 hostile, all by the prover, plus 3 candidates beyond
    `3/2`.
  * `j = 27`: 25 hostile, plus 7 candidates beyond `3/2`. 24 are
    certified by the prover. The prover stops at `10^8` nodes on
    `-11330911601771/7625597484987`, but that point is a canonical
    Theorem-P step satisfying (A) and (B) from a certified parent, so
    Theorem P proves it.

  The per-`j` counts stay in the same range (at most 29).
* **Outside the window.** No `-p/3^j` with `1/4 <= |x| < 1` and no
  positive `p/3^j <= 2` (`j <= 25`) lies in any class of `Bad_61`.
* **Coverage.** Every class of the projected set `Bad_61 mod 2^L` at
  `L = 20` (147 classes) and `L = 30` (270 classes) contains a `Z[1/3]`
  point with `j <= 25`: a certified hostile one (132 and 233) or a
  slow-descender candidate (15 and 37). At these levels nothing in the
  data calls for a hostile point outside `Z[1/3]`.

**The census is a tree of canonical perturbations of bounded depth
(FINITE-EXACT, `..._dim_census_tree.py`).** Give each census point
`x != -1` its **parent**: the census point `h` with smaller exponent
`j_h < j_x` that agrees with `x` to the greatest 2-adic depth
`i = v_2(x - h)`. Then, for all 342 points, without exception:

* **Scale.** `i >= 1.585 j_x - 3`; the deficit `1.585 j_x - i` is 1, 2
  or 3 (43, 269 and 30 points). Each child's scale also exceeds its
  parent's own scale, so the parent relation is a genuine tree rooted
  at `-1`.
* **Shape.** `x - h = ±2^i/3^(alpha_h(i)+e)` with `e in {-1, 0, 1}`.
  Measured by the relative cost `rho = (|x| - |h|) R_min^h(i)`:

  | `rho` | `+1` | `-1` | `+1/3` | `-1/3` | `+3` | `-3` |
  |---|---|---|---|---|---|---|
  | points | 312 | 15 | 6 | 4 | 4 | 1 |

  So 91% are exactly the canonical Theorem-P step
  `x = h - 2^i/3^(alpha_h(i))`, which satisfies (B). For all but one of
  these 312 steps, condition (A) also holds. The exception is `-13/9`,
  the equality case.

  The other 30 steps shift the exponent by one or perturb upward, and
  lie outside Theorem P. **All 30 have `j <= 14`**, spread over
  `j = 4..14`. For `15 <= j <= 25`, every one of the 197 census points
  is a canonical step, satisfying both (A) and (B), from its nearest
  lower-exponent census point.

  Every census point outside the closure (157) descends from `-13/9` or
  from one of these 30 early steps. For example, `-781553/531441` is a
  canonical step from the non-closure point `-216827/177147`.
* **Generations** (tree depth):

  | depth | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  |---|---|---|---|---|---|---|---|---|
  | points | 1 | 37 | 58 | 147 | 62 | 27 | 10 | 1 |

  **The maximal depth per `j` grows to 6 at `j = 12`, hits 7 once at
  `j = 14`, and stays at 5–6 for all `j = 15..25`.** Chains of depth
  about 12 would fit in this range of `j` if the real-value budget
  allowed them. Only 48 of the 343 nodes have children; `-1` has 37.

### 4.6 Where the credit construction fails, and a conditional countability statement

The heuristic said: a stretch near `-1` banks credit `>= 3/2`, so the
next segment only needs to avoid losing more than that, and concatenated
blocks give uncountably many hostile points. Its first half is right:
Lemma S is exactly this bookkeeping. **The step that fails is the
independence of the blocks.**

* For a start `x in Z[1/3]_(<0)`, Lemma C says the only credit that
  matters at every integer state is the real number `beta - delta`.
* On the cheapest branches, where all multiplications come first,
  `beta = (1 - 3^-a)/2 < 1/2`. So **the total credit is the fixed
  slack `1/2 - delta`**. "Credit `3/2` at `-1`" is the same number, read
  as `1/(1 - beta)` at `beta = 1/3`.
* A block (a 2-adic perturbation at scale `2^-i`) changes `delta` by its
  real size `2^i c/3^j`. Every block draws on the same slack, which is
  condition (A) of Theorem P.
* A block costs about `1/R_min^h(i)`. In all 342 census steps it costs
  at least `1/(3 R_min^h(i))`, since `|rho| >= 1/3`. For `h = -1`,
  `R_min(i) < 6` at every tested scale, so a first block costs more than
  `1/6`.
* The number of blocks is therefore limited. The census tree has depth
  at most 7, and at most 6 for `j >= 15` (§4.5). The Theorem-P closure
  has depth 3.
* Infinitely many blocks with summable costs would need
  `rho -> 0`, i.e. `j_k - i_k log_3 2 -> inf`: ever longer non-integer
  phases. In those phases halvings can drive values into `(-1/3, 0)`,
  and the orbit escapes. No census point has `|rho| < 1/3`.

**Conjecture G (finite seed).** Every hostile point `x in Z[1/3]` with
3-adic exponent `j >= 15` is a canonical Theorem-P child
`h - 2^i/3^(alpha_h(i))` of a hostile point `h` of smaller exponent,
with (A) and (B) both holding. Equivalently, `Bad_inf(E) ∩ Z[1/3]` is
the Theorem-P closure of the finite seed of 146 hostile points with
`j <= 14`. This is FINITE-EXACT-verified for `15 <= j <= 27`. For
`j = 26, 27`, all 39 points are canonical steps satisfying (A) and (B),
with depth deficit 1–2 (`..._dim_conjG_check.py`). If it
holds, the whole rational part of `Bad_inf` is governed by Theorem P.
Its depth is then bounded by the budget argument, provided
`R_min^h(i)` stays bounded.

**Conditional statement (HEURISTIC; hypotheses FINITE-EXACT-supported).**
Assume:

* **(i)** `Bad_inf(E)` is the closure of its points in `Z[1/3]`
  (coverage at `L = 20, 30`);
* **(ii)** the number of hostile points with 3-adic denominator exponent
  exactly `j` is bounded, `<= C` (bounded to `j = 25`);
* **(iii)** every hostile point of exponent `j` agrees 2-adically, to
  within `2^-(1.585j - O(1))`, with a hostile point of smaller exponent.
  This holds with `O(1) = 3` for all 342 census points. For a Theorem-P
  point `h - 2^i c/3^j` it follows from an upper bound on `R_min^h(i)`,
  which is `< 6` for `h = -1`.

Then the number of classes mod `2^L` meeting `Bad_inf` is
`O(C L)`. So `dim_B Bad_inf = dim_H Bad_inf = 0`. With a finite
generation depth, `Bad_inf` is also countable, of finite Cantor–
Bendixson rank. In the data the splitting tree has Strahler number 5
and the census tree has depth at most 7. None of (i)–(iii) is proved.
Proving (ii) needs uniform control of `alpha_h(i)`, which is a
"loops through `-1` with bounded ratio" statement, the mirror of
HYP-9122.

### 4.7 The backward side (PROVED lemmas, FINITE-EXACT census)

* **Lemma S_b.** Positive integers are closed under legal reverse moves.
  The move `k = 0` from `1` gives `0`, which is illegal. From an integer
  `y >= 1` with multiplier `R = 2^K/3^s`, every continuation has
  `R' >= 1/y`, since `2^K' y >= 3^s' + B'`. So `R > y` makes the branch
  safe.
* **Lemma C_b.** Write `B = sum_t 3^(t-1) 2^(K - K_t)` and
  `beta = sum_t 3^(t-1)/2^(K_t)`. Then `y = R(x - beta)`. In particular
  **a start `x <= 1` is safe at every integer state**, so for dyadic
  `x <= 1` only the non-integer phase matters.
* **Census.** Every projected class of `Bad_41` at `r = 10, 15, 20`
  (`26, 44, 72` classes) contains a positive dyadic point `p/2^e`
  (`e <= 40`) with value in `[1/2, 1.54]`; at `r = 25` the figure is
  98 of 102. `Bad_33` gives the same list. There are 98 such points in all, including `1` and `1/2`.
  Their per-`e` counts stay small (at most 15), with peaks at
  `e = 17, 25, 28, 33, 36`.
* **Tree shape.** The backward projected counts and Strahler number
  (constant `3`) match the forward picture.
* A backward version of Theorem P (perturbing `1` or `1/2` by
  `3^i c/2^e`) should follow from S_b and C_b in the same way. It is
  not written out here.
* **Observation.** Three numerators occur on both sides:
  `1`, `793585` (forward over `3^12`, backward over `2^19`) and
  `419868489953` (over `3^24` and `2^38`). They come from the convergents
  `19/12` and `38/24` of `log_2 3`. This is a near-coincidence, not a
  duality.

## 5. Consequences for E-SCC and next steps

* **Q2.** Verified below `7.87*10^17` (§2).
* **Escape lemmas.** The forward escape obligations are not "finitely
  many rationals plus `-1`". They are a comb of combs: generation-1
  points `x_i`, one per scale, and their descendants. Their number up
  to scale `L` grows polynomially, so an Applegate–Lagarias see-saw
  needs **parametrized** escape families, indexed by the generation
  tree. This is the refined form of the atlas's warning.
* **Next probes.**
  1. Prove (ii) for generation 1, i.e. an upper bound on
     `R_min^h(i)` for `h = x_i`, from the structure of loops through
     `-1`. Test Conjecture G beyond `j = 25`.
  2. Prove that `-p/3^j` with `|x| > 3/2` is never hostile, from a
     reachability lemma for `-1` in the negative E-graph with small
     final carry.
  3. Push the census to `j = 30` (needs the `m ~ 64` dump and a C
     prover).
  4. Write out backward Theorem P and generation 1 near `1/2`.

## 6. Reproduction

```bash
bash 04-computation/experiments/collatz_procgen_20260922_dim_run.sh           # QUICK, ~10 min, <0.3 GB
FULL=1 bash 04-computation/experiments/collatz_procgen_20260922_dim_run.sh    # hours; forward 0.8 GB, backward 0.4 GB
```

The quick mode runs:

* all validations of §1.3;
* forward to `m = 50` (`P0 = 26`) and backward to `r = 33`;
* the `alpha` table and the Theorem-P closure to scale 56;
* prover examples, and the census to `j = 18`;
* the analysis of §3.

The quick mode took 6 min on the mac-mini, and every check passed.
`FULL=1` runs the production parameters of this note. They took 3.0 h of
CPU forward (`m <= 64`) and 1.3 h backward (`r <= 41`), with other jobs
running.

The census extension of §4.5 is two separate commands (about 40 min, in
Python):

```bash
JMIN=26 python3 .../collatz_procgen_20260922_dim_z13scan.py <dumpdir>/fwd_bad_m63.txt 63 27 100000000
python3 .../collatz_procgen_20260922_dim_conjG_check.py <dumpdir>/fwd_bad_m63.txt 63 <census j<=25> <dim_alpha binary> 26,27
```

The `.out` file concatenates:

* the validation;
* the production logs, with per-level counts and minimal
  representatives;
* the analysis;
* the `alpha` table, the closure, the prover certificates and descents;
* the census, its tree statistics, and the list of 343 certified points;
* the TT comparison and the `j = 26, 27` extension.
