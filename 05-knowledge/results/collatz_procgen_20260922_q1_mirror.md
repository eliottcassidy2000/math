# The Q1 mirror: loops through -1 live on the lower best approximations of log_2 3, every exit from -1 costs more than 1, the endgame is the low binary digits of 3^A u, and the shared numerators are a straddle, not a duality

**Status.**
- **PROVED** (hand proofs below; each also checked by code where marked):
  - the loop equation through `-1` (legality automatic; the first run of `x3`-moves and the last run of halvings have odd length);
  - the sharp ratio bound. Every loop through `-1` has `2^(K+1) <= 3^a + 1`, with equality only for the basic loop `-1 -> -2 -> -1`. So `a >= 2` forces `3^a > 2^(K+1)`: the loop ratio exceeds `2`, and `a >= a0(K) := ceil((K+1) log_3 2)`. The bound `2*3^a/(3^a+1)` proposed in the task is correct but not sharp; it is exactly the cost of the pure climb `M^a`;
  - the cost formula, the splicing lemma, and the reduction of the mirror of HYP-9122 to the records of `eta(n) = ceil(n log_3 2) - n log_3 2`. These records are the **lower** best approximations of `log_2 3` (`19/12, 84/53, 569/359, 1054/665, ...`), where Q2's are the upper ones;
  - the height lower bound;
  - the exceptions: for `K <= 9` the only loops through `-1` are `(MH)^K` (exhaustive, with the loop equation as a certificate). So `amin(K) = K > a0(K)` for `K = 5..9`;
  - the exit lemma near `-1`. Every word consuming at most `m` bits of `n = 2^m u - 1` acts as at `-1`, so its multiplier is `> 1` (`w(-1) = -lambda + beta <= -1`). No path descends within `m+1` halvings. The exact escape price is `P(m) = 3^(A*(m))/2^(m+1) > 1`, and it equals `3^(eta(m+1))`, which lies in `(1,3)`, iff an `a0`-loop with `m` halvings exists. **So there is no direct escape from `-1` at any precision `m >= 2`: the Q1 base point behaves like Q2's `1/2`, not like Q2's `1`**;
  - the landing law: landings are uniform over all classes, and the landing depth is `v_2(A - lambda(u)) + 1`, with `lambda(u)` the 2-adic logarithm `log_3(-1/u)`. Hence the chain recursion;
  - **the endgame reduction:** a chain of hostile exits continues exactly as dictated by the low binary digits of `3^A u`;
  - the straddle identity for shared numerators;
  - the finite criterion for backward hostility below `3/2`;
  - no letter-count-linear word map preserves non-descent (using Gelfond–Schneider).
- **FINITE-EXACT:**
  - `amin(K) = a0(K)` for every `10 <= K <= 4000` (layered DP, cap `3*10^7`, with exact minimal heights);
  - explicit loops, each verified by an independent exact checker, for every `K in {2,3,4} ∪ [10,4000]` (3994 values). They are spliced from 5 base loops (`B_11, B_19, B_84, B_569, B_1054`) and 7 record cycles (`C_1, ..., C_1054`); the record loop `B_1054` (height `19,012,579`) and `C_1054` were reconstructed by divide and conquer;
  - `alpha(i) = a0(i)` for `10 <= i <= 4000`;
  - exit prices for `m <= 40` (cross-checked against the dimension lane's program);
  - the level-24 certificate landscape;
  - chain statistics for every `n = 7 mod 8` below `2^32`;
  - 117 backward dyadic points (every one of value `< 3/2` in the backward census, `e <= 45`) certified hostile. The 4 census points above `3/2` descend, with explicit certificates;
  - 11 of 11 backward generation-1 points at the lower clocks (`e <= 200`) and 25 of 72 forward candidates above `3/2` descend, with explicit certificates.
- **CITED:** Dupuy–Weirich (J. Number Theory 158 (2016) 268–280), for averaged equidistribution of the low binary digits of `3^n`. The bibliographic data were verified on the author's page and in arXiv:2601.12753. The original is paywalled and was **not read**. Inherited and not re-proved: the dimension lane's Lemma S, Lemma C, Lemma A and Theorems P and F; the loops lane's Q2 results; and the three-adic lane's cycle gate and sign law. Gelfond–Schneider (`log_3 2` is transcendental) is used only in section 6.6.
- **REFUTED:** a forward/backward duality sending hostile points to hostile points. The two nontrivial numerators the dimension lane found on both sides, `793585` and `419868489953`, pair forward hostile points with backward points that **descend**, with explicit certificates at depths 53 and 212. Within the census ranges (`j <= 27`, `e <= 45`) the only numerator shared by two hostile points is `1`.
- **OPEN:** the mirror of HYP-9122 for all `K`; whether any point above `3/2` is hostile, on either side; the pointwise digit statement of the endgame; termination of hostile chains for all integers; Q1, Q2 and Collatz.

Session `collatz-procgen-20260922`, Q1-mirror sub-lane (mac-mini). Scripts:
`04-computation/experiments/collatz_procgen_20260922_mirror_*`. Output:
[collatz_procgen_20260922_mirror.out](collatz_procgen_20260922_mirror.out). Reproduction is in
section 10. Parent notes: [loops and escapes](collatz_procgen_20260922_loops_and_escapes.md) (Q2),
[exceptional dimension](collatz_procgen_20260922_exceptional_dimension.md) (Q1 threads),
[synthesis](collatz_procgen_20260922_synthesis.md).

## 0. Setting and notation

Graph `E`: `x -> x/2` (x even) and `x -> 3x+1` (every x); on `Z_2` "even" means `v_2(x) >= 1`. A word
`w` over `M` (`x -> 3x+1`) and `H` (`x -> x/2`) with `a` letters `M` and `b` letters `H` acts affinely,
`w(x) = (3^a x + B)/2^b`, where `B = sum_{i<a} 3^(a-1-i) 2^(b_i)` and `b_i` is the number of `H` before
the `i`-th `M`. Its multiplier is `lambda_w = 3^a/2^b`, and `beta_w = B/2^b`. Legality of `w` at `x`
depends only on `x mod 2^b`.

Notation:
- `theta = log_3 2` and `L = log_2 3`;
- `a0(K) = ceil((K+1) theta)`, the least `a` with `3^a > 2^(K+1)`;
- `eta(n) = ceil(n theta) - n theta`, which lies in `(0,1)`.

Then `3^(a0(K))/2^K = 2*3^(eta(K+1))`, which lies in `(2,6)`.

The negative integers are closed under both moves (dimension lane, Lemma S). A **loop through `-1`**
is a legal closed walk `-1 -> ... -> -1` in this negative graph.

## 1. Loops through -1 (PROVED)

**Proposition 1.1 (loop equation; legality automatic).** Let `0 <= K_0 <= ... <= K_(a-1) <= K` and
let `w` be the word with `a` letters `M` in which the `i`-th `M` is preceded by `K_i` letters `H`, and
`K` letters `H` in all. Then `w` is a legal E-loop through `-1` iff

`2^K + B = 3^a`, where `B = sum_{i<a} 3^(a-1-i) 2^(K_i)`.

In that case every state is a negative integer, `K_0 = 0`, the initial run of `M` has odd length, and
the final run of `H` has odd length.

*Proof.* (⇒) The endpoint is `(-3^a + B)/2^K = -1`.

(⇐) Compute the states `x_t` formally in `Q`:
- forwards from `x_0 = -1`, each `x_t` lies in `Z[1/2]`;
- backwards from `x_T = -1`, via `x -> (x-1)/3` and `x -> 2x`, each `x_t` lies in `Z[1/3]`.

The identity says the forward computation ends at `-1`, so both computations give the same `x_t`, and
`x_t` lies in `Z[1/2] ∩ Z[1/3] = Z`. Every `H` is applied to `x_(t-1) = 2x_t`, which is an even
integer. Negativity follows from closure.

Parity structure:
- Mod 2 (`K >= 1`): `#{i : K_i = 0}` is odd. So `K_0 = 0`, and the initial `M`-run is odd.
- Mod 3: `2^K + 2^(K_(a-1)) = 0 (mod 3)`. So `K - K_(a-1)` is odd. ∎

(`..._mirror_loops_verify.py`: all 10,400,574 words with `K, a <= 12` satisfy "identity ⟺ legal loop
through -1", and every solution has the parity structure.)

**Proposition 1.2 (sharp ratio bound).** Every loop through `-1` satisfies `2^(K+1) <= 3^a + 1`. Equality
holds only for `(a,K) = (1,1)`, the basic loop `MH: -1 -> -2 -> -1` with ratio `3/2`. Hence for
`a >= 2`:
- `3^a > 2^(K+1)`;
- the ratio `3^a/2^K` exceeds `2`;
- `a >= a0(K)`.

The least conceivable ratio at `K` halvings is `2*3^(eta(K+1))`.

*Proof.* `B >= sum 3^(a-1-i) = (3^a-1)/2`, with equality iff all `K_i = 0`. So
`3^a - 2^K >= (3^a-1)/2`, i.e. `2^(K+1) <= 3^a + 1`.

Equality needs `3^a + 1` to be a power of 2:
- `a` odd gives `3^a + 1 = 4 (mod 8)`, so `(a,K) = (1,1)`;
- `a` even gives `3^a + 1 = 2 (mod 8)`, so `K = 0`.

For `a >= 2` we get `2^(K+1) <= 3^a`, and since `3^a` is odd, `2^(K+1) < 3^a`. ∎

*On the bound proposed in the task.* `2^(K+1) <= 3^a + 1` gives `3^a/2^K >= 2*3^a/(3^a+1)`, which is
correct. But the right-hand side is the cost of the pure climb `M^a` (Proposition 1.3), which closes
into a loop only for `a = 1`. The sharp statement is: ratio `3/2` for `MH`, and `> 2` otherwise.

**Proposition 1.3 (cost formula).** `3^a/2^K = prod_v 3|v|/(3|v|-1) = prod_v (1 + 1/(3|v|-1))` over the
`a` multiplication points `v`. The climb `-1 -> -2 -> -5 -> -14 -> ...` (points `-(3^t+1)/2`, `t < n`)
costs exactly `2*3^n/(3^n+1)`.

*Proof.* Around the loop the product of the step ratios is 1. An `M` at `v` contributes `(3v+1)/v`, an
`H` contributes `1/2`. The climb product telescopes:
`prod_{t<n} (3^(t+1)+3)/(3^(t+1)+1) = 2*3^n/(3^n+1)`. ∎

**Corollary 1.3a.** Every loop through `-1` has `a <= K`, with equality only for `(MH)^K`.

*Proof.* Each factor `1 + 1/(3|v|-1)` is at most `3/2`, with equality iff `v = -1`. So
`3^a/2^K <= (3/2)^a`, i.e. `2^a <= 2^K`. ∎

**Mirror of HYP-9122.** *For every `K >= 10` there is a loop through `-1` with exactly `a0(K)`
multiplications*, i.e. with the least conceivable ratio `2*3^(eta(K+1))`.
- It holds also at `K = 2, 3, 4` (`(MH)^K`).
- It fails exactly at `K = 1` (the basic loop, `a = 1 < a0 = 2`) and `K = 5..9`.

**Proposition 1.4 (the exceptions, PROVED by exhaustion).** For `K <= 9` the only loops through `-1` are
`(MH)^K`. Hence `amin(K) = K > a0(K)` for `K = 5..9`. For `K = 10, 11, 12` there are 2, 4 and 6 loops
with `a0(K)` multiplications (plus `(MH)^K`).

*Proof.* By Corollary 1.3a every loop has `a <= K`. By Proposition 1.1 loops are exactly the
nondecreasing `(K_i)` in `[0,K]` solving the identity. So for `K <= 12` it suffices to enumerate
`a <= 12`: 10,400,574 words, all checked. ∎

*Remark.* The first primitive loop occurs at `K = 10`: `M3 H1 M1 H1 M2 H1 M1 H7`, visiting
`-1,-2,-5,-14,-7,-20,-10,-29,-86,-43,-128,...,-2,-1`. Q2 has only one exception (`s = 1`). Q1 has
five more, because the `-1` neighbourhood is rigid: every primitive loop must climb to `-14` and
come back through a landing `-(2^j+1)/3` with `j >= 5`.

**Structure of every loop (PROVED).**
- The last multiplication lands on `-2^j` with `j` odd, from `-(2^j+1)/3`. These are
  `-1, -11, -43, -683, -2731, ...`; the value `-3` and all `j = 3 (mod 6)` are excluded, since those
  values are multiples of 3, which are unreachable from `-1`.
- Every primitive loop starts `-1 -> -2 -> -5 -> -14`. So it is the basic loop with a closed walk at
  `-2` spliced in. That walk has ratio `(4/3)*3^(eta(K+1))` for an `a0`-loop.
- This is the mirror of Q2's "trivial loop plus a closed walk at 4, ratio `(9/8)*2^(eps)`" and of its
  landings on `R_j = (4^j-1)/3`.

## 2. Records, splicing, hubs (PROVED reduction; FINITE-EXACT instance)

**Records.** Call `n` a record if `eta(n) < eta(n')` for all `n' < n`. These are exactly the `n` for
which `3^(ceil(n theta))/2^n` is a new minimum above 1, i.e. the **lower** best approximations
`n/a < log_2 3`. Exact list for `n <= 60000`:

`1/1, 3/2, 11/7, 19/12, 84/53, 569/359, 1054/665, 25781/16266, 50508/31867` (`n/a`)

with `eta = 0.369, 0.107, 0.0598, 0.0123, 1.90e-3, 9.70e-4, 3.97e-5, 2.32e-5, 6.61e-6`. Q2's records are
the **upper** best approximations `2/1, 5/3, 8/5, 27/17, 46/29, 65/41, 149/94, ..., 485/306, ...`
(`p/n`). The hard cases of the two problems therefore sit on opposite sides of `log_2 3`.

**Lemma 2.1 (splicing).** Suppose a loop `(K,a)` visits the value `h`, and `C` is a closed walk at `h`
with `p` halvings and `q` multiplications. Splicing `C` in gives a loop `(K+p, a+q)`, and the ratio is
multiplied by `3^q/2^p`.

If `a = a0(K)` and `q = ceil(p theta)`, the spliced loop is an `a0`-loop iff
`eta(K+1) + eta(p) < 1`. Then `eta(K+p+1) = eta(K+1) + eta(p)`.

*Proof.* `a0(K) + ceil(p theta) = (K+p+1) theta + eta(K+1) + eta(p)`, and this is an integer. ∎

**Theorem 2.2 (reduction to records; mirror of the loops note's Theorem 4.2).** Suppose that:
- for every record `n >= 11` there is an `a0`-loop `B_n` with `n-1` halvings;
- for every record `n` there is a cycle `C_n` with `(p,q) = (n, ceil(n theta))` through a value `h_n`;
- every `B_n` visits `h_(n')` for every record `n' <= n`.

Then every `K` with `K+1 >= 11` has an `a0`-loop.

*Proof.* Let `n*` be the largest record `<= K+1`. If `K+1 = n*`, take `B_(n*)`.

Otherwise put `m = K+1-n*`. We have `eta(K+1) > eta(n*)`. Also `eta(m) = eta(K+1) - eta(n*)`: both
sides are `= -m theta (mod 1)` and lie in `(0,1)`.

Greedy decomposition: while `m` is not a record, subtract the largest record `n' < m`; then
`eta(m - n') = eta(m) - eta(n')`. This writes `m = n_2 + ... + n_t` with records `n_j <= n*` and
`sum eta(n_j) = eta(m)`, hence `sum ceil(n_j theta) = ceil(m theta)`.

Splice the cycles `C_(n_j)` into `B_(n*)`. The hubs remain visited after each splice. The result has
`K` halvings and `ceil(n* theta) + ceil(m theta) = ceil((K+1) theta)` multiplications. ∎

The Q2 arguments carry over verbatim with `eta` in place of `eps`:
- the record base loops `B_n` are not splices of anything smaller;
- `B_n` is the basic loop plus a closed walk at `-2` whose ratio exceeds `4/3` by the factor
  `3^(eta(n))`.

The hub condition is automatic for the first three records. Every primitive loop visits `-1`, `-5`
and `-14`, and these are the hubs used below for `C_1`, `C_3`, `C_11`. So only `C_19, C_84, ...`
impose a real condition, namely that the base loops climb high enough.

**FINITE-EXACT instance.**
- Base loops `B_11, B_19, B_84, B_569, B_1054` (heights 43, 547, 23,947, 458,899 and 19,012,579),
  plus `(MH)^2`.
- Cycles:
  - `C_1 = {-1,-2}` at `-1`;
  - `C_3` at `-5`: the negative Collatz 2-cycle `{-5,-7}`;
  - `C_11` at `-14`: an E-cycle with an `E`-only arrow; the negative Collatz 7-cycle also realizes
    `(11,7)` through `-41`;
  - `C_19` at `-41`, `C_84` at `-365`, `C_569` at `-1094`;
  - `C_1054` at `-2391485 = -(3^14+1)/2`, height `29,732,963`.
- The closure produces loops for **every `K in {2,3,4} ∪ [10,4000]`** (3994 values). Each is verified
  independently: E-simulation from `-1`, the identity (Horner), `a = a0(K)` by exact powers, and the
  cost formula as the integer identity `prod v * 2^K = prod (3v+1)` (`..._mirror_loops_family.py`).
- The family fails exactly at `K = 5..9`, because `(MH)^2` does not visit `-5`, matching
  Proposition 1.4.
- `B_1054` and `C_1054` are too tall for checkpointing within 0.9 GB. `..._mirror_mitm.c`
  reconstructs them in `O(N)` memory: it recomputes forward and backward layers and splits at a
  midpoint state, Hirschberg-style (0.23 GB and 7.5 min for `B_1054`; 0.36 GB for `C_1054`).
  - `B_1054` climbs 15 steps, to `-(3^15+1)/2`.
  - It lands on `-(2^25+1)/3 = -11184811` and then falls through `2^25`.
  - The mirror of Q2's landings on `R_5 = 341` and `R_7 = 5461`.

**Hub cycles (FINITE-EXACT, `mirror_recon1 ... scan`).** Successive minima of `eta(p)` over cycles
`(p,q)` through the climb points `-(3^j+1)/2`. The scans run to `p <= 100` with cap `3*10^5`, and to
`p <= 600` with cap `2*10^6` at `-3281` and `-9842`. The `(569,359)` cycle at `-1094` comes from the
reconstruction of `C_569`, and `(1054,665)` at `-2391485` from section 9 of the `.out`:

| hub | minima `(p,q)` |
|---|---|
| `-1` | `(1,1)` |
| `-5` | `(3,2)`, `(68,43)` |
| `-14` | `(11,7)`, `(57,36)` |
| `-41` | `(11,7)`, `(19,12)` |
| `-122` | `(19,12)` |
| `-365` | `(19,12)`, `(84,53)` |
| `-1094` | `(19,12)`, `(84,53)`, `(569,359)` |
| `-3281`, `-9842` | `(57,36)`, `(84,53)`, `(569,359)` |
| `-2391485 = -(3^14+1)/2` | `(1054,665)`, ratio `1.000044` |

As in Q2, high enough hubs realize every record.

*The first three records are the negative Collatz cycles.* The records `(1,1), (3,2), (11,7)` are the
parameters of the three negative Collatz cycles `{-1}`, `{-5,-7}`, `{-17,...,-91}`. From `19/12` on,
the record cycles need `E`-only arrows (`3n+1` at even `n`).

## 3. The DP and the heights (FINITE-EXACT; bound PROVED)

`mirror_loops_dp.c` is a layered DP by halvings. It keeps, for every state `w` right after a halving,
the minimal `a` and then the minimal height `H` (the largest `|v|` over multiplication points). The cap
is `|v| <= 3*10^7`. States are pruned when `A_L(w) >= (L+1) theta + log_3|w| + 3`. This pruning is
exact for loops with `a <= a0 + 2`: continuing from `w` needs `3^(a'')|w| >= 2^(K-L)`.

The lexicographic DP is exact wherever `amin = a0`: an `a0`-loop is `a`-minimal at every state, else
Proposition 1.2 would be violated.

**Result:** `amin(K) = a0(K)` for every `10 <= K <= 4000`. The run takes 288 s and 0.54 GB. Every
minimal height lies below the cap, so the heights are exact.

Cross-checks:
- exhaustive enumeration (`K <= 12`);
- 79 minimal-height reconstructions (`K <= 80`) and 10 record objects, each rechecked by the Python
  verifier;
- the independent splicing family (`K <= 4000`), which realizes every `a0(K)` found by the DP with
  explicit, independently verified words.

**Corollary 3.0.** `alpha(i) = a0(i)` for `10 <= i <= 4000`, extending the dimension lane's
`9 <= i <= 56`. We have `alpha(i) <= amin(i) = a0(i) <= alpha(i)` by Lemma A. So Theorem F's points
`x_i = -1 - 2^i/3^(a0(i))` are hostile, with this explicit exponent, for all these `i`.

**Proposition 3.1 (height bound, mirror of Prop 5.1).** Let an `a0`-loop with `K+1 = n` have ratio
`< 9/4` (so it is primitive) and all multiplication points `|v| <= H`. Then

`H > (a0(K) - 2 - log_3(2H)) / (3 ln3 * eta(n))`.

*Proof.* The loop starts with a climb `M^(a_1)`, `a_1` odd and `>= 3`. By Proposition 1.3 the climb
together with the forced next point `-(3^(a_1)+1)/4` costs `2*3^(a_1+1)/(3^(a_1+1)-1) > 2`. Every
other point costs `>= exp(1/(3H))`. So `3^(eta)` exceeds `exp((a - a_1 - 1)/(3H))`. Finally
`(3^(a_1-1)+1)/2 <= H`. ∎

Hard cases (`eta < 0.006`; full table in `.out` section 1):

| `n = K+1` | `a0` | ratio | `eta(n)` | `H(K)` | proved bound | `H*eta*3ln3/a0` |
|---|---|---|---|---|---|---|
| 84 | 53 | 2.004181 | 1.90e-3 | 23,947 | 8,040 | 2.83 |
| 569 | 359 | 2.002133 | 9.70e-4 | 458,899 | 111,448 | 4.09 |
| 1054 | 665 | 2.000087 | 3.97e-5 | 19,012,579 | 5,062,252 | 3.74 |
| 1623 | 1024 | 2.002220 | 1.01e-3 | 1,229,387 | 306,828 | 4.00 |
| 2108 | 1330 | 2.000175 | 7.95e-5 | 20,008,879 | 5,065,926 | 3.94 |
| 2677 | 1689 | 2.002308 | 1.05e-3 | 1,877,195 | 487,440 | 3.85 |
| 3162 | 1995 | 2.000262 | 1.19e-4 | 20,008,879 | 5,071,399 | 3.94 |

The last column is 3.7 to 4.1 at every hard case with `n >= 569`. Q2 had 3.85 to 3.99. In both
problems the constant is about 1.3 in `H ≈ 1.3 * (#multiplications)/(ln(prime) * defect)`. Median
heights are 655, 2,467, 4,111 and 4,691 on `K` in `(0,1000]`, `(1000,2000]`, `(2000,3000]` and
`(3000,4000]`.

Minimal-height loops (`.out` section 2):
- they climb 3 to 9 steps;
- they wander;
- they mostly land on `-43 = -(2^7+1)/3`, then fall through `-128, ..., -2`. The other landings are
  `-11` and `-683`.

Example: `K=18` is `M7 H1 M1 H3 M1 H1 M1 H3 M1 H3 M1 H7`, height 547.

## 4. The escape near -1 (PROVED; prices FINITE-EXACT)

**Lemma 4.1 (words inside the thread).** Let `n = 2^m u - 1`, with `u` odd and `n > 0`, and let `w` be
a word with at most `m` halvings. Then:
- `w` is legal at `n` iff it is legal at `-1`;
- `w(n) = w(-1) + lambda_w (n+1)`, where `w(-1) = -lambda_w + beta_w <= -1` by Lemma S;
- hence `lambda_w >= 1 + beta_w > 1` for every nonempty word, and `w(n) = lambda_w n + beta_w > n`.

*Proof.* The parity at the `j`-th halving depends only on `x mod 2^j`, and `w` is affine. The `-1`-state
is a negative integer. `beta_w > 0` once `a >= 1`, and every nonempty word from odd `n` begins with
`M`. ∎

**Theorem 4.2 (exit and exact price).** Let `m >= 2` and `n = 2^m u - 1` with `u` odd.
- (i) Every E-path from `n` with at most `m` halvings stays above `n`.
- (ii) The state right after the `(m+1)`-th halving is `y = (v + 3^A u)/2`. Here `(v, A)` runs over the
  endpoints of `-1`-paths with `m` halvings and trailing multiplications that end at an **odd**
  `v <= -1`, and `A` counts all their multiplications. All of them satisfy `3^A >= 2^(m+1)|v| + 1`.
  Hence `y > n`: **no path from `n` descends within `m+1` halvings.**
- (iii) The least `y/n` tends to `P(m) := 3^(A*(m))/2^(m+1)` as `u -> inf`, where `A*(m) = min A`.
  We have `A*(m) >= a0(m)`, with equality iff an `a0`-loop with `m` halvings exists. In that case the
  unique optimal endpoint is `v = -1`, the exit is "loop, then one halving", `y = (3^(a0) u - 1)/2`,
  and `P(m) = 3^(eta(m+1))`, which lies in `(1,3)`.
- (iv) So `P(m) > 1` for every `m >= 2`, and a see-saw is needed. For `m = 1` the price is `3/4`.

*Proof.* (i) is Lemma 4.1. Values grow under `M`, so minima are attained right after halvings.

(ii) The prefix `w'` up to the `(m+1)`-th halving has `m` halvings, so
`w'(n) = w'(-1) + 3^A (n+1)/2^m = v + 3^A u`. This is even iff `v` is odd.

`M^A H^m` is illegal from `-1` for `m >= 2`, because `v_2((3^A+1)/2) <= 1`. So some multiplication
follows a halving, and `B >= (3^A+1)/2`. Hence `2^m|v| = 3^A - B <= (3^A-1)/2`. Then
`2(y - n) = u(3^A - 2^(m+1)) - |v| + 2`, which is `>= u+1` if `|v| = 1` and
`>= (2^(m+1)-1)(|v|-1)+2` otherwise.

(iii) At `A = a0(m)` we have `3^(a0) < 3*2^(m+1)`, so `2^(m+1)|v| + 1 <= 3^(a0)` forces `|v| < 3`,
i.e. `v = -1`. ∎

Exact values (`..._mirror_escape.py`, exact layered search from `-1` with safe pruning, cross-checked
against the dimension lane's `dim_alpha` as `A*(m) = f_(m+1)(2^m - 1)`, 39 of 39 agree):

| `m` | 2 | 3 | 4 | **5** | **6** | **7** | **8** | **9** | 10 | 12 | 13 | 18 | 29 | 37 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `A*` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 8 | 7 | 9 | 9 | 12 | 19 | 24 | 26 |
| `P(m)` | 1.125 | 1.688 | 2.531 | **3.797** | **5.695** | **8.543** | **12.814** | **6.407** | 1.068 | 2.403 | 1.201 | 1.0136 | 1.0824 | 1.0275 | 1.1559 |

In bold are the rigid precisions `m = 5..9`, where no `a0`-loop exists. There `A*` is attained by
`(MH)^m` together with non-loop endpoints (`v = -5` for `m = 6`; `-5, -7` for `m = 7`;
`-5, -7, -11` for `m = 8`; only `v = -5` for `m = 9`). The least landing uses the most negative `v`,
but the price `3^(A*)/2^(m+1)` is the same. For `m in {2,3,4} ∪ [10,4000]`,
`P(m) = 3^(eta(m+1))`.

**The mirror of Q2's prices.**
- Q2's 1-escape factor is `c(k-1)/3 = 2^(eps(k))/2`, in `(1/2,1)`: direct descent.
- Q2's 1/2-price is `2c(k-1)/3 = 2^(eps(k))`, in `(1,2)`.
- Q1's `-1` price is `3^(eta(m+1))`, in `(1,3)`.

Structurally, the Q1 exit mirrors the 1-escape: a loop through the base point, then the one move that
is illegal at the base point (a halving at an odd value, versus the move 0 at 1). Its price, however,
mirrors the 1/2-price. The reason is the Catalan constant `1/2` in `B >= (3^a-1)/2`:
- in Q1 it is cancelled exactly by the exit's factor `1/2`;
- in Q2 the exit's factor `1/3` beats it.

**Proposition 4.3 (landings; PROVED).**
- (a) For the loop exit `y = (3^A u - 1)/2`, the map `u (odd, mod 2^(r+1)) -> y (mod 2^r)` is a
  bijection. So landings are uniform over all classes.
- (b) The landing lies on the `-1` thread with precision `m' = v_2(y+1) = v_2(3^A u + 1) - 1`:
  - if `u = 1, 3 (mod 8)`, then `m' <= 1`;
  - if `u = 5, 7 (mod 8)`, let `lambda(u)` be the element of `Z_2` with `3^lambda = -1/u`. It exists
    because `-1/u = 1, 3 (mod 8)`, which is the closure of `<3>`. Then `m' = 0` if
    `A - lambda(u)` is odd, and `m' = v_2(A - lambda(u)) + 1` if it is even.
- (c) (Chain recursion.) If `m' >= 2`, then `y = 2^(m') u' - 1` with
  `u' = (3^A u + 1)/2^(m'+1)`. So repeated hostile landings obey
  `u_(t+1) = (3^(A*(m_t)) u_t + 1)/2^(m_(t+1)+1)`, where `m_(t+1) = v_2(3^(A*(m_t)) u_t + 1) - 1`.
  This mirrors Q2's `w_(t+1) = (2^(K_t+3) w_t - 1)/3^(j_(t+1))`. Values strictly increase along a
  chain (Theorem 4.2(ii)), so there is no cycle of hostile exits.

*Proof.* (b) `v_2(3^A u + 1) = v_2(3^(A-lambda) - 1)`, and 2-adic LTE gives `v_2(9^t - 1) = v_2(t) + 3`.
(118,499 pairs `(u,A)` checked.) ∎

After a landing with price `P`:
- if `y` is even, halving gives factor `P/2`, which is `< 1` iff `P < 2`;
- if `y = 0 (mod 4)`, the factor is `P/4`;
- if `y = 1 (mod 4)`, the path `M, H, H` gives `3P/4`, which is `< 1` iff `P < 4/3`;
- if `y = 3 (mod 4)`, the landing is hostile again.

In general `n` descends as soon as the class of `y` has a certificate of factor `f < 1/P(m)`.

**Level-24 landscape (FINITE-EXACT, `mirror_seesaw.c`; mirror of Q2's level 12).**
- There are **369 exceptional classes mod `2^24`**, equal to the dimension lane's `|Bad_24|` (an
  independent check).
- The worst certified factor is `3^5/2^8 = 0.949219`. Only 19 certified classes have factor `>= 0.9`,
  and 10,713 have `>= 1/2`.
- The worst factor by distance `d` (shared low bits with the exceptional set) stays in
  `[0.63, 0.95]` for all `d >= 4`, and equals `0.855` for most `d >= 16`.
- Q2's thread of 1 forces the worst factor at level `r` to be at least `max_k c(k-1)/3`, whose
  limsup is 1. The `-1` thread forces no such bound, because its exit costs more than 1.
- The one-round see-saw fails on the following fraction of landing classes mod `2^24`:
  - `2.2e-5` to `9.8e-4` for every `m in {2,3,4} ∪ [10,40]`;
  - `3.1e-3, 7.4e-3, 1.1e-2, 2.9e-2, 8.6e-3` at the rigid `m = 5..9`.

  Q2's figure was `<= 2.7e-4`.

**Chain statistics on integers (FINITE-EXACT).**
- *Canonical chains* (always take the exit of Theorem 4.2):
  - for 200,000 random `n = 7 mod 8` below `2^40`, `P(length >= k) = 0.2505, 0.0625, 0.0151, 0.0037,
    0.0010` for `k = 2..6`. This is exactly the `4^(1-k)` of the random model, since `m' >= 2` needs
    `u = 5,7 mod 8` and the right parity;
  - for all `n = 7 mod 8` below `2^22`, the longest chain has 9 links. The largest cumulative
    canonical price is 4494 at `n = 3563631` (precisions `4,2,5,7,7,6`), i.e. `n^0.557`. The
    expensive rigid precisions dominate it.
  - The canonical price per consumed bit is at most `ln P(8)/9 = 0.283` nats. For `m >= 10` it is
    at most `0.068`. Q2's figure was `0.114`.
- *Shortest descents* (`mirror_chain.c`, the mirror of `half_chain.c`): iterative deepening on the
  number of halvings, excursions `<= 1000n`, pruning by the exact oracle `f_22`. Over **all 536,870,912
  integers `n = 7 mod 8` below `2^32`**:
  - every one descends within `D <= 43` halvings; the maximum is at `n = 780669439`;
  - `D >= m+2` always, as Theorem 4.2 requires;
  - excess `D-(m+1) = 1` occurs exactly when `P(m) < 2` and the landing is even (e.g. half of all
    `m = 10`);
  - only **2.20%** of the shortest paths re-land on the `-1` thread at precision `>= 3` after the exit
    (4.89% at precision `>= 2`).

  These fractions are identical at `3*10^7` and `2^32`. Q2's `31.4%` counts any visit to the classes
  `1, 14 mod 27`, including the trivial transfer through 1, so the two numbers measure different
  things. The refined Q1 statistic shows that shortest descents rarely chain.

## 5. The endgame: low binary digits of 3^A u (reduction PROVED; digit statement OPEN; averaged version CITED)

**Proposition 5.1.** Along a chain of hostile exits (Proposition 4.3(c)), the next landing depth and the
next landing class are functions of the lowest binary digits of `3^(A*(m_t)) u_t + 1`:
- `m_(t+1) = v_2(3^A u_t + 1) - 1`, i.e. the agreement between the binary digits of `A = A*(m_t)` and
  those of the 2-adic number `lambda(u_t)`;
- `y mod 2^r` is determined by `3^A u_t mod 2^(r+1)`.

Once the digits of the starting integer are spent, each further link reads the next block of low binary
digits of `3^A u_t + 1`. So **the Q1 endgame is a statement about the low binary digits of `3^A u`**,
with `A = A*(m)` running over the Beatty-type exponents `a0(m) = ceil((m+1) log_3 2)` (except at the
rigid `m = 5..9`). This is the exact mirror of
Q2's "3-adic digits of `2^K w`". There, the 1/2-thread landing condition is
`K + 1 = mu(w) (mod 2*3^(j-1))`, with `2^(mu) = 1/w` in `Z_3^x`, and `v_3(4^t - 1) = v_3(t) + 1`.

*Proof.* Proposition 4.3, applied link by link. ∎

**What is and is not known.**
- (1) The low digits of `3^A` are periodic in `A`. The chain question is therefore whether the Beatty
  integers `a0(m_t)` keep agreeing, to many binary digits, with the 2-adic logarithms `lambda(u_t)`
  of the current cofactors. It is a transversality statement between an archimedean clock
  (`log_3 2`) and 2-adic logarithms. Q2's version reads the same way with 2 and 3 exchanged.
- (2) CITED: Dupuy–Weirich prove that the lowest `q`-adic digits of `p^n` (their title case: bits of
  `3^n`) are equidistributed **on average over `n`**, via a nonexistence theorem for "higher Wieferich
  primes". A non-averaged version would imply Erdős's conjecture on the ternary digits of `2^n`. The
  endgame needs a **pointwise** statement along the specific exponents `A*(m_t)`: OPEN.
- (3) Refinement of the Q2 wording. The synthesis calls Q2's endgame "Erdős ternary type". In both
  problems, what the chain actually reads are the **low** `p`-adic digits, where the averaged theory
  (Dupuy–Weirich) lives, not the full expansion that Erdős's problem concerns. Neither version is
  proved. This agrees with the parallel
  [transversality lane](collatz_procgen_20260922_transversality_foundry.md), which has the Q2 form of
  Proposition 4.3(b) (its "kappa formula" for `v_3(2^K - r)`) and separates the low digits (Q2) from the
  top digits (Erdős).

*Why Q1's chains are structurally a Collatz map with memory.* The chain
`(m,u) -> (m',u')`, `3^(a0(m)) u + 1 = 2^(m'+1) u'`, is an accelerated `3x+1`-type map whose multiplier
exponent depends on the previous valuation. Its termination for all positive integers is OPEN; the
data give `P(len >= k) = 4^(1-k)`.

## 6. The duality question (REFUTED as a hostile-to-hostile duality; straddle PROVED; coincidence list FINITE-EXACT)

**6.1 Backward finite criterion (PROVED).** Let `y0` be a positive dyadic 3-adic unit with `y0 < 3/2`,
and let `S* = max{S : 3^(-S) > 3 - 2y0}`. Every path state satisfies `R_s = (x_s + beta'_s)/y0` with
`beta'_s = B_s/3^s >= (1 - 3^(-s))/2`. Hence:
- integer states `z >= 2` have `R > (2 + 1/3)/y0 > 1`;
- the state 1 at depth `s > S*` has `R > 1`;
- positive integers are closed under legal moves (Lemma S_b);
- non-integer states have `K < e` and depth `< e theta` while `R > 1`, so they form a finite tree;
- an integer `z` at depth `s` can reach 1 by depth `S*` only if `z <= (3^(S*-s+1)-1)/2`.

So hostility below `3/2` is decided by a finite search. `..._mirror_bwd_prover.py` memoises each state
with the least multiplier seen; a larger multiplier on the same state is dominated. A reachable
negative state, whose descendants are all negative, is followed greedily; the search then reports a
descent or "undecided". It never occurred for the points certified below. The prover reproduces the
loops lane's seven dyadic hostile points.

**6.2 The census (FINITE-EXACT).** The backward census consists of every dyadic `p/2^e` (`e <= 45`)
with value in `[0.3,2]` in a class of `Bad_41`. That is 121 points; the dimension lane's 98 for
`e <= 40` are among them.
- All **117 below `3/2` are PROVED hostile**. Since every hostile point lies in `Bad_41`, these are
  exactly the hostile `p/2^e`, `e <= 45`, with value in `[0.3, 3/2)`.
- All **4 above `3/2` descend**. Each has an explicit certificate: first move `k_1`, then the greedy
  map, replayed exactly:

  | point | value | certificate |
  |---|---|---|
  | `793585/2^19` | 1.5136 | `2^84 < 3^53` |
  | `1675672562339/2^40` | 1.5240 | `2^103 < 3^65` |
  | `419868489953/2^38` | 1.5275 | `2^336 < 3^212` |
  | `1690578594467/2^40` | 1.5376 | `2^168 < 3^106` |

  **No hostile dyadic point with value in `(3/2, 2]` and `e <= 45` exists.**

**6.3 Shared numerators (FINITE-EXACT, complete for the stated ranges).** Compare the forward hostile
census (382 points, `j <= 27`, `|x|` in `[1, 3/2)`) and the forward candidates above `3/2` (72) with the
backward census (121). The shared numerators are exactly:

| `N` | forward point | status | backward point | status |
|---|---|---|---|---|
| 1 | `-1` | hostile | `1`, `1/2` | hostile |
| `793585 = 3^12 + 2^18` | `-1 - 2^18/3^12` (1.4933) | hostile (Theorem F) | `1/2 + 3^12/2^19` (1.5136) | **descends**, depth 53 |
| `419868489953 = 3^24 + 2^37` | `-1 - 2^37/3^24` (1.4866) | hostile | `1/2 + 3^24/2^38` (1.5275) | **descends**, depth 212 |
| `196249027 = 3^17 + 2^26` | `-1 - 2^26/3^17` (1.5197) | **descends**: `M^48`, then negative Collatz, `3^253 < 2^401` | `1/2 + 3^17/2^27` (1.4622) | hostile |

The last row is new: the dimension lane compared only hostile censuses.

**6.4 The straddle (PROVED).** Every nontrivial shared numerator is of generation-1 form
`N = 3^j + 2^(e-1)`. For such `N`, put `x = -N/3^j = -1 - 2^(e-1)/3^j`, `y = N/2^e = 1/2 + 3^j/2^e`
and `rho = 3^j/2^e`. Then

`(|x| - 3/2)(y - 3/2) = -(rho-1)^2/(2 rho) < 0`,

so **exactly one of the two points lies above 3/2**. *Proof:* `|x| = 1 + 1/(2 rho)` and
`y = 1/2 + rho`. ∎

- The numerator is shared because the roles of `3^j` and `2^(e-1)` are exchanged:
  - forward: base `-1 = -3^j/3^j`, perturbation `2^(e-1)/3^j`;
  - backward: base `1/2 = 2^(e-1)/2^e`, perturbation `3^j/2^e`.
- Both points belong to the generation-1 families at the same clock `(e,j)`. On the lower side
  (`rho > 1`, so `j = a0(e-1)`) the forward point is Theorem F's `x_(e-1)`: the 2-adic preimage of
  `-2` under an `a0`-loop with `e-1` halvings, since `w(x) = w(-1) + (3^j/2^(e-1))(x+1) = -2`. On the
  upper side, `j = a0(e-1) - 1` and the forward point lies above `3/2`.
- The regime is `|rho - 1|` small, i.e. `e/j` near a convergent of `log_2 3`:
  - `rho > 1` (lower side: `19/12`, `38/24`, `84/53`, ...): the forward point is hostile and the
    backward one lies above `3/2`;
  - `rho < 1` (upper side: `8/5`, `27/17`, `46/29`, ...): the other way round.

**6.5 The clocks, beyond the census (FINITE-EXACT).**
- Lower side, `1 < rho < 1.06`, `e <= 200`: **all 11 backward points descend.** The clocks are
  `(19,12), (38,24), (57,36), (76,48), (84,53), (103,65), (122,77), (141,89), (160,101), (168,106),
  (187,118)`. The certificates are `2^K < 3^s` at lower near-convergents: `84/53, 336/212, 271/171,
  420/265, 569/359, 4300/2713, 420/265, 569/359, 504/318, 2192/1383, 905/571`. Their forward partners
  are hostile by Theorem F with `alpha(e-1) = a0(e-1)` (Corollary 3.0).
- Upper side, `0.94 < rho < 1`, `e <= 70`:
  - backward hostile at `(8,5)` (`371/256`), `(27,17)` and `(46,29)`; undecided at `(65,41)` (state cap);
  - the forward partners `-371/243` and `-196249027/3^17` descend, via `3^147 < 2^233` and
    `3^253 < 2^401`. The forward partners at `(46,29)` and `(65,41)`: no route found.
- Forward candidates above `3/2` (`j <= 27`): the route `M^a` then negative Collatz (`a <= j+400`)
  gives explicit descents for 25 of the 72. The certificates are at upper near-convergents
  `149/94, 233/147, 298/188, 317/200, 401/253, 485/306, 970/612, 1203/759, 1371/865, 1539/971,
  2910/1836, 2994/1889`, mostly Q2-record clocks or their multiples. Together with the dimension
  lane's three extra routes (`-9905/6561`, `-89641/59049`, `-9953/6561`), 28 of the 72 are PROVED
  non-hostile; the rest are undetermined.
- (`-371/243` and `-43/27`, which the choice-ladder note reported as showing "no descent" within 2M
  nodes, are both outside `Bad_61`, so both descend.)

**The mirror of slow descenders.** Forward slow descents certify at **upper** near-convergents
(`3^a` just below `2^b`: `65/41`, `233/147`, `485/306`, ...). Backward slow descents certify at
**lower** ones (`2^K` just below `3^s`: `84/53`, `569/359`, ...). The same two sides separate the
loop records (section 2).

**6.6 No word duality (PROVED).**
- **Lemma.** No map sending forward words to backward words with letter counts
  `(s,K) = (gamma a + delta b, alpha a + beta b)` preserves the non-descent condition for all `(a,b)`.

  *Proof.* Non-descent is `3^a > 2^b` forward and `2^K > 3^s` backward. Both conditions are
  half-planes through the origin, and the boundary `a log3 = b log2` has irrational slope, so
  lattice points of the quadrant lie arbitrarily close to it on both sides. Agreement on all of
  them therefore forces proportional forms, `K log2 - s log3 = c (a log3 - b log2)` with `c > 0`.
  That is, `c = alpha theta - gamma = delta L - beta` with `L = 1/theta`. Multiplying by `theta`
  gives `alpha theta^2 + (beta - gamma) theta - delta = 0`.

  `theta = log_3 2` is transcendental (Gelfond–Schneider, classical, CITED). So `alpha = delta = 0`
  and `beta = gamma`, and then `c = -gamma <= 0`: a contradiction. ∎

  Concretely, the basic loop `MH` (ratio `3/2`) goes to a word of ratio `2/3` under the letter swap
  and under reversal alike.
- **Relation to the inherited cycle gate** ([three-adic G map](collatz_mod6_20260917_three_adic_g_map.md),
  Theorem 4.2; not re-proved). The identity `m_0(2^J - 3^L) = bB'(w)`, `B'(w) = B(reversed w)`, relates
  one and the same cycle read forwards and backwards. It sends the forward loops through `-1` to
  backward loops through `-1`, not to loops through `1` or `1/2`.
  - By the sign law `sign(m_0) = sign(b) sign(2^J - 3^L)`, forward `-1`-loops have `3^a > 2^K`, while
    backward `1`-loops have `2^K > 3^s`. They lie on opposite sides of `log_2 3`.
  - The shared pairs are distinct rationals with prime-power denominators `3^j` and `2^e`, not
    cycle-gate points `B/(2^J - 3^L)`.
  - Hence **the shared-numerator phenomenon is neither a word-reversal nor a cycle-gate identity at
    `19/12` and `38/24`**. It is the exchange-and-straddle coincidence of 6.4.
  - The joint carry identity S5 ([joint carry](collatz_mod6_20260921_g_negatives_joint_carry.md)) is
    not needed.

**Classification.** A duality sending hostile points to hostile points is **REFUTED**: the two
nontrivial census coincidences pair a hostile point with a descending one, with explicit certificates.
The coincidence list `{1, 793585, 419868489953, 196249027}` is **FINITE-EXACT** and complete for
`j <= 27`, `e <= 45`. The straddle identity (6.4) is **PROVED** and fixes the regime: at each
convergent-type clock at most one partner can lie below `3/2`. Whether the other partner is always
non-hostile is the OPEN question "no hostile point above `3/2`".

## 7. The mirror dictionary

| object | Q2 (backward, 3-adic) | Q1 (forward, 2-adic) | status of the Q1 row |
|---|---|---|---|
| hostile base points | `1`, `1/2` (positive dyadic) | `-1` only (negative, minus sheet) | PROVED (inherited Lemma S) |
| loop equation | `2^K = 3^s + sum 3^(i-1) 2^(K-K_i)`, legality automatic | `2^K + B = 3^a`, legality automatic, odd first and last runs | PROVED |
| ratio bound | `2^(K+1) > 3^(s+1)` (`s >= 2`), ratio `> 3/2` | `2^(K+1) < 3^a` (`a >= 2`), ratio `> 2` | PROVED |
| least ratio | `c(s) = (3/2) 2^(eps(s+1))` | `2*3^(eta(K+1))` | PROVED |
| loop HYP | HYP-9122, `s <= 6000` | `amin = a0` for `K >= 10`, `K <= 4000`; false for `K = 5..9` | FINITE-EXACT; exceptions PROVED |
| records | upper best approx. `2/1, 5/3, 8/5, 27/17, ..., 485/306` | lower best approx. `1/1, 3/2, 11/7, 19/12, 84/53, 569/359, 1054/665` | PROVED (arithmetic) |
| hubs | climb `4, 13, 40, 121, 1093` | climb `-1, -2, -5, -14, -41, -122, -365, -1094` | FINITE-EXACT |
| record cycles | `C_1, C_3, C_5, C_17, ...` | `C_1, C_3, C_11` = negative Collatz cycles; `C_19, C_84, C_569` use `E`-only arrows | FINITE-EXACT |
| reduction to records | Theorem 4.2 | Theorem 2.2 | PROVED |
| explicit family | `s <= 2500` | `K in {2,3,4} ∪ [10,4000]` | FINITE-EXACT |
| heights | `H >= (s-...)/(3 ln2 eps)`; ratio 3.85–3.99 | `H > (a-2-log_3 2H)/(3 ln3 eta)`; ratio 3.7–4.1 | bound PROVED; law conjectural |
| landing | `R_j = (4^j-1)/3`, fall `4^j` | `-(2^j+1)/3`, `j` odd `>= 5`, fall `2^j` | PROVED |
| base-point escape | 1: factor `2^(eps(k))/2 < 1` (direct) | none: every exit costs `> 1` | PROVED |
| price | 1/2: `2^(eps(k))`, in `(1,2)` | `-1`: `3^(eta(m+1))`, in `(1,3)`; rigid `m = 5..9`: 3.8–12.8 | PROVED + FINITE-EXACT |
| no-descent window | `k` moves from `1+3^k u` | `m+1` halvings from `2^m u - 1` | PROVED |
| landings | `2^K w`, all unit classes | `(3^A u - 1)/2`, all classes (bijection) | PROVED |
| two-move exits | `c < 9/4`, ... | `P < 2` (even), `P < 4/3` (`1 mod 4`) | PROVED |
| certificate landscape | level 12: 30 exceptional classes, worst 0.98654, fail `<= 2.7e-4` | level 24: 369 exceptional classes, worst 0.949, fail `2.2e-5`–`9.8e-4` (rigid `m`: up to `2.9e-2`) | FINITE-EXACT |
| chain recursion | `w' = (2^(K+3) w - 1)/3^j` | `u' = (3^(A*(m)) u + 1)/2^(m'+1)` | PROVED |
| chain statistics | 31.4% visit a hostile class mod 27 (crude) | 2.20% re-land (precision `>= 3`); `P(len >= k) = 4^(1-k)` | FINITE-EXACT |
| price per digit | 0.114 nats | 0.283 (rigid `m = 8`); `<= 0.068` for `m >= 10` | FINITE-EXACT |
| endgame | low ternary digits of `2^K w` | low binary digits of `3^A u` (`v_2(A - lambda(u))`) | reduction PROVED; statement OPEN; averaged CITED |
| perturbation census | 117 hostile dyadics (`e <= 45`), none above `3/2` | 382 hostile `-p/3^j` (`j <= 27`, dim lane) | FINITE-EXACT |
| slow descenders | above `3/2`, at lower clocks `84/53`, `569/359` | above `3/2`, at upper clocks `65/41`, `233/147` | FINITE-EXACT |
| duality | — | straddle, not duality | REFUTED / PROVED / FINITE-EXACT |

## 8. Hypothesis candidates (no HYP files created)

- **(M1) Mirror of HYP-9122.**
  - *Statement:* for every `K >= 10` there is a loop through `-1` with `a0(K)` multiplications.
  - *Record form:* for every eta-record `n >= 11` there are a base loop `B_n` and a cycle `C_n`,
    compatible at hubs.
  - *Evidence:* the statement holds for `K <= 4000`, with explicit record objects up to `n = 1054`
    (`B_1054`, `C_1054`).
  - *Next records:* `25781/16266`, `50508/31867`.
- **(M2) Height law.** `H(n-1) ≈ 1.3 * a0(n-1)/(ln3 * eta(n))` at the eta-records. The constant is
  1.25–1.36 for `n >= 569`, the same as Q2's 1.3.
- **(M3) No hostile point above 3/2, on either side.**
  - Backward evidence: 4 of 4 census points and 11 of 11 lower clocks with `e <= 200` descend.
  - Forward evidence: 28 of 72 candidates descend; the rest are undetermined.
  - Consequence: with the straddle identity, (M3) implies that the only numerator shared by two
    hostile generation-1 points is 1.
- **(M4) Chain termination.** Every positive integer's canonical exit chain is finite. The random
  model `P(len >= k) = 4^(1-k)` fits the data.
- **(M5) Endgame transversality (the digit form of M4).** For every positive integer, the canonical
  chain eventually reaches a link at which `u_t = 1, 3 (mod 8)`, or `A*(m_t) - lambda(u_t)` is odd,
  or `v_2(A*(m_t) - lambda(u_t)) = 0`. That is, the Beatty exponent disagrees with the 2-adic
  logarithm of the current cofactor in its lowest binary digits.
  - This is the pointwise low-digit statement whose average over exponents is the Dupuy–Weirich
    theorem.
  - The Q2 mirror concerns `K0(k)` and `mu(w) = log_2(1/w)` in `Z_3`.

## 9. What changed relative to the parent notes

- The dimension lane's "near-coincidence, not a duality" is now precise: it is a straddle, and the
  backward partners of `793585` and `419868489953` **descend**. The observation compared a hostile
  census with a census that still contained slow descenders.
- `alpha(i) = ceil((i+1) log_3 2)` extends from `i <= 56` to `10 <= i <= 4000`.
- The Q1 mirror of HYP-9122 holds with explicit, independently verified loops for every
  `10 <= K <= 4000`. The Q2 lane has explicit loops for `s <= 2500` and DP values for `s <= 6000`.
- The backward dyadic census is certified: 117 hostile points, 4 slow descenders.
- The choice-ladder note's open points `-371/243` and `-43/27` both descend.
- Synthesis section 4 ("the endgame is Erdős ternary type") is refined: in both problems the chains
  read low `p`-adic digits (section 5.3).

## 10. Reproduction

```
bash 04-computation/experiments/collatz_procgen_20260922_mirror_run.sh > 05-knowledge/results/collatz_procgen_20260922_mirror.out
```

- FULL (default): 38 minutes on one core (measured); peak memory 0.55 GB (the loop DP). The backward
  prover stays below 0.55 GB (state cap 2.5M); the divide-and-conquer reconstructions use 0.23 and
  0.36 GB.
- `QUICK=1`: `K <= 1200`, chains to `3*10^7`, and no section 9 (`B_1054`, `C_1054`, the family to `K <= 4000`).
- The duality section reads the dimension lane's dumps (`DIM`, default `scratch/procgen_dim`: `fwd30/`
  and `bwd17/`), which `OUTDIR=<dir> FULL=1 bash ..._dim_run.sh` regenerates.

| program | what it does |
|---|---|
| `mirror_loops_dp.c` | layered min-`a`/min-height DP for loops through `-1`, with pruning; optional reconstruction |
| `mirror_recon1.c` | checkpointed reconstruction of one loop or hub cycle; hub scans |
| `mirror_mitm.c` | divide-and-conquer (Hirschberg) reconstruction of the tall record objects `B_1054`, `C_1054` |
| `mirror_loops_verify.py` | exhaustive loop equation for `K <= 12`; independent word verifier |
| `mirror_loops_family.py` | records, the splicing closure, verification |
| `mirror_escape.py` | `A*(m)` and `P(m)`, exit lemma on integers, landing law, canonical chains |
| `mirror_seesaw.c` | level-`r` forward certificates, landing-failure counts |
| `mirror_chain.c` | shortest descents from the `-1` thread (mirror of `half_chain.c`) |
| `mirror_bwd_prover.py` | backward finite criterion below `3/2`; slow-descent search above |
| `mirror_duality.py` | censuses, shared numerators, straddle, clocks, forward routes |
| `mirror_run.sh` | runs everything into the `.out` |
