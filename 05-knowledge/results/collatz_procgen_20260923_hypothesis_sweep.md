# Hypothesis sweep HYP-9120–9127, 9132: HYP-9121 is refuted, the loop hypotheses become FINITE-EXACT far beyond the old ranges, and the rest stay open with exact placements

**Status.**
* **REFUTED: HYP-9121.** The number of exceptional classes of `E_{6 mod 8}` grows exponentially. For every `m >= 120` it is at least `2^(0.0536 m - 6)` (Theorem R, PROVED by hand; its two finite lemmas are checked by code). The converged covering numbers grow like `2^(0.17 m)` on `m = 25..45` (FINITE-EXACT).
* **FINITE-EXACT (extended), still OPEN:**
  * **HYP-9122** now holds for every `2 <= s <= 190535` (was `s <= 6000`). The proof uses the splicing theorem (loops note, Theorem 4.2), applied to explicit record objects `B_q`, `C_q` for every record `q <= 79335` of `log_2 3`. Each object is independently verified.
  * **HYP-9125** now holds for every `10 <= K <= 176249` (was `K <= 4000`), by the mirror theorem (Theorem 2.2) with record objects to `n <= 50508`.
* **PROVED (new):**
  * The **tightness lemmas** T1–T3: on every record object, the halving count of every prefix is forced by its endpoint. This yields an exact, 1-bit, sparse DP whose layers hold `O(q)` states instead of `O(height)`.
  * The **carry-mod-D L^2 obstruction.** Every Fourier argument that bounds the frequencies one by one fails to prove HYP-9122. The chain-sum problem sits in the sparse regime `#chains < D^2`.
  * The **ladder of implications** `HYP-9126 ⟹ HYP-9124 ⟹ Q2`.
  * The **cost criterion** for descents through the node `1`.
  * The **canonical-route dichotomy** at the generation-1 points above `3/2` on both sides. At the "easy" clocks (`c(i-1) < 9/4`, respectively `3^(a0(i)) < 2^(i+2)`) these points descend, conditionally on HYP-9122(`i-1`) respectively HYP-9125(`i`), hence unconditionally inside the FINITE-EXACT ranges. At the other clocks, every descent through the canonical landing is **impossible** (Propositions C3–C6).
  * `|Bad_m(E_S)| >= m - 2` for **every** `S`.
  * In `E_{6 mod 8}` the analogue of HYP-9125 fails for every `K >= 35`. Every loop through `-1` there has `3a >= 2K+1` (§2.4).
* **FINITE-EXACT (new):**
  * `E_{6 mod 8}` thread counts to `m = 64` and converged projections to `m = 59`. The code reproduces the lane counts at `m <= 27` and the dimension lane's E counts and projections exactly.
  * The generation-1 sweep of §4.4: all 589 generation-1 upper partners with `i <= 300` descend.
  * The tightness check on 49 record objects; independent re-verification of every stored object; hub checks.
* **CITED:** only as recorded in the parent notes (Applegate–Lagarias 2006; Lagarias 1985 and 2009; Bernstein–Lagarias 1996; Dupuy–Weirich; Nesterenko via Bertrand; BDGP 1996; Chudnovsky–Bertrand). No new literature was read in this sweep.
* **OPEN (unchanged status):**
  * HYP-9120: contains E-SCC (Q1 ∧ Q2).
  * HYP-9123: part of Lagarias's Periodicity Conjecture.
  * HYP-9124: implies Q2.
  * HYP-9126: implies HYP-9124.
  * HYP-9127: the 2-adic cubic theta value.
  * HYP-9132: p-adic Nesterenko.
  * Collatz.

Session `collatz-procgen-20260922`, hypothesis-sweep sub-lane (mac-mini, 2026-09-23). No HYP or THM file was created or edited; status changes are proposed in §8.

* **Scripts** (`04-computation/experiments/`):
  * `collatz_procgen_20260923_sweep_tightdp.c`: sparse tight DP, checkpointed reconstruction, `ES_FILTER` option;
  * `..._sweep_records.py`: record driver and independent verifier;
  * `..._sweep_splice.py`: explicit splicing closure and hub check;
  * `..._sweep_gen1.py`: generation-1 points above `3/2`;
  * `..._sweep_es_forward.c`: `E_S` thread DFS, adapted from the dimension lane's `dim_forward.c`;
  * `..._sweep_es_refute.py`: Lemmas F and G, Theorem R, cross-checks;
  * `..._sweep_run.sh`.
* **Output:** [`collatz_procgen_20260923_sweep.out`](collatz_procgen_20260923_sweep.out).

## 1. Triage (one line per hypothesis)

| HYP | relation to a known open problem | most plausible proof route | cheapest decisive test | status after this sweep |
|---|---|---|---|---|
| 9120 | Contains E-SCC = Q1 ∧ Q2. **Collatz ⟹ Q1.** Q2 ⟸ X_min (HYP-9124) ⟸ HYP-9126. The dimension clause (c) is implied by the box-dimension form of HYP-9121, which is now false, so that route is gone (§2.5). | Applegate–Lagarias see-saw with parametrised escape families. It is blocked by the sharp price per digit `c* = 0.1144` and by a low-digit (Dupuy–Weirich-type) endgame. | An integer in `Bad_inf` (none below `10^18` backward). For (c): converged E-projections. | OPEN; the thread of 1 is settled to depth `190536` (§3). |
| 9121 | None; it is a statement about the relaxation `E_{6 mod 8}`. As stated, it would imply the dimension clause of HYP-9120(c). | — | Thread counts at `m = 30..64` (done), and the forward orbit of `-1` in `E_S` (finite!). | **REFUTED** (§2) |
| 9122 | Not a known problem. It is the *existence* counterpart of the Collatz cycle equation: `D = 2^p - 3^q` must divide a monotone 2–3 chain sum. The famous problem is the *non-existence* of Collatz cycles; here about `2^(0.9q)` solutions are expected. | Equidistribution of chain sums mod `D`. The L^2 obstruction (§3.5) rules out frequency-by-frequency bounds. What is needed is either cross-frequency cancellation or a constructive recursion along the continued fraction. | The next record objects (done to `q = 79335`). | OPEN; FINITE-EXACT for `s <= 190535` |
| 9123 | **Part of Lagarias's Periodicity Conjecture:** `PC ⟹ C2 ⟹ HYP-9127`. It is the bounded-discrepancy supercritical slice of PC. | A new irrationality mechanism for words with `Dio <= eta` and no functional equation. The coupled Z-number form is Proposition M of the hard-class note. | HYP-9127, a single explicit word. | OPEN |
| 9124 | **⟹ Q2** (endgame Theorem 6.1). **⟸ HYP-9126** (Proposition C1). Under HYP-9122 it is equivalent to `S_1/2` (endgame Theorem 6.2). | An amortised see-saw that pays for chains of 1/2-landings. It is blocked by the 3/2 wall. | An integer `m >= 2` in `Bad_inf`; none below `10^18`. | OPEN; FINITE-EXACT `<= 10^18` (inherited) |
| 9125 | The Q1-side mirror of HYP-9122 (lower best approximations). It governs `alpha(i) = a0(i)`, and so it is the quantity that keeps E's hostile tree finite (§2.5). | Same as HYP-9122. | The next record objects (done to `n = 50508`). | OPEN; FINITE-EXACT for `K <= 176249` |
| 9126 | **⟹ HYP-9124 ⟹ Q2** (Proposition C1: an integer `m >= 2` is a dyadic point above `3/2`). The forward half contains "multiplicative Q1 on the minus sheet". | Near-optimal climbs, i.e. near-solutions of `2^A p - 2^B ~ 3^C`, at the hard clocks. At the easy clocks the proof is reduced to HYP-9122 and HYP-9125 (§4). | A generation-1 upper partner at a hard clock with no descent. | OPEN; partial results in §4 |
| 9127 | **⟸ C2 ⟸ PC.** It is equivalent to the 2-adic irrationality of `sum_k (2^10/3^9)^(k^3)`, the 2-adic cousin of the open irrationality problem for non-unit cubic theta values. | A 2-adic zero estimate (HC-CT1) or `U(1/3)` (cube-theta note). | None finite: heights `<= 2^4999999` are already excluded. | OPEN |
| 9132 | **⟸ a p-adic Nesterenko theorem** at `q = rho^2`. Only Chudnovsky–Bertrand's `tr.deg >= 2` is known p-adically. Corollary M (not both `theta_3(rho)` and `psi(rho^2)` algebraic) is the proved part. | A p-adic Nesterenko theorem. Rationality of the value is already excluded (Theorem Y). | None finite. | OPEN |

## 2. HYP-9121 is false: the forward `E_{6 mod 8}`-orbit of `-1` is a finite trap of rate `2/3`

`E_S`, `S = {6 mod 8}`: odd `x -> 3x+1` (forced); even `x -> x/2`; and, only at even `x = 6 mod 8`, the
excursion `x -> 3x+1 -> 9x+4`. `Bad_m(E_S)` is the set of classes mod `2^m` with no multiplicatively
descending certificate of precision `m` (choice-ladder note §1). HYP-9121 asserts `|Bad_m(E_S)| = 2^(o(m))`.

### 2.1 Two finite lemmas (PROVED; finite checks in `..._es_refute.py`)

**Lemma F (finite trap).** The set

`F = {-1, -2, -4, -5, -7, -8, -10, -14, -16, -20, -29, -32, -43, -64, -86, -128}`

is closed under the `E_S` moves, so every `E_S`-path from `-1` stays in `F`.

*Proof (check).* Among the evens of `F`, only `-2` and `-10` are `= 6 mod 8`; the others are `2, 4, 0 mod 8`
and can only halve. The edges are `-1 -> -2`; `-2 -> -1, -5`; `-5 -> -14 -> -7 -> -20 -> -10`;
`-10 -> -5, -29`; `-29 -> -86 -> -43 -> -128 -> -64 -> -32 -> -16 -> -8 -> -4 -> -2`. ∎

(In the full graph `E`, `-14 = 2 mod 8` and `-20 = 4 mod 8` may also multiply, and the orbit of `-1` is
infinite; that is where E's record loops through `-1` live.)

**Lemma G (rate `2/3`).** Weight each `x3`-move by `+3` and each halving by `-2`. The simple cycles of the
graph on `F` are `-1 <-> -2` (weight `+1`), the Collatz 2-cycle `{-5,-14,-7,-20,-10}` (`2` moves up, `3`
down, weight `0`) and the 15-cycle through `-2, -10, -29, -43, -128` (`6` up, `9` down, weight `0`). No cycle
has negative weight, and the least weight of a path from `-1` is `0` (the empty path). Hence **every prefix
of every `E_S`-path from `-1` satisfies `3a >= 2b`** (`a` multiplications, `b` halvings), and

`alpha_S(b) := min #x3-moves at the b-th halving over E_S-paths from -1  >=  ceil(2b/3)`.

(Exact values from the DP on `F`: `alpha_S(b) - 2b/3` lies in `[0, 1]` for `b <= 400`. Compare E, where
`alpha(b) = a0(b) = ceil((b+1) log_3 2)` for `9 <= b <= 4000`.)

### 2.2 The perturbation theorem for `E_S` (PROVED)

**Theorem P_S.** Let `i >= 4`, `j >= 1`, `c >= 1` with
* **(B_S)** `j <= alpha_S(i-3)`;
* **(A)** `2^i c/3^j < (1 - 3^(-j))/2`.

Then `x = -1 - 2^i c/3^j` lies in `Bad_inf(E_S)`: every `E_S`-path from `x` has all values `< 0` and
multiplier `> 1` at every halving.

*Proof.* We use the dimension lane's Lemmas S and C, stated for E-paths; every `E_S`-path is an E-path.
1. `x = -1 (mod 2^i)`. Run an `E_S`-path from `x` with states `v_t`, and apply the same move sequence to `-1`, with states `u_t`. Then `v_t - u_t = -3^(a_t) 2^(i-b_t) c/3^j`.
2. While `b_t <= i-3` this difference lies in `8 Z_2`. So `v_t = u_t (mod 8)`, and the parity and the `6 mod 8` test agree: the same moves are `E_S`-legal at both points.
3. By (B_S), each such path has made its `j`-th multiplication, at time `tau`, no later than its `(i-3)`-th halving.
4. Before `tau`:
   * the multipliers at halvings are those of an `E_S`-path from `-1`, hence `> 1` (Lemma S at `v = -1`, or Lemma G);
   * the values `v_t = u_t - (positive non-integer)` are `< 0`.
5. At `tau`, `v_tau = u_tau - 2^(i - b_tau) c`, an integer `<= -1 - 8c`, since `u_tau` lies in `F`.
6. The carry satisfies `beta_tau >= (1 - 3^(-j))/2 > 2^i c/3^j = |x| - 1` by (A). So Lemma C gives `R_tau > |v_tau|`.
7. Lemma S then shows that every continuation stays at integers `<= -1`, with running multiplier `> R_tau/|v_tau| > 1` at each halving. ∎

**Theorem R (HYP-9121 is false).** Put `c_max(i) = floor((3^(alpha_S(i-3)) - 2)/2^(i+1))`.
* **(Count.)** For every `m`, `|Bad_m(E_S)| >= max_(4<=i<m) min(c_max(i), 2^(m-i))`.
* **(Growth.)** By Lemma G, `c_max(i) >= lambda^i/18 - 1` with `lambda = 3^(2/3)/2 = 1.040042`. Hence

  `|Bad_m(E_S)| >= 2^(0.0536 m - 6)` for every `m >= 120`,

  and `liminf log2|Bad_m(E_S)|/m >= log2(lambda)/(1 + log2(lambda)) = 0.0536`.
* **(Dimension.)** The lower box dimension of `Bad_inf(E_S)` is `>= 0.0536`.

*Proof.*
1. For fixed `i` and `j = alpha_S(i-3)`, Theorem P_S applies to each `c = 1..c_max(i)`.
2. The points `-1 - 2^i c/3^j` are pairwise distinct mod `2^m` as long as the `c` are distinct mod `2^(m-i)`.
3. A class that contains a hostile point has no certificate of any precision. This gives the count.
4. For the growth, Lemma G gives `c_max(i) >= lambda^i/18 - 1`. Take `i = ceil((m + log2 18)/(1 + log2 lambda))`.
5. Then `min(c_max(i), 2^(m-i)) >= 2^(m-i) - 1 >= 2^(0.0536 m - 4.95) - 1`, which is `>= 2^(0.0536 m - 6)` for `m >= 120`. ∎

(Machine check of the exact bound: it exceeds `2^(0.0536 m - 6)` by at least 3.39 bits for every
`200 <= m <= 2000`, and `alpha_S(b) >= ceil(2b/3)` holds for every `b <= 2100`.)

Explicit values (`..._es_refute.py`):

| `i` | 60 | 80 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| `alpha_S(i-3)` | 39 | 52 | 65 | 99 | 132 |
| `c_max(i)` | 1 | 2 | 4 | 60 | 297 |

The lower bound `max_i min(c_max(i), 2^(m-i))` equals `278`, `12059` and `522659` at `m = 200, 300, 400`.

**Cross-checks (FINITE-EXACT).**
* Every family member with `i < 52` (3 points) lies in a class of `Bad_55` as computed by the independent thread DFS.
* The exact `E_S` prover (a DFS over all paths, with the Lemma S/C safety test) certifies the members at `i = 10` and `i = 60` hostile.
* At `i >= 80` the full path tree from `-1` exceeds the prover's cap of `5*10^6` nodes. Theorem P_S covers these members; the prover is only a check.

### 2.3 The counts (FINITE-EXACT; `..._sweep_es_forward.c`)

The thread DFS of the dimension lane, with the excursion allowed only at `6 mod 8`, reproduces three sets
of reference values exactly:
* the lane's `E_S` counts for `m <= 27`: `327, 1619, 2109, 2025` at the drop levels `16, 21, 24, 27`;
* in its `E`-limit (`J=1, S={0}`), the dimension lane's E counts (`1052, 1197, 561, 454, 1030, 809, 518` at `m = 30, 31, 32, 35, 40, 43, 46`);
* its E projections (`P_16..P_44 = 91, 147, 190, 250, 300, 354, 416, 468`).

New:

| drop level `m` | 27 | 32 | 35 | 40 | 43 | 46 | 51 | 54 | 59 | 62 |
|---|---|---|---|---|---|---|---|---|---|---|
| `E_S`: `|Bad_m|` | 2,025 | 5,292 | 5,545 | 15,924 | 16,908 | 15,627 | 44,688 | 51,017 | 165,171 | 199,337 |
| E: `|Bad_m|` (dimension lane) | 255 | 561 | 454 | 1,030 | 809 | 518 | 988 | 758 | 1,601 | 1,198 |

The converged projections `P_m(M) = |Bad_M mod 2^m|`, here with `M = 59`. `P_m(55)` and `P_m(59)` agree to within 1% for `m <= 40`, and `P_m(50)` is within 5% for `m <= 35`:

| `m` | 20 | 25 | 30 | 35 | 40 | 45 |
|---|---|---|---|---|---|---|
| `E_S` | 554 | 1,034 | 1,821 | 3,333 | 6,075 | 12,038 |
| E | 147 | 200 | 270 | 343 | 416 | 517 |

`E_S` grows by a factor of about `1.8` per 5 levels, i.e. `2^(0.17 m)`, with local power-law exponents
rising (`3.1, 3.9, 4.5`). E grows like `m^1.6`, with local exponents falling. The drop-level counts at
`m = 21..27` (`1619, 2109, 2025`), which looked flat and motivated HYP-9121, are a local plateau. Those counts
depend on `frac(m log_3 2)`.

**Where the hostile classes are.** Rational reconstruction of the `27,859` classes of `Bad_55 mod 2^50`
finds `5,371` points `-p/3^j` (`j <= 25`, `|x| <= 2.5`), against about `156` chance matches. They include
`-1, -13/9, -35/27, -97/81, ...` (the E census), all negative, with values in `[-2.4, -1.0]` and a mode
near `-1.65`. The counts per `j` grow (`84, 179, 354, 842` at `j = 10, 15, 20, 25`), whereas E's stay bounded
(11–29). Theorem R explains the part below `3/2`.

### 2.4 What exactly failed in the rationale of HYP-9121, and the mirror of HYP-9125 in `E_S` (PROVED)

The choice-ladder argument was that the `6 mod 8` excursion "interrupts every rising run". It does. But
near `-1` the interrupted run lands in the Collatz 2-cycle `{-5,-7}` and in the 15-cycle of `F`, whose rate
`2/3` exceeds `log_3 2 = 0.6309`. E reaches the record rates `12/19, 53/84, ...` (HYP-9125) through the
E-only arrows at the evens `= 2, 4 mod 8` (at `-14`, `-20`, ...), and `E_S` removes exactly those.

Precisely: every closed walk through `-1` in `F` leaves by `-1 -> -2` (`x3`) and returns by `-2 -> -1` (halving),
which together have weight `+1`. So **every `E_S`-loop through `-1` with `K` halvings has `3a >= 2K + 1`.**

Since `a0(K) = ceil((K+1) log_3 2) < (2K+1)/3` for `K >= 37`, the `E_S`-analogue of HYP-9125 fails for
every `K >= 37`. The exact DP on `F` shows that it holds exactly for

`K in {2, 3, 4, 10, ..., 17, 19, 20, 22, 23, 25, 26, 28, 31, 34}` (all `K <= 80` checked).

The tight DP with `ES_FILTER=1` agrees on every tested `K`. So the dichotomy is sharp:
* E: loops through `-1` of least ratio exist for all `10 <= K <= 176249` (§3.4). This keeps E's hostile tree finite.
* `E_S`: they stop at `K = 34`, and the hostile set explodes (Theorem R).

### 2.5 Consequence for HYP-9120(c) and E itself

* The evidence that `dim Bad_inf(E) = 0` is fragile in a precise sense. Theorem P_S works verbatim for E with `alpha` in place of `alpha_S`.
* The E-tree stays finite because `2^i/3^(alpha(i)) = 1/R_min(i)` is bounded below (by `1/6`). That bound holds because loops through `-1` with `a0(i)` multiplications exist, i.e. HYP-9125.
* Any excess `alpha(i) - i log_3 2 -> infinity` along a positive fraction of scales would already give exponentially many E-hostile points.
* So **dimension 0 for E needs the HYP-9125 loops**, now FINITE-EXACT to `K <= 176249`. This is the "uniform control of `alpha_h(i)`" that the dimension lane named as its missing step, here realised as a mechanism.

### 2.6 Which partial choices escape the mechanism

Theorem R needs only two things: that the forward `E_S`-orbit of `-1` is finite, and that its cycles have
rate `> log_3 2`. A computation on the orbit of `-1` (values `<= 10^6`) gives:
* **The orbit is finite** for `S = {6 mod 8}` (the 16-point trap `F`, rate `2/3`). For `S = {2 mod 8}`, `{0 mod 4}`, `{0 mod 8}` and `{54 mod 64}` it is `{-1, -2}`, of rate 1, and Theorem P_S then gives `c_max(i) ~ (3/2)^i`. **Theorem R's argument applies verbatim to all of these**, with even larger families when the rate is 1. This matches the choice ladder's finding that they gain "little" or "nothing measurable".
* **The orbit grows with the cap** for `S = {2 mod 4}` (984 values below `10^4`, 113,972 below `10^6`) and for E itself. The multiplication at `-14 = 2 mod 8` opens the trap. These are the two systems whose data look subexponential. For them the question stays OPEN, and it is tied to loops through `-1` of least ratio (§2.5).

## 3. HYP-9122 and HYP-9125: tightness, an exact sparse DP, and the extended ranges

### 3.1 Setting

Reverse moves from a start `h`:
* `x -> (2^k x - 1)/3` (Q2 side, sign `+`);
* `u -> (2^k u + 1)/3` (Q1 side, with `u = -v`).

Write `K_i` for the halvings of the first `i` moves and `B_i = sum_(t<=i) 3^(t-1) 2^(K_i - K_t)`. Then
`2^(K_i) h = 3^i x_i + B_i` (sign `+`), respectively `3^i u_i = 2^(K_i) h + B_i` (sign `-`), and
`B_i >= (3^i - 1)/2`.

Put `E_i = B_i - (3^i-1)/2 >= 0` and `G_i = E_i/(2^(K_i) h)`. From `B_(i+1) = 2^k B_i + 3^i`:

**(R)** `G_(i+1) = G_i + (1 - 2^(-k)) nu_i`, where `nu_i = (3^i - 1)/(2^(K_i+1) h)`. So `G` is nondecreasing along every path.

Moreover
* `G = 1 - (3^i(2x+1) - 1)/(2^(K+1) h)` (sign `+`);
* `G = (3^i(2u-1) + 1)/(2^(K+1) h) - 1` (sign `-`).

### 3.2 Tightness lemmas (PROVED)

**T1 (record cycles, Q2 side).** Let `p = ceil(q log_2 3)` and let `h -> ... -> h` be a cycle with `q` moves
and `p` halvings. Then every prefix has `K_i = kappa(i, x_i) := min{m : 2^m h >= 3^i x_i + (3^i-1)/2}`. So
`0 <= G_i < 1/2`, and `G_i <= G_q = 1 - (3^q(2h+1) - 1)/(2^(p+1) h)`.

*Proof.* `K_i >= kappa` from the prefix identity. The suffix `x_i -> h` gives `2^(p-K_i) x_i >= 3^(q-i) h + (3^(q-i)-1)/2`.
If `K_i >= kappa + 1`, multiplying the two bounds gives `2^p h x_i > 2 * 3^q h x_i`, contradicting `2^(p-1) < 3^q`. ∎

**T2 (record loops through 1).** Let `s >= 2` and `K = K0(s) = floor((s+1) log_2 3)`. On every loop through 1,
`K_i = kappa(i, x_i)` (with `h = 1`) and `G_i <= G_s = 1 - (3^(s+1) - 1)/2^(K+1)`.

*Proof.* The suffix gives `2^J x_i >= (3^(r+1)-1)/2` and `2 x_i <= 3^(r+1) - 1` (from `B' <= 2^J (3^r-1)/2`).
If `K_i >= kappa + 1`, then
`2^K >= (2*3^i x_i + 3^i - 1)(3^(r+1)-1)/(2x_i) >= 3^(s+1) - 3^i + (3^i - 1) = 3^(s+1) - 1`.
But `2^(K0(s)) < 3^(s+1)` and `3^(s+1) - 1` is not a power of 2 for `s >= 2`. ∎

**T3 (Q1 side).**
* **(a)** Loops through `-1` with `a0(K)` multiplications, `n = K+1`, `eta(n) < log_3(3/2) = 0.369`: every prefix has the largest conceivable `K_i`.
* **(b)** Cycles `(p,q) = (n, ceil(n log_3 2))` with `eta(n) < log_3 2`: the same.

*Proof.* The same two-sided squeeze, now from above: the prefix gives `2^(K_i) h <= 3^i u_i - (3^i-1)/2` and the
suffix gives `2^(J) u_i <= 3^r h - (3^r-1)/2` (`J` = halvings of the suffix, `r` = its moves).
* **(a)** A deficit of one halving gives `2^K <= 3^i (3^r + 1)/4`. The loop has `2^K = 3^a/(2*3^eta)`, so `2/3^eta <= 1 + 3^(-r) <= 4/3` for `r >= 1`, i.e. `eta >= log_3(3/2)`.
* **(b)** A deficit of one halving gives `2^p h u_i <= (3^i u_i)(3^r h)/2`, i.e. `2^(p+1) <= 3^q`. This contradicts `3^q = 2^p 3^eta < 2^(p+1)` (`eta < log_3 2`). ∎

**Corollary (exact 1-bit DP).**
* A record object exists iff its target is reachable in the digraph of tight states `(i, x)` with `G(i,x) <= G_fin`, where `K` is a function of `(i,x)`.
* Every intermediate state lies in a window of relative width `G_fin`. So a layer holds at most `2 G_fin X/(1 - G_fin) + O(i + log X)` values below a cap `X`.
* With `X = lambda * q/G_fin` this is `O(lambda q)` states per layer, **independent of the height** (`~q/eps`). The old DP stored every value up to the height.
* The cap makes the search one-sided: found ⟹ exists.

Machine check (`..._sweep_splice.py tightcheck`, exact integers): all 49 record objects with record `<= 6000`
are tight at every prefix. These are 37 on the Q2 side and 12 on the Q1 side, including the hand-fixed `C_1`, `B_3` and `C_3`.

### 3.3 Implementation and verification

`..._sweep_tightdp.c` keeps each layer as a sorted array. It updates `G` and `nu` incrementally in double precision (errors `<= 1e-11` relative, against a pruning margin of `1e-9 G_fin`), deduplicates with a radix sort, writes checkpoints every 64 layers, and reconstructs by recomputing segments from the checkpoints with a looser margin; the recomputed segments are supersets of the true layers, and every state in them is genuine. Every reconstructed object is re-verified by `verify()` in `..._sweep_records.py`, a separate code path that checks:
* legality of every move;
* the return to `h`;
* the halving count;
* a forward E-simulation;
* the carry identity `2^K h = 3^s h + B`, respectively `3^s h = 2^K h + B`, in exact integers.

Base loops are required to contain the climb `1 -> 4 -> ... -> y_N` (respectively `-1 -> -2 -> -5 -> ...`),
and `B_hubs_visited_max` records how far. Cycles are searched at the lowest climb hub `y_n` (`w_n`) that
leaves `G_fin >= 0.25 eps ln 2`. The hub hypothesis of Theorems 4.2 and 2.2 is checked explicitly
(`..._sweep_splice.py hubcheck`). Spliced loops are built explicitly at sample lengths and verified.

### 3.4 Results (FINITE-EXACT)

**The reformulation used throughout (PROVED, from Proposition 1.1 of the loops note).** Put
`R_s = {sum_(j<s) 3^j 2^(f_j) : f_0 >= ... >= f_(s-1) >= 0}`. Then `R_s = U_(a>=0) 2^a (3^(s-1) + R_(s-1))`, and
**HYP-9122(s) ⟺ `2^(K0(s)) - 3^s in R_s`**. (The target is odd, which forces `a = 0` and `e_(s-1) = e_s = 0`.)

**Q2 side (HYP-9122).**
* **Objects.** For every record `q` of `eps` with `3 <= q <= 79335` (the 35 records `3, 5, 17, 29, 41, 94, 147, 200, 253, 306`, `971 + 665j` for `j = 0..22` (i.e. `971, ..., 14936, 15601`), `47468, 79335`), the tight DP found and reconstructed:
  * a base loop `B_q` (length `q-1`, `K0(q-1)` halvings) that climbs at least to `y_2` (`q = 5`), `y_8` (`q = 17`), `y_10` (`q = 29`), and `y_14` or more for every `q >= 41`, `y_20` for `47468` and `y_22` for `79335`;
  * a cycle `C_q = (q, ceil(q log_2 3))` through a climb hub.
* **Hubs.** `4` (`q = 3`); `13` (`5`); `40` (`17, 29`); `121` (`41, 94, 147`); `364` (`200, 253`); `1093` (`306` to `6291`); `3280` (`6956` to `12276`); `9841` (`12941` to `14271`); `29524` (`14936`); `88573` (`15601, 47468`); `265720` (`79335`).
* **Checks.**
  * `C_1` (trivial) and `B_3` (trivial twice) are fixed by hand.
  * Every object passes the independent exact check (`reverify`).
  * Every `B_q` visits the hub of every `C_(q')` with `q' <= q` (`hubcheck`, 0 failures).
* **Conclusion.** By Theorem 4.2 of the loops note (PROVED), **HYP-9122 holds for every `2 <= s <= 190535`**, since `s+1` stays below the next record `190537`.
* **Explicit spliced samples.** Exact spliced loops, built by the rule in `..._sweep_splice.py` and each verified independently, at `s = 2, 3, 40, 305, 6001, 6290, 12345, 15600, 20000, 47466, 47467, 79333, 79334, 120000, 190535`. The largest has `190535` multiplications, `301992` halvings and height `1.33*10^11`; it is `B_79335` with nine record cycles spliced in, among them `C_79335` and `C_15601`.

**Q1 side (HYP-9125).**
* **Objects.** Base loops `B_n` for the eta-records `n = 11, 19, 84, 569, 1054, 25781, 50508`. The largest, `B_50508`, has `31867` multiplications and climbs at least to `-(3^20+1)/2`.
* **Cycles.** `C_n` at the hubs `-14` (`n = 11`), `-122` (`19`), `-365` (`84`), `-1094` (`569`), `-29525` (`1054, 25781`) and `-265721` (`50508`). `C_1 = {-1,-2}` and `C_3 = {-5,-7}` are fixed by hand.
* **Checks.** Every object is independently verified; `hubcheck` reports 0 failures.
* **Conclusion.** By Theorem 2.2 of the mirror note, **HYP-9125 holds for every `10 <= K <= 176249`**. The next eta-record is `176251`.
* **Explicit spliced samples** are verified at `K = 10, 18, 4001, 25779, 25780, 40000, 50506, 50507, 100000, 176249`. The last has `111202` multiplications.

**Consequences.**
* The 1-escape exists at every depth `k <= 190536`, and so does the canonical escape `Psi` of the endgame note.
* `alpha(i) = a0(i)` for `10 <= i <= 176249` (mirror note, Corollary 3.0), so Theorem F's hostile points `-1 - 2^i/3^(a0(i))` carry an explicit exponent in that range.

**Cost** (one core, other lanes running). Every Q2-side record `<= 47468` took about 12 min in total, including reconstruction and verification; `79335` took 16.5 min (`B` 7.4, `C` 9.0). The Q1 side to `50508` took 4 min. The largest layer had `106,160` states (`C_79335`) and the peak RSS was below 150 MB. For comparison, the old DP needed 361 MB and 8 min for `s <= 6000` and cannot reach heights `~10^10`. The object heights reported in the `.out` are those of the reconstructed paths, which run close to the caps. They are not minimal heights.

| record `q` | 971 | 5626 | 14936 | 15601 | 47468 | 79335 |
|---|---|---|---|---|---|---|
| `eps(q)` | `1.41e-3` | `9.71e-4` | `8.92e-5` | `2.62e-5` | `1.58e-5` | `5.29e-6` |
| cap of `B_q` | `1.0e7` | `2.2e7` | `6.3e8` | `2.2e9` | `1.1e10` | `5.6e10` |
| max states per layer of `B_q` | 4,319 | 6,512 | 17,269 | 18,088 | 55,101 | 92,013 |
| hub of `C_q` | 1093 | 1093 | 29524 | 88573 | 88573 | 265720 |
| max states per layer of `C_q` | 1,258 | 7,538 | 19,873 | 20,812 | 63,517 | 106,160 |

### 3.5 Why the existence proof is hard: the sparse regime (PROVED remark)

* **The counting problem.** For a record `q` let `N = C(p+q-1, q-1) ~ 2^(2.489 q)` be the number of monotone chains and `D = 2^p - 3^q`. A record cycle exists iff `D` divides some chain sum `B`. Let `S(r) = sum_chains e(rB/D)`.
* **Why frequency-by-frequency bounds fail.** Any argument of the form `#{B = 0 (mod D)} >= N/D - (1/D) sum_(r != 0) |S(r)|` needs `sum_(r != 0)|S(r)| < N`. Parseval gives `sum_r |S(r)|^2 >= D N`. Even with the uniform-model collision count, Cauchy–Schwarz leaves an error `~ sqrt(N)`, which beats the main term `N/D` unless `N > D^2`.
* **We are in the sparse regime.** Here `N ~ 2^(2.49 q) < D^2 ~ 2^(3.17 q)` for all `q >= 60`. Even square-root cancellation at every single frequency (`|S(r)| ~ sqrt N`) is not enough, since `sqrt N > N/D`.
* **What pointwise bounds would need.** Bounding `max_(r != 0)|S(r)|` alone needs `|S(r)| < N/D` at every nonzero `r`: saving a factor `D = 2^(Theta(q))` at every frequency. Tao-type 3-adic Fourier decay for Syracuse variables (as recorded in the parent notes) is polynomial, and it is proved on average; it is nowhere near this.
* **What a proof needs.** Either cancellation *across* frequencies (the multiplicative structure `r -> 2r, 3r` of the problem is the obvious lever), or a constructive argument. The tightness lemmas turn the constructive route into a reachability question in a window digraph of width `O(q)` per layer. A two-sided "expansion" statement would suffice by pigeonhole. It would say that at some middle layer, both
the states reachable from the start and the states from which the target is reachable fill more than half of
the window. It is OPEN.
* **The same for the loop form.** The representation form of HYP-9122(s), `2^(K0(s)) - 3^s in R_s = {sum_(j<s) 3^j 2^(f_j) : f nonincreasing}`, sits in the same regime.
* **Numerics.** `R_s` covers all but an exponentially thin set of the non-multiples of 3 in `[0.5, 3] 3^s`. At `s = 16` there are `3,757` uncovered values, of relative density `~2^(-s)`, concentrated at the lower edge `1/2`. The record targets `2^(K0(q-1)) - 3^(q-1)`, of ratio `1/2 + O(eps)`, sit exactly there. So a covering induction cannot prove HYP-9122 either; it can at best prove it for `s` with `c(s)` bounded away from `3/2`, and even that is not uniform (section S8 of the `.out`).

## 4. HYP-9126: the 3/2 wall

### 4.1 The ladder (PROVED)

**Proposition C1.** HYP-9126 (backward half) implies X_min (HYP-9124), which implies Q2.

*Proof.* An integer `m >= 2` prime to 3 is a positive dyadic 3-adic unit (`e = 0`) of value `> 3/2`. So
HYP-9126 says it is not hostile, which is X_min. Then apply the endgame note, Theorem 6.1. ∎

The forward half, read on the integers `-m <= -2`, is "every negative integer has an E-path with `3^a < 2^b` at
some halving" (multiplicative Q1 on the minus sheet). So **HYP-9126 is at least as hard as X_min.** Its new
content lies in the non-integer points above `3/2`.

### 4.2 The cost criterion (PROVED)

Let `x > 0` be a 3-adic unit in `Z[1/2]` and let `x = x_0 -> x_1 -> ... -> x_s = 1` be a legal reverse path. Then

`2^K/3^s = C/x`, where `C = prod_(t=1..s) (1 + 1/(3 x_t))`

(telescope `x_(t-1) = (3/2^(k_t)) x_t (1 + 1/(3x_t))`). So **a path to `1` descends iff `C < x`**. Every such
`C` exceeds `4/3` (the node `1` contributes `4/3`). By the carry bound, `C >= 3/2 - 3^(-r)/2`, where `r` is the
number of moves after the first integer node. Above `3/2`, descent therefore needs a near-optimal climb: the
path, read forwards from `1`, must stay high until its last `x3`-step.

### 4.3 The canonical-route dichotomy at the generation-1 points (PROVED)

Backward points. Let `T_i^(-1)(2) = 1/2 + 3^i/2^a`, `a = K0(i-1)`. Its value is `1/2 + 3/c(i-1)`, in `(3/2, 5/2)`.
It is `Psi`-mapped to `2` with price `2c(i-1)/3`.

**C3 (easy clocks).** If `c(i-1) < 9/4` (equivalently `frac(i log_2 3) > log_2(4/3)`, density `0.585`) and
HYP-9122(`i-1`) holds, then the path `(k_1+1, k_2, ..., k_(i-1), 0, 1)` descends, with multiplier `4c/9 < 1`. This is
the 1/2-transfer (loops note, Lemma 7.1, with `w = 2^(1-a)`, 3-adically), followed by `2 -> 1`.

**C4 (obstruction).** If a reverse path reaches the node `2` with multiplier `R_0 >= 3/2` (in particular after
the canonical transfer when `c(i-1) > 9/4`), then no continuation reaches total multiplier `< 1`.

*Proof.* A continuation node `y` satisfies `y = M(2 - beta) >= 1`. For `R_0 M < 1` we need `M < 2/3`, hence
`y < 4/3`, i.e. `y = 1`. Then `M = C/2 >= (4/3)/2 = 2/3`. ∎

Forward points. Let `X_i = -1 - 2^i/3^(A-1)`, `A = a0(i)`. Its absolute value `1 + 3 * 2^i/3^A` lies in `(3/2, 5/2)`. These
are the upper straddle partners of Theorem F.

**C5 (easy clocks).** If `3^A < 2^(i+2)` (equivalently `eta(i+1) < log_3 2`, density `0.631`) and HYP-9125(`i`)
holds, then the loop word followed by `H H` is legal at `X_i`. It ends `-4 -> -2 -> -1`, with multiplier
`3^A/2^(i+2) < 1` (mirror note, Lemma 4.1: `w(X_i) = -1 + (3^A/2^i)(X_i + 1) = -4`).

**C6 (obstruction).** If `3^A > 2^(i+2)`, no descent passes through that landing. By Lemma S at `v = -2` with
multiplier `3^A/2^(i+1)`, every continuation keeps multiplier `>= 3^A/2^(i+2) > 1`. ∎

So the canonical escapes decide exactly the easy clocks, and fail exactly where the census's slow descenders
live (`c` near 3). At the hard clocks every descent must leave the canonical landing, which is why the
certificates found so far start with first moves `k_1 = 33..1545`.

### 4.4 Sweep (FINITE-EXACT; `..._sweep_gen1.py`)

All generation-1 points of §4.3 with `i <= 300` were checked on both sides.

**Backward: `T_i^(-1)(2)`, `3 <= i <= 300` (298 points).**
* **Easy clocks (175 points).** All descend by the explicit route of C3. The `K0(i-1)`-loop is spliced from the record objects of §3.4, and the whole path is replayed in exact dyadic arithmetic.
* **Hard clocks (123 points).** For 122 of them a certificate was found: first move `k_1 <= 3000`, then the greedy map `G`, replayed exactly. The `k_1` range from 9 to 1545; the lanes' 11 lower-clock points (`e <= 200`) are among them, with the same certificates.
* **Hardest case, `i = 159`** (`c(158) = 2.9813`, value `1.50628`). No certificate exists with `k_1 <= 3000`. The deeper search (`gen1.py deep-back`, same route) finds `k_1 = 13097`: a 22,663-move path with `2^35920 < 3^22663`, a margin of `0.005` bits, replayed exactly.
* **So all 298 backward generation-1 upper partners with `3 <= i <= 300` descend.**

**Forward: `X_i = -1 - 2^i/3^(a0(i)-1)`, `10 <= i <= 300` (291 points).**
* **Easy clocks (184 points).** All descend by C5. The loop through `-1` is spliced from the record objects, followed by `H H`, and replayed exactly on the 2-adic rational.
* **Hard clocks (107 points).** 100 descend by the route "`M^a`, then the negative Collatz map" (`a <= j + 600`).
* **Remaining 7 (`i = 45, 64, 148, 213, 232, 278, 297`, all with `|x|` in `(1.5036, 1.5165)`).** These resist the plain route but descend when one extra E-only `x3`-move is allowed at an even value `< 10^7` in absolute value (`gen1.py deep-fwd`, `a <= j + 1500`, exact integer arithmetic). The certificates are `3^306 < 2^485`, `3^1224 < 2^1940`, `3^5879 < 2^9318`, `3^612 < 2^970` (twice), `3^1118 < 2^1772` and `3^7768 < 2^12312`.
* **So all 291 forward generation-1 upper partners with `10 <= i <= 300` descend.**

**Reading.**
* At the easy clocks the 3/2 wall is crossed canonically. Their density is `0.585` (backward) and `0.631` (forward), and the result holds for every `i` in the FINITE-EXACT ranges, by C3 and C5.
* At the hard clocks it takes a near-optimal climb with a long first move. The points that needed the deeper searches are exactly those closest to `3/2` (`|x| - 3/2 < 0.017`), where the needed carry precision is highest. The first move grows accordingly, up to `k_1 = 13097`.
* **No generation-1 point above `3/2` with `i <= 300` is hostile, on either side (589 points, all certified).**

## 5. HYP-9124 and HYP-9120: what the sweep adds

* **HYP-9124 (X_min).** Unchanged: OPEN, FINITE-EXACT to `10^18`.
  * New: X_min is implied by HYP-9126 (Proposition C1).
  * Theorem 6.2 of the endgame note (under HYP-9122, X_min ⟺ `S_1/2` ⟸ `X_T`) is now unconditional for every `m` whose thread precision is `<= 190536`.
* **HYP-9120.** OPEN. Its parts:
  * (a) is Q1; Collatz implies Q1.
  * (b) is Q2, implied by X_min.
  * (c) is structure. Its dimension clause cannot come from partial choice any more (§2), and it requires the HYP-9125 loops (§2.5).
  * The thread of `1` is settled to depth `190536` (1-escape at every precision `k <= 190536`).
  * The fragment "`-1-2^i/3^alpha(i)` hostile" now has `alpha(i) = a0(i)` for `10 <= i <= 176249`.

## 6. HYP-9123, 9127, 9132: exact placement, no new partial result

* **HYP-9123 (C2).**
  * **Placement.** The Periodicity Conjecture (every rational has an eventually periodic parity vector; Lagarias 1985 §2.8, Bernstein–Lagarias 1996) implies C2. C2 is its restriction to parity vectors of bounded discrepancy around a supercritical slope.
  * **C2 implies HYP-9127.** An eventually-`Y3` vector has width 1 at slope `9/10`.
  * **Coupled Z-number form.** By Proposition M of the hard-class note, C2 for integers is exactly this coupled Z-number statement.
  * **No overlap with Mahler.** It is logically independent of Mahler's Z-number problem, which concerns a different map.
  * This sweep proves nothing new here. The one mechanism beyond `Dio > eta` (q-difference Padé) needs a functional equation; generic strip words have none, and have `Dio = 1` almost surely.
* **HYP-9127.** Equivalent to `sum_k (2^10/3^9)^(k^3)` being irrational in `Q_2` (block identity, re-verified modulo `2^20000` in the cube-theta lane). It sits between C2 and the open real problem for non-unit cubic theta values. The proved conditional routes are HC-CT1 and `U(1/3)`. No new result.
* **HYP-9132.** Equivalent to the transcendence of `theta_3(rho)`, `rho = 2^10/3^9`, in `Q_2`.
  * Two results are proved: `theta_3(rho)` is irrational (Theorem Y), and `theta_3(rho)` and `psi(rho^2)` are not both algebraic (Corollary M, modulo BDGP 1996).
  * A p-adic analogue of Nesterenko's `tr.deg >= 3` would give the full statement; only Chudnovsky–Bertrand's `tr.deg >= 2` is known p-adically.
  * The Padé forms of Lemma P approximate `theta_3(rho)` only to exponent about `1.05`, far below the p-adic Roth threshold `2`, so Roth–Ridout cannot give transcendence either.

## 7. Reproduction

```
bash 04-computation/experiments/collatz_procgen_20260923_sweep_run.sh > 05-knowledge/results/collatz_procgen_20260923_sweep.out
```

The run script regenerates every record object with reconstruction (`FB=2`, `FC=3`, checkpoints every 64
layers). With `REUSE=1` it re-verifies the stored objects instead of recomputing them. The `.out` concatenates:
* the record logs (`objects/{pos,neg}_log.jsonl`);
* `reverify`, `tightcheck`, `hubcheck` and the spliced samples;
* the `E_S` counts and projections;
* `..._es_refute.py`;
* the generation-1 sweep and the deep searches;
* the `R_s` coverage table.

Measured on the mac-mini, with other lanes running:
* Q2-side records `<= 47468`: about 12 min, including reconstruction; `79335`: 16.5 min.
* Q1-side records `<= 50508`: 4 min.
* `E_S` DFS to `m = 64`: 2.6 CPU-min.
* Generation-1 sweep: about 12 min.

Peak memory `< 150 MB` per process. Transient checkpoint files are written under `scratch/procgen_sweep/ck`
(up to about 1 GB) and deleted after each object.

## 8. Proposed status changes (for the integrator; no HYP file edited)

* **HYP-9121: OPEN → REFUTED.**
  * Theorem R: `|Bad_m(E_{6 mod 8})| >= 2^(0.0536 m - 6)` for every `m >= 120`, so the lower box dimension of `Bad_inf(E_S)` is `>= 0.0536`.
  * The mechanism: the forward `E_S`-orbit of `-1` is the finite set `F` (16 values), whose cycles have rate `>= 2/3 > log_3 2`. So `alpha_S(b) >= 2b/3`, and the Theorem-P perturbations of `-1` cost only `(2/3^(2/3))^i`.
  * FINITE-EXACT: converged covering numbers `~2^(0.17 m)` on `m = 25..45`; raw counts to `m = 64`.
  * Suggested replacement question (OPEN): determine `dim_B Bad_inf(E_{6 mod 8})`. It is PROVED to lie in `[0.0536, h(log_3 2) = 0.95]`, and the data suggest about `0.17`.
* **HYP-9122: OPEN**, FINITE-EXACT range `s <= 6000` → **`s <= 190535`**. Record objects for every record `q <= 79335` are explicit and verified. Add Lemma T2 and the sparse-regime remark (§3.5).
* **HYP-9125: OPEN**, FINITE-EXACT range `K <= 4000` → **`K <= 176249`**. Record objects for every record `n <= 50508` are explicit and verified. The `E_S` analogue is false for all `K >= 35` (§2.4).
* **HYP-9126: OPEN.** Add:
  * `HYP-9126 ⟹ HYP-9124 ⟹ Q2` (C1);
  * the cost criterion;
  * C3 and C5, the easy clocks, PROVED given the loop hypotheses and hence unconditional for `i <= 190536` (backward) and `i <= 176249` (forward);
  * the obstructions C4 and C6;
  * FINITE-EXACT: every generation-1 upper partner with `i <= 300` descends (298 backward, 291 forward; §4.4).
* **HYP-9124: OPEN.** Add `⟸ HYP-9126`. Theorem 6.2 of the endgame note is now unconditional for thread precisions `<= 190536`.
* **HYP-9120: OPEN.** Add:
  * the route "HYP-9121 ⟹ dimension clause of (c)" is void, since HYP-9121 is false;
  * E's dimension-0 evidence rests on the HYP-9125 loops (§2.5);
  * the thread of `1` is settled to depth `190536`;
  * `alpha(i) = a0(i)` for `10 <= i <= 176249`.
* **HYP-9123, HYP-9127, HYP-9132: OPEN, unchanged.** Placement: `PC ⟹ C2 ⟹ HYP-9127`; HYP-9132 ⟸ p-adic Nesterenko.

Nothing in the sweep proves Collatz, E-SCC, Q1, Q2, X_min or PC, or brings any of them closer by a proved
reduction beyond those stated.
