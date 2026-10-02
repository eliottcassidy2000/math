# Syracuse drop multiplicities, a supplement to THM-4530: m(d) = O(log d / sqrt(log log d)), the Dirichlet series, 30-digit densities, and an independent cross-check

**Provenance.** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-01/02;
temporary identity, no `.machine-id`. This thread continued the "cdiff" lane listed in the wave-24 letter of
collatz-procgen-20260922 without knowing that the same lane was finishing on mac-mini. When this thread synced
with `origin/main` (05791bbd), [THM-4530](../../01-canon/theorems/THM-4530-syracuse-drop-multiplicities-two-copies-over-z-densities-and-pair-injectivity.md)
and its note [`procgen_cdiff_20261001_drop_multiplicities.md`](procgen_cdiff_20261001_drop_multiplicities.md)
were already there, independently audited. Most of this thread's derivation duplicates them, so only two things
are recorded here: what this thread adds, and the independent cross-check.
**Builds on:** THM-4530 (above), [THM-4527](../../01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md),
[HYP-9171](../hypotheses/HYP-9171-drop-multiplicity-maximal-order-and-pair-injectivity.md).
**Runner:** `04-computation/experiments/drop_multiplicities_supplement_20261002_run.py` (written from scratch, did not
read the THM-4530 code), stdout in `05-knowledge/results/drop_multiplicities_supplement_20261002.out`, 92 checks,
ending ALL CHECKS PASSED.
**Nothing here bears on the Collatz conjecture.** THM-4530 §4.7 explains why drop statistics cannot.

Notation as in THM-4530: `M_v = 2^v - 3`, `m(d) = #{v >= 2 : M_v | 6d+1}` for `d >= 0` (the number of odd `A >= 1`
with `(A - S(A))/2 = d`), `N(n) = #{v >= 3 : M_v | n}`, and `L_T = lcm{M_v : v in T}`.

## Status table

| # | Claim | Label |
|---|---|---|
| 1 | **Upper bound on the maximal order.** For `d >= 1`, with `n = 6d+1`, `k = m(d) - 1` and `tau` the divisor function: `k(k+3)/2 <= log2 n + (log2 n + 1) log2 tau(n)`. Hence `m(d) = O(log d / sqrt(log log d))`. This improves the upper side of THM-4530 Prop. 2.5(a), which is `O(log d)`, toward HYP-9171(a) | PROVED (Theorem 2) |
| 2 | Gap bound: `m(d) <= log2(1 + sqrt(6d+5))`, slightly sharper than `1 + (1/2) log2(6d+1)` | PROVED (Theorem 1); checked for every `d <= 10^8` |
| 3 | **Explicit tail.** For `k >= 2`, the density of `{d : m(d) >= k+1}` is at most `sum_{|S| = k} 1/L_S <= (4k+5)/(8 * 3^k)`. THM-4530 Lemma 2.2 gives `O(u^(-k))` for every `u < 4`, without a constant | PROVED (Lemma 3) |
| 4 | **Dirichlet series.** `F(s) = sum over n = 1 mod 6 of m'(n) n^(-s) = (L(s,chi_0) Q(s) + L(s,chi) Q^-(s))/2`, with `Q(s) = sum_{v>=2} M_v^(-s)` and `Q^-(s) = sum_{v>=2} (-1)^v M_v^(-s)`. `F` is meromorphic on C; in `Re s > -1` its only pole is `s = 1` (simple, residue `Q(1)/6`); `F(0) = 1/6` | PROVED (Theorem 3; standard facts about zeta and L(s, chi_-3) CITED); VERIFIED numerically |
| 5 | The densities `delta_1, ..., delta_11` of `{d : m(d) = k}` to 30 significant digits, from exact rationals at truncation `v <= 100` with error `< 7.9e-31` | PROVED enclosures |
| 6 | Independent cross-check of THM-4530: every printed digit of its densities `delta_1..delta_10`, its mean `1.34367343318176901854448...`, its record values for 2 to 11 copies, its "two copies over Z" count and its 5x+1 law (here to 12 digits, `P(m_5 = 0) = 0.592675448437`) agree with this thread's code | FINITE-EXACT (independent code; for the densities the same method family, plus method-independent checks) |
| 7 | Dictionary with the repo's cycle gates: the gate parameter `q(w) = \|Delta_w\|/gcd(c(w), \|Delta_w\|)` is the additive order of the word's drop class modulo `\|2^p - 3^k\|`, so `q = 1` means that the class is 0 | PROVED (restatement); checked for `k <= 6`, `p <= 2k+2` |
| 8 | Two small wording points: THM-4527 Statement item 2 lacked the range "m >= 0" (repaired in the same commit as this note); the upstream note's §2.5(c) says v = 19 and v = 16 enter "from k = 10 on", but v = 19 enters at k = 9 and v = 16 only at k = 12 | FINITE-EXACT (recomputed) |
| 9 | `max_{d <= X} m(d) = O(sqrt(log X))` (HYP-9171(a)) | OPEN. Data: the primitive part of `M_v` exceeds `2^(v/2)` for every `10 <= v <= 200` |

## 1. Three lemmas on the moduli M_v

**Lemma 1 (2-adic gap principle).** Let `n > 0` and `2 <= v < w` with `M_v | n` and `M_w | n`. Then
`n/M_v - n/M_w` is a positive multiple of `2^v`. In particular `n >= (2^v + 1) M_v`.

*Proof.* Put `a = n/M_v` and `b = n/M_w`. These are positive integers, and `a > b` because `M_v < M_w`. From
`a(2^v - 3) = n` we get `3a = -n (mod 2^v)`, and from `b(2^w - 3) = n` we get `3b = -n (mod 2^w)`, hence also
mod `2^v`. As 3 is invertible modulo `2^v`, `a = b (mod 2^v)`. So `a - b >= 2^v`, `a >= 2^v + 1`, and
`n = a M_v >= (2^v + 1) M_v`. QED

(Equivalently `gcd(M_v, M_w) = gcd(M_v, 2^(w-v) - 1)`; THM-4530 Prop. 2.3(a) has `gcd(M_v, M_w) | 2^(w-v) - 1`.)

**Lemma 2 (arithmetic progressions).** Let `Q = p^e > 1` be a prime power dividing some `M_v` with `v >= 2`, and let
`o = ord_Q(2)`. Then `{v >= 2 : Q | M_v} = {r + jo : j >= 0}` for some `r` with `3 <= r < o`. Moreover
`Q | 3 * 2^(o-r) - 1`, so `Q < 2^(o-r+2)`; and `Q < 2^r`. Hence `o > log2 Q + 1` and `o > 2 log2 Q - 2`.

*Proof.* `M_v` is prime to 6, so `o` is defined. `Q | M_v` iff `2^v = 3 (mod Q)`, so the set is a single residue
class modulo `o` intersected with `v >= 2`; let `r` be its least element. `r != 2` since `M_2 = 1`. If `r >= o`,
then `2^(r-o) = 3 (mod Q)` with `r - o >= 0`. Here `r - o = 0` would give `Q | 2`, `r - o = 1` would give `Q | 1`,
and `r - o >= 2` contradicts the minimality of `r`. So `r < o`. Then
`Q | 2^o - 1 - 2^(o-r) (2^r - 3) = 3 * 2^(o-r) - 1`, and `Q <= M_r < 2^r`. Multiplying the two bounds gives
`Q^2 < 2^(o+2)`. QED (Checked for all 110 distinct prime powers dividing some `M_v` with `v <= 64`; these form
150 pairs `(v, Q)`.)

**Lemma 3 (explicit tail).** For `k >= 2`, `sigma_k := sum over S subset {3, 4, 5, ...} with |S| = k of 1/L_S`
satisfies

```
sigma_k <= sum_{w >= k+1} C(w-3, k-2) (w+1) 2^(1-2w) = (4k+5)/(8 * 3^k).
```

Since `m(d) >= k+1` means that `L_S | 6d+1` for some `k`-set `S`, the density of `{d : m(d) >= k+1}` is at most
`sigma_k`.

*Proof.* Write `S = {v_1 < ... < v_k}`, `w = v_(k-1)` and `V = v_k`. By Lemma 1 applied to `n = L_S`,
`L_S >= max(M_V, (2^w + 1) M_w)`. If `w < V < 2w`, the term is at most `1/((2^w + 1) M_w) < 2^(1-2w)` (as `w >= 3`),
and there are `w - 1` such `V`. If `V >= 2w`, the term is at most `1/M_V <= 2^(1-V)`, and these terms add up to at
most `2^(2-2w)`. So the sum over `V` is at most `(w+1) 2^(1-2w)`. The elements `v_1, ..., v_(k-2)` lie in
`[3, w-1]`, which gives `C(w-3, k-2)` choices, and `w >= k+1`. The closed form follows from
`sum_j C(j,r) x^j = x^r/(1-x)^(r+1)` and `sum_j j C(j,r) x^j = x^r (r+x)/(1-x)^(r+2)` at `x = 1/4`, `r = k-2`,
`j = w-3`. QED

The true decay is like `2^(-k^2/2)` (§4), so the bound is weak; its point is the explicit constant.

## 2. The maximal order

**Theorem 1 (gap bound).** For every `d >= 1`, `m(d) <= log2(1 + sqrt(6d+5))`.

*Proof.* Let `n = 6d+1` and let `v_1 < ... < v_m` be the `v >= 2` with `M_v | n` (so `v_1 = 2`). If `m = 1` there is
nothing to prove. Otherwise Lemma 1 for `v_(m-1) < v_m` gives `n >= (2^(v_(m-1)) + 1)(2^(v_(m-1)) - 3)`. The `v_i`
are distinct integers `>= 2`, so `v_(m-1) >= m`, and `x -> (2^x + 1)(2^x - 3)` increases for `x >= 1`. So
`n >= (2^m + 1)(2^m - 3) = (2^m - 1)^2 - 4`, that is, `(2^m - 1)^2 <= 6d + 5`. QED
(Checked for every `d <= 10^8`. It is sharper than THM-4530 Prop. 2.5(a), `m(d) <= 1 + (1/2) log2(6d+1)`, for every
`d >= 1`, because `1 + sqrt(6d+5) < 2 sqrt(6d+1)` there; for large `d` the gain is about 1.)

**Theorem 2 (AP bound).** Let `d >= 1`, `n = 6d+1`, `k = m(d) - 1`, and let `tau(n)` be the number of divisors
of `n`. Then

```
k(k+3)/2  <=  log2 n + (log2 n + 1) * log2 tau(n).
```

Since `log2 tau(n) = O(log n / log log n)`, it follows that `m(d) = O(log d / sqrt(log log d))`.

*Proof.* Let `S = {v >= 3 : M_v | n}`, so `|S| = k`; if `k = 0` there is nothing to prove. Let `V = max S`. As
`2^(V-1) <= M_V <= n`, `V <= log2 n + 1`. For each `v`, `log2 M_v` is the sum of `log2 p` over the prime powers
`Q = p^e` (`e >= 1`) dividing `M_v`, and all of these divide `n`. So

```
sum_{v in S} log2 M_v = sum_{Q | n} c_Q log2 p,    c_Q = #{v in S : Q | M_v}.
```

By Lemma 2, the `v` counted by `c_Q` lie in one progression of difference `o_Q > log2 Q` inside `[3, V]`, so
`c_Q <= 1 + V/log2 Q`. For `Q = p^e`, `log2 p / log2 Q = 1/e`. Hence the sum is at most
`sum_{Q|n} log2 p + V sum_{p^a || n} (1 + 1/2 + ... + 1/a) <= log2 n + V log2 tau(n)`. Here we used
`1 + 1/2 + ... + 1/a <= log2(a + 1)`: it holds for `a = 1`, and the left side grows by `1/(a+1)` while the right
grows by `log2(1 + 1/(a+1)) > 1/(a+1)` (as `1/(a+1) <= 1/2`). Finally `prod (a+1) = tau(n)`.

On the other side, `S` consists of `k` distinct integers `>= 3`, and `log2 M_v >= v - 1` for `v >= 3`. So the sum
is at least `sum_{j=3}^{k+2} (j - 1) = k(k+3)/2`.

The divisor bound `log tau(n) = O(log n / log log n)` is classical (Wigert). An elementary version: for `y >= 2`,
each prime `p >= y` with `p^a || n` contributes `a + 1 <= 2^a <= p^(a ln 2/ln y)`, and each prime `p < y`
contributes at most `1 + log2 n`; so `ln tau(n) <= y ln(1 + log2 n) + (ln 2) ln n / ln y`, and `y = ln n/(ln ln n)^2`
gives the claim. QED

*Remarks.* Theorem 2 beats Theorem 1 only for astronomically large `d`; its point is that `m(d) = o(log d)`. It
does not reach HYP-9171(a), `max_{d <= X} m(d) = (1 + o(1)) sqrt(2 log2 X)`. That would follow from a lower bound
`log2 L_S >= c * sum_{v in S} v` for all finite `S`. The obstacle is the shared-prime "loss"
`log2(prod_{v in S} M_v / L_S)`. Lemmas 1 and 2 bound it, but for `S = [3, V]` both bounds still allow a loss of
order `V^2`. A primitive-divisor theorem would suffice: if the primitive part
`P_v = M_v / gcd(M_v, lcm_{u<v} M_u)` satisfies `P_v >= 2^(cv)`, then `L_S >= prod_{v in S} P_v`. The data
(EMPIRICAL): `P_v > 2^(v/2)` for every `10 <= v <= 200`, the minimum of `log2(P_v)/v` being 0.640 at `v = 19`. For
initial segments the loss is `0.24, 0.31, 0.40, 0.46` times `V log2 V` at `V = 20, 50, 100, 200`, far below `V^2`.

## 3. The Dirichlet series

Let `m'(n) = #{v >= 2 : M_v | n}` for `n = 1 (mod 6)`, so `m'(6d+1) = m(d)`. Put

```
F(s) = sum_{n = 1 mod 6} m'(n) n^(-s),   Q(s) = sum_{v>=2} M_v^(-s),   Q^-(s) = sum_{v>=2} (-1)^v M_v^(-s).
```

Let `chi_0` and `chi` be the principal and the non-principal character modulo 6, so that
`L(s, chi_0) = (1 - 2^-s)(1 - 3^-s) zeta(s)` and `L(s, chi) = (1 + 2^-s) L(s, chi_-3)`.

**Theorem 3.** For `Re s > 1`, `F(s) = (L(s,chi_0) Q(s) + L(s,chi) Q^-(s))/2`. The functions `Q` and `Q^-` extend
meromorphically to C:

```
Q(s)   = 1 + sum_{j>=0} ((s)_j / j!) 3^j 2^(-3(s+j)) / (1 - 2^-(s+j)),
Q^-(s) = 1 - sum_{j>=0} ((s)_j / j!) 3^j 2^(-3(s+j)) / (1 + 2^-(s+j)),
```

where `(s)_j = s(s+1)...(s+j-1)`. `Q` has simple poles exactly at `s = -j + 2 pi i k/ln 2` (`j >= 0`, `k` in Z), and
`Q^-` has simple poles exactly at `s = -j + (2k+1) pi i/ln 2`; `Q^-(0) = 1/2`. Consequently `F` is meromorphic on
C. In the half-plane `Re s > -1` its only pole is `s = 1`, which is simple with residue `Q(1)/6`, and `F(0) = 1/6`.

*Proof.* `m'(n) = sum_{v>=2} [M_v | n]`, so `F(s) = sum_v M_v^(-s) sum_{t >= 1, M_v t = 1 (mod 6)} t^(-s)`, absolutely
convergent for `Re s > 1`. For `t` prime to 6, `M_v t = 1 (mod 6)` iff `chi(t) = chi(M_v)`, and `chi(M_v) = (-1)^v`
(`M_v = 1 (mod 6)` for even `v`, `5 (mod 6)` for odd `v`). Since the sum of `t^(-s)` over `t` prime to 6 with
`chi(t) = e` is `(L(s,chi_0) + e L(s,chi))/2`, the formula follows.
Continuation: for `v >= 3`, `M_v^(-s) = 2^(-vs) (1 - 3 * 2^(-v))^(-s) = sum_j ((s)_j/j!) 3^j 2^(-v(s+j))`; summing the
geometric series in `v >= 3` gives the displayed formulas. For `j > |s| + 1` the `j`-th term is `O(j^|s| (3/8)^j)`
uniformly on compact sets, so both series are meromorphic, and their poles are the zeros of `1 - 2^-(s+j)`, resp.
`1 + 2^-(s+j)`, all simple. The residue of `Q` at `s_0 = -j + 2 pi i k/ln 2` is `((s_0)_j/j!) 3^j/ln 2`, which is
nonzero: for `k = 0` it is `(-3)^j/ln 2`, and for `k != 0`, `s_0` is not an integer. The same holds for `Q^-`, whose
poles are never real. At `s = 0` only the `j = 0` term of `Q^-` survives: `Q^-(0) = 1 - 1/2 = 1/2`.
Poles of `F`: in `L(s,chi_0) Q(s)`, the `j = 0` poles of `Q` (on the imaginary axis) are cancelled by the zeros of
`1 - 2^-s`; at `s = 0` the double zero of `(1 - 2^-s)(1 - 3^-s)` even leaves a zero. In `L(s,chi) Q^-(s)`, the
`j = 0` poles of `Q^-` are cancelled by the zeros of `1 + 2^-s`, and `L(s, chi)` is entire. All other poles of `Q`
and `Q^-` have real part `<= -1`. What remains in `Re s > -1` is the pole of `zeta` at `s = 1`, with residue
`(1/2)(1 - 1/2)(1 - 1/3) Q(1) = Q(1)/6`. Finally `F(0) = (1/2)(0 + L(0,chi) Q^-(0)) = (1/2)(2 * 1/3)(1/2) = 1/6`,
using `L(0, chi_-3) = 1/3` (CITED, standard). QED

On the line `Re s = -1` and further left there are genuine poles, at `s = -j + i pi l/ln 2` with `j >= 1`, except at
`s = -2, -4, ...` where `zeta` vanishes; this uses the standard fact that `zeta` and `L(s, chi_-3)` have no
non-real zeros with negative real part (CITED).
*Numerical checks (runner section D):* `F(2) = 1.04553040592...` and `F(3) = 1.00422049393706...` from the sieve
to `d = 10^8` with a rigorous tail bound agree with the closed form; the continuation agrees with the series at
`s = 1, 2`; the residues at `s = 0, -1`, `Q^-(0) = 1/2` and `F(0) = 1/6` are reproduced. `Q^-(1) = 0.85347503151184735...`.
*Remark.* The elementary count `sum_{d <= X} m(d) = sum_v (floor((X - e_v)/M_v) + 1)` (THM-4530 §3) gives a better
error term than any contour argument; Theorem 3 records the structure: the moduli enter only through `Q` and `Q^-`,
the residue classes mod 6 only through `L(s, chi_0)` and `L(s, chi)`.

## 4. Densities to 30 digits

The runner computes the truncated densities `delta_k^(100)` (truncation `v <= 100`) exactly, as rationals with
denominator `L_{[3,100]}` (1439 digits). Method: by the Chinese remainder theorem, for `n` uniform modulo the lcm,
the capped valuations at the 13 primes shared by two or more of `M_3, ..., M_100` (5, 11, 13, 19, 23, 29, 37, 47,
53, 61, 71, 97, 431) and the divisibility events for the pairwise coprime private parts are independent. So the
probability generating function of `N` factors: each private part contributes a linear factor, and each connected
component of the shared primes (two shared primes are joined when they divide a common `M_v`) contributes a sum,
over the capped valuation levels of its primes, of products of linear factors (4614 level configurations over
all components). The truncation error is `< 7.9e-31` (as in THM-4530 Thm 2.1(c)). Method-independent
checks: a brute-force count over a full period at truncation 7 (period `2874625 = 5^3 * 13 * 29 * 61`), the
inclusion-exclusion formula over all `2^14` subsets at truncation 16, the exact mean, and a sieve over `d <= 10^8`
(`{1: 69617491, 2: 26623796, 3: 3538961, 4: 213471, 5: 6192, 6: 89}`).

| k | `delta_k` (within one unit of the last digit) |
|---|---|
| 1 | 0.6961749034401325952680065255 |
| 2 | 0.2662380366777383738915135713 |
| 3 | 0.03538955883931588244282135104 |
| 4 | 0.00213458515130198837708089211 |
| 5 | 6.20614225702492134536025e-5 |
| 6 | 8.4914929170526775682603e-7 |
| 7 | 5.3028197674900875448e-9 |
| 8 | 1.680074489300358645e-11 |
| 9 | 2.866616892777891e-14 |
| 10 | 2.6972975032442e-17 |
| 11 | 1.436873267e-20 |

These agree with every digit THM-4530 prints (`delta_1 = 0.696174903440132595`, ..., `delta_10 = 2.70e-17`).
Both computations condition on the shared primes, so for the densities the agreement is between independent
codes, not independent methods; the full-period count and the sieve are the method-independent checks.

## 5. Independent cross-check of THM-4530 (FINITE-EXACT)

This thread's runner, written without reading the THM-4530 code, also reproduces (runner sections A, C, E, H):
* the drop identity for odd `A < 4 * 10^6` and the multiplicity theorem for `|d| <= 2 * 10^5`;
* Prop. 2.6 (two copies over Z) for `|d| <= 10^4`, by brute force over both signs;
* the mean `sum_{v>=2} 1/M_v = 1.34367343318176901854448283338124120618880718...` (enclosure of width `< 1e-61`);
* the least `d` with `m(d) >= k` for `k = 2, ..., 11` (exhaustive branch-and-bound over sets of moduli, sieve to
  `10^8` for `k <= 6`): 2, 24, 314, 7854, 479104, 121213354, 49576261854, 50617363353104, 115043252496711229,
  117459160799142164979, identical to THM-4530 Prop. 2.5(c);
* the 5x+1 law (THM-4530 §5.3, there to 5 digits). For `Syr_5(A) = (5A+1)/2^v` on odd `A >= 1`,
  `10d + 1 = (2^v - 5) Syr_5(A)`; `d = 0` never occurs; every `d <= -1` occurs exactly `1 + [d = 2 (mod 3)]`
  times (the unit `2^2 - 5 = -1` and the branch `2^1 - 5 = -3`); every `d >= 1` occurs
  `m_5(d) = #{v >= 3 : (2^v - 5) | 10d+1}` times. The densities of `m_5 = 0, 1, 2, 3, 4` on `d >= 1` are
  `0.592675448437, 0.327617684117, 0.072727956895, 0.006719211654, 0.000255113375`, each within `8.3e-25`
  (exact generating function at truncation `v <= 80`, 11 shared primes);
* the k-step identity, the unit words and the 3x-1 control (THM-4530 §§4.1, 4.4, 4.7), on finite ranges.

*Wording point in the upstream note §2.5(c).* It says the records take `v = 19` and `v = 16` "from k = 10 on".
Recomputing the sets `{v >= 3 : M_v | n_k}`: `n_9` uses `{3, ..., 9, 11, 19}`, `n_10` uses `{3, ..., 11, 19}`,
`n_11` uses `{3, ..., 12, 19}`, and `n_12` uses `{3, ..., 12, 16, 19}`. So `v = 19` enters at `k = 9` and `v = 16`
only at `k = 12`. The sentence about `n_12` there is correct.

## 6. Remark: drop classes and the repo's cycle gates

By THM-4530 Theorem 4.1, the points with valuation word `w` (length `k`, total valuation `p`, carry `c(w)`) have
k-step drops `D` in a half-progression with step `|Delta_w|`, `Delta_w = 2^p - 3^k`, inside the class
`D_w = -c(w) (2 * 3^k)^(-1) (mod |Delta_w|)`; here `2 * 3^k` is a unit because `Delta_w` is prime to 6. The repo's
minimal gate parameter is `q(w) = |Delta_w|/gcd(B(w), |Delta_w|)` with `B(w) = c(w)`
([`arithmetic_braids2_20260917_signed_cycles.md`](arithmetic_braids2_20260917_signed_cycles.md) §3, eq. (8)). So
`q(w)` is the additive order of `D_w` in `Z/|Delta_w|`, and the cycle gate `q(w) = 1` says `D_w = 0`: the word
can carry the drop 0, which is a cycle, subject to the sign condition of THM-4530 Prop. 4.2(a). For the map
`3A + b` the drop class of `w` is `b D_w`, so "the word is a cycle word at parameter `b` iff `q(w) | b`" reads
"`b D_w = 0`". (Checked for `k <= 6`, `p <= 2k + 2`, runner section I.)

## 7. THM-4527: one wording repair (applied)

THM-4527 (independently audited 2026-10-01, wording fixes applied) read, in Statement item 2: "So `m` occurs as a
descent exactly `#{v >= 2 : (2^v - 3) | 6m+1}` times, and every negative `m` occurs once, as an ascent." For `m < 0`
the displayed set is never empty (`v = 2` always counts; `m = -1` gives `{2, 3}`), while a negative `m` never occurs
as a descent. The commit that adds this note changes the sentence to "So every `m >= 0` occurs as a descent exactly
... times" and logs the slip in `01-canon/MISTAKES.md`. The title (with "quotient > 0") and the note's Theorem 2(b)
were already correct. The other repair this thread had found (the trunk labels `(4^j + 2)/6`) had already been
applied upstream.

## 8. Reproduction

`nice -n 10 python3 04-computation/experiments/drop_multiplicities_supplement_20261002_run.py > 05-knowledge/results/drop_multiplicities_supplement_20261002.out`
Single process, 44 s on one core; elapsed times go to stderr. Needs numpy, sympy and mpmath. Sections: A (drop identity, multiplicity
theorem), B (THM-4527 items 1-2), C (exact densities, Lemma 3), D (generating function, Dirichlet series), E (maximal
order: Lemmas 1-2, Theorems 1-2, records, primitive parts), F (label map), G (k-step drops), H (controls 3A-1,
5A+1, both signs), I (2-adic and 3-adic structure). Every check prints `PASS <claim>`.
