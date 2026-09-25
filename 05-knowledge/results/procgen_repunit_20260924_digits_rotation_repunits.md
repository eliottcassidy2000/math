# Digits, rotations and repunits: the consecutive-prime digit bias against Collatz's exact digit chains, circular primes as rotation orbits, and the repunits inside Collatz

**Status.**
* **PROVED** (hand proofs below; every statement is also checked by the scripts):
  1. **Rotation theorem.** Rotating the `k` base-`b` digits of `n` is multiplication by `b` modulo `b^k - 1`. Its fixed points are exactly the repdigits `c R_k`. A prime with `k >= 2` digits lies on a rotation orbit of exactly `k` numbers unless it is a repunit, which is fixed. So the repunit primes are exactly the multidigit prime fixed points of rotation, and `{337, 373, 733}` is a free orbit of size 3. Every repunit `R_k/(b^k - 1) = 1/(b-1)` is the same fixed point of the circle map `x -> bx`, seen at period `k`.
  2. **Exact Haar laws of Collatz's digit chains.** Under Haar measure, consecutive odd Syracuse terms `U(n)` form an exact Markov chain modulo 10, 3, 9 and 8, and the shortcut map `T` does too (explicit matrices in §2). Examples:
     * mod 3, every row is `(0, 1/3, 2/3)`: the classes are i.i.d.;
     * mod 10, `3 -> 5` has probability 1 and `9 -> 9` has probability `8/15`;
     * mod 9 the stationary law is `(8, 16, 11, 4, 2, 22)/63` on `(1, 2, 4, 5, 7, 8)`.
  3. **The Collatz repunit dictionary.**
     * The odd branch of `T` is the dilation `x + 1 -> (3/2)(x + 1)` about the fixed point `-1`, so the rising run of `n` has length `v_2(n+1)`.
     * `T^j(2^k - 1) = 3^j 2^(k-j) - 1` for `j <= k`, and `T^(k+1)(2^k - 1) = (3^k - 1)/2`.
     * Exact stopping times and heights for the trunk `(4^k-1)/3`, the minus-trunk `(2^k+1)/3` and the base-3 repunits.
     * The **shadowing lemma** for every even base `b`: `T^j(R_k^(b)) = T^j(xi_b) + 3^(a_j) b^k/((b-1) 2^j)` for `j <= k v_2(b)`, where `xi_b = 1/(1-b)`.
     * The pairing `tau(2^(2j) - 1) = tau(2^(2j-1) - 1) + 1`.
     * 2-adic continuity: the parity bits of `2^k - 1` beyond its forced prefix depend on `k` only 2-adically.
  4. **The 2-adic macrocosm.**
     * `xi_2 = -1` is fixed.
     * `xi_4 = -1/3 -> 0`.
     * `xi_(-2) = 1/3 -> 1 -> 2`.
     * `xi_10 = -1/9 -> 1/3 -> 1`.
     * `xi_6 = -1/5 -> 1/5 -> 4/5 -> 2/5 -> 1/5`.
     * `T` acts on the base by the Möbius map `M(b) = -(b+2)/(b-4)`.
  5. **The Bernstein–Lagarias conjugacy.** Its known odd fixed points `-1` and `1/3` are the 2-adic repunits `...1111` in bases `2` and `-2`. Its known 2-cycle `{1, -1/3}` joins the two ends of the base-4 trunk. A **second odd 2-cycle `{-1/5, 5/7}`** exists: `Q(-1/5) = 5/7` and `Q(5/7) = -1/5`, where `Q = Phi^(-1)` is the parity-vector map. Its point `-1/5` is the base-6 repunit limit. Bernstein–Lagarias 1996 knew only `{1, -1/3}`; whether later work records the new cycle is UNVERIFIED.
* **FINITE-EXACT:**
  * **Lemke Oliver–Soundararajan tables reproduced exactly:** the first `10^6` primes `> 3` mod 3, the first `10^8` primes `>= 11` mod 10, the four exact mod-8 counts at `10^11`, and 34 four-digit table entries.
  * Consecutive-prime transitions to `10^11` (`4,118,054,813` primes).
  * Circular primes complete below `10^9` (sieve) and for every length `<= 16` (necklace search): exactly OEIS A068652/A016114.
  * Permutable primes below `10^9` are exactly A003459, and every multiset of 4 to 40 digits other than the all-ones one has a composite permutation.
  * `R_n` is prime for `n <= 1100` exactly when `n = 2, 19, 23, 317, 1031` (BPSW and Miller–Rabin agree).
  * Among odd rationals `p/q` with `|p| <= 201` and odd `q < 100`, the fixed points of the Bernstein–Lagarias map `Q` are exactly `-1` and `1/3`, and its 2-cycles are exactly `{1, -1/3}` and `{-1/5, 5/7}`.
  * Integral Collatz cycles with a period `<= 18` are the 5 known ones.
  * **Real orbits follow the exact Haar digit law.** Beyond the random bits, deterministic orbits above `2^40` agree with Theorem H over `3*10^9` transitions. The apparent deviations below `2^40` and in full-orbit tallies come from orbits merging: counting each integer once removes them.
* **CITED** (read): Lemke Oliver–Soundararajan (arXiv:1603.03720v4 full text; PNAS 113(31) (2016) E4446–E4454 via Europe PMC); OEIS A016114, A068652, A293663, A003459, A004023, A000043, A000978, A028491; P. De Geest's circular-prime table; Wikipedia "Circular prime" and "Permutable prime"; Bernstein–Lagarias, Canad. J. Math. 48 (1996) 1154–1169; Lagarias's annotated bibliographies I and II. **Not read:** Richert 1951 and Johnson 1977 (quoted via Wikipedia); Hotzel 2003.
* **VERDICTS on the owner's text:**
  * "The odds from one of `{1,3,7,9}` to the next are not uniform" is **REAL** (Lemke Oliver–Soundararajan, reproduced), but it is a **vanishing** second-order effect.
  * "That generates the finite family `{337,373,733}`" is **NOT SUPPORTED**. Circular primes are a rotation phenomenon; the consecutive-prime bias plays no part in them.
  * "Brouwer: fixed point versus orbit" is **REAL** (PROVED).
  * "The length of all-one primes is hidden structure" is, for Collatz, **NUMEROLOGY** (tested). What is real is that the length of the all-ones binary *tail* is the length of the rising run.
  * The Collatz analogue of the digit bias is **exact and first-order**.
* **CORRECTED leads:**
  * `-1/3` is *not* a fixed point of `T` (`T(-1/3) = 0`). It is the fixed point of the ladder map `p -> 4p+1` and half of a Bernstein–Lagarias 2-cycle.
  * The minus-trunk numbers are *all* `(2^k+1)/3` with `k` odd (`1, 3, 11, 43, 171, 683, ...`). The list `3, 11, 43, 683` gives only the prime ones (Wagstaff primes).
* **OPEN:**
  * Lemke Oliver–Soundararajan's Main Conjecture.
  * Infinitely many repunit primes.
  * Finitely many non-repunit circular primes.
  * Bernstein–Lagarias's Fixed Point Conjecture.
  * Collatz.

Session `collatz-procgen-20260922` (mac-mini), lane `procgen_repunit_20260924`, 2026-09-24. No HYP or THM file was created.

* **Scripts:** `04-computation/experiments/procgen_repunit_20260924_{primes,collatz,repunits,dictionary}.py`, with C helpers `..._{primes,collatz,circular,collatz_repunits}.c` (the last needs GMP). The runner is `..._run.py`.
* **Output:** [`procgen_repunit_20260924.out`](procgen_repunit_20260924.out). It takes about 10 minutes and peaks below 700 MB, running one process at a time, and every check raises on failure.

## 0. The request, the leads, and what became of each

The owner wrote: "Think Brouwer's fixed point theorem, and how in base 10 primes must end in one of {1,3,7,9}. When put in numerical order, the odds of going from any one of those 4 nodes to the next are not uniform, and that generates the finite family of numbers like {337,373,733} which are prime as you cycle through their digits, along with the infinite family of repunit primes." He added: "Think of the length of all-one primes being important hidden structure, micro/macrocosm."

| lead | verdict | where |
|---|---|---|
| (a) the digit bias of consecutive primes is Lemke Oliver–Soundararajan's, a Hardy–Littlewood second-order term of size `log log x/log x` | **TRUE, CITED and reproduced exactly**. The paper's Main Conjecture has `c1 = 1/2 - (phi(q)/2)[a=b]` plus a constant `c2/log x`. At `10^11` the two-term prediction explains 86% of the diagonal deficit | §1 |
| (b) rotation is multiplication by `b` mod `b^k - 1`; its fixed points are the repdigits; repunit primes are the multidigit prime fixed points, and `{337,373,733}` is a prime orbit | **TRUE, PROVED**. Sharpened: every prime with `k >= 2` digits other than a repunit has an orbit of exactly `k` elements | §3.1 |
| (c) base 2: `T^k(2^k-1) = 3^k - 1`, and `...1111 = -1` is the hostile fixed point | **TRUE, PROVED**, with the exact run lemma (`v_2(n+1)` odd steps for every `n`) | §4.1 |
| (c) base 4: the trunk `(4^k-1)/3` tends to `-1/3 = E(0)`, the limit of every sibling ladder | **TRUE**. The ladder map `p -> 4p+1` *is* the base-4 repunit generator | §4.2 |
| (c) "the repunits of bases 2 and 4 are exactly Collatz's 2-adic fixed points `-1` and `-1/3`" | **CORRECTED**. `-1` is `T`-fixed, but `T(-1/3) = 0`. The fixed-point reading is right one level up: `-1` and `1/3` (bases 2 and `-2`) are the known odd fixed points of the Bernstein–Lagarias conjugacy, and `-1/3` (base 4) lies on its known 2-cycle | §4.4, §4.5 |
| (c) base `-2`: `(2^k+1)/3`, `k` odd, are the Jacobsthal/Wagstaff minus-trunk numbers `3, 11, 43, 683, ...` | **TRUE with a correction**: all odd `k` give minus-trunk numbers (`1, 3, 11, 43, 171, 683, ...`), `171 = 3^2*19` included; the listed ones are the primes. The even-length base `-2` repunits `-(4^j-1)/3` are the negative trunk of `3n-1` | §4.2 |
| "the dictionary: which bases and signs give which Collatz objects" | **COMPLETE for even `b`**: shadowing lemma plus the orbit of `xi_b` for every even `|b| <= 64` on both sheets | §4.3, §4.4 |

## 1. Consecutive primes: the last-digit and mod-3 bias (CITED; FINITE-EXACT to `10^11`)

**1.1 Exact reproduction.** A segmented sieve records every consecutive pair `(p_n, p_(n+1))` mod 120, so every modulus dividing 120 is an aggregation. Up to `10^11` there are `4,118,054,813` primes.
* The first `10^6` primes greater than 3 (pairs `n = 3..10^6+2`), taken mod 3, give `215873, 283957, 283957, 216213`. This is the paper's p. 2 table exactly.
* The first `10^8` primes, taken mod 10, give the paper's 16-entry table exactly. The window is `n = 5..10^8+4`, i.e. the primes from 11 on; with the window `n = 4..10^8+3` two cells move by one.
* At `10^11` mod 8, counting pairs with both primes `<= x` gives `278676326, 278696997, 278692843, 278681776`: the Conjecture 1.6 example exactly. Our convention `p_n <= x` adds the boundary pair `(99999999977, 100000000003) = (1,3) mod 8`.
* All 34 four-significant-digit "Actual" entries of the paper's Section 5 (`q = 3, 4, 8` at `10^9, 10^10, 10^11`; `q = 5` at `10^9`) agree.

**1.2 The transition matrices** `P(next last digit = b | current = a)`:

| `x` | `1->1` | `1->3` | `1->7` | `1->9` | `9->1` | `9->9` | `3->3` | mod 3: `1->1` | max row TV from uniform (mod 10) |
|---|---|---|---|---|---|---|---|---|---|
| `10^6` | 0.1638 | 0.3211 | 0.3339 | 0.1812 | 0.3546 | 0.1610 | 0.1430 | 0.4179 | 0.155 |
| `10^8` | 0.1770 | 0.3039 | 0.3103 | 0.2087 | 0.3308 | 0.1768 | 0.1667 | 0.4388 | 0.114 |
| `10^9` | 0.1832 | 0.2986 | 0.3023 | 0.2159 | 0.3220 | 0.1832 | 0.1754 | 0.4451 | 0.101 |
| `10^10` | 0.1883 | 0.2943 | 0.2962 | 0.2212 | 0.3152 | 0.1883 | 0.1823 | 0.4500 | 0.091 |
| `10^11` | 0.1927 | 0.2908 | 0.2914 | 0.2251 | 0.3096 | 0.1927 | 0.1879 | 0.4540 | 0.082 |

* The diagonal is suppressed.
* `9 -> 1` is the most favoured transition and `1 -> 9` the least, exactly the extremes predicted by the `c2` constants (`+3.93` and `-3.70`).
* The ordering of the 12 off-diagonal cells matches the prediction with Spearman correlation 0.94.
* The stationary law of the empirical matrix, which is the Brouwer/Perron fixed point of this four-node chain, is uniform to `2*10^-6` at `10^11`, as Dirichlet requires. **The bias lives in the transitions, never in the stationary law.**

**1.3 The decay against the predicted second-order term.** Write `delta(a,b;x) = phi(q)^2 pi(x;q,(a,b))/li(x) - 1`. The prediction is `c1 LL/L + c2/L`, with `L = log x` and `LL = log log x`.
* For `q = 3`, `c2 = -/+ (1/2) log(2 pi/3)`: formula (1.1) and the `q`-prime case of (2.23).
* For the last digit, `c2` is taken from the explicit `q = 5` formula of Section 5, since primes `> 5` mod 10 are primes mod 5. The constants involved are `L(0,chi) = (3+i)/5`, `L(1,chi) = 0.8648+0.2042i` and `A_(5,chi) = 1.8910+1.5590i`.
* Checks: every row and column of `c2` sums to 0; relation (1.2) holds; the symmetry `c2(a,b) = c2(-b,-a)` holds; `A_(12,chi) = 1.0356` against the paper's "`~ 1.036`"; and the paper's (5.1) predictions are reproduced to four digits.

| `x` | `LL/L` | mod 10 diagonal: actual | predicted | ratio | mod 3 diagonal: actual | predicted | ratio |
|---|---|---|---|---|---|---|---|
| `10^6` | 0.190 | -0.387 | -0.310 | 1.25 | -0.165 | -0.122 | 1.36 |
| `10^8` | 0.158 | -0.313 | -0.256 | 1.22 | -0.122 | -0.099 | 1.23 |
| `10^10` | 0.136 | -0.259 | -0.219 | 1.18 | -0.100 | -0.084 | 1.19 |
| `10^11` | 0.128 | -0.239 | -0.205 | 1.17 | -0.092 | -0.078 | 1.17 |

* **What the table shows.** The bias decays, the two-term formula has the right shape, and the ratio actual/predicted falls monotonically toward 1: from 1.25 to 1.17 (mod 10) and from 1.36 to 1.17 (mod 3).
* **The residual.** It is close to a constant times `(LL/L)^2`: `-2.1 (LL/L)^2` on the mod-10 diagonal and `-0.85 (LL/L)^2` mod 3. These are the "non-oscillating lower order terms of size `(log log x/ log x)^2`" that the paper anticipates. The paper also reports that its unsimplified integral (2.16) fits better than the two-term formula; that integral was not recomputed here.
* **Summary measures of the bias.**
  * The largest row distance from uniform falls from 0.155 to 0.082 (mod 10) and from 0.082 to 0.046 (mod 3) between `10^6` and `10^11`.
  * The second eigenvalue of the empirical matrix falls from 0.17 to 0.098 (mod 10) and from 0.16 to 0.092 (mod 3).
  * **The prime chain forgets its state faster as `x` grows, and in the limit it is memoryless and uniform.**

## 2. The Collatz analogue: exact, first-order digit chains (PROVED; FINITE-EXACT)

`U(n)` = the odd part of `3n+1` (consecutive odd terms); `T(n) = n/2` or `(3n+1)/2`.

**2.1 Theorem H (exact Haar laws).** Let `n` be Haar-random in `Z_2 x Z_m`, `m` odd, with `n` odd for `U`. Then `v = v_2(3n+1)` is geometric, `P(v = j) = 2^-j`, independent of `n mod m`, and `U(n) = (3n+1) 2^-v (mod m)`. Consequently:

* **mod 3:** `U(n) = 2^-v = (-1)^v (mod 3)`, so **every row is `(0, 1/3, 2/3)`**. The mod-3 classes of consecutive odd terms are i.i.d., and 0 never occurs.
* **mod 9:** `U(n) = (3(n mod 3) + 1) 5^v (mod 9)` (`5 = 2^-1`). It depends on `n` only mod 3 and on `v` mod 6, where `P(v = r mod 6) = 2^(6-r)/63`. The matrix, in 63rds:

  | from `n mod 3` | to 1 | 2 | 4 | 5 | 7 | 8 |
  |---|---|---|---|---|---|---|
  | 0 | 1 | 2 | 4 | 32 | 16 | 8 |
  | 1 | 16 | 32 | 1 | 8 | 4 | 2 |
  | 2 | 4 | 8 | 16 | 2 | 1 | 32 |

  The stationary law is `(8, 16, 11, 4, 2, 22)/63` on the classes `(1, 2, 4, 5, 7, 8)`. It is *not* uniform: `8 = -1 mod 9` gets 35%, `7` gets 3%. `P^2` has identical rows, so the chain forgets its state in two steps.
* **mod 10 (last digit):**

  | from | to 1 | 3 | 5 | 7 | 9 |
  |---|---|---|---|---|---|
  | 1 | 4/15 | 2/15 | 0 | 8/15 | 1/15 |
  | 3 | 0 | 0 | **1** | 0 | 0 |
  | 5 | 1/15 | 8/15 | 0 | 2/15 | 4/15 |
  | 7 | 8/15 | 4/15 | 0 | 1/15 | 2/15 |
  | 9 | 2/15 | 1/15 | 0 | 4/15 | **8/15** |

  * `3 -> 5` is forced, and `5` is entered only from `3`.
  * The stationary law is uniform, `1/5` each, **including the digit 5**.
  * The characteristic polynomial is `(lambda - 1)(15 lambda^4 + 2 lambda^3 - 8 lambda^2 + 1)/15`, with `|lambda_2| = 0.7007`.
* **mod 8:** rows 1 and 5 are uniform; `3 -> {1,5}` and `7 -> {3,7}` with `1/2` each; `P^2` is uniform.
* **`T` mod `M`:** `T(x) mod M` depends on `x mod 2M`, and the two lifts are equally likely, so every row has two entries `1/2`.
  * mod 10: `d` goes to the two digits `= 3d (mod 5)` if `d` is even, and `= 3(3d+1) (mod 5)` if `d` is odd; the next parity bit is fresh. `|lambda_2| = (1+sqrt5)/4 = cos 36°`.
  * mod 3: `0 -> {0,2}`, `1 -> 2`, `2 -> {1,2}`, with stationary law `(0, 1/3, 2/3)`.
  * mod 9: the stationary law equals `U`'s `(8,16,11,4,2,22)/63`, because the backward gaps between odd steps are again `Geom(1/2)`.
  * mod 8: the window chain of 3 i.i.d. bits, mixed exactly in 3 steps.
* **Markov property.** `(U^t n mod M)_t` and `(T^t n mod M)_t` are Markov chains with these matrices. For `T` the state mod `2^k` is the window of the next `k` i.i.d. parity bits. For `U` it is the induced chain on windows that start with an odd bit, and the odd-prime part is updated by a fresh `v`.

*Proof.* Haar measure on `Z_2 x Z_m` is a product, and by Terras the parity vector of a Haar-random 2-adic integer is i.i.d. fair. So `v` is geometric and independent of the `m`-adic coordinate, and `U(n) mod m` is the stated function of `(n mod m, v)`. Given `v = j`, `n` lies in one class mod `2^(j+1)`, and `U(n) = (3n+1)/2^j` maps that class affinely onto the odd 2-adic integers, so `U(n) mod 2^k` is uniform on the odd classes. The matrices follow by summing geometric series over `v mod ord_m(2)`. A second code path re-derives the `U` laws by brute force over all odd `n mod 2^20 m`, agreeing to within the undetermined mass (`<= 3*10^-5`). ∎

**2.2 Empirical confirmation.** The starts are random in `[2^100, 2^101)`, whose low 100 bits are uniform.

The sample is `10^7` starts, each followed for at most 200 `U`-steps and 500 `T`-steps, stopping below `10^6`. The transitions fall into three regimes.

| regime | `U`-transitions | `T`-transitions | what it tests | result (chi-square p-values, mod 10 / 3 / 9 / 8 / 72) |
|---|---|---|---|---|
| **Terras**: determined by the 100 random low bits | `3.1*10^8` | `9.1*10^8` | Theorem H itself (exact) | `U`: 0.12 / 0.021 / 0.27 / 0.87 / 0.22; `T`: 0.93 / 0.99 / 0.65 / 0.73 / 0.83 |
| **post-high**: deterministic orbit, value `>= 2^40` | `1.2*10^9` | `2.1*10^9` | the heuristic that real orbits follow the Haar law | `U`: 0.077 / 0.58 / 0.82 / 0.95 / 0.91; `T`: 0.75 / 0.37 / 0.32 / 0.75 / 0.12 |
| **post-low**: `10^6 <= ` value `< 2^40` | `3.1*10^8` | `8.7*10^8` | the same, where the `10^7` orbits funnel through few integers | `U`: 0.004 / 0.19 / 0.043 / 0.50 / `< 10^-3`; `T`: 0.12 / 0.096 / 0.42 / 0.41 / `< 10^-3` |

**The deviation in the low band is a sampling effect: orbits merge.** Take the `U`-transitions from odd values in `[10^6, 10^8)`. There are `7.6*10^7` visits to only `1.87*10^7` distinct integers, a mean multiplicity of 4.1.
* Counted per visit, they deviate: mod 72, chi-square is 532 on 168 degrees of freedom; mod 10, `p < 10^-3`.
* Counted once per integer, they agree with Theorem H: p = 0.81, 0.28 and 0.998 mod 10, 3 and 9, and chi-square 47 on 168 mod 72.

The one-step law does not change. Merged paths are counted with multiplicity, which is a tree-weighted sample, exactly as in the full-orbit tallies below.

**No deviation from Theorem H is detected in deterministic orbits above `2^40`**, over `3*10^9` transitions: the Haar law describes real orbits there, as the heuristic assumes.

* **Terras regime.** These transitions are determined by the random bits, so they are an **exact** test of Theorem H.
  * `3.1*10^8` `U`-steps and `9.1*10^8` `T`-steps; the chi-square p-values are in the table above (mod 5 adds 0.12 for `U` and 0.52 for `T`).
  * There are **zero** counts in cells of probability 0.
  * An independent histogram of `v` over `3.1*10^9` Terras-regime steps gives `P(v even) = 0.3333215` (exact `1/3`, `z = -1.40`) and `P(v=1..6) = 0.500013, 0.249989, 0.125000, 0.062500, 0.031245, 0.015627`.
  * The two-step Markov checks give p = 0.11, 0.94, 0.85, 0.71 and 0.73.
  * The bulk sample, as row frequencies times 15, reads `1: 4.000 1.999 0 8.001 1.000`, `9: 2.000 1.001 0 4.000 7.999`, and `3 -> 5` in 100% of `6.2*10^7` cases.
* **Full orbits** of every `n <= 10^7`, followed to 1. The tallies weight each integer by the number of starts whose orbit passes through it; orbits merge, and "highways" such as the trunk into 1 are counted millions of times. These tallies are therefore far from Haar (chi-square in the millions), even when restricted to values `> 10^4`.
  * Example: `1 -> 9` falls from `1/15` to `0.82/15`, and `5 -> 1` rises from `1/15` to `2.5/15`, because the trunk `(4^k-1)/3 -> 1` ends in `...1 -> 1` or `...5 -> 1`.
  * Even so, **no transition ever lands in a probability-0 cell**, in any regime: the support of the law is exact for every integer, since `U(n) mod M` is always some legal branch.

**2.3 The 10-adic reading of the forced cells.**
* `9 -> 9` has probability `8/15`: when `v = 1 mod 4`, `n = -1 mod 5` stays `-1`. This is the mod-5 shadow of `T(-1) = -1`, and `-1 = ...9999` in `Z_10`: the hostile fixed point is the 10-adic repdigit 9.
* `3 -> 5` has probability 1: `n = 3 = -1/3 (mod 5)` makes `3n + 1 = 0 (mod 5)`. This is the shadow of `3(-1/3) + 1 = 0`, i.e. of `T(-1/3) = 0`, the base-4 repunit limit.
* So the rigid cells of the Collatz last-digit chain are the 5-adic shadows of the same two 2-adic points that §4 finds as repunit limits.

**2.4 Why Collatz's bias is exact and first-order while the primes' vanishes.** In both cases the next residue is `f(current residue, gap)`.
* **Primes.** The gap `g = p_(n+1) - p_n` has scale `log x` (Poisson/Gallagher), so `g mod q` equidistributes as `x -> infinity` and the matrix tends to the uniform rank-one matrix.
  * The departures come from the *correlations* between the events "`p + h` is prime" for different `h`, encoded in the singular series `S(H)`.
  * Montgomery–Soundararajan's `-(1/2) log H` correction to sums of `S_0({0,h})` gives the `c1 log log x/log x` term; the arithmetic of short gaps gives `c2/log x`.
  * These are sieve correlations that average out, so the bias is second order and vanishes.
* **Collatz.** The "gap" is `v = v_2(3n+1)`, and under Haar measure its law is **exactly** `Geom(1/2)` at every scale. Parity bits are i.i.d. (Terras) and independent of every odd-prime residue, because `Z_hat` is a product.
  * The multiplier `2^-v` therefore has a fixed, non-uniform law on the cyclic group `<2>` in `(Z/m)^*`, and it never spreads.
  * The residue map `n -> 3n+1` is deterministic, and not a bijection mod 9, since `3n+1 mod 9` sees only `n mod 3`.
  * So the transition law is first order (row distance from stationarity 0.8 mod 10), exact, and **independent of the size of `n`**.
* **In one line.** Primes have weak correlations with a growing scale; Collatz has exact independence with a fixed scale. The prime chain tends to "memoryless and uniform"; the Collatz chain is the same fixed matrix forever.

## 3. Rotation, circular primes, permutable primes, repunit primes

**3.1 Theorem R (PROVED; checked on 36,000 random `(b, k, n)` with `b in {2,3,4,10,12,16}`, `k <= 30`).** Let `b >= 2`, `0 <= n <= b^k - 1`, and let `rho(n)` move the leading digit `d` of the `k`-digit base-`b` string of `n` (leading zeros allowed) to the end.
1. `rho(n) = bn - d(b^k - 1)`. So `rho(n) = bn (mod b^k - 1)`, and `rho` is multiplication by `b` on `Z/(b^k - 1)` (the all-`(b-1)` string represents 0).
2. `rho(n) = n` iff `n = d R_k`, i.e. `n` is a repdigit. There are `b` such strings, and `b - 1` classes in `Z/(b^k-1)`.
3. Let `n` be a prime with exactly `k >= 2` digits. Its rotation orbit has exactly `k` elements unless `n = R_k`, which is fixed. In particular **the repunit primes are exactly the multidigit prime fixed points of rotation**, and `{337, 373, 733}`, `{113, 131, 311}` and `{1193, 1931, 9311, 3119}` are free orbits.
4. A circular prime with `k >= 2` digits uses only `1, 3, 7, 9`: a digit `0, 2, 4, 5, 6, 8` rotated into the last place gives a multiple of 2 or 5.

*Proof.*
* Item 1: `bn = d b^k + (n - d b^(k-1)) b = rho(n) + d(b^k - 1)`.
* Item 2: `rho(n) = n` iff `(b-1)n = d(b^k-1)`.
* Item 3: if `rho^j(n) = n` with `j < k` minimal, then `j | k`, and the string is `j`-periodic with block `m`. So `n = m (b^k-1)/(b^j-1)`, and the factor `F = 1 + b^j + ... + b^(k-j)` satisfies `1 < F`. Primality forces `m = 1`. The block `0...01` then gives a string with leading digit 0 unless `j = 1`, and `j = 1` means `n = R_k`. ∎

Checked: among the primes `10 < p < 10^7`, the only one whose orbit is smaller than its number of digits is `11`.

**The Brouwer/Lefschetz reading.**
* The circle map `x -> bx mod 1` has degree `b` and Lefschetz number `1 - b`. Its `b - 1` fixed points `c/(b-1)` are the repdigit expansions `0.ccc..._b`.
* Its points of period dividing `k` are `n/(b^k - 1)`, whose expansions repeat the `k`-digit block of `n`, and on these points the map *is* rotation.
* **Every repunit is the same fixed point**: `R_k/(b^k-1) = 1/(b-1) = 0.111..._b` for every `k`. The finite repunits are that one fixed point written at period `k`: this is the exact "microcosm = macrocosm" statement.
* The lengths are organised by divisibility: `R_d | R_k` for `d | k`, and `gcd(R_m, R_n) = R_gcd(m,n)`. Hence `R_k` prime forces `k` prime.
* For prime `k` the `b^k` strings split into the `b` fixed repdigits and `(b^k - b)/k` free orbits: Fermat's little theorem in necklace form.
* This is the only sense in which the "length of all-one primes" is structure: a prime length is forced, and which prime lengths succeed is, as far as anyone knows, random (§3.5, §5).

**3.2 Circular primes (FINITE-EXACT; C helper `..._circular.c`).** There are two independent code paths.
* A bitset sieve to `10^9` checks every rotation of every prime.
* An FKM necklace enumeration visits every rotation orbit of strings over `{1,3,7,9}` of length `<= 16`, skipping orbits with digit sum `= 0 mod 3`. Each orbit is tested by trial division and a base-2 strong-pseudoprime filter, and survivors are certified by a deterministic 12-base Miller–Rabin.

The two paths agree on every length `<= 9`.
* **The circular primes below `10^9`** are exactly 55 numbers, all below `10^6`: `2 3 5 7 11 13 17 31 37 71 73 79 97 113 131 197 199 311 337 373 719 733 919 971 991 1193 1931 3119 3779 7793 7937 9311 9377 11939 19391 19937 37199 39119 71993 91193 93719 93911 99371 193939 199933 319993 331999 391939 393919 919393 933199 939193 939391 993319 999331`.
  * This is OEIS A068652, whose JSON data field ends at 319993.
  * The orbit representatives are exactly A016114 below `10^9`: `2, 3, 5, 7, 11, 13, 17, 37, 79, 113, 197, 199, 337, 1193, 3779, 11939, 19937, 193939, 199933`.
  * The 54 non-repunit ones match A293663's comment ("54 terms that have 6 or fewer decimal digits, the largest of which is 999331").
* **No circular orbit of length 7 to 16 exists** (268 million length-16 necklaces). So the list is complete below `10^16`. De Geest's table (read) reports none of lengths 7–25 other than `R19` and `R23`; lengths 24 and 25 are due to S. Castelli (2024). Lengths 26 and 27 are marked "unknown" there.
* **Eligible primes** (all digits in `{1,3,7,9}`) number 757, 2709, 9177 and 33191 for lengths 6–9, equal to De Geest's counts.
* **Near misses** (orbits with exactly one composite rotation) exist only at lengths 2–9 and 12 within `<= 16`:
  * two at length 7, where `3999131 = 17*235243` and `7919777 = 83*95419` are the single composite rotations;
  * one at length 8, with `37177739 = 29*683*1877`;
  * one at length 9, with `391331191 = 29*131*239*431`;
  * one at length 12, with `373379311139 = 23*53*306299681`: eleven of its twelve rotations are prime.

  The near misses of lengths 7–12 agree with De Geest's "near miss" entries.

**3.3 Permutable primes (FINITE-EXACT).**
* Below `10^9` they are exactly A003459's 22 terms: `2, 3, 5, 7, 11, 13, 17, 31, 37, 71, 73, 79, 97, 113, 131, 199, 311, 337, 373, 733, 919, 991`.
* Every one of the 135,716 digit multisets of length 4–40 over `{1,3,7,9}`, except the all-ones multiset, has an explicit composite permutation. So the permutable primes with 4 to 40 digits are exactly `R19` and `R23`.
* **Richert's mechanism**, reproduced. Take a near-repdigit `x...xy` of length `4 <= n <= 20000`. Its permutation `x R_n + (y-x) 10^j` is divisible by a small prime `q <= 131` for some position `j`, because `10^j` runs through a full coset of `<10>` mod `q`. With Johnson 1977 (via Wikipedia, not read: no permutable prime uses three distinct digits, or two digits twice each), this is the no-non-repunit statement to 20000 digits. Richert 1951 (via Wikipedia and De Geest, not read) states it for `3 < n < 6*10^175`.

**3.4 Repunit primes, `n <= 1100` (FINITE-EXACT).**
* For composite `n`, `R_d | R_n` with `1 < R_d < R_n`.
* For the 184 primes `n <= 1100`, gmpy2's BPSW (`is_bpsw_prp`) and 25-round Miller–Rabin agree: `R_n` is prime exactly for **`n = 2, 19, 23, 317, 1031`**. This is OEIS A004023, whose first five terms are proved primes.
* Every prime factor `q < 10^5` of `R_p` with `5 <= p < 200` satisfies `q = 1 mod 2p` (`ord_q(10) = p`).

**3.5 Why "finitely many non-repunit circular primes" is conjectured.** A `k`-digit circular prime must avoid `0, 2, 4, 5, 6, 8`, so there are only `4^k` candidate strings, and **all `k` rotations must be prime**. Model the rotations as independent with the local corrections:
* each rotation is prime with probability `3.75/log N` (coprime to 30);
* at 3, all rotations share the digit sum;
* at every `q | 10^k - 1`, rotation multiplies by 10 mod `q`, so either all rotations or none are divisible by `q`.

The expected count then decays like `(c/k)^k` with `c = 4*3.75/log 10 = 6.5`.

| length `k` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 12 | 16 | 18 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| expected (numbers) | 12.1 | 11.1 | 9.1 | 3.8 | 8.7 | 0.56 | 0.43 | 0.06 | 0.034 | `8.5e-6` | `1.1e-5` | `4.5e-15` |
| actual | 8 | 12 | 8 | 10 | 12 | 0 | 0 | 0 | 0 | 0 | 0 (De Geest) | 0 (De Geest) |

* For `k = 2..6` the model predicts 44.8 numbers against 50 actual.
* For all `k >= 7` together it predicts **1.12 numbers (about 0.15 orbits)**, and none are known.
* Since the sum over `k` converges, finiteness is the natural conjecture. It is OEIS A293663's "Conjecture: The sequence is finite", and it is unproved.
* The spikes at `k = 12` and `k = 18` come from the many small primes dividing `10^k - 1`.
* **The repunits escape the decay** because a repunit is a single string (the fixed point), not `k` independent events. For prime `k`, its primality probability is of order `(log k)/k` (the Mersenne-type heuristic, whose factor comes from `q = 1 mod 2k`). The sum of these over prime `k` diverges, so infinitely many repunit primes are expected; this too is unproved.
* **The owner's dichotomy, finite orbit family versus infinite fixed-point family, is exactly this: `(c/k)^k` versus `(log k)/k`.**

## 4. Collatz's repunits: the dictionary

**4.1 The rising run is a repunit (PROVED; checked for every `n <= 2*10^5`).** For odd `x`, `T(x) + 1 = (3/2)(x + 1)`: **the odd branch is the dilation by `3/2` about the fixed point `-1`**. Hence, if `r = v_2(n+1)`:
* `T^j(n) + 1 = (3/2)^j (n+1)` for `j <= r`;
* `T^j(n)` is odd for `j < r`, and `T^r(n)` is even.

So **the initial run of odd steps of any `n` has length exactly `v_2(n+1)`**, the length of the all-ones tail of `n` in binary. Every integer carries, at its low end, a finite repunit that shadows the fixed point `-1 = ...1111_2`. The run lemma rests on 2-adic valuations, not on archimedean size, which is why it holds for every integer.

For `n = 2^k - 1` (checked `k <= 400` exactly, and by GMP iteration for `k <= 4000`):
* `T^j(2^k-1) = 3^j 2^(k-j) - 1` for `j <= k`;
* `T^k = 3^k - 1`;
* `T^(k+1) = (3^k-1)/2`, the base-3 repunit;
* `v_2(3^k - 1) = 1` for `k` odd and `v_2(k) + 2` for `k` even (lifting the exponent).

The longest rising run turns the `k`-th 2-adic approximant of `-1` into the `k`-th base-3 repunit.

**4.2 The families (PROVED unless marked; the ranges are checked exhaustively).**

| base `b` | numbers | forced prefix (exact) | stopping time `sigma` | total stopping time `tau` | height |
|---|---|---|---|---|---|
| `2` | `M_k = 2^k - 1` (Mersenne numbers) | `k` odd steps, then `v_2(k)+2` halvings (`k` even) or 1 (`k` odd) | `>= k+2`; least-squares fit `3.838k - 43` (the model predicts `3.819k`) to `k <= 12000` | `k + 1 + tau((3^k-1)/2)`; fit `8.630k + 57` (model `8.638k`) | `>= 3^k - 1`, with equality for 9551 of 11999 `k` (79.6%) |
| `3` | `(3^k-1)/2` | none: it is `T^(k+1)(M_k)` | **1 (`k` even), 2 (`k` odd)** | `tau(M_k) - k - 1` | random |
| `4` | `(4^k-1)/3` (trunk) | `1 0^(2k-1)`: `T = 2^(2k-1)` | **2** (`k >= 2`) | **`2k`** | **`2^(2k-1)`** |
| `-2`, odd `k` | `W_k = (2^k+1)/3` (Jacobsthal; prime: Wagstaff) | `T(W) = 2^(k-1)+1`, `T^(1+2i)(W) = 3^i 2^(k-1-2i)+1`, `T^k(W) = 3^((k-1)/2)+1` | **5** (`k >= 5`; 4 at `k = 3`) | random tail from `3^((k-1)/2)+1` | `3*2^(k-2) + 2`, except `k = 9` (FINITE-EXACT, odd `k <= 12000`) |
| `-2`, odd `k`, **3n-1 sheet** | same `W_k` | `T_-(W) = 2^(k-1)`: the minus trunk into 1 | 2 | `k` | `2^(k-1)` |
| `-2`, even `k = 2j` | `-(4^j-1)/3` | `T_-(x) = -2^(2j-1)`: the negative trunk of 3n-1 (the `x -> -x` image of the plus trunk) | | | |
| `10` | `R_k = (10^k-1)/9` | `T = (10^k+2)/6`, `T^(2+2i) = 3^i 5^k 2^(k-2-2i) + 1` (`2+2i <= k`): the word `1 1 (1 0)^*` | **8** (`k >= 8`) | random tail from `~ 8.66^k` | `3*5^k*2^(k-3) + 2`, except `k in {5,6,9,13,21,29,37}` (FINITE-EXACT, `k <= 3000`) |

*Proofs.*
* **Trunk.** `3R + 1 = 4^k`.
* **Base 3.**
  * `k` even: `(3^k-1)/2` is even, so one halving already goes below it.
  * `k` odd: `T((3^k-1)/2) = (3^(k+1)-1)/4`, which is even since `8 | 3^(k+1) - 1`, and `T^2 = (3^(k+1)-1)/8 < (3^k-1)/2`.
* **`W_k`.** A value `A 2^m + 1` (`A` odd, `m >= 2`) goes to `3A 2^(m-1) + 2` and then `3A 2^(m-2) + 1`. So `T^5 = 9*2^(k-5)+1 < W` while `T^1..T^4 >= W`.
* **`R_k`.** The same two-step rule applies with `A = 3^i 5^k`. Comparing `T^j` with `(10^k-1)/9`, the first value below it is `T^8 = (27/256)10^k + 1`.
* **Mersenne `sigma`.** The values up to index `k+1` are at least `(3^k-1)/2 > 2^k - 1`.
* **Pairing (PROVED).** For odd `x`, `T(x) = T(3x+1)`. With `x = (3^(2j-1)-1)/2` one has `3x + 1 = (3^(2j)-1)/2`, so **the orbits of `2^(2j-1)-1` and `2^(2j)-1` merge**, and `tau(2^(2j)-1) = tau(2^(2j-1)-1) + 1`. Checked for all 5999 odd `k < 12000`.
* The fitted slopes agree with the drift model: the tail must fall by `log(3^k/2)` at `log(2/sqrt3) = 0.1438` per step.

**4.3 The shadowing lemma (PROVED; checked for 18 even bases and `k <= 60`, 1080 cases).** Let `b` be even and `xi_b = 1/(1-b)`, the 2-adic limit of `R_k^(b) = (b^k-1)/(b-1)`. For `j <= k v_2(b)`, the first `j` parity bits of `R_k^(b)` are those of `xi_b`, and

`T^j(R_k^(b)) = T^j(xi_b) + 3^(a_j) b^k / ((b-1) 2^j)`,

where `a_j` is the number of odd steps of `xi_b` among its first `j` steps.

*Proof.* `T^j(x) = (3^a x + c)/2^j` with `c` depending only on the first `j` parity bits, and these agree for `x = y mod 2^j`. Here `R_k - xi_b = b^k/(b-1)` has valuation `k v_2(b)`. ∎

So each finite repunit is the 2-adic repunit (its "macrocosm") followed for exactly `k v_2(b)` steps, plus an explicit archimedean correction. The families of §4.2 are the cases `xi_2 = -1` (fixed), `xi_4 = -1/3`, `xi_(-2) = 1/3` and `xi_10 = -1/9`.

**4.4 The macrocosm: Collatz orbits of the 2-adic repunits (FINITE-EXACT for every even `|b| <= 64`, both sheets).**

| `b` | `xi_b` | `T_+`-orbit | parity word | `T_-`-orbit |
|---|---|---|---|---|
| 2 | `-1 = ...1111_2` | **fixed** | `1^inf` | `-1 -> -2 -> -1` |
| 4 | `-1/3 = ...1111_4 = E(0)` | `-1/3 -> 0` (fixed) | `1 0^inf` | `-> -1 -> -2` |
| `-2` | `1/3 = ...1111_(-2) = E_+(1) = E_-(0)` | `1/3 -> 1 -> 2 -> 1` | `1 (10)^inf` | `-> 0` |
| 10 | `-1/9 = ...1111_10` | `-1/9 -> 1/3 -> 1 -> 2` | `1 1 (10)^inf` | `-> -1 -> -2` after 3 steps |
| 6 | `-1/5` | `-> 1/5 -> 4/5 -> 2/5 -> 1/5` (the 3x+5 cycle `{1,4,2}`, scaled) | `1 (100)^inf` | periodic `{-1/5, -4/5, -2/5}` |
| 8 | `-1/7` | `-> 2/7 -> 1/7 -> 5/7 -> {5,11,20,10}/7` | `101 (1100)^inf` | `-> {-5/7, ...}` |
| `-4` | `1/5` | **periodic** `{1/5, 4/5, 2/5}` | `(100)^inf` | `-> {-1/5, ...}` |
| `4 - 2^p` | `1/(2^p - 3)` | **periodic**, word `1 0^(p-1)` (checked `p <= 11`) | | |

* **Periodic limits.** Among even `b in [-64, 64]`, `xi_b` is `T_+`-periodic exactly for `b = 2, -4, -10, -12, -16, -28, -40, -42, -54, -58, -60`.
* **Eventual periodicity.** Every `xi_b` computed is eventually periodic on both sheets, as Lagarias's Periodicity Conjecture predicts for all rationals with odd denominator (open in general).
* **`T` acts on the base by a Möbius map** (PROVED): `T(xi_b) = xi_(M(b))` with `M(b) = -(b+2)/(b-4)`. This gives `10 -> -2 -> 0` (then the cycle `{1,2}`, since `xi_0 = 1`), `6 -> -4` (periodic), `2 -> 2`, and `4 -> infinity` (`xi = 0`).
* **The coordinator's statement, corrected.** Of the four limits, only `-1` is a fixed point. `-1/3` is preperiodic to the fixed point 0; its fixed-point role is for the ladder map `p -> 4p+1`, the base-4 repunit generator. `1/3` and `-1/9` fall into the trivial cycle. The fixed-point reading becomes exact one level up, in the conjugacy of §4.5.

**4.5 The Bernstein–Lagarias conjugacy: its known fixed points and 2-cycles are repunit limits (PROVED values; FINITE-EXACT search).**
* **The setting.** `Phi` conjugates the 2-adic shift `S` to `T`: `Phi S Phi^(-1) = T` (Bernstein–Lagarias 1996, read). Its inverse is the parity-vector map `Q(x) = sum_i x_i 2^i`, where `x_i` is the parity of `T^i(x)`. So `Q(x) = x` means that **the binary digits of `x` are its own parity vector**.
* **What Bernstein–Lagarias found.** Their Fixed Point Conjecture says `Phi` has exactly two odd fixed points; they "immediately found two: `x = -1` and `x = 1/3`", and they knew one odd 2-cycle, `{1, -1/3}`.
* **Repunit readings.**
  * `-1 = ...1111` in base 2, and its parity vector is `1^inf`.
  * `1/3 = ...0101011_2 = ...1111` in base `-2`, and its parity vector `1,1,0,1,0,1,...` spells `1/3` again.
  * The 2-cycle joins `1 = R_1^(4)` and `-1/3 = R_inf^(4)`: the two ends of the base-4 trunk.
* **Our search.** Over odd `x = p/q` with `|p| <= 201` and odd `q < 100`, `Q` has exactly the fixed points `{-1, 1/3}` and exactly two 2-cycles: `{1, -1/3}` and **`{-1/5, 5/7}`**. A scratch run outside the deliverable output (`|p| <= 151`, periods up to 6, image denominators capped at `10^7`) found no cycles of period 3 to 6.
* **Checking the new cycle.**
  * `-1/5 -> 1/5 -> 4/5 -> 2/5 -> 1/5` has word `1 (100)^inf`, so `Q(-1/5) = 1 + 2/(1-8) = 5/7`.
  * `5/7 -> 11/7 -> 20/7 -> 10/7 -> 5/7` has word `(1100)^inf`, so `Q(5/7) = 3/(1-16) = -1/5`.
  * An independent computation of the first 80 parity bits from residues mod `2^160` confirms both, and the four known values.
* **Standing of the new cycle.**
  * It gives `F_1 >= 2` in Bernstein–Lagarias's notation, where they wrote `F_1 >= 1`.
  * Lagarias's annotated bibliography I (read) lists the same two fixed points and the single 2-cycle.
  * Hotzel's 2003 thesis (Chapter 7: periodic points of the conjugacy map; not read) may contain it, so **novelty is UNVERIFIED**.
  * Its point `-1/5 = ...1111_6` is again a repunit limit, of base 6.
* **Among all repunit limits** `xi_b` (`b` even, `|b| <= 64`, plus `xi_0 = 1`), exactly `b = 2, -2` are `Q`-fixed and `b = 0, 4, 6` lie on `Q`-2-cycles.

This is the exact content of the owner's "Brouwer" instinct on the Collatz side. `T`'s own fixed points, `0` and `-1`, are the images of the binary repdigit blocks `0^k` and `1^k` under the conjugacy (§4.6). The fixed points of the conjugacy itself are the base `+-2` repunits. The Fixed Point Conjecture (only these two) is **OPEN**.

**4.6 Rotation versus Collatz cycles.**
* **Rotation.** The shift `S` (drop the last binary digit) acts on the period-`k` points `-B/(2^k-1)`, where `B` is the `k`-bit block, by **rotating the block**: this is base-2 digit rotation (§3.1).
* **Cycles.** Since `T = Phi S Phi^(-1)`, each rational `T`-cycle is the image of the rotation orbit of its parity word `w`. Its points are `x_w = c_w/(2^p - 3^a)`, and `x_(w^m) = x_w`: a repeated word is the same point, just as `R_k/(10^k-1) = 1/9`.
* **Fixed points.** Rotation fixes the blocks `0^k` and `1^k`, and these correspond to `T`'s fixed points `0` and `-1`. For prime `p` there are `(2^p - 2)/p` cycles of exact period `p`.
* **Circular primes versus integral cycles (corrected 2026-09-25).** Circular primality, meaning every rotation is prime, **is rotation-invariant** by definition. Ordinary primality of one member is not: 19 is prime while its rotation91=7*13 is composite. On a finite Collatz cycle, one integral vertex forces every vertex integral because the integer map preserves integrality. This compares predicate propagation; it does not identify integral cycles with circular-prime orbits. See the [digit audit](ternary_digits_20260925.md) and MISTAKES for the correction lineage.
* **The integral cycles found.** The integral periodic points with primitive words of length `<= 18` are exactly the five known cycles:
  * `0` (word `0`);
  * `-1` (word `1`);
  * `{1,2}` (word `10`);
  * `{-5,-7,-10}` (word `110`);
  * `{-17,...}` (word `11110111000`).
* **Expected counts.** The expected number of integral `p`-cycles, `sum_a C(p,a)/(p |2^p - 3^a|)`, decays slowly: `0.54, 0.15, 0.060, 0.053, 0.028, 0.018` at `p = 5, 10, 20, 30, 40, 60`. It behaves like `2^(-0.05p)` with spikes at convergents of `log_2 3`. Both "finitely many" statements are therefore heuristically natural, but they decay at very different speeds: Collatz exponentially and slowly, circular primes like `(6.5/k)^k`.

**4.7 The 3-adic side.**
* **The base-3 repunits.** `(3^k-1)/2 = 0 -> 1 -> 4 -> 13 -> 40 -> 121 -> ...` under the unaccelerated odd step `x -> 3x+1` (the E-graph climb). They converge in `Z_3` to `-1/2`, the fixed point of that 3-adic contraction.
* **Base 2 to base 3.** The identity `T^(k+1)(2^k-1) = (3^k-1)/2` carries the `k`-th 2-adic approximant of `-1` (distance `2^k`) to the `k`-th 3-adic approximant of `-1/2` (distance `3^k/2`).
* **Every repunit family as an affine orbit.** Every base-`b` repunit sequence is the orbit of 0 under `x -> bx + 1`, whose fixed point is `-1/(b-1)`.
  * `b = 4` is the sibling ladder `S(p) = 4p+1` of the inverse tree, with fixed point `-1/3 = E(0)`.
  * `b = 3` is the E-graph climb.
  * `b = 2` is the all-ones binary tail, whose fixed point `-1` is `T`'s hostile fixed point.

## 5. "The length of all-one primes": what is real and what is numerology

**5.1 What is exactly true about lengths (PROVED).**
1. **The all-ones tail sets the rising run.** The length of the all-ones *tail* of any `n` is its rising-run length, `v_2(n+1)` (§4.1). The all-ones repunit `2^k - 1` is the smallest integer with a run of `k`, and the run density is `2^-k` (the Haar measure of the ball `-1 + 2^k Z_2`).
2. **Repunit primes have prime length.** A repunit prime has prime length in every base, because `R_d | R_k` (§3.1).
3. **Consecutive Mersenne lengths pair up.** The orbits of `2^(2j-1) - 1` and `2^(2j) - 1` merge after `2j+1` and `2j+2` steps (§4.2), so all Collatz statistics of `2^k - 1` come in pairs `(odd k, k+1)`.
4. **The exponent acts only 2-adically.** Beyond the forced prefix `1^k 0`, the first `m` parity bits of `2^k - 1` (the bits of `(3^k-1)/2`) depend only on `k mod 2^(m-1)`, since 3 has order `2^(m-1)` mod `2^(m+1)`; this is checked for `m <= 13`, `k < 3000`. So at every finite resolution the orbit data of `2^k - 1` are a 2-adically continuous function of the exponent `k`. The only exact `k`-effect is `v_2(3^k - 1) = v_2(k) + 2` halvings after the prefix for even `k` (mean 3 over even `k`).

**5.2 The test.** The families of lengths tested are:
* Mersenne exponents (A000043);
* Wagstaff exponents (A000978);
* base-3 repunit exponents (A028491);
* base-10 repunit exponents (A004023).

The reference is all odd primes `64 <= k <= 12000` (exponents below 64 sit in octaves with too few non-family primes, and there smallness dominates).

The statistics of `n = 2^k - 1` are:
* its stopping time;
* the total stopping time of its tail (the orbit of `(3^k-1)/2`);
* its height excess over `3^k - 1`;
* the odd fraction of its tail;
* the tail's parity-word discrepancy `max|a_j - j a/tau|/sqrt(tau)`;
* the tail's log-bridge supremum.

Linear trends in `k` are removed from the hitting times, and each statistic is standardized within its octave by the non-family primes. Each family is then compared with random primes of the same octaves (a stratified randomization test, 20000 draws).

A naive global comparison of whole-orbit discrepancies gives spurious "effects": `-0.77` sd for the Mersenne family in the output, kept for the record, and down to `-1.5` sd at `p < 0.001` in development runs with `k <= 2000`. These are **artefacts**: the forced prefix `1^k` makes the whole-orbit discrepancy grow like `sqrt(k)`, and the exponent families crowd the small `k`. The stratified test on the tail removes them.

| statistic of `2^k - 1` (mean z, p) | Mersenne (`n = 14`) | Wagstaff (`n = 16`) | base 3 (`n = 9`) | base 10 (`n = 2`) |
|---|---|---|---|---|
| stopping time `sigma` | +0.40, 0.085 | -0.20, 0.37 | -0.56, 0.085 | +0.17, 0.81 |
| `tau` of the tail | +0.12, 0.54 | -0.38, 0.12 | -0.16, 0.66 | +1.21, 0.11 |
| height excess over `3^k - 1` | -0.22, 0.52 | -0.01, 0.89 | -0.35, 0.41 | +0.78, 0.13 |
| odd fraction of the tail | +0.16, 0.47 | -0.34, 0.14 | -0.19, 0.56 | +1.29, 0.083 |
| word discrepancy of the tail | +0.11, 0.66 | +0.39, 0.15 | +0.34, 0.36 | +0.09, 0.92 |
| bridge supremum of the tail | +0.11, 0.66 | +0.39, 0.15 | +0.34, 0.35 | +0.10, 0.91 |

* **Result.** Of the 24 tests, none has `p < 0.05` (1.2 would be expected by chance). The smallest is 0.083, far from the Bonferroni level 0.0021.
* **Each family on the total stopping time of its own repunit:** Wagstaff exponents on `(2^k+1)/3`, `p = 0.52`; base-3 exponents on `(3^k-1)/2`, `p = 0.99`; base-10 exponents on `R_k`, `p = 0.93`.
* **Power.** With 14 Mersenne exponents, a mean shift of 0.75 sd would be seen at the 5% level with 80% power; the two base-10 exponents in range can only reveal shifts above 2 sd.

**5.3 A tested coincidence: the convergents of `log_2 3`.**
* Three of the five base-10 exponents `<= 1100` (`2, 19, 317`) are numerators of convergents or semiconvergents of `log_2 3` = `[1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, ...]`: `19/12` is a convergent and `317/200` a semiconvergent. The hypergeometric chance is `1.1*10^-3`, but the target set was chosen after seeing the data.
* The Mersenne and Wagstaff exponents `<= 1100` score *more* hits (`2, 3, 5, 19` and `3, 5, 11, 19`), because the numerator set `{1, 2, 3, 5, 8, 11, 19, 27, 46, 65, 84, 149, 233, 317, 401, 485, 569, 1054}` is dense among small numbers.
* Above 20, the only hit in all four families is `317`, and none of the six later base-10 exponents (49081, ..., 8177207) is a numerator.
* **NUMEROLOGY.**

**5.4 Verdict.**
* **REAL:** the length of the all-ones *tail* in binary is the Collatz rising run (§4.1). Prime length is forced for repunit primes by the divisibility lattice. The Collatz data of `2^k - 1` depend on `k` 2-adically and in pairs.
* **NUMEROLOGY:** whether the length is a Mersenne, Wagstaff, base-3 or base-10 repunit exponent. No statistic distinguishes those lengths from other primes of the same size; the power is enough to see a `0.75`-sd shift for the Mersenne family.
* **Mechanism.** The Collatz orbit of `2^k - 1` sees `k` through its 2-adic digits (and through `3^k`), while the primality of `2^k - 1` or `R_k` is decided by the multiplicative orders of 2 or 10 modulo primes `q = 1 mod 2k`. No known mechanism couples the two, and the data show no coupling.

## 6. Verdicts

| claim | verdict | why |
|---|---|---|
| primes `> 5` end in `1, 3, 7, 9`, and the transitions between them are not uniform | **REAL** (CITED; reproduced exactly) | Lemke Oliver–Soundararajan. The diagonal deficit is `-0.24` at `10^11` against `-0.20` predicted, and the ratio falls toward 1 |
| that non-uniformity "generates" the finite family `{337, 373, 733}` | **NOT SUPPORTED** | circular primes need all `k` rotations prime; rotations are not consecutive primes. Their scarcity is the `(c/k)^k` decay of `k` coincidences (§3.5). The consecutive-prime bias vanishes as `x` grows and enters nowhere |
| Brouwer: repunits are to rotation what fixed points are to a map | **REAL, PROVED** | Theorem R. Rotation is `x -> bx` on `Z/(b^k-1)` and on the circle; its fixed points are the repdigits, and the only multidigit prime ones are repunits. Every repunit is the single fixed point `1/(b-1)` |
| a finite family of orbits plus an infinite family of fixed points | **REAL as a heuristic dichotomy; both halves OPEN** | `k` coincidences decay like `(c/k)^k`; one coincidence per prime length (probability of order `(log k)/k`) diverges |
| "the length of all-one primes" is hidden Collatz structure | **NUMEROLOGY** | no statistic of `2^k - 1` separates the exponent families (§5.2). The convergent coincidence fails out of sample (§5.3) |
| micro/macrocosm | **REAL, two exact forms** | rotation: `R_k/(b^k-1) = 1/(b-1)` for all `k`. Collatz: the shadowing lemma (a finite repunit follows its 2-adic limit for exactly `k v_2(b)` steps), and the rising run of every `n` is a microcosm of `-1` |
| the Collatz analogue of the digit bias | **REAL, exact, first-order** | Theorem H: `3 -> 5` forced, `9 -> 9` with probability `8/15`, mod-3 classes i.i.d. `(1/3, 2/3)`, mod-9 stationary law `(8,16,11,4,2,22)/63`. It never decays, because the gap `v` has the same law at every scale |
| the Collatz fixed points are 2-adic repunits | **REAL** | `-1 = ...1111_2` is `T`-fixed. One level up, the Bernstein–Lagarias conjugacy's known odd fixed points are the base `+-2` repunits, its 2-cycles contain the base 0, 4 and 6 repunit limits, and the 2-cycle `{-1/5, 5/7}` is new relative to B–L 1996 |
| lead (c) as stated ("bases 2 and 4 give Collatz's 2-adic fixed points `-1` and `-1/3`") | **CORRECTED** | `T(-1/3) = 0`. `-1/3` is the fixed point of the ladder map and lies on a B–L 2-cycle |
| real orbits obey the Haar digit law | **CONFIRMED** (FINITE-EXACT) | deterministic orbits above `2^40` match it over `3*10^9` transitions; the apparent low-band deviations vanish on deduplication (§2.2) |

## 7. Reproduction

`python3 04-computation/experiments/procgen_repunit_20260924_run.py` writes [`procgen_repunit_20260924.out`](procgen_repunit_20260924.out). The whole run takes about 10 minutes on the mac-mini (105 + 87 + 281 + 119 s), runs one process at a time, and peaks below 700 MB (the largest allocation is the 63 MB bitset of Part 3).
* **Part 1:** `..._primes.py` and `.c`. The segmented sieve to `10^11` takes about 105 s CPU and under 2 MB.
* **Part 2:** `..._collatz.py` and `.c`. This covers `10^7` random starts, full orbits of `n <= 10^7`, and a `v`-histogram over `10^8` starts.
* **Part 3:** `..._repunits.py` and `..._circular.c`. The sieve to `10^9` is followed by the necklace search to length 16 (about 4 minutes).
* **Part 4:** `..._dictionary.py` and `..._collatz_repunits.c` (GMP, from `/opt/homebrew`).

Environment variables (`PR_*`) shrink every range for quick checks, and every check raises on failure. The random samples use `splitmix64` with fixed seeds, so every number in the output is deterministic.

## 8. Sources

**Read this session:**
* R. J. Lemke Oliver, K. Soundararajan, "Unexpected biases in the distribution of consecutive primes": arXiv:1603.03720v4 full text. Journal data (PNAS 113(31) (2016) E4446–E4454, PMID 27418603) via Europe PMC and Crossref (DOI 10.1073/pnas.1605366113).
* D. J. Bernstein, J. C. Lagarias, "The 3x + 1 conjugacy map", Canad. J. Math. 48 (1996) 1154–1169 (retypeset PDF at cr.yp.to): abstract, §2 Fixed Point Conjecture, §6.
* J. C. Lagarias, "The 3x+1 problem: an annotated bibliography" I (arXiv:math/0309224) and II (arXiv:math/0608208): entries on Bernstein–Lagarias, Akin, Hotzel, López–Stoll and Monks–Yazinski.
* OEIS (JSON): A016114, A068652, A293663, A003459, A004023, A000043, A000978, A028491, A001045.
* P. De Geest, "Circular Primes" (worldofnumbers.com): the table of lengths 1–27 and the near misses.
* Wikipedia: "Circular prime" and "Permutable prime" (wikitext).

**Not read (cited through the sources above):**
* H.-E. Richert (1951), Norsk Mat. Tidsskr. 33, 50–54.
* A. W. Johnson (1977), Math. Mag. 50, 100–103.
* W. Hotzel (2003), Hamburg dissertation.
* Montgomery–Soundararajan (2004), as used inside Lemke Oliver–Soundararajan.
* Lagarias (1985), Theorem L, via Bernstein–Lagarias.

**Repository parents:**
* the [inverse tree mod 192](collatz_procgen_20260924_inverse_tree_mod192.md) (sibling ladder `S(p) = 4p+1`, limit `-1/3`, trunk);
* the [pairings and transitions](procgen_brackets_20260924_pairings_transitions.md) note (`{2,3,11}`, Wagstaff and minus trunk);
* the [barrier atlas](collatz_procgen_20260922_barrier_atlas.md) (Bernstein–Lagarias as a reformulation);
* the synthesis, §2f.
