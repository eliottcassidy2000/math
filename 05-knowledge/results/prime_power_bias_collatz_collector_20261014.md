# Prime powers, not primes: the square terms behind Chebyshev's race, Skewes' crossing and Pólya's conjecture (Mertens has none), computed to 10^10; the Collatz collector (odd steps land in 2 mod 3, trajectories visit 2 : 1; for px+1 the doubling subgroup, which is the squares at p = 7); the two sheets as one pair {3k+1, 3k+2}, the carry sign and why 3x−1 has three cycles; and a three-way basin race for 3x−1

2026-10-14, session opus-2026-10-14-S32 (worktree `codex/session-prime-power-bias-collector-20261014`).

Owner's prompt (pasted, verbatim core): "merge in related to skewe's number, chebyshev bias, the prime race and mertens disproved conjecture all feature a large microcosmic region from 0 to a large number, and then a macrocosmic region up to infinity afterward. it has nothing to do with the primes themselves. It has to do with prime POWERS. Every odd prime, when squared, lands in residue 1 mod 4. [...] To keep the natural weighted count balanced, the actual primes ≡ 1 mod 4 have to come up short by ~√x / log x. [...] think 26861 and 1473. consider the structure of (3N+1)/2 accelerated collatz as opposed to its negative version (3N-1)/2 where the former maps odds 1,3,5,7,9 and 11 mod 12 all to values 2 mod 3, and the latter maps those to 1 mod 3, specifically one less than what the former maps to. you can think of each instance of 1,2 mod 3 as coming in a pair, and understand why the negative version produces 3 cycles while the positive only produces one, yet they both produce full integer coverage with just the 4 cycles. in collatz mod 3 we have 1 and 0 pointing to 2 mod 3, which cycles to itself, and in negative collatz we have 0 and 2 pointing to 1, which cycles to itself. note that primes squared mod 3 have 0 and 2 pointing to 1, which cycles to itself, and mod 4, 3 points to 1 which cycles to itself, and 2 points to 0, which cycles to itself."

**Status.**

- **PROVED (elementary):**
  - the prime-power identity of section 1 (the exact square term in every prime race and in `π − li`), and `λ = μ * 1_squares` for Pólya versus Mertens;
  - Proposition 2: the collector law. The exact law of `T_s^k(n) mod p` for every map `n/2, (pn + s)/2`, its limit on the doubling coset `s<2>`, zero mass elsewhere; the Syracuse corollary;
  - Proposition 3: the two sheets share one pair `{3k+1, 3k+2}`; the carry sign; the carry budget for a cycle's least odd element.
- **FINITE-EXACT** (scripts [prime_power_bias_races_20261014.py](../../04-computation/experiments/prime_power_bias_races_20261014.py) and [collatz_collector_basins_20261014.py](../../04-computation/experiments/collatz_collector_basins_20261014.py); outputs [.out](prime_power_bias_races_20261014.out) and [.out](collatz_collector_basins_20261014.out), both ending `ALL CHECKS PASSED`):
  - `π(26861; 4, 1) = 1473 > 1472 = π(26861; 4, 3)`;
  - all excursions of the mod-4 race to `10^10`; no lead for `1 mod 3` to `10^10`; `π(p) < li(p)` at every prime `p <= 10^10`;
  - Pólya's first excursion: `906150257 ... 906488079`, maximum `L = 829` at `906316571`;
  - the collector law against brute force for 8 primes and both signs;
  - every `n <= 10^9` reaches one of the three cycles of `3x − 1`, with the basin counts and every change in their order;
  - the integer-cycle census for all words of length `<= 24`.
- **NUMERICAL:**
  - log-grid averages and sign-change counts of the normalised races to `10^10`, and of `L(x)/√x`, `M(x)/√x` to `10^9`;
  - the 3-adic visit law mod `3^m`, `m <= 9`;
  - the backward-tree balance to depth 50.
- **CITED:**
  - Chebyshev (1853), Littlewood (1914), Skewes (1933, 1955), Leech (1957), Haselgrove (1958), Lehman (1960);
  - Bays–Hudson (1978, 2000), Tanaka (1980), Odlyzko–te Riele (1985), Rubinstein–Sarnak (1994), Ng (2004), Granville–Martin (2006), Humphries (2013);
  - Saouter–Trudgian–Demichel (2015), Büthe (2018), Hurst (2018), Kim–Nguyen (2024), Bařina (2025);
  - THM-4631 (doubling on units, the Paley primes) and THM-4641 row 2b (the Gersonides cycles) of this repository.
- **ANALOGY:** the microcosm/macrocosm reading of section 5.
- **NUMEROLOGY:** section 6.
- **Not yet audited** (section 8).

## 0. Answers in brief

1. **"It has nothing to do with the primes themselves. It has to do with prime POWERS."** For the *bias* this is exactly right and classical.
   - **The identity.** For a real character `χ` modulo `q` (here `χ_(−4)`, `χ_(−3)`), the prime count `π(x; χ) = Σ_(p<=x) χ(p)` is the prime-power count `Π(x; χ) = Σ_(p^k<=x) χ(p^k)/k` minus `½ #{odd p <= √x}` (the squares, all quadratic residues), minus smaller terms.
   - **What it implies.** Under GRH, `Π` has no bias, so `π` inherits a deficit of about `√x / log x` in the residue class (section 1).
   - **Measured on `[10^3, 10^10]` (log-averages, normalised by `log x / √x`):**

     | quantity | prime count only | prime powers, weight `1/k` |
     |---|---|---|
     | race mod 4 | `+1.175`, share `> 0`: `0.9986`, 15 sign changes | `+0.020`, share `0.53`, 1070 sign changes |
     | race mod 3 | `+1.186`, never negative | `+0.040`, share `0.54`, 1106 sign changes |
     | `π − li` | `−1.415`, never positive | `−0.023`, share `0.46`, 1671 sign changes |

     The prime-power term alone accounts for the whole difference (`+1.154`, `+1.145`, `+1.392`: it is `1 + 2/log x + ...` at finite `x`).
   - **What it does not mean.** The *crossings* do have to do with the primes themselves: they are made by the zeros of `L(s, χ)` and `ζ(s)`, which set where the microcosm ends (26861, 616841, ...).
2. **"Think 26861 and 1473."**
   - At `x = 26861`, `π(x; 4, 1) = 1473` and `π(x; 4, 3) = 1472`: the first lead of the residue class (Leech 1957).
   - Below it, class 1 never leads; it ties at `p = 5, 17, 41, 461, 26833, 26849`.
   - The class has by then received the squares of all 37 odd primes up to 163 (`163^2 = 26569 < 26861 < 167^2`).
   - At the crossing the square term is worth `+1.15` in normalised units, so the zeros had to supply about `−1.21`.
   - Up to `10^10` the race has exactly five excursions, starting at 26861, 616841, 12306137, 951784481 and 6309280709 (lead 1374 at most).
3. **The four "microcosms" share one mechanism, except Mertens.**

   | phenomenon | first failure | where the bias comes from |
   |---|---|---|
   | Chebyshev mod 4 | 26861 | prime squares |
   | Chebyshev mod 3 | 608981813029 (Bays–Hudson) | prime squares |
   | Skewes | below `e^727.95` (Saouter–Trudgian–Demichel), above `10^19` (Büthe) | prime squares |
   | Pólya | 906150257 (Tanaka; recomputed here) | squares through `ζ(2s)`: `Σ λ(n) n^(−s) = ζ(2s)/ζ(s)`; log-mean of `L(x)/√x` is `1/ζ(1/2) = −0.6848` (measured `−0.6813`) |
   | **Mertens** | none known (no counterexample below `10^16`; one below `e^(1.96·10^19)`) | **none**: `Σ μ(n) n^(−s) = 1/ζ(s)`, `μ` kills every square and `M(x)/√x` has log-mean 0 (measured `−0.002`) |

   The identity `λ = μ * 1_squares` is the whole difference between the two.
   - **Pólya.** Its long microcosm is a square bias: `L(x)` leads negative.
   - **Mertens.** Its microcosm is the thin tail of an *unbiased* law: `max |M(x)|/√x = 0.464` on `[10^4, 10^9]`.
4. **Collatz mod 3, as stated, is exact.**
   - **The odd steps.** `(3n+1)/2 ≡ 2` and `(3n−1)/2 ≡ 1 (mod 3)` for every odd `n`, and the two images are consecutive integers. For `n = 2k + 1` they are `3k + 2` and `3k + 1`, the "pair": the two non-multiples of 3 in the block `[3k, 3k+2]`. The positive sheet takes the upper member and the negative sheet the lower.
   - **The even steps** multiply the residue by `2^(−1) = 2`, i.e. swap 1 and 2 and fix 0.
   - **The collector theorem (Proposition 2).** Along trajectories class 2 is visited with limiting frequency exactly `2/3` and class 1 with `1/3` (the reverse for `3x − 1`). Class 0 is never visited after the first odd step.
     - Measured on all trajectories up to `10^7`: `0.6725 / 0.3275`.
     - The odd-only Syracuse map gives the same `2/3`, from the parity of the halving exponent.
   - **The general form:** for every `px ± 1` map the odd step is the constant `±2^(−1) (mod p)`, so the trajectory's residues collect into the coset `±<2>` of the doubling subgroup, with weights `2^(d−m)/(2^d − 1)`.
5. **Your squaring comparison, typed.**
   - Mod 4 it is correct: `3 → 1 → 1`, `2 → 0 → 0`.
   - Mod 3, squaring sends `1, 2 → 1` but fixes 0 (`3^2 = 9 ≡ 0`). For primes `p ≠ 3`, `p^2 ≡ 1 (mod 3)` always.
   - The structural parallel that survives at every modulus: squaring collects into the squares; the Collatz odd step collects into the doubling subgroup `<2>`.
   - These coincide exactly when `ord_p(2) = (p − 1)/2`. For the `7x+1` map the residues mod 7 end up on the quadratic residues `{1, 2, 4}`, with weights `1/7, 2/7, 4/7`: the same numbers as the trivial cycle's code in `ramanujan_heegner7_trivial_cycle_20261006.md`.
   - **Two differences.** The Collatz bias is first order (odd steps are half of all steps), while the prime-square bias is second order (`√x` of `x`). And it goes the other way: the collector class is *over*-visited, while the prime-square class is *under*-supplied with primes.
6. **Three cycles versus one.**
   - **The carry sign.** Each odd step adds a carry of `+1/2` (positive sheet) or `−1/2` (negative sheet), and a cycle has to pay it back.
     - A positive cycle with `l` odd steps and `L` total steps needs `2^L > 3^l` (net contraction).
     - A cycle of the negative sheet needs `3^l > 2^L`.
   - **The carry budget.** Exactly how small the least odd element must be: `n_min <= 1/(3(e^(|δ|/l) − 1))` above the line and `1/(3(1 − e^(−|δ|/l)))` below, with `δ = L log 2 − l log 3`.
   - **The four cycles are the four lattice points that paid** (the census of every word of length `<= 24` finds nothing else):
     - `(L, l) = (2, 1)` above, i.e. `{1, 2}`;
     - `(1, 1)`, `(3, 2)`, `(11, 7)` below, i.e. `{1}`, `{5, 7, 10}`, `{17, ..., 136}` for `3x − 1`.
   - **Three are the Gersonides near-equalities** `3 − 2`, `4 − 3`, `9 − 8` (`|2^L − 3^l| = 1`; THM-4641 row 2b); the fourth has `3^7 − 2^11 = 139` dividing its numerator.
   - **So "three versus one" is a microcosm statement about where the small coincidences between powers of 2 and 3 fall.** Two of the three Gersonides pairs, and the semiconvergent `11/7`, lie below the line. It is not a mod-3 effect.
7. **"Full integer coverage with just the 4 cycles."**
   - For `3x + 1` on the positives this is the Collatz conjecture, verified to `2^71` (Bařina 2025).
   - For `3x − 1`, every `n <= 10^9` reaches one of its three cycles (computed here).
   - Both sheets have the same drift, `½ log(3/4)` per step, so trajectories fall into the microcosm, where the cycles that exist there catch them.
   - The three basins of `3x − 1` give a genuine **Collatz race** (section 4):
     - densities at `10^9`: `0.32687`, `0.32495`, `0.34818`;
     - the lead changes 25 times, the last at `N = 1389`;
     - second place changes 204 times, the last at `N = 30613`;
     - the gap between the basins of `{1}` and `{5, 7, 10}` is still shrinking (`0.0043` at `10^6`, `0.0019` at `10^9`).
8. **Verdict.**
   - **Prime races.** The owner's mechanism is the classical one (Rubinstein–Sarnak; Granville–Martin's exposition) for Chebyshev, Skewes and Pólya. Mertens is the exception.
   - **Collatz.** The analogue is exact at the level of the collector: odd steps collect into one residue class as squaring collects into the squares, with an exact `2 : 1` visit law and an exactly balanced tree. The cycle count is decided by the powers `2^L`, `3^l`, not by the integers: "prime powers, not primes" in its Collatz form.
   - **No consequence for the Collatz conjecture**, whose difficulty is unchanged.

## 1. Prime powers carry the bias (PROVED identity, CITED theory, FINITE-EXACT and NUMERICAL data)

**The identity.** Let `χ` be a real non-principal character modulo `q`, `Π(x; χ) = Σ_(n<=x) Λ(n) χ(n)/log n = Σ_(k>=1) (1/k) Σ_(p^k <= x) χ(p)^k`. Since `χ(p)^2 = 1` for `p ∤ q`,

`π(x; χ) = Π(x; χ) − ½ #{p <= √x : p ∤ q} − Σ_(k>=3) (1/k) Σ_(p^k<=x) χ(p)^k`,

and in the same way `π(x) = Π(x) − ½ π(√x) − ...` for Riemann's count `Π(x) = Σ_(p^k<=x) 1/k`.

**What the identity gives.** It is the bookkeeping behind the theorem of Rubinstein and Sarnak (1994).
- The explicit formula expresses `Π(x; χ)` (and `Π(x) − li(x)`) through the zeros. Under GRH it has logarithmic mean 0, given `L(½, χ) ≠ 0`, which holds for `χ_(−4)`, `χ_(−3)`.
- So `(log x/√x)(π(x;4,3) − π(x;4,1))` and `(log x/√x)(li(x) − π(x))` have logarithmic mean `1`, the normalised `½ π(√x)`.
- Under GRH and the linear independence of the zeros, the logarithmic densities of the leads are `0.9959` (mod 4), `0.9990` (mod 3) and `1 − 2.6×10^(−7)` (`li > π`).
- The general rule: a non-residue class beats a residue class, because the squares fall on residues. The bias is `(c(q, b) − c(q, a))/2`, where `c(q, a)` is the number of square roots of `a`. This is the owner's "the square of any number is always a quadratic residue".

**Computed** (Part A sweeps every prime up to `10^10` in order; Part B sieves `λ` and `μ` up to `10^9`):

* **The mod-4 race.**
  * `π(10^10) = 455052511`, `π(10^10; 4, 1) = 227523275`, `π(10^10; 4, 3) = 227529235`.
  * The excursions where `π(x; 4, 1) > π(x; 4, 3)` (runs of leads merged when closer than 1 % of `x`):

    | starts at | ends at | runs | max lead |
    |---|---|---|---|
    | 26861 | 26863 | 1 | 1 |
    | 616841 | 633799 | 52 | 8 |
    | 12306137 | 12382327 | 75 | 24 |
    | 951784481 | 952223491 | 182 | 48 |
    | 6309280709 | 6403150199 | 3117 | 1374 |

  * In all, 3427 maximal runs below `10^10`.
* **The mod-3 race.** No `x <= 10^10` has `π(x;3,1) > π(x;3,2)`, consistent with Bays–Hudson's first crossing at `608981813029`.
* **Skewes.** `π(p) < li(p)` at every prime `p <= 10^10`. The largest value of `π − li` is `−0.165`, at `p = 3`; below `10^5` this was checked with exact `li`.
* **Pólya.**
  * The first `n >= 2` with `L(n) > 0` is `906150257` (Tanaka 1980).
  * The excursion lasts to `906488079` (136 runs), with maximum `L = 829` at `n = 906316571`.
  * No other `n <= 10^9` has `L(n) > 0`.
* **Mertens.** `max |M(x)|/√x = 0.4642` on `[10^4, 10^9]` (at `x = 179919749`).
* **Log-averages on `[10^4, 10^9]`.** `L(x)/√x`: `−0.6813` (prediction `1/ζ(½) = −0.6848`, Humphries 2013); `M(x)/√x`: `−0.0021` (prediction 0).
* **The log-grid table of section 0** (step `10^(−3)` in `log x`).

**Mertens is different.**
- `Σ λ(n) n^(−s) = ζ(2s)/ζ(s)`, i.e. `λ(n) = Σ_(d^2 | n) μ(n/d^2)`. The squares are inside `λ`, and the pole of `ζ(2s)` at `½` is Pólya's bias.
- `Σ μ(n) n^(−s) = 1/ζ(s)` has no such pole: `μ` vanishes on every multiple of a square.
- Under RH, LI and a growth condition, `M(x)/√x` has a limiting logarithmic distribution of mean 0 with very thin tails (Ng 2004). The Mertens inequality `|M(x)| < √x` fails on a set of positive but tiny logarithmic density.
- The conjecture was disproved without an example (Odlyzko–te Riele 1985). Today a counterexample is known to lie below `e^(1.96·10^19)` (Kim–Nguyen), and there is none below `10^16` (Hurst).
- **Typing.** The owner's "prime powers" mechanism explains Chebyshev, Skewes and Pólya, and explicitly does not apply to Mertens.

## 2. The Collatz collector (PROVED + FINITE-EXACT + NUMERICAL)

Write `T_s(n) = n/2` for even `n` and `(pn + s)/2` for odd `n` (`p` an odd prime, `s = ±1`). Collatz is `p = 3, s = +1` (the shortcut map `T_+`); the negative version is `s = −1` (`T_−`, equivalently `3x+1` on the negative integers, since `T_−(n) = −T_+(−n)`).

**Proposition 2 (the collector law).** Let `d = ord_p(2)`.

1. For odd `n`, `T_s(n) ≡ s · 2^(−1) (mod p)`, a constant. For even `n`, `T_s(n) ≡ 2^(−1) n`.
2. **Exact finite law.** For `n` uniform modulo `p 2^k`,

   `P(T_s^k(n) ≡ r) = Σ_(j=0)^(k−1) 2^(−(j+1)) [r ≡ s 2^(−(j+1))] + 2^(−k)/p`.

3. **Limit.** As `k → ∞` the law of `T_s^k(n) mod p` converges to the measure with mass `2^(d−m)/(2^d − 1)` on `s 2^(−m)`, `m = 1, …, d`. It gives zero mass to `0` and to every unit outside the coset `s<2>`, which after the first odd step is never visited at all.

*Proof.*
- Item 1 holds because `pn ≡ 0`.
- **The parity vector.** The parity vector of the first `k` steps depends only on `n mod 2^k` and runs through `{0,1}^k` bijectively (Terras's argument, valid for every odd multiplier). By CRT it is independent of `n mod p`.
- **The residue after `k` steps.** If the last odd step is followed by `j` even steps, item 1 gives `T^k(n) ≡ s 2^(−(j+1))`; this happens with probability `2^(−(j+1))`. With no odd step at all (probability `2^(−k)`), `T^k(n) ≡ 2^(−k) n` is uniform.
- **The limit.** Summing the geometric weights over `j ≡ m − 1 (mod d)` gives the limit. ∎

*Checked:*
- the formula equals brute-force counts for `p = 3, 5, 7, 11, 13, 17, 23, 31`, both signs, `k <= 12` (`.out` §1);
- the limits are printed there.

**Readings.**

* **`p = 3` (the owner's case).**
  * `T_+` visits classes `(1, 2)` with frequencies `(1/3, 2/3)` and `T_−` with `(2/3, 1/3)`. Class 0 is transient.
  * On actual trajectories up to `10^7`, after the first odd step:

    | map | visits | class 0 | class 1 | class 2 |
    |---|---|---|---|---|
    | `T_+` | `1.017·10^9` | `0` | `0.3275` | `0.6725` |
    | `T_−` | `0.809·10^9` | `0` | `0.6669` | `0.3331` |

  * The excess over `2/3` for `T_+` comes from the deterministic descent `… → 8 → 4 → 2 → 1` at the end.
* **Syracuse.** For the odd-to-odd map `n ↦ (3n+1)/2^v`, the residue is `2^(−v) ≡ 2^v (mod 3)`, so `P(class 2) = P(v odd) = 2/3` exactly (`.out` §2: `4095 : 8190`, etc.).
* **Squares versus the doubling subgroup.**
  * The support `s<2>` equals the squares exactly when `ord_p(2) = (p − 1)/2` and `s` is a square: `p = 7` (`s = +1`), `17` (both signs), `23` (`s = +1`), ...
  * For the `7x+1` map the limit law mod 7 is `4/7, 2/7, 1/7` on the residues `4, 2, 1`.
  * When 2 is a primitive root (`p = 3, 5, 11, 13`), the support is every unit, with geometric weights.
  * This is the stationary side of THM-4631, whose fast modes are the doubling cycles on the same units.
* **The 3-adic visit law.**
  * The stationary law of the fair-coin chain mod `3^m` is exact for `m = 1, 2`. Mod 9 it is: `1: 8/63`, `2: 16/63`, `4: 11/63`, `5: 4/63`, `7: 2/63`, `8: 22/63`.
  * Up to `m = 9` it is computed by iteration. The new digit's conditional entropy rises: `0.58`, `0.85`, `0.90`, `0.93`, `0.94`, `0.96`, `0.96`, `0.97`, `0.974` (`log 3 = 1`).
  * The mass-weighted deviation of the top digit from uniform falls from `0.67` to `0.19`, but the worst residue stays near `0.41`.
  * **The 3-adic law is biased at the coarse scale and approaches uniformity at finer scales only on average.** Tao's fine-scale mixing (2022) is a statement of the same averaged kind, about a different object: the Syracuse random variable at level `n`. The first draft of this section said "uniform at the fine scale"; the computation refuted that before the note was written.
* **The tree is balanced.**
  * The backward tree of `T_+` from 1 (preimages `2n`, and `(2n−1)/3` when `n ≡ 2 mod 3`) has layer shares `0.3332, 0.3333, 0.3334` for classes `0, 1, 2` at depth 50 (`1185305` nodes).
  * Its layer growth ratio is `4/3` (to 4 digits).
  * A two-class transfer model gives eigenvalues `4/3` and `−1`. So the class-1/class-2 imbalance of layer `k` decays like `(−3/4)^k` relative to the layer size; the data oscillate in sign and shrink, as predicted.
  * **The complete count (the tree, i.e. all integers) is balanced; the trajectory visits are biased.** This is the Collatz form of the race statement "prime powers balanced, primes biased".

**The owner's squaring remark, typed.**
- Squaring mod 4: `0, 2 → 0` and `1, 3 → 1`. Correct.
- Squaring mod 3: `1, 2 → 1` and `0 → 0` (`3^2 = 9`). "0 points to 1" is not right for squaring, though it is right for `T_−`'s odd step.
- On the units mod 3, `T_−`'s odd step coincides with squaring (both constant 1) and `T_+`'s with minus squaring.
- **What carries over** to every modulus is Proposition 2's coset `s<2>`, not the squares.

## 3. Two sheets, one pair, and the carry sign (PROVED + FINITE-EXACT + CITED)

**Proposition 3.**

1. **The pair.** For odd `n = 2k + 1`, `T_+(n) = 3k + 2` and `T_−(n) = 3k + 1`. The odd steps of the two sheets together hit every non-multiple of 3 exactly once: the positive sheet the class 2, the negative sheet the class 1.
2. **Mirror.** `T_−(n) = −T_+(−n)`. The cycles of `3x − 1` on the positives are the negatives of the negative cycles of `3x + 1`.
3. **The carry sign.** A rational cycle of `T_+` with `l` odd and `L` total steps (shortcut form) has `x_0 (2^L − 3^l) = c_w > 0`, with `c_w = Σ_k w_k 3^(#ones after k) 2^k`. So its points are positive iff `2^L > 3^l`, i.e. `L/l > log_2 3`. For `T_−` the sign is reversed: positive cycles need `3^l > 2^L`.
4. **The carry budget.** Along a cycle, `Π (1 + s/(3 n_i)) = 2^L/3^l` over its odd elements `n_i`. If `δ = L log 2 − l log 3`, the least odd element satisfies

   `n_min <= 1/(3(e^(δ/l) − 1))` for `s = +1`, and `n_min <= 1/(3(1 − e^(−|δ|/l)))` for `s = −1`.

*Proof.*
1. `(3(2k+1) ± 1)/2 = 3k + 2` or `3k + 1`.
2. Direct.
3. Unroll `x_(k+1) = (3^(w_k) x_k + w_k)/2`.
4. Take logarithms of the product of the per-step ratios `(3/2)(1 + s/(3n))` and `1/2`, and bound each term by its value at `n_min`. ∎

**The four cycles.**
- **Census** (all `2^25` words of length `L <= 24`): the integer cycles of `T_+` on Z are exactly `{0}`, `{1, 2}`, `{−1}`, `{−5, −7, −10}` and `{−17, −25, −37, −55, −82, −41, −61, −91, −136, −68, −34}`. The repo's earlier censuses (the eleventh S15 note to `A = 22`; THM-4641's `|n| <= 100001`) agree.
- **Where they sit:**
  - `(L, l) = (2, 1)`, above the line;
  - `(1, 1)`, `(3, 2)`, `(11, 7)`, below the line.
- **The budget table** (`.out` §8):
  - it is `1, 1, 5.83, 35.7` at those four points, against actual least odd elements `1, 1, 5, 17`;
  - the large budgets that stayed empty are `(5, 8)` (`31.8`, above) and `(12, 19)` (`295`, below).

**Why three against one.**
- The three unit-denominator cycles are Gersonides' near-equalities `3 − 2 = 1`, `4 − 3 = 1`, `9 − 8 = 1` (THM-4641 row 2b: `{2, 4, 8}` index them). Two of them, `(1, 1)` and `(3, 2)`, lie below the line.
- The fourth cycle uses the semiconvergent `11/7 < log_2 3` with `3^7 − 2^11 = 139` and `c_w = 139 · 17`.
- Above the line the only comparable chance, `8/5` with `2^8 − 3^5 = 13`, produced no integer point.
- So the asymmetry is decided entirely in the microcosm, exponents `<= 11`, by which side of `log_2 3` the lucky small coincidences of `2^L` and `3^l` fall on.
- The negative sheet's budget is larger by exactly `1/3` for the same `|δ|/l` (`1/(1 − e^(−t)) = 1 + 1/(e^t − 1)`): rounding down gives slightly more room. That is negligible next to the Diophantine luck.
- **The owner's mod-3 pairing explains the two collector classes, not the cycle count.**

## 4. The three basins of 3x − 1: a Collatz race (FINITE-EXACT; densities NUMERICAL)

Every `n <= 10^9` reaches `{1}`, `{5, 7, 10}` or the 11-cycle `{17, …, 136}` under `T_−` (largest value met: `1.12·10^17`). The basin shares:

| N | `{1}` | `{5,7,10}` | `{17,…}` |
|---|---|---|---|
| `10^3` | 0.349 | 0.306 | 0.345 |
| `10^4` | 0.3244 | 0.3213 | 0.3543 |
| `10^5` | 0.3303 | 0.3210 | 0.3487 |
| `10^6` | 0.32768 | 0.32335 | 0.34897 |
| `10^7` | 0.32738 | 0.32450 | 0.34812 |
| `10^8` | 0.32697 | 0.32471 | 0.34832 |
| `10^9` | 0.326869 | 0.324952 | 0.348179 |

**The race.**
- **Leader.** The order of the three counts changes 205 times. The leader changes 25 times, the last at `N = 1389`, after which the basin of the 11-cycle leads to `10^9`.
- **Second place** changes 204 times, the last at `N = 30613`, after which `{1}` stays ahead of `{5, 7, 10}`.
- **The gap.** The gap between them shrinks slowly: `0.0043` (`10^6`), `0.0029` (`10^7`), `0.0023` (`10^8`), `0.0019` (`10^9`).

**Status.**
- **Whether the densities exist is OPEN.** Their existence would follow from a strong form of the conjecture for `3x − 1`.
- **Whether those of `{1}` and `{5, 7, 10}` differ is OPEN.**
- If they were equal, this would be a Chebyshev-type race, decided by lower-order terms, with crossings possibly recurring far out. If they differ, the race has ended and 30613 is its last crossing.
- We found no published basin densities for `3x − 1`.

## 5. Microcosm and macrocosm (ANALOGY, typed)

**In the prime races.**
- **The macrocosm law is an exact constant set by prime powers:** the bias `+1` from squares, or `0` for Mertens.
- **The macrocosm is scale-invariant in `log x`:** bias and oscillation are both of size `√x / log x` at every scale.
- **The microcosm ends at a rare alignment of zeros.** Its length is a large-deviation statement whose rarity grows with bias over spread:
  - log densities `0.0041` (mod 4), `0.0010` (mod 3), `2.6·10^(−7)` (Skewes);
  - first failures `2.7·10^4`, `6.1·10^11`, `~10^316`.

**In Collatz.**
- **The macrocosm law is also exact, but first order:** contraction `½ log(3/4)` per step, collector frequencies `1/3 : 2/3`, the balanced tree.
- **The microcosm is where the carry `±1/(3n)` is not negligible**, i.e. where cycles can close. Proposition 3's budget makes this precise: a cycle with `l` odd steps needs `n_min <~ l/(3|δ|)`.
- **It is tiny.** All cycle elements are at most 136.
- **The coincidences that fill it are between prime powers `2^L` and `3^l`.** Gersonides and 139; Catalan–Mihăilescu says `8, 9` is the last pair of consecutive powers.
- **The macrocosm's absence of cycles is the Diophantine statement** that `|2^L − 3^l|` cannot be small relative to `3^l/l^c` (Baker, Rhin), combined with verification to `2^71`.

**The shared slogan holds in both: the long-run law is decided by powers, the microcosm by small coincidences.**
- For the primes the coincidences are alignments of zeros.
- For Collatz they are near-equalities of powers of 2 and 3.
- This is an analogy of shape, not a transfer of mechanism.

## 6. Hostile controls and numerology

1. **"Nothing to do with the primes themselves."** True of the bias, false of the crossings, which the zeros produce. Typed above.
2. **`26861` and `163^2 = 26569`.** That the first lead falls just past the square of Ramanujan's Heegner number is NUMEROLOGY: the first lead is a zero alignment, and the square count only sets the bias.
3. **`1473`** is a count (`π(26861; 4, 1)`), with no further structure.
4. **Mod-3 squaring.** "0 points to 1" is a slip (`0 → 0`); corrected in section 2.
5. **Bias directions differ.** Prime squares make the residue class *short of primes*; the Collatz collector makes its class *over-visited*. The analogue of "primes" in Collatz is the trajectory visits, and of "prime powers" the whole tree. Both are typed, not merged.
6. **First versus second order.** The Collatz collector bias is a density-one effect; prime races have equal first-order densities (Dirichlet). There is no Collatz phenomenon here that needs GRH or zeros. The closest analogue of a genuine race is section 4's basin race, if two basin densities turn out equal.
7. **`7x+1` and the S17 code.** The weights `1/7, 2/7, 4/7` equal the trivial cycle's code `{1, 2, 4}/7`. The identity is exact: a geometric law wrapped mod `d` is the binary code of the word with one `1`. Its appearance at `p = 7` is the Paley coincidence `<2> = QR_7` (THM-4631).

## 7. Open threads

1. **The `3x − 1` basin densities (FINITE, extendable; OPEN as theorems).**
   - Push to `10^11` with a two-pass method to see whether the `{1}`/`{5, 7, 10}` gap keeps shrinking.
   - Model the densities as "harmonic measure from infinity" of the microcosm.
2. **A Collatz race with a provable square-type term (OPEN, speculative).** Look for a Collatz counting function whose first-order densities are equal by symmetry, so that a collector term decides the race at second order, as prime squares do.
3. **The 3-adic visit law (OPEN).** Its fine-scale structure: is the conditional entropy of the new digit tending to 1, and at what rate? Compare with the repo's Syracuse-law Fourier profile (S20).

## 8. Audit

Pending: an independent adversarial audit will be recorded here.

## 9. Reproduction and references

**Reproduction.**
- `python3 04-computation/experiments/prime_power_bias_races_20261014.py` (about 2 minutes, numba): races to `10^10`, `λ` and `μ` to `10^9`.
- `python3 04-computation/experiments/collatz_collector_basins_20261014.py` (about 30 seconds plus 1 GB of memory): the collector law, the tree, the basins to `10^9`, the census.
- Outputs `.out` and `.json` beside this note.

**References.**
- **Chebyshev's bias and prime races:**
  - P. L. Chebyshev, letter (1853);
  - J. Leech (1957), the first lead at 26861;
  - C. Bays and R. H. Hudson, *Math. Comp.* (1978), the mod-3 crossing 608981813029;
  - M. Rubinstein and P. Sarnak, *Chebyshev's bias*, *Experiment. Math.* 3 (1994);
  - A. Granville and G. Martin, *Prime number races*, *Amer. Math. Monthly* 113 (2006), arXiv:math/0408319.
- **Skewes:**
  - J. E. Littlewood (1914);
  - S. Skewes (1933, 1955);
  - C. Bays and R. H. Hudson, *Math. Comp.* 69 (2000), crossings near `1.39822·10^316`;
  - Y. Saouter, T. Trudgian and P. Demichel, *Math. Comp.* (2015), the first crossing below `e^727.951336105`;
  - J. Büthe (2018), `π(x) < li(x)` for `x <= 10^19` (from secondary sources);
  - D. Platt and T. Trudgian (2016), to `1.39·10^17`.
- **Pólya:**
  - C. B. Haselgrove, *Mathematika* 5 (1958);
  - R. S. Lehman (1960), the counterexample 906180359;
  - M. Tanaka, *Tokyo J. Math.* 3 (1980), 906150257;
  - P. Humphries, *Acta Arith.* (2013), the distribution of `L(x)/√x` (from memory, not re-fetched).
- **Mertens:**
  - A. Odlyzko and H. te Riele, *J. reine angew. Math.* 357 (1985);
  - J. Pintz (1987);
  - T. Kotnik and H. te Riele (2006);
  - N. Ng, *Math. Comp.* 73 (2004);
  - G. Hurst, *Experiment. Math.* (2018), `M(x)` to `10^16` and `limsup M(x)/√x > 1.826`, `liminf < −1.837`;
  - J. Rozmarynowycz and S. Kim (2023);
  - S. Kim and P. Q. Nguyen (2024), the bound `e^(1.96·10^19)`.
- **Collatz:**
  - R. Terras (1976), the parity-vector bijection;
  - D. Bařina, *J. Supercomputing* 81 (2025), verification to `2^71`;
  - T. Tao, *Forum Math. Pi* (2022), almost all orbits attain almost bounded values (fine-scale mixing of Syracuse variables);
  - THM-4631 and THM-4641 (this repository; [graph_theory_as_common_carrier_20261013.md](graph_theory_as_common_carrier_20261013.md), [repair_algebra_minus_one_and_wilson_dictionary_20261014.md](repair_algebra_minus_one_and_wilson_dictionary_20261014.md));
  - the Paley and S17 notes ([collatz_paley_bridge_20261001.md](collatz_paley_bridge_20261001.md), [ramanujan_heegner7_trivial_cycle_20261006.md](ramanujan_heegner7_trivial_cycle_20261006.md)).
