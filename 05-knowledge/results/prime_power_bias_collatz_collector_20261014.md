# Prime powers, not primes: the square terms behind Chebyshev's race, Skewes' crossing and Pólya's conjecture (Mertens has none), computed to 10^10; the Collatz collector (odd steps land in 2 mod 3; the 2-adic ensemble visits 2 : 1; for px+1 the doubling coset, which is the squares at p = 7); the two sheets as one pair {3k+1, 3k+2}, the carry sign and the three cycles of 3x−1; and a three-way basin race for 3x−1 to 10^10

2026-10-14, session opus-2026-10-14-S32 (worktree `codex/session-prime-power-bias-collector-20261014`).

Owner's prompt (pasted, verbatim core): "merge in related to skewe's number, chebyshev bias, the prime race and mertens disproved conjecture all feature a large microcosmic region from 0 to a large number, and then a macrocosmic region up to infinity afterward. it has nothing to do with the primes themselves. It has to do with prime POWERS. Every odd prime, when squared, lands in residue 1 mod 4. [...] To keep the natural weighted count balanced, the actual primes ≡ 1 mod 4 have to come up short by ~√x / log x. [...] think 26861 and 1473. consider the structure of (3N+1)/2 accelerated collatz as opposed to its negative version (3N-1)/2 where the former maps odds 1,3,5,7,9 and 11 mod 12 all to values 2 mod 3, and the latter maps those to 1 mod 3, specifically one less than what the former maps to. you can think of each instance of 1,2 mod 3 as coming in a pair, and understand why the negative version produces 3 cycles while the positive only produces one, yet they both produce full integer coverage with just the 4 cycles. in collatz mod 3 we have 1 and 0 pointing to 2 mod 3, which cycles to itself, and in negative collatz we have 0 and 2 pointing to 1, which cycles to itself. note that primes squared mod 3 have 0 and 2 pointing to 1, which cycles to itself, and mod 4, 3 points to 1 which cycles to itself, and 2 points to 0, which cycles to itself."

**Status.**

- **PROVED (elementary):**
  - the prime-power identity of section 1, which gives the exact square term in every prime race and in `π − li`; and `λ = μ * 1_squares`, which separates Pólya from Mertens;
  - Proposition 2, the collector law: the exact law of `T_s^k(n) mod p` over the 2-adic ensemble, for every map `n/2, (pn + s)/2`; its limit on the doubling coset `s<2>`, with zero mass elsewhere; the Syracuse corollary;
  - Proposition 3: the two sheets share one pair `{3k+1, 3k+2}`; the carry sign; the carry budget.
- **FINITE-EXACT:**
  - `π(26861; 4, 1) = 1473 > 1472 = π(26861; 4, 3)`, with ties at `5, 17, 41, 461, 26833, 26849`;
  - all runs and excursions of the mod-4 race to `10^10`; no lead for `1 mod 3` to `10^10`;
  - Pólya's first excursion `906150257 ... 906488079`, maximum `829` at `906316571`;
  - `L(10^9)` and `M(10^9)`, `max |M|/√x` (my float32-assisted sieve, confirmed by the audit's exact integer sieve);
  - the collector law against brute force for every `k <= 12`, 8 primes, both signs;
  - every `n <= 10^10` reaches one of the three cycles of `3x − 1`, with the basin counts and every change in their order;
  - the integer-cycle census for all words of length `<= 24`.
- **VERIFIED (floating point with a large margin):**
  - `π(p) < li(p)` at every prime `p <= 10^10`. `li` is accumulated by Simpson's rule and re-anchored by a series each segment; below `10^5` it is exact (mpmath). The audit re-verified it with its own series.
- **NUMERICAL:**
  - log-grid averages and sign-change counts of the normalised races to `10^10`, and of `L(x)/√x` and `M(x)/√x` to `10^9`;
  - the 3-adic stationary law mod `3^m`, `m <= 9`;
  - the backward-tree balance to depth 50;
  - visit shares along actual trajectories to `10^7`;
  - basin densities.
- **CITED:**
  - Chebyshev (1853), Littlewood (1914), Skewes (1933, 1955), Leech (1957), Haselgrove (1958), Lehman (1960);
  - Bays–Hudson (1978, 1979, 2000), Tanaka (1980), Odlyzko–te Riele (1985), Pintz (1987), Rubinstein–Sarnak (1994);
  - Ng (2004), Kotnik–te Riele (2006), Granville–Martin (2006), Humphries (2013), Saouter–Trudgian–Demichel (2015), Platt–Trudgian (2016), Büthe (2018), Hurst (2018);
  - Rozmarynowycz–Kim (2023), Kim–Nguyen (2024), Bařina (2025);
  - Terras (1976), Tao (2022);
  - from this repository: THM-4631, THM-4641 row 2b, HYP-9258 and THM-4629.
- **CONDITIONAL (cited):** the logarithmic-mean statements require GRH (or RH) and further conditions, stated where used.
- **ANALOGY:** the microcosm/macrocosm reading (section 5) and the "Collatz form" sentences.
- **NUMEROLOGY:** section 6, items 2 and 7.
- **AUDITED.** An independent adversarial audit (2026-10-14) recomputed everything with its own code. It found eleven errors and a list of typing and citation issues, all corrected in place (section 8; MISTAKE-607).

## 0. Answers in brief

1. **"It has nothing to do with the primes themselves. It has to do with prime POWERS."**
   - **This is the classical explanation of the bias (Rubinstein–Sarnak).** It is a statement about logarithmic means, conditional on GRH and `L(½, χ) ≠ 0`. The latter holds for the characters here: `L(½, χ_(−4)) = 0.6677`, `L(½, χ_(−3)) = 0.4809`.
   - **The identity.** For a real character `χ` modulo `q`, the prime count `π(x; χ) = Σ_(p<=x) χ(p)` equals the prime-power count `Π(x; χ) = Σ_(p^k<=x) χ(p^k)/k`, minus `½ #{p <= √x : p ∤ q}` (the squares, all quadratic residues), minus smaller terms. `Π` has log-mean 0 under GRH, so `π` inherits, *on average*, a deficit of about `√x / log x` in the residue class.
   - **Measured on `[10^3, 10^10]`** (log-averages of the races normalised by `log x / √x`):

     | race | prime count only | prime powers, weight `1/k` |
     |---|---|---|
     | mod 4 | `+1.175`; share `> 0`: `0.9986`; 15 sign changes | `+0.020`; share `0.53`; 1070 sign changes |
     | mod 3 | `+1.186`; never negative | `+0.040`; share `0.54`; 1106 sign changes |
     | `π − li` | `−1.415`; never positive | `−0.023`; share `0.46`; 1671 sign changes |

     The difference between the two columns is the prime-power term by construction. The finding is that the corrected column is near 0, with about a thousand sign changes.
   - **The crossings are another matter.** They come from the zeros of `L(s, χ)` and `ζ(s)`, i.e. from the primes themselves, and the zeros decide where the microcosm ends.
2. **"Think 26861 and 1473."**
   - At `x = 26861`, `π(x; 4, 1) = 1473 > 1472 = π(x; 4, 3)`: the first lead of the residue class (Leech 1957).
   - Before it, class 1 never leads and only ties, at `p = 5, 17, 41, 461, 26833, 26849`.
   - By `x = 26861` the class has received the squares of the 37 odd primes up to 163 (`163^2 = 26569 < 26861 < 167^2`).
   - At the crossing, in normalised units, the race is `−0.06` and the prime-power term is `+1.19` (`+1.15` from the squares alone). So the zeros supplied about `−1.25` (`−1.21` measured against the squares alone).
   - Up to `10^10` the race has five excursions, starting at 26861, 616841, 12306137, 951784481 and 6309280709 (strict leads; maximum lead 1374).
3. **The four "microcosms" share one mechanism, except Mertens.**

   | phenomenon | first failure | where the bias comes from |
   |---|---|---|
   | Chebyshev mod 4 | 26861 | prime squares |
   | Chebyshev mod 3 | 608981813029 (Bays–Hudson) | prime squares |
   | Skewes | none known; the first known region of positivity is near `1.398·10^316`; none below `10^19` (Büthe) | prime squares |
   | Pólya | 906150257 (Tanaka; recomputed here) | squares through `Σ λ(n) n^(−s) = ζ(2s)/ζ(s)` |
   | **Mertens** | none known (none below `10^16`; one below `e^(1.96·10^19)`) | **none**: `Σ μ(n) n^(−s) = 1/ζ(s)`, and `μ` kills every square |

   - **Log-means.** Under RH and further conditions (Humphries 2013; Ng 2004):
     - `L(x)/√x` has log-mean `1/ζ(1/2) = −0.6848`; measured `−0.6813`;
     - `M(x)/√x` has log-mean 0; measured `−0.0021`.
   - **The difference between the two.** The identity `λ = μ * 1_squares` is the whole difference.
     - Pólya's long microcosm is a square bias.
     - Mertens' microcosm is the thin tail of an unbiased law: `max |M(x)|/√x = 0.464` on `[10^4, 10^9]`.
4. **Collatz mod 3, as stated, is exact.**
   - **The odd step.** `(3n+1)/2 ≡ 2` and `(3n−1)/2 ≡ 1 (mod 3)` for every odd `n`, and the two images are consecutive. For `n = 2k + 1` they are `3k + 2` and `3k + 1`, the "pair": the two non-multiples of 3 in the block `[3k, 3k+2]`. The positive sheet takes the upper member and the negative sheet the lower.
   - **The even step** multiplies the residue by `2^(−1) = 2`, i.e. swaps 1 and 2 and fixes 0.
   - **The collector theorem (Proposition 2) is an ensemble statement.** For `n` uniform modulo `3 · 2^k`, `T_+^k(n)` lies in class 2 with probability tending to `2/3`, and in class 1 with probability tending to `1/3`. For `3x − 1` it is the reverse. Class 0 is never visited after the first odd step.
   - **It is not a time average along one orbit.** For trajectories known to reach `{1, 2}`, that cycle visits the two classes `1 : 1`; universal entry remains open (MISTAKE-607's scope distinction).
     - **Measured.** Pooling all `T_+` trajectories up to `10^7` gives `0.6725 : 0.3275`.
     - **The excess over 2/3 is a length bias** (audit measurement). The pool is dominated by long trajectories. In their middle part the odd-step frequency is about `0.57`, and the class-2 share is `0.698`, close to the fair-coin chain's `1/(2 − q) = 0.700`. The final descent contributes a deficit.
   - **The general form.** For every `px ± 1` map the odd step is the constant `±2^(−1) (mod p)`, so the ensemble's residues collect on the coset `±<2>`, with weights `2^(d−m)/(2^d − 1)`.
5. **Your squaring comparison, typed.**
   - Mod 4 it is correct: `3 → 1 → 1`, `2 → 0 → 0`.
   - Mod 3, squaring sends `1, 2 → 1` but fixes 0 (`3^2 = 9 ≡ 0`). For primes `p ≠ 3`, `p^2 ≡ 1 (mod 3)` always.
   - The structural parallel that survives at every modulus is that squaring collects into the squares and the Collatz odd step collects into the doubling subgroup `<2>`. They coincide exactly when `ord_p(2) = (p − 1)/2`.
   - For the `7x+1` map the ensemble's residues mod 7 end up on the quadratic residues `{1, 2, 4}`, with weights `1/7, 2/7, 4/7`.
   - **Where those weights come from.** Their denominator is 7 because `7 = 2^3 − 1` is a Mersenne prime: the weights equal `r/p` exactly when `p = 2^d − 1`. Their support is `{1, 2, 4}` because `<2> = QR_7`. Matching them to the trivial cycle's code `{1, 2, 4}/7` is NUMEROLOGY; that code belongs to the unaccelerated cycle `1 → 4 → 2`.
   - **Two differences from the prime races:**
     - the collector bias is first order (odd steps are about half of all steps), while the prime-square bias is second order;
     - the collector class is *over*-visited, while the prime-square class is *under*-supplied with primes.
6. **Three cycles versus one.**
   - **The carry sign.**
     - A positive cycle of `3x + 1` with `l` odd steps among `L` steps needs `2^L > 3^l`, i.e. net contraction to pay the `+½` carries.
     - A positive cycle of `3x − 1` needs `3^l > 2^L`.
   - **The carry budget** bounds how small the least odd element must be: `n_min <= 1/(3(e^(|δ|/l) − 1))` above the line and `1/(3(1 − e^(−|δ|/l)))` below it, with `δ = L log 2 − l log 3`.
   - **The census.** The census of all words of length `<= 24` finds exactly the four known nonzero cycles:
     - `(L, l) = (2, 1)` above the line: `{1, 2}`;
     - `(1, 1)`, `(3, 2)`, `(11, 7)` below the line: `{1}`, `{5, 7, 10}` and `{17, …, 136}` for `3x − 1`.
   - **Where they come from.**
     - Three are Gersonides' near-equalities `3 − 2`, `4 − 3` and `9 − 8` (`|2^L − 3^l| = 1`; THM-4641 row 2b).
     - The fourth has `3^7 − 2^11 = 139` dividing its numerator.
   - **The asymmetry.** Among the known cycles, "three against one" is decided by the small exponents, by which side of `log_2 3` the small coincidences of `2^L` and `3^l` fall on. Whether other cycles exist, of either sign, is OPEN. It is not a mod-3 effect.
7. **"Full integer coverage with just the 4 cycles."**
   - For `3x + 1` on the positives this is the Collatz conjecture, verified to `2^71` (Bařina 2025).
   - For `3x − 1`, every `n <= 10^10` reaches one of its three cycles (computed here).
   - Both sheets have the same drift, `½ log(3/4)` per step.
   - The three basins of `3x − 1` give a genuine **Collatz race** (section 4):
     - densities at `10^10`: `0.32701`, `0.32498`, `0.34801`;
     - the order of the counts changes 204 times;
     - the lead changes 25 times, the last at `N = 1389`;
     - second place is last tied at `N = 30613`, and `{1}` leads `{5, 7, 10}` strictly from `30614` to `10^10`;
     - the gap between those two basins has settled near `0.0020` between `10^9` and `10^10`.
8. **Verdict.**
   - **Prime races.** The owner's mechanism is the classical one for Chebyshev, Skewes and Pólya (Rubinstein–Sarnak; Granville–Martin), stated for logarithmic means under GRH. Mertens does not fit it.
   - **Collatz.** The analogue is an exact ensemble law (the collector) and a numerically balanced backward tree. In its Collatz form, "prime powers, not primes" says that the known cycles are decided by the powers `2^L` and `3^l`. That is an ANALOGY: in the races "prime powers" means `p^k` for every prime `p`. Nothing here changes the difficulty of the Collatz conjecture.

## 1. Prime powers carry the bias (PROVED identity, CITED theory, FINITE-EXACT and NUMERICAL data)

**The identity.** Let `χ` be a real non-principal character modulo `q`, and `Π(x; χ) = Σ_(n<=x) Λ(n) χ(n)/log n = Σ_(k>=1) (1/k) Σ_(p^k<=x) χ(p)^k`. Since `χ(p)^2 = 1` for `p ∤ q`,

`π(x; χ) = Π(x; χ) − ½ #{p <= √x : p ∤ q} − Σ_(k>=3) (1/k) Σ_(p^k<=x) χ(p)^k`,

and likewise `π(x) = Π(x) − ½ π(√x) − ...` for Riemann's count `Π(x) = Σ_(p^k<=x) 1/k`. The audit checked the identity by enumerating prime powers directly.

**What it gives** (Rubinstein and Sarnak 1994):
- Under GRH, with `L(½, χ) ≠ 0`, the explicit formula gives `Π(x; χ)` and `Π(x) − li(x)` logarithmic mean 0.
- So `(log x/√x)(π(x;4,3) − π(x;4,1))` and `(log x/√x)(li(x) − π(x))` have logarithmic mean `1`.
- Under GRH and the linear independence of the zeros, the logarithmic densities of the leads are `0.9959` (mod 4), `0.9990` (mod 3) and `1 − 2.6×10^(−7)` (`li > π`).
- **The general rule:** for the race `π(x; q, a) − π(x; q, b)`, the normalised log-mean is `(c(q, b) − c(q, a))/φ(q)`, where `c(q, a)` counts the square roots of `a` modulo `q`. For `q = 3, 4` this is `(2 − 0)/2 = 1`. Non-residue classes beat residue classes because the squares fall on residues: the owner's "the square of any number is always a quadratic residue".

**Computed.** Part A sweeps every prime up to `10^10` in order; Part B sieves `λ` and `μ` up to `10^9`.

* **The mod-4 race.**
  * `π(10^10) = 455052511`, `π(10^10; 4, 1) = 227523275`, `π(10^10; 4, 3) = 227529235`.
  * Excursions where `π(x; 4, 1) > π(x; 4, 3)` strictly (runs merged when closer than 1 % of `x`; the gaps between excursions are factors of at least 6.6):

    | first prime of a strict lead | first prime where the strict lead is lost | runs | max lead |
    |---|---|---|---|
    | 26861 | 26863 | 1 | 1 |
    | 616841 | 633799 | 52 | 8 |
    | 12306137 | 12382327 | 75 | 24 |
    | 951784481 | 952223491 | 182 | 48 |
    | 6309280709 | 6403150199 | 3117 | 1374 (at 6345026777) |

  * 3427 runs in all.
  * Counting ties as Bays and Hudson do (as quoted by Granville–Martin), the regions read `... 952223506` and `6309280697 ... 6403150362`. The race is tied at `6309280697` and `6403150357`.
* **The mod-3 race.** `π(x; 3, 2) > π(x; 3, 1)` strictly for all `2 <= x <= 10^10`. This is consistent with Bays and Hudson's first crossing at `608981813029`.
* **Skewes (VERIFIED).** `π(p) < li(p)` at every prime `p <= 10^10`. The largest value of `π(p) − li(p)` is `−0.0452`, at `p = 2`; at `p = 3` it is `−0.1636` (exact `li` below `10^5`).
* **Pólya.**
  * The first `n >= 2` with `L(n) > 0` is `906150257`.
  * All 136 runs with `L > 0` lie in `[906150257, 906488079]`, with maximum `L = 829` at `906316571`.
  * No other `n <= 10^9` has `L(n) > 0`. Also `L(906180359) = +1` (Lehman's example).
* **Mertens.** `max |M(x)|/√x = 0.464162` on `[10^4, 10^9]`, at `x = 179919749`; `L(10^9) = −25216`, `M(10^9) = −222`.
* **Log-averages on `[10^4, 10^9]`:** `L(x)/√x = −0.6813` and `M(x)/√x = −0.0021`.

**Mertens is different.**
- `Σ λ(n) n^(−s) = ζ(2s)/ζ(s)`, i.e. `λ(n) = Σ_(d^2 | n) μ(n/d^2)`. The squares are inside `λ`, and the pole of `ζ(2s)` at `½` is Pólya's bias.
- `Σ μ(n) n^(−s) = 1/ζ(s)` has no such pole.
- Under RH and the Gonek–Hejhal condition, `M(x)/√x` has a limiting logarithmic distribution of mean 0 (Ng 2004). With LI its tails are very thin, and `|M(x)| < √x` fails on a set of positive but tiny logarithmic density.
- No counterexample is known. There is none below `10^16` (Hurst 2018), and one exists below `e^(1.96·10^19)` (Kim–Nguyen 2024; earlier `e^(1.017·10^29)`, Rozmarynowycz–Kim 2023).
- **Typing.** The owner's mechanism explains Chebyshev, Skewes and Pólya; it does not apply to Mertens.

## 2. The Collatz collector (PROVED + FINITE-EXACT + NUMERICAL)

Write `T_s(n) = n/2` for even `n` and `(pn + s)/2` for odd `n` (`p` an odd prime, `s = ±1`). Collatz is `p = 3, s = +1` (`T_+`); the negative version is `s = −1` (`T_−`, equivalently `3x+1` on the negative integers, since `T_−(n) = −T_+(−n)`).

**Proposition 2 (the collector law; an ensemble statement).** Let `d = ord_p(2)`.

1. For odd `n`, `T_s(n) ≡ s · 2^(−1) (mod p)`, a constant. For even `n`, `T_s(n) ≡ 2^(−1) n`.
2. **Exact finite law.** For `n` uniform modulo `p 2^k`,

   `P(T_s^k(n) ≡ r) = Σ_(j=0)^(k−1) 2^(−(j+1)) [r ≡ s 2^(−(j+1))] + 2^(−k)/p`.

3. **Limit.** As `k → ∞` the law converges to mass `2^(d−m)/(2^d − 1)` on `s 2^(−m)` (`m = 1, …, d`). It gives zero mass to `0` and to every unit outside the coset `s<2>`, which a given orbit never visits after its first odd step.

*Proof.*
- Item 1 holds because `pn ≡ 0`.
- **The parity vector.** The parity vector of the first `k` steps depends only on `n mod 2^k` and runs through `{0,1}^k` bijectively (Terras's argument, valid for every odd multiplier). By CRT it is independent of `n mod p`.
- **The residue after `k` steps.** If the last odd step is followed by `j` even steps, item 1 gives `T^k(n) ≡ s 2^(−(j+1))`; this happens with probability `2^(−(j+1))`. With no odd step at all (probability `2^(−k)`), `T^k(n) ≡ 2^(−k) n` is uniform.
- **The limit.** Summing the geometric weights over `j ≡ m − 1 (mod d)` gives the limit. ∎

*Checked:* brute force equals the formula for every `k <= 12`, for `p = 3, 5, 7, 11, 13, 17, 23, 31` and both signs (`.out` §1). The audit repeated this for 14 primes.

**Readings.**

* **`p = 3` (the owner's case).**
  * In the ensemble, `T_+` puts mass `(1/3, 2/3)` on classes `(1, 2)` and `T_−` puts `(2/3, 1/3)`.
  * Along trajectories known to reach the listed cycles, the statement is different: `{1, 2}` for `T_+` visits `1 : 1`; the listed `T_−` cycles visit class 1 with frequencies `1`, `2/3`, `8/11`. This does not assert that every orbit reaches those cycles.
  * Pooling all trajectories up to `10^7`, values after the first odd step:

    | map | visits | class 0 | class 1 | class 2 |
    |---|---|---|---|---|
    | `T_+` | `1.017·10^9` | `0` | `0.3275` | `0.6725` |
    | `T_−` | `0.809·10^9` | `0` | `0.6669` | `0.3331` |

  * **The excess over `2/3` is a length bias** (audit measurement).
    * Long trajectories dominate the pool.
    * In their middle the odd-step frequency is `q ≈ 0.572`, and the class-2 share `0.698` matches the fair-coin chain's `1/(2 − q) = 0.700`.
    * The final descent contributes a deficit, not the excess (the first draft said the opposite).
* **Syracuse.** For the odd-to-odd map `n ↦ (3n+1)/2^v`, the residue is `2^(−v) ≡ 2^v (mod 3)`, so in the ensemble `P(class 2) = P(v odd) = 2/3` exactly (`.out` §2: `4095 : 8190`, etc.).
* **Squares versus the doubling subgroup.**
  * The support `s<2>` equals the squares exactly when `ord_p(2) = (p − 1)/2` and `s` is a square: `p = 7, 17, 23, 41, 47, …` for `s = +1`, and `17, 41, 97, …` for `s = −1`.
  * When 2 is a primitive root (`p = 3, 5, 11, 13`), the support is every unit.
  * This is the stationary side of THM-4631, whose fast modes are the doubling cycles on the same units.
* **The 3-adic visit law.**
  * The stationary law of the fair-coin chain mod `3^m` is exact for `m = 1, 2`. Mod 9 it is: `1: 8/63`, `2: 16/63`, `4: 11/63`, `5: 4/63`, `7: 2/63`, `8: 22/63`.
  * For larger `m` it is computed by iteration. The conditional entropy of each new digit rises to `0.974 log 3` at `m = 9`.
  * The mass-weighted deviation of the top digit from uniform falls from `0.67` to `0.19`, while the worst residue stays near `0.41`.
  * So the law is biased at the coarse scale and approaches uniformity at finer scales only on average. Tao's fine-scale mixing (2022) is a statement of the same averaged kind, about the Syracuse random variable at level `n`, a different object.
* **The backward tree is balanced to the level of fluctuations (NUMERICAL).**
  * The backward tree of `T_+` from 1 (preimages `2n`, and `(2n−1)/3` when `n ≡ 2 mod 3`) has layer shares `0.3332, 0.3333, 0.3334` at depth 50 (`1185305` nodes). Its layer growth ratio is `4/3` to about 3 digits.
  * **The model.** A mean-field transfer model, which assumes the odd preimage of a class-2 node is uniform mod 3, has matrix `[[1, 0, 1/3], [0, 0, 4/3], [0, 1, 1/3]]` (rows: to classes 0, 1, 2; columns: from) and eigenvalues `4/3, 1, −1`. Its `4/3`-eigenvector `(1, 1, 1)` is the balance.
  * **The data.** The observed class-1/class-2 imbalance is fluctuation noise. It decays by a fitted factor `0.845` per layer relative to the layer size, `|c_1 − c_2|/√size` stays between `0.02` and `0.54` over depths 20–50, and its sign flips between consecutive layers in 19 of 30 cases. The first draft's "`(−3/4)^k` decay, as predicted" was not supported.
  * The tree contains every integer only if the Collatz conjecture holds.

**The owner's squaring remark, typed.**
- Squaring mod 4: `0, 2 → 0` and `1, 3 → 1`. Correct.
- Squaring mod 3: `1, 2 → 1` and `0 → 0` (`3^2 = 9`). "0 points to 1" is not right for squaring, though it is right for `T_−`'s odd step.
- On the units mod 3, `T_−`'s odd step coincides with squaring (both constant 1) and `T_+`'s with minus squaring.
- What carries over to every modulus is Proposition 2's coset `s<2>`, not the squares.

## 3. Two sheets, one pair, and the carry sign (PROVED + FINITE-EXACT + CITED)

**Proposition 3.**

1. **The pair.** For odd `n = 2k + 1`, `T_+(n) = 3k + 2` and `T_−(n) = 3k + 1`. Together the odd steps of the two sheets hit every non-multiple of 3 exactly once: the positive sheet the class 2, the negative sheet the class 1.
2. **Mirror.** `T_−(n) = −T_+(−n)`.
3. **The carry sign.** A rational cycle of `T_+` with `l` odd and `L` total steps (shortcut form) has `x_0 (2^L − 3^l) = c_w > 0`, with `c_w = Σ_k w_k 3^(#ones after k) 2^k`. So its points are positive iff `2^L > 3^l`. For `T_−` the sign is reversed.
4. **The carry budget.** Along a cycle, `Π (1 + s/(3 n_i)) = 2^L/3^l` over the odd elements `n_i`. With `δ = L log 2 − l log 3`, the least odd element satisfies

   `n_min <= 1/(3(e^(δ/l) − 1))` for `s = +1`, and `n_min <= 1/(3(1 − e^(−|δ|/l)))` for `s = −1`.

*Proof.*
1. `(3(2k+1) ± 1)/2 = 3k + 2` or `3k + 1`.
2. Direct.
3. Unroll `x_(k+1) = (3^(w_k) x_k + w_k)/2`.
4. Take logarithms of the per-step ratios `(3/2)(1 + s/(3n))` and `1/2`, and bound each term by its value at `n_min`. ∎

The audit verified the product identity exactly on all four cycles.

**The census.**
- All `2^25` words of length `L <= 24`: the integer cycles of `T_+` on Z are exactly `{0}`, `{1, 2}`, `{−1}`, `{−5, −7, −10}` and `{−17, −25, −37, −55, −82, −41, −61, −91, −136, −68, −34}`.
- They sit at `(L, l) = (2, 1)` above the line, and `(1, 1)`, `(3, 2)`, `(11, 7)` below.
- The budget table (`.out` §8):
  - at the four cycles the budget is `1, 1, 5.83, 35.7`, against actual least odd elements `1, 1, 5, 17`;
  - the large budgets that stayed empty include `(L, l) = (8, 5)` and `(16, 10)` (both `31.8`, above) and `(19, 12)` (`295`, below).

**Why three against one, among the known cycles.**
- The three unit-denominator cycles are Gersonides' near-equalities `3 − 2 = 1`, `4 − 3 = 1`, `9 − 8 = 1` (THM-4641 row 2b: `{2, 4, 8}` index them). Two of them, `(L, l) = (1, 1)` and `(3, 2)`, lie below the line.
- The fourth uses the semiconvergent `11/7 < log_2 3` with `3^7 − 2^11 = 139` and `c_w = 139 · 17`.
- Above the line the comparable chance `8/5` (`2^8 − 3^5 = 13`) gave no integer point.
- The negative sheet's budget is larger by exactly `1/3` for the same `|δ|/l`, since `1/(1 − e^(−t)) = 1 + 1/(e^t − 1)`. This is negligible next to which side the coincidences fall on.
- Whether further cycles exist, of either sign, is OPEN. Baker- and Rhin-type bounds, together with verification to `2^71`, force any further positive cycle to be very long.
- The owner's mod-3 pairing explains the two collector classes, not the cycle count.

## 4. The three basins of 3x − 1: a Collatz race (FINITE-EXACT; densities NUMERICAL)

**The computation.**
- Every `n <= 10^10` reaches `{1}`, `{5, 7, 10}` or the 11-cycle `{17, …, 136}` under `T_−` (`.out` and `collatz_collector_basins_1e10_20261014.out`).
- In the `10^10` run, six starting values had orbits passing `3·10^18` (largest value `1.85·10^19`, beyond 64-bit range). They were finished with exact Python integers. The first `10^10` attempt without this guard overflowed.

**The basin shares:**

| N | `{1}` | `{5,7,10}` | `{17,…}` | gap `{1} − {5,7,10}` |
|---|---|---|---|---|
| `10^3` | 0.349 | 0.306 | 0.345 | 0.043 |
| `10^4` | 0.3244 | 0.3213 | 0.3543 | 0.0031 |
| `10^5` | 0.3303 | 0.3210 | 0.3487 | 0.0093 |
| `10^6` | 0.32768 | 0.32335 | 0.34897 | 0.0043 |
| `10^7` | 0.32738 | 0.32450 | 0.34812 | 0.0029 |
| `10^8` | 0.32697 | 0.32471 | 0.34832 | 0.0023 |
| `10^9` | 0.326869 | 0.324952 | 0.348179 | 0.0019 |
| `2·10^9` | 0.326822 | 0.324982 | 0.348197 | 0.0018 |
| `5·10^9` | 0.327033 | 0.324993 | 0.347974 | 0.0020 |
| `10^10` | 0.327009 | 0.324985 | 0.348007 | 0.0020 |

**The race.**
- The order of the three counts changes 204 times.
- The leader changes 25 times, the last at `N = 1389`. After that the basin of the 11-cycle leads strictly, to `10^10`.
- Second place is last tied at `N = 30613` (`{5, 7, 10}` last strictly ahead at `30612`). `{1}` is strictly ahead for every `N` in `[30614, 10^10]`.
- The gap fell from `0.0043` (`10^6`) to `0.0019` (`10^9`) and then held near `0.0020` up to `10^10`.

**Status.**
- **Whether the densities exist is OPEN.** We know of no conditional result that would give them; the `3x − 1` conjecture alone does not.
- **Whether `{1}` and `{5, 7, 10}` have different densities is OPEN.** The data between `10^9` and `10^10` suggest a gap near `0.002`.
- If `d({1}) > d({5,7,10})`, `{1}` eventually leads for good. If the order is reversed, a further crossing must occur. Even with `d({1})` larger, a crossing above `10^10` is not excluded.
- We found no published basin densities for `3x − 1`.

### 4b. The repository's own Collatz prime race: the parity of the stopping time (CITED from the repository; ANALOGY)

The closest Collatz analogue of a prime race was already in the repository: the race between even and odd total stopping times `τ(n)` (HYP-9258, S25; THM-4629, S26).

- **Equal first-order densities.** The ROOT-supported parity function has ordinary mean zero (HYP-9258, PROVED from THM-4581). This is the analogue of Dirichlet's equidistribution.
- **An "explicit formula" with lines.** By THM-4629, every character `exp(2πi(ατ + βω)/m)` resonates at frequencies set by the near-returns `2^P ≈ 3^Q`. The Collatz analogue of the zeros is therefore indexed by coincidences of `2^P` and `3^Q`: the owner's slogan in its Collatz form, at the level of oscillations.
- **No square-type term.** The measured finite-`N` mean of `(−1)^τ` (`+0.00145` at `2^34`) is the Cesàro shadow of the lines (model decay `N^(−0.0018)`), not a fixed bias.
- **Classification.** The stopping-time parity race is of Mertens type (unbiased, oscillating), not Chebyshev type (a fixed square term). Section 4's basin race is a first-order race whose densities appear to differ.

## 5. Microcosm and macrocosm (ANALOGY, typed)

**In the prime races.**
- **The macrocosm law is a constant set by prime powers**, under GRH: the bias `+1` from squares, or `0` for Mertens.
- **It is scale-invariant in `log x`.** Bias and oscillation have the same size at every scale: `√x / log x` for the prime counts, `√x` for `L` and `M`.
- **The microcosm ends at a rare alignment of zeros.** Its length grows with bias over spread:
  - log densities `0.0041` (mod 4), `0.0010` (mod 3), `2.6·10^(−7)` (Skewes, conditional);
  - first failures `2.7·10^4` and `6.1·10^11`, and the first known positivity region of `π − li` near `10^316`.

**In Collatz.**
- **The macrocosm laws are first order:** contraction `½ log(3/4)` per step and the ensemble collector `1/3 : 2/3`.
- **The microcosm is where the carry `±1/(3n)` is not negligible**, i.e. where cycles can close. Proposition 3's budget gives `n_min <~ l/(3|δ|)`.
- **It is tiny:** every element of a known cycle has absolute value at most 136.
- **The known coincidences that fill it are between `2^L` and `3^l`:** Gersonides and 139; by Catalan–Mihăilescu, 8 and 9 are the last consecutive perfect powers.
- **Whether the microcosm holds all the cycles is OPEN.** Baker/Rhin-type bounds and verification to `2^71` only force any other positive cycle to be long.

**The shared shape.** A long-run law decided by powers, and a microcosm decided by small coincidences: alignments of zeros for the primes, near-equalities of powers of 2 and 3 for Collatz. This is an analogy of shape, not a transfer of mechanism.

## 6. Hostile controls and numerology

1. **(Correction.)** "Nothing to do with the primes themselves": true of the mean bias, false of the crossings, which the zeros produce.
2. **(NUMEROLOGY.)** `26861` just past `163^2 = 26569`, the square of Ramanujan's Heegner number. The first lead is a zero alignment; the squares only set the bias.
3. **(Count.)** `1473` is `π(26861; 4, 1)`, with no further structure.
4. **(Correction.)** Mod-3 squaring: "0 points to 1" is a slip (`0 → 0`); see section 2.
5. **(Typing.)** The bias directions differ. Prime squares leave the residue class *short of primes*; the Collatz collector makes its class *over-represented* in the ensemble. The analogue of "primes" is the ensemble of iterates, and of "prime powers" the tree.
6. **(Typing.)** The collector bias is first order; prime races have equal first-order densities (Dirichlet). No Collatz phenomenon here needs zeros. The closest genuine races are the stopping-time parity race (4b) and the basin race (4).
7. **(NUMEROLOGY.)** The `7x+1` weights `1/7, 2/7, 4/7` versus the S17 code `{1, 2, 4}/7`.
   - The denominator is the Mersenne coincidence `7 = 2^3 − 1`: the weights equal `r/p` iff `p = 2^d − 1`, checked at `p = 3, 7, 31`.
   - The support is the Paley coincidence `<2> = QR_7`. At the Paley primes 23 and 47 the support is QR, but the denominators are `2^11 − 1` and `2^23 − 1`.
   - The S17 code belongs to the unaccelerated cycle `1 → 4 → 2`. Under the shortcut map the trivial cycle `{1, 2}` has code `{1, 2}/3`, which matches the `p = 3` weights `1/3, 2/3`.

## 7. Open threads

1. **The `3x − 1` basin densities (FINITE, extendable; OPEN as theorems).**
   - Push to `10^11` with a two-pass method.
   - Model the densities as "harmonic measure from infinity" of the microcosm.
2. **A Collatz race with a provable square-type term (OPEN, speculative).** Find a Collatz counting function whose first-order densities are equal by symmetry, so that a collector term decides the race at second order.
3. **The 3-adic visit law (OPEN).** Does the conditional entropy of the new digit tend to 1, and at what rate? Compare with the repository's Syracuse-law Fourier profile (S20).

## 8. Audit

An independent adversarial audit (2026-10-14, blind subagent, its own code: parallel sieve with exact log-averages, an exact integer `λ`/`μ` sieve, a CRT-state brute force for Proposition 2, an independent census and basin run with explicit tie handling) reproduced every number. It found eleven errors, all corrected above (MISTAKE-607):

1. The Skewes maximum was misstated; it is `−0.0452` at `p = 2`.
2. The collector law was stated as a time average along trajectories; it is an ensemble law.
3. The excess over `2/3` was attributed to the final descent; it is a length bias, and the final descent gives a deficit.
4. The tree's `(−3/4)^k` claim was unsupported; the imbalance is fluctuation noise.
5. "Exactly balanced tree" overstated a numerical observation.
6. The order-change count was 205; it is 204.
7. The race logic at 30613 was wrong: 30613 is a tie, and a further crossing is not excluded.
8. The bias formula's `/2` should be `/φ(q)`.
9. The 7 in the `7x+1` weights is the Mersenne coincidence, not the Paley one.
10. "The macrocosm's absence of cycles" overstated what is known; existence of further cycles is OPEN.
11. Minor: the `(l, L)` ordering, the `−1.21` figure, and the counted set `{p ∤ q}`.

**Typing changes:**
- conditional statements are marked;
- puns are retyped as ANALOGY;
- float-based checks are typed VERIFIED;
- the brute-force range now covers every `k <= 12`;
- the transfer model is stated with its matrix and assumption.

**Citation venues corrected:** Hurst, *Math. Comp.* 87; Ng, *Proc. London Math. Soc.* 89; Humphries, *J. Number Theory* 133. Büthe's `10^19` is in his abstract.

**What survives:** every computation, Propositions 2 and 3, the prime-power identity, and the classification of Pólya versus Mertens.

## 9. Reproduction and references

**Reproduction.**
- `python3 04-computation/experiments/prime_power_bias_races_20261014.py` (numba, a few minutes): races to `10^10`, `λ` and `μ` to `10^9`.
- `python3 04-computation/experiments/collatz_collector_basins_20261014.py [N]` (`N = 10^9` default: about 30 s and 1 GB; `N = 10^10`: about 3 minutes and 10 GB).
- Outputs `.out` and `.json` beside this note.

**References.**
- **Chebyshev's bias and prime races:**
  - P. L. Chebyshev (1853);
  - J. Leech (1957);
  - C. Bays and R. H. Hudson, *Math. Comp.* (1978) and their tables of the mod-4 and mod-8 regions (1979);
  - M. Rubinstein and P. Sarnak, *Chebyshev's bias*, *Experiment. Math.* 3 (1994);
  - A. Granville and G. Martin, *Prime number races*, *Amer. Math. Monthly* 113 (2006), arXiv:math/0408319.
- **Skewes:**
  - J. E. Littlewood (1914);
  - S. Skewes (1933, 1955);
  - C. Bays and R. H. Hudson, *Math. Comp.* 69 (2000);
  - Y. Saouter, T. Trudgian and P. Demichel, *Math. Comp.* 84 (2015) (the exponent `727.951336105` from secondary sources);
  - D. Platt and T. Trudgian (2016), `θ(x) < x` below `1.39·10^17`;
  - J. Büthe, *Math. Comp.* 87 (2018), `li(x) > π(x)` for `2 <= x <= 10^19`.
- **Pólya:**
  - C. B. Haselgrove, *Mathematika* 5 (1958);
  - R. S. Lehman (1960), 906180359;
  - M. Tanaka, *Tokyo J. Math.* 3 (1980), 906150257;
  - P. Humphries, *J. Number Theory* 133 (2013) 545–582.
- **Mertens:**
  - A. Odlyzko and H. te Riele, *J. reine angew. Math.* 357 (1985);
  - J. Pintz (1987);
  - T. Kotnik and H. te Riele (2006);
  - N. Ng, *Proc. London Math. Soc.* (3) 89 (2004) 361–389;
  - G. Hurst, *Math. Comp.* 87 (2018) 1013–1028;
  - J. Rozmarynowycz and S. Kim, arXiv:2305.00345 (2023);
  - S. Kim and P. Q. Nguyen, ANTS XVI, arXiv:2502.21021 (2024).
- **Collatz:**
  - R. Terras (1976);
  - D. Bařina, *J. Supercomputing* 81 (2025);
  - T. Tao, *Forum Math. Pi* (2022);
  - THM-4631 and THM-4641 ([graph_theory_as_common_carrier_20261013.md](graph_theory_as_common_carrier_20261013.md), [repair_algebra_minus_one_and_wilson_dictionary_20261014.md](repair_algebra_minus_one_and_wilson_dictionary_20261014.md));
  - HYP-9258 and THM-4629;
  - the Paley and S17 notes ([collatz_paley_bridge_20261001.md](collatz_paley_bridge_20261001.md), [ramanujan_heegner7_trivial_cycle_20261006.md](ramanujan_heegner7_trivial_cycle_20261006.md)).
