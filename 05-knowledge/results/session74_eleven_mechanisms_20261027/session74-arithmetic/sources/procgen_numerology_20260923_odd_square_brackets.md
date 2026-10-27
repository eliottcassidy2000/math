# The owner's odd-square prime brackets: {2,3,11} is a theorem, "every 5th" is a local-density effect, the rest is small-number coincidence

**Status.**
* **PROVED (elementary):** the primes `p` with `2p` in the same odd-square bracket as `p` are exactly `{2, 3, 11}`.
* **FINITE-EXACT** (brackets `n <= 30000`, i.e. odd squares up to `3.6*10^9`): the owner's distance sequences are reproduced exactly.
* **Explained (heuristic, standard):** the owner's "something happens every 5th occurrence" is a real `13σ` signal. It is the Hardy–Littlewood local-density (singular-series) effect of the odd square's residue mod 5. The same effect appears mod 3 and mod 7.
* **NUMEROLOGY:** the "n+3" prime counts, the "doubled 10 / tripled 20" and "fractal patterns into the macrocosm". No structure survives beyond the local densities.

Session `collatz-procgen-20260923`, orchestrator. Script: [`procgen_numerology_20260923_odd_square_brackets.py`](../../04-computation/experiments/procgen_numerology_20260923_odd_square_brackets.py). Output: [`procgen_numerology_20260923_odd_square_brackets.out`](procgen_numerology_20260923_odd_square_brackets.out). It is a vectorized segmented sieve and runs in 12 s.

## 1. Setting

The bracket is `B_n = ((2n-1)^2, (2n+1)^2]` for `n >= 1`. The owner's lists are the primes of `B_n`, with `2, 3, 11` set aside, followed by the square `(2n+1)^2`. Their "lengths" are (prime count) + 1.

## 2. The {2,3,11} theorem (PROVED)

`p in B_n` and `2p in B_n` require `(2n-1)^2 < p <= (2n+1)^2/2`. That needs `2(2n-1)^2 < (2n+1)^2`, i.e. `n(sqrt2 - 1)*2 < 1 + sqrt2`, i.e. `n <= 2`.
* In `B_1 = (1, 9]` the qualifying primes are `2` (with `4`) and `3` (with `6`).
* In `B_2 = (9, 25]` the only one is `11` (with `22`); `13` fails since `26 > 25`.

An exhaustive check over `p < 10^6` returns exactly `[2, 3, 11]`. So "11 escapes its primality within its bracket" is true. It is true because brackets grow more slowly than doubling from the third bracket on.

## 3. The sequences (FINITE-EXACT)

* **Prime counts, `n = 1..24`:** `4, 5, 6, 7, 8, 9, 9, 13, 11, 13, 14, 15, 15, 17, 16, 19, 19, 19, 21, 23, 20, 23, 23, 28`.
  * The owner's "5,6,7,8,9,10,10,14" is these counts plus 1.
  * The count equals `n+3` for `n <= 6` and then only sporadically (`n = 10, 11, 12, 14, 16, 20`). The true size is about `4n/log(2n+1)`, and the ratio has mean `1.0000` with no residue-class effect. The early `n+3` run is a small-number coincidence.
* **Largest prime to the next odd square, `n = 1..24`:** `2,2,2,2,8,2,2,6,2,2,6,6,2,2,8,2,2,2,10,12,2,8,2,2`. This matches the owner's list. Over `n <= 30000` the value `2` occurs `18.8%` of the time, because `(2n+1)^2 - 2` is often prime.
* **Previous odd square to the smallest prime:** `1,2,4,4,2,6,4,2,...`. With `2, 3, 11` removed it becomes `4,4,4,4,2,6,4,2,...`, which matches the owner.

## 4. "Every 5th occurrence" is real: the local density of the square's residue

Mean distance from `(2n-1)^2` to the next prime, by class, over `n <= 30000` (standard error about `0.22` per class):

| square `(2n-1)^2` mod 5 | classes `n mod 5` | mean |
|---|---|---|
| 0 (every 5th bracket) | 3 | **17.06** |
| 4 | 2, 4 | 19.12 |
| 1 | 0, 1 | 20.08 |

**Mechanism.** The candidates are `(2n-1)^2 + j` with `j = 2, 4, 6, 8, ...`.
* If `5` divides the square, then `j = 2, 4, 6, 8` are all prime to 5, so primes come sooner.
* If the square is `1 mod 5`, then `j = 4` is killed.
* If it is `4 mod 5`, then `j = 6` is killed.
* `j = 2` and `j = 8` are never killed, because `3` and `2` are not squares mod 5.

The ordering of the three means is exactly this ranking. The same effect appears mod 3 (mean 16.2 when `3` divides the square, against 20.5 otherwise) and mod 7 (17.97 against about 19.3). It is the Hardy–Littlewood singular series for primes of the form `x^2 + j`. The owner saw its mod-5 shadow.

For the distance *down* from `(2n+1)^2` the mod-5 means (17.66, 17.85, 17.54) are within noise, but the mod-3 and mod-7 effects are strong. The difference comes from which offsets are quadratic non-residues.

## 5. Relation to known problems and to this session

* Every bracket in range is nonempty; the minimum count for `n >= 2` is 5. Legendre's conjecture (a prime between consecutive squares) would force at least 2 primes per bracket, and Oppermann's refinement at least 4. Both are OPEN.
* No link to the Collatz/PC work was found. The owner's `{2, 3, 11}` recurs in the session only as numerology: `11/8 = 1 + 3/8` in THM-4467, and the level-11 eta product of `M_24`. We record the coincidence but claim no mechanism.
