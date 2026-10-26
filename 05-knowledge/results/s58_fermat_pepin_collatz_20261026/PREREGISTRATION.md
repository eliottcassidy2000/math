# S58 pre-registration: primes on the Collatz Fermat and Mersenne lines

opus-2026-10-26-S58. These predictions were computed by
`04-computation/fermat_pepin_s58/collatz_line_primes_20261026.py predict` and committed **before** the
`test` mode, which runs the primality tests, was executed. The prediction output is
`collatz_line_primes_prediction.json`.

**Why this test.** The only quantitative argument that `{3, 5, 17, 257, 65537}` is the complete list of Fermat primes
is the random model. Boklan–Conway (arXiv:1605.01371) bound `P(F_n prime)` by `4/2^n`, which sums to `2^-30` over
`n >= 33`; the cruder model `2/ln F_n` gives `6.7e-10`. The test applies the same model to integers that Collatz
itself attaches to `F_n` and `M_p`. These are the exact Syracuse runs:
- the valuation-two Fermat line `x_j = 3^j 2^(m-2j) + 1` (THM-4793 (F));
- the valuation-one Mersenne line `y_k = 3^k 2^(p-k) - 1` (THM-4556).

**Model.** Every value is sieved exactly by all primes `q <= 2^16`, using the congruences
`(3/4)^j = -2^-m` and `(3/2)^k = 2^-p (mod q)`. A surviving value `x >= 2^32` is prime with probability
`prod_{q <= 2^16} (1 - 1/q)^-1 / ln x`. Values below `2^32` are classified exactly and are not part of the test.

| id | family | predicted count of probable primes among survivors `>= 2^32` | window (E ± 2.5 sd) |
|---|---|---|---|
| L-P1 | Fermat-type lines of `2^m + 1`, `m = 8..1200` | `E = 2707.11`, `sd = 49.51` | `[2583.3, 2830.9]` |
| L-P2 | Mersenne-type lines of `2^p - 1`, `p = 8..400` | `E = 1087.58`, `sd = 30.11` | `[1012.3, 1162.9]` |
| L-P3 | genuine Fermat lines `m = 2^e`, `e = 3..12` (up to 4096-bit values) | `E = 22.91`, `sd = 4.43` | `[11.8, 34.0]` |

**Reading.**
- A pass says the Collatz-line integers behave like the random model after exact small-prime sieving. That is the
  assumption under every heuristic that no sixth Fermat prime exists.
- A failure in either direction would show a Collatz-line bias the model misses, and would weaken that heuristic's
  calibration.
- Probable primes are counted with BPSW (`sympy.isprime`), typed PRP.
