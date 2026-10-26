# The Legendre tournament on primes congruent to 3 modulo 4

2026-10-25, klein. **EXACT** definition (quadratic reciprocity); **FINITE-EXACT** census
(`qr_tournament.py`, receipts `qr-normal.jsonl` = `qr-optimized.jsonl`); quasirandomness is
**HYP-9375** (OPEN). Inheritance: Paley tournaments THM-640 (prime scope only, MISTAKE-232), the
QR-decomposition of A000568 modulo `p` in THM-307, and MISTAKE-214 (node repetition is not score
repetition).

## Definition

For primes `p != q` both `≡ 3 (mod 4)`, quadratic reciprocity gives `(p/q) = -(q/p)`; put `p -> q`
iff `(p/q) = +1`. This is a tournament `L` on `P_3 = {3, 7, 11, 19, 23, 31, 43, 47, ...}`. On pairs
with a prime `≡ 1 (mod 4)` the relation is symmetric and is a graph, not a tournament.

## Census

| universe | triangle density | TT4 | source+C3 | C3+sink | strong |
|---|---:|---:|---:|---:|---:|
| random tournament | 1/4 | 3/8 | 1/8 | 1/8 | 3/8 |
| first 120 primes of `P_3` (exact, all `C(120,4)` 4-sets) | 0.25339 | 0.3678 | 0.1254 | 0.1255 | 0.3813 |
| first 2,000 primes (300,000 sampled 4-sets) | 0.2507 | 0.3750 | 0.1240 | 0.1252 | 0.3759 |
| first 20,000 primes (300,000 sampled 4-sets) | 0.2497 | 0.3764 | 0.1254 | 0.1243 | 0.3739 |

All four isomorphism classes of 4-vertex tournaments occur; scores on the first 120 primes range
from 47 to 68 (mean 59.5), so `L` is not regular (unlike Paley tournaments). The displayed sampled frequencies
are within0.002 of the random values; this is not a certified confidence interval.

## What quadratic reciprocity does and does not say

It says exactly that `L` is a tournament (antisymmetry on `P_3`) and that the relation on the
remaining pairs is a graph. It does not say anything about `A000568(4) = 4`; the user's sentence
"reciprocity describes that a 4-vertex tournament has four isomorphism classes" is a numeral
coincidence (four residue cases versus four orbit classes), recorded here so that it is not
re-promoted. The exact content salvaged is the census above.

## Recovered stronger structural theorem

The [known CRT/Dirichlet extension theorem](../root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md)
shows that L is the countable random tournament, hence contains every finite
tournament. This goes beyond the finite census but does not settle its
increasing-prime-prefix frequencies. The [2026-10-26 continuation](../orientation_monodromy_rule30_20261026/README.md)
gives a graded exchange algebra and a four-object converse-equivariant
bijection, explicitly retaining their losses and scope.

## HYP-9375 (OPEN)

`L` is quasirandom in the sense of Chung–Graham for tournaments: for every fixed `k`, the induced
`k`-vertex isomorphism classes occur with the random frequencies `|Aut|^{-1} k!/2^{C(k,2)}` as the
prime bound grows. Evidence: the table. Likely route: large-sieve bounds for sums of products of
Legendre symbols over primes; no literature search beyond one web query found a statement either way.
A refutation would be a persistent bias in some `k`-class; the `k = 3, 4` data show none.
