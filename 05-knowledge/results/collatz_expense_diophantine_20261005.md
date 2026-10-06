# The expense is Diophantine: weak resets are the upper convergents of log2 3, the exact segment-level tail, and what a bank buys

2026-10-05, opus session `opus-2026-10-05-S9` (expense-Diophantine). Owner's
seed: "work next admissible targets as they emerge". The targets left by the
excess-rate family note were (i) a segment-adaptive deadline, (ii) the decay
law of the expense tail, (iii) the family relative to a bank of certified
excursions. Nothing new arrived from Codex between the two sessions.

**Status.** PROVED (elementary + CITED classical one-sided best approximation):
Theorem Q (the expense depends only on the segment type, its records are the
upper best approximations of `log_2 3`), Theorem T (the exact segment-level
tail as a cylinder series), Proposition B (the two-tier and bank-relative
families and their deadline). FINITE-EXACT: the DP over first-descent types
to length 200, the type frequencies against the census, the record list to
length 320. VERIFIED: the census below `2^20` (float excess, exact at the
reported types). **OPEN:** universal entry. Not a canon promotion.

## 0. Answers

1. **The expense is a Diophantine quantity.** The least admissible rate of a
   first-descent segment depends only on its type `(l, A)`:
   `q(l,A) = ceil(l/(A - l log_2 3))`. For each length the worst type is the
   minimal valuation `A_min(l) = ceil(l log_2 3)`, and the lengths at which
   `q_max(l) = q(l, A_min(l))` sets a record are exactly the denominators of
   the best upper rational approximations of `log_2 3`:
   `2/1, 5/3, 8/5, 27/17, 46/29, 65/41, 149/94, 233/147, 317/200, 401/253,
   485/306`, with `q_max = 3, 13, 67, 306, 804, 2480, 6951, 13984, 26668,
   56382, 207489`. The hardest source below `2^20` (`q = 6951`, a 94-step
   segment) is the `149/94` approximation. The staircase of the coverage
   curve is the continued fraction of `log_2 3`.
2. **The decay law is exact and is a staircase.** The proportion of odd
   sources whose own first segment has rate above `Q` is
   `sum N(l,A)/2^A` over first-descent types with `q(l,A) > Q`, where
   `N(l,A)` counts words of length `l` and valuation `A` whose proper prefixes
   all rise. The DP gives `0.630, 0.259, 0.187, 0.082, 0.055, 0.0070,
   0.0065, 0.00093, 0.00078, 0.00018` at `Q = 2, 3, 7, 13, 31, 67, 104, 310,
   800, 2000`; the census below `2^20` agrees to four decimals, and the type
   frequencies with `A <= 16` agree to `1.4e-6`. The steps sit at the record
   types: the `8/5` type `(5,8)` alone carries `2.7%` of all sources (seven
   words), the `27/17` type `0.23%`, `46/29` `0.054%`, `65/41` `0.017%`,
   `149/94` `2.6e-6`. The decay is a stretched exponential in `Q`: the hard
   types at scale `Q` have length about `sqrt(Q)` (successive denominators
   multiply), and each costs `2^(-(1-h*) A)` with `A ~ 1.585 l`.
3. **Inheritance is the dominant effect.** Below `2^20`, 55% of sources have
   chain rate above 30, but only 5.5% have a hard own segment; the rest
   inherit through running minima. The most popular hard excursions are
   `47 -> 23` (34 steps, `q = 31`, a running minimum of 9.5% of all sources),
   `31 -> 23` (35 steps, `q = 67`, 9.3%), `91` (28 steps, `q = 46`, 8.8%), and
   then the `(5,8)` type at `95, 379, 847, 455, 335, 1243, 1711` (each `q = 67`).
4. **What a bank buys is coverage, not deadline.** With a table of all odd
   sources below `Y` and rate `q_hi` above it, the deadline is still
   `floor(1.051 q_hi log_2 n) + max tau(odd < Y')`, because the last
   admissible segment may land far below `Y'`; the coverage among sources at
   or above `Y = 2^16` is `49%` at `q_hi = 3`, `62%` at 8, `80%` at 16, `86%`
   at 32. A bank of the ten most popular hard excursions raises rate-3
   coverage from `6.1%` to `9.6%`, a thousand to `16.6%`, ten thousand to
   `23.3%`. The segment-adaptive deadline of target (i) is therefore not
   available as a source-only quantity: it is the chain rate itself.

## 1. Inheritance and board

Closest mechanisms: the excess-rate family `F_q`
([weak-reset family](collatz_weak_reset_family_20261005.md), Theorem D and
Proposition R), the rising-prefix count and first-descent census of the
measurement-independence note, the Beatty structure of primitive rising cones
(lengths with `{l log_2 3} < log_2(3/2)`), and, in the Collatz thread's
earlier notes, the role of the convergents of `log_2 3` in the torsion census
(`2^11` against `3^7`) and in the mediant tree of rational cycles. Canonical
hostile: the `(5,8)` type, seven words of density `7/256`, every one a
near-cycle of the positive rational point `c_w/13`. Corrected near miss: a
two-tier deadline written with `log_2(n/Y')`; the drop telescopes to the full
`log_2 n` because the last admissible segment can land anywhere below the
threshold. Least-used sidecar: the count `N(l, A)` of first-descent words,
which turns the empirical staircase into an exact series. Board: **segment
type (l,A) / minimal valuation / upper semiconvergent / first-descent word
count / inheritance through running minima / bank threshold.**

## 2. Theorem Q: the expense and the continued fraction of log2 3

**Theorem Q (PROVED).** (a) For a first-descent segment of type `(l, A)` the
least admissible rate is `q(l,A) = ceil(l/(A - l log_2 3))`; it depends on the
word only through its type, and admissibility is monotone in `q`. (b) At fixed
length the largest rate is attained at the minimal descending valuation
`A_min(l) = ceil(l log_2 3)` (the bit length of `3^l`). (c) The lengths at
which `q_max(l) = ceil(l/(A_min(l) - l log_2 3))` exceeds all earlier values
are the denominators of the best upper approximations `p/l > log_2 3`, i.e.
the upper semiconvergents of its continued fraction.

*Proof.* (a) `2^(qA-l) >= 3^(ql)` is `q(A - l log_2 3) >= l`; the excess is
positive and irrational, so the least integer `q` is the ceiling, and the
inequality persists for larger `q` because `2^A > 3^l`. (b) Larger `A` has
larger excess at the same `l`. (c) `q_max(l)` is a record iff
`A_min(l) - l log_2 3` is a record minimum among `l' <= l`, i.e. `A_min(l)/l`
is a best upper approximation of `log_2 3` in the one-sided sense; by the
classical theorem these are the intermediate fractions of the continued
fraction on the upper side (CITED: Khinchin, Continued Fractions, the theory of
best one-sided approximations). The program verifies the identity of the two
lists for `l <= 320` by exact integer comparison of `2^A/3^l`. QED.

Table of records (`2^A - 3^l` is the denominator of the positive rational
cycle point of the shadow theorem's dual):

| `l` | `A_min` | `q_max` | `2^A - 3^l` | fraction |
|---:|---:|---:|---:|---:|
| 1 | 2 | 3 | 1 | 2/1 |
| 3 | 5 | 13 | 5 | 5/3 |
| 5 | 8 | 67 | 13 | 8/5 |
| 17 | 27 | 306 | 5,077,565 | 27/17 |
| 29 | 46 | 804 | `1.7e12` | 46/29 |
| 41 | 65 | 2,480 | `4.2e17` | 65/41 |
| 94 | 149 | 6,951 | `6.7e42` | 149/94 |
| 147 | 233 | 13,984 | `1.0e68` | 233/147 |
| 200 | 317 | 26,668 | `1.4e93` | 317/200 |
| 253 | 401 | 56,382 | `1.6e118` | 401/253 |
| 306 | 485 | 207,489 | `1.0e143` | 485/306 |

For `l = 94` the exact integer confirmation of `q = 6951` is performed; beyond
`qA ~ 2e6` bits the values are 60-digit evaluations of the ceiling, which are
exact unless `l/(A - l log_2 3)` lies within `10^-50` of an integer.

## 3. Theorem T: the exact segment-level tail

Let `N(l, A)` be the number of valuation words `(a_1, ..., a_l)` with
`2^(A_j) < 3^j` for every `j < l` (rising proper prefixes) and `2^A > 3^l`
(first descent at step `l`).

**Theorem T (PROVED).** The natural density, among odd integers, of the
sources whose own first-descent segment has type `(l, A)` is `N(l,A)/2^A`,
and the density of those with own rate above `Q` is

    sum_{(l,A): q(l,A) > Q} N(l,A) / 2^A.

*Proof.* A source follows a word of total valuation `A` iff it lies in a
residue class modulo `2^(A+1)` (the shadow theorem's forward half), which has
density `2^-A` among odd integers; it descends at the end of the word iff it
exceeds the positive rational cycle point of the word (Proposition R), which
excludes at most finitely many sources per word, of density zero; and the
first descent is at step `l` iff the proper prefixes rise. QED.

The DP enumerates all types with `l <= 200` and terminal valuation at most
`A_min + 40` (8,200 types) with total density `1 - 5.5e-8`. Values:

| `Q` | 2 | 3 | 7 | 8 | 13 | 31 | 67 | 104 | 310 | 800 | 2000 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| exact tail | 0.63006 | 0.25903 | 0.18743 | 0.17631 | 0.08182 | 0.05549 | 0.00696 | 0.00646 | 0.00093 | 0.00078 | 0.00018 |
| census `2^20` | 0.63009 | 0.25900 | 0.18711 | 0.17598 | 0.08138 | 0.05526 | 0.00686 | 0.00643 | 0.00088 | 0.00072 | 0.00015 |

The record types carry the steps: `N(5,8) = 7`, `N(17,27) = 312,455`,
`N(29,46) = 3.8e10`, `N(41,65) = 6.1e15`, `N(94,149) = 1.9e39`, with
densities `2.7e-2, 2.3e-3, 5.4e-4, 1.7e-4, 2.6e-6`. Between records the tail
is flat; the law is a staircase whose steps are at `q_max` of the successive
upper semiconvergents and whose heights decay like `2^(-(1-h*) A_min(l_k))`
with `l_k` the denominators, i.e. roughly `exp(-c sqrt(Q))` since successive
denominators multiply to about `q_max`.

## 4. Inheritance and Proposition B

Every running minimum `y` of `n` is a start of `n`'s chain, so `q(n) >= q(y)`
and a popular hard excursion is inherited by every source whose chain passes
through it. Below `2^20`: 289,666 sources (55.3%) have `q(n) > 30`, of which
29,072 have a hard own segment. The twelve most inherited hard starts:

| start | type | `q` | inheritors | share |
|---:|---:|---:|---:|---:|
| 47 | (34,55) | 31 | 49,927 | 0.095 |
| 31 | (35,56) | 67 | 48,898 | 0.093 |
| 91 | (28,45) | 46 | 45,951 | 0.088 |
| 95 | (5,8) | 67 | 27,694 | 0.053 |
| 379 | (5,8) | 67 | 11,384 | 0.022 |
| 103 | (26,42) | 33 | 11,201 | 0.021 |
| 847 | (5,8) | 67 | 6,664 | 0.013 |
| 455 | (5,8) | 67 | 5,349 | 0.010 |
| 71 | (32,51) | 114 | 4,479 | 0.009 |
| 335, 1243, 1711 | (5,8) | 67 | 4,025 / 3,766 / 3,321 | 0.008 / 0.007 / 0.006 |

**Proposition B (two-tier and bank-relative families; PROVED).** Fix a
threshold `Y` and a rate `q_hi`, and let `Y' = max(Y, 10 q_hi)`. Let
`F(q_hi, Y)` be the sources whose chain segments with start at least `Y` are
`q_hi`-admissible (segments starting below `Y` are unconstrained, their odd
steps being bounded by the finite table of sources below `Y'`). Then every
member satisfies `tau(n) <= floor(1.051 q_hi log_2 n) + max{tau(y): y odd,
y < Y'}`, hence the compiled floor of the weak-reset note. The same holds when
the table is replaced by a finite bank of certified starts at which the chain
stops. *Proof.* Lemma D1 of the weak-reset note applies to every segment with
start at least `Y'`; their lengths telescope against the full drop from
`log_2 n`, because the last such segment may land anywhere below `Y'`; the
remainder is the table. QED. (Checked for every member below `2^20` at `Y`
in `2^6..2^16` and `q_hi` in `3, 8, 16, 32`.)

Coverage among sources at or above `Y`:

| `Y` | `q_hi = 3` | 8 | 16 | 32 |
|---:|---:|---:|---:|---:|
| `2^6` | 0.112 | 0.216 | 0.474 | 0.571 |
| `2^10` | 0.204 | 0.346 | 0.605 | 0.700 |
| `2^14` | 0.373 | 0.515 | 0.731 | 0.804 |
| `2^16` | 0.494 | 0.623 | 0.801 | 0.857 |

Bank of the most inherited hard excursions, all sources below `2^20`:

| bank size | rate 3 | rate 8 |
|---:|---:|---:|
| 0 | 0.061 | 0.149 |
| 10 | 0.096 | 0.227 |
| 100 | 0.123 | 0.272 |
| 1,000 | 0.166 | 0.331 |
| 10,000 | 0.233 | 0.406 |

The returns are logarithmic in the bank size: hard excursions are not rare
exceptions but a positive-density population (Theorem T), renewed at every
scale by the `(5,8)` and longer record types.

## 5. What this settles and what it leaves

- Target (ii) is settled: the tail is an exact cylinder series, a staircase
  on the upper semiconvergents of `log_2 3`. Weak resets are not an accident
  of small numbers; they are the Diophantine shadow of `2^A` being barely
  above `3^l`, and their density at the `k`-th record is about
  `2^(-0.079 l_k)`.
- Target (iii) is settled quantitatively: a bank buys coverage with
  logarithmic returns and buys no deadline. The two-tier family is the right
  object when a finite table is accepted.
- Target (i) is settled negatively: a segment-adaptive deadline is the chain
  rate, which is not source-only; the bank threshold is the only adaptive
  parameter available, and it does not shorten the deadline.
- The frontier therefore moves to the inheritance structure: the proportion
  of sources whose chain avoids every hard excursion. Since the record types
  recur at every scale with positive density, no family with a bounded rate
  and a finite bank has coverage tending to one as the scale grows; the exact
  asymptotic coverage of `F(q_hi, Y)` as `n -> infinity` with `Y` fixed is a
  renewal question on the chain of running minima, left open here with the
  census as its finite evidence.

## 6. Reproduction

[Script](../../04-computation/experiments/collatz_expense_diophantine_20261005.py),
[output](collatz_expense_diophantine_20261005.out),
[JSON](collatz_expense_diophantine_20261005.json):

```text
python3 04-computation/experiments/collatz_expense_diophantine_20261005.py --census-bits 20 --json 05-knowledge/results/collatz_expense_diophantine_20261005.json
python3 -O 04-computation/experiments/collatz_expense_diophantine_20261005.py --census-bits 16
```

Exact integer arithmetic for the record identity (`l <= 320`), the DP counts
and the type densities; 60-digit ceilings with exact confirmation up to
`qA = 2e6` bits; numba for the census, inheritance counts and bank coverage.
Hostiles: the `(l, A)` formula is checked against exact bisection for
`l <= 120`; the two-tier deadline is checked at every member; the census type
frequencies must match the DP to `1e-4` on fully represented cylinders.
