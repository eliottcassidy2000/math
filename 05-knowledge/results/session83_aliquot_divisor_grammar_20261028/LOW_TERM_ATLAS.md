# A complete low-term atlas: at most nine proper divisors

Date: 2026-10-28. Status: **PROVED** elementary classification;
**FINITE-EXACT** independent verification. No literature-priority claim.

## Result

The only perfect numbers with **at most nine proper divisors** are
**6,28,496**. The complete abundant/deficient classification in these layers is:

| Number of proper-divisor terms | All source shapes | Abundant sources | Perfect sources |
|---:|---|---|---|
|0|1|none|none|
|1|p|none|none|
|2|p²|none|none|
|3|p³ or pq|none|6|
|4|p⁴|none|none|
|5|p⁵ or p²q|12,18,20|28|
|6|p⁶|none|none|
|7|p⁷ or p³q or pqr|24,40,54,56,70,88,104; and6q for every prime q≥5|none|
|8|p⁸ or p²q²|36,100,196|none|
|9|p⁹ or p⁴q|48,80,112,162,176,208,272,304,368,464|496|

Different prime letters denote distinct primes. Every unlisted source in each
layer is deficient. The exceptional unit1 has no proper divisors, sum0, and
is deficient. This is a theorem about a restricted divisor-count universe,
**not** a classification of all perfect numbers or all abundant numbers.

## 1. Dependencies and what is added here

[DIVISOR_TERM_CLASSES.md](DIVISOR_TERM_CLASSES.md) proves the source-shape
classification, deficiency of prime powers, and the three- and five-term rows.
[SEVEN_TERM_ABUNDANCE.md](SEVEN_TERM_ABUNDANCE.md) proves the complete seven-term
row, including its infinite part and the correction70≠210. The only new
case splits needed to close this atlas are the eight- and nine-term rows.
These are elementary derived corollaries, with no novelty claim.

## 2. Eight terms: the mixed-square layer

The shapes are p⁸ and p²q². Prime powers are deficient. For odd distinct primes
p,q, the largest possible abundancy ratio of p²q² is at p=3,q=5:

`σ(p²q²)/(p²q²) ≤ (1+1/3+1/9)(1+1/5+1/25)=403/225<2`.

Here each factor decreases as its prime increases. If one prime is2, write
n=4q² with q odd prime. Then

`s(4q²)−4q² = −q²+7q+7`.

This is positive for q=3,5,7, giving36,100,196 with excesses19,17,7. For
q≥11 it is negative: it equals−37 at11 and decreases thereafter. No odd prime
lies strictly between7 and11, and none of the tested primes gives zero.
Therefore no eight-term source is perfect, and the abundance list is complete.

## 3. Nine terms: one more two-prime layer

The shapes are p⁹ and p⁴q. Again prime powers are deficient. Put

`A=p⁴+p³+p²+p+1`, `B=p⁴−p³−p²−p−1`.

For p⁴q, its excess is A−Bq.

* For p=2, the excess is `31−q`. The odd primes below31 are
  3,5,7,11,13,17,19,23,29. Multiplying these by16 gives exactly the nine
  abundant sources in the table other than162. At q=31 the source is496 and
  is perfect. Every larger prime gives deficiency.
* For p=3, the excess is `121−41q`. Only q=2 is positive, giving162 with
  excess39. Every other allowed prime q is at least5 and gives deficiency;
  equality has no integral q.
* For p≥5, q≥2 gives the upper bound
  `−p⁴+3(p³+p²+p+1)<−p⁴+3p⁴/(p−1)≤−p⁴/4<0`.

This proves the nine-term row. Combining all rows proves the atlas and the
exact perfect list6,28,496. In particular an odd integer with at most nine
proper divisors is always deficient; no assertion about larger divisor counts
or the existence of odd perfect numbers is made.

## 4. An exact compressed description

Let

```
E = {12,18,20,24,36,40,48,54,56,70,80,88,100,104,
     112,162,176,196,208,272,304,368,464}.
```

The abundant sources with at most nine proper divisors form exactly the
disjoint union **`E ∪ {6q:q≥5 prime}`**. The finite set E has23 members, none
equal to6 times a prime≥5. Consequently their counting function is

`A_≤9(X)=max(π(X/6)−2,0)+#{e∈E:e≤X}`,

and for **X≥464**, it is **`π(X/6)+21`**. This uses only the exact classification,
not a prime-density estimate. It counts sources rather than distinct targets.

The term-count model exposes6,28,496 at layers3,5,9, respectively, but it must
not be extrapolated to a new perfect-number recurrence. The next even perfect
number8128 has13 proper divisors; this is a direct verified example outside
the atlas, not a theorem about all later layers.
The separate [PERFECT_TERM_LAYERS.md](PERFECT_TERM_LAYERS.md) gives a guarded
general perfect-layer result; this finite-layer atlas does not substitute for
its hypotheses or proof.

## 5. Reproduction and stop boundary

Standalone Python3 standard library:

```
python3 check_low_term_atlas.py > low_term_results.json
python3 -O check_low_term_atlas.py > low_term_results.optimized.json
cmp low_term_results.json low_term_results.optimized.json
```

The checker independently sieves proper-divisor sums and counts over every
source1..100000, compares every abundant/perfect source in each layer0..9
against the table, verifies the counting formula at every integer cutoff, and
tests named excluded sources210 and8128. A second divisor-pair oracle checks
36,100,196,162,496 and all named boundary controls. The existing companion
checkers independently verify the factorization-based source shapes and the
seven-term transport. Explicit exceptions remain active under optimization;
generated JSON is a deterministic local receipt, not a dependency.

Global conclusions follow from the all-prime proofs. The computation has an
explicit finite universe and does not promote a finite absence to a general
perfect-number theorem. This atlas stops at nine proper divisors.
