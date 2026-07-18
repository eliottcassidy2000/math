---
id: THM-983
title: Prime common-scale Hamming-six scalar obstruction above seventeen
status: PROVED STRUCTURAL + FINITE-EXACT — the residue-block recurrence and exact B6/B5 bounds exclude every prime common scale p>=19 at scalar owner capacity, with independent exhaustive 924-support checks for the only numerical exceptions p=23 and p=29
source: codex-2026-07-17-S66 prime-scale residue-capacity synthesis
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-981]
related: [THM-962, THM-974, THM-980, THM-982, HYP-6820]
---

# THM-983 — every prime common scale above seventeen is scalar-impossible

Let `p>=19` be prime.  The primitive proper AP-centred common-scale-`p`
Hamming-six sheet bank is empty already at unit-independent scalar owner
capacity.

The only effective orders are `1,p`.  For an order-`p` provider/owner ratio
`r in F_13^*`, unit choice does not change the number of covered sheets.  The
cardinality has the exact residue-class description

```text
a_p(r)=#{x in Z : -p < x <= p and x == p*r (mod 13)}.
```

Write `p=13q+s`, `1<=s<=12`.  Each `a_p(r)` is `2q` plus a bonus depending
only on `(s,r)`.  If `B_6(s)` and `B_5(s)` are the sums of the six and five
largest bonuses, their complete twelve-row table is

```text
s        1  2  3  4  5  6  7  8  9 10 11 12
B_6(s)   1  3  5  6  6  6  7  9 11 12 12 12
B_5(s)   1  3  5  5  5  5  6  8 10 10 10 10.
```

An order-one provider covers all `p` sheets only at its matching owner and
zero sheets at the other five owners.  At any different owner the remaining
five providers cover at most `10q+B_5(s)<13q+s=p`; hence a scalar row must be
all-order-`p`.

For an all-order-`p` row, any owner has capacity at most
`12q+B_6(s)`.  Since `B_6(s)-s<=2`, this is strictly below `p` whenever
`q>=3`.  The primes below 39 are `17,19,23,29,31,37`.  THM-981 classifies the
exceptional scale seventeen.  The same table immediately excludes
`19,31,37`.  Only `23` and `29` reach the numerical threshold:

- at `p=23`, an owner must see at least five of the seven high-cardinality
  ratios; the exact twelve-vertex ratio digraph has no six-vertex support with
  that condition at every owner;
- at `p=29`, every owner must see all five high-cardinality ratios.  The
  resulting multiplicative closure cannot fit in six labels; equivalently the
  exact 924-support check is empty.

Thus no prime `p>=19` reaches even the scalar owner gate.  The proof is uniform
in `p`, not a cutoff scan: adding thirteen to the scale adds two complete
residue blocks, and `B_6(s)-s<=2<q` closes every `q>=3`.  Within THM-860's
finite primitive range this simultaneously removes all 139 primes from 19
through 839.  The primary certificate additionally scans all 129,360 supports
at those primes plus the `p=17` anchor and finds only the two expected
quadratic rows at seventeen.

## Tournament and carrier audit

There are two useful but inequivalent tournament shadows.  Ordering the twelve
ratio weights `a_s(r)` and breaking ties by ratio label gives a transitive
tournament for every residue `s`, useful for checking the bonus order
statistics.  Orienting each reciprocal provider pair by the two directed
capacities instead gives, at both exceptional scales, score histogram
`2,3,3,3,5,5,6,6,8,8,8,9`, 40 directed triangles, one SCC of size twelve,
124,961 Hamiltonian paths, 42 ties, and 12 flips.  Neither shadow retains the
six-support incidence together with absolute owner sums.  The faithful carrier
is therefore the labelled induced ratio-capacity digraph with its weight
vector and owner thresholds; signed projective classes alone do not preserve
the `p=23` condition.

## Certificates

The primary C++ implementation proves the complete-residue-block lemma,
reconstructs the bonus rows by literal integer counts, maximizes all five- and
six-subsets, scans the exceptional banks, and audits every THM-860 prime.  A
separately written Python referee derives the same recurrence and tables,
independently scans `p=17,23,29`, and supplies the reciprocal-capacity
tournament fingerprint.  C++ `-O3`, `-O0`, and ASan/UBSan outputs agree
byte-for-byte; normal and optimized Python outputs agree.

| artifact | SHA-256 |
|---|---|
| `04-computation/lrc13_prime_scale_hamming_six_scalar_obstruction_codex_prime.cpp` | `86b34f9d404e10831ef9a86003dbf9f6fbc062054a9dfa884e20e853947911e8` |
| `05-knowledge/results/lrc13_prime_scale_hamming_six_scalar_obstruction_codex_prime.out` | `87eb95da5e554ab7ffa960ed8a485ca97ae5ffc75b78df7444119810b01620cc` |
| `04-computation/lrc13_prime_scale_hamming_six_scalar_referee_codex_prime.py` | `4ffd04a51c9c4584ec4c2e43806cd22703765f7921af42436ed41718b57c9883` |
| `05-knowledge/results/lrc13_prime_scale_hamming_six_scalar_referee_codex_prime.out` | `2570e3b5b285fb2b901e89113734ffdcbe8b7586865976c2e5b69c9d0bf19112` |

The theorem concerns only prime common-scale Hamming-six faces.  Composite
scales, H5 ramification, non-AP/deep sheets, and global sporadic emptiness
remain open.
