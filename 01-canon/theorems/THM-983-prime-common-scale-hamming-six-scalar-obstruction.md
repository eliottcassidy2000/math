---
id: THM-983
title: Prime common-scale Hamming-six scalar obstruction above seventeen
status: CLAIMED STRUCTURAL — an exact residue-class cardinality formula reduces every prime common scale p>=19 to a six-largest-capacity bound, except for the finite support checks p=23 and p=29; both exceptional 924-support banks are empty in an independent scratch reconstruction, and frozen referees/formalization are in progress
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
`q>=3`.  The primes below 39 are `17,19,23,29,31,37`.  THM-981 handles the
exceptional scale seventeen.  The same table immediately excludes
`19,31,37`.  Only `23` and `29` reach the numerical threshold:

- at `p=23`, an owner must see at least five of the seven high-cardinality
  ratios; the exact twelve-vertex ratio digraph has no six-vertex support with
  that condition at every owner;
- at `p=29`, every owner must see all five high-cardinality ratios.  The
  resulting multiplicative closure cannot fit in six labels; equivalently the
  exact 924-support check is empty.

Thus no prime `p>=19` reaches even the scalar owner gate.  Within THM-860's
finite primitive range this simultaneously removes every prime scale after
seventeen.

## Tournament and carrier audit

The natural pair observable is whether the provider/owner ratio lies in the
high-cardinality residue set.  This produces a directed ratio graph, not a
tournament: reciprocal ratios can tie or have different bonuses.  A switch
can orient ties along the multiplicative label order, but the resulting
tournament forgets the absolute owner-capacity sum.  The faithful carrier is
the labelled six-vertex induced ratio digraph with vertex-wise capacity
thresholds.  Signed projective classes simplify some residue rows but do not
preserve the `p=23` threshold without multiplicities.

Promotion requires a frozen implementation of the residue-class formula,
the complete bonus table, the two exceptional support scans, and an
independently written referee.  The theorem concerns only prime common-scale
Hamming-six faces.  Composite scales, H5 ramification, non-AP/deep sheets, and
global sporadic emptiness remain open.
