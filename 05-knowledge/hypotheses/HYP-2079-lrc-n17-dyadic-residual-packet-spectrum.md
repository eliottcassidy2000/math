---
id: HYP-2079
status: SUPPORTED - exact S566 audit; proof formula open
source: codex-2026-06-02-S566
related:
  - HYP-1844
  - HYP-1866
  - HYP-1868
  - HYP-2064
  - HYP-2069
  - HYP-2073
  - HYP-2075
  - HYP-2076
  - THM-369
  - THM-391
---

# HYP-2079: n=17's dyadic surprise is a full skip-spectrum lift, with scalar/product gauge split

## Statement

For the prime LRC denominator `n=17`, the S562 dyadic residual tier is not an
accident of the `skip 8` half-gate packet.  Every one-gate packet

```text
{1} union {17*q : 1 <= q <= 16, q != skip}
```

has an exact dyadic lift law:

```text
gap_ratio(h, skip) = gap_ratio(0, skip) / 2^h
boundary(h, skip)  = boundary(0, skip)  * 2^h
product(h, skip)   = gap_ratio(h, skip) * boundary(h, skip)
                   = product(0, skip).
```

The half-gate packet `skip=8` is selected by the scalar visible-gap gauge, but
the lower conserved product belongs to `skip=6`.  Thus the n=17 surprise points
back to the n=14 predecessor-of-apex channel in a product/debt gauge, even
though the scalar gap initially selected the half-gate.

## Evidence

`lrc_n17_dyadic_packet_surprise_s566.py` audits all skips `1..16` through
depths `h=0,1,2` at scale `17*2^h`.  Every skip conserved product.

Key rows:

```text
skip 6:
  gap0 = 1/204, boundary0 = 330, product = 55/34

skip 8:
  gap0 = 1/272, boundary0 = 450, product = 225/136

skip 16:
  gap0 = 1/255, boundary0 = 450, product = 30/17
```

The scalar-gap ordering picks `skip=8`:

```text
scalar-gap winner: skip 8
```

The conserved product ordering picks `skip=6`:

```text
product-debt winner: skip 6
```

There are `15` edge flips between the scalar-first and product-first
tournament gauges on the sixteen skipped labels.

## Prime Control

The same half-gate dyadic lift law appears in neighboring odd primes:

```text
p=5,  skip=2:  gap0=1/20,  boundary0=18,   product=9/10
p=7,  skip=3:  gap0=1/42,  boundary0=60,   product=10/7
p=11, skip=5:  gap0=1/110, boundary0=216,  product=108/55
p=13, skip=6:  gap0=1/156, boundary0=198,  product=33/26
p=17, skip=8:  gap0=1/272, boundary0=450,  product=225/136
p=19, skip=9:  gap0=1/342, boundary0=680,  product=340/171
p=23, skip=11: gap0=1/506, boundary0=1260, product=630/253
```

For these half-gate packets, the audited formula is:

```text
gap0 = 1 / (p*(p-1)).
```

Boundary and product constants remain arithmetically irregular in the audit;
they are the next target for a closed form.

## Breaker Parity Collapse

For the `n=17`, `skip=8` family with raw breaker `r in {1,...,16}`, even
breakers do not create new packets.  Normalization drops them to lower dyadic
depth:

```text
h=0: {0:16}
h=1: {0:8, 1:8}
h=2: {0:4, 1:4, 2:8}
h=3: {0:2, 1:2, 2:4, 3:8}
h=4: {0:1, 1:1, 2:2, 3:4, 4:8}
h=5: {1:1, 2:1, 3:2, 4:4, 5:8}
```

All rows keep product `225/136`.  The primitive new packets at each depth are
the odd-breaker rows.

## Interpretation

S559 found `skip=8` because it ranked by visible gap.  S560 found `n=14` wants
`skip=6`, the predecessor of the apex `7`, because the apex must stay as a
shield.  S566 shows these are not contradictory:

```text
skip 8 = scalar-gap face of n=17
skip 6 = lower-product/debt face of n=17
```

So the n=17 prime control already contains both messages.  The scalar view
points to the half-gate.  The recursive frontier-product view points to the
predecessor channel that later becomes visible in `14=2*7`.

Incoming HYP-2075/HYP-2076 sharpen the role of this split.  If pair-sum moduli
are the apex-free multi-sieve primitive, then the n=17 skip spectrum is a small
laboratory for choosing which packet gauge should feed the recursive
lower-bound ledger: scalar witness extraction sees `8`, while conserved
frontier debt sees `6`.

## Tournament Analysis

Vertices were skipped quotient labels `q=1..16`, each carrying the whole
`n=17` residual packet.

Assumption challenged: the tournament vertices need not be runners.  I tested
the runner-level instinct against skipped labels, dyadic depth events, endpoint
owners, boundary residues, and proof-obligation packet states.  The selected
skip-label quotient preserves the LRC predicate "this packet has a positive
lonely gap and a stable dyadic lift law" together with the scalar/product debt
ledger.  It destroys endpoint ownership, witness location, wall-crossing order,
and the signs that would be needed for a cyclic obstruction.

Pairwise observable:

```text
(gap/th at scale 17, conserved gap*boundary, boundary, skip)
```

Two switches were compared:

```text
scalar switch:  smaller gap wins, then smaller product
product switch: smaller product wins, then smaller gap
```

Both tournaments are transitive with no directed 3-cycles and singleton SCCs,
but their top vertices differ:

```text
scalar top = skip 8
product top = skip 6
edge flips = 15
```

Thus the current quotient is a ledger, not a cyclic obstruction.  The
important signal is the gauge disagreement.

## Open Work

1. Prove the dyadic lift formula for all odd prime `p` and all skipped labels:
   gap halves, boundary doubles, forbidden length is stable.
2. Find the closed form for `boundary0(p, skip)` and for the half-gate product
   constants.
3. Re-run the n=14 transfer using product-first ranking; predict whether the
   predecessor-of-apex channel appears before the explicit apex shield is added.
4. Lift HYP-2073's product-building norm so the scalar/product gauge split is
   weighted by endpoint owners, not raw boundary count only.

## Files

- `04-computation/lrc_n17_dyadic_packet_surprise_s566.py`
- `05-knowledge/results/lrc_n17_dyadic_packet_surprise_s566.out`
- `07-reflections/lrc-n17-dyadic-residual-packet-surprise-s566.md`
