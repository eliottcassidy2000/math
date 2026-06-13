---
source: codex-2026-06-02-S566
status: reflection + exact audit
tags: [LRC, n17, CRT, dyadic, residual-packets, endpoint-debt, Tournament-Analysis]
---

# The n=17 dyadic residual packet surprise

The S562 surprise was phrased this way:

```text
n=17 has a dyadic residual tier even though 2 is not a base factor of 17.
```

That remains true, but S566 sharpens it.  The dyadic tier is not unique to the
`skip 8` packet.  It is a whole skip-spectrum phenomenon.

For every skipped label `s=1..16`, the packet

```text
{1} union {17*q : 1 <= q <= 16, q != s}
```

obeys:

```text
scale 17 -> 34 -> 68
gap halves
boundary debt doubles
gap*boundary is conserved
forbidden length is stable
```

So the dyadic coordinate is being introduced by the endpoint/residual
denominator.  It is not inherited from the base CRT factorization, because
`17` has no dyadic factor.

## The real surprise: two gauges

S559 selected `skip 8` for a good reason:

```text
skip 8 has the smallest visible gap.
```

At scale `17`, the scalar-gap winner is:

```text
skip 8: gap/th=1/272, boundary=450, product=225/136.
```

But the product/debt gauge selects:

```text
skip 6: gap/th=1/204, boundary=330, product=55/34.
```

That is the new hinge.  `skip 8` is the scalar face.  `skip 6` is the
frontier-product face.

This makes the S560 carryover less mysterious.  When the n=17 idea goes back
to `n=14=2*7`, the visible structure picks the predecessor-of-apex `skip 6`
because the apex `7` must stay as a shield.  S566 says that predecessor channel
was already latent in the n=17 product gauge.

## Neighboring primes

The half-gate dyadic mechanism appears for primes:

```text
p=5,7,11,13,17,19,23.
```

For the half-gate skip `(p-1)/2`, the audit finds:

```text
gap0 = 1/(p*(p-1))
gap(h) = gap0/2^h
boundary(h) = boundary0*2^h
```

The product constants are not monotone and do not immediately reveal a simple
closed form:

```text
9/10, 10/7, 108/55, 33/26, 225/136, 340/171, 630/253.
```

That irregularity is probably not noise.  It is where endpoint ownership or
residue arithmetic is hiding.

## Breakers do not add a second mystery

For `n=17`, `skip=8`, and raw breaker `r=1..16`, even breakers normalize down
to lower dyadic depth.  All rows keep product `225/136`.

The depth histograms are exactly what the 2-adic valuation predicts:

```text
h=0: {0:16}
h=1: {0:8, 1:8}
h=2: {0:4, 1:4, 2:8}
h=3: {0:2, 1:2, 2:4, 3:8}
h=4: {0:1, 1:1, 2:2, 3:4, 4:8}
h=5: {1:1, 2:1, 3:2, 4:4, 5:8}
```

The odd breakers are the genuinely primitive new rows at each dyadic depth.

## Tournament Analysis

Vertices:

```text
skipped quotient labels 1..16, with whole packet data attached.
```

This deliberately challenges the default runner-vertex assumption.  Alternative
vertices considered were raw runners, skipped labels, endpoint residues,
boundary owners, dyadic depth events, and proof-obligation packet states.  The
skip-label quotient preserves the positive-gap LRC predicate and the dyadic
product ledger; it destroys endpoint ownership, witness location, and
wall-crossing order.  That is why this tournament records gauge disagreement
but not the cyclic obstruction itself.

Observable:

```text
(gap/th, conserved product, boundary count, skipped label)
```

Two switches:

```text
scalar-first:  smaller gap, then smaller product
product-first: smaller product, then smaller gap
```

Both tournaments are transitive:

```text
score_hist={0:1,1:1,...,15:1}
directed_3_cycles=0
sccs=singletons
```

But the top vertex changes:

```text
scalar top = 8
product top = 6
edge flips between gauges = 15
```

This is exactly the kind of ranker/ledger split the repository has been
circling.  The cyclic object is not in the bare skip spectrum.  It should enter
when product debt is decorated by endpoint owners, pressure leaves, or
cross-prime signs.

The incoming HYP-2075/HYP-2076 multi-sieve work gives this a sharper job.  The
pair-sum sieve wants fast witness extraction; the recursive lower-bound ledger
wants stable product debt.  The n=17 skip spectrum is a small place where those
two requirements visibly choose different labels.

## What I now believe

The phrase "n=17 dyadic surprise" should now mean:

```text
prime residual packets acquire dyadic valuation tiers after the gate,
and different proof gauges select different skipped labels.
```

The scalar gauge sees the half-gate.  The product gauge sees the predecessor.
The n=14 transfer is the moment when the predecessor gauge becomes visible in
the actual composite proof route.

**Artifacts:** `04-computation/lrc_n17_dyadic_packet_surprise_s566.py`,
`05-knowledge/results/lrc_n17_dyadic_packet_surprise_s566.out`, HYP-2079.
