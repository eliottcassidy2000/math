# Endpoint Half-Filter Enumeration, S671

The nicest thing that happened this session is that the old three-state
automaton finally found a very concrete job.

For a tournament `T`, delete each vertex and remember only:

```text
the class of T-v,
and whether v has low, middle, or high outdegree.
```

Low/middle/high means below, equal to, or above the half-score threshold.  This
is not the exact deleted score.  It is just the upper/lower filter side plus a
middle fixed state.

That tiny trace separates every tournament isomorphism class through `n=8`.

This is exactly HYP-2245 becoming algorithmic.  The incident edges of a new
endpoint form a Boolean cube, so there are literal principal filters before
quotienting.  The quotient leaks unless we carry an address.  The address here
is almost absurdly small:

```text
(deleted card, L/M/U side)
```

It also gives a credible A000568 enumerator architecture.  Generate children
from parent classes by `Aut(parent)`-orbits of incident bit patterns.  Then
compute the half-filter deletion trace of each raw child.  If a trace bucket is
pure, keep one representative and skip full child canonicalization.  Through
`7 -> 8`, every half-filter bucket is pure.

The raw scale numbers are good.  To build `n=8`, endpoint orbits give `54,256`
candidates instead of `2,097,152` fixed-path tilings or `268,435,456` labelled
tournaments.  The `n=8 -> 9` estimate is `1,716,608` endpoint-orbit candidates
for the known `191,536` classes.  That is still a lot, but it is a sane lot,
and the trace suggests most expensive canonicalization can be delayed or
avoided.

The contrast is sharp.  The cheap scalar `(score histogram, c3, SCC)` has only
`167` buckets for `6880` classes at `n=8`, with a largest bucket of `577`.
The unpaired deletion-card multiset is almost enough but has two collision
buckets at `n=8`.  The three-state side bit kills those collisions.

The collision anatomy is cleaner than I expected.  Both bad `n=8` card-deck
buckets are regular strong tournaments with identical scalar payloads: same
score histogram, same `c3=20`, and same Hamiltonian-path count inside each
pair (`633` for one pair, `645` for the other).  The missing information is not
more scalar weight.  It is which side of the half-filter owns the deleted
card.  The repair is literally an `L/U` address coordinate.

So the speculative slogan is:

```text
A000568 may be recursively enumerable by feasible half-filter decks.
```

Not proved.  But through `n=8`, it is behaving like the right sufficient
statistic.

The next move is obvious and useful: inspect the two `n=8` card-deck collision
pairs and see why the `L/M/U` side separates them.  That is the local lemma.
If that mechanism generalizes, the half-filter trace becomes a reconstruction
principle rather than a lucky invariant.

For LRC14, the translation is also direct.  The owner/carry problem has been
begging for a small state automaton.  Here is one: owner below threshold,
owner exactly balanced, owner above threshold.  The same left/middle/right
state can label cover arcs, deletion owners, and carry branches before we ask
for exact numeric scores.
