# Endpoint 591: a free layer with a sharp fragile core

**FINITE-EXACT reflection on THM-4314. LRC(14) remains open.**

## Anchor

The THM-4313 exchange was designed only to repair endpoint 592 while
preserving 467 rows. Its unchanged 3,925-mask carrier unexpectedly closes the
next 13 endpoint-591 rows with no additional mask. Two independent complete
programs scan 185,992,950 bodies and find zero failures, advancing the typed
frontier to endpoint 590.

## Niche: low multiplicity, not only failures

The most informative object was not another response cover but the complete
low-witness quotient. It exposes a very small fragile core:

```text
two one-witness obligations;
26 two-witness obligations;
23 distinct two-witness pairs;
all other obligations have at least three witnesses.
```

That quotient converts immediately into an exact theorem for every deletion
set of size at most two. The two singleton-critical witnesses are joint masks,
which explains why the earlier protected-joint census first looked more
robust: retaining the joint sidecar hides precisely the two arbitrary-single
obstructions.

## Wildcard: obstruction hypergraphs by deletion radius

For a fixed layer, witness sets form a hypergraph on carrier masks. Deletion
robustness through radius `k` depends only on hyperedges of size at most `k`.
This suggests a cleaner recursive object than repeated singleton scans:

1. enumerate the low-witness hypergraph through a chosen radius;
2. quotient repeated hyperedges and record joint/nonjoint types;
3. derive all unsafe deletion sets of that radius by upward closure;
4. optimize the next exchange against those exact obstructions;
5. retain a hostile all-witness check so protected-sidecar claims are not
   misreported as arbitrary-deletion claims.

At endpoint 591, radius two is tiny enough to classify completely. Endpoint
590 should first be tested unchanged; if it fails, the low-witness hypergraph
offers the right cost model for additions and deletions. The central missing
bridge remains semantic: a long finite carrier staircase is not a physical
lonely-runner entry theorem.
