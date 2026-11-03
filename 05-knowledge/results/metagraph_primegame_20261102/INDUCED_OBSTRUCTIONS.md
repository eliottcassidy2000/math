# Actual source containment is stronger than a count recurrence

**PROVED all-order reductions with FINITE-EXACT complete core tables.**
All H values below are raw Hamilton-path counts. The unrestricted spectrum
conjecture remains OPEN; every nontrivial target in the obstruction table
is attained.

## Core reduction for a prime source

Let p be an attained odd prime, and let p<q<3p. There exists an induced
embedding from some tournament of count p into some tournament of count q
iff some **strong** tournament of count q has an induced subtournament of
count p.

Proof: SCCs of a tournament form an ordered join and their H counts multiply.
A count-p source therefore has one strong component of count p, with the
other components singleton vertices. Under an induced embedding, that strong
component lies entirely in one target SCC C. Induced path monotonicity gives
H(C)>=p. Also H(C) divides q. Every proper divisor of odd q is at most q/3<p,
so H(C)=q. Conversely, a count-p subset of a strong count-q tournament itself
provides the required source and target. QED.

For q<125 the inherited strong minimum f(10)=125 puts every possible strong
target at order at most9. The complete census and all its subsets therefore
decide this question **at every order**, including arbitrary transitive
padding of the source or target. This does not bound arbitrary intermediate
rewiring: it addresses induced, old-edge-preserving extensions.

## Complete small prime binary branches

For every attained odd prime p<=61 and q=2p+-1, the complete test has these
nontrivial negative pairs:

| p | q | number of strong target classes checked |
|---:|---:|---:|
|13|27|5|
|29|57|7|
|31|63|2|
|43|85|10|
|43|87|11|
|47|93|15|
|53|105|18|

The additional failures3->7 and11->21 have forbidden targets. The formal
p=7 rows of the receipt have an unavailable source and are excluded from
this classification. Every other branch in this declared range has an
explicit **one-vertex** extension witness in INDUCED_BRANCHES.json.

Thus no operation that only adds vertices and retains every old edge can
repair, for example,31->63. Allowing more outgoing ports does repair31->61,
but cannot repair31->63. These are different failure boundaries.

Every nontrivial negative target above is composite. This motivates testing
prime-target constructions separately, but is not an all-order theorem.
Even a successful prime-to-prime binary rule would not settle the spectrum:
prime targets13 and43 have forbidden odd half-predecessors7 and21, and
infinitely many construction seeds cannot simply be assumed available.

## A precise failure of composing count-level arrows

Write a~>b when some count-a tournament embeds inducedly into some count-b
tournament. This is the image of a composable embedding relation under the
map T|->H(T). It is **not transitive**:

    47 ~> 73,       73 ~> 117,       but 47 !~> 117.       (1)

For reproducibility, rows below are outgoing-neighbor bitmasks with labels
0,...,n-1. In both positive examples the induced source is the first n-1
vertices, so there is no search over relabellings in the certificate.

    H=73 on8 vertices:
      [134,252,248,241,225,193,1,64]
      first7 vertices have H=47.

    H=117 on9 vertices:
      [298,508,505,496,481,448,385,1,128]
      first8 vertices have H=73.

All27 strong H117 classes have been checked: none contains a count47
subset. Since117<3*47 and117<125, the core reduction proves the negative
part of(1) at all orders. The two middle H73 witnesses cannot be identified
or substituted: one contains a count47 source and the other does not.

This is an explicit obstruction to composing scalar recurrence arrows.
The actual graph embeddings compose perfectly; the count quotient discards
which middle object they use. A correct grammar must keep a compatible
witness family or provide a legal, costed transport between its middle
representatives. Choosing an arbitrary new representative is not such a
transport. PRIMEGAME's formal-return/executable-return distinction is an
analogous loss with a different exact predicate; it is not an asserted
isomorphism of the two systems.

## Audits and scope

`induced_branches.cpp` computes all subset counts by dynamic programming.
`audit_induced.py` independently enumerates every directed simple path,
without a subset-count recurrence. They check153914 nonempty proper subsets
of all687 strong sources below125, with713827 explicitly enumerated paths.
The audit additionally supplies the two witnesses in(1).

`padded_ears.cpp` tests a separate repair boundary: all cuts on the nine
one-port-dead source counts, with transitive prefix and suffix lengths
summing to at most2. Its297236 cuts add no binary branch or lower-port repair
beyond the unpadded examples. This last assertion is **only FINITE-EXACT**;
the prefix length is a declared bound, not an all-padding theorem for ears.

All-order proof dependencies: inherited SCC factorization, induced path
monotonicity (also proved in TWO_BIRTH_THEOREM.md), and the Busch floor in
Session142 `minimal_gap/PROOF.md`. Exact class-bank pins are in INPUTS.json
and SOURCES.json. No universal realization statement follows from this table.
