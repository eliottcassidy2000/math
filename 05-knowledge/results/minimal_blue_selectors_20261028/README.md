# Exact minimum for singleton composition; one-port separation and its limit

2026-10-28. **PROVED:** M1--M5 in [THEORY.md](THEORY.md), for the explicitly
specified singleton/native-complement/fresh-re-marking language. The result
settles the previous packet's minimality question for that language.
It does not settle arbitrary strong-atom substitution or named-guard closure.

The smallest behavioral carrier is

    K_min(T) = (endpoint-rooted type of T,
                endpoint-rooted type of native complement(T),
                B(T)).

The minimal blue interface B has exactly4n-8 states for n>=5, compared with
the previous O(n)-bit boundary record. It needs a signed offset, one polarity
bit and the exceptional flat/center/dead cases: O(log n) bits. This bound
does not include the two retained rooted relations.

A four-singleton context recovers the ordered endpoints from an ordinary
unlabelled output class, with no strongness assumption. That proves the two
rooted relations are observable and makes the sufficient carrier minimal.
The all-future quotient at size n is determined by contexts up to order
max(n+4,2n-3). At7 it has exactly32272 cells among32768 presentations.

For the seven-vertex pair1072/4176, **one extra interior port is sufficient
and minimum in added-port arity**. Marking the current second vertex works
after complement. A corresponding start-conditioned reset has23 versus21
paths; the ordinary event "partner H=29" has probabilities1/23 and0. A single
Boolean named legality guard also works: the word(0,1,3,4,2,5,6). Its three
moved positions and span3 are minimal among permutation guards on this pair.
These are separate, explicitly costed notions of weakest selector.

**FINITE-EXACT:** a successor port in both rooted types separates all32768
seven-vertex presentations. **PROVED limitation:** a fixed boundary width w
still fails on a strong pair of order2w+5. In particular, the successor-port
enhancement already fails at9. The proof pads the original pair and extends
its two incompatible endpoint maps by the identity on every added vertex.

[CONNECTIONS.md](CONNECTIONS.md) recovers the published Borwein thresholds
{15,113,80499}, tournament/even-graph equinumerosity, the exact Euler-mark
criterion for excess-minus-one graph carriers, and older oddness results.
The synthesis distinguishes support compression, symmetry-compatible maps,
and averaging over odd fibers. It proposes four concrete next tests, including
the least endpoint-decoding cost and a moving-selector model.

The81/27 control passes unchanged: same67 labelled cyclic triples and462915
Hamilton paths, but176/512 binary attachment orbits and27/0 transporters.
Paley7 retains seven anti-involutions, a centralizer of order3 inside Aut21,
and the explicit Frobenius phase i->2i. The inherited F64/F256 orientation
guards and LCA witness constraints remain necessary for richer interfaces.

## Inheritance and scope

Anchor: minimizing K and separating the7-pair. Niche: odd fibers versus
equivariant choices. Wildcard: support thresholds and marked graph carriers.
The live board is boundary offsets, endpoint decoding, selector costs,
joint transport, symmetry marks, and delayed support crossings.

Closest proved mechanism: [S2--S4](../singleton_pair_reconstruction_20261028/THEORY.md)
and the [strong-core sandwich/reset lemmas](../unnamed_carrier_dissection_20261026/README.md).
Canonical hostiles: the7-pair and81/27. Corrected near miss: a sufficient
boundary profile need not be minimal, and separation through7 after adding
one port need not extend to9. Least-used sidecar: the common polarity of
boundary runs, forced by reflected cross-pairs.

The old [classification and enumeration packet](../singleton_pair_reconstruction_20261028/README.md)
retains routes to4p thresholds, modular atoms, rank charts, LCA, and the
A000568350-prefix/order10 generation work. Its minimality OPEN statement is
superseded here. Counts and runtime bounds for those enumeration jobs are
not extended in this session. LRC(14) remains OPEN.

Read [AUDIT.md](AUDIT.md) for proof gates, finite universes, source pins and
reproduction. No independent-agent audit or global theorem-ID promotion
was performed. The normal and optimized verifiers use active exception
checks rather than removable Python assertions.
