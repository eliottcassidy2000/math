# Independent review of the scaled-pair and one-subdivision packet

2026-10-19. **PASS.** This staged packet is separate from the frozen two-subdivision/forest packet. No mathematical checkout or Git writes were performed by this reviewer.

## Main conclusions checked

- A single common square scaling of any finite positive integer seed permits the simultaneous addition of one independent new vertex for each selected pair, with exactly that pair as its old neighborhood.
- Every finite two-degenerate graph consequently embeds induced in the positive square-sum graph. The construction has log_2 height at most (16/3)(4^(n-1)-1), rather than a claimed polynomial bound.
- Simultaneous free placement gives the exact one-subdivision of every n-vertex graph at height below 4096n^4.
- The universal square-scaled prescribed-neighborhood theorem fails for three neighbors: the particular seed {1,2,3} is excluded by its inherited gap-one/congruent-number obstruction. This does not exclude any graph from a different placement.
- Countably infinite two-degenerate universality is not proved by iteration: previously assigned labels change under further scaling.

The coordinating reviewer independently confirmed the main proof in [review-root/PROOF_REVIEW.md](review-root/PROOF_REVIEW.md).

## Arithmetic checks that matter

The common scale 4lambda, with lambda a power of four, is an actual integer square; therefore old nonedges as well as edges persist. Scaled originals being 0 modulo four alone would not establish this. New vertices are 1 modulo four, which does establish all new-new nonedges.

For a putative third original incidence, the square root has the same odd parity as the designated root, so the root difference t is even. From lambda|t(t-2), the valid conclusion is t=0 or2 modulo lambda/2, not modulo lambda. The selected scale and the bound |t|<5H leave only the intended representatives despite this loss of one power of two.

Distinct difference classes have root-square separation at least 3lambda^2-2lambda, exceeding all changes in the subtracted seed term. Equal difference classes have distinct subtracted seed labels. New labels exceed the scaled original maximum.

The fictitious endpoint H+1 handles one earlier neighbor without inserting an actual extra vertex. Its distance from a genuine base label is at most H, so neither the incidence estimate nor the height estimate changes. A zero-neighbor vertex labeled 2 is safe modulo four.

For the optimized one-subdivision formula, t is bounded by 4n, while the root-class modulus is at least 8n. The height constant 4096 is valid, including original labels and the trivial n=1 case.

## Independent constructions and scope

The separately derived dyadic/mod-seven proof in [provenance/ONE_SUBDIVISION_DYADIC_GUARD.md](provenance/ONE_SUBDIVISION_DYADIC_GUARD.md) also passes. Its original family is 3 modulo seven and its new family 5 modulo seven; the power-of-two condition selects the two correct endpoints. Its 50176n^4 constant is superseded by the square-scaling formula, and it is retained only as independent provenance.

The odd-prime construction in [ODD_PRIME_ONE_SUBDIVISION.md](ODD_PRIME_ONE_SUBDIVISION.md) passes for any fixed odd prime p with 2 a quadratic nonresidue. Its odd power p^a makes original-pair valuations odd, while new-pair sums are 2 modulo p. Since p does not divide both t and t-2, the entire p^a divides one, enforcing the endpoints. It gives O_p(n^4), including the independent p=3 construction with constant 26244.

The optional minimal-mod-seven note is correct only for the stated three-distinct-residue-family pattern over prime fields: p=2,3,5 fail and p=7 has the unique unordered triple {3,5,6}. It is not a lower bound for the one-subdivision method or for arbitrary composite-modulus constructions.

## Exact controls

All scripts are Python standard-library programs. Their explicit checks remain active under -O, and each normal receipt agrees byte-for-byte with its optimized receipt.

1. `verify_one_subdivision.py`: **6,509,625 checks**. Three independent constructions on every labeled simple graph through order five; complete graphs through order thirty-two; one hundred deterministic randomized graph edge lists on at most eighty original vertices; simultaneous pair extensions of 130 arbitrary seeds, including every full interval seed through height thirty. Every output pair is checked against exact adjacency. The four small-prime residue searches are exhaustive.
2. `verify_two_degenerate.py`: **28,747 checks**. Every ordered at-most-two-earlier-parent graph through order five (counts 1,2,8,56,616); one hundred order-six/seven samples; per-step range and full adjacency checks; complete enumeration of all common-neighbor candidates for the scaled first two anchors of {1,2,3}, for square factors S^2 with 1<=S<=500, checking exclusion of the third. The latter covers 2,384 candidates and is a finite control, not the all-S Fermat proof.
3. `verify_prime_family.py`: **1,453,434 checks**. Complete graphs through order twenty-four for p in {3,5,11,13,19,29}, with full adjacency and height checks.

Total: **7,991,806 always-active checks**. These tests validate formulas and finite boundary cases; the universal conclusions use the proofs.

## Inheritance and attribution

The AP-anchor criterion in [ap_triangle_common_neighbors_20261019.md](../ap_triangle_common_neighbors_20261019.md) already characterizes the existence of a common neighbor for three AP anchors by the congruent-number class of the gap. Its square-gap specialization and the cited Fermat input are also present in [congruent_square_gap_classes_20261019.md](../congruent_square_gap_classes_20261019.md). This packet uses that inherited boundary; it does not claim to rediscover the AP-anchor classification, the previously known scaled K4 example, or Fermat's theorem.

The constructive contribution in this packet is the universal induced scaled-pair extension with all nonedges controlled, its finite two-degenerate closure, and the simultaneous quartic one-subdivision formulas. No external novelty claim is made. No outstanding mathematical blocker was found.
