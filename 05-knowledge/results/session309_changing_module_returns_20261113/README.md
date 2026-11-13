# Session 309: changing-module repairs and real clock sites

**Proved scoped results; Chan–Pak remains open.** This session begins at the determinant and boundary-refill equations of the [target-transition packet](../target_collatz_transitions_20261113/README.md), and separates a neutral refill, an unbounded physical move, and a complete target selector.

## A prime target really can refill its clock

[Physical refill](PHYSICAL_REFILL.md) gives a neutral rebase at the prime **83,405,299**. The original right module is strong and has no detachable trailing singleton. The replacement right module uses a different strong core and adds six actual trailing singleton sites, and the new left module is transitive. The seven exterior labels and all exterior-to-exterior arcs may be retained. Six successive one-arc transfers exhaust the new sites and end at the prime **950,033,299**. The intermediate transfer is H↦(3H−499)/2, with the new physical core furnishing the anchor 499. This is not ordinary Collatz, and we have not supplied another neutral return at the terminal prime.

A second [exact rebase at 7p](SEVEN_REFILL.md) changes two actual tail sites to five at **6,064,590,511 = 7·866,370,073**, with strong nontransitive bases on both sides. All five new transfers are literal. The new offset 1711 differs from the old module-count product, as a neutral refill must.

[Whole-module migration and recharge](MIGRATING_CLOCK.md) is a separate all-order rule: move an entire final ordinal factor across the separator and recharge with any actual module. This can continue indefinitely with nontransitive modules and gives an explicit family whose 7-divisible subsequence has cofactors running isometrically through every 7-adic residue level. This is a repeatable growth rule, not a repeatable **neutral** return, and makes no prime-cofactor infinitude claim.

## Why a return needs genuinely new modules and composite states

- [Fixed-bank obstructions](FIXED_BANKS.md): no fixed finite bank of strong components, even with unbounded ordinal repetition and interleaving their order, can produce infinitely many nontransitive determinant ±1 pairs. The same restriction permits only finitely many changed-intercept, target-neutral clock collisions at primes or 7p. These are applications of the S-unit theorem, with the degenerate cancellations treated separately. With fixed intercept, a dyadic valuation increase is impossible even with arbitrary module changes. A determinant ±1 output cannot keep its original product of module counts upon entering this clock, by a mod-12 mismatch.
- [Prime-lane Collatz obstruction](PRIME_COLLISION.md): for every fixed positive depth, infinitely many prime inputs have no prime or 7p descendant anywhere in their bounded-depth mixed-sign odd Collatz tree. This holds also from 7 times such primes. The argument uses a single dyadic cylinder to fix all valuations on the finite tree, then CRT and Dirichlet to force large independent factors at every descendant. Thus a bounded reverse-Collatz compiler cannot use only prime/7p states, even before one asks for tournaments. Certified composite carriers are a mathematical requirement of that particular strategy, not failed prime targets.

## Structural reading

[Hypergraph seams](HYPERGRAPH_SEAMS.md) follows the older Lonely Runner moving-clock work and the nearby historical hypergraph and tournament threads. Its same-sector refill analogy tracks a binding resource rather than aggregate counts. The tournament transfer uses the 3-uniform hypergraph of directed triangles: a disjoint packing of r triples forces H≥3^r, so the union of a maximal packing is a triangle transversal of size at most 3 floor(log₃ H). The proof isolates a bounded exceptional set without assuming that a finite set of strong component types suffices. See also [the synthesis map](THREE_SOURCES.md).

The next obligation is to grow a stock of **new** strong modules with proved port and tail behavior, then connect the target-directed determinant and clock charts through those modules and through composites. Neither a one-time exact refill, an indefinitely repeatable growth rule at selected values, nor a residue-class cover establishes that universal return.

Run `python verify_physical_refill.py`, `python verify_seven_refill.py`, and `python verify_migrating.py` beside the proof. All pass also with `python -O`. The verifiers check literal module responses and arc moves; the inherited universal frame formula supplies the large exact counts.
