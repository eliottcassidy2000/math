# A signed move bank across a fixed tournament seam

Research session, 2026-11-16. Exact identities and literal constructions below do not settle Chan–Pak. The initial selection really was random: the seed and 24-entry sample are recorded in [the archive trail](ARCHIVE.md), supplemented by other independent seeded passes.

## Three consecutive counts at one growing cut

In a transitive tape with one head, insert an unmarked spacer of length t in one fixed cut. The count is **P·2^t+C**, where **P is a product of two selectable factors** and C is an exact bit-reversal/interaction remainder. Three choices move both factors while keeping P=812130 and making C=−391,−389,−387. The resulting strong tournaments give **1,623,869 → 1,623,871 → 1,623,873** at t=1, and the same two successive +2 repairs at every spacer length. Only head contacts change. This is stronger than matching one residue or producing two independently amplified moves; it gives two distinct transitions on the same growing cut.

[Product identity and exact divisor decoder](PRODUCT.md) extends the factorization to two arbitrary marked cores. It also proves why a common factor in the ordinary core counts rules out a unit on a fixed product once the spacer is nonempty. The earlier adjacent 13/13 repair used a term absent at every positive spacer length. A second construction uses two fixed **nontransitive** cores, with 55 and 31 paths. They really trade factors (1107,880)→(1188,820) without changing the high product; the low terms cancel as −9238+9240. It gives +2 across an arbitrarily long inner spacer with just four contact flips. [Fixed strong factor trade](INNER_SWAP.md).

## The rest of the bank

- **Fixed nontransitive seam, distinct signed repairs:** two fixed converse cores give literal differences **2, 6 and 8** at every symmetric guard height, while a pivot/separator variant gives **2, 8 and 10**. These are alternate edges, not independent commuting operations. [Mirror bank and interaction rank](seam_moving.md).
- **Genuine carry with minimal edit support:** when both guards are at least three, a nonzero carry on the stated transitive tape cannot be paid by one or two contact flips. An explicit family achieves **+2 with three flips**, while its actual high displacement grows like 3·2^i. There is also a finite **seven-count** window on one fixed transitive core with actual nonzero carry on all six steps; it is not known to persist under padding. [Carry theorem and fixtures](seam_bezout.md).
- **Another adjacent two-strong-core repair:** physical alternatives on the 13/13 seam and on a 9/13 seam let large contacts move. A fixed nine-vertex core supplies a strict three-odd lane for independently arbitrary exterior guard lengths. [Two-core checks and small-port obstruction](seam_cube.md).
- **A complete but bounded symmetric sector:** the exact trailing-one carry tariff gives arbitrarily many different +2 edges on one symmetric transitive tape as its width grows, but the longest consecutive window stays at three odds for even width and two for odd width. [Dyadic carry graph](DYADIC.md).
- **Different signed directions really are necessary:** at a prime p≡3 mod 4 dividing h, a repeated modular drift word must have length at least p; four successive steps cannot use only two directions. For a full modular cycle and p≥7, at least four directions are necessary. The local norm cycle still has to lift to physical masks and pay the absolute carry. [Direction-clock proof](DRIFTS.md).

Each explicit witness has an independent arithmetic or literal-path replay. The fixed-cut search is bounded; none of the finite negative scans is an impossibility theorem.

## Next exact target

Use the product decoder in reverse. Choose a changing **common product** with several admissible factor pairs, then force their digital remainders to hit the specified consecutive odd values. For adjacent nontransitive seams an alternative is to lift a curved drift through the actual two-contact geometry, without asking the same direction to repeat or presuming a modular cycle pays its absolute carry.
