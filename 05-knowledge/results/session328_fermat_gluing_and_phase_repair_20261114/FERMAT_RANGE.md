# Three genuine Fermat blocks cover every odd target from 69 to 79,115

**Finite-exact theorem.** Every odd integer 69≤N≤79115 is the Hamiltonian-path count of a strong tournament obtained by one ordered overlap of a strong source A of order at most nine with one of three blocks R₁,R₂,R₄. The outer block's own count is respectively 3,5,17, the first three genuine Fermat numbers. The final tournament has order at most fourteen. This strengthens the proposed end of 60,000 and does not rely on arbitrary exponents j, products of counts, or a claim about all larger targets.

## The operation and exact witness formula

R_j is the near-transitive order-(j+2) tournament: put vertices in order and reverse only the extreme edge. It shares its first vertex v with a marked source vertex v; all arcs between A−v and R_j−v go from A to R_j. Both factors being strong implies the result is strong. The block has H=2^j+1 by splitting the interior vertices on the two sides of its reversed edge, plus the direct forward path.

For the source keep H_A and the exact relaxation L_A. Remove v and split the other vertices into an ordered before piece and after piece, each a directed path, with empty pieces allowed. Require only that the after piece be empty or start at an out-neighbor of v. The resulting number of choices is L_A. As shown by the [ordered seam formula](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md), the exact output is

    N=H_A+2^j L_A.

Thus using j=1,2,4 does not mean multiplying H_A by 3,5,17. It reads a marked endpoint state and adds a scaled response. That distinction is necessary for primes and is made explicit in [the historical comparison](FERMAT_GLUE_PROVENANCE.md).

## Self-contained finite certificate

Run

    g++ -O3 -std=c++17 replay_fermat_range.cpp -o replay_range
    ./replay_range

This has no external witness files or number library. It enumerates every labelled source through order six, then uses an explicit unsigned 64-bit xorshift stream (seed 325, shifts 13,7,17) to generate a deterministic set of 200,000 candidate masks at each of orders seven, eight and nine. This is a reproducible way to **supply positive witnesses**; it is not claimed to exhaust those orders. For each strong candidate the verifier computes all subset Hamilton-path/end tables directly, using the converse to compute starts. For each root, summing over the after subset gives L exactly. The verifier constructs a witness for each of the 39,524 prescribed odd integers and fails if any is not present. It also independently assembles full graphs and counts paths on the selected boundary, historical-prime and block-type controls without using the transfer formula for those checks. The finite cover result is literal because every accepted pair has a recovered graph and a legal seam.

The cached receipt may report a first unwitnessed target for a restricted choice of blocks. That is only relative to this witness generator, **not** a proof of universal failure for that restriction. The next paragraph gives a separate universal bound that does not depend on sampling.

## A genuine obstruction to removing the third block

For any root of a source of order at most nine, H_A≤L_A≤Q_(A−v)(2). The last is the optional two-slot count with no root guard. The exact maximum of Q(2) over all tournaments of order eight is 8848, proved by complete enumeration of the 6880 isomorphism classes in the supplied [small verifier](replay_fermat_bound.py). Smaller orders only reduce the maximum, since singleton ordinal padding doubles Q(2) per vertex. Consequently blocks restricted to j≤2 have N≤(1+4)·8848=44240, even permitting nonstrong sources. They cannot cover the proposed range. Allowing j=4 suffices on the entire certified interval.

## Scope at infinity

For a fixed finite bank of order-at-most-nine sources and roots, each rooted pair contributes O(log log X) values through X if j is restricted to genuine Fermat exponents (powers of two), and O(log X) even if all j are allowed. The strong finite coverage here cannot simply be extrapolated. An all-scale proof needs a growing source grammar or another type of exact cut. The prior non-Fermat selected-ear charts supply a different finite construction and do not depend on this certificate.

The inherited complete H-census for order nine has maximum 3357, so combining that independently with the Q8 maximum sharpens the cited obstruction to N≤3357+4·8848=38749 for j≤2. The 44240 argument above does not rely on that separate census.

To print a literal source, root and length for an individual target after checking the entire range, run `./replay_range 30013` (or pass several target integers). The frozen positive certificate with j≤4 produced digest 4952078363158291100; the three preferred exponent counts were 536,6280,32708. The separate conditional maxima above do not rely on this digest or on the witness generator's choice of preferred source.
