# Independent audit of the six-name modulo-four response bank

**PASS.** Every one of the 15,859,712 canonical residue bytes was independently reconstructed, and all 32,768 complete named vectors are distinct. This is finite exact injectivity modulo four at six old names, not an all-order reconstruction theorem or a new maximum-order exact-integer record.

## Independent complete reconstruction

The producer's Hamilton-word algebra is correct: old-old edges contribute either an old-edge variable or its complement; old-port incidences contribute the appropriately named port variable or its complement; an adjacent y,x pair is impossible because the fixed arc is x->y. Expanding old-edge factors followed by the Boolean subset zeta transform evaluates the polynomial at each labelled tournament. No edge can be required both ways by a literal Hamilton word. The per-port degree is at most two.

The independent `audit.cpp` avoids that old-edge polynomial expansion and zeta transform. For each literal word it records the required old-edge orientations and explicitly enumerates every compatible assignment of the remaining old edges. Its port finite differences have mandatory incoming variables and optional outgoing variables with alternating signs. Contributions are accumulated directly at physical tournament codes. This covers all 35,280 permitted words and 98,795,520 compatible word/tournament pairs, with 344,309,760 direct residue updates. The resulting complete byte array equals the supplied bank byte for byte.

Both implementations retain all named coordinates, including zero, cubic and quartic coordinates. The independent basis is constructed from the Cartesian product of the empty set, six singletons and fifteen pairs for each port, then ordered by the combined physical mask. Its 22 squared coordinates agree with the declared 484-coordinate ordering. The independent old-edge index formula agrees with lexicographic physical pairs; bit one means the smaller endpoint beats the larger. All fixed x=6, y=7 incidences were checked against actual integer queries.

## Schema, exact equality and literal checks

The bank is exactly 484 by 32,768 bytes, coefficient-major. Each byte is checked to lie in 0,1,2,3. There are no multi-byte integers, so there is no native-endian or alignment convention to interpret. The audit rejects a short file, extra bytes, and noncanonical residues. Complete 484-byte vectors are inserted into an ordered exact container, not compared by checksums. All 32,768 insertions are unique.

The Boolean coefficient/value transform is triangular with diagonal one and has integer Möbius inverse, so it remains invertible modulo four. Distinct retained complete coefficient vectors therefore mean distinct actual response functions modulo four; this step does not divide by two or assume that the coefficient ring is a field.

A separate backward-ending integer subset recurrence evaluates every one of the 484 lower-closed interpolation queries on thirty source tournaments, using a control list different from the producer validator. Integer finite differences independently recover and compare all 14,520 coefficients. It also evaluates all 4,096 physical queries on each of codes 408 and 15462, retaining dense attachments that exercise every polynomial degree. Their fixed query (6,1) gives exactly 455 and 477.

The audit independently counts complete ordered-cover profiles of every old tournament by literal permutations and cut choices: every backwards adjacency forces a cut, and each forward adjacency is optional. All 196,608 coordinates are counted. Their q1 values agree modulo four with the all-zero two-port coordinate, and are positive odd integers. The universal q5=1,800 and q6=720 controls also hold. This does not infer universal profile log-concavity.

## Independent truncation and module-family diagnostic

`diagnostic.py` imports no producer code. It sorts exact byte vectors at each degree cutoff and independently reproduces every degree-class histogram and all 240 quadratic collision pairs. It constructs original modules from their physical exterior-incidence definition and exhaustively enumerates subsets of those modules, rather than using the producer's linear-algebra reduction. Exactly 120 pairs have their unique representation as fully covered laminar module XOR families; the other 120 lie outside the original-module span. Every saved family agrees. The first pair's complete modular contrast has exactly 48 cubic and 48 quartic nonzero coefficients and no lower-degree part. Degree at most three distinguishes this finite universe; the degree-four bank is still retained in full.

## Reproduction and pins

The independent C++ audit passes with both -O2 and -O0, with identical receipts. The independent Python diagnostic passes normally and with -O, also identically. The source C++ producer and validator were rebuilt with -O2 and -O0 under -Wall -Wextra -Werror. Their receipts match the saved records, and both produced banks match the frozen bank exactly. Both source Python checkers were replayed normally and with -O and reproduce their saved records. These are local exact executions, not a GitHub CI claim.

`source_pins.json` pins the reviewed source files, bank, and independent audit code. Restore the sibling bank with its supplied lossless reassembler before running:

    g++ -std=c++17 -O2 -Wall -Wextra -Werror audit.cpp -o audit
    ./audit ../covered_laminar/mod4_six_vectors.bin
    python diagnostic.py

No mathematical claim in this review uses a hash as an equality oracle. Integrity hashes accompany the complete data and exact comparisons.
