# Reconstructive support, first-birth multiplicities, and primitive seams

**Anchor: the arcs you can strip away.** [Support rigidity](SUPPORT-RIGIDITY.md) removes every arc unused by all Hamiltonian paths. The support still determines the tournament: on each missing pair the orientation is the direction with odd relative-order count. The opposite order occurs a multiple of four times. This gives a recursive, target-pruned grammar of path words over partial orientations rather than a search over arbitrary missing edges. [The path-conflict transform](PATH-CONFLICT.md) is a second grammar, using compatible path atoms and Burnside to give both the scalar and the marked one-arc class polynomials.

**Niche: why minimum birth need not be unique.** [Fixed-count size series](FIXED-COUNT-SERIES.md) extends the inherited rational SCC factorization with exact transient thresholds, converse-parity terms, first-witness recursions, and uniform asymptotics on a growing finite-prime-bank fiber. Minimum size optimizes total strong order; eventual growth optimizes the number of strong factors instead. [The complete first-birth census through order nine](BIRTH-MULTIPLICITIES.md) explains actual collisions: the two smallest H=9 classes are both self-converse with identical scores, but have different deletion profiles and different stripped supports. The H=1467 and H=1503 minimum fibers each have 257 classes, one self-converse and 128 converse pairs. The [full 1520-value table](first_birth_fibers.tsv) and a large-order class-distribution theorem accompany the examples.

**Wildcard: the joint perfect center.** [Primitive seam grammar](PRIMITIVE-SEAMS.md) proves that after removing common transitive padding, every literal +2 pair with lower count at least 5 has two strong cores. Its exact size series is the finite primitive polynomial divided by \((1-z)^2\); the flipped arc must be used by Hamiltonian paths on both sides. [Exact perfect centers](PERFECT-BIRTH-SEAMS.md) gives minimum orders 6, 8 and 10 at exponents 3, 5 and 7, and a strong 24-vertex witness at 13 without claiming minimality. At exponent 5, 22 and 33 endpoint classes first occur at order eight, but only 18 simultaneous pair classes connect them by one arc there. The earlier [exponent-31 product cube](../chanpak_e31_physical_20261122/PRODUCT-CUBE.md) provides a complementary multiplicative mechanism: component ordering distinguishes classes even when only the number of reflected choices survives in the scalar count. Its exterior multiplication is exactly the scale ruled out for a genuine +2 edge.

Exponent 2 cannot straddle because 7 is absent. An unbounded literal straddling construction and Chan–Pak still remain **OPEN**.

## Replay

    python3 check_core_grammar.py
    python3 check_fixed_count_series.py
    g++ -O3 -std=c++17 birth_classes.cpp -o /tmp/birth-classes
    /tmp/birth-classes 9 > /tmp/birth-classes.tsv
    python3 check_birth_multiplicities.py /tmp/birth-classes.tsv
    python3 check_perfect_birth_seams.py /tmp/birth-classes.tsv

Compile `check_support.cpp` and run it from the repository root to use the stored class decks. The census tests and finite masks are independent of the all-size proofs.
