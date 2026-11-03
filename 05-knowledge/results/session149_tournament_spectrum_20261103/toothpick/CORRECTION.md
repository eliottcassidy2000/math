# A scoped correction to the toothpick–tournament correspondence

This note preserves the earlier record in `chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213.md`. Its common involution identity, |X| = |Fix(g)| + 2·(# two-element orbits), is correct. But that identity alone does not give an orbit-for-orbit correspondence between two actions. The later reflection `the-self-similar-object-is-the-operation-congruence-not-the-count-codex-S16.md`, §1, already cautions that early toothpick and tournament census coincidences are not a bijection.

## Convention and symmetry

Start with one vertical unit toothpick centered at the origin. At each subsequent generation attach a perpendicular unit toothpick centered at each exposed endpoint of the preceding generation, using the standard simultaneous update. Occupied endpoints do not support a new toothpick. Equivalently, scale the lattice by two, so toothpicks have length two and endpoints have integer coordinates.

Horizontal and vertical reflection, and hence half-turn rotation, preserve the seed and the update. By induction every finite stage has these symmetries. The quarter-turn does not preserve the initial seed: it sends it to a horizontal toothpick centered at the origin. Nor can this missing segment appear later, because the origin is occupied by the interior of the initial toothpick and a new toothpick is centered at an exposed endpoint. Thus the standard seeded configuration has the four-element D2 symmetry group, not D4 symmetry. This conclusion concerns the geometric pattern including segment orientations, not just the set of intersection points.

The half-turn fixes a segment exactly when its midpoint is the origin. Thus it fixes exactly the seed at every stage. All other segments come in two-element orbits. In particular every cumulative stage has odd cardinality. For either coordinate reflection the number of fixed segments is also odd, by the same orbit-counting identity; the seed need not be the only reflection-fixed segment.

## A concrete obstruction

For every n ≥ 3, at least two distinct isomorphism classes of n-vertex tournaments are fixed by converse. One is the transitive tournament Tr_n. For odd n = 2k+3 another is the ordered sum Tr_k ⊕ C3 ⊕ Tr_k. For even n = 2k+4 another is Tr_k ⊕ Q4 ⊕ Tr_k, where Q4 has arcs

    0→1, 0→2, 1→2, 1→3, 2→3, 3→0.

The reversal permutation (0 3)(1 2) identifies Q4 with its converse, and C3 is self-converse. Reversing block order and identifying the two transitive blocks proves each displayed ordered sum self-converse. Each contains a directed triangle and is therefore distinct from Tr_n.

An equivariant injection must send distinct fixed points to distinct fixed points. Consequently, for any n ≥ 3 and any toothpick stage, there is **no equivariant injection**, let alone an equivariant bijection, from the set of all tournament classes with converse action to the set of toothpicks with half-turn action. This is an explicit obstruction to the suggested literal orbit-for-orbit dictionary for that chosen geometric involution. For reflection instead of half-turn, its odd fixed-point count already rules out an equivariant bijection at n = 3: both tournament classes are self-converse. It does not rule out all possible weaker or differently defined constructions.

The verified count identity and the signed-cycle-index result `THM-587` remain separate valid statements. Nothing here proves or disproves a general parity theorem for tournament counts, a connection between their spectra and toothpick dynamics, or any nonstandard encoding with multiple geometric objects per tournament. Those require additional definitions and proofs.
