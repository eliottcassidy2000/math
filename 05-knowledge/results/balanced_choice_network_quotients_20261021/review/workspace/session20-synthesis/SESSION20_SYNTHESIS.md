# Session 20: what the 4–2–1 connection preserves

Status: synthesis of four proof packets, each separately checked. This is progress toward tournament-metagraph classification, not a complete classification or a claim that all occurrences of 4,2,1 are the same structure.

## Incoming precursor and what this session adds

While this session was running, mainline added `metagraph_fibre_network_20261020k/README.md` (commit118eec39549). It already proves the restricted idoneal classification, marked regular-action fibre identifications, and one-step choice electrical, held-Markov, and entropy formulas. Those are inherited overlap, not newly claimed discoveries here. The present additions are the canonical all-dyadic conductor maps and blue-character marking with a nonliftability boundary, the literal all-depth two-seam defect theorem, exact choice critical groups and Kron comparison, and the actual reflection-paired metagraph parity/tree invariants. The prior packet's missing representative-port compatibility remains an open obligation.

## 1. The strongest positive connection

The balanced family's literal central-pair deletion has two commuting insertion toggles. Its fibres are Klein-four orbits with stabilizers 0, the left axis, the right axis, or the whole group. This explains 4,2,1 and, importantly, distinguishes two different types of size-two fibre.

In the Eisenstein conductor tower, the class groups at conductors 2,4,8 have sizes 1,2,4. These correspond to the idoneal integers 3,12,48. The inherited classification of this family does not depend on asserting a complete global list of 65 idoneal numbers.

There is an actual marked local map, not merely an analogy. Above a finite blue parent, label its four lifts by insertion bits (e,h). In the conductor-eight class group write A=w, B=j, C=wj. The map

    (e,h) -> B^e C^h

sends the two path toggles to B and C. Actual ideal extension to conductor four detects e XOR h. Its kernel {1,A} consists exactly of the two blue lifts. The genus characters recover the individual bits.

This requires a finite blue parent, an origin in its fibre, and a shore marking. For a black parent all four children remain black. On the exceptional size-two strata one must retain which shore is inactive. Neither individual-shore stabilizer is the diagonal conductor kernel.

## 2. A proved stopping point, rather than an indefinitely extended analogy

For a>=2 the dyadic class group is C_(2^(a−2)) × C2, with explicit residue generators and canonical extension maps. At conductor16 it is C4×C2, not an elementary Boolean group. Every lift of A has order four, so the map to conductor eight has no group section.

There is an even sharper symmetry obstruction. Path shore exchange swaps B,C and fixes A. It is a valid automorphism of the three-level 8→4→2 group tower, but cannot lift through16→8: the image of the conductor-sixteen two-torsion is {1,B}, which the swap does not preserve. Complex conjugation is inversion and is trivial at conductor eight; it is not this shore exchange.

The exact local bridge is useful precisely because its boundary is now proved.

## 3. Information lost by path deletion is concentrated at two seams

For any valid deleted Hamiltonian path p, deleting after reversing the original path arcs differs from reversing p by precisely the shortcut edges spanning internal deleted blocks. In the balanced family, arbitrarily many consecutive centre-flanking deletions still produce at most two such edges.

This is a constant-size incidence interface: retain the two actual endpoint pairs and whether each seam is occupied. It reconstructs the deleted target exactly. It does not make the undecorated deletion map a metagraph homomorphism, nor does it assert that the reconstructed target is reached by the same parent path.

An explicit pair of lifts in one fibre produces different deleted isomorphism types. Conversely the displayed target that is not a labelled one-step neighbour can be an unlabelled neighbour by another path. Both controls are necessary to state the obstruction correctly.

The exact defect counts also distinguish two growth regimes. Write q for the number of parent choices, D=2^d for the lifts of each finite coordinate, and Q=D(q−1)+1. For a uniformly chosen child path, the correction-size probabilities (0,1,2) tend, at fixed q and d→infinity, to

    (1, 2(q−2), (q−2)^2)/(q−1)^2.

At fixed d and q→infinity they instead tend to

    (0, 2/D^2, 1−2/D^2).

Proof: divide the all-depth packet's exact N0,N1,N2 formulas by Q² and take the leading coefficients in D² or q² respectively. Thus a bounded defect is not a rare defect. These are statements about the labelled reversal/deletion square, not adjacency after quotienting by isomorphism.

## 4. Three different meanings of an electrical quotient

The literal merged metagraph has network G=B+2Q: each black reflection pair contributes two parallel physical edges. On each connected component, its colour tree polynomial has all positive coefficients between powers c_B−1 and v−c_Q, and all nonzero roots are negative. Its value at2 is the physical spanning-tree count.

Its exact critical-group two-rank is c_B−1 plus the blue bicycle dimension. Thus non-self-converse vertices carry parity obstructions to integer chip firing. These are stronger than mass and degree counts, but still do not determine wiring. An explicit degree- and colour-preserving abstract switch changes the tree count from41096 to39232; realizability of that switch as tournament ownership is not claimed.

For the balanced **choice graph**, separately, contracting choice fibres of weights s,...,s,1 gives conductances w_x w_y. If m is the number of finite parent choices and Q=sm+1, the resulting critical group is

    C_s × (C_(sQ))^(m−1).

For s=2^d with d>=1 its two-primary part is (C_(2^d))^m. This is an exact full formula for that auxiliary network, not for the tournament metagraph. Its smallest V4 coincidence with the arithmetic group fails at the next rank or depth.

Shorting these fibres, eliminating interior vertices by Kron reduction, and lumping the complete-graph random walk are distinct operations. The first preserves constant-fibre energy; the second preserves terminal response; the third preserves the specified transition law only when holding probabilities are retained. None automatically preserves the integer firing lattice or the path-flip incidence.

## 5. A useful classification discipline

A proposed simplification should specify four data:

1. Carrier: literal path, tiling, ordinary class, converse-merged owner, choice state, or ideal class.
2. Operation: insertion, deletion, shore reflection, free-arc flip, shorting, Kron elimination, or ideal extension.
3. Preserved observables: blue predicate, incidence, energy, transition law, or integral parity.
4. Required retained data: inactive shore, parent coincidence, two seam endpoints, loop mass, harmonic extension, or marked ideal generator.

The actual metagraph's order-five resistance example makes the need concrete: terminal resistance is6092/5137, with odd denominator, while an internal voltage is11327/10274, with even denominator. Kron reduction is correct electrically yet hides the full-vector parity witness unless harmonic extension is retained.

## 6. Next proof obligations

- Extend the two-seam description beyond this balanced family using an explicit structural hypothesis, rather than arbitrary vertex deletion.
- Determine which local corrections remain legal marked Hamiltonian moves and how their owner labels transport.
- Classify the all-order blue bicycle space and the remaining two-primary critical-group data; finite values alone do not justify a recurrence.
- Test realizability criteria on full labelled incidence, not only counts, parity, or the colour tree polynomial.
- Keep black-core connectivity, all bridge classification, and unrestricted wiring classification open until proofs address those carriers directly.

The arithmetic and electrical calculations supply additional constraints and counterexamples. They do not substitute for these outstanding incidence obligations.
