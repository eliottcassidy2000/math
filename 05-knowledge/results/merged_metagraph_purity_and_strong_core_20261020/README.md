# Pure-blue points, the strong core, and exact complement cuts

This packet concerns the converse-merged tournament metagraph whose literal lines pair fixed-Hamiltonian-path tilings by complementing every nonpath arc. It does not concern the single-arc-flip isomorphism graph.

## Main results

* The unique transitive blue pipe lands in a mixed node for every n>=4, with an explicit black normalization in that target fibre.
* Every nontrivial strong pure-blue tournament is odd, regular, and circulant along every directed Hamiltonian cycle. The proof uses all Hamiltonian-path reversal permutations, Camion's cycle theorem, and Lu's published crossing theorem.
* Every even pure-blue tournament is self-converse and rigid. At even n>=4 its merged node is a blue leaf attached to a mixed strong node.
* Every complement partner of a nonstrong tournament is strong. The nonstrong locus is an independent periphery, with no loops.
* Exact recurrences for strong and strong-blue tilings force the total and colour-resolved line counts across the strong/periphery cut and inside the core.
* An exact n5 two-edge surgery preserves all node colour-degrees and total counts but violates the blue triangle-current law. The regular pure-blue n5 node also refutes the historical unrestricted pendant claim.

## Reading order

1. FOUNDATIONS_AND_MIXED_LANDING.md: definitions, inherited constraints, hostile rewiring, correction lineage, and the transitive-pipe proof.
2. PURE_BLUE_STRONG_CYCLE_REDUCTION.md: cycle-cut reduction, even rigidity, core/periphery support, and Hamiltonian-cycle injection.
3. STRONG_CORE_COUNT_RECURSIONS.md: exact count-to-arrangement formulas and finite checks.
4. LU_CROSSING_SOURCE.md: primary-source statement and precise access scope.
5. AUDIT.md: proof-review and computation scope, including inherited claims not reaudited.

Source snapshot: repository commit 129154b04024f1246570592ac8df4db249712bc4. Earlier theorem references are mathematical dependencies, not claims of fresh proof. In particular the through62 rigid-atom consequence is explicitly conditional on that inherited theorem; this packet makes no independent all-order rigid-atom or finiteness classification claim.

## Reproduction

Python3, standard library only. From any working directory:

    python3 /path/to/this/packet/verify.py

The verifier runs all three portable scripts and compares their stdout byte-for-byte with the included outputs. They cover all64 tilings at n5 with exact isomorphism canonicalization; the explicit transitive landing family n3..100; and every fixed-path tiling n3..7 with independent reachability and coordinate-mirror tests. Recurrence-generated values through n12 are explicitly distinguished from direct enumeration.
