# Exact self-converse completion by transitive cloning

This reviewed Session197 packet extends the one-root clone-length theorem to every positive weight vector on a fixed canonical transitive-module quotient.

- [All-order proof](PROOF.md)
- [Independent audit](audit/REPORT.md)
- [Standalone source verifier](verify.cpp)
- [Standalone independent checker](audit/check.cpp)

Self-converse clone lengths satisfy a finite union of explicit block-weight equality systems. The fewest added vertices needed for self-converseness is exactly the minimum, over reversing involutions of the quotient, of the sum of absolute weight differences on exchanged block pairs. Immutable blocks add explicit feasibility conditions. This lower bound covers arbitrary repeated cloning and changing roots, and the proof supplies an optimal physical construction whenever one exists.

A quotient that is not self-converse can never be repaired by this cloning language. When the seed is strong, an attainable completion remains strong. The operation does not preserve Hamilton-path counts.

The packet also gives exact ordinary and self-converse class cycle indices, labelled counts and asymptotic symmetry frequencies for all blow-ups of one fixed quotient. For the Paley seven-vertex quotient, the self-converse fraction is asymptotic to105/n^3. The singleton quotient, even-order parity and nontrivial anti-involution centralizers are treated explicitly.

These are fixed-language structural and enumeration results, not unrestricted spectrum coverage or prime infinitude. The separate C3 completion-cost distribution is not part of this checkpoint.

## Reproduce the controls

Both programs are standalone C++17 with no input-file dependencies.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug
    g++ -O3 -std=c++17 audit/check.cpp -o independent
    ./independent
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 audit/check.cpp -o independent-debug
    ./independent-debug

Source outputs match receipt.txt and receipt-debug.txt. Independent outputs match audit/normal.json and audit/debug.json. Source controls check3938 actual blow-ups and368286 finite completion problems. The independent checker uses reverse breadth-first search for the completion minima, separate physical graph-isomorphism backtracking, full subset module tests, and disjoint-set orbit enumeration. Finite checks supplement the all-order proofs.

## Unchanged mathematical dependencies

The source and review retain their frozen session-local names. Their release aliases are dependencies/CONVERSE_CLONE_LENGTHS.md for the canonical-block theorem and dependencies/CONVERSE_MASS.md for the odd-order group and anti-involution results. The other six included dependency proofs are the previously reviewed arithmetic prerequisites of the canonical-block packet: UNIT_CLOCK_REDUCTION.md, MARKED_BINARY.md, MINIMUM_TERNARY.md, NEGATIVE_BASE_RIGIDITY.md, SYMMETRIC_SQUARE.md and EVEN_CONVERSE.md.

Ordinary and twisted Burnside counting is proved directly here and is standard group-action machinery; no literature-priority claim is made. RELEASE_FILES.sha256 binds every other included file.
