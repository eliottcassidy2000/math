HAMILTON-PATH RELOCATION CONNECTIVITY
2026-10-20

PROVED. For every finite nonempty tournament, its labelled directed Hamiltonian words form a connected graph under legal single-vertex relocation.

The proof combines front promotion inside the initial strong component with vertex-insertion fibres and edge lifting. It explicitly proves why deleting the moved vertex from two distinct adjacent words leaves a directed path. No assertion is made that arbitrary vertex deletion preserves a path.

Consequences at every order:
- Auxiliary relocation components in the fixed-path coordinate cube are exactly the ordinary HP(T)/Aut(T) fibres.
- Adding reflection gives exactly the converse-merged fibres.
- The previously sufficient local relocation/reflection forest criterion is also necessary, with the original coordinate, complement, colour, and exact node-count conditions retained.

Read THEOREM.txt for the complete proof and quotient conventions, LOCAL_FOREST_ADDENDUM.txt for the exact realizability consequence, AUDIT.txt for proof review and computational scope, and PROVENANCE.txt for sources and boundaries.

EXACT_CHECK_RESULTS.txt records an exhaustive independent labelled check through n=6. It covers 33,867 labelled tournaments and 745,167 labelled Hamilton words across n=1,...,6, and verifies every relevant front-promotion orbit and 4,216,560 edge-lift configurations.

Reproduce from this directory with a C++17 compiler:
 g++ -std=c++17 -O2 -Wall -Wextra -pedantic check_labelled_connectivity.cpp -o check_labelled_connectivity
 ./check_labelled_connectivity 6

Do not disable assertions. The source rejects compilation with NDEBUG. No external packages, tournament atlas, network access, or precomputed isomorphism classes are used. The program intentionally accepts only exhaustive limits 1 through 6.

These auxiliary move edges are not complement-line metagraph edges. The triangle has three connected labelled words but just one Aut(T)-quotient vertex, so nontrivial labelled moves may collapse to quotient loops. No positive nonloop degree is claimed for every quotient fibre.

The result is a connectivity theorem, not a Gray-code, diameter, polynomial-runtime, or certificate-discovery theorem. The finite checks are supporting evidence, not the all-order proof. No literature novelty claim is made. Earlier result packets remain unchanged; this is an additive result.
