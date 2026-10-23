# Conference Gram/theta collapse: independent audit

PASS after the producer applied the explicit n>2 qualifier to the universal-switch theorem heading, assumptions, conclusion, and final summary. The corrected final wording was inspected. The n=2 class has only two switches, each containing both a source and a sink; it does not have four distinct universal switches.

## General algebra

For positive even n and C^2=-(n-1)I, diagonal-sign switching gives C_d^2=-(n-1)I entrywise. Thus Q=nI exactly, not only up to orthogonal or integral equivalence. Every stated finite heat, spectral zeta, theta and Epstein equality follows directly. Determinant normalization sends Q to I. The rational resolvent expression is proved away from its poles and the adjugate identity then extends polynomially, so no singular-point division is hidden.

Set x=2lambda+1 and a=n-1. The determinant lemma gives det(lambda I-A)=2^{-n}(x^2+a)^(n/2)(1-nx/(x^2+a)). The numerator x^2+a-nx=(x-1)(x-(n-1)). Replacing x proves exactly lambda(lambda-(n-2)/2)(lambda^2+lambda+n/4)^(n/2-1). Thus the assertion of adjacency cospectrality is valid, independently of the finite data.

For each prescribed vertex and uniform sign there is one sign-switch class modulo global sign. Two sources or two sinks are impossible in any tournament. A source and sink have row inner product -(n-2), incompatible with orthogonality when n>2. This proves distinctness and exactly2n such switches for n>2. Deleting the universal vertex gives B1=0 and B^2=J-(n-1)I by multiplying blocks, the skew-matrix characterization of a doubly regular tournament. The theorem is all-order, not an extrapolation from order8.

## Independent order8 certificate

The alternate script checks Q=8I and the exact adjacency characteristic polynomial on EACH of the128 switches. It independently classifies them using recursive partial-isomorphism extension with edge-consistency pruning. There are four score sequences, and every switch with each score sequence is explicitly found isomorphic to its chosen representative. Because the four score sequences differ, the four classes are distinct. Counts are8,56,56,8.

Hamiltonian paths on representatives were counted by direct enumeration of all8! vertex permutations, independently of the producer's subset DP. Directed cycles were counted by anchored subset DP, independently of the producer's cyclic permutation enumeration. Pair packings were then enumerated by disjoint vertex masks. Results:
- universal types: H189, cycles/packings (14,42,24,7,0), Aut21;
- nonuniversal types: H333, cycles/packings (14,42,48,19,12), Aut3.
Thus the inherited odd-cycle formula yields the exact144 difference as48 from7-cycles,48 from disjoint triangle pairs,48 from triangle-pentagon pairs. These counts do not follow merely from equality of adjacency spectrum.

Source masks are exactly {0,23,46,57,75,92,101,114}. Complemented sink masks are exactly {0,29,39,58,78,83,105,116}. Both are independently verified as size8 XOR-linear codes with seven weight4 nonzero vectors. The complement supports have size3 and cover each unordered pair once. They are different labelled codes; the fixed-gauge caution is correct.

## Inheritance review

Read the cited local sources. THM-472 explicitly records four skew-conference extremal classes at order8. THM-450 section6 and THM-454 record the189/333 doubling split. THM-505's order-eight extension supplies the disjoint (3,3) and (3,5) packing carriers. Therefore the note correctly avoids claiming these inherited facts or the OCF as new. Its stated emphasis is the explicit full entrywise Gram/theta collapse and accompanying labelled certificate.

The test universe is exactly the128 sign switches of the supplied bordered Paley matrix, not the entire labelled order8 universe. No reconstruction, RH, or word-algebra completeness claim is made. The supplied output was reproduced byte-for-byte. Alternate tests are check_new_additions.py/.out in this audit directory.

## Final artifact SHA256 pins

- session35-zeta/CONFERENCE_COLLAPSE.md: 11530a605669f749a6c4e2eef68409a203180b3678f4835ba8000d46ffe51953
- session35-zeta/CONFERENCE_COLLAPSE.py: 62a87ff7672d3eb50d01f649d274a91ae19dd0e36ebb287f1553322a93d3f08d
- session35-zeta/CONFERENCE_COLLAPSE.out: 7f0f6609a58c96099bbb22fd4aa0ba06998e4bd1421efcd9ddfba4e730e34c14

No Git operations or repository edits by this auditor. The initial optional NetworkX-based alternate test could not run because NetworkX was absent; it was replaced with a self-contained exact backtracking implementation, which passed. No test failure or dependency remains unresolved.
