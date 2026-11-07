# First-occurrence fibers through eight vertices

Status: FINITE-EXACT for complete orders 1 through 8. PROVED refutation of the proposed one-converse-orbit first-occurrence law.

Define b(H)=min{n>=1: some n-vertex tournament has H Hamiltonian paths}; B_H is the set of isomorphism classes on b(H) vertices with that count. Converse reversal preserves H and induces an involution on B_H. Its complete cardinality profile is (s(H),p(H)): s self-converse classes and p two-element converse pairs. Thus |B_H|=s+2p and the number of converse orbits is s+p. The proposed law is exactly s+p=1. The only general law supplied here is involutivity, not transitivity of its action.

## Smallest counterexample

H=9 is born at n=5 with (s,p)=(2,0). On vertices 0,...,4, two witnesses have outgoing-neighbor bitmasks

A=(14,28,24,16,1), B=(10,28,25,16,1).

Both have exactly nine Hamiltonian paths (permutation enumeration suffices). Both are self-converse: anti-isomorphism maps are respectively (4,3,2,1,0) and (0,4,3,2,1). They are nonisomorphic: the unique outdegree-2 vertex of A has outneighbor-degree multiset {1,1}, while that of B has {1,3}. This also proves they are not converses. The complete n<=4 spectrum is {1,3,5}; H=1,3,5 have one class at birth. Accordingly H=9 is the smallest realized count contradicting the law, and n=5 is its earliest failing order. No claim about unrealized H=7 is needed for this conclusion.

H=15 is also born at n=5 with (2,0). H=33 is born at n=6 with (1,2), i.e. five ordinary classes and three converse orbits.

## Exact census and verification

birth_fibers.csv records all 320 H values attained through order 8; their b(H) entries are globally valid, because every smaller order is exhaustively checked. For values not in this table no assertion of birth/nonexistence is made. Of the 320 values, 285 have more than one converse orbit at birth. This is not a statistical sample of all odd numbers.

Universe: ALL tournaments, with no strong/prime/module-free filter. H means actual Hamiltonian-path count, never divided by automorphism order. Positive controls: class counts 1,1,2,4,12,56,456,6880; self-converse counts 1,1,2,2,8,12,88,176; H=9 and H=33 are hostile controls for uniqueness.

The verifier independently enumerates all labelled tournaments at n<=5. At n=6,7 it reads existing representatives from census.json, and at n=8 from n8_classes.csv, then recomputes Hamiltonian paths by subset dynamic programming, canonical forms by all permutations within outdegree cells, automorphism orders, self-converseness, converse closure, and distinctness. Degree cells are ordered by increasing degree; any isomorphism respects them, so this is a genuine canonical form. Counting permutations giving the minimum gives the automorphism group size. For each n, distinct canonical classes satisfy sum n!/|Aut(T)| = 2^(n choose 2). This independently certifies completeness: they are disjoint labelled orbits whose total size equals the entire labelled universe. The Hamiltonian presentation identity sum H(T)/|Aut(T)|=2^((n-1)(n-2)/2) supplies an additional independent channel. Checks use explicit exceptions and remain active under python -O.

A separate C++ rebuild of all n=8 representatives from the 456 n=7 classes agrees byte-for-byte with the inherited n8_classes.csv. Every order-8 tournament has an order-7 deletion, and every orientation of the seven new arcs is generated, proving augmentation exhaustiveness.

## Reproduction

Inputs reside in ../session24-multiplicities (or pass that directory explicitly): census.json, n8_classes.csv, n7_representatives.txt. Their hashes are in receipt.json. No source file is edited.

python3 verify.py /path/to/session24-multiplicities

g++ -O2 -std=c++17 replay_n8.cpp -o /tmp/replay_n8

/tmp/replay_n8 /path/to/session24-multiplicities/n7_representatives.txt > /tmp/n8.csv

cmp /tmp/n8.csv /path/to/session24-multiplicities/n8_classes.csv

The replay_n8.cpp is the inherited session24 augment8.cpp copied without mathematical changes; verify.py is an independent Python verifier and birth-ledger extractor. The minimal H9 witness does not depend on the larger atlas.
