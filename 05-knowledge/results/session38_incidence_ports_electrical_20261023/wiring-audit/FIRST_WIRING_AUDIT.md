# Session38 independent audit

23 October 2026. **PASS / FINITE-EXACT** for the complete order-eight class,
port, colored multigraph and dyadic computations listed here. **PROOF REVIEW
PASS** for the scoped incidence and electrical arguments. No Git files were
modified by this audit. The inherited input is pinned to repository commit
`2130e7633805521d4e24355c37596bd2cd79f044`.

## 1. Independence and exhaustive coverage

This is not a sample and does not use the producer's Hamiltonian-path map
as a classification oracle. The audit independently:

1. Decodes the 6,880 inherited big-endian upper-triangle representatives.
   It minimizes over every ascending-score-block permutation, checks each
   automorphism count, checks distinct canonical words, and verifies the
   orbit-size sum is exactly 268,435,456 = 2^28. There are 403,193 checked
   score-block permutations. Independent subset DP verifies all 6,880
   Hamiltonian-path counts. Converse labels are separately canonicalized.
2. Computes exact E signatures directly from cyclic triples and proves them
   distinct on that independently certified complete class universe.
3. Builds each of the 2,097,152 fixed-path masks directly and classifies it by
   E. This differs from the producer's enumeration and normalization of
   Hamiltonian paths of representatives. Complete class coverage plus
   invariant injectivity justifies every classification, not a hash guess.
   C++ hash-table equality compares the complete exact signatures.
4. Rebuilds every colored ordinary endpoint-port and merged edge entry,
   including loops. The independently produced complete u16 class map is
   byte-for-byte identical to the Hamiltonian-path producer:

   `636f0ae8f69455d4d39143a7c612260cfbb2a9736c91ca67fe2ac3ea65449be1`.

The C++ audit checks signed-root-curl reconstruction on every one of the
2^21 masks. It checks every endpoint transport and every reflection/converse
relation. Reflection is position reversal plus converse, not free-arc
complement. Blue means a fixed actual presentation under that reflection,
not merely a self-converse class label. All black twin multiplicities are
checked even, for directed and undirected supports, including loops.

## 2. Exact wiring results and entrywise comparison

- 6,880 ordinary classes; 176 self-converse classes; 3,528 merged nodes.
- 4,096 blue presentations; 2,048 blue lines.
- 1,046,528 literal black lines; 523,264 paired-black lines.
- 956,339 ordinary endpoint/color/target entries, all compared exactly.
- 98,580 endpoint-capacity entries, all compared exactly.
- 460,699 directed merged endpoint pairs, or 460,816 color-separated cells.
- 438,216 undirected endpoint pairs, or 438,374 color-separated cells.
- 618 literal black loops on 232 supports, and four blue loops.

Every class-field and colored merged multiplicity comparison passes.
Independent recursive graph traversal/bridge computation agrees with the
producer's iterative implementation on all five inventories: literal full,
blue, literal black, paired black, and fundamental blue-plus-paired-black.
The full support is connected. Black support has one 3,525-node component
and three isolated nodes; blue support has one 176-node component and 3,352
isolated nodes. The full graph has eight support bridges but only three
literal multigraph bridges. All degrees are checked against actual port
fibers; loops contribute twice to degree. These are finite n8 results,
not all-order connectivity or bridge classifications.

## Reproduction scope of this first checkpoint

This checkpoint freezes only the independent full-wiring audit above. Its exact
E-injectivity check is part of the C++ replay. Later incidence minimality,
all-order proof, and complete electrical audits are separate additions. Read
FIRST_README.md for commands and FIRST_FREEZE_RECEIPT.json for explicit pins.
