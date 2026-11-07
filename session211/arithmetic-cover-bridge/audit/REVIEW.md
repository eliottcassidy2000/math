# Independent arithmetic-cover bridge review

**Verdict: PASS.** Final proof SHA-256: `d87ac9b71303f3114b9ddd9c15b27cbbf59d24984d4cd9a1583bc6c12c6601f8`.

## Reviewed scope

1. The word-cut formula and arc-complement reciprocity hold for every finite loopless digraph, allowing missing and reciprocal arcs. The empty core is treated separately where necessary. Tournament parity is a specialization, not a recognition criterion for arbitrary digraphs.
2. The square-zero subset-series proof correctly extends the gate operator to every such digraph, including the empty one. The path/cycle bijection requires a nonempty core. For the empty core, `L(1)=q^2` remains valid but `HC(Gamma(empty))=0` differs from `H(empty)=1`.
3. Missing-pair and reciprocal-pair counts are conserved by every gate iteration. The gate does not complete or orient the old comparison alphabet. This obstruction is independent of counts and strongness.
4. The literal increasing arithmetic graph on `1,...,N` with steps `+2` and doubling has no Hamiltonian path for `N>=3`, but exactly two unordered two-path covers. The largest-odd endpoint argument completes the proof that those covers are the only ones. Consequently its strong gate has exactly four Hamiltonian paths and no Hamiltonian cycle. Symmetrization and transitive completion have the different stated path/cycle behavior.
5. The order-three hostile with arcs `0->1,1->0,1->2` has the same full cover polynomial `q^3` as a transitive tournament while retaining one missing and one reciprocal pair. All iterated cover/path observations coincide; positive-depth cycle observations coincide too. The order-two comparison polynomials prove minimality.
6. The final permutohedron warning is explicitly local to one-dimensional faces. Missing or reciprocal comparisons cannot be represented by the same local sink predicate as a single ordinary edge orientation.

The statement distinguishing literal `BA=A^2B` from commuting row-scaled odd-core operations is correctly typed. The square-sum values and tournament omissions are inherited from the separately audited companion packet. No general spectrum completeness, arithmetic Hamiltonicity threshold, or novelty claim is made.

## Reproduction and independence

The final verifier was copied and executed under normal Python and `python3 -O -B`; both outputs match frozen `RESULT.json` byte for byte. This exhausts all 4,166 labeled loopless digraphs of orders zero through four, checks arithmetic graphs through order ten, and checks four nested responses of the explicit order-three hostile. The all-order conclusions rest on the reviewed proofs, not the finite census.

An independent pure-permutation checker, `independent_common_gate.py`, uses no producer functions. It verifies complement reciprocity and the gate polynomial on all 69 labeled loopless digraphs of orders one through three, including reciprocal/missing arcs, at sufficiently many positive and negative integer arguments to determine each tested polynomial. Rooted explicit cycle enumeration checks the nonempty-core bijection. Separate controls record the correct empty-core polynomial and excluded cycle-bijection boundary. Its normal and optimized outputs agree exactly.

All four entries in `SOURCES.json` were independently verified against the immutable reference commit, using both Git blob IDs and SHA-256 values.

## Artifact hashes

- `verify.py`: `3b1f1c4423f1833c0b2f180162512476f1160e5aa9fd6200794a34303d174567`
- `RESULT.json`: `5d7d0279e5edeae5486df560b298e309292231d666a1c8709fac608b88e572e1`
- `SOURCES.json`: `f3a29de6cfcdeb334c03302a5f41b325d433e9c5b691576c8339c651199f4dd6`

The independent checker and its normal/optimized outputs are approved for inclusion beside this review.
