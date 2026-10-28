# Four vertices, exact compression, and faster tournament enumeration

Session 87, 2026-10-28. **Proved mechanisms, independently checked finite examples, and measured implementation improvements.** Classical Burnside counting, partial-spread theory and canonical construction paths are credited as inherited mathematics. This packet does not claim a new formula for A000568 or an unrestricted prime-class atlas.

## Main results

1. **A genuine K4 realization of the F64 six-to-four query gap.** The three binary planes of the nine-port alphabet can be assigned to three adjacent edge pairs of K4. Three vertex-boundary bits plus total edge parity separate the nine ports. Exactly eight edge-pair partitions work, each a root and a directed triangle. Opposite-edge pairings fail. The map merges different tournament classes on the full 64-orientation universe.
2. **A general compression theorem with a closure boundary.** One-active-block linear identification is precisely a partial-spread problem. Three two-bit blocks need six fixed local queries, five adaptive local queries, or four unrestricted linear queries. Five such blocks give fifteen states and still need four unrestricted queries. Allowing more simultaneous active blocks changes the required code distance; full additive closure restores the full information cost.
3. **Faster scalar counting against the real prior frontier.** Degree buckets, shared transfer metadata and a terminal unit-tail fold preserve the October degree-local response quotient. Five-run median CPU gains against its telemetry-light implementation were 1.54x, 1.50x and 1.60x for complete prefixes through 100, 150 and 200. This is an implementation improvement, not a smaller state quotient. A reusable kernel helps sufficiently large palette batches but loses on one explicitly reported smaller workload.
4. **Actual class generation without building the full metagraph.** An invariant root filter followed by an exact aligned deletion-key minimum selects one pointed orbit per child. Both selectors reproduce all 6,880 classes at order eight. The refined selector uses 693 competitor-key requests instead of 434,048 full-Q entries; these are work counts, not a speedup ratio. A separate same-cache benchmark measured 0.932 versus 8.933 seconds, with the important distinction that the faster path does not output Q itself.
5. **A single lost four-vertex type bit changes larger attachment symmetry.** The successful wedge scaffold forgets whether its root is a source or sink. In a three-block cyclic composition, changing that bit in one block gives the 81/27 automorphism pair and 176/512 binary attachment-orbit counts. The same labelled cyclic triples and all 4,096 induced Hamilton-path counts survive. Thus the local decoder cannot replace ordinary child types and transporters in the composition grammar.
6. **The historical comparison found two executable bugs.** A legacy odd-part generator gives A000568(5)=10 instead of12. A hypergraph enumerator's silent composition cap gives92 instead of93 subset orbits and even depends on the order of identical cycle-block data. Narrow, tested source repairs and a repository-relative regression checker accompany the research. General legacy CRT and large-input C arithmetic remain outside the repair claim.

## Read the parts

- [K4 bridge, earlier four-vertex routes, Hamming/Fano codes and arithmetic limits](K4_BRIDGE.md).
- [Block compression, partial spreads, adaptive queries and closure](BLOCK_COMPRESSION.md).
- [Faster exact scalar and colored counts, proofs and benchmark tables](FAST_COUNTS.md).
- [Exact class-generation selector, 81/27 hostile and composition grammar](EXACT_GENERATION.md).
- [March-to-October history, related sequences, marking identities and legacy repairs](HISTORY_AND_TRANSFERS.md).

## Three different exact quotients

The useful unification is a rule about what a quotient is asked to preserve, not an identification of all these objects.

**Identification on a promised source.** A linear map L is injective on a set S exactly when its kernel avoids every nonzero difference of two points in S. The K4 bicycle kernel passes for a suitably marked nine-port alphabet and fails for its additive closure. Enlarge the permitted source or operations and recheck the difference set.

**Aggregation for future counting.** Two partial cycle histories may be merged when every still-possible continuation has the same response to them. Their weights are summed, not individually reconstructed. The future-GCD profile provides that exact equivalence for the finite-horizon Burnside dynamic program. Degree, stage, arithmetic scale and observer payload must remain. Colored or marked counts can reuse this machinery; full prime cycle indices require more than a scalar payload.

**Forgetting a mark while generating classes.** Pointed objects are quotiented by actual root transport. An invariant filter chooses which roots to consider, but only the complete aligned deletion key resolves the remaining root orbits. Matching degree signatures, cyclic triples or Hamilton counts does not suffice. The 81/27 example is the concrete obstruction.

This separates three legitimate savings: a code tailored to a restricted alphabet; a sufficient response state for a weighted sum; and an exact orbit selector avoiding redundant generation. None creates missing information. Each avoids paying for information that its particular output does not require.

## What transfers to other problems

- **Square sums and graceful labels:** incidence signs and binary boundary structure can share a K4 scaffold, but actual integer labels, positivity, distinctness, difference resources, endpoint names and connectivity are extra data. The notes give explicit failures of discarding them. No square-sum Hamiltonicity or graceful-tree conjecture is solved.
- **Related counting series:** the colored count C_n(c) gives rooted tournaments via C_(n-1)(2), and ordered k marked vertices via 2^binom(k,2) C_(n-k)(2^k). These compute attachment workload and, using the inherited root-switch theorem, even-order binary homology dimensions without constructing Q. Ordinary graphs use all partition types; the tournament odd-cycle filter cannot be copied unchanged.
- **Royle-even graphs:** their count matches A000568, but they are not graphs whose vertex degrees are all even. Literal small-graph controls retain this distinction.
- **Other relational structures:** the root-selection proof transfers to graphs, digraphs, hypergraphs and colored structures with an exact deletion/attachment encoding. It does not authorize quotienting fixed arithmetic labels by an abstract shape symmetry.

## Verification and use

All new production/checker Python uses the standard library. The optional historical C comparisons and legacy C regression require a C compiler and GMP. Existing baseline files are loaded from their repository locations with explicit hashes; historical snapshots are not duplicated here.

```sh
python3 check_connections.py
python3 check_block_compression.py
node check_block_compression_independent.js
python3 verify.py --receipt verification.json
python3 check_pointed_observers.py
python3 check_generation.py --max-order 8 --selector neighbor_degrees
python3 check_repository_legacy_enumeration.py
python3 benchmark.py --receipt benchmarks.json
python3 benchmark_generation.py
```

Repeat Python checks under `python3 -O`. Exact receipts are identical; benchmark timings are expected to change. Receipts are generated locally rather than committed as raw audit data. The full labelled tournament universe was independently enumerated through order six. A separate auditor checked the generator through seven and inspected the producer's full order-eight checks; the order-eight check is not advertised as a second independent implementation run.

Scalar verification includes complete-prefix agreement with the prior counter through200, independent Fraction sums through50, and literal permutation/action controls. K4 and block-code results have separate independent finite implementations. The legacy regression deliberately states what it does not validate.

## Remaining boundary

The all-order count formula was already classical, and the repository already had faster counting work before Session86. The new contribution is the exact marked K4 bridge, its larger-composition obstruction, a usable compression generalization, tested execution improvements, and recovered implementation repairs. Individual unrestricted prime families, a universal smallest closed attachment carrier, and polynomial-time all-order class generation remain open. Aggregate symmetry counts do not supply missing class representatives.
