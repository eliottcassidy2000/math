# Freeze review: mixed-(1,2) witness migration

Date: 2026-10-19. Both full receipts replayed under Python 3.12.14 with `-O`, byte-for-byte identical to the originals. Mathematical verdict: **PASS for the stated c≤3 maximum-span-root exclusion and its order dual.** This is a finite rooted-core certificate with an analytic arbitrary-exterior lifting proof. It is not a proof for c=4, all seven-thin posets, or the general 1/3–2/3 conjecture.

## Structural proof audit

1. **Majority order.** With no probability in the closed interval [1/3,2/3], each pair has a strict >2/3 direction. A directed triangle would force the sum of its three comparison probabilities above 2, although at most two of its directed comparisons hold in any extension. Therefore the tournament has no directed triangle and is transitive. The resulting total order extends the poset.
2. **Maximum-span criticality and root orientation.** A violation of D(x)⊆D(y) or U(y)⊆U(x) produces an incomparable pair spanning farther than the ordered roots. In the stipulated one-lower-slave morphology, y before x is itself impossible for a maximum-span root: b<y would put b before y before x, and b∥x would have greater span. Thus the x-before-y filter in the verifier is necessary, not an extra assumption. A critical-pair swapping injection is an alternative justification.
3. **Allowed core relations.** The listed optional directions are exhaustive. The prohibited reverse directions force a forbidden root/common comparability or x<y. Natural labels within Z and within the upper group are choices of names, not additional order relations. Enumerating increasing internal edges and taking transitive closure therefore covers every required core, allowing harmless duplication of isomorphism types.
4. **Projection.** Every exterior point is comparable to both roots, hence is either below both or above both; mixed directions would compare the roots. No exterior point is between the roots in the majority order. Projection preserves the root span and can only reduce other core spans, so every actual maximum-span majority projection appears in the enumerated set.
5. **Fibre product.** All lower exterior points must occur before the first root and all upper exterior points after the last root. Every lower exterior point is below every upper exterior point. For a fixed core extension, the remaining freedom is exactly an interleaving of A with its pre-root prefix and of its post-root suffix with B. Their normalized counts multiply, proving the fibre formula for arbitrary actual exteriors.
6. **Nonnegative increments.** A prefix letter cannot lie below any element of A, or it would lie below both roots and fall outside this core. Appending it after all A elements preserves each old interleaving. The dual argument prepends a suffix letter before all B elements. Thus the stated increment differences are nonnegative, including after any actual exterior constraints. No independence or unrestricted realizability of these increments is needed.
7. **Generator interpretation and orientation.** Telescoping along the two word trees gives a generator for each compatible initial word and terminal word. Its entries count exactly the core extensions satisfying those endpoints. Reversing an incomparable comparison changes the row 3N−2E to E−3N, as implemented. A nonpositive positive-weight sum of these rows in every column yields the ambient weighted ≤2/3 inequality and contradicts strict majority on all selected forward pairs.

## Exceptional inequality

The exceptional six-point core has 37 extensions. Its four comparison counts are 27,24,22,23, giving weighted numerator

    3·27 + 4·24 + 22 + 23 = 222 = 6·37.

The 3×9 deficit matrix is the one displayed in both proof notes; every coefficient is nonnegative. The matrix is for 6E minus the weighted numerator. The verifier's row convention is three times the negative of that deficit, so the normalization is consistent. The label-swapped exceptional case has the same proof. Bare-core equality establishes sharpness of this inequality.

The independent verifier checks recursive versus permutation enumeration of this core and exact whole-poset subset-DP evaluations for 200 actual exteriors of order at most 12. These finite exterior checks audit the analytic lifting argument; they do not substitute for it.

## Code and completeness audit

The main search generates all optional edge subsets and their closures; recursively enumerates every core extension; constructs all compatible generator columns; filters only the necessary maximum-span condition; and tests each selected comparison multiset against every column. The two exceptional weights are explicitly checked. The expected family counts and weight histograms are checked separately. The c=3 receipt stores per-core certificate-stream hashes rather than individual witnesses, but a full replay reconstructs and checks every witness before checking its receipt; a hash alone is never treated as proof.

All eleven main-script assertions and thirteen independent-script assertions were replaced with explicit conditional RuntimeError checks. No mathematical expression, enumeration order, witness choice, or serialization was changed. The replay driver checks structural equivalence of the original and normalized ASTs after this transformation; generated diagnostic-message wording is ignored for Python-version portability. Both normalized verifiers have zero Assert nodes, and deliberately cyclic input still triggers their guards under `python -O`.

## Exact family totals

| c | Cores | Bare core extensions | Admissible projected orders | Weight 1 | Weight 2 | Weight 3 | Weight 9 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 7 | 92 | 60 | 60 | 0 | 0 | 0 |
| 1 | 40 | 2,123 | 722 | 710 | 3 | 7 | 2 |
| 2 | 357 | 79,681 | 15,899 | 15,672 | 209 | 18 | 0 |
| 3 | 4,824 | 4,753,521 | 609,424 | 601,284 | 8,008 | 132 | 0 |

Together these are 5,228 cores, 4,835,417 bare extensions, and 626,105 admissible projected orders. Counts of naturally labeled cores are not counts of isomorphism classes.

## Preservation and reproduction

The original source files, proof notes, and receipts in the research directory are unchanged. The 10,496,542-byte main receipt is preserved byte-for-byte inside a 371,398-byte deterministic gzip. Original verifier and note copies are included as provenance. The frozen copy of the main note adds only reproduction guidance and two proof clarifications: root orientation and the increment injections.

Run `python replay_frozen_mixed12_20261019.py --check-only` for integrity and normalization checks. Run `python replay_frozen_mixed12_20261019.py` to recompute both certificates under `-O` in isolation and require exact byte equality with both preserved receipts. See `freeze_replay_receipt_20261019.json` for the completed replay record and `manifest.json` for file hashes.

The paper cited for the related balanced-triple and slave methods is Marcin Peczarski, [Comments on the Golden Partition Conjecture](https://cdm.ucalgary.ca/article/download/62349/46785/178291), 2017. That source was checked directly. This review does not claim a complete comparison with historical configuration tables or any priority result.
