# Independent review of the finite named marked-cover reduction

**PASS.** Frozen proof SHA256 `01b1f9b904d3db87f3e1882367ff9dcf045c22cd7c5ae35750507ce26424b5d4`; verifier `9a22fd3641a20bd386710a3ba1cc0340a0ed1a9baf5bb439f9786a6c688d7251`. No mathematical correction was needed.

## Free-slot action and exact physical substitution

With s distinct prescribed slots, permuting the remaining h-s slots preserves every endpoint prescription. The action is free because cover paths are nonempty and vertex-disjoint. Thus the exact count is divisible by `(h-s)!`, independently of tournament parity. Since eight divides 4!, the claimed cutoff `h<=s+3` follows. The same argument gives the stated general factorial-divisibility modulus cutoff.

A Hamilton word with two ports has at most four old contact roles. Splitting its old segments into maximal child runs and reading in global word order yields an ordered cover in each child, with no additional slot factor. A child run touching both ports can carry two endpoint roles; in particular the same physical name at both endpoints forces that run to be a singleton. Repeating a name in different slots is impossible. The marked tables retain these constraints rather than factorizing joint endpoint counts.

Fixing the run skeleton and all contact names makes every inter-child and port incidence explicit. Summing products of the corresponding marked child counts therefore gives the full physical response. At most four slots in total need endpoint prescriptions, giving the seven-path per-child and `3b+4` total-run sufficient bounds modulo eight.

## Coefficient refinement and literal exceptional-name cutoff

For a coefficient, expand outgoing contact factors before extraction. Choosing a constant term removes that endpoint prescription. Summing the unspecified contact name is exact: a fixed realized path has only one endpoint at that contact, so the summation adds no multiplicity. Slots carrying no selected physical name are again freely permutable. A physical name cannot belong to two different child runs; shared-index port roles either occupy the same compatible run or contribute zero. Thus the sharper number of distinguished slots is at most the number of distinct supported physical names in that child, as claimed.

For a literal query, partition covers by the slots containing exceptional names. There are at most s such slots. Permuting all other child-run slots preserves their port incidences because every vertex in those paths has the same joint two-port pattern; it also preserves all quotient joins. If fewer than s slots contain exceptions, the larger factorial still implies the displayed weaker divisibility. This validates the cutoff even when an exceptional name is internal rather than at a prescribed endpoint.

The result is a sufficient finite-data reduction retaining the complete bidegree-(2,2) polynomial. It neither says the number of labelled table entries is independent of child order nor proves that the retained data force noncancellation. The proof correctly makes no reconstruction or exact-collision claim.

## Independent computation

The producer enumerates permutations and valid cut sets. The independent audit instead enumerates unordered vertex-set partitions, chooses a directed path on each part, and assigns the designated slots by ordered selection of distinct component paths. It then multiplies by the factorial of the remaining slots. This reconstructs every nonzero marked-table entry and exactly matches the producer's complete table hashes/statistics for all 78 supplied children.

Additional explicit controls verify repeated-name prescriptions in different slots are zero, and a slot with identical start and end u has exactly the ordered-cover count of the vertex-deleted child in the remaining slots. These directly challenge the shared-name convention.

A separate C++ predecessor recurrence counts the number of first-child runs on each complete actual augmented Hamilton word. It agrees with every run coordinate of all 96 queries on six full four-variable faces, including constant nonzero backgrounds and shared physical names. All face coefficients are re-extracted, and every individually forbidden high-run contribution vanishes modulo eight. Positive high-run counts are retained, so the cutoff is actual divisibility rather than absence of contributing words.

Normal and optimized audit receipts are byte-identical. The audit imports no producer implementation, uses only the standard Python library and available C++17 compiler, and preserves the final proof/source pins. Its finite controls supplement the all-order slot-action and substitution arguments rather than replacing them.
