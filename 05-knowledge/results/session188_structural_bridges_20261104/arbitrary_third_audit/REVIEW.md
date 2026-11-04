# Independent audit of arbitrary third-block separators and root families

**PASS on proof and execution review. No mathematical correction found.** The source packets are `arbitrary_third_block` and `arbitrary_third_root_family`. Their precise claims and ambient restrictions are retained. The receipt records the completed independent and producer execution checks.

## Six finite-support identities

The inherited cyclic old-run matrix and two-port word formula have already been independently audited. Here the finite-support extraction has an additional independent check: `audit.py` directly enumerates quotient-run words with fixed multiplicities, treating the two ports as distinct single occurrences. It sums against the exact TT2 and TT3 ordered-cover profiles to recover all six coefficient vectors, without importing the producer's rational polynomial algebra.

The tail bound is valid: deleting the two ports leaves at most three old walk segments; in each directed cyclic segment, the C count exceeds the A count by at most one. Thus total C runs are at most total A runs plus three, hence at most m+3. No longer C profile coordinate can affect the enumerated fields. The endpoint zero in each displayed vector is checked, but the proof does not rely on that zero to establish the tail bound.

## Universal separation, signs, and cubes

For every nontrivial C, q2/2 is the sum of odd products over all unordered nonempty bipartitions. There are 2^(|C|-1)-1 such partitions, which is odd, so q2 is two modulo four. q3 is divisible by six because permuting ordered path slots acts freely. Consequently -8H-14q2+4q3 is four modulo eight. The singleton boundary is separately -8 and is not incorrectly assigned the same residue.

For |C|>=3, cutting any Hamilton path at an internal edge injects the path/cut pairs into ordered two-path covers: concatenate the resulting two ordered paths to recover the source path and its cut. This proves q2>=(|C|-1)H>=2H. Order two is handled directly. The TT3 root sum is therefore at most -324H-552q3 and is strictly negative for every nontrivial C, whereas its singleton value is +242. No log-concavity claim is needed.

All root queries are uniform on every child, hence independent of the child-converse bits over the integers. Existing lifted child isolation decodes each nontrivial interior bit while all other children remain homogeneous to the ports. The adaptive/fixed counts are therefore one root query plus r/2r child queries for TT2, or two root queries plus r/2r for TT3. A singleton third child has no bit. The stated 16-state and 8-state counts and query budgets are correct. The exact zero of a single TT3/TT3/C3 probe is a legitimate hostile control and does not contradict the nonzero two-query sum.

## Root-containing arbitrary original-module families

The C3 child quotient forces every proper original module into one child, even when a child is itself nontrivial or has mixed internal incidences. A module meeting two children must contain every vertex of the third, because each omitted third vertex sees the first two parts oppositely. This then forces every remaining vertex of the first two children as well. Hence only the root spans children.

A selected root with odd parity therefore reverses the quotient once. Proper selected modules only change the child interiors. The original-module decomposition theorem gives hereditary Hamilton/profile domination. If the full scalar count increases, query (0,0) separates. If it stays equal, the equality characterization gives equality of all complete child q profiles on their fixed physical subsets. Every uniform root query is consequently the same as if the original child interiors had been retained and only the quotient conversed. The already proved fixed root battery then separates.

This argument allows arbitrary depth, number and overlap among original proper modules. It does not assert that crossing selected modules remain valid in sequential processing, nor all-pairs injectivity among arbitrary targets. The scalar query can be omitted in the known laminar class because each reversal remains legal and preserves scalar H. It cannot be omitted merely by assuming arbitrary original-module XOR is legal.

The frozen packets concern the whole cyclic core and separately justified ordered SCC padding. They do not assume transfer of their exact integer gaps through an arbitrary strongly connected exterior. Any later modulo-eight transfer corollary is an additional theorem, not an implicit premise of these proofs.

## Independent finite controls

The audit imports no producer helpers. Besides the literal run-word calculation, it reconstructs all physical cube states and every original-module XOR target from the saved masks. It independently enumerates each source's proper modules, checks the complete root-containing-family universe, and recomputes every declared query by endpoint subset dynamic programming. It checks the exact nontrivial-child query counts, root independence, distinct cube signatures, and scalar-equality/strict censuses. Both source verifiers are replayed normally and optimized, with byte and saved-record comparisons. All checks survive optimized Python; receipts bind both source proofs, both verifiers, and the independent audit script.

Final independent counts: 222 literal run-word terms recover all six coefficient forms; 126 actual query controls, 64 cube states with 440 query values, and all 476 root-containing families with 1,144 query values pass. Both source normal/optimized replays and both independent normal/optimized receipts are identical. The crossing TT4-third-child census is 96 equality and 160 strict cases; the crossing TT3-first-two-child census is 72 equality and 56 strict cases. Exact frozen source hashes are in `receipt.json`.
