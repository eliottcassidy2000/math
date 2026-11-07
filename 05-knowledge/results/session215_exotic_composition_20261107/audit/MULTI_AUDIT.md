# Independent audit: arbitrary-apex positive block compiler

**Verdict: PASS.** Mathematical proofs, deterministic receipt, independent permutation controls and source pins verified. No substantive mathematical corrections required. Audit signed by the independent Session215 mathematical auditor.

The checked artifacts are named and hashed in MULTI.sha256. Source manifest SOURCES.json was checked against all eight actual Git blobs at repository commit 0f8750677dd20b594a0e032a3b56aa23f58a950c, including SHA256 digests.

## Proof review

The apex-order/run bijection correctly delays the final run-edge test, permits repair of a temporarily invalid last edge, and handles empty internal runs using actual apex arcs. The dimension k!·2·3^k is the stated redundant full space, not a reachable-state minimum. At k=0 the unique empty Hamilton word and the chain count one agree. At empty chain length and positive k, accepted apex orders count precisely H(A).

For each fixed letter, slot tests and writes are coordinate-local. Every slot operator is idempotent, including the zero image of a forbidden empty slot. Distinct slot operators commute. Expanding words by their used-slot subsets proves the nonnegative onto/Stirling formula, including m=0 with 0^0=1. Simultaneous image/kernel decomposition proves the squarefree annihilator with roots 0,...,k+1, and excludes polynomial factors in each homogeneous-block exponent. This does not assert commutativity for different letters.

The two-block EGF has the correct initial/final factors and internal descent factors. It is correctly credited as a specialization of existing substitution/profile theory. The cyclic-three-apex coefficient vector (87,-189,72,36) is correct for m>=1; m=0 must be separate and is. The all-k highest coefficient (k!)² follows by assigning the initial k vertices one each to slots 1,...,k and then allowing the final block all k+1 slots. Any missing initial slot removes the highest eigenmode. This separates parameterized sequences, not individual numerical realizability.

## Computation review

- Normal and Python -O executions of verify.py are byte-identical to RESULT.json. The final receipt records 7543 exact subset-DP controls, 384 block-expansion controls, 43670 idempotence/commutation checks, 1200 EGF controls, 20 explicit-formula controls and 75 highest-clock controls.
- check_multi_independent.py uses complete permutation enumeration on 75 deterministic eight-vertex tournaments with 0,...,4 apexes. It does not use subset DP or the run representation for actual H. All match the transfer.
- check_multi_leading.py independently extracts the highest mode for 60 deterministic random apex tournaments at k=1,...,6, using the EGF and annihilating lower roots. All match (k!)².
- The independent scripts find the sibling multi-exit packet when published; workspace fallbacks support local replay.

No all-odd spectrum, all-k inverse synthesis, or literature-priority claim is certified. The finite checks support the written all-k proofs rather than replacing them.

Reproduce from this audit folder:

    python3 check_multi_independent.py
    python3 check_multi_leading.py
    python3 ../multi-exit/verify.py
    python3 -O ../multi-exit/verify.py
