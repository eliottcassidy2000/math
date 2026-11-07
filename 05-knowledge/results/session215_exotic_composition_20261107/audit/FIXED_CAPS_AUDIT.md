# Independent audit: fixed cyclic caps

**Verdict: PASS WITH PRIORITY CORRECTION, applied.** Signed by the independent Session215 mathematical auditor. Exact checked files are recorded in FIXED_CAPS.sha256.

## Correction lineage

The initial presentation treated actual positive-context profile recovery as a new general contribution. Pinned prior `05-knowledge/results/proof_principles_20261103/profiles/PROOF.md`, R3, already proves recovery by actual tournament probes and compressed known-order recovery, crediting Session176. The final PROOF.md explicitly cites both sources and narrows the contribution to the alternative symmetric cyclic-cap family, its (m!)² diagonal, and explicit arbitrary-seed applications. Generic contextual equivalence, positive probes and parity compression are inherited. This correction affects attribution, not the formulas.

## Proof review

The inherited directed-triangle kernel has coefficient three on equal positive run counts, one on positive counts differing by at most one but not equal, and zero otherwise. Fixed cap orders bound the middle run count by one plus the smaller order. The finite cover functional and its Newton-transform coefficients therefore follow exactly. Ordinal multiplication of Q then yields mixed-seed products and repeated-seed exponential sums. These concern nonempty blocks and r>=1, not balanced growing caps or empty-block limits.

The singleton-cap formula has coefficient 3H(A)-Q_A(2), not 4H(A)-Q_A(2). The five-vertex profile (11,78,198,240,120) gives Q(2)=100, so the physical strong family has H=100^r+11^r and order 5r+2. The printed recurrence is correct. Strongness follows by detours around the nonempty cyclic quotient.

For two TT_m caps, middle run count m+1 forces both cap counts to be m, giving diagonal (m!)². Triangular inversion establishes the claimed observational equivalence, including different orders after zero extension. Known-order compression follows from monic order parity and Vandermonde interpolation in j²; the n=1,2 degenerate cases and n=3,4 scalar cases are sound. No minimum number of finite-set queries or isomorphism reconstruction is claimed.

## Reproduction and scope

Final normal and Python -O verifications reproduce receipt.json byte-for-byte. The verifier uses independent permutation/cut counts for all 75 labeled orders 1–4 and the five-vertex seed; independent forward subset DP verifies 160 physical compositions, with 76 full-profile and 79 compressed-recovery cases. The coefficient bridge (87,-189,72,36) agrees with the independent multi-exit EGF. All eight source manifests were checked by Git blob and SHA256 against pinned main 0f8750677dd20b594a0e032a3b56aa23f58a950c.

The result is an explicit compositional family and observer specialization, not a cofinite-spectrum theorem or a new generic profile-recovery theory. No numerical first-occurrence or literature-priority assertion is approved.
