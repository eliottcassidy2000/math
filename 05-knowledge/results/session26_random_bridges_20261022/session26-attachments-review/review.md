# Independent hostile review: fair attachments to a transitive base

Verdict: PASS. No mathematical correction required in the frozen proof. This verdict covers its self-contained statements, not novelty, canon provenance, or repository publication.

## Binding and inspected sources

I read proof.md, audit.py, and audit.out in full from /workspace/shared/session26-attachments. The reviewed final proof SHA256 is `8d4f75f0bf851b9544e046f6dc8e2b11dda7fedc71e1980d9dc716eede24587d`. Audit source SHA256 is `dd95ac36eba1ea12e8c786516585b0367ba56a4eab87b411a02b50361efcf047`; output SHA256 is `e0c7a19fc451457f21dd1c7774045e979f5b4cc79d750ae29ceb7083972f03aa`. These hashes bind the reflection and negative-Hankel additions, not an earlier proof draft. Producer normal and optimized replays both match the frozen audit.out byte for byte.

## Mathematical review

1. The maximal-core-run decomposition is bijective. An ordered r-path cover has r-1 compulsory base slots and two optional end slots. Surjections count assignments of distinct base vertices to occupied slots; the increasing order within each slot is unique. Every candidate uses exactly 2r-2+u distinct cross edges, so its probability is 2^(-2r+2-u). No independence between different words is assumed. This establishes the exact transform including m=0.
2. The EGF uses the surjection EGF (exp(x)-1)^q. Its apparent division by y cancels because every cover is nonempty. Polynomial substitution gives the stated finite spectrum and rational OGF. The zero channel contributes at m=0 only. A tail recurrence without the zero root is not automatically an all-index recurrence; the proof explicitly respects this boundary.
3. At observation m, only r<=m+1 contributes, and the r=m+1 diagonal is m!/4^m. The cut dictionary for the backward-adjacency distribution is triangular. Reversing a word interchanges backward and forward adjacent pairs, giving B_b=B_(k-1-b). Thus ceil(k/2) initial means do suffice. The conclusion is profile recovery, not tournament reconstruction or a minimal sample claim.
4. Top profile formulas correctly count respectively singleton blocks, one pair, and either one triple or two disjoint pairs. The factor ordering the blocks is essential and present. The leading coefficient is k!/4^k; the next is (k+1)k!/4^k. Consequently the other lane's sole input k!/4^k*(k+1)^m+O_k(k^m) is sound for fixed k>=1. For k=1 the direct expression is (2^m+2+0^m)/4. There is no need for a connectivity hypothesis.
5. The first nonuniversal channel is the triangle term. The largest differing profile index s contributes degree s+1 with nonzero coefficient delta(p_s)/4^s; lower coordinates cannot cancel it. This proves eventual sign, not necessarily sign at every m.
6. Writing P(t)/t in the B basis makes t -> -1-t a signed exchange of complementary exponents. Substitution t=(z-1)/4 gives reflection about z=-1 with parity (-1)^(k+1). The extra even-k root multiplicity follows from parity plus the already present square factor. This is a polynomial symmetry, not dynamics or origin-reflection of the spectral support.
7. T4's coefficient vector is (-5,-9,6,22,15,3)/32. Its nonzero negative atom at j=1 persists after dropping m=0. The shifted Hankel matrix equals V diag(j*a_j) V^T and has determinant -16238475/1024. A nonnegative measure on [0,infinity) reproducing all positive-index moments would make this matrix positive semidefinite via integral x*f(x)^2 dmu, which is impossible. The proof does not confuse signed spectral weights with negative expectation values.

## Independent controls

The independent.py test does not import producer code. It computes path counts on every vertex subset and convolves over ordered subset partitions to get profiles, instead of using the producer's permutation-cut method. It expands the polynomial directly, instead of the producer's surjection formula. Expected Hamiltonian counts are independently obtained by enumerating all full vertex words, rejecting illegal deterministic adjacencies and assigning weight 2^(-number of cross adjacencies). This uses neither slot/run counts nor orientation-mask DP.

Passed controls:
- All 1,099 labelled cores of orders 1 through 5: subset-partition profile computation, leading channels, triangle coordinate, reflection, and half-prefix reconstruction.
- 2,348 direct weighted-word mean comparisons: all cores k<=4 with m=0..3, and all k=5 cores with m=0..1.
- T4 determinant independently calculated by the 120-term Leibniz expansion, not Gaussian elimination or Vandermonde evaluation.
- All eight C3-plus-one-vertex masks: H=3 occurs twice and H=5 occurs six times. The fair unconditional mean is 9/2. Strong-conditioned mean is 5, furnishing an explicit hostile to silently adding strong conditioning.
- The producer's nonisomorphic ordered-SCC collision is mathematically valid: K1▷C3 and C3▷K1 have different sorted degree lists but identical ordered profile (3,18,36,24). Hence tomography cannot be promoted to isomorphism classification.

The tests use exact Fraction arithmetic and explicit exceptions, not assertions or floating point. The finite tests supplement, not replace, the proof. No external search or novelty claim was undertaken. No Git or repository writes were made.

Final replay: independent normal and python -O outputs byte-match. Independent source SHA256: `fb45b22d0bdf12f85d123b22601409340ad57f1135468b3e831c7e4e4dbfde97`; independent output SHA256: `9c7bc30a098fa6298e3ac8f7624c7032f49773559a7423a9028e1c8144e3826a`. A final source rehash confirms the reviewed proof remained unchanged.

## Final publication status-only rebind

Verified final proof SHA256: `1af17f5b1954af08c9e4fb65d71af3b3e764e0b1b0dee926e47b2e347a8b43a2`. The only change from the reviewed proof is replacement of its stale candidate-status paragraph by the PROVED/research-packet status paragraph. Reversing that exact replacement reproduces the previously reviewed SHA256 `8d4f75f0bf851b9544e046f6dc8e2b11dda7fedc71e1980d9dc716eede24587d`, establishing byte-for-byte identity of all mathematical content. Audit source and frozen audit output hashes remain unchanged. PASS carries forward to this final publication version; no new research or mathematical change was involved. Earlier independent output intentionally records the original reviewed hash.

## Final Markdown-only rebind

Verified final proof SHA256: `82776d149b480cbe6068664864e2bea807236de424f758d8c3bb00b8ab37673c`. The only change from the preceding reviewed publication version is replacing `[2(k-2)!/4^{k-2}]` by `(2(k-2)!/4^{k-2})` to prevent accidental Markdown link interpretation. Reversing that single replacement reproduces SHA256 `1af17f5b1954af08c9e4fb65d71af3b3e764e0b1b0dee926e47b2e347a8b43a2` exactly. This is a grouping-delimiter typography repair with no mathematical change. PASS carries forward; no computation rerun was needed.
