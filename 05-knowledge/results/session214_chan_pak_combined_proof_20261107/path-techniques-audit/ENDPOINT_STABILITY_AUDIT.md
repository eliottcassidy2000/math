# Session214 sole-exit endpoint stability audit

**Verdict: PASS.** The fixed-defect finite-skeleton proof, finite union of affine-dyadic rays, bounded-defect sparsity barrier and exact defect-three classification are valid under the explicitly stated positive endpoint and sole-exit hypotheses. No mathematical correction was required. This is not a theorem for arbitrary output endpoints or unrestricted total path counts.

## Frozen author bytes

Staging packet: `session214-path-techniques`.

- PROOF.md: `8dc84768893ecb8f0fcfcce4ed3ef89e16e98a37d0b38b6f22a4e5ea5b46c420`
- verify.py: `7588cc0496a0f5b13cc128a65a6362a7a2265b58db89cf8ac59137545eac47a3`
- RESULT.json and RESULT.optimized.json: `43942d03fcde1ec8e0df50ee4ad8c28f0770e21e5505550366bb78ffea38f9e9`
- RESULT.classes.json: `0e678ed7ce8f835eb5085e9a24e5100e5e60fd56c6d9629b6eaee569b8f29458`
- README.md: `050f02f9f1b04b1f475f4160e461b286deec9752dca699b85d3abfac4735edd0`
- WITNESSES.json: `e95271e15b8ad315f8d0627cef1f753ea97c9eb0ca52b31c45c45e4a7e621234`
- SOURCES.json: `970a32734abf22346cf8c4b077488d268e0728a4d707842b863ccc28a1de85bc`

Normal and optimized receipts agree. The class-representative receipt is supplementary inherited-input evidence; this audit's primary finite check regenerates all labeled tournaments through6 without relying on those representatives.

## Independent proof review

For the sole exit a->r and U=T-{a,r}, paths ending at a are exactly the valid words P,r,Q,a. Their failed end(P)->r junctions are exactly K, giving e_a=L-K. Moving a immediately before r injects these into paths ending in U; the missed words have precisely the same K obstruction. Paths ending at r are the H(U) words followed by a,r. Thus H-2e_a=H(U)+K.

The recalled exposed-block two-cover argument has at least two insertions for every old Hamilton word. Cutting at the exposed step does not require that internal step to be an arc. Removing its two boundary vertices recovers the old word, which verifies the required injection and the capacity floor. In particular K=0 iff at most one exit of r enters U.

For strong T, r must reach U's initial SCC directly. A singleton SCC after every exit can never belong to the K-cover's first path; deletion is bijective there. In an L-cover it can be placed in either labeled path, including an empty first path, so each such singleton doubles L. Removing it cannot empty the second path, whose start is an exit. Nontrivial SCCs are retained, and singleton deletion neither splits them nor merges distinct SCCs.

The skeleton order bound correctly separates nontrivial SCC vertices (at most H(U)<=d) from singleton components. Between the first and last exit-containing SCC, independent singleton assignments give at least 2^b distinct K-covers. Since K<=d-1, their number is logarithmically bounded. The j=1 case and d=1 logarithm convention are explicitly valid. Reverse insertion proves that the finite list of physical skeletons gives exact rays rather than merely necessary numerical forms.

For a nonstrong tournament, e_a>0 puts a in the final SCC. Its sole exit remains in that SCC. The external multiplier M divides d, and only physically attained multipliers are retained. The finite prefix test and logarithmic spectrum conclusion therefore do not assume arbitrary odd integers are attainable.

The defect-three proof exhausts (H(U),K)=(1,2),(3,0). The chain pair weights, short endpoint cases, triangle two-cover count12, and external multiplier3 all check. The exact endpoint spectrum is

    {1} union {3*2^q:q>=0} union {3*2^q-2:q>=1} union {5*2^q-2:q>=1}.

The output sole-exit restriction is essential: a sole-exit vertex in an input cycle compiler does not automatically imply this restriction on the remaining marked endpoint. The packet makes that distinction explicitly and leaves Chan–Pak completeness open.

## Independent exhaustive controls

`endpoint.cpp` enumerates every directed simple vertex word in every labeled tournament through order6. It obtains cover counts by partitioning the actual word lists, not by importing the author's endpoint/cover routines.

- 33,867 labeled tournaments.
- 32,110 sole-exit roots encountered; the two order2 roots have empty U and zero endpoint, outside the theorem's scope.
- Exact defect and endpoint identities checked at the remaining32,108 roots.
- 29,998 positive endpoint roots.
- 22,654 strong sole-exit roots checked for compression, preserved strongness and skeleton bound.
- 4,824 nontrivial compressions.
- 6,024 defect-three roots, all in the displayed exact spectrum.

Reproduce:

    g++ -O2 -std=c++17 -Wall -Wextra -Werror -pedantic endpoint.cpp -o /tmp/session214-endpoint-audit
    /tmp/session214-endpoint-audit > /tmp/session214-endpoint-audit.json
    cmp endpoint.json /tmp/session214-endpoint-audit.json

All validation gates remain active independently of optimization flags. No external literature-priority claim or multi-exit extension is certified by this audit.

## Final metadata and literal witness check

The final README explicitly restricts the displayed skeleton bound to strong T and separately invokes SCC factorization for arbitrary T. This fixes a summary-scope ambiguity without changing the proved theorem. The extended source manifest and all author checksums were rehashed. `endpoint_witness_audit.py` independently enumerates every permutation for the six supplied witnesses, confirming H=23,47,95,19,39,79 and the claimed sole-exit endpoint data.

    python3 endpoint_witness_audit.py ../path-techniques/WITNESSES.json > /tmp/session214-endpoint-witness-audit.json
    cmp endpoint_witness_audit.json /tmp/session214-endpoint-witness-audit.json
