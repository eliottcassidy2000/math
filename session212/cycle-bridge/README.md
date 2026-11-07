# Complete tournament Hamiltonian-cycle count support

**PROVED, independently audited:** Every nonnegative integer occurs as a tournament Hamiltonian-cycle count. Every positive integer occurs even in a strong tournament with a sole-exit vertex and no degree-defined gate. Cycles are counted modulo cyclic rotation. No literature-priority claim.

The coherent unified family is a transitive chain with two marked extra vertices. Binary orientations encode q; a root-orientation bit supplies tau; reversing one specified apex arc supplies a separate increment:

    HC(B-) = q+tau,       HC(B+) = q+tau+1.

At each N>=5 this realizes the complete interval 0 through 2^(N-4)+1 across the two variants, with all positive values strong. The guaranteed nongated positive interval is 1 through 2^(N-4), not the extra top value. Literal marked insertion yields C->2C and C->2C+1 inside the same strong nongated family. For every count at least two, these witnesses have feedback vertex number exactly two. This is an all-order construction, not a finite search extrapolation.

The endpoint companion shows that every positive integer is also a specified-start Hamiltonian-path count of a strong tournament, and a specified-start-and-end count of a strong tournament. These remain marked counts. The total Hamiltonian-path spectrum, including its conjectural completeness apart from 7 and 21, is not solved here.

## Files

- `PROOF.md`: Unified binary family, exact cycle partition, all-order support, order bound, nongated family, and marked growth grammar.
- `ENDPOINT_BRIDGE.md`: Inherited exact endpoint/sole-exit dictionary and independently derived rooted binary code, including strong cores and both fixed endpoints.
- `verify_binary_comb.py`, `RESULT.json`: 2,040 family instances through order twelve, 248 growth transitions, complete permutation checks through order seven, 75 inherited-gate controls, and explicit target witnesses.
- `verify_endpoint.py`, `ENDPOINT_RESULT.json`: Separate labeling and dense dynamic program, 511 endpoint instances, 1,022 cycle instances, 62 full-permutation controls through order eight, the m=0 boundary, the strong-core criterion, and specified endpoint checks.
- `SOURCES.json`, `sources/`: Exact inherited source pins from live main revision c137c698e632cb76e232a79b48133987755e080a.
- `audit/`: Independent audit supplied separately by the session reviewer.

## Replay

From this directory:

    python3 verify_binary_comb.py
    cmp RESULT.json RESULT.normal.json
    python3 -O verify_binary_comb.py
    cmp RESULT.json RESULT.normal.json
    python3 verify_endpoint.py > /tmp/endpoint-replay.json
    cmp ENDPOINT_RESULT.json /tmp/endpoint-replay.json
    python3 -O verify_endpoint.py > /tmp/endpoint-replay-O.json
    cmp ENDPOINT_RESULT.json /tmp/endpoint-replay-O.json

The primary verifier regenerates RESULT.json. Both verifiers use explicit exceptions rather than assertions. The finite universes test code and conventions; the all-order theorem follows from the displayed cycle bijection.

## Scope

This classifies possible **numerical counts over all orders**, not all realizing tournaments, optimal sizes, or full fixed-order spectra. Its sole-exit constructor is inherited from THM-4862; the stronger two-vertex gate correspondence is inherited from Session210; HC=7/21 witnesses already existed in Session211. Two independently developed binary endpoint/comb constructions in this session were identified and merged, rather than counted as separate results. No priority claim is made against the literature.
