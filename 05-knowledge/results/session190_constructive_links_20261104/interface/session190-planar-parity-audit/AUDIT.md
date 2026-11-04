# Independent audit: six-port Hamiltonian parity gadget

**Verdict: PASS.** Audited the final revision 2 of `session190-interface` on 4 November 2026. No mathematical blocker found. The universal assertions are supported by explicit proofs, not extrapolation from the finite cases. This audit makes no literature-priority claim.

## Exact subject and necessary hypotheses

The count is the number of **unoriented spanning Hamiltonian cycles**. Every interior vertex is mandatory. A and B are simple connected graphs with the same six distinct boundary vertices in the prescribed order. Each exterior intersects the patch exactly in those vertices, has an edge set disjoint from the patch, and has at least one interior vertex. The planar integer identity further requires an opposite-disk embedding with that boundary order.

These qualifications matter. Ordinary nonspanning rim cycles do not count. The hub is not optional. The formal difference between the profiles is not claimed to be the profile of another graph. All-five-unit-response equality modulo two is not integer equality and does not establish Hamiltonian existence.

## Independent complete enumeration and geometry

`audit.py` was written independently, with no producer imports. It enumerates every degree-valid edge subset by branching over all graph edges and pruning only impossible degree bounds. It does not begin by assuming a selected boundary degree pattern or forcing the subdivisions in code. It then removes cyclic components and records every remaining subset and its pairing.

The complete output is:

- B has 18 vertices, 18 edges, and exactly two partial Hamiltonian forests. Its entire profile is `(1,0,0,1,0)` in the producer's matching order.
- A has 19 vertices, 24 edges, and exactly nine partial Hamiltonian forests. Its entire profile is `(3,1,1,3,1)`.
- Every one of these forests has boundary degrees `(1,1,1,1,1,1)`. No other degree sector occurs. There are no cycle-containing degree-valid candidates.

The elementary explanation was checked separately. Each mandatory degree-two subdivision forces its two spoke edges. The ring then has exactly its two perfect matchings. In A, the hub chooses an even-indexed and an odd-indexed rim vertex. There are nine such pairs, each leaving a unique rim perfect matching. Six adjacent choices give three copies each of M0 and M3; the three opposite choices give M1, M2, and M4 once each. These choices exhaust all selected subgraphs.

The producer's integer-coordinate drawing was checked independently with exact segment-parameter and polygon-halfplane tests, including incident-edge overlap. A has 276 edge-pair checks; B has 153. All nonboundary vertices lie strictly in the boundary hexagon, and no forbidden intersection occurs. Graph simplicity and connectedness also follow directly from the specified ring and spokes.

## Matching matrix and Smith invariants

In the stated order, direct matching union gives

    K = 0 0 0 1 0
        0 0 1 0 1
        0 1 0 0 1
        1 0 0 0 0
        0 1 1 0 0.

Independent calculations give determinant −2 and rank four over F2. The leading 4-by-4 minor has determinant one. Therefore the fourth determinantal divisor is one and the Smith invariants are exactly `(1,1,1,1,2)`, not merely a rank-based guess.

The full six-port alphabet was generated independently: 365 live states. The other 360 states form a permutation block; there are 362 rows of degree one and three rows of degree two. Its F2 rank is 364, and its integral Smith invariants are 364 copies of one and one copy of two.

## Universal planar identity and physical witnesses

The full integer profile difference is

    d = (2,1,1,2,1),       Kd = (2,2,2,2,2).

For any permitted planar exterior X, restriction of a spanning Hamiltonian cycle to either side is a spanning partial forest. A closed restriction could not cover the other side's interior. Forced boundary degree one means that only X's all-ones sector contributes. Conversely, compatible forests whose pairings form one cycle glue to exactly one unoriented spanning Hamiltonian cycle. The inverse restriction map prevents orientation or multiplicity factors.

Consequently, with S(X) the sum of all five all-ones coordinates,

    h(X union A) − h(X union B) = 2 S(X).

This proof applies to every permitted planar exterior, including those having additional boundary degree sectors. It is not an empirical generalization from the test set.

An independent vertex-order Hamiltonian DFS, which does not use the gluing formula or the forest enumeration, verifies the five unit exterior A/B counts:

    M0: 3/1; M1: 2/0; M2: 2/0; M3: 3/1; M4: 2/0.

Thus the profiles differ modulo two, every exterior sees the same parity, and an explicit planar exterior distinguishes Hamiltonian existence. Forty additional planar exteriors with boundary edges, optional hubs, and deleted edges were checked by the same vertex-order DFS. Every one satisfies the integer identity. Their all-ones masses range from zero to 63.

The producer verifier was separately replayed in normal and optimized modes and reproduced its frozen output byte-for-byte. It covers all 365 live planar unit states, ten additional crossing pairing probes, and exterior copies of A and B: 377 exterior cases, each glued to both gadgets.

## Arbitrary-exterior strengthening

The stronger parity statement is also correct. An arbitrary exterior can induce any of the 15 perfect matchings, but still only the all-ones degree sector contributes. Independent enumeration and independent physical vertex-order tests give these exact differences:

    01|23|45: 2    01|24|35: 4    01|25|34: 2
    02|13|45: 4    02|14|35: 6    02|15|34: 4
    03|12|45: 2    03|14|25: 4    03|15|24: 6
    04|12|35: 4    04|13|25: 6    04|15|23: 4
    05|12|34: 2    05|13|24: 4    05|14|23: 2.

The revision-2 nonenumerative proof was independently checked. M1, M2, and M4 are the three odd bijections from `{0,2,4}` to `{1,3,5}`. A bipartite exterior matching has either zero or two connected partners among them, by relative permutation parity. A nonbipartite matching has exactly one even-even edge, one odd-odd edge, and one mixed edge. Exactly one of the three special matchings shares that mixed edge. The other two have no common edge and therefore form a single alternating six-cycle. The parity contribution is always even.

Hence parity equality holds for every arbitrary exterior satisfying the intersection, edge-disjointness, and nonempty-interior assumptions. The planar formula with a constant coefficient two must not be extended to arbitrary exteriors: the displayed differences are not constant. Revision 2 keeps this distinction explicit.

## Recovery and general integer-observability statements

All 1,024 profile vectors modulo four were enumerated independently. Every response residue determines a unique profile parity, and the stated half-sum formulas recover it. There are 512 distinct response residues. Modulo two, the lost direction is exactly `(0,1,1,0,1)`. The specified divided difference recovers the parity of the five-coordinate total, which resolves that missing bit; the other 360 coordinates already have permutation probes.

The general exact sequence

    0 -> ker(R)/p ker(R) -> ker(R mod p) -> coker(R)[p] -> 0

was checked directly, including well-definedness, surjectivity, and the integral-kernel term. The full-column-rank precision threshold `1 + max v_p(d_i)` follows correctly in unimodular Smith coordinates. The theorem is about arbitrary integer vectors, not an unproved physical realizability assertion.

For the order-three Walsh comparison, the polynomial was independently evaluated on all eight labelled tournaments: six values are one and two values are minus three. The matrix `[[1,1],[1,-3]]` has Smith invariants `(1,4)`. All 64 coefficient vectors modulo eight confirm coefficient-parity recovery. The vectors `(0,0)` and `(1,-1)` have identical responses modulo four and different coefficient parities. The claimed contrast with the six-port factor two is correct. The all-order Walsh discussion is conditional on the explicitly cited inherited invariant-basis theorem; this audit does not reprove that separate theorem. Its Fourier-denominator argument is consistent with that basis convention.

## Frozen producer hashes

All entries in the producer's revision-2 manifest were checked:

- `PROOF.md`: `f33466652c3d8a88b437185930da9034869e95327ba480c157d1360e4f3c7566`
- `verify.py`: `1831f5eec0655ac178e024a654888893f65a152a3aa9342f9336d7d0b958e327`
- `RESULTS.json`: `c8a3337459b6484ea61d6b99e6661be56824f45f06f1dedb8de344189b53164c`
- `RESULTS.optimized.json`: `c8a3337459b6484ea61d6b99e6661be56824f45f06f1dedb8de344189b53164c`
- `README.md`: `0fc1208fea19bcd0b209770460ec73d9005590b427102a1e8d90761da4886758`

## Independent audit artifacts and reproduction

Run from the workspace root:

    python3 session190-planar-parity-audit/audit.py > /tmp/parity-audit.json
    python3 -O session190-planar-parity-audit/audit.py > /tmp/parity-audit-optimized.json
    cmp /tmp/parity-audit.json /tmp/parity-audit-optimized.json
    cmp /tmp/parity-audit.json session190-planar-parity-audit/RESULTS.json

The script uses only the Python standard library. All test failures use explicit exceptions and survive optimization.

- `audit.py`: `6affd3549996e067584746e48a90f288d3c58e43b6845998a21928930e4303e4`
- `RESULTS.json` and `RESULTS.optimized.json`: `3929723a3d410b0f14634ee808ad7eb3f77a3f31b7897cede6a401660c97d189`
- `PRODUCER_REPLAY.json` and `PRODUCER_REPLAY.optimized.json`: `c8a3337459b6484ea61d6b99e6661be56824f45f06f1dedb8de344189b53164c`

The accompanying `FROZEN_MANIFEST.json` hashes this audit report and all listed audit artifacts.
