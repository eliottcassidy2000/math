# Independent review: exact one-port fibres and binary decoder

**Verdict: APPROVED in the stated fixed-label, fresh-query scope.** This
includes the complete one-port fibre over every paired-black quotient leaf,
not only a constructed subfamily, and the sharp second-stage modulo-four
query complexity. No unresolved mathematical defect remains.

## Proof and provenance

The C3/C5 binary-color reflection argument is valid in every uniform ambient
module. The chosen reflection may depend on the probe, which is consistent
with equality of the entire fixed-name response function. It does not make
the core arc arrays equal. Disjoint reversals are genuinely distinct named
arrays but remain isomorphic as unmarked tournaments.

The local11/13 and89/91 counts, dominance factorization and odd exterior
multiplier give exactly one independently readable bit per localized query.
Modulo-four answers are binary by tournament Hamilton-path parity; the
adaptive worst-case lower bound is therefore valid, even with more ports.
It is not a lower bound for full-integer answers or for computing H.

The exact-fibre extension has the necessary stronger input: current C4 in
`proof_principles_20261103/continuation/PROOF.md` proves equality of the
actual ordered named SCC partitions from equality of full one-port tables.
The reviewer fetched and read its positive-circulation proof at the pinned
commit286f5f4ddbdf25a74ecbfe74e8ef9a735a142393. Mutual refinement of the two
partitions justifies the claimed equality. Hence only the two triangle
orientations can vary on three-vertex SCCs. The24 labelled strong S4
centers have distinct one-port tables; localization divides out the same
positive exterior factor3^r and fixes the central labelled core. This
proves that the whole fibre has exactly2^r members. It does not prove
universal arbitrary-core two-port reconstruction.

Review caught and resolved a provenance issue before approval: Session27's
historical open module question was superseded by later repository work.
The final packet explicitly credits the stronger module invisibility,
visibility, SCC-boundary theorem and local C5 constants. Its exact-fibre
specialization and sharp oracle statement make no literature-priority claim.

## Independent verification

`independent_port_audit.py` counts Hamilton paths by literal backtracking,
without dynamic programming or producer imports. It checks:

-4,392 complete one-port mask comparisons for a C3 module in every ambient
  tournament with at most three exterior vertices;
- the exact local11/13 two-port pair;
- all32 leaf-word flip cubes through order10, with192 localized two-port
  counts and injective second-bit signatures;
- all24 labelled strong four-vertex cores and their16-entry one-port tables,
  independently recovering24 distinct tables.

Run `python3 independent_port_audit.py` and
`python3 -O independent_port_audit.py` in this audit directory. The outputs
are byte-identical to `independent_port_audit.json`.

The final complete producer checker was separately run in an isolated copy
under normal and optimized Python. Both JSON receipts and stdout streams
agree, and its RESULT.json matches the frozen producer receipt. This replay
includes677,320 explicit C3/C5 isomorphism checks, both local constants by
permutations, the nonmodule hostile, all stated SCC-word variants and the
new complete S4 sweep.

Reviewed SHA256 values:

- PROOF.md: 0ff8d999e25614d52290a9abf5228dd3741d6d46e14008152723c50bb0249e09
- verify.py: 2c5204bbfc7dacda1ea9262803a4d008ab236c28bc2c937308a9f06de0f11bff
- RESULT.json: e808dbe52dc4644918dc8f6e587cfce7cacd1ed3fbb5385a5c03d4a3a47360d3

## Boundaries

Each probe attaches fresh ports to the original core. Retaining a port is
a different, multiport experiment. Complete first-stage preprocessing can
cost2^n queries; the sharp r-query result concerns the second stage. The
2^r fibre consists of fixed-label presentations of one unmarked leaf, not
2^r metagraph nodes. No result recovers arbitrary wiring or settles the
global Hamilton-count spectrum or its eventual threshold.
