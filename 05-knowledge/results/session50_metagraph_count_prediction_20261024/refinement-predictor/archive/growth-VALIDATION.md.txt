# Exact validation and replay boundary

Status: FINITE-EXACT. Python standard library only (Python3.10 or later).
No randomized tests, external packages, Burnside fixed-tournament sum, or
child-class count input is used. All validation is explicit exceptions,
not Python assertions.

## Universes and outputs

The predictor starts from the unique empty parent and builds the complete
pointed extension inventory at every successive order. Each child-class
representative used at the next stage is selected from a distinct root-transfer
row, not obtained by canonicalizing the current child order. The runtime
records the largest canonicalized order during prediction; it is exactly
N-1 at each order N. Child isomorphism canonicalization is used only by the
separate audit routine.

| N | full pointed states | endpoint pointed states | fixed-path words | A000568 | self-converse | merged | symmetry defect |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 |
| 2 | 2 | 1 | 1 | 1 | 1 | 1 | 0 |
| 3 | 4 | 2 | 2 | 2 | 2 | 2 | 2 |
| 4 | 12 | 7 | 8 | 4 | 2 | 3 | 4 |
| 5 | 48 | 33 | 64 | 12 | 8 | 10 | 12 |
| 6 | 296 | 236 | 1024 | 56 | 12 | 34 | 40 |
| 7 | 3040 | 2681 | 32768 | 456 | 88 | 272 | 152 |
| 8 | 54256 | 50709 | 2097152 | 6880 | 176 | 3528 | 784 |

The order-eight fixed-path-word value is recovered by summing the proved
parent endpoint fibres, not by an independent scan of all2,097,152 words.
The independent full cube scan stops at order7. The rows are outputs;
known A000568 and self-converse prefixes are checked afterward, never supplied
to the predictor.

## Exact controls

All marked states through order8 are checked for:

- Parent subset orbit-stabilizer, endpoint lift divisibility, and deletion-atlas membership.
- D row sums, D squared=N D, equal rows within a root-transfer block.
- Stabilizer/isotropy recovery and detailed balance.
- Pointed converse involution and its commutation with D.
- Endpoint projector R squared=R, K squared=R, and RK=KR=K.
- Full and endpoint trace counts, converse trace, and merged integrality.
- Individual observer endpoint weights L/h and global word/right-lift masses.
- Directed-refinement upper bounds, odd root-orbit sizes, and explicit minimal hostile controls.

Independent paths in the same executable are:

- Every fixed-path word through order7, verifying each endpoint-pointed fibre
  and the unpointed class set via separate child-order canonicalization.
- Every pointed extension through order6 under all child vertex permutations,
  directly computing Hamiltonian end counts and child automorphism groups.
- Labelled orbit-mass checksum through order6.
- Invariant-refinement child canonicalization through order7, compared with
  the predicted transfer partition.

The separate `refinement_hostile.py` imports no transfer implementation. It
checks all5040 permutations of the regular-seven witness and explicitly
verifies its automorphism group, vertex orbits, and trace excess. It also
checks clone degree cells, neighbor counts, and canonical homogeneous-pair
chains for inherited inflations at multipliers1 through6 (orders7 through42).
Its normal and optimized stdout agree byte-for-byte.

These controls certify the stated universes. They do not independently
exhaust all labelled order-eight tournaments, do not prove a compact
metagraph recurrence, and are not performance benchmarks.

## Recorded normal/optimized agreement

The final source is `transfer.py` SHA256
`fd3db75c8425ff4fc2946835e7e4e55ffba905ae4d66fb1a6efd5dda61c9c4ba`.

- Default maximum order7: 172,308 explicit checks.
- Extended maximum order8: 2,641,571 explicit checks.
- `normal.out` equals `optimized.out` byte-for-byte.
- `receipt.json` equals `optimized-receipt.json` byte-for-byte.
- `order8.out` equals `optimized-order8.out` byte-for-byte.
- `receipt-order8.json` equals `optimized-receipt-order8.json` byte-for-byte.

Receipts contain deterministic complete row summaries and exact kernel
hashes. They intentionally omit machine-dependent timings and temporary paths.

## Manifest-only relocation

The manifest uses literal-byte SHA256 hashes. `source-pins.json` authenticates
all inherited Markdown copies preserved unchanged as `.txt` against the
pinned source commit. No repository files are imported or read during replay.

    python replay.py --manifest-only --max-order 7
    python replay.py --manifest-only --max-order 8

The replay verifies every manifest file, copies **only** those files plus the
manifest into a fresh temporary directory, runs normal and optimized Python
from that directory, compares stdout and JSON receipts against the sealed
outputs, and runs the independent regular-seven/inflation hostile checker.
Its terminal JSON attestation includes the exact manifest hash. Replay logs
are not themselves inside the manifest, avoiding a self-referential hash.

## Scope of compression

At order7, stable refinement supplies exact root-orbit sizes for3037/3040
states, but its trace upper bound3206 differs from the true3192. At order8,
it is exact for54250/54256 states and certifies53882 singleton roots; its
trace upper bound55068 differs from55040. The infinite inflation corollary
explains why this tempting simplification cannot be exact at every order.
The retained deletion-overlap transport resolves the discrepancy without
assuming a pre-existing child classification.
