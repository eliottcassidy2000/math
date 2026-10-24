# Session50: parent-only marked-incidence growth counts

Read `THEOREM.md` for all-size proofs and the non-closure boundary.
`transfer.py` constructs order-N extension states, root transfer, endpoint
projection, and converse using only order-(N-1) isomorphism tests during the
prediction stage. It has no Burnside fixed-tournament summation and does not
use known counts as inputs. Known prefixes and child canonicalization occur
only in separate controls.

Run from any directory:

    python /path/to/transfer.py --max-order 7 --receipt /tmp/growth-normal.json
    python -O /path/to/transfer.py --max-order 7 --receipt /tmp/growth-optimized.json

The larger order-eight run uses the same prediction algorithm and all lower
controls; it does not claim an independent full order-eight word-cube census.
See `VALIDATION.md`, `receipt.json`, and `receipt-order8.json` for exact scope.
The inherited Markdown sources are archived byte-for-byte as `.txt`; original
commit/path/hash pins are in `source-pins.json`.

A closed recurrence for A000568 from bare metagraph invariants remains open.
The result here is an exact, auditable parent-incidence transfer that locates
precisely the gluing and stabilizer information such a compression must retain.
