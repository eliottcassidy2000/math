# Reproduce the source-cut and rearrangement packet

Run from this directory:

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

Only Python's standard library is used. The verifier writes RESULT.json as
well. All proof checks use exceptions, never removable assertions. Physical
source words are independently constructed and counted by subset dynamic
programming; these are word-universe checks, not an isomorphism-class census.

The independent audit and its separate implementation are under audit/.
Its frozen read receipt identifies the final proof and verifier hashes.
Run its checker with `python3 audit/independent.py` and
`python3 -O audit/independent.py`; its receipt is written beside the checker.

The all-order conclusions are proved in PROOF.md, not inferred from the finite
checks. Coverage of every sufficiently large odd integer remains open.

The six dependency documents are frozen source snapshots, not new claims or
fresh audits of every inherited theorem. SOURCE_MANIFEST.json records their
origin and exact bytes. Where an inherited snapshot refers to sibling files,
use its named original published packet. No such external file is needed to
run the new verifier.
