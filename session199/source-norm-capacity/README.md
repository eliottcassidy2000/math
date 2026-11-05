# Source norm and capacity

Sharp all-order source-cover minima, a quadratic C2 norm chart, and a
strictly scoped density-zero theorem for finite common-source C2/C3 families.
The unrestricted odd spectrum remains open.

- `PROOF.md`: complete proofs and scope boundaries.
- `verify.py`: standalone standard-library verification.
- `NORMAL.json`, `OPTIMIZED.json`: matching exact receipts.
- `dependencies/`: unchanged inherited proof snapshots.
- `audit/`: independent proof review and exact reconstruction.

Run `python3 verify.py` and `python3 -O verify.py` from this directory.
The independently reviewed proof and verifier hashes are pinned in
`audit/INDEPENDENT.json`; `MANIFEST.sha256` pins the original release packet.
