# Session124: three-ear transfer and fixed radix-16 obstructions

This packet proves an **architecture restriction**, not a new forbidden
Hamiltonian-path count. In tournaments built from three fixed marked vertices
and a transitive spine, no fixed repeatable bank of length-two templates at
an actual seed H>1 can realize two distinct balanced radix-16 digit actions.
In particular a complete sixteen-digit bank cannot exist. Phase-dependent
banks that cannot repeat a physical template indefinitely are not excluded.

The complementary-template family also has a uniform exponential-polynomial
formula with a forced cubic term `(4/27)k^3`.

- Proof and exact scope: `THREE_EAR_RADIX16_OBSTRUCTION.md`.
- Independent mathematical audit: `audit/AUDIT.md`.
- Exact rational certificate banks and producer checkers: `code/`.
- Independently reconstructed transfer, spectrum, and example: `audit/`.
- Matched normal/optimized outputs: `receipts/`.
- Claim provenance and superseded intermediate conclusions: `CORRECTION_LINEAGE.md`.

## Replay

```sh
python /absolute/path/to/this/packet/verify_all.py
```

All five replay jobs require only Python's standard library (tested with
Python3.11). Each runs normally and under `-O`; outputs must match both each
other and their bundled receipts. No tests depend on removable `assert`
statements. Default replay does not write files or bytecode. All paths are
relative to source files, not to the invocation directory.

LP discovery, exploratory enumeration, SymPy, SciPy, and the original workspace
are not needed. The independent audit imports none of the producer modules.
The `--write-receipts` option intentionally regenerates receipt files and
should not be used to validate an unchanged packet.

The whitelist and SHA-256 manifest record curated payload files. Publication
transport batches are prepared separately from this portable payload.
