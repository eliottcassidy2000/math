# Balanced-source normal form

This packet proves a canonical graph inverse for balanced sources, an all-order
physical rerooting identity, and an exact residual family for C1 graphs that
have no balanced root. It does not prove universal odd coverage. The proposed
cofinite arithmetic lemma is explicitly open.

Run from this directory:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json
    sha256sum -c MANIFEST.sha256

The verifier is standalone, uses only Python's standard library, and does not
read a checkout or contact the network. The large exact computation compares
all projected profiles, not merely scalar supports. Its 42 balanced gaps
through50,001 and40 repairs are finite results only.

`PROOF.md` gives domains, physical constructions, the root-choice boundary,
explicit graph isomorphism, residual scalar formula, and the open arithmetic
statement. `REPAIRS.json` contains the small exceptional source certificates.
`DEPENDENCIES.md` records pinned provenance. The independent audit is published
alongside this packet, rather than mixed into its verifier's input files.
