# Binary cycle-to-path arithmetic audit

This packet applies known path-count mechanisms to Session212's universal
cycle compiler. Its endpoint core has total count
`H=2h+D(h-1)`, where `D` is the total count of a residual one-apex tournament.
The hoped-for value `2h+1` occurs in this construction exactly for powers
of two `h`. Additive-vertex repairs are blocked by induced monotonicity.
The general odd Hamiltonian-path spectrum remains open.

- `PROOF.md`: all-order identities, operation obstructions, provenance and scope.
- `verify.py`: standalone exact checks with explicit failures under normal or optimized Python.
- `RESULT.json`: deterministic verification receipt.
- `SOURCES.json`: exact source paths and SHA-256 digests.

Reproduce with `python3 verify.py`, then with `python3 -O verify.py`.
Both outputs must equal `RESULT.json`.
