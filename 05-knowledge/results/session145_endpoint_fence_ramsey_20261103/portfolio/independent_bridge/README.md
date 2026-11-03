# Independent Session145 bridge audit

**PASS for the final publication scope.** Read [AUDIT.md](AUDIT.md) and the resolved [scope correction lineage](CORRECTION_LINEAGE.md). The three checkers are self-contained and import no discovery implementations.

Requirements: Python 3.12 and SymPy 1.14.0 for the magic/fence and theta checkers. The Ramsey/knot checker uses only the Python standard library.

Run from this directory:

```sh
python3 verify_package.py
mkdir -p ../session145_audit_replay
python3 check_magic_fence.py ../session145_audit_replay/magic_fence.normal.json
python3 -O check_magic_fence.py ../session145_audit_replay/magic_fence.optimized.json
cmp ../session145_audit_replay/magic_fence.normal.json magic_fence.normal.json
cmp ../session145_audit_replay/magic_fence.optimized.json magic_fence.optimized.json
python3 check_theta.py ../session145_audit_replay/theta.normal.json
python3 -O check_theta.py ../session145_audit_replay/theta.optimized.json
cmp ../session145_audit_replay/theta.normal.json theta.normal.json
cmp ../session145_audit_replay/theta.optimized.json theta.optimized.json
python3 check_ramsey_knots.py > ../session145_audit_replay/ramsey.normal.json
python3 -O check_ramsey_knots.py > ../session145_audit_replay/ramsey.optimized.json
cmp ../session145_audit_replay/ramsey.normal.json ramsey.normal.json
cmp ../session145_audit_replay/ramsey.optimized.json ramsey.optimized.json
```

Replay outputs are written to the named sibling directory, leaving the sealed package unchanged. The integrity verifier enforces this exact artifact whitelist:

- `README.md`, `AUDIT.md`, `CORRECTION_LINEAGE.md`
- `check_magic_fence.py`, `check_theta.py`, `check_ramsey_knots.py`
- `magic_fence.normal.json`, `magic_fence.optimized.json`
- `theta.normal.json`, `theta.optimized.json`
- `ramsey.normal.json`, `ramsey.optimized.json`
- `SOURCE_EVIDENCE.json`, `AUDIT_RECEIPTS.json`
- `verify_package.py`, `MANIFEST.json`

Independent coverage includes the full 9! magic census, 2,025 rational parameter pairs, exact radical geometry, bounded-theta identities and cyclic witness, matchings through fourteen ports, 2,208,196 linear-forest pairs, Fox diagram constraints, direct Paley orbit/cycle counts, and 496 closed-component sidecar cases.

The source hashes document what was reviewed. They do not re-prove inherited elliptic-curve geometry, the classical fixed-side polygon theorem, or the external computer-assisted Ramsey bound. Those boundaries remain explicit in the audit.
