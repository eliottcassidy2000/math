# Endpoint/deletion response for square-sum arrivals

Read `PROOF.md` for the exact all-component transfer, deletion-card update, sharp finite-horizon truncation, four-port K3 response inversion, and bounded square-sum consequences.

Canonical reproducibility files:
- `verify.py`: exact transfer, independent edge-subset controls, finite experiments.
- `direct.py`: independent vertex-order path enumeration.
- `RESULTS.json` and `RESULTS.optimized.json`: deterministic matching outputs.
- `verify_structural.py`, `STRUCTURAL.json`, and `STRUCTURAL.optimized.json`: human obstruction adjacency checks and positive witness.
- `FROZEN_MANIFEST.json`: canonical file hashes.

Run `python3 verify.py > RESULTS.json`, repeat with `python3 -O`, and compare. Run `python3 verify_structural.py > STRUCTURAL.json` for the structural certificate. No external packages or network are used. `prototype.py` and `explore.py`, if present in the working directory, are development scratch files and are not part of the canonical packet.

New concrete findings: Q29's 19 Hamiltonian paths all use the new vertex internally despite Q28 having 52 paths; Q25 minus 3 is traceable but Q26 minus 3 is not; full old endpoint counts cannot determine next-step traceability. Infinite endpoint/deletion statements remain open here.
