# Packaging signoff

Status: author packaging checks PASS; independent review recorded separately in REVIEW.md. No publication performed by the author.

- README.md contains the corrected proof, explicit boundary-input handling, finite algorithm, exact scope and inherited attribution.
- check.py uses only the Python standard library. It exhaustively enumerates labelled tournaments through order five, counts paths by permutations, computes strong connectivity independently by reachability, and checks the SCC recurrence.
- check.out is the frozen complete standard output. Replayed `python3 check.py` output matches it byte for byte.
- The proof's all-order statements are mathematical deductions; the verifier supplies finite controls and does not purport to computationally prove an all-order theorem.
- No stale working-tree edits, secrets, dependencies, generated binaries or external data are included.
- `manifest.json` records SHA-256 hashes of all release files other than itself.

Reproduce from the release directory:

    python3 check.py > replay.out
    cmp check.out replay.out

The successful normal Python run executes assertions. Optimized Python disables assertions and must not be used as the mathematical validation command.

Intended repository destination: `session207/inverse-size-spectrum/`. Publication awaits the parent coordinator's assigned lane and final acceptance.
