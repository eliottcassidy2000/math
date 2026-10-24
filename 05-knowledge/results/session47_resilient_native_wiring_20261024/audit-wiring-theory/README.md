# Session 47 independent wiring-theory audit

PASS: minimax/cut/MST and neutral-budget quotient proofs; helper supply-to-route
interface; explicit suspension class-descent obstruction; all-size triangle
identities; and the two-port determinant/closed-walk refinement.

Read `PROOF_AUDIT.md` for the proof review, exact scopes, and qualifications.
No native guard transfer, new helper supply, full metagraph enumeration, or
universal numeric move threshold is inferred.

Independent results in each normal and optimized run:

- `independent_check.py`: 2,340,879 checks.
- `independent_generating.py`: 6,740 checks.
- Total: 2,347,619 checks, with matching mode-independent results.

Producer replays separately reproduce 2,386,018 / 12,118 / 3,019 checks.
The independent code imports no producer code. It uses literal simple paths,
all labelled small posets, product-group charges, full owner-permutation
canonicalization, induced-degree triangle counts, literal directed cycles,
and Leibniz determinant expansion as different verification routes.

`INPUT_PINS.json` and `GENERATING_INPUT_PINS.json` authenticate the reviewed
snapshots under `inputs/`. `MANIFEST.json` authenticates this delivered packet.
Run `python verify_manifest.py` or `python -O verify_manifest.py` to check it.

The mathematical producer files and inherited source bytes remain unchanged.
The final portability refresh uses .txt archival names and updated provenance;
see PORTABILITY_UPDATE.json. Both independent normal/-O receipts and the
bottleneck replays were regenerated. The global-observer audit was untouched.
No Git operations or publication actions were performed in this lane.
