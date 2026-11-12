# Session285 release receipt

2026-11-12. Final mathematical review: **PASS WITH MINOR SCOPE/TERMINOLOGY
REPAIRS**, all applied. See `audit/REVIEW.md` for the independently reviewed
universe and limitations. No claim of a complete Chan–Pak proof.

## Replays

The final `replay.py` run passed all ten normal/optimized pairs, each matching
its saved receipt byte-for-byte. Four are independent audit controls. All
compressed and raw catalog hashes passed. Environment versions and exact
receipt hashes are in `REPLAY.json`;52 source/receipt SHA256 identities are
in `audit/SOURCE_HASHES.json` (excluding itself and this release receipt).

The lexical search has34,864 read paths among34,887 candidates;7,509 seed
hits and10,069 nonseed expansion paths. Omitted23 paths are explicit.
Reachable commit history is shallow;4,541 commit subjects yielded157 hits.
These are search leads, not counts of deeply read sessions. Source-body
consultation and a deterministic lateral sample are separately identified.

## Repairs during review

- Monodromy group is abelian (nilpotency class≤1), with square-zero deviation
  algebra; the earlier 'two-step nilpotent' wording was unnecessarily ambiguous.
- Midpoint stationarity in the transitive attachment example is specific to
  the displayed symmetric direction, not asserted for every direction.
- Existing renewable six/nine banks are retained: what remains missing is a
  universal absolute-target section or renewable interval growth sufficient
  for all allowed odd targets, not renewable banks altogether.
- Unequal pair-orbit energies, not merely multiple orbits, generate variation.
- 'Five-state pentagon' avoids misreading it as a five-tail-label object.
- Frozen historical subject count is4,541, not the draft4,542.

No substantive theorem failure was found. Classical stress/polytope results
are credited; historical repository claims were not all independently
re-proved. New physical unit-switch witnesses were independently counted.

## Repository-wide check

`python3 agents/check_docs.py` was run and retains two pre-existing warnings:
`META-PATTERNS.md` has401 maintained lines against400, and startup surface
has180862 bytes against180000. This session did not edit the startup surface
or the warning-bearing file. These are not reported as a full aggregate pass.

## Publication discipline

The session publishes only its packet and `00-navigation/S285.md`.
Unrelated inbox changes in the base checkout were not staged or altered.
Initial and mathematical checkpoints used GitHub's fast-forward connector
because shell push lacked credentials. Final tree is assembled from the
current remote base and publication content, then fetched and compared by
Git blob identity before reporting completion. No force-push is used.
