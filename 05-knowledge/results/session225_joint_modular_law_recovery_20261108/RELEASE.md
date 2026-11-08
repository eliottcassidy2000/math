# Recovery release and verification ledger

This packet is reconstructed, not restored. It makes no claim that unavailable original Session225 files or logs survived.

## Checkpoint 1

Commit: 0adbabe5845f44668da6b3e42fa07488957b113c
Parent: 6872236362fbb7077314e25f3a66a7d0070c5402
Tree: 5e7cad0e0d9dd90284de9003d4a3303c2f56270b

Exact immutable file readbacks were compared with the prepared content, not merely inspected by path:
- PROOF.md: de148345d85c5653d379cff64be485e489d6c3a5
- check.py: c6aff4f5c200b39cd3a603c0b7aaa1f366c2cc58
- README.md: b5d71c15e8ca653c057b95345e2496862a7cb633
- VERIFICATION.md: 34c0002edeed5783f6561815b16a39ee8ed49de9

Compare against the parent confirmed exactly four added files, all inside this recovery folder, with no deletions. Small individual blobs and SHA-only tree entries were used.

## Checkpoint 2

The second checkpoint adds the independent actual-artifact review and explicit prime-power kernels, this ledger, and separately scoped Session225 navigation. It preserves the existing graft recovery and leaves root descent under separate review. Its own immutable commit ID is provided by Git history rather than a self-referential placeholder in this file.

Before each branch movement, main is checked again and updated non-forcibly with expected_sha. Exact immutable readbacks are required before publication. No unrelated paths are overwritten and no original artifact path is reused.

## Verification limits

The fresh mathematical runs and failed optional launch are recorded in VERIFICATION.md. Repository-wide scripts and a full CI run were not executed. Empty GitHub combined-status/workflow responses are not a CI pass. This recovery creates no local repository/worktree; therefore no local dirty state or unpublished local commit is left by it.

Scope remains fixed odd q and fixed-length two-J words. Height laws, scalar coverage, bank completeness, priority, and the Chan-Pak cutoff are not claimed.
