# Reviewed release and scope

The physical compiler, complete target-dependent inverse, finite normal-menu
reduction, graph-order bound, and conditional implication are independently
reviewed. The all-target arithmetic occupancy lemma remains OPEN.

- `PROOF.md` is the promoted statement. `originals/` and `PROMOTION.diff`
  preserve its exact unreviewed predecessor and status/session-label changes.
- `audit/REPORT.md` reviews the mathematics and independent checks, not just
  the finite occupancy result. `audit/PROMOTION.md` verifies release bytes.
- `audit/REPLAY_100M_ADDENDUM.md` separately records the completed replay of
  the reviewed code through100,000,001. It tests19,856,206,889 index pairs and
  leaves29 length-eight gaps, largest2647. It is bounded evidence only.
- `explorations/` preserves author-only discarded routes and the predecessor
  scan, with its original limitations. Those older sources are not promoted
  merely because the reviewed successor independently reproduces the same
  numerical100m result. The new compact12/26-template companion is absent.
- No existing repository indexes or unrelated checkout files were changed.

Routine reproduction commands are in `README.md`. All file checksums below
this directory are in `MANIFEST.sha256`, excluding that manifest itself.
