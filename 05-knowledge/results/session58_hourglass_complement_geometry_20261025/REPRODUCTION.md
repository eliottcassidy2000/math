# Reproduction, mode awareness and provenance

First Session 58 checkpoint; session ongoing. All supplied executables use the Python standard library. Run in a disposable copy because producer and audit scripts write their JSON certificates beside the script. Do not treat regenerated certificates as an automatic publication update.

## Final producer checks

Run each command both normally and with `python3 -O` substituted for `python3`:

    python3 hourglasses/check.py
    python3 hourglasses/check_geometry.py
    python3 hourglass-fiber-product/check.py
    python3 complement-tournaments/verify.py
    python3 absolute-geometry/verify_bridges.py

The final producer checks use exception gates, not removable assertions. Normal and optimized outputs are compared in the preparation replay receipt. Enumerations are finite corroboration of the proofs, not proofs outside their ranges. These scripts may take several seconds, especially the independent exhaustive censuses.

## Independent checks and mutation controls

    python3 audit-hourglasses/check.py
    python3 audit-hourglasses/projective_check.py
    python3 audit-hourglasses/mitre_check.py
    python3 audit-complement-zeta/independent_audit.py
    python3 audit-complement-zeta/negative_controls.py
    python3 audit-complement-zeta/corrected_negative_controls.py

Run these in both modes too. The two negative-control drivers themselves run both normal and optimized child executions. Their expected result includes rejection of deliberately false checks, not a clean execution of the corrupted child.

IMPORTANT: `audit-complement-zeta/replays/absolute/verify_bridges.py` is intentionally the **historical assertion-based snapshot**. Its false assertion disappears under optimization; `negative_controls.py` demonstrates this old defect. It is not the current verifier. The corrected snapshot is `replays/absolute/corrected.py`, and the live producer is `absolute-geometry/verify_bridges.py`. The corrected negative control rejects the injected false parity check in both modes. This lineage is retained so the repair is auditable, not to endorse the obsolete snapshot.

The original hourglass audit harness referred to temporary absolute certificate locations. The publication edition reads the sibling producer certificates and writes beside itself. `PUBLICATION_ADAPTATIONS.json` records the exact before/after hashes. The historical prose in `audit-hourglasses/AUDIT.md` still describes the original harness. No mathematical test was changed by this portability adaptation.

## Provenance and integrity

`PROVENANCE.json` explicitly inventories every selected input with original path, size and hash. Temporary absolute paths inside source-pin files and audit reports are historical provenance, not links that consumers must resolve. Producer source pins identify the inherited repository baseline and external primary sources. Original source manifests are named `SOURCE_SHA256SUMS.txt`; they describe original bytes and locations, not the final portable package.

The publication map contains source/destination/bytes/sha256 for every delivered file. It is an external packaging artifact rather than a self-hashing member of its own inventory. Preparation and replay receipts likewise remain external review artifacts. All mathematical source packets are independently audited at their recorded hashes; the root synthesis still requires parent review before publication. No remote mutation is performed by this preparation.

Only selected UTF-8 text is included. No full articles, images, binaries, caches, worker conversations, private attachment identifiers, or user attachment text are published. All cited external sources are primary publications, authors' records, or authoritative primary mathematical reference entries; inherited attribution is identified where it was not freshly checked.
