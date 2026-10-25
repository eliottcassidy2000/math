# Reproduction, mode awareness and provenance

Second Session 58 checkpoint; session ongoing. The first-checkpoint instructions below remain valid. All supplied executables use the Python standard library. Run in a disposable copy because producer and audit scripts write their JSON certificates beside the script. Do not treat regenerated certificates as an automatic publication update.

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


## Second-checkpoint checks

Run each normally and with `python3 -O`:

    python3 hourglass-prime-support/check.py
    python3 audit-prime-support/audit.py
    python3 tropical-gap-triples/check.py
    python3 audit-tropical-triples/check.py
    python3 observer-profile-algebra/analyze.py 6
    python3 observer-profile-algebra/profile_collision.py
    python3 observer-profile-algebra/direct_check.py
    python3 observer-profile-algebra/family_check.py
    python3 observer-profile-algebra/quadratic_check.py
    python3 observer-profile-audit/audit.py
    python3 observer-profile-audit/audit_walsh.py
    python3 observer-profile-audit/audit_quadratic.py
    python3 audit-profile-algebra/audit.py
    python3 audit-profile-algebra/check_producer.py
    python3 audit-profile-algebra/family_check.py
    python3 audit-profile-algebra/negative_controls.py
    python3 audit-profile-algebra/quadratic_audit.py
    python3 diagrams/check.py

Run the observer source checks in this order so generated input certificates are available. The original SVGs can be regenerated with `python3 diagrams/draw.py`; no raster renderer is needed for the exact validation. Published SVGs are UTF-8 text, and the PNG review previews are excluded.

Second-checkpoint checks are replayed in an isolated packet copy in both modes. `SECOND_PROVENANCE.json` pins the explicit input allowlists; `SECOND_PUBLICATION_ADAPTATIONS.json` records all original/publication hash pairs for path-only edits and the prime-support source-copy clarification. The primary PDF and extracted full text, fetched metadata, caches and previews do not ship.

The first observer audit preserves two historical findings: a graph-count placeholder and an external helper dependency. Both were corrected in the final producer before this freeze: the unused graph-placeholder field is removed and `edge_actions.py` is bundled. That report is a correction history, not a live portability blocker. The audit scripts' remaining temporary producer paths are adapted to sibling packet paths, without changing mathematical checks. Source SHA256 manifests remain historical input manifests, not hashes of the adapted edition.

The first checkpoint's original provenance/adaptation files and all packet files remain untouched. Only README.md, SYNTHESIS.md and REPRODUCTION.md are replaced; second-checkpoint provenance files use new names. The new map and replay/preparation receipts remain external review artifacts. No remote change is performed by this preparation. The analytic Dirichlet-series wave is excluded until separately audited.


### Explicit audit hardening and source-boundary controls

The first observer audit originally used removable assertions. Its publication edition replaces those exact AST predicates with RuntimeError guards, so they remain active under optimization. Run any of its three scripts with `--negative-check` to inject an expected guard failure; do this in normal and optimized modes. Before/after hashes and transformation scope are recorded in the adaptations receipt.

The prime-support audit originally required the excluded source PDF and fetched metadata. Its publication edition retains only authored DOI/date/preprint/checksum facts in `audit-prime-support/SOURCE_STATUS.json`. All arithmetic checks still run. It explicitly reports `source_pdf_rechecked: false` unless the user supplies `--source-dir /path/to/local/source` containing MagicHourGlassPrimes.pdf. The original independent source verification remains historical audit evidence, not a claim of fresh verification by the portable replay.
