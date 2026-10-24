# Blue internal-port existence probe

This standalone packet separates three statuses:

- **PROVED:** anti-involution-preserving regular completion and the exact
  symmetric-cycle/downhill-path bijection; elementary existence of an
  anti-involution in any self-converse tournament.
- **FINITE-EXACT:** complete fixed-anti-involution controls at even orders
  2, 4 and 6, including realization of every paired arc as a central arc.
- **CONDITIONAL:** the all-order existence conclusion depends on an original
  1980 Salvi-Zagaglia theorem not inspected here. A precise secondary survey
  statement supplies the bibliographic lead. No new Hamiltonicity theorem
  or primary-source verification is claimed.

Read PROBE.md for the proof and exact stopping boundary.

## Reproduce

From this directory:

    python3 verify_manifest.py
    python3 check.py > replay-normal.json
    python3 -O check.py > replay-optimized.json
    cmp normal.json replay-normal.json
    cmp normal.json replay-optimized.json

Python standard library only; no network, atlas, package, sibling folder,
working-directory-specific path, or literature theorem enters check.py.
The two generated replay files are disposable, not manifest payloads.

## Frozen context

sources/balanced-metagraph-theorem.md preserves the inherited extension and
endpoint context. sources/anti-involution-context.md is a located canonical
reference, not an imported classification dependency. Both mathematical
steps needed here are re-proved in PROBE.md. SOURCE_PINS.json records their
original paths and hashes. Source paths in these contextual documents refer
to their original repository locations; the local computation does not read
them or require those transitive references.

The original paper was not downloaded. The secondary URL, printed/PDF page,
formal paraphrase, short quotation and limitations are in PROBE.md. This
packet is not a replacement primary-literature certificate.
