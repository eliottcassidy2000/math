# Reproduce the frozen checkpoints

The newest `FILE_MANIFEST.checkpoint*.json` is the complete current inventory.
Earlier checkpoint manifests, native manifests, code, mathematical notes and
receipts are preserved byte-for-byte; later changes use additive files.
The aggregate manifest itself is committed and externally identified by its
Git blob and publication receipt; it cannot hash its own bytes.

From this directory run:

```bash
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
python3 -I -B -O verify_provenance.py
python3 -I -B test_packet_corruption.py .
python3 -I -B test_provenance_corruption.py .
```

The aggregate checks every declared byte count, SHA256, mode and relative path,
rejects symlinks and unlisted payloads, then copies the packet into a temporary
directory. Only the declared current scripts execute, both normally and under
`-O`. Outputs must match their frozen receipts byte-for-byte. Replays may not
alter any frozen file or add payloads. Repository source programs are inert
snapshots, not replay targets. Source origin checks bind each such snapshot to
its exact commit, original path, SHA256 and Git blob; audit bindings additionally
retain the exact reviewed native bytes. Corruption tests exercise independent
rejection gates rather than accepting a digest-only happy path.

These are exact finite checks and proof-document provenance checks, not a formal
proof assistant build, exhaustive all-order computation, or CI pass. Consult each
native report for its precise universe and open boundaries. Remote publication
receipts separately establish whole-tree identity, exact changed paths and bytes,
ancestry, expected-head non-force updates, and local clean 0/0 synchronization.
No GitHub checks or workflow runs is reported as **none**, never as passing CI.
