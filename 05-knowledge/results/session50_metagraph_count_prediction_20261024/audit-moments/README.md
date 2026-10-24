# Session50 independent finite-moment reconstruction audit

Verdict: **PASS**, subject to the explicit inherited-atlas scope in `AUDIT.md`.

- `check.py`: independent atlas H/Aut/converse verification, fresh fixed-word enumeration n=1,…,6, polynomial arithmetic and adversarial controls.
- `reconstruct.py`: independent exact decoder receiving only n and supplied moments; optimized rational bounds are proved by coefficient identities.
- `simple_certificate.py`: shorter standalone eight-moment decoder using a fixed three-coefficient upper polynomial. It needs neither matrix inversion nor a lower bound.
- `normal.json`, `optimized.json`: deterministic audit results, byte-identical.
- `sources.txt`, `SOURCE_PINS.json`: source paths, pinned commit and SHA-256 provenance. Source snapshots are inert data, never imported code.
- `MANIFEST.json`: every packet file except the manifest itself is byte-pinned. The receipt is deterministic and itself included after the first replay.
- `verify_packet.py`: actual manifest-only fresh-directory normal and optimized replay. `REPLAY.json` records the deterministic result.

Run from this directory:

    python3 check.py --output fresh-normal.json
    python3 -O check.py --output fresh-optimized.json
    cmp fresh-normal.json fresh-optimized.json
    python3 verify_packet.py

All checks use standard-library Python. No network or Git writes are needed. The independent checker never imports producer code; pinned producer code is replayed separately in isolated directories. The source quotient completeness at order eight is inherited; every stored row and all numeric normalizations are independently checked. The complete fresh quotient census here is through order six.

This is statistical reconstruction from already available metagraph masses. The moments in this packet came from a full atlas, so the packet makes no claim to predict an unseen term from n alone or to outperform enumeration.
