# Independent regular-eleven audit

Read AUDIT.md for findings, proof scope, algorithm independence, and controls.
The frozen production packet is under producer/. No publication is performed.

From any working directory:

    python3 /path/to/packet/verify_audit.py
    python3 -O /path/to/packet/verify_audit.py
    python3 /path/to/packet/replay.py --output /tmp/regular-eleven-audit-replay.json

The replay compiles default, O2, and O3/NDEBUG audit executables and both
production modes, using only files physically copied from MANIFEST.json.
No network, original workspace, target count oracle, or target-eleven
canonicalizer is required. Reports produced by the independent algorithm
and literal degree-sequence recurrence are retained beside their source.
