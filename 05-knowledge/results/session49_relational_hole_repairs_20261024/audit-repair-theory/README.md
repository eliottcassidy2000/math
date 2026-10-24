# Independent Session49 relational repair audit

This packet audits the three new theory packets, not the native 82-slot
numerical receiver census. `AUDIT.md` contains the proof and scope verdict.

Run from any working directory:

    python -I -B verify_manifest.py
    python -I -B independent_check.py --output /tmp/repair-normal.json
    python -I -O -B independent_check.py --output /tmp/repair-optimized.json
    python -I -B replay_producers.py
    python -I -B portable_check.py

Use absolute script paths when outside this directory. Compare generated
receipts with the corresponding frozen `normal.json` and `optimized.json`.
The independent checker imports no producer or repository code and uses no
removable assertions. All scripts use the Python standard library only.

`INPUT_PINS.json` identifies the exact reviewed producer manifests and the
two inherited source snapshots. Entire final producer packets are archived
under `inputs/`, so their executable checks can also be replayed without
the repository. This producer replay is distinct from independent evidence.

`MANIFEST.json` closes every shipped regular file except itself, using raw
byte sizes and SHA256. `portable_check.py` creates a genuinely new temporary
packet from only those listed files, uses an unrelated working directory,
verifies its closure, runs independent normal/optimized checks and all three
producer checks in both modes, compares frozen receipts, and checks closure
again. It is an executed relocation, not a proposed command.

No Git operations were performed by this audit lane. The parent owns
integration and publication.
