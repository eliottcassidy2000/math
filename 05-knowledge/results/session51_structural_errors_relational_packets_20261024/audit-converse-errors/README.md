# Independent converse-error audit

Read `REPORT.md` for the mathematical verdict, scope, inherited completeness,
and independent mechanisms. This packet does not publish anything.

Run from any directory:

    python3 /path/to/packet/independent.py > /tmp/audit.json
    python3 -O /path/to/packet/independent.py > /tmp/audit-O.json
    cmp /tmp/audit.json /tmp/audit-O.json
    python3 /path/to/packet/prime_necklaces.py
    python3 /path/to/packet/ports_independent.py
    python3 /path/to/packet/validate.py
    python3 /path/to/packet/validate_candidate.py /path/to/session51-converse-errors
    python3 /path/to/packet/verify_manifest.py

`SOURCE_PINS.json` pins copied mathematical inputs; the independent checker
verifies them before use. `MANIFEST.json` pins every artifact in this audit
packet except itself and Python caches. `candidate-validation.json` pins the
producer's frozen packet and records literal-copy normal/optimized replays,
source stability, mutation rejection, and complete independent matching.
