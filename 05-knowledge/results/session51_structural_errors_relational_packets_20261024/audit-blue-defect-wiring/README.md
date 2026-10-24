# Independent actual blue-wiring audit

Read `REPORT.md` for the finite-exact PASS, inherited dependencies, and scope.
The packet verifies only the induced blue graph on the 24 positive-error
self-converse order-ten nodes; it does not identify boundary target classes
or recover black wiring.

Reproduce from any current directory:

    python3 /path/to/packet/independent.py > /tmp/blue-audit.json
    python3 -O /path/to/packet/independent.py > /tmp/blue-audit-O.json
    cmp /tmp/blue-audit.json /tmp/blue-audit-O.json
    python3 /path/to/packet/validate.py /path/to/session51-blue-defect-wiring
    python3 /path/to/packet/verify_manifest.py

`SOURCE_PINS.json` authenticates copied inputs; `MANIFEST.json` freezes every
packet file except itself and Python caches. The validator performs actual
copied-directory normal/-O replays and semantic/source mutation tests.
