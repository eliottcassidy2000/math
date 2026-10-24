# Independent nonrigid regular-eleven audit

AUDIT.md records the complete proof audit, genuinely independent small-family
computation, class-by-class cross-bijection, and subsequent complete-atlas holdout.
The frozen producer is under producer/. The correction audit remains separate.

From any working directory:

    python3 /path/to/packet/verify_audit.py
    python3 -O /path/to/packet/verify_audit.py
    python3 /path/to/packet/replay.py --output /tmp/nonrigid-audit-replay.json

The independent structural generator never reads the complete parent atlas.
Only atlas_control.cpp subsequently uses that atlas as a holdout. No full
regular-eleven census or target count is an input to the structural theorem.
