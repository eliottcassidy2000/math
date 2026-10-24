# Reproduce the independent audit

Python standard library only, no network or external files required.
Run from this folder:

    python3 verify_manifest.py
    python3 audit.py > /tmp/blue-audit-normal.json
    python3 -O audit.py > /tmp/blue-audit-optimized.json
    cmp normal.json /tmp/blue-audit-normal.json
    cmp normal.json /tmp/blue-audit-optimized.json
    python3 producer/check.py > /tmp/blue-producer-normal.json
    python3 -O producer/check.py > /tmp/blue-producer-optimized.json
    cmp producer/normal.json /tmp/blue-producer-normal.json
    cmp producer/normal.json /tmp/blue-producer-optimized.json

The audit verifies full-matrix and full-permutation controls independently.
AUDIT.md reviews proof and source scope. All-order existence is CONDITIONAL,
not established by these finite checks. The producer snapshot is immutable.
