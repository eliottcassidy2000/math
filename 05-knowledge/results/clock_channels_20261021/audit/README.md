# Independent clock-channel audit

Result: PASS. Read `AUDIT.md` for the mathematical review and scope.

`independent_clock_channels.py` accepts `--packet PATH` pointing to the frozen
producer directory. With the original sibling layout, run:

    python3 independent_clock_channels.py > independent.normal.json
    python3 -O independent_clock_channels.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json

Both retained independent receipts contain 46,083 successful checks. The two
retained producer reruns exactly match its existing 263-check receipts.
No producer code is imported; no external packages or network are required.

Verify this audit's retained files with `sha256sum -c SHA256SUMS`.
