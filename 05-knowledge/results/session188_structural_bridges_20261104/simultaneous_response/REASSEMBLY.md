# Restore all exact response vectors and reversal records

The large records are stored losslessly as bounded gzip/base64 chunks. Nothing is omitted: both11,356,449-byte anchor receipts, both complete search receipts, and all five complete response banks are restored. The small summary files are indexes, not substitutes for these records. Every cubic and quartic coordinate, including zeros in each named basis, remains present.

From this directory run:

    python3 reassemble.py
    sha256sum -c manifest.sha256
    python3 search.py > replay.search.json
    cmp receipt.json replay.search.json
    python3 anchored.py > replay.anchored.json
    cmp anchored.json replay.anchored.json

Repeat the two checkers under `python3 -O` for the optimization control. The checkers use the already-published sibling `../response/verify.py`. The separate direct-count audit is `../anchored_audit/check.py` and does not import producer code or read polynomial banks.

`REASSEMBLY.json` binds every chunk and restored file by exact size and SHA-256. The reassembler verifies all groups and any existing outputs before writing anything, rejects corrupted chunks, and refuses to overwrite different bytes. Gzip/base64 is storage packaging only; it changes no finite universe or theorem scope.
