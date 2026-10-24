# Replay

Requires Python 3 with NumPy.

- python audit.py producer
- python -O audit.py producer
- python verify_topology.py producer/normal.json
- python -O verify_topology.py producer/normal.json
- python hostile.py producer
- python -O hostile.py producer

`producer/` is a frozen manifest-listed copy, including every pinned source dependency. `MANIFEST.json` lists the exact standalone closure. `verify_manifest.py` verifies its hashes. Receipts document normal/optimized independent runs and hostile controls. See `AUDIT.md` for mathematical scope and the 455-edge topology certificate.
