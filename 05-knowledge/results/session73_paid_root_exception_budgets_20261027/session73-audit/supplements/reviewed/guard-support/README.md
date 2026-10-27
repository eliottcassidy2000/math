# Session73 guard-support audit

Read [REPORT.md](REPORT.md): HYP-9419 clause 2's new support-rank equality is
refuted at order four literally, and at order five even after downward
closure. The separate Session72 sharp-depth conjecture remains OPEN.

- `check.py`: dependency-free exact checker, valid under normal Python and `-O`.
- `check.normal.json`, `check.optimized.json`: identical verified outputs.
- `SOURCE_PINS.json`, `sources/`: immutable incoming source evidence.
- `VERIFICATION.json`, `MANIFEST.json`: execution receipt and content hashes.

Reproduction: `python3 check.py --verify-manifest` (also run with `python3 -O`).
